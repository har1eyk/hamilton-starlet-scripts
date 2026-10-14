#!/usr/bin/env python3
"""Exp1034 single-chain experiment draft and offline, mass-aware dry run.

This standalone entrypoint is offline. The Exp1034 notebook executes the same
plan through a simulated or explicitly approved STARlet session. See README.md.
All volumes are microlitres; concentrations are ng/microlitre; heights are mm.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import inspect
import json
import math
import platform
import subprocess
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pylabrobot
from pylabrobot.molecular_devices.filtermax import FilterMaxF5
from pylabrobot.resources.volume_tracker import VolumeTracker

from exp1034_deck import build_deck, choose_plate_tip

HERE = Path(__file__).resolve().parent
SOURCE_PAGE = "https://chatgpt.com/space/page_2505a363d96081918edc799c4d38a292"
VOLUMES = (200, 100, 75, 50)
CONCENTRATIONS = {"A31": 0, "A30": 16, "A29": 8, "A28": 4, "A27": 2, "A26": 1, "A25": 0.5}
INITIAL = {"A32": 4000, "A31": 4000, "A30": 3000, "A29": 0, "A28": 0, "A27": 0, "A26": 0, "A25": 0}
DNA_SOURCES = ("A30", "A29", "A28", "A27", "A26", "A25")
# Explicit mapping documented by Exp1031's physical_to_site, NOT inferred geometry.
SITE_CANDIDATES = {"A32": 0, "A31": 1, "A30": 2, "A29": 3, "A28": 4, "A27": 5, "A26": 6, "A25": 7}


def require(condition: bool, message: str) -> None:
  """Raise a validation error even when Python is run with optimization."""
  if not condition:
    raise ValueError(message)


def write_json(path: Path, data: Any) -> None:
  """Write strict JSON, refusing NaN and infinity."""
  path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def write_csv(path: Path, records: list[dict]) -> None:
  """Write a nonempty sequence of records with a header."""
  with path.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(records[0]))
    writer.writeheader()
    writer.writerows(records)


@dataclass(frozen=True)
class WellRecord:
  """One immutable well assignment, with preparation identity explicit."""

  well: str
  source_tube: str
  final_concentration: float
  final_volume: int
  dna_type: str
  replicate: int
  replicate_type: str
  preparation_id: str


def manifest() -> list[WellRecord]:
  """Expand the supplied map without changing its well or replicate order."""
  records = []
  for volume, concentration, source, wells in json.loads(
    (HERE / "well_assignments.json").read_text()
  ):
    for replicate, well in enumerate(wells, 1):
      records.append(
        WellRecord(
          well,
          source,
          concentration,
          volume,
          "plasmid DNA" if concentration else "buffer only",
          replicate,
          "dispensing",
          "single_dilution_chain" if concentration else "matching_buffer",
        )
      )
  expected = {f"{row}{column}" for row in "ABCDEFGH" for column in range(1, 13)}
  require(
    len(records) == 96 and {r.well for r in records} == expected,
    "The manifest must contain every well exactly once.",
  )
  for volume in VOLUMES:
    group = [r for r in records if r.final_volume == volume]
    for source, concentration in CONCENTRATIONS.items():
      subgroup = [r for r in group if r.source_tube == source]
      count = 6 if concentration == 0 else 3
      require(len(subgroup) == count, f"Wrong count for {volume}/{source}.")
      require(
        [r.replicate for r in subgroup] == list(range(1, count + 1)), "Replicate order changed."
      )
      require(all(r.final_concentration == concentration for r in subgroup), "Source mismatch.")
    require(Counter(r.well[0] for r in group) == {r: 3 for r in "ABCDEFGH"}, "Row balance changed.")
    require(
      Counter(int(r.well[1:]) for r in group) == {c: 2 for c in range(1, 13)},
      "Column balance changed.",
    )
  return records


def dispensing_order(records: list[WellRecord]) -> list[WellRecord]:
  """Group wells by source and rotate volumes within each group, preserving replicate labels."""
  ordered = []
  for source_index, source in enumerate(CONCENTRATIONS):
    for replicate in range(1, 7 if source == "A31" else 4):
      shift = (replicate - 1 + source_index) % 4
      volumes = VOLUMES[shift:] + VOLUMES[:shift]
      for volume in volumes:
        ordered.extend(
          r
          for r in records
          if r.source_tube == source and r.replicate == replicate and r.final_volume == volume
        )
  require(
    len(ordered) == len(records) and len({r.well for r in ordered}) == len(records),
    "Dispense order lost or duplicated wells.",
  )
  return ordered


def environment() -> dict:
  """Record the actual interpreter, package source and source revision."""
  checkout = Path(pylabrobot.__file__).resolve().parent.parent
  revision = subprocess.run(
    ["git", "-C", str(checkout), "rev-parse", "HEAD"], capture_output=True, text=True, check=False
  )
  return {
    "python": platform.python_version(),
    "pylabrobot_version": pylabrobot.__version__,
    "pylabrobot_file": str(pylabrobot.__file__),
    "pylabrobot_commit": revision.stdout.strip(),
    "filtermax_backend": inspect.getfile(FilterMaxF5),
    "filtermax_read_signature": str(inspect.signature(FilterMaxF5.read_absorbance)),
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "assignment_sha256": hashlib.sha256((HERE / "well_assignments.json").read_bytes()).hexdigest(),
  }


def inspect_resources() -> dict:
  """Construct the corrected deck without connecting to hardware."""
  return build_deck().audit


def hardware_blockers(resources: dict) -> list[str]:
  """List remaining method inputs, without reopening the user's geometry corrections."""
  return [
    "Enter measured A30 concentration (16 ng/uL in 3000 uL total).",
    "Record method validation for immersion, noncontact delivery, mixing, tip reach and source residual/loss limits.",
    "Approve each physical run in the notebook before setup; dry-run success is not physical validation.",
  ]


class DryRun:
  """Execute the full plan on PLR volume trackers, with a separate DNA mass ledger.

  A residual floor is dead volume plus a conservative TOTAL loss allowance per
  source. Nominal transfers conserve mass; losses are reserved, not fabricated.
  Empty receiving tubes may start below that floor; every aspiration must end above it.
  """

  def __init__(
    self,
    concentration: float,
    resources: dict,
    limits: dict,
    mix_volume: float = 750,
    mix_cycles: int = 8,
    air_gap: float = 0,
    immersion: float = 2,
    rim_margin: float = 2,
  ) -> None:
    """Initialize an entirely in-memory experiment."""
    require(
      math.isfinite(concentration) and concentration > 0, "Starting concentration must be positive."
    )
    require(
      math.isclose(concentration, 16, rel_tol=0, abs_tol=1e-9),
      f"A30 is {concentration:g}, not 16 ng/uL. Actual chain would be "
      f"{[concentration / 2**i for i in range(6)]}; targets must not be silently relabeled.",
    )
    require(
      math.isfinite(mix_volume) and mix_volume > 0 and mix_cycles > 0, "Invalid mixing inputs."
    )
    require(math.isfinite(air_gap) and air_gap >= 0, "Invalid air-gap capacity reserve.")
    require(set(limits) == set(INITIAL), "Supply limits for exactly A32 through A25.")
    self.floor = {}
    self.limit_status = {}
    for source, item in limits.items():
      dead, loss = item["dead_volume"], item["total_loss_allowance"]
      for value in (dead, loss):
        require(value is None or (math.isfinite(value) and value >= 0), "Invalid source allowance.")
      self.floor[source] = (dead or 0) + (loss or 0)
      self.limit_status[source] = (
        "supplied_assumption" if dead is not None and loss is not None else "unverified"
      )
    self.trackers = {source: VolumeTracker(thing=source, max_volume=5000) for source in INITIAL}
    for source, volume in INITIAL.items():
      self.trackers[source].set_volume(volume)
    self.mass = {source: concentration * 3000 if source == "A30" else 0.0 for source in INITIAL}
    self.wells: dict[str, dict] = {}
    self.events: list[dict] = []
    self.phase = "preparation"
    self.snapshots: list[dict] = []
    self.tips = Counter()
    self.tip: tuple[int, str] | None = None
    self.contacted_source: str | None = None
    self.resources = resources
    self.mix_volume = mix_volume
    self.mix_cycles = mix_cycles
    self.air_gap = air_gap
    self.immersion = immersion
    self.rim_margin = rim_margin
    self.loss = {source: item["total_loss_allowance"] or 0 for source, item in limits.items()}

  def volumes(self) -> dict[str, float]:
    """Return a fresh volume snapshot."""
    return {source: tracker.get_used_volume() for source, tracker in self.trackers.items()}

  def log(self, action: str, **details: Any) -> None:
    """Record a virtual event; its index is not an actual execution timestamp."""
    self.events.append(
      {"step": len(self.events) + 1, "action": action, "channel": 0, "phase": self.phase, **details}
    )

  def snapshot(self, stage: str) -> None:
    """Retain a dilution-stage ledger row."""
    self.snapshots.append({"stage": stage, **self.volumes()})

  def pickup(self, size: int) -> None:
    """Allocate one fresh tip in column-major rack order."""
    require(self.tip is None, "A tip is already attached.")
    index = self.tips[size]
    require(index < self.resources["tips"][str(size)]["available"], "Tip inventory exhausted.")
    spot = f"{'ABCDEFGH'[index % 8]}{index // 8 + 1}"
    self.tip = (size, spot)
    self.tips[size] += 1
    self.contacted_source = None
    self.log("pick_up_tip", tip_size=size, rack_spot=spot)

  def discard(self) -> None:
    """Discard the current tip; never return it to the rack."""
    require(self.tip is not None, "No tip to discard.")
    self.log("discard_tip", tip_size=self.tip[0], rack_spot=self.tip[1])
    self.tip = None
    self.contacted_source = None

  def capacity(self, volume: float) -> None:
    """Check commanded liquid plus air against nominal and defined tip capacity."""
    require(self.tip is not None, "No attached tip.")
    tip = self.resources["tips"][str(self.tip[0])]
    require(
      0 < volume and volume + self.air_gap <= min(tip["nominal"], tip["maximum"]),
      f"Liquid {volume} + air {self.air_gap} exceeds tip capacity.",
    )

  def withdraw(self, source: str, volume: float) -> float:
    """Withdraw nominal liquid, checking residual plus the reserved loss budget."""
    self.capacity(volume)
    require(
      self.contacted_source in (None, source), "Destination-contacted tip returned to source."
    )
    before = self.trackers[source].get_used_volume()
    require(
      before - volume >= self.floor[source] - 1e-9,
      f"{source}: {before - volume:g} uL residual is below dead+loss floor {self.floor[source]:g}.",
    )
    concentration = self.mass[source] / before
    self.trackers[source].remove_liquid(volume)
    self.trackers[source].commit()
    self.mass[source] -= concentration * volume
    # Buffer-only reuse is allowed. DNA exposure binds this tip to its source.
    if concentration > 0:
      self.contacted_source = source
    return concentration * volume

  def transfer(
    self, source: str, destination: str, volume: float, *, contact: bool = False
  ) -> None:
    """Deliver noncontact, except for contact prefills into DNA-free dilution tubes."""
    if contact:
      require(
        self.phase == "preparation"
        and source in ("A32", "A31")
        and destination in DNA_SOURCES[1:]
        and self.mass[source] == 0
        and self.mass[destination] == 0
        and self.contacted_source is None,
        "Contact buffer prefills require a clean tip and DNA-free source and destination.",
      )
    before = self.volumes()
    mass = self.withdraw(source, volume)
    if destination in self.trackers:
      self.trackers[destination].add_liquid(volume)
      self.trackers[destination].commit()
      self.mass[destination] += mass
      after_destination = self.trackers[destination].get_used_volume()
    else:
      require(destination not in self.wells, f"Well {destination} filled twice.")
      self.wells[destination] = {"volume": volume, "mass": mass, "source": source}
      after_destination = volume
    self.log(
      "transfer",
      source=source,
      destination=destination,
      volume=volume,
      source_before=before[source],
      source_after=self.volumes()[source],
      destination_before=before.get(destination, 0),
      destination_after=after_destination,
      tip_size=self.tip[0],
      noncontact=not contact,
      buffer_prefill=contact,
      air_capacity_reserve=self.air_gap,
      aspiration_height=None,
      dispense_height=None,
      geometry_status="unverified; no coordinates or hardware commands emitted",
    )

  def mix(self, source: str) -> None:
    """Model repeated partial-volume mixing without assuming proven homogeneity."""
    self.capacity(self.mix_volume)
    last_destination = self.events[-1].get("destination") if self.events else None
    require(
      self.contacted_source in (None, source) or last_destination == source,
      "Mixing would introduce a contaminated tip into an unrelated source.",
    )
    self.contacted_source = source
    before = self.volumes()[source]
    for cycle in range(self.mix_cycles):
      mass = self.withdraw(source, self.mix_volume)
      self.trackers[source].add_liquid(self.mix_volume)
      self.trackers[source].commit()
      self.mass[source] += mass
      self.log(
        "mix_cycle",
        source=source,
        volume=self.mix_volume,
        cycle=cycle + 1,
        transient_source_volume=before - self.mix_volume,
        source_after=before,
        aspiration_height=None,
        dispense_height=None,
        validation="homogeneity pending",
      )

  def refill(self, source: str, acknowledged: bool) -> None:
    """Change tracking only after a simulated acknowledgement; require a clean tip."""
    require(
      source in ("A32", "A31") and self.mass[source] == 0, "Only clean buffer may be refilled."
    )
    require(self.contacted_source is None, "DNA-exposed tip present at refill.")
    before = self.volumes()[source]
    self.log(
      "operator_refill_checkpoint",
      source=source,
      before=before,
      target_total=4000,
      nominal_addition=4000 - before,
      acknowledgement="SIMULATED" if acknowledged else "NOT_ACKNOWLEDGED",
    )
    require(acknowledged, "Refill not acknowledged; tracking was not updated.")
    self.trackers[source].set_volume(4000)

  def plate(self, records: list[WellRecord]) -> None:
    """Fill the supplied wells with one tip per source, sized for the group's lowest level."""
    self.phase = "plate"
    ordered = dispensing_order(records)
    for source in CONCENTRATIONS:
      group = [record for record in ordered if record.source_tube == source]
      if not group:
        continue
      residual = (
        self.volumes()[source] - sum(record.final_volume for record in group) - self.loss[source]
      )
      self.pickup(choose_plate_tip(residual, self.resources, self.immersion, self.rim_margin))
      for record in group:
        self.transfer(source, record.well, record.final_volume)
      self.discard()

  def run_prepared(self, source_volumes: dict[str, float]) -> dict:
    """Plan all wells from operator-supplied prepared volumes, without dilution or refill."""
    require(
      set(source_volumes) == set(CONCENTRATIONS),
      "Enter current total volumes for exactly A31, A30, A29, A28, A27, A26 and A25.",
    )
    records = manifest()
    for source, volume in source_volumes.items():
      require(
        isinstance(volume, (int, float))
        and not isinstance(volume, bool)
        and math.isfinite(volume)
        and 0 < volume <= 5000,
        f"Enter the current total volume in {source}, in uL (greater than 0, at most 5000).",
      )
      needed = sum(r.final_volume for r in records if r.source_tube == source)
      require(
        volume >= needed + self.floor[source],
        f"{source} needs at least {needed + self.floor[source]:g} uL TOTAL "
        f"({needed:g} for the plate plus {self.floor[source]:g} reserved); entered {volume:g}.",
      )
    for source in INITIAL:
      volume = source_volumes.get(source, 0)  # A32 is unused and untracked in this mode.
      self.trackers[source].set_volume(volume)
      self.mass[source] = CONCENTRATIONS.get(source, 0) * volume
    initial_mass = sum(self.mass.values())
    self.snapshot("operator_supplied_prepared_tubes")
    self.plate(records)
    self.snapshot("plate_complete")
    for record in records:
      well = self.wells[record.well]
      require(
        well["volume"] == record.final_volume
        and well["source"] == record.source_tube
        and math.isclose(well["mass"] / well["volume"], record.final_concentration),
        f"Incorrect prepared-source plan for {record.well}.",
      )
    require(self.tips == {300: 7}, "Expected seven 300 uL plate tips.")
    require(
      math.isclose(
        sum(self.mass.values()) + sum(w["mass"] for w in self.wells.values()), initial_mass
      ),
      "Prepared-source DNA mass is not conserved.",
    )
    return {
      "status": "PASS: nominal offline checks; NOT physical validation",
      "hardware_ready": False,
      "entry_phase": "plate",
      "initial_source_volumes": dict(source_volumes),
      "prepared_concentrations": dict(CONCENTRATIONS),
      "unused_sources": ["A32"],
      "wells": len(self.wells),
      "whole_plate_volume": sum(w["volume"] for w in self.wells.values()),
      "tip_counts": dict(self.tips),
      "plate_tip_policy": "one_tip_per_source",
      "source_residuals": self.volumes(),
      "dead_plus_loss_floors": self.floor,
      "limit_status": self.limit_status,
      "refills": [],
      "mixing": {"cycles": 0, "flow": None, "homogeneity_validated": False},
      "preparation": "Supplied by operator; not performed in this session.",
    }

  def run(self) -> dict:
    """Simulate prefilling, five serial links, refill and all 96 final fills."""
    records = manifest()
    self.snapshot("initial")
    self.pickup(1000)
    for destination in DNA_SOURCES[1:]:
      for buffer in ("A32", "A31"):
        if self.volumes()[buffer] - 750 < self.floor[buffer]:
          self.refill(buffer, acknowledged=True)
        self.transfer(buffer, destination, 750, contact=True)
    self.snapshot("buffer_prefill_complete")
    self.mix("A30")
    for index, (source, destination) in enumerate(zip(DNA_SOURCES, DNA_SOURCES[1:])):
      if index:
        self.pickup(1000)
      self.transfer(source, destination, 750)
      self.transfer(source, destination, 750)
      self.mix(destination)
      self.discard()
      self.snapshot(f"mixed_{destination}")
    require(
      self.volumes()
      == {
        "A32": self.volumes()["A32"],
        "A31": self.volumes()["A31"],
        **{s: 1500 for s in DNA_SOURCES[:-1]},
        "A25": 3000,
      },
      "Incorrect chain residuals.",
    )
    self.phase = "refill"
    self.refill("A31", acknowledged=True)
    self.snapshot("A31_refilled_to_4000_total")
    self.plate(records)
    self.snapshot("plate_complete")
    for record in records:
      measured = self.wells[record.well]
      require(measured["volume"] == record.final_volume, "Well volume mismatch.")
      require(measured["source"] == record.source_tube, "Well source mismatch.")
      require(
        math.isclose(measured["mass"] / measured["volume"], record.final_concentration),
        "Well concentration mismatch.",
      )
    require(
      sum(self.tips.values()) == 12, "Expected five preparation tips plus seven source-group tips."
    )
    require(
      math.isclose(sum(self.mass.values()) + sum(w["mass"] for w in self.wells.values()), 48000),
      "DNA mass is not conserved.",
    )
    return {
      "status": "PASS: nominal offline checks; NOT physical validation",
      "hardware_ready": False,
      "wells": len(self.wells),
      "whole_plate_volume": sum(w["volume"] for w in self.wells.values()),
      "tip_counts": dict(self.tips),
      "plate_tip_policy": "one_tip_per_source",
      "source_residuals": self.volumes(),
      "dead_plus_loss_floors": self.floor,
      "limit_status": self.limit_status,
      "nominal_loss": 0,
      "loss_interpretation": "Supplied total losses are reserved, not simulated as exact losses.",
      "reagent_budget": {
        "initial_buffer": 8000,
        "dilution_buffer": 7500,
        "blank_buffer": 2550,
        "total_buffer_delivered": 10050,
        "minimum_initial_buffer_shortfall": 2050,
        "nominal_refill_buffer": sum(
          e["nominal_addition"] for e in self.events if e["action"] == "operator_refill_checkpoint"
        ),
        "initial_DNA_solution": 3000,
        "initial_DNA_mass_ng": 48000,
        "plate_DNA_solution": 7650,
        "per_concentration_plate_solution": 1275,
      },
      "refills": [e for e in self.events if e["action"] == "operator_refill_checkpoint"],
      "mixing": {
        "volume": self.mix_volume,
        "cycles": self.mix_cycles,
        "flow": None,
        "homogeneity_validated": False,
      },
      "air_gap_capacity_reserve": self.air_gap,
      "checks": [
        "96 unique wells",
        "exact concentration/source mapping",
        "replicate labels",
        "row/column balance",
        "tip capacity including air reserve",
        "tip inventory",
        "one plate tip per source; blanks and DNA concentrations never share a tip",
        "no source overdraft",
        "supplied residual floors on every aspiration",
        "contact buffer prefills only into DNA-free dilution tubes",
        "noncontact DNA split transfers before mixing and noncontact plate dispensing",
        "single-channel operations",
        "DNA mass conservation",
        "refill acknowledgement before tracking update",
      ],
      "hardware_blockers": hardware_blockers(self.resources),
    }


def create_bundle(args: argparse.Namespace) -> None:
  """Write an immutable offline verification bundle without any device setup."""
  resources = inspect_resources()
  limits = (
    json.loads(args.limits.read_text())
    if args.limits
    else {source: {"dead_volume": None, "total_loss_allowance": None} for source in INITIAL}
  )
  experiment = DryRun(
    args.a30_concentration, resources, limits, args.mix_volume, args.mix_cycles, args.air_gap
  )
  report = experiment.run()
  args.output.mkdir(parents=True, exist_ok=False)
  records = manifest()
  write_csv(args.output / "well_manifest.csv", [asdict(r) for r in records])
  write_csv(args.output / "dilution_ledger.csv", experiment.snapshots)
  write_json(args.output / "operation_plan.json", experiment.events)
  write_json(args.output / "verification.json", report)
  write_json(args.output / "resource_audit.json", resources)
  write_json(
    args.output / "metadata.json",
    {
      "experiment": "Exp1034",
      "mode": "OFFLINE_DRY_RUN",
      "source_page": SOURCE_PAGE,
      "generated_at": datetime.now(timezone.utc).isoformat(),
      "environment": environment(),
      "A30_concentration_input": args.a30_concentration,
      "A30_total_volume": 3000,
      "input_status": "planning inputs, not independently measured",
      "source_limits": limits,
      "preparation": "one serial dilution chain; all sample replicates are dispensing replicates",
      "units": {"volume": "uL", "concentration": "ng/uL", "mass": "ng", "height": "mm"},
      "actual_run_started_at": None,
      "plate_completed_at": None,
      "reader_settings": None,
      "reader_instrument_identity": None,
      "pathlength_correction_applied": None,
    },
  )
  write_csv(
    args.output / "raw_a260_template.csv",
    [
      {
        "well": r.well,
        "a260": "",
        "valid": "",
        "exclusion_reason": "",
        "read_id": "",
        "timestamp": "",
      }
      for r in records
    ],
  )
  print(
    json.dumps(
      {
        "output": str(args.output.resolve()),
        "status": report["status"],
        "tips": report["tip_counts"],
        "residuals": report["source_residuals"],
        "hardware_blocker_count": len(report["hardware_blockers"]),
      },
      indent=2,
    )
  )


def main() -> None:
  """Run only offline planning; supplying a concentration is mandatory."""
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument(
    "--a30-concentration", type=float, required=True, help="Explicit ng/uL input; must equal 16."
  )
  parser.add_argument(
    "--output", type=Path, required=True, help="New output directory; never overwrites a run."
  )
  parser.add_argument(
    "--limits", type=Path, help="Per-source dead volume and total loss allowances JSON."
  )
  parser.add_argument(
    "--mix-volume", type=float, default=750, help="Draft mixing volume, uL; validation pending."
  )
  parser.add_argument(
    "--mix-cycles", type=int, default=8, help="Draft cycle count; homogeneity validation pending."
  )
  parser.add_argument(
    "--air-gap", type=float, default=0, help="Capacity reserve only; not a hardware setting."
  )
  create_bundle(parser.parse_args())


if __name__ == "__main__":
  main()
