"""STARlet execution of the Exp1034 plan, with the same path used in simulation."""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import math
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from pylabrobot.legacy.liquid_handling import LiquidHandler
from pylabrobot.legacy.liquid_handling.backends import LiquidHandlerChatterboxBackend, STARBackend
from pylabrobot.legacy.liquid_handling.liquid_classes.hamilton import get_star_liquid_class
from pylabrobot.resources import Container, Liquid, set_tip_tracking, set_volume_tracking

from exp1034 import (
  CONCENTRATIONS,
  DryRun,
  INITIAL,
  environment,
  manifest,
  require,
  write_csv,
  write_json,
)
from exp1034_deck import DeckSetup, build_deck, tip_insertion_limit


class RefillNotAcknowledged(ValueError):
  """Pause at an unchanged manual checkpoint without latching a hardware failure."""


@dataclass(frozen=True)
class Settings:
  """Run inputs in PLR units, with optional prior method notes for exploratory runs."""

  simulate: bool = True
  channel: int = 0
  a30_concentration: float = 16.0
  measured_a30_concentration: float | None = None
  operator: str = ""
  method_validation_reference: str = ""  # Optional prior method record; blank for a first run.
  source_dead_volume: float = 200.0
  source_total_loss_allowance: float = 25.0
  mix_volume: float = 750.0
  mix_cycles: int = 8
  mix_height: float = 6.0
  mix_flow: float = 400.0
  preparation_aspirate_flow: float | None = 250.0
  preparation_dispense_flow: float | None = 250.0
  preparation_aspiration_lowering: float = 0.0
  preparation_minimum_tip_height: float = 4.3
  buffer_prefill_submersion: float = 1.0
  aspirate_flow: float = 100.0
  dispense_flow: float = 120.0
  plate_aspiration_lowering: float = 0.0
  plate_dispense_lowering: float = 0.0
  plate_dispense_height: float | None = None  # Fixed inner-bottom height overrides lowering.
  immersion: float = 2.0
  minimum_tip_height: float = 1.0
  noncontact_clearance: float = 3.0
  tip_rim_margin: float = 2.0
  traverse_height: float = 247.6
  settling_time: float = 1.0


def validate_settings(settings: Settings, physical: bool = False) -> None:
  """Check candidate parameters, measured starting concentration and the physical operator."""
  require(type(settings.simulate) is bool, "simulate must be a boolean.")
  require(type(settings.channel) is int and 0 <= settings.channel < 8, "Invalid channel index.")
  require(
    type(settings.mix_cycles) is int and 1 <= settings.mix_cycles < 10,
    "number_of_tube_mixes must be an integer from 1 to 9.",
  )
  positive = [
    settings.a30_concentration,
    settings.mix_volume,
    settings.mix_height,
    settings.mix_flow,
    settings.buffer_prefill_submersion,
    settings.preparation_minimum_tip_height,
    settings.aspirate_flow,
    settings.dispense_flow,
    settings.immersion,
    settings.minimum_tip_height,
    settings.noncontact_clearance,
    settings.tip_rim_margin,
    settings.traverse_height,
  ]
  nonnegative = [
    settings.source_dead_volume,
    settings.source_total_loss_allowance,
    settings.settling_time,
    settings.preparation_aspiration_lowering,
    settings.plate_aspiration_lowering,
    settings.plate_dispense_lowering,
  ]
  require(
    all(math.isfinite(x) and x > 0 for x in positive), "Nonfinite/nonpositive method parameter."
  )
  require(
    all(math.isfinite(x) and x >= 0 for x in nonnegative), "Invalid residual/loss/time setting."
  )
  require(
    settings.plate_dispense_height is None
    or (
      math.isfinite(settings.plate_dispense_height)
      and settings.plate_dispense_height >= settings.minimum_tip_height
    ),
    "Fixed plate dispense height must be finite and at least the minimum bottom clearance.",
  )
  require(
    all(
      value is None or (math.isfinite(value) and value > 0)
      for value in (settings.preparation_aspirate_flow, settings.preparation_dispense_flow)
    ),
    "Preparation flow rates must be None for liquid-class defaults or finite positive values.",
  )
  if physical:
    issues = []
    if settings.measured_a30_concentration is None or not math.isclose(
      settings.measured_a30_concentration, 16, abs_tol=1e-9, rel_tol=0
    ):
      issues.append(
        "Physical run requires a measured A30 concentration of 16 ng/uL. "
        f"Current measured_a30_concentration={settings.measured_a30_concentration!r}; "
        f"a30_concentration={settings.a30_concentration!r} is only the design target. "
        "Enter the actual measurement in measured_a30_concentration."
      )
    if not settings.operator.strip():
      issues.append("operator is blank. Enter the operator for the physical preparation record.")
    require(
      not issues,
      "Physical run inputs need attention:\n- "
      + "\n- ".join(issues)
      + "\nAfter editing, rerun the full Settings cell to recreate session before calling start().",
    )


def resolve_flow_rates(event: dict, audit: dict, settings: Settings) -> tuple[float, float]:
  """Resolve explicit preparation rates, using only liquid-class flow defaults when set to None."""
  if event["action"] == "mix_cycle":
    return settings.mix_flow, settings.mix_flow
  if event["phase"] != "preparation":
    return settings.aspirate_flow, settings.dispense_flow
  tip = audit["tips"][str(event["tip_size"])]
  aspiration = settings.preparation_aspirate_flow
  dispense = settings.preparation_dispense_flow
  if aspiration is None:
    liquid_class = get_star_liquid_class(
      tip_volume=tip["maximum"],
      is_core=False,
      is_tip=True,
      has_filter=tip["filtered"],
      liquid=Liquid.WATER,
      jet=False,
      blow_out=False,
    )
    require(liquid_class is not None, "No preparation aspiration default for the loaded tip.")
    aspiration = liquid_class.aspiration_flow_rate
  if dispense is None:
    liquid_class = get_star_liquid_class(
      tip_volume=tip["maximum"],
      is_core=False,
      is_tip=True,
      has_filter=tip["filtered"],
      liquid=Liquid.WATER,
      jet=event["noncontact"],
      blow_out=True,
    )
    require(liquid_class is not None, "No preparation dispense default for the loaded tip/mode.")
    dispense = liquid_class.dispense_flow_rate
  require(
    all(math.isfinite(value) and value > 0 for value in (aspiration, dispense)),
    "Invalid preparation liquid-class flow defaults.",
  )
  return aspiration, dispense


def geometry_for_plan(setup: DeckSetup, plan: list[dict], settings: Settings) -> list[dict]:
  """Resolve command heights, tip reach and flow rates before any hardware session begins."""
  resolved = []
  tip_size = None
  for event in plan:
    item = dict(event)
    item["channel"] = settings.channel
    if item["action"] == "pick_up_tip":
      tip_size = item["tip_size"]
    elif item["action"] == "discard_tip":
      tip_size = None
    elif item["action"] in ("transfer", "mix_cycle"):
      tube = setup.tubes[item["source"]]
      loss = settings.source_total_loss_allowance
      if item["action"] == "transfer":
        before, after = item["source_before"], item["source_after"]
        start = tube.compute_height_from_volume(before - loss) - settings.immersion
        end = tube.compute_height_from_volume(after - loss) - settings.immersion
        if item["phase"] == "plate" or item["destination"] in setup.tubes:
          lowering = (
            settings.plate_aspiration_lowering
            if item["phase"] == "plate"
            else settings.preparation_aspiration_lowering
          )
          minimum = max(
            settings.minimum_tip_height,
            tube.get_size_z()
            - tube.material_z_thickness
            - tip_insertion_limit(tip_size, setup.audit)
            + settings.tip_rim_margin,
          )
          require(
            end >= minimum,
            "Unadjusted aspiration exceeds the bottom-clearance or checked tip insertion limit.",
          )
          if item["phase"] != "plate":
            minimum = max(minimum, settings.preparation_minimum_tip_height)
            require(
              minimum < tube.compute_height_from_volume(after - loss),
              "Preparation aspiration floor must remain below the modeled final liquid surface.",
            )
          requested_start = start - lowering
          requested_end = end - lowering
          lowered_start = max(minimum, requested_start)
          lowered_end = max(minimum, requested_end)
          item.update(
            aspiration_lowering_requested=lowering,
            aspiration_start_lowering_applied=start - lowered_start,
            aspiration_end_lowering_applied=end - lowered_end,
            aspiration_minimum_height=minimum,
            aspiration_height_limited=requested_end < minimum,
          )
          start, end = lowered_start, lowered_end
        item["surface_following_distance"] = start - end
        if item["destination"] in setup.tubes:
          destination: Container = setup.tubes[item["destination"]]
        else:
          destination = setup.plate.get_item(item["destination"])
        final_surface = destination.compute_height_from_volume(item["destination_after"])
        if item.get("buffer_prefill", False):
          dispense = final_surface - settings.buffer_prefill_submersion
        else:
          dispense = final_surface + settings.noncontact_clearance
        if item["phase"] == "plate" and item["destination"] not in setup.tubes:
          if settings.plate_dispense_height is not None:
            dispense = settings.plate_dispense_height
            item["dispense_height_mode"] = "fixed_above_inner_bottom"
            item["dispense_lowering_applied"] = 0.0
          else:
            dispense -= settings.plate_dispense_lowering
            item["dispense_height_mode"] = "final_surface_with_offset"
            item["dispense_lowering_applied"] = settings.plate_dispense_lowering
          item["noncontact"] = dispense > final_surface + 1e-9
        require(
          settings.minimum_tip_height
          <= dispense
          < destination.get_size_z() - destination.material_z_thickness,
          f"Dispense height exits {destination.name} or violates bottom clearance.",
        )
        item["destination_final_surface"] = final_surface
        item["dispense_height"] = dispense
        destination_bottom = (
          destination.get_location_wrt(setup.deck).z + destination.material_z_thickness
        )
        item["destination_inner_bottom_z"] = destination_bottom
        item["dispense_tip_z"] = destination_bottom + dispense
        require(
          destination_bottom + dispense + settings.tip_rim_margin < settings.traverse_height,
          "Dispense exceeds the configured traverse clearance.",
        )
      else:
        before, after = item["source_after"], item["transient_source_volume"]
        start = end = settings.mix_height
        require(
          end + settings.immersion <= tube.compute_height_from_volume(after - loss),
          "Mixing tip would lose immersion at its transient low volume.",
        )
        item["surface_following_distance"] = 0.0
        item["dispense_height"] = settings.mix_height
      require(
        end >= settings.minimum_tip_height, "Aspiration would cross the bottom-clearance limit."
      )
      tip = setup.audit["tips"][str(tip_size)]
      insertion_limit = tip_insertion_limit(tip_size, setup.audit)
      penetration = tube.get_size_z() - tube.material_z_thickness - end
      require(
        penetration + settings.tip_rim_margin <= insertion_limit,
        f"{tip_size} uL tip exceeds the checked insertion limit in {item['source']}.",
      )
      require(item["volume"] <= min(tip["nominal"], tip["maximum"]), "Tip volume exceeded.")
      item.update(
        {
          "tip_size": tip_size,
          "aspiration_height": start,
          "aspiration_end_height": end,
          "source_surface_before": tube.compute_height_from_volume(before),
          "source_surface_after": tube.compute_height_from_volume(after),
          "tip_penetration_at_end": penetration,
          "geometry_status": "resolved from resource curve and configured allowances; method validation required",
        }
      )
      aspiration_flow, dispense_flow = resolve_flow_rates(item, setup.audit, settings)
      item.update(
        aspirate_flow_rate=aspiration_flow,
        dispense_flow_rate=dispense_flow,
        # Height corrections do not change the programmed jet/blow-out mode.
        dispense_jet=item["action"] == "transfer"
        and event.get("dispense_jet", event["noncontact"]),
        dispense_blow_out=item["action"] == "transfer",
      )
    resolved.append(item)
  tallest = max(t.get_absolute_location(z="t").z for t in setup.tubes.values())
  require(
    tallest + settings.tip_rim_margin < settings.traverse_height, "Tube rim exceeds safe traverse."
  )
  return resolved


class ExperimentSession:
  """Run once, stop on error, and log intent/completion around each liquid-handling action."""

  def __init__(
    self,
    settings: Settings,
    output: Path,
    confirm: Callable[[str], str] = input,
    *,
    prepared_source_volumes: dict[str, float] | None = None,
  ) -> None:
    """Plan offline; supplied prepared volumes select plate-only execution on an empty plate."""
    validate_settings(settings)
    self.settings = settings
    self.output = output
    self.confirm = confirm
    self.setup = build_deck()
    self.limits = {
      s: {
        "dead_volume": settings.source_dead_volume,
        "total_loss_allowance": settings.source_total_loss_allowance,
      }
      for s in INITIAL
    }
    dry = DryRun(
      settings.a30_concentration,
      self.setup.audit,
      self.limits,
      settings.mix_volume,
      settings.mix_cycles,
      immersion=settings.immersion,
      rim_margin=settings.tip_rim_margin,
    )
    self.entry_phase = "preparation" if prepared_source_volumes is None else "plate"
    if prepared_source_volumes is None:
      self.initial_volumes = dict(INITIAL)
      self.report = dry.run()
      self.report["mixing"]["flow"] = settings.mix_flow
    else:
      self.report = dry.run_prepared(prepared_source_volumes)
      self.initial_volumes = {"A32": 0, **prepared_source_volumes}
    self.plan = geometry_for_plan(self.setup, dry.events, settings)
    self.ledger = dry.snapshots
    self.lh: LiquidHandler | None = None
    self.started = False
    self.failed = False
    self.stopped = False
    self.next_step = 0
    self.attached_tip: tuple[int, str] | None = None
    self.completed_phases: list[str] = []
    self.events: list[dict] = []
    self.metadata: dict = {}

  def record(self, event: dict) -> None:
    """Persist each attempt and outcome so interrupted runs are not silently replayed."""
    row = {"timestamp": datetime.now(timezone.utc).isoformat(), **event}
    self.events.append(row)
    with (self.output / "execution.jsonl").open("a") as stream:
      stream.write(json.dumps(row, allow_nan=False) + "\n")

  def volumes(self) -> dict:
    """Read the actual PLR trackers used by the LiquidHandler calls."""
    return {source: tube.tracker.get_used_volume() for source, tube in self.setup.tubes.items()}

  def persist(self) -> None:
    """Record current tracked volumes and lifecycle state."""
    self.metadata.update(
      {
        "completed_phases": self.completed_phases,
        "source_volumes": self.volumes(),
        "failed": self.failed,
        "attached_tip": self.attached_tip,
        "next_step_index": self.next_step,
      }
    )
    write_json(self.output / "preparation_record.json", self.metadata)

  async def start(self) -> None:
    """Ask for this physical run's approval before constructing or setting up STARBackend."""
    require(not self.started and self.lh is None, "Session already attempted; do not replay setup.")
    validate_settings(self.settings, physical=not self.settings.simulate)
    require(
      not self.output.exists(), "Use a fresh output directory; existing run records are protected."
    )
    method_reference = self.settings.method_validation_reference.strip()
    method_note = (
      f"Operator-supplied method reference: {method_reference}"
      if method_reference
      else "Exploratory run; no prior method validation supplied"
    )
    if not self.settings.simulate:
      if self.entry_phase == "plate":
        prepared = "; ".join(
          f"{source}: {self.initial_volumes[source]:g} uL at {concentration:g} ng/uL"
          for source, concentration in CONCENTRATIONS.items()
        )
        deck_note = (
          f"Verify prepared, mixed sources at these current TOTAL volumes: {prepared}. "
          "A31 must be clean matching buffer; A32 is unused. Verify the corrected deck, "
          "a fresh empty 3635 plate, fresh 300 uL tips in A1-G1 and no attached old tip. "
          "Only plate filling will run; no dilution or buffer refill. "
        )
      else:
        deck_note = (
          "Verify the corrected deck, 3000 uL TOTAL at measured 16 ng/uL in A30, "
          "4000 uL clean matching buffer each in A31/A32, empty A29-A25, replenished tip racks, "
          "no attached old tip and an empty 3635 plate. "
        )
      answer = self.confirm(
        deck_note + f"{method_note}. Type APPROVE EXP1034 RUN to authorize setup and this run: "
      )
      require(
        answer.strip() == "APPROVE EXP1034 RUN", "Not approved; no hardware backend constructed."
      )
    self.output.mkdir(parents=True, exist_ok=False)
    self.started = True
    self.metadata = {
      "experiment": "Exp1034",
      "mode": "OFFLINE_SIMULATION" if self.settings.simulate else "PHYSICAL_PREPARATION",
      "entry_phase": self.entry_phase,
      "initial_source_volumes": {
        source: volume
        for source, volume in self.initial_volumes.items()
        if self.entry_phase != "plate" or source in CONCENTRATIONS
      },
      "prepared_concentrations": dict(CONCENTRATIONS) if self.entry_phase == "plate" else None,
      "unused_sources": ["A32"] if self.entry_phase == "plate" else [],
      "settings": asdict(self.settings),
      "method_validation": {
        "status": "reference_provided" if method_reference else "not_provided",
        "reference": method_reference or None,
        "note": method_note,
      },
      "environment": environment(),
      "a30_measured_concentration": self.settings.measured_a30_concentration,
      "a30_total_volume": self.initial_volumes["A30"],
      "plate_catalog": "3635",
      "source_limits": self.limits,
      "operator": self.settings.operator,
      "started_at": datetime.now(timezone.utc).isoformat(),
      "plate_completed_at": None,
      "per_well_dispensing_timestamps": {},
      "refill_acknowledgements": [],
      "nominal_loss_model": "No programmed over-aspiration. Total source losses reserved in the floor/height checks.",
      "backend_version": "PLR legacy LiquidHandler + STARBackend or ChatterboxBackend",
      "protocol_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      "tube_definition_sha256": hashlib.sha256(
        (Path(__file__).parent / "mckesson_tubes.py").read_bytes()
      ).hexdigest(),
    }
    set_tip_tracking(True)
    set_volume_tracking(True)
    for source, volume in self.initial_volumes.items():
      self.setup.tubes[source].set_volume(volume)
    for record in manifest():
      self.setup.plate.get_item(record.well).tracker.set_volume(0)
    write_json(self.output / "resource_audit.json", self.setup.audit)
    write_json(self.output / "planned_operations.json", self.plan)
    write_json(self.output / "verification.json", self.report)
    write_csv(self.output / "well_manifest.csv", [asdict(r) for r in manifest()])
    write_csv(self.output / "dilution_ledger.csv", self.ledger)
    self.metadata["well_manifest_sha256"] = hashlib.sha256(
      (self.output / "well_manifest.csv").read_bytes()
    ).hexdigest()
    backend = (
      LiquidHandlerChatterboxBackend(num_channels=8) if self.settings.simulate else STARBackend()
    )
    self.lh = LiquidHandler(backend=backend, deck=self.setup.deck)
    self.persist()
    try:
      with (
        contextlib.redirect_stdout(io.StringIO())
        if self.settings.simulate
        else contextlib.nullcontext()
      ):
        if self.settings.simulate:
          await self.lh.setup()
        else:
          await self.lh.setup(skip_autoload=True)
          self.metadata["instrument_serial_and_firmware"] = {
            "serial": await backend.request_device_serial_number(),
            "firmware": await backend.request_firmware_version(),
          }
          targets = [*self.setup.tubes.values(), *self.setup.plate.get_all_items()]
          for rack in self.setup.racks.values():
            targets.extend(rack.get_all_items())
          for target in targets:
            require(
              backend.can_reach_position(
                self.settings.channel, target.get_location_wrt(self.setup.deck, x="c", y="c")
              ),
              f"Configured channel cannot reach {target.name}; no pipetting was started.",
            )
      self.record({"action": "setup", "status": "completed", "simulation": self.settings.simulate})
      self.persist()
    except BaseException:
      self.failed = True
      self.persist()
      raise

  async def aspirate(self, event: dict, mixing: bool) -> None:
    """Use explicit settings so no hidden liquid-class correction or air overhead is introduced."""
    require(self.lh is not None, "Start the session first.")
    source = self.setup.tubes[event["source"]]
    after = source.tracker.get_used_volume() - event["volume"]
    require(
      after >= self.settings.source_dead_volume + self.settings.source_total_loss_allowance,
      "Runtime source volume no longer covers the validated floor.",
    )
    await self.lh.aspirate(
      [source],
      vols=[event["volume"]],
      use_channels=[self.settings.channel],
      flow_rates=[
        event.get(
          "aspirate_flow_rate", self.settings.mix_flow if mixing else self.settings.aspirate_flow
        )
      ],
      liquid_height=[event["aspiration_height"]],
      surface_following_distance=[event["surface_following_distance"]],
      immersion_depth=[0.0],
      lld_mode=[STARBackend.LLDMode.OFF],
      transport_air_volume=[0.0],
      pre_wetting_volume=[0.0],
      blow_out_air_volume=[0.0],
      hamilton_liquid_classes=[None],
      disable_volume_correction=[True],
      swap_speed=[10.0],
      settling_time=[self.settings.settling_time],
      minimum_traverse_height_at_beginning_of_a_command=self.settings.traverse_height,
      min_z_endpos=self.settings.traverse_height,
    )

  async def dispense(self, event: dict, mixing: bool) -> None:
    """Use the resolved dispense heights, modes and flows, including plate height corrections."""
    require(self.lh is not None, "Start the session first.")
    destination_name = event["source"] if mixing else event["destination"]
    destination = (
      self.setup.tubes[destination_name]
      if destination_name in self.setup.tubes
      else self.setup.plate.get_item(destination_name)
    )
    await self.lh.dispense(
      [destination],
      vols=[event["volume"]],
      use_channels=[self.settings.channel],
      flow_rates=[
        event.get(
          "dispense_flow_rate", self.settings.mix_flow if mixing else self.settings.dispense_flow
        )
      ],
      liquid_height=[event["dispense_height"]],
      surface_following_distance=[0.0],
      immersion_depth=[0.0],
      lld_mode=[STARBackend.LLDMode.OFF],
      jet=[event.get("dispense_jet", not mixing and event.get("noncontact", True))],
      blow_out=[event.get("dispense_blow_out", not mixing)],
      transport_air_volume=[0.0],
      blow_out_air_volume=[0.0],
      stop_back_volume=[0.0],
      hamilton_liquid_classes=[None],
      disable_volume_correction=[True],
      swap_speed=[10.0],
      settling_time=[self.settings.settling_time],
      minimum_traverse_height_at_beginning_of_a_command=self.settings.traverse_height,
      min_z_endpos=self.settings.traverse_height,
    )

  def acknowledge_refill(self, event: dict) -> None:
    """Confirm the named buffer total before updating tracking, without issuing robot commands."""
    require(
      math.isclose(self.volumes()[event["source"]], event["before"], abs_tol=1e-6),
      "Source tracking differs at the refill checkpoint.",
    )
    if not self.settings.simulate:
      response = self.confirm(
        f"Channels have returned to traverse height. Confirm {event['source']} contains "
        f"4000 uL TOTAL of identical clean buffer (nominal addition "
        f"{event['nominal_addition']:g} uL). If already refilled, do not add more. "
        f"Type refilled or REFILLED {event['source']} to acknowledge, or cancel to pause: "
      )
      normalized = " ".join(response.casefold().split())
      if normalized not in ("refilled", f"refilled {event['source']}".casefold()):
        raise RefillNotAcknowledged(
          "Refill cancelled; no tracking update. Rerun this checkpoint when ready and type "
          f"refilled or REFILLED {event['source']}. Do not refill twice."
        )
    self.setup.tubes[event["source"]].set_volume(4000)
    self.metadata["refill_acknowledgements"].append(
      {
        "source": event["source"],
        "nominal_addition": event["nominal_addition"],
        "target_total": 4000,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "simulated": self.settings.simulate,
      }
    )

  async def execute_event(self, event: dict) -> None:
    """Dispatch only the validated plan's small set of operations."""
    require(self.lh is not None, "Start the session first.")
    action = event["action"]
    if action == "pick_up_tip":
      require(self.attached_tip is None, "Tip already attached.")
      rack = self.setup.racks[event["tip_size"]]
      await self.lh.pick_up_tips(
        [rack.get_item(event["rack_spot"])],
        use_channels=[self.settings.channel],
        minimum_traverse_height_at_beginning_of_a_command=self.settings.traverse_height,
      )
      self.attached_tip = (event["tip_size"], event["rack_spot"])
    elif action == "discard_tip":
      await self.lh.discard_tips(
        use_channels=[self.settings.channel],
        allow_nonzero_volume=False,
        minimum_traverse_height_at_beginning_of_a_command=self.settings.traverse_height,
        z_position_at_end_of_a_command=self.settings.traverse_height,
      )
      self.attached_tip = None
    elif action in ("transfer", "mix_cycle"):
      mixing = action == "mix_cycle"
      expected = event["source_after"] if mixing else event["source_before"]
      require(
        math.isclose(self.volumes()[event["source"]], expected, abs_tol=1e-6),
        "Runtime source tracking diverged from the plan; inspect before recovery.",
      )
      self.record(
        {
          "action": "aspirate",
          "step": event["step"],
          "status": "attempt",
          "volume": event["volume"],
        }
      )
      await self.aspirate(event, mixing)
      self.record({"action": "aspirate", "step": event["step"], "status": "completed"})
      self.record(
        {
          "action": "dispense",
          "step": event["step"],
          "status": "attempt",
          "volume": event["volume"],
        }
      )
      await self.dispense(event, mixing)
      self.record({"action": "dispense", "step": event["step"], "status": "completed"})
      if not mixing and event["destination"] not in self.setup.tubes:
        self.metadata["per_well_dispensing_timestamps"][event["destination"]] = datetime.now(
          timezone.utc
        ).isoformat()
    elif action == "operator_refill_checkpoint":
      self.acknowledge_refill(event)
    else:
      raise ValueError(f"Unsupported action: {action}")

  async def execute_phase(self, phase: str) -> None:
    """Run the next phase once; unacknowledged refills pause, hardware failures latch."""
    require(
      self.started and not self.failed and not self.stopped,
      "Session is not ready, failed or stopped. "
      "For a new empty plate with prepared tubes, create a new ExperimentSession with "
      "prepared_source_volumes and a fresh output directory, then await session.start(). "
      "Do not reset a completed or failed session's flags.",
    )
    require(
      phase not in self.completed_phases, "Phase already completed; never replay on a filled deck."
    )
    require(
      self.next_step < len(self.plan) and self.plan[self.next_step]["phase"] == phase,
      "Run the next planned phase; plate-only sessions skip preparation and refill.",
    )
    try:
      with (
        (self.output / "simulation_backend.log").open("a")
        if self.settings.simulate
        else contextlib.nullcontext() as log
      ):
        with (
          contextlib.redirect_stdout(log) if self.settings.simulate else contextlib.nullcontext()
        ):
          while self.next_step < len(self.plan) and self.plan[self.next_step]["phase"] == phase:
            event = self.plan[self.next_step]
            self.record({"action": "plan_step", "status": "attempt", "operation": event})
            if event["action"] == "operator_refill_checkpoint":
              # Use the current checkpoint handler with an existing live session as well.
              ExperimentSession.acknowledge_refill(self, event)
            else:
              await self.execute_event(event)
            self.next_step += 1
            self.record(
              {
                "action": "plan_step",
                "status": "completed",
                "step": event["step"],
                "source_volumes": self.volumes(),
              }
            )
            self.persist()
      if phase == "plate":
        for record in manifest():
          require(
            math.isclose(
              self.setup.plate.get_item(record.well).tracker.get_used_volume(), record.final_volume
            ),
            f"Incorrect final tracked volume in {record.well}.",
          )
        require(
          self.volumes() == self.report["source_residuals"], "Incorrect final source volumes."
        )
        self.metadata["plate_completed_at"] = datetime.now(timezone.utc).isoformat()
      if phase == "preparation":
        self.metadata["dilution_completed_at"] = datetime.now(timezone.utc).isoformat()
      self.completed_phases.append(phase)
      self.persist()
      print(f"{phase}: complete ({'SIMULATION' if self.settings.simulate else 'PHYSICAL RUN'}).")
    except RefillNotAcknowledged as error:
      self.record(
        {
          "action": "refill_paused",
          "error": str(error),
          "next_step_index": self.next_step,
        }
      )
      self.persist()
      raise
    except BaseException as error:
      self.failed = True
      self.record(
        {
          "action": "failure",
          "error": repr(error),
          "next_step_index": self.next_step,
          "recovery": "Inspect the deck and attached tip; no automatic retry, discard or movement.",
        }
      )
      self.persist()
      raise

  async def stop(self) -> None:
    """Close the backend connection; this does not perform automatic tip recovery."""
    if self.lh is not None and not self.stopped:
      await self.lh.stop()
      self.stopped = True
      self.record({"action": "stop", "status": "completed"})


async def resume_preparation(session: ExperimentSession, *, no_liquid_aspirated: bool) -> None:
  """Prepare a retry of the first A25 buffer aspiration; issue no hardware commands.

  Requires the operator's confirmation that the failed A32 aspiration took no liquid
  and the original empty tip and deck are still in place. Only this pre-mixing
  checkpoint is supported. The caller then executes the remaining preparation.
  """
  require(
    no_liquid_aspirated is True,
    "Resume requires confirmation that no liquid was aspirated and the attached tip is empty.",
  )
  require(
    session.started and session.lh is not None and not session.stopped and session.failed,
    "Keep the original failed live session; do not restart the kernel or create a new session.",
  )
  require(
    session.entry_phase == "preparation"
    and session.completed_phases == []
    and session.next_step == 9
    and len(session.plan) > 10,
    "Recovery is limited to the failed A32-to-A25 buffer aspiration at step 10.",
  )
  event = session.plan[session.next_step]
  expected = {
    "step": 10,
    "action": "transfer",
    "phase": "preparation",
    "source": "A32",
    "destination": "A25",
    "volume": 750,
    "source_before": 1000,
    "source_after": 250,
    "destination_before": 0,
    "destination_after": 750,
    "tip_size": 1000,
    "buffer_prefill": True,
  }
  require(
    all(event.get(key) == value for key, value in expected.items()),
    "The pending operation is not the supported A32-to-A25 aspiration.",
  )
  require(len(session.events) >= 3, "Missing aspiration failure history.")
  attempt, aspiration, failure = session.events[-3:]
  require(
    attempt.get("action") == "plan_step"
    and attempt.get("status") == "attempt"
    and attempt.get("operation") == event
    and aspiration.get("action") == "aspirate"
    and aspiration.get("status") == "attempt"
    and aspiration.get("step") == 10
    and failure.get("action") == "failure"
    and failure.get("error") == "ChannelizedError()"
    and failure.get("next_step_index") == 9,
    "Recovery requires an aspiration failure with no completed aspiration or dispense attempt.",
  )
  completed = [
    row["step"]
    for row in session.events
    if row.get("action") == "plan_step" and row.get("status") == "completed"
  ]
  require(completed == list(range(1, 10)), "Completed steps differ; no transfers may be replayed.")
  expected_volumes = {
    "A32": 1000,
    "A31": 1000,
    "A30": 3000,
    "A29": 1500,
    "A28": 1500,
    "A27": 1500,
    "A26": 1500,
    "A25": 0,
  }
  require(
    session.volumes() == expected_volumes, "Source trackers changed; inspect before recovery."
  )
  require(
    not session.metadata["per_well_dispensing_timestamps"]
    and all(w.tracker.get_used_volume() == 0 for w in session.setup.plate.get_all_items()),
    "Plate filling has begun; this preparation recovery does not apply.",
  )
  head = session.lh.head[session.settings.channel]
  spot = session.setup.racks[1000].get_item("A1")
  require(
    session.attached_tip == (1000, "A1")
    and head.has_tip
    and not spot.has_tip()
    and head.get_tip_origin() is spot
    and head.get_tip().tracker.get_used_volume() == 0,
    "The original empty 1000 uL A1 tip must remain attached on the same channel.",
  )
  original_settings = asdict(session.settings)
  settings = Settings(**{**original_settings, "preparation_minimum_tip_height": 4.3})
  validate_settings(settings, physical=not settings.simulate)
  resolved = geometry_for_plan(session.setup, session.plan, settings)
  plan = session.plan[: session.next_step] + resolved[session.next_step :]
  archive = session.output / "planned_operations.before_preparation_resume.json"
  require(
    not archive.exists(), "A recovery archive already exists; do not prepare the retry twice."
  )
  require(
    json.loads((session.output / "planned_operations.json").read_text()) == session.plan,
    "The saved plan differs from this live session; no recovery was applied.",
  )
  revision = {
    "action": "preparation_resume",
    "next_step_index": session.next_step,
    "operator_confirmed_no_liquid_aspirated": True,
    "original_tip_reused": [1000, "A1"],
    "source_volumes": session.volumes(),
    "previous_plan": archive.name,
    "previous_settings": original_settings,
    "settings": asdict(settings),
    "aspiration_height": plan[9]["aspiration_height"],
    "aspiration_end_height": plan[9]["aspiration_end_height"],
    "protocol_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  }
  write_json(archive, session.plan)
  write_json(session.output / "planned_operations.json", plan)
  session.record(revision)
  session.settings = settings
  session.plan = plan
  session.metadata["settings"] = asdict(settings)
  session.metadata["preparation_resume"] = revision
  try:
    session.failed = False
    session.persist()
  except BaseException:
    session.failed = True
    raise
  print("Recovery prepared at step 10: A32 → A25, aspiration 9.08 → 4.30 mm.")
  print("Same empty tip and completed transfers preserved; no robot commands issued.")
  print('Continue with await session.execute_phase("preparation").')


async def resume_refill(session: ExperimentSession) -> None:
  """Resume only the final A31 acknowledgement on the same live session, without robot motion."""
  require(
    session.started and session.lh is not None and not session.stopped,
    "Keep the existing live session; a stopped or unstarted session cannot resume here.",
  )
  require(
    session.completed_phases == ["preparation"] and session.attached_tip is None,
    "Resume requires completed preparation, no attached tip and no plate dispensing.",
  )
  require(0 <= session.next_step < len(session.plan) - 1, "No pending refill checkpoint.")
  event = session.plan[session.next_step]
  require(
    event["action"] == "operator_refill_checkpoint"
    and event["phase"] == "refill"
    and event["source"] == "A31"
    and session.plan[session.next_step + 1]["phase"] == "plate",
    "Resume is limited to the final A31 refill immediately before plate dispensing.",
  )
  require(
    all(well.tracker.get_used_volume() == 0 for well in session.setup.plate.get_all_items())
    and not session.metadata["per_well_dispensing_timestamps"],
    "Plate dispensing has already begun; do not replay the refill.",
  )
  prepared = next(row for row in session.ledger if row["stage"] == "mixed_A25")
  require(
    all(
      math.isclose(volume, prepared[source], abs_tol=1e-6)
      for source, volume in session.volumes().items()
    ),
    "Tracked sources differ from completed preparation; inspect before recovery.",
  )
  if session.failed:
    require(len(session.events) >= 2, "Missing failure record; cannot resume.")
    failure, attempt = session.events[-1], session.events[-2]
    require(
      failure.get("action") == "failure"
      and failure.get("error") == "ValueError('Refill cancelled; no tracking update.')"
      and failure.get("next_step_index") == session.next_step
      and attempt.get("action") == "plan_step"
      and attempt.get("status") == "attempt"
      and attempt.get("operation") == event,
      "Resume can clear only the refill acknowledgement error, not an instrument failure.",
    )
  session.record(
    {
      "action": "refill_resume",
      "next_step_index": session.next_step,
      "cleared_acknowledgement_failure": session.failed,
      "protocol_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
  )
  session.failed = False
  session.persist()
  await ExperimentSession.execute_phase(session, "refill")


async def update_plate_tip_reuse(session: ExperimentSession) -> None:
  """Replan unfilled wells on a live session, verifying an interrupted pickup before recovery.

  Query physical tip presence when needed; issue no motion or liquid commands. A mounted,
  empty tip from an interrupted pickup is scheduled for discard before the grouped transfers.
  A completed A25 diagnostic is accounted for separately and leaves execution paused for resume.
  """
  require(
    session.started and session.lh is not None and not session.stopped,
    "Keep the existing live session; do not rerun setup or preparation.",
  )
  require(
    session.completed_phases == ["preparation", "refill"],
    "Tip reuse can be applied only after refill and before the plate phase completes.",
  )
  require(
    "plate_resume" not in session.metadata, "A resume plan is already saved; do not replan it."
  )
  diagnostic = session.metadata.get("a25_bottom_300ul_test")
  diagnostic_withdrawal = 0.0
  if diagnostic is not None:
    require(
      diagnostic["status"] == "completed"
      and diagnostic["stage"] == "aspirate_300ul_and_retract"
      and diagnostic["source"] == "A25"
      and diagnostic["tip_choice"] == "attached_B1"
      and diagnostic["volume"] == diagnostic["tip_volume_after"] == 300
      and diagnostic["source_before"] == 3000
      and diagnostic["source_after"] == 2700,
      "The A25 diagnostic is incomplete or differs from the confirmed 300 uL test.",
    )
    diagnostic_withdrawal = 300.0
  require(0 <= session.next_step < len(session.plan), "No remaining plate steps.")
  interrupted = session.plan[session.next_step]
  require(
    interrupted["phase"] == "plate" and interrupted["action"] == "pick_up_tip",
    "Replanning requires a plate tip-pickup boundary, not an aspiration or dispense.",
  )
  if session.failed:
    require(len(session.events) >= 2, "No cancellation record to reconcile.")
    failures = [i for i, row in enumerate(session.events) if row.get("action") == "failure"]
    require(bool(failures) and failures[-1] > 0, "No cancellation record to reconcile.")
    failure_index = failures[-1]
    failure, attempt = session.events[failure_index], session.events[failure_index - 1]
    require(
      diagnostic is not None or failure_index == len(session.events) - 1,
      "Additional actions followed the cancellation; inspect before replanning.",
    )
    require(
      failure.get("action") == "failure"
      and failure.get("error") == "CancelledError()"
      and failure.get("next_step_index") == session.next_step
      and attempt.get("action") == "plan_step"
      and attempt.get("status") == "attempt"
      and attempt.get("operation") == interrupted,
      "Only a cancelled tip pickup can be reconciled here; inspect other failures separately.",
    )

  records = manifest()
  completed = set(session.metadata["per_well_dispensing_timestamps"])
  require(completed <= {r.well for r in records}, "Unknown completed wells in the run record.")
  remaining = [r for r in records if r.well not in completed]
  require(bool(remaining), "All wells are already filled.")
  for record in records:
    volume = session.setup.plate.get_item(record.well).tracker.get_used_volume()
    expected = record.final_volume if record.well in completed else 0
    require(
      math.isclose(volume, expected, abs_tol=1e-6), f"Inspect {record.well}: tracking differs."
    )
  refilled = next(row for row in session.ledger if row["stage"] == "A31_refilled_to_4000_total")
  for source, volume in session.volumes().items():
    expected = (
      refilled[source]
      - sum(r.final_volume for r in records if r.well in completed and r.source_tube == source)
      - (diagnostic_withdrawal if source == "A25" else 0)
    )
    require(
      math.isclose(volume, expected, abs_tol=1e-6),
      f"Inspect {source}: volume does not match completed wells; no plan was changed.",
    )

  head = session.lh.head[session.settings.channel]
  if not session.settings.simulate and diagnostic is None:
    presence = await session.lh.backend.request_tip_presence()
    require(
      len(presence) > session.settings.channel
      and presence[session.settings.channel] is head.has_tip,
      "Hardware tip presence differs from tracking or is unknown; inspect the cancelled pickup.",
    )
  mounted = None
  if head.has_tip:
    spot = session.setup.racks[interrupted["tip_size"]].get_item(interrupted["rack_spot"])
    require(
      session.failed
      and head.get_tip_origin() is spot
      and not spot.has_tip()
      and head.get_tip().tracker.get_used_volume() == diagnostic_withdrawal,
      "The attached tip does not match the cancelled pickup and diagnostic record.",
    )
    mounted = (interrupted["tip_size"], interrupted["rack_spot"])
  require(
    session.attached_tip in (None, mounted), "Recorded attached tip differs from the pickup state."
  )
  require(
    diagnostic is None or mounted == (300, "B1"),
    "The diagnostic B1 tip must still be attached; inspect any manual changes before resume.",
  )

  # Resource dimensions are unchanged; add the operator-confirmed insertion evidence to old sessions.
  session.setup.audit["tip_reach_checks"] = build_deck().audit["tip_reach_checks"]

  dry = DryRun(
    session.settings.a30_concentration,
    session.setup.audit,
    session.limits,
    session.settings.mix_volume,
    session.settings.mix_cycles,
    immersion=session.settings.immersion,
    rim_margin=session.settings.tip_rim_margin,
  )
  for source, volume in session.volumes().items():
    dry.trackers[source].set_volume(volume)
    dry.mass[source] = CONCENTRATIONS.get(source, 0) * volume
  dry.plate(remaining)
  residuals = {
    source: refilled[source]
    - sum(r.final_volume for r in records if r.source_tube == source)
    - (diagnostic_withdrawal if source == "A25" else 0)
    for source in session.volumes()
  }
  require(dry.volumes() == residuals, "Remaining plan changed final volumes.")
  available = {
    size: [
      f"{'ABCDEFGH'[index % 8]}{index // 8 + 1}"
      for index, spot in enumerate(rack.get_all_items())
      if spot.has_tip()
    ]
    for size, rack in session.setup.racks.items()
  }
  used = {
    size: session.setup.racks[size].num_items - len(spots) for size, spots in available.items()
  }
  pending = []
  if mounted is not None and diagnostic is None:
    pending.append(
      {"action": "discard_tip", "phase": "plate", "tip_size": mounted[0], "rack_spot": mounted[1]}
    )
  current_spot = None
  for event in dry.events:
    item = dict(event)
    if item["action"] == "pick_up_tip":
      size = item["tip_size"]
      require(bool(available[size]), f"Insufficient unused {size} uL tips.")
      current_spot = available[size].pop(0)
      item["rack_spot"] = current_spot
    elif item["action"] == "discard_tip":
      item["rack_spot"] = current_spot
    item["step"] = session.next_step + len(pending) + 1
    pending.append(item)
  if mounted is not None and diagnostic is None:
    pending[0]["step"] = session.next_step + 1
  plan = session.plan[: session.next_step] + geometry_for_plan(
    session.setup, pending, session.settings
  )
  report = {
    **session.report,
    "plate_tip_policy": "one_tip_per_source",
    "tip_counts": {
      size: used[size] + dry.tips[size]
      for size in session.setup.racks
      if used[size] + dry.tips[size]
    },
    "completed_wells_preserved": sorted(completed),
    "remaining_plate_tips": dict(dry.tips),
    "source_residuals": residuals,
    "diagnostic_withdrawal": {"A25": diagnostic_withdrawal},
  }
  revisions = session.metadata.get("plate_plan_revisions", [])
  archive = session.output / f"planned_operations.before_tip_reuse_{len(revisions) + 1}.json"
  require(not archive.exists(), "Plan archive already exists; inspect before retrying.")
  write_json(archive, session.plan)
  write_json(session.output / "planned_operations.json", plan)
  write_json(session.output / "verification.json", report)
  revision = {
    "action": "plate_plan_updated",
    "tip_policy": "one_tip_per_source",
    "previous_plan": archive.name,
    "completed_wells_preserved": sorted(completed),
    "remaining_wells": len(remaining),
    "remaining_plate_tips": dict(dry.tips),
    "cancelled_pickup_tip_attached": mounted,
    "cleared_cancelled_pickup": session.failed and diagnostic is None,
    "paused_for_resume": diagnostic is not None,
    "protocol_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  }
  session.record(revision)
  session.plan = plan
  session.report = report
  session.attached_tip = mounted
  session.failed = diagnostic is not None
  session.metadata["plate_plan_revisions"] = [*revisions, revision]
  if diagnostic is not None:
    session.metadata["plate_resume"] = {
      "status": "prepared_not_authorized",
      "source_volumes": session.volumes(),
      "completed_wells": sorted(completed),
      "next_step_index": session.next_step,
      "test_tip": {"tip_size": 300, "rack_spot": "B1", "volume": 300, "source": "A25"},
      "remaining_wells": len(remaining),
      "remaining_volume": sum(r.final_volume for r in remaining),
      "hardware_tip_verification": "deferred_until_resume",
    }
    write_json(session.output / "resource_audit.plate_resume.json", session.setup.audit)
    write_csv(
      session.output / "remaining_plate_transfers.csv",
      [e for e in plan[session.next_step :] if e["action"] == "transfer"],
    )
    write_json(
      session.output / "plate_resume_checkpoint.json",
      {
        "resume": session.metadata["plate_resume"],
        "settings": asdict(session.settings),
        "source_residuals": residuals,
        "remaining_plate_tips": dict(dry.tips),
        "remaining_operations": plan[session.next_step :],
        "available_tip_spots": {
          str(size): [
            f"{'ABCDEFGH'[i % 8]}{i // 8 + 1}"
            for i, spot in enumerate(rack.get_all_items())
            if spot.has_tip()
          ]
          for size, rack in session.setup.racks.items()
        },
      },
    )
  session.persist()
  print(
    f"Plate plan updated: {len(completed)} filled wells preserved; {len(remaining)} wells remain. "
    + (
      "Paused. When ready, run await notebook_protocol.resume_plate_after_test(session)."
      if diagnostic is not None
      else "Continue with await session.execute_phase('plate')."
    )
  )


async def resume_plate_after_test(session: ExperimentSession) -> None:
  """After fresh approval, discard the diagnostic aliquot and execute the saved remaining plate."""
  require(
    session.started and session.lh is not None and not session.stopped,
    "Use the existing live session. If the kernel was restarted, recover the saved checkpoint first.",
  )
  resume = session.metadata.get("plate_resume")
  require(
    resume is not None and resume["status"] == "prepared_not_authorized",
    "No unstarted resume plan; do not repeat a partially executed resume.",
  )
  require(session.next_step == resume["next_step_index"], "The plan advanced after it was saved.")
  require(session.volumes() == resume["source_volumes"], "Sources changed after the pause.")
  require(
    sorted(session.metadata["per_well_dispensing_timestamps"]) == resume["completed_wells"],
    "The plate changed after the pause.",
  )
  for record in manifest():
    expected = record.final_volume if record.well in resume["completed_wells"] else 0
    require(
      session.setup.plate.get_item(record.well).tracker.get_used_volume() == expected,
      f"Inspect {record.well}: volume changed after the pause.",
    )
  head = session.lh.head[session.settings.channel]
  require(
    head.has_tip
    and head.get_tip_origin() is session.setup.racks[300].get_item("B1")
    and head.get_tip().tracker.get_used_volume() == 300
    and session.attached_tip == (300, "B1"),
    "The diagnostic tip changed during the pause; inspect before resuming.",
  )
  if not session.settings.simulate:
    answer = session.confirm(
      "Check the returned samples, plate and racks against the saved deck: A9 already has "
      "200 uL buffer; B1 is attached with the diagnostic aliquot. Resume will discard B1 "
      "and its liquid, then fill 95 remaining wells using seven fresh 300 uL tips. "
      "Type yes to approve this plate run, or anything else to leave it paused: "
    )
    require(answer.strip().casefold() in {"yes", "resume"}, "Still paused; no robot command sent.")
    presence = await session.lh.backend.request_tip_presence()
    require(
      presence == [i == session.settings.channel for i in range(session.lh.backend.num_channels)],
      "Actual tip occupancy differs from the saved state; no motion was started.",
    )
  resume["status"] = "started"
  session.record({"action": "plate_resume_approved", "simulation": session.settings.simulate})
  session.persist()
  try:
    if not session.settings.simulate:
      backend = session.lh.backend
      # Restore normal Z drive settings after the deliberately slow diagnostic approach.
      current_z = await backend.request_probe_z_position(session.settings.channel)
      await backend.move_channel_stop_disk_z(
        channel_idx=session.settings.channel,
        z=current_z,
        speed=125.0,
        acceleration=800.0,
        current_limit=3,
      )
    session.record({"action": "discard_diagnostic_tip", "status": "attempt", **resume["test_tip"]})
    await session.lh.discard_tips(
      use_channels=[session.settings.channel],
      allow_nonzero_volume=True,
      minimum_traverse_height_at_beginning_of_a_command=session.settings.traverse_height,
      z_position_at_end_of_a_command=session.settings.traverse_height,
    )
    session.attached_tip = None
    session.record(
      {"action": "discard_diagnostic_tip", "status": "completed", **resume["test_tip"]}
    )
    session.failed = False
    session.persist()
    await ExperimentSession.execute_phase(session, "plate")
    resume["status"] = "completed"
    session.persist()
  except BaseException as error:
    session.failed = True
    resume["status"] = "failed"
    session.record({"action": "plate_resume_failed", "error": repr(error)})
    session.persist()
    raise
