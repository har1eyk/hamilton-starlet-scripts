"""Offline regression checks. No test opens USB/serial or drives hardware."""

import asyncio
import csv
import json
import math
import statistics
import tempfile
import unittest
from pathlib import Path
from unittest.mock import AsyncMock, patch

from pylabrobot.molecular_devices.filtermax import (
  AbsorbanceResult,
  FilterSelection,
  FilterMaxStatus,
  InstalledFilterSlides,
  InstrumentInfo,
  PlateGeometry,
)

import read_filtermax
from analyze_a260 import analyze, load_readings
from exp1034 import INITIAL, DryRun, dispensing_order, inspect_resources, manifest


def limits(dead: float = 0, loss: float = 0) -> dict:
  """Return explicit synthetic limits, not measurements of the physical tubes."""
  return {s: {"dead_volume": dead, "total_loss_allowance": loss} for s in INITIAL}


def synthetic_readings() -> dict:
  """Generate known linear responses with nonzero intercept and variable blanks."""
  readings = {}
  for record in manifest():
    value = (
      0.1 + (record.replicate - 3.5) * 0.001
      if record.final_concentration == 0
      else 0.3 + record.final_concentration * record.final_volume / 10000
    )
    readings[record.well] = {"a260": value, "valid": True, "exclusion_reason": ""}
  return readings


def synthetic_ranges() -> dict:
  """Return a reviewed test-only full range."""
  return {
    str(v): {
      "concentrations": [0.5, 1, 2, 4, 8, 16],
      "reviewed": True,
      "reason": "synthetic linear fixture",
    }
    for v in (200, 100, 75, 50)
  }


class ExperimentTests(unittest.TestCase):
  """Check experiment invariants and failure paths."""

  @classmethod
  def setUpClass(cls) -> None:
    """Inspect installed resources without setting up a device."""
    cls.resources = inspect_resources()

  def test_nominal_budget_and_chain(self) -> None:
    """Verify complete plate, source balances, five serial tips and one refill."""
    run = DryRun(16, self.resources, limits())
    report = run.run()
    self.assertEqual(report["whole_plate_volume"], 10200)
    self.assertEqual(report["tip_counts"], {1000: 5, 300: 7})
    self.assertEqual(
      report["source_residuals"],
      {
        "A32": 250,
        "A31": 1450,
        "A30": 225,
        "A29": 225,
        "A28": 225,
        "A27": 225,
        "A26": 225,
        "A25": 1725,
      },
    )
    self.assertEqual(len(report["refills"]), 1)
    self.assertEqual(report["refills"][0]["nominal_addition"], 3750)
    self.assertEqual(len([e for e in run.events if e["action"] == "mix_cycle"]), 48)

  def test_manifest_labels_and_order(self) -> None:
    """Check dispersed positions and label preservation in execution order."""
    records = manifest()
    mapped = {r.well: r for r in records}
    self.assertEqual((mapped["H2"].source_tube, mapped["H2"].replicate), ("A30", 1))
    self.assertEqual((mapped["A1"].final_volume, mapped["A1"].replicate), (200, 2))
    self.assertEqual((mapped["C10"].final_concentration, mapped["C10"].replicate), (0.5, 3))
    self.assertEqual(set(records), set(dispensing_order(records)))

  def test_one_plate_tip_per_concentration(self) -> None:
    """A plate tip serves one entire source group and never crosses concentrations."""
    run = DryRun(16, self.resources, limits(200, 25))
    run.run()
    sources = []
    group = []
    for event in (e for e in run.events if e["phase"] == "plate"):
      if event["action"] == "pick_up_tip":
        self.assertEqual(group, [])
        self.assertEqual(event["tip_size"], 300)
      elif event["action"] == "transfer":
        group.append(event)
        self.assertTrue(event["noncontact"])
      elif event["action"] == "discard_tip":
        self.assertEqual(len({e["source"] for e in group}), 1)
        source = group[0]["source"]
        sources.append(source)
        self.assertEqual(len(group), 24 if source == "A31" else 12)
        self.assertEqual({e["volume"] for e in group}, {50, 75, 100, 200})
        group = []
    self.assertEqual(sources, ["A31", "A30", "A29", "A28", "A27", "A26", "A25"])
    self.assertEqual(group, [])

  def test_wrong_or_nonfinite_concentration_rejected(self) -> None:
    """Do not silently relabel a dilution from the wrong stock."""
    for concentration in (15, 32, 0, math.nan, math.inf):
      with self.subTest(concentration=concentration), self.assertRaises(ValueError):
        DryRun(concentration, self.resources, limits())

  def test_dna_dead_volume_and_loss_allowance(self) -> None:
    """A reserve exceeding the 225 uL allowance must fail before physical use."""
    for dead, loss in ((226, 0), (200, 26)):
      with (
        self.subTest(dead=dead, loss=loss),
        self.assertRaisesRegex(ValueError, "below dead\\+loss"),
      ):
        DryRun(16, self.resources, limits(dead, loss)).run()
    self.assertEqual(DryRun(16, self.resources, limits(200, 25)).run()["wells"], 96)

  def test_earlier_buffer_refills(self) -> None:
    """Higher buffer residual limits trigger top-ups while the tip is clean."""
    reserve = limits(200)
    reserve["A32"]["dead_volume"] = reserve["A31"]["dead_volume"] = 300
    result = DryRun(16, self.resources, reserve).run()
    self.assertEqual([e["nominal_addition"] for e in result["refills"]], [3000, 3000, 750])
    self.assertEqual(result["source_residuals"]["A32"], 3250)
    self.assertEqual(result["source_residuals"]["A31"], 1450)

  def test_unacknowledged_refill_does_not_change_tracking(self) -> None:
    """Cancellation must never set the source to 4000 uL."""
    run = DryRun(16, self.resources, limits())
    run.trackers["A31"].set_volume(250)
    with self.assertRaisesRegex(ValueError, "not acknowledged"):
      run.refill("A31", acknowledged=False)
    self.assertEqual(run.volumes()["A31"], 250)

  def test_capacity_and_air_reserve(self) -> None:
    """Mixing and split aspirations include air in the capacity limit."""
    for mixing, air in ((1001, 0), (750, 251)):
      with self.subTest(mixing=mixing, air=air), self.assertRaisesRegex(ValueError, "capacity"):
        DryRun(16, self.resources, limits(), mix_volume=mixing, air_gap=air).run()

  def test_destination_contact_prevents_source_return(self) -> None:
    """A mixing tip may not re-enter the higher-concentration source."""
    run = DryRun(16, self.resources, limits())
    run.pickup(1000)
    run.transfer("A32", "A29", 750)
    run.transfer("A30", "A29", 750)
    run.mix("A29")
    with self.assertRaisesRegex(ValueError, "Destination-contacted"):
      run.transfer("A30", "A29", 750)

  def test_contact_prefill_rejects_dna_before_changing_tracking(self) -> None:
    """Neither a DNA-containing destination nor a DNA-exposed tip may contact clean buffer."""
    for contaminated in ("destination", "tip"):
      with self.subTest(contaminated=contaminated):
        run = DryRun(16, self.resources, limits())
        run.pickup(1000)
        if contaminated == "destination":
          run.trackers["A29"].set_volume(100)
          run.mass["A29"] = 100
        else:
          run.mix("A30")
        volumes, mass, events = run.volumes(), dict(run.mass), list(run.events)
        with self.assertRaisesRegex(ValueError, "Contact buffer prefills require"):
          run.transfer("A32", "A29", 750, contact=True)
        self.assertEqual(run.volumes(), volumes)
        self.assertEqual(run.mass, mass)
        self.assertEqual(run.events, events)

  def test_corrected_geometry_is_used(self) -> None:
    """Honor the user's carrier, track and plate corrections."""
    self.assertEqual(self.resources["deck_tracks"], 30)
    self.assertEqual(self.resources["tip_carrier_track"], 25)
    self.assertEqual(self.resources["plate_carrier_model"], "MFX_CAR_L5_base")
    self.assertTrue(self.resources["plate_reference"]["physical_3635_correspondence_verified"])
    self.assertEqual(self.resources["tube_coordinates"]["A32"]["site"], 0)
    self.assertEqual(self.resources["tube_coordinates"]["A25"]["site"], 7)

  def test_correct_blank_sd_intercept_and_mass_lod(self) -> None:
    """Verify sample SD (not SEM), free intercept and concentration-to-mass units."""
    results, corrected = analyze(synthetic_readings(), synthetic_ranges())
    self.assertEqual(len(corrected), 72)
    for row in results:
      slope = row["volume"] / 10000
      sd = statistics.stdev(row["blank_values"])
      self.assertAlmostEqual(row["slope_A260_per_ng_uL"], slope)
      self.assertAlmostEqual(row["intercept_A260"], 0.2)
      self.assertAlmostEqual(row["preliminary_concentration_lod_ng_uL"], 3.3 * sd / slope)
      self.assertAlmostEqual(row["preliminary_mass_lod_ng"], 3.3 * sd / slope * row["volume"])

  def test_invalid_blank_and_unreviewed_range_rejected(self) -> None:
    """Do not reduce blank count or choose a calibration range silently."""
    readings = synthetic_readings()
    readings["A9"]["valid"] = False
    with self.assertRaisesRegex(ValueError, "six valid blanks"):
      analyze(readings, synthetic_ranges())
    ranges = synthetic_ranges()
    ranges["200"]["reviewed"] = False
    with self.assertRaisesRegex(ValueError, "Inspect the measured range"):
      analyze(synthetic_readings(), ranges)

  def test_invalid_sample_is_retained_but_excluded(self) -> None:
    """Reject a flagged sample from regression while preserving its raw value."""
    readings = synthetic_readings()
    readings["H2"].update(a260=999, valid=False, exclusion_reason="synthetic saturation")
    results, corrected = analyze(readings, synthetic_ranges())
    self.assertEqual(results[0]["n"], 17)
    sample = next(r for r in corrected if r["well"] == "H2")
    self.assertEqual(sample["raw_a260"], 999)
    self.assertFalse(sample["used_in_fit"])

  def test_duplicate_raw_well_rejected(self) -> None:
    """Protect against incomplete/duplicated manual reader exports."""
    with tempfile.TemporaryDirectory() as directory:
      path = Path(directory) / "raw.csv"
      with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["well"])
        writer.writeheader()
        writer.writerows([{"well": "A1"}] * 96)
      with self.assertRaisesRegex(ValueError, "exactly once"):
        load_readings(path)

  def test_reader_requires_verified_profile(self) -> None:
    """A draft profile cannot initiate an instrument session."""
    with tempfile.TemporaryDirectory() as directory:
      path = Path(directory) / "profile.json"
      path.write_text(json.dumps({"physical_plate_catalog": "3635", "geometry_verified": False}))
      with self.assertRaisesRegex(ValueError, "verified 3635"):
        read_filtermax.load_profile(path)

  def test_reader_adapter_uses_all_wells_and_saves_raw_data(self) -> None:
    """Exercise acquisition with mocked I/O and check row-major well mapping."""
    info = InstrumentInfo("F5", "fixture", 1, "fixture", 1, 2, 57855, True, "fixture", 1, "raw")
    result = AbsorbanceResult(
      [[r * 12 + c for c in range(12)] for r in range(8)], (260,), False, 22, 1700000000
    )
    reader = AsyncMock()
    reader.get_instrument_info.return_value = info
    reader.get_status.return_value = FilterMaxStatus(True, "READY", ())
    reader.get_filter_slides.return_value = InstalledFilterSlides(1, 2)
    reader.read_absorbance.return_value = [result]
    with tempfile.TemporaryDirectory() as directory:
      root = Path(directory)
      prep = root / "prep.json"
      prep.write_text(
        json.dumps(
          {
            "experiment": "Exp1034",
            "mode": "PHYSICAL_PREPARATION",
            "a30_measured_concentration": 16,
            "a30_total_volume": 3000,
            "plate_completed_at": "2023-11-14T00:00:00+00:00",
          }
        )
      )
      with (
        patch.object(read_filtermax, "FilterMaxF5", return_value=reader),
        patch("builtins.input", return_value="LOADED"),
        patch.object(
          read_filtermax.FilterSlideCatalog,
          "resolve",
          return_value=FilterSelection("excitation", 1, 1, 260, 10),
        ),
      ):
        catalog = read_filtermax.FilterSlideCatalog(())
        asyncio.run(
          read_filtermax.acquire(
            {"expected_serial": "fixture"},
            PlateGeometry.costar_96_clear_landscape(),
            catalog,
            "MOCK",
            root / "output",
            prep,
          )
        )
      self.assertEqual(len(reader.read_absorbance.call_args.kwargs["wells"]), 96)
      raw = json.loads((root / "output/raw_backend_result.json").read_text())
      self.assertEqual(raw["data"][7][11], 95)
      with (root / "output/raw_a260.csv").open() as stream:
        rows = {r["well"]: r for r in csv.DictReader(stream)}
      self.assertEqual(rows["A12"]["a260"], "11")
      self.assertEqual(rows["H1"]["a260"], "84")
      reader.stop.assert_awaited_once()


if __name__ == "__main__":
  unittest.main()
