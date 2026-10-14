"""Offline checks for Orange G acquisition, well mapping and approval boundaries."""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import AsyncMock, Mock, patch

from pylabrobot.molecular_devices.filtermax import (
  AbsorbanceResult,
  FilterMaxF5,
  FilterMaxStatus,
  FilterSlideCatalog,
  InstalledFilterSlides,
  InstrumentInfo,
)

import read_orange_g as script


def measurement() -> AbsorbanceResult:
  """Give each row and column a distinct value to expose transposed exports."""
  return AbsorbanceResult(
    data=[[row + column / 100 for column in range(12)] for row in range(8)],
    wavelengths=(450,),
    reference_subtracted=False,
    temperature=23.5,
    timestamp=1791375773,
  )


def fake_reader() -> Mock:
  """Mock only public reader calls, without constructing a serial transport."""
  reader = Mock(spec=FilterMaxF5)
  reader.setup = AsyncMock()
  reader.stop = AsyncMock()
  reader.get_instrument_info = AsyncMock(
    return_value=InstrumentInfo(
      model="Anthos Fluoro",
      firmware="test",
      device_number=141,
      serial_number="1191",
      excitation_slide_id=2,
      emission_slide_id=1,
      device_code=57855,
      plate_control=True,
      pic_firmware="test",
      cpld_version=6,
      raw="synthetic test identity",
    )
  )
  reader.get_status = AsyncMock(return_value=FilterMaxStatus(True, "READY", ()))
  reader.get_filter_slides = AsyncMock(return_value=InstalledFilterSlides(2, 1))
  reader.move_plate_tray_out = AsyncMock()
  reader.read_absorbance = AsyncMock(return_value=[measurement()])
  return reader


class OfflineTests(unittest.TestCase):
  """Check the default offline path, plate identities and optical data processing."""

  def test_default_does_not_construct_hardware_or_generate_data(self) -> None:
    """No flags must be safe even when no adapter is present."""
    with tempfile.TemporaryDirectory() as folder:
      output = Path(folder) / "unused"
      with (
        patch.object(script, "FilterMaxF5", side_effect=AssertionError("hardware forbidden")),
        contextlib.redirect_stdout(io.StringIO()),
      ):
        self.assertEqual(script.main(["--output", str(output)]), 0)
      self.assertFalse(output.exists())

  def test_fluorescence_only_filter_is_rejected(self) -> None:
    """485 nm exists in the saved catalog but cannot be selected for absorbance."""
    with (
      patch.object(script, "FilterMaxF5", side_effect=AssertionError("hardware forbidden")),
      contextlib.redirect_stderr(io.StringIO()),
    ):
      self.assertEqual(script.main(["--wavelength", "485"]), 1)

  def test_mapping_and_matched_volume_blanks(self) -> None:
    """All 96 coordinates, sample IDs and four independent blank means survive export."""
    records = script.load_manifest()
    self.assertEqual(len(records), 96)
    self.assertEqual(sum(record["source_tube"] == "A31" for record in records), 24)
    rows, blanks = script.result_rows(measurement(), records, 450, True)
    by_well = {row["well"]: row for row in rows}
    self.assertEqual(by_well["A1"]["raw_od"], 0.0)
    self.assertEqual(by_well["H12"]["raw_od"], 7.11)
    self.assertEqual(by_well["H2"]["sample_id"], "A30_200uL_R1")
    self.assertNotIn("final_concentration", by_well["H2"])
    for blank in blanks:
      expected = sum(by_well[well]["raw_od"] for well in blank["blank_wells"].split(";")) / 6
      self.assertAlmostEqual(blank["blank_mean_od"], expected)
      for row in rows:
        if row["planned_volume_uL"] == blank["planned_volume_uL"]:
          self.assertAlmostEqual(row["blank_corrected_od"], row["raw_od"] - expected)
    self.assertLess(by_well["A1"]["blank_corrected_od"], 0)

  def test_bad_blank_invalidates_only_its_volume_group(self) -> None:
    """Do not quietly calculate a six-blank mean from five finite values."""
    result = measurement()
    result.data[0][8] = float("nan")  # A9 is a 200 uL A31 blank.
    rows, blanks = script.result_rows(result, script.load_manifest(), 450, True)
    self.assertEqual(
      next(b for b in blanks if b["planned_volume_uL"] == 200)["finite_blank_count"], 5
    )
    for row in rows:
      if row["planned_volume_uL"] == 200:
        self.assertIsNone(row["blank_corrected_od"])
        self.assertEqual(row["correction_status"], "incomplete_blank_group")
      else:
        self.assertIsNotNone(row["blank_corrected_od"])
    self.assertEqual(script.json_safe(result.data)[0][8], "nan")

  def test_wrong_shape_or_optical_mode_is_rejected(self) -> None:
    """Malformed data must not silently acquire correct-looking well labels."""
    for result in (
      replace(measurement(), data=[[0.0] * 8 for _ in range(12)]),
      replace(measurement(), wavelengths=(260,)),
      replace(measurement(), reference_subtracted=True),
    ):
      with self.subTest(result=result), self.assertRaises(ValueError):
        script.result_rows(result, script.load_manifest(), 450, False)


class AcquisitionTests(unittest.IsolatedAsyncioTestCase):
  """Exercise the physical workflow against mock reader methods only."""

  def setUp(self) -> None:
    """Create an isolated output destination and offline inputs."""
    self.temp = tempfile.TemporaryDirectory()
    self.addCleanup(self.temp.cleanup)
    self.args = argparse.Namespace(
      output=Path(self.temp.name) / "run",
      port="FAKE",
      expected_serial="1191",
      wavelength=450,
      blank_correct=False,
      geometry=None,
    )
    self.catalog = FilterSlideCatalog.from_softmax_xml(script.DEFAULT_CATALOG)
    self.records = script.load_manifest()
    self.geometry = script.load_geometry(None)
    self.reader = fake_reader()

  async def execute(self, answers: list[str]) -> None:
    """Run with supplied approvals and a mocked constructor."""
    with (
      patch.object(script, "FilterMaxF5", return_value=self.reader),
      contextlib.redirect_stdout(io.StringIO()),
    ):
      await script.acquire(
        self.args, self.catalog, self.geometry, self.records, {}, Mock(side_effect=answers)
      )

  async def test_rejected_approval_prevents_construction(self) -> None:
    """Rejection creates neither a reader nor a run directory."""
    with (
      patch.object(script, "FilterMaxF5") as factory,
      contextlib.redirect_stdout(io.StringIO()),
      self.assertRaisesRegex(ValueError, "not approved"),
    ):
      await script.acquire(self.args, self.catalog, self.geometry, self.records, {}, lambda _: "NO")
    factory.assert_not_called()
    self.assertFalse(self.args.output.exists())

  async def test_existing_run_is_preserved(self) -> None:
    """Refuse replay into an existing result directory before approval or I/O."""
    self.args.output.mkdir()
    with (
      patch.object(script, "FilterMaxF5") as factory,
      self.assertRaisesRegex(ValueError, "already exists"),
    ):
      await script.acquire(
        self.args, self.catalog, self.geometry, self.records, {}, lambda _: script.APPROVAL
      )
    factory.assert_not_called()

  async def test_identity_and_filter_failures_prevent_motion(self) -> None:
    """A different reader or fluorescence slide is rejected before tray movement."""
    for failure in ("serial", "slide"):
      with self.subTest(failure=failure):
        self.args.output = Path(self.temp.name) / failure
        self.reader = fake_reader()
        if failure == "serial":
          self.args.expected_serial = "different"
        else:
          self.args.expected_serial = "1191"
          self.reader.get_filter_slides.return_value = InstalledFilterSlides(1, 1)
        with self.assertRaises(Exception):
          await self.execute([script.APPROVAL])
        self.reader.move_plate_tray_out.assert_not_awaited()
        self.reader.read_absorbance.assert_not_awaited()
        self.reader.stop.assert_awaited_once()

  async def test_loading_cancel_does_not_measure(self) -> None:
    """Cancelling at the open tray leaves measurement unrequested."""
    with self.assertRaisesRegex(ValueError, "Loading cancelled"):
      await self.execute([script.APPROVAL, "CANCEL"])
    self.reader.read_absorbance.assert_not_awaited()
    self.reader.stop.assert_awaited_once()

  async def test_one_endpoint_uses_current_api_and_saves_96_wells(self) -> None:
    """No fluorescence, reference, kinetic, shaking or obsolete keyword is sent."""
    await self.execute([script.APPROVAL, "LOADED"])
    self.reader.read_absorbance.assert_awaited_once_with(
      self.geometry,
      wavelength=450,
      wells=script.WELLS,
      slide_id=2,
      filter_position=4,
    )
    self.reader.stop.assert_awaited_once()
    with (self.args.output / "absorbance.csv").open() as stream:
      rows = list(csv.DictReader(stream))
    self.assertEqual(len(rows), 96)
    self.assertEqual(rows[-1]["well"], "H12")
    self.assertFalse((self.args.output / "blank_summary.csv").exists())
    metadata = json.loads((self.args.output / "reader_metadata.json").read_text())
    self.assertEqual(metadata["status"], "complete")
    self.assertEqual(metadata["mode"], "PHYSICAL_ORANGE_G_READ")

  async def test_raw_payload_is_saved_before_invalid_shape_is_rejected(self) -> None:
    """Keep the original driver result available when export validation fails."""
    self.reader.read_absorbance.return_value = [replace(measurement(), data=[[1.0]])]
    with self.assertRaisesRegex(ValueError, "8 by 12"):
      await self.execute([script.APPROVAL, "LOADED"])
    self.assertTrue((self.args.output / "raw_backend_results.json").exists())
    metadata = json.loads((self.args.output / "reader_metadata.json").read_text())
    self.assertEqual(metadata["status"], "failed_or_cancelled")
    self.reader.stop.assert_awaited_once()

  async def test_read_error_survives_shutdown_error(self) -> None:
    """A cleanup failure cannot hide the acquisition failure or trigger a retry."""
    self.reader.read_absorbance.side_effect = RuntimeError("read failed")
    self.reader.stop.side_effect = RuntimeError("stop failed")
    with (
      contextlib.redirect_stderr(io.StringIO()),
      self.assertRaisesRegex(RuntimeError, "read failed"),
    ):
      await self.execute([script.APPROVAL, "LOADED"])
    self.reader.read_absorbance.assert_awaited_once()
    metadata = json.loads((self.args.output / "reader_metadata.json").read_text())
    self.assertIn("read failed", metadata["error"])
    self.assertIn("stop failed", metadata["shutdown_error"])


if __name__ == "__main__":
  unittest.main()
