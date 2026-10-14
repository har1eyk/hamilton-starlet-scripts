"""Offline tests for DNA A260 acquisition, exports and reader failure handling."""

from __future__ import annotations

import argparse
import asyncio
import contextlib
import csv
import io
import json
import math
import tempfile
import unittest
from dataclasses import asdict, replace
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

import analyze_a260
import read_dna_a260 as script
import read_orange_g


def measurement() -> AbsorbanceResult:
  """Give every well a distinct value, exposing transposition and blank-group mistakes."""
  return AbsorbanceResult(
    data=[[row + column / 100 for column in range(12)] for row in range(8)],
    wavelengths=(260,),
    reference_subtracted=False,
    temperature=23.5,
    timestamp=1791375773,
  )


def fake_reader() -> Mock:
  """Mock public reader methods without constructing a serial transport."""
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
  """Check the no-hardware entrypoint and numerical data handling."""

  def test_default_and_dry_run_never_construct_hardware_or_write_outputs(self) -> None:
    """Preview works without the adapter and resolves the exact A260 absorbance filter."""
    for flags in ([], ["--dry-run"]):
      with self.subTest(flags=flags), tempfile.TemporaryDirectory() as folder:
        output = Path(folder) / "unused"
        captured = io.StringIO()
        with (
          patch.object(script, "FilterMaxF5", side_effect=AssertionError("hardware forbidden")),
          contextlib.redirect_stdout(captured),
        ):
          self.assertEqual(
            script.main([*flags, "--output", str(output), "--port", "/dev/nonexistent_f5"]), 0
          )
        self.assertFalse(output.exists())
        plan = json.JSONDecoder().raw_decode(captured.getvalue())[0]
        self.assertEqual(plan["wavelength"], 260)
        self.assertEqual(plan["well_count"], 96)
        self.assertEqual(
          plan["catalog_filter"],
          {
            "kind": "excitation",
            "slide_id": 2,
            "position": 1,
            "wavelength": 260,
            "bandwidth": 12,
          },
        )
        self.assertFalse(plan["pathlength_correction_applied"])
        self.assertFalse(plan["reference_subtracted"])
        self.assertFalse(plan["raw_csv_blank_subtracted"])

  def test_geometry_and_cli_overrides(self) -> None:
    """Use the working portrait geometry, allowing explicit connection and profile overrides."""
    expected = {
      "rows": 8,
      "columns": 12,
      "length": 127.7,
      "width": 85.7,
      "height": 14.27,
      "well_depth": 10.69,
      "left_column_offset": 14.05,
      "top_row_offset": 11.18,
      "column_spacing": 9.02,
      "row_spacing": 9.0,
      "well_size_x": 6.4,
      "well_size_y": 6.4,
      "absorbance_z": 10.27,
      "name": "96 Well Costar clear [Portrait]",
      "orientation": "portrait",
    }
    self.assertEqual(asdict(script.load_geometry(None)), expected)
    self.assertEqual(script.load_geometry(None), read_orange_g.load_geometry(None))
    self.assertIn(
      "usb-Prolific_Technology_Inc._USB-Serial_Controller_BFADb2A8J14", script.DEFAULT_PORT
    )
    with tempfile.TemporaryDirectory() as folder:
      geometry_path = Path(folder) / "geometry.json"
      geometry_path.write_text(json.dumps({**expected, "orientation": "landscape"}))
      captured = io.StringIO()
      with (
        patch.object(script, "FilterMaxF5", side_effect=AssertionError("hardware forbidden")),
        contextlib.redirect_stdout(captured),
      ):
        self.assertEqual(
          script.main(
            [
              "--port",
              "FAKE",
              "--expected-serial",
              "another-reader",
              "--geometry",
              str(geometry_path),
              "--catalog",
              str(script.DEFAULT_CATALOG),
            ]
          ),
          0,
        )
      plan = json.JSONDecoder().raw_decode(captured.getvalue())[0]
      self.assertEqual(plan["port"], "FAKE")
      self.assertEqual(plan["expected_serial"], "another-reader")
      self.assertEqual(plan["plate_geometry"]["orientation"], "landscape")
      self.assertEqual(plan["geometry_source"], str(geometry_path))
      for change in ({"absorbance_z": 100}, {"rows": 12}, {"column_spacing": 99}):
        geometry_path.write_text(json.dumps({**expected, **change}))
        with self.subTest(change=change), self.assertRaises(ValueError):
          script.load_geometry(geometry_path)

  def test_missing_a260_filter_is_rejected_offline(self) -> None:
    """Do not substitute a visible filter when 260 nm is absent from the catalog."""
    with tempfile.TemporaryDirectory() as folder:
      catalog = Path(folder) / "without-a260.xml"
      catalog.write_text(script.DEFAULT_CATALOG.read_text().replace("<WL>260</WL>", "<WL>261</WL>"))
      with (
        patch.object(script, "FilterMaxF5", side_effect=AssertionError("hardware forbidden")),
        contextlib.redirect_stderr(io.StringIO()),
      ):
        self.assertEqual(script.main(["--catalog", str(catalog)]), 1)

  def test_manifest_and_volume_matched_mean_and_sample_sd(self) -> None:
    """Preserve every identity and compute independent means and n-minus-one blank SDs."""
    records = script.load_manifest()
    original = read_orange_g.load_manifest()
    self.assertEqual(len(records), 96)
    for record, prior in zip(records, original):
      for key in ("well", "sample_id", "source_tube", "planned_volume_uL", "replicate"):
        self.assertEqual(record[key], prior[key])
      self.assertEqual(
        record["planned_dna_concentration_ng_uL"], script.CONCENTRATIONS[record["source_tube"]]
      )
    self.assertEqual(sum(r["role"] == "dna_sample" for r in records), 72)
    raw, corrected, blanks = script.result_rows(measurement(), records)
    by_well = {row["well"]: row for row in raw}
    self.assertEqual(by_well["A1"]["a260"], 0)
    self.assertEqual(by_well["H12"]["a260"], 7.11)
    self.assertEqual(by_well["H2"]["sample_id"], "A30_200uL_R1")
    self.assertTrue(all(row["valid"] == row["exclusion_reason"] == "" for row in raw))
    self.assertTrue(all("blank_corrected_a260" not in row for row in raw))
    self.assertEqual(len(blanks), 4)
    for blank in blanks:
      values = [by_well[well]["a260"] for well in blank["blank_wells"].split(";")]
      self.assertEqual(len(values), 6)
      expected = sum(values) / 6
      self.assertAlmostEqual(blank["blank_mean_a260"], expected)
      self.assertAlmostEqual(
        blank["blank_sample_sd_a260"], math.sqrt(sum((v - expected) ** 2 for v in values) / 5)
      )
      for row in corrected:
        if row["planned_volume_uL"] == blank["planned_volume_uL"]:
          self.assertAlmostEqual(
            row["blank_corrected_a260"], by_well[row["well"]]["a260"] - expected
          )
    self.assertLess(corrected[0]["blank_corrected_a260"], 0)

  def test_missing_sample_and_blank_preserve_other_volumes(self) -> None:
    """One nonfinite blank withholds only its volume's correction; raw data remains intact."""
    for value in (None, float("nan"), float("inf"), float("-inf")):
      with self.subTest(value=value):
        result = measurement()
        result.data[0][8] = value  # A9: 200 uL buffer blank.
        result.data[0][0] = value  # A1: 50 uL DNA.
        raw, corrected, blanks = script.result_rows(result, script.load_manifest())
        blank = next(b for b in blanks if b["planned_volume_uL"] == 200)
        self.assertEqual(blank["finite_blank_count"], 5)
        self.assertIsNone(blank["blank_mean_a260"])
        self.assertIsNone(blank["blank_sample_sd_a260"])
        self.assertIsNone(raw[0]["a260"])
        self.assertEqual(raw[0]["data_status"], "missing_or_nonfinite")
        for row in corrected:
          if row["planned_volume_uL"] == 200:
            self.assertIsNone(row["blank_corrected_a260"])
            self.assertEqual(row["correction_status"], "incomplete_blank_group")
          elif row["well"] == "A1":
            self.assertIsNone(row["blank_corrected_a260"])
            self.assertIsNotNone(row["blank_mean_a260"])
          else:
            self.assertIsNotNone(row["blank_corrected_a260"])
        self.assertEqual(script.json_safe(result.data)[0][8], None if value is None else str(value))

  def test_unexpected_matrix_or_optical_mode_is_rejected(self) -> None:
    """Never apply well labels to a transposed matrix or a different optical measurement."""
    for result in (
      replace(measurement(), data=[[0.0] * 8 for _ in range(12)]),
      replace(measurement(), wavelengths=(450,)),
      replace(measurement(), reference_subtracted=True),
    ):
      with self.subTest(result=result), self.assertRaises(ValueError):
        script.result_rows(result, script.load_manifest())


class AcquisitionTests(unittest.IsolatedAsyncioTestCase):
  """Run the physical workflow exclusively against mocked reader methods."""

  def setUp(self) -> None:
    """Create a temporary output and the actual offline catalog and plate map."""
    self.temp = tempfile.TemporaryDirectory()
    self.addCleanup(self.temp.cleanup)
    self.args = argparse.Namespace(
      output=Path(self.temp.name) / "run",
      port="FAKE",
      expected_serial="1191",
      catalog=script.DEFAULT_CATALOG,
      geometry=None,
    )
    self.catalog = FilterSlideCatalog.from_softmax_xml(script.DEFAULT_CATALOG)
    self.records = script.load_manifest()
    self.geometry = script.load_geometry(None)
    self.reader = fake_reader()

  async def execute(self, answers: list[str]) -> None:
    """Call acquisition with a mocked reader constructor and explicit test responses."""
    with (
      patch.object(script, "FilterMaxF5", return_value=self.reader),
      contextlib.redirect_stdout(io.StringIO()),
    ):
      plan = script.build_plan(self.args, self.catalog, self.geometry, self.records)
      await script.acquire(
        self.args, self.catalog, self.geometry, self.records, plan, Mock(side_effect=answers)
      )

  def metadata(self) -> dict:
    """Read this test's acquisition metadata."""
    return json.loads((self.args.output / "reader_metadata.json").read_text())

  async def test_rejected_approval_prevents_construction(self) -> None:
    """No approval means no device object or output directory."""
    with (
      patch.object(script, "FilterMaxF5") as factory,
      contextlib.redirect_stdout(io.StringIO()),
      self.assertRaisesRegex(ValueError, "not approved"),
    ):
      await script.acquire(self.args, self.catalog, self.geometry, self.records, {}, lambda _: "NO")
    factory.assert_not_called()
    self.assertFalse(self.args.output.exists())

  async def test_existing_output_is_preserved(self) -> None:
    """Refuse an existing run before prompting or opening hardware."""
    self.args.output.mkdir()
    confirm = Mock(side_effect=AssertionError("Unexpected prompt"))
    with (
      patch.object(script, "FilterMaxF5") as factory,
      self.assertRaisesRegex(ValueError, "already exists"),
    ):
      await script.acquire(self.args, self.catalog, self.geometry, self.records, {}, confirm)
    factory.assert_not_called()
    confirm.assert_not_called()

  async def test_identity_readiness_and_filter_failure_prevent_motion(self) -> None:
    """Wrong serial, unready reader or wrong installed slide never reaches tray movement."""
    for failure in ("serial", "ready", "slide"):
      with self.subTest(failure=failure):
        self.args.output = Path(self.temp.name) / failure
        self.reader = fake_reader()
        if failure == "serial":
          self.reader.get_instrument_info.return_value = replace(
            self.reader.get_instrument_info.return_value, serial_number="different"
          )
        elif failure == "ready":
          self.reader.get_status.return_value = FilterMaxStatus(False, "BUSY", ())
        else:
          self.reader.get_filter_slides.return_value = InstalledFilterSlides(1, 1)
        with self.assertRaises(Exception):
          await self.execute([script.APPROVAL])
        self.reader.move_plate_tray_out.assert_not_awaited()
        self.reader.read_absorbance.assert_not_awaited()
        self.reader.stop.assert_awaited_once()
        self.assertEqual(self.metadata()["status"], "failed_or_cancelled")

  async def test_loading_cancel_does_not_measure(self) -> None:
    """Stop after a cancelled loading acknowledgement without requesting a read."""
    with self.assertRaisesRegex(ValueError, "Loading cancelled"):
      await self.execute([script.APPROVAL, "CANCEL"])
    self.reader.read_absorbance.assert_not_awaited()
    self.reader.stop.assert_awaited_once()

  async def test_one_endpoint_and_analyzer_compatible_exports(self) -> None:
    """Export all 96 raw readings and separate corrections without changing optical modes."""
    await self.execute([script.APPROVAL, "LOADED"])
    self.reader.read_absorbance.assert_awaited_once_with(
      self.geometry,
      wavelength=260,
      wells=script.WELLS,
      slide_id=2,
      filter_position=1,
    )
    self.reader.stop.assert_awaited_once()
    expected_files = {
      "raw_a260.csv",
      "blank_corrected_a260.csv",
      "blank_summary.csv",
      "raw_backend_results.json",
      "well_manifest.csv",
      "reader_metadata.json",
    }
    self.assertEqual({p.name for p in self.args.output.iterdir()}, expected_files)
    raw_path = self.args.output / "raw_a260.csv"
    original_bytes = raw_path.read_bytes()
    with raw_path.open(newline="") as stream:
      rows = list(csv.DictReader(stream))
    self.assertEqual(len(rows), 96)
    self.assertEqual(rows[-1]["well"], "H12")
    self.assertTrue(all(row["valid"] == row["exclusion_reason"] == "" for row in rows))
    # The existing analyzer requires review; simulate that review on a separate copy.
    with self.assertRaisesRegex(ValueError, "valid=true or false"):
      analyze_a260.load_readings(raw_path)
    reviewed = Path(self.temp.name) / "reviewed.csv"
    script.write_csv(reviewed, [{**row, "valid": "true"} for row in rows])
    loaded = analyze_a260.load_readings(reviewed)
    self.assertEqual(loaded["H12"]["a260"], 7.11)
    self.assertEqual(raw_path.read_bytes(), original_bytes)
    metadata = self.metadata()
    self.assertEqual(metadata["status"], "complete")
    self.assertEqual(metadata["mode"], "PHYSICAL_DNA_A260_READ")
    self.assertEqual(metadata["wavelength"], 260)
    self.assertFalse(metadata["pathlength_correction_applied"])
    self.assertFalse(metadata["reference_subtracted"])
    self.assertEqual(metadata["validity_review"], "pending")
    self.assertEqual(metadata["incomplete_blank_volumes"], [])

  async def test_incomplete_blanks_and_nonfinite_payload_are_preserved(self) -> None:
    """Finish with explicit missing-data status and lossless raw JSON."""
    self.reader.read_absorbance.return_value[0].data[0][8] = float("nan")
    await self.execute([script.APPROVAL, "LOADED"])
    metadata = self.metadata()
    self.assertEqual(metadata["status"], "complete_with_missing_values")
    self.assertEqual(metadata["incomplete_blank_volumes"], [200])
    self.assertEqual(metadata["missing_well_count"], 1)
    raw = json.loads((self.args.output / "raw_backend_results.json").read_text())
    self.assertEqual(raw[0]["data"][0][8], "nan")

  async def test_raw_payload_saved_before_invalid_result_rejection(self) -> None:
    """Capture unexpected results before failing their optical or shape checks."""
    for label, results in (
      ("shape", [replace(measurement(), data=[[1.0]])]),
      ("wavelength", [replace(measurement(), wavelengths=(450,))]),
      ("reference", [replace(measurement(), reference_subtracted=True)]),
      ("count", [measurement(), measurement()]),
    ):
      with self.subTest(label=label):
        self.args.output = Path(self.temp.name) / label
        self.reader = fake_reader()
        self.reader.read_absorbance.return_value = results
        with self.assertRaises(ValueError):
          await self.execute([script.APPROVAL, "LOADED"])
        self.assertTrue((self.args.output / "raw_backend_results.json").exists())
        self.assertFalse((self.args.output / "raw_a260.csv").exists())
        self.assertEqual(self.metadata()["status"], "failed_or_cancelled")
        self.reader.stop.assert_awaited_once()

  async def test_setup_failure_is_closed_and_recorded(self) -> None:
    """A failed connection is closed without reaching motion or measurement."""
    self.reader.setup.side_effect = RuntimeError("setup failed")
    with self.assertRaisesRegex(RuntimeError, "setup failed"):
      await self.execute([script.APPROVAL])
    self.reader.stop.assert_awaited_once()
    self.reader.move_plate_tray_out.assert_not_awaited()
    self.assertIn("setup failed", self.metadata()["error"])

  async def test_primary_error_survives_cleanup_and_metadata_failures(self) -> None:
    """No retry and no masking of the acquisition error when shutdown or final logging fail."""
    self.reader.read_absorbance.side_effect = RuntimeError("read failed")
    self.reader.stop.side_effect = RuntimeError("stop failed")
    with (
      contextlib.redirect_stderr(io.StringIO()),
      self.assertRaisesRegex(RuntimeError, "read failed"),
    ):
      await self.execute([script.APPROVAL, "LOADED"])
    self.reader.read_absorbance.assert_awaited_once()
    metadata = self.metadata()
    self.assertIn("read failed", metadata["error"])
    self.assertIn("stop failed", metadata["shutdown_error"])
    self.args.output = Path(self.temp.name) / "save-failure"
    self.reader = fake_reader()
    self.reader.read_absorbance.side_effect = RuntimeError("read failed")
    write = script.write_json

    def fail_final_save(path: Path, value: dict) -> None:
      """Fail only final metadata persistence, after the primary read error."""
      if path.name == "reader_metadata.json" and "finished_at" in value:
        raise OSError("final save failed")
      write(path, value)

    with (
      patch.object(script, "write_json", side_effect=fail_final_save),
      contextlib.redirect_stderr(io.StringIO()),
      self.assertRaisesRegex(RuntimeError, "read failed"),
    ):
      await self.execute([script.APPROVAL, "LOADED"])
    self.reader.stop.assert_awaited_once()

  async def test_read_cancellation_is_recorded_without_retry(self) -> None:
    """Cancellation is propagated after connection cleanup."""
    self.reader.read_absorbance.side_effect = asyncio.CancelledError()
    with self.assertRaises(asyncio.CancelledError):
      await self.execute([script.APPROVAL, "LOADED"])
    self.reader.read_absorbance.assert_awaited_once()
    self.reader.stop.assert_awaited_once()
    self.assertIn("CancelledError", self.metadata()["error"])


if __name__ == "__main__":
  unittest.main()
