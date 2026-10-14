"""Offline regression tests for selected-well F5 acquisition and concentration analysis."""

from __future__ import annotations

import asyncio
import contextlib
import csv
import io
import json
import tempfile
import unittest
from dataclasses import asdict, replace
from pathlib import Path
from unittest.mock import AsyncMock, Mock, patch

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import nbformat
import pandas as pd

from pylabrobot.molecular_devices.filtermax import (
  AbsorbanceResult,
  FilterMaxF5,
  FilterMaxFilterNotFoundError,
  FilterMaxStatus,
  InstalledFilterSlides,
  InstrumentInfo,
)

import read_dna_a260_selected_wells as a260

NOTEBOOK = Path(__file__).with_name("read_dna_a260_selected_wells.ipynb")


def example_readings() -> pd.DataFrame:
  """Return the user's two samples and minimum-absorbance blank."""
  return pd.DataFrame({"well": ["A9", "A10", "A12"], "a260": [0.08558, 0.0874, 0.05665]})


def fake_reader() -> Mock:
  """Supply public reader methods without constructing any hardware transport."""
  data = [[None for _ in range(12)] for _ in range(8)]
  data[0][8], data[0][9], data[0][11] = 0.08558, 0.0874, 0.05665
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
  reader.read_absorbance = AsyncMock(
    return_value=[AbsorbanceResult(data, (260,), False, 23.5, 1791375773)]
  )
  return reader


class PlanningTests(unittest.TestCase):
  """Verify offline planning, selection validation, and the notebook defaults."""

  def test_well_selection(self) -> None:
    """Normalize case and spacing, retaining order and rejecting ambiguous entries."""
    self.assertEqual(a260.parse_wells(" a10, A9 "), ("A10", "A9"))
    for text in ("", " ", "A8,", "A8,,A9", "A8, a8", "A0", "A13", "I1", "A8:A9"):
      with self.subTest(text=text), self.assertRaises(ValueError):
        a260.parse_wells(text)

  def test_preview_never_constructs_hardware_or_creates_results(self) -> None:
    """Resolve local filters, geometry, and the inclusive row span entirely offline."""
    with (
      tempfile.TemporaryDirectory() as folder,
      patch.object(
        a260, "FilterMaxF5", side_effect=AssertionError("hardware forbidden")
      ) as factory,
    ):
      plan = a260.prepare_read("a8, C1", results_root=Path(folder), notebook_path=NOTEBOOK)
      self.assertEqual(plan.wells, ("A8", "C1"))
      self.assertEqual(plan.metadata["hardware_rows"], ["A", "B", "C"])
      self.assertEqual(plan.metadata["well_count"], 2)
      self.assertEqual(plan.metadata["catalog_filter"]["wavelength"], 260)
      self.assertFalse(plan.output.exists())
      self.assertFalse(plan.metadata["pathlength_correction_applied"])
      self.assertNotIn("assignment_sha256", plan.metadata)
      factory.assert_not_called()

  def test_invalid_geometry_and_missing_filter_fail_offline(self) -> None:
    """Retain plate bounds and exact-wavelength validation."""
    with tempfile.TemporaryDirectory() as folder:
      path = Path(folder) / "geometry.json"
      valid = asdict(a260.load_geometry(None))
      for change in ({"rows": 12}, {"absorbance_z": 100}, {"column_spacing": 99}):
        path.write_text(json.dumps({**valid, **change}))
        with self.subTest(change=change), self.assertRaises(ValueError):
          a260.load_geometry(path)
      catalog = Path(folder) / "filters.xml"
      catalog.write_text(a260.DEFAULT_CATALOG.read_text().replace("<WL>260</WL>", "<WL>261</WL>"))
      with self.assertRaises(FilterMaxFilterNotFoundError):
        a260.prepare_read("A9, A12", catalog_path=catalog)

  def test_notebook_structure_and_safe_default(self) -> None:
    """Ship a valid unexecuted notebook with a separate table and plot cell."""
    notebook = nbformat.read(NOTEBOOK, as_version=4)
    nbformat.validate(notebook)
    cells = {cell.id: cell for cell in notebook.cells}
    self.assertIn("RUN_READER = False", cells["run-read"].source)
    self.assertIn("save_analysis", cells["show-results"].source)
    self.assertIn("plot_concentrations", cells["plot-results"].source)
    for cell in notebook.cells:
      if cell.cell_type == "code":
        self.assertIsNone(cell.execution_count)
        self.assertEqual(cell.outputs, [])


class AnalysisTests(unittest.TestCase):
  """Verify the supplied numerical example and the plotted calibration."""

  def test_example_and_display_rounding(self) -> None:
    """Round only for display and apply the dilution factor after inversion."""
    table, metadata = a260.analyze_readings(example_readings())
    self.assertEqual(metadata["blank_well"], "A12")
    self.assertEqual(metadata["blank_a260"], 0.05665)
    self.assertEqual(metadata["slope"], 0.00338)
    self.assertEqual(metadata["intercept"], 0.0022)
    self.assertEqual(metadata["dilution_factor"], 10)
    self.assertFalse(metadata["additional_volume_correction_applied"])
    by_well = table.set_index("well")
    for well, corrected, concentration in (
      ("A9", 0.02893, 7.908284023668639),
      ("A10", 0.03075, 8.446745562130177),
    ):
      self.assertAlmostEqual(by_well.loc[well, "blank_corrected_a260"], corrected)
      self.assertAlmostEqual(by_well.loc[well, "calculated_concentration_ng_uL"], concentration)
      self.assertAlmostEqual(by_well.loc[well, "reported_concentration_ng_uL"], concentration * 10)
    self.assertEqual(by_well.loc["A12", "blank_corrected_a260"], 0)
    self.assertTrue(pd.isna(by_well.loc["A12", "calculated_concentration_ng_uL"]))
    displayed = a260.concentration_display_table(table)
    self.assertEqual(displayed.loc["A9", "Calculated concentration (ng/ul)"], "7.91")
    self.assertEqual(displayed.loc["A9", "Reported concentration ×10 (ng/ul)"], "79.1")
    self.assertEqual(displayed.loc["A10", "Reported concentration ×10 (ng/ul)"], "84.5")
    self.assertEqual(displayed.loc["A12", "Reported concentration ×10 (ng/ul)"], "—")

  def test_tie_uses_entered_order_and_retains_signed_estimates(self) -> None:
    """Do not sort tied blanks by well or clip concentrations below the intercept."""
    table, metadata = a260.analyze_readings(
      pd.DataFrame({"well": ["B1", "A1", "A2"], "a260": [0.05, 0.05, 0.051]})
    )
    self.assertEqual(metadata["blank_well"], "B1")
    self.assertEqual(table["well"].tolist(), ["B1", "A1", "A2"])
    self.assertTrue((table.loc[table.role == "sample", "calculated_concentration_ng_uL"] < 0).all())

  def test_missing_values_do_not_become_blank_or_disappear(self) -> None:
    """Keep invalid rows available for review while using only finite values."""
    raw = pd.DataFrame(
      {
        "well": ["A1", "A2", "A3", "A4", "A5", "A6"],
        "a260": [None, float("nan"), float("inf"), float("-inf"), 0.04, 0.06],
      }
    )
    table, metadata = a260.analyze_readings(raw)
    self.assertEqual(metadata["blank_well"], "A5")
    self.assertEqual(len(table), 6)
    self.assertEqual(table.role.tolist()[:4], ["unavailable"] * 4)
    self.assertTrue(table["reported_concentration_ng_uL"].iloc[:4].isna().all())

  def test_insufficient_finite_readings_are_rejected(self) -> None:
    """A blank alone cannot produce a sample concentration."""
    for values in ([None, float("inf")], [0.02, None], [0.02]):
      raw = pd.DataFrame({"well": [f"A{i + 1}" for i in range(len(values))], "a260": values})
      with self.subTest(values=values), self.assertRaisesRegex(ValueError, "two finite"):
        a260.analyze_readings(raw)

  def test_saved_analysis_keeps_raw_file_and_precision(self) -> None:
    """Save derived data and its blank/calibration provenance separately."""
    with tempfile.TemporaryDirectory() as folder:
      path = Path(folder) / "raw_a260.csv"
      example_readings().to_csv(path, index=False)
      original = path.read_bytes()
      table, metadata = a260.save_analysis(path)
      self.assertEqual(path.read_bytes(), original)
      saved = pd.read_csv(path.parent / "concentrations.csv")
      self.assertAlmostEqual(saved.loc[0, "reported_concentration_ng_uL"], 79.08284023668639)
      self.assertEqual(len(saved), len(table))
      persisted = json.loads((path.parent / "analysis_metadata.json").read_text())
      self.assertEqual(persisted["blank_well"], "A12")
      self.assertEqual(persisted["raw_csv_sha256"], metadata["raw_csv_sha256"])

  def test_plot_uses_calibration_x_and_reported_labels(self) -> None:
    """Place points on the supplied line and label only the samples with 10x values."""
    table, _ = a260.analyze_readings(example_readings())
    fig = a260.plot_concentrations(table)
    self.addCleanup(plt.close, fig)
    ax = fig.axes[0]
    self.assertEqual(ax.get_title(), a260.PLOT_TITLE)
    self.assertEqual(ax.get_xlabel(), "Concentration (ng/ul)")
    self.assertEqual(ax.get_ylabel(), "Blank-corrected A260 (OD)")
    points = ax.collections[0].get_offsets()
    self.assertEqual(len(points), 2)
    self.assertAlmostEqual(points[0][0], 7.908284023668639)
    self.assertAlmostEqual(points[0][1], 0.02893)
    labels = [text.get_text() for text in ax.texts]
    self.assertIn("A9\n79.1ng/ul", labels)
    self.assertIn("A10\n84.5ng/ul", labels)
    self.assertFalse(any("A12" in label for label in labels))
    self.assertIn("y = 0.00338x + 0.0022\nSlope = 0.00338\ny-intercept = 0.0022", labels)
    for x, y in zip(ax.lines[0].get_xdata(), ax.lines[0].get_ydata()):
      self.assertAlmostEqual(y, a260.SLOPE * x + a260.INTERCEPT)
    with tempfile.TemporaryDirectory() as folder:
      path = Path(folder) / "plot.png"
      fig.savefig(path)
      self.assertGreater(path.stat().st_size, 1000)


class AcquisitionTests(unittest.IsolatedAsyncioTestCase):
  """Exercise the read lifecycle only against mocked reader methods."""

  def setUp(self) -> None:
    """Prepare a fresh temporary output and intercept every reader construction."""
    self.temp = tempfile.TemporaryDirectory()
    self.addCleanup(self.temp.cleanup)
    self.root = Path(self.temp.name)
    self.plan = a260.prepare_read("A10, A9, A12", port="FAKE", results_root=self.root)
    self.reader = fake_reader()
    self.factory_patch = patch.object(a260, "FilterMaxF5", return_value=self.reader)
    self.factory = self.factory_patch.start()
    self.addCleanup(self.factory_patch.stop)
    self.confirm = Mock(return_value="LOADED")

  async def acquire(self) -> Path:
    """Execute one mocked read without printing its operator instructions."""
    with contextlib.redirect_stdout(io.StringIO()):
      return await a260.acquire_selected(self.plan, self.confirm)

  def metadata(self) -> dict:
    """Read this test's acquisition record."""
    return json.loads((self.plan.output / "reader_metadata.json").read_text())

  async def test_only_loaded_prompt_and_correct_well_mapping(self) -> None:
    """Run on request, with exactly one confirmation after tray opening."""

    def loaded(prompt: str) -> str:
      """Check that the loading prompt follows tray motion."""
      self.reader.move_plate_tray_out.assert_awaited_once()
      self.reader.read_absorbance.assert_not_awaited()
      self.assertIn("LOADED", prompt)
      self.assertNotIn("APPROVE", prompt)
      return "LOADED"

    self.confirm.side_effect = loaded
    output = await self.acquire()
    self.confirm.assert_called_once()
    self.reader.read_absorbance.assert_awaited_once_with(
      self.plan.geometry, wavelength=260, wells=("A10", "A9", "A12"), slide_id=2, filter_position=1
    )
    raw = pd.read_csv(output / "raw_a260.csv")
    self.assertEqual(raw.well.tolist(), ["A10", "A9", "A12"])
    self.assertEqual(raw.a260.tolist(), [0.0874, 0.08558, 0.05665])
    self.assertEqual(self.metadata()["status"], "complete")
    self.assertEqual(self.metadata()["missing_well_count"], 0)
    self.assertEqual(self.metadata()["well_count"], 3)
    self.assertNotIn("source_tube", raw.columns)
    self.reader.stop.assert_awaited_once()

  async def test_cancelled_loading_does_not_measure(self) -> None:
    """A non-LOADED response cancels measurement and closes the connection."""
    self.confirm.return_value = "CANCEL"
    with self.assertRaisesRegex(ValueError, "Loading cancelled"):
      await self.acquire()
    self.reader.read_absorbance.assert_not_awaited()
    self.reader.stop.assert_awaited_once()
    self.assertEqual(self.metadata()["status"], "failed_or_cancelled")

  async def test_existing_output_never_connects_or_prompts(self) -> None:
    """Preserve an existing run instead of overwriting it."""
    self.plan.output.mkdir()
    with self.assertRaisesRegex(ValueError, "already exists"):
      await self.acquire()
    self.factory.assert_not_called()
    self.confirm.assert_not_called()

  async def test_identity_readiness_and_filter_checks_prevent_motion(self) -> None:
    """Reject wrong serials, busy readers, and incompatible slides before tray movement."""
    for failure in ("serial", "ready", "slide"):
      with self.subTest(failure=failure):
        self.plan = a260.prepare_read("A9, A12", port="FAKE", results_root=self.root)
        self.reader = fake_reader()
        self.factory.return_value = self.reader
        if failure == "serial":
          self.reader.get_instrument_info.return_value = replace(
            self.reader.get_instrument_info.return_value, serial_number="wrong-reader"
          )
        elif failure == "ready":
          self.reader.get_status.return_value = FilterMaxStatus(False, "BUSY", ())
        else:
          self.reader.get_filter_slides.return_value = InstalledFilterSlides(1, 1)
        with self.assertRaises(Exception):
          await self.acquire()
        self.reader.move_plate_tray_out.assert_not_awaited()
        self.reader.read_absorbance.assert_not_awaited()
        self.reader.stop.assert_awaited_once()
        self.confirm.assert_not_called()

  async def test_nonfinite_payload_is_saved_before_analysis(self) -> None:
    """Retain raw backend values even when a selected well is missing."""
    self.reader.read_absorbance.return_value[0].data[0][8] = float("nan")
    output = await self.acquire()
    self.assertEqual(self.metadata()["missing_well_count"], 1)
    self.assertEqual(self.metadata()["status"], "complete_with_missing_values")
    raw = json.loads((output / "raw_backend_results.json").read_text())
    self.assertEqual(raw[0]["data"][0][8], "nan")
    with (output / "raw_a260.csv").open(newline="") as stream:
      self.assertEqual(list(csv.DictReader(stream))[1]["a260"], "")

  async def test_invalid_result_is_preserved_before_validation(self) -> None:
    """Keep the payload when its shape, wavelength, mode, or count is invalid."""
    result = self.reader.read_absorbance.return_value[0]
    for results in (
      [replace(result, data=[[1.0]])],
      [replace(result, wavelengths=(450,))],
      [replace(result, reference_subtracted=True)],
      [result, result],
    ):
      with self.subTest(results=results):
        self.plan = a260.prepare_read("A9, A12", port="FAKE", results_root=self.root)
        self.reader.read_absorbance.return_value = results
        with self.assertRaises(ValueError):
          await self.acquire()
        self.assertTrue((self.plan.output / "raw_backend_results.json").exists())
        self.assertFalse((self.plan.output / "raw_a260.csv").exists())
        self.assertEqual(self.metadata()["status"], "failed_or_cancelled")

  async def test_setup_failure_is_closed_without_motion(self) -> None:
    """Disconnect after a failed setup and record its cause."""
    self.reader.setup.side_effect = RuntimeError("setup failure")
    with self.assertRaisesRegex(RuntimeError, "setup failure"):
      await self.acquire()
    self.reader.move_plate_tray_out.assert_not_awaited()
    self.reader.stop.assert_awaited_once()
    self.assertIn("setup failure", self.metadata()["error"])

  async def test_primary_error_survives_cleanup_failure(self) -> None:
    """Do not retry or replace the acquisition error with a shutdown error."""
    self.reader.read_absorbance.side_effect = RuntimeError("read failure")
    self.reader.stop.side_effect = RuntimeError("stop failure")
    with (
      contextlib.redirect_stderr(io.StringIO()),
      self.assertRaisesRegex(RuntimeError, "read failure"),
    ):
      await self.acquire()
    self.reader.read_absorbance.assert_awaited_once()
    self.reader.stop.assert_awaited_once()
    self.assertIn("read failure", self.metadata()["error"])
    self.assertIn("stop failure", self.metadata()["shutdown_error"])

  async def test_cancellation_is_recorded_and_propagated(self) -> None:
    """Cancellation closes the connection without retrying a read."""
    self.reader.read_absorbance.side_effect = asyncio.CancelledError()
    with self.assertRaises(asyncio.CancelledError):
      await self.acquire()
    self.reader.read_absorbance.assert_awaited_once()
    self.reader.stop.assert_awaited_once()
    self.assertIn("CancelledError", self.metadata()["error"])


if __name__ == "__main__":
  unittest.main()
