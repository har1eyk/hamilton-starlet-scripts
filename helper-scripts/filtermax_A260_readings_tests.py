"""Offline tests for the interactive FilterMax UV script."""

import asyncio
import contextlib
import importlib.util
import io
import math
import sys
import unittest
from dataclasses import replace
from pathlib import Path
from typing import Iterator
from unittest.mock import AsyncMock, MagicMock, patch

from pylabrobot.molecular_devices.filtermax import (
  AbsorbanceResult,
  FilterDefinition,
  FilterMaxCatalogError,
  FilterMaxF5,
  FilterMaxStatus,
  FilterMaxTimeoutError,
  FilterSlideCatalog,
  InstalledFilterSlides,
  InstrumentInfo,
)
from pylabrobot.molecular_devices.filtermax.protocol import ACK, ENQ, EOT, build_frame

SPEC = importlib.util.spec_from_file_location(
  "filtermax_uv_script", Path(__file__).with_name("filtermax-A260-readings.py")
)
assert SPEC is not None and SPEC.loader is not None
script = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(script)


def catalog_fixture() -> FilterSlideCatalog:
  """Construct a synthetic catalog with A260 and no other UV absorbance wavelengths."""
  primary = FilterDefinition(
    kind="excitation",
    slide_id=7,
    slide_name="Synthetic UV",
    position=1,
    wavelength=260,
    bandwidth=10,
    apply_to=1,
    techniques=frozenset({"absorbance"}),
  )
  return FilterSlideCatalog(
    (
      primary,
      replace(primary, kind="emission", slide_id=9, slide_name="Synthetic emission"),
      replace(primary, slide_id=8, slide_name="Uninstalled set"),
    )
  )


def result_fixture(wavelength: int, base: float = 0.1) -> AbsorbanceResult:
  """Make canonical plate data with only row A populated."""
  return AbsorbanceResult(
    data=[[base + column / 100 for column in range(12)]] + [[None] * 12 for _ in range(7)],
    wavelengths=(wavelength,),
    reference_subtracted=False,
    temperature=24.0,
    timestamp=0.0,
  )


def reader_fixture() -> MagicMock:
  """Make a ready F5 with mocked read-only queries and no open serial connection."""
  reader = MagicMock(spec=FilterMaxF5)
  reader.port = "/dev/ttyUSB0"
  reader.io = FilterMaxF5("offline", reader.port).io
  reader.get_instrument_info.return_value = InstrumentInfo(
    model="Anthos Fluoro",
    firmware="V1.2 b13 11.12.2014",
    device_number=141,
    serial_number="1191",
    excitation_slide_id=7,
    emission_slide_id=9,
    device_code=57855,
    plate_control=True,
    pic_firmware="V1.0",
    cpld_version=6,
    raw="",
  )
  reader.get_status.return_value = FilterMaxStatus(ready=True, text="READY", state_bits=())
  reader.get_filter_slides.return_value = InstalledFilterSlides(7, 9)
  reader.get_temperature.return_value = 25.6
  return reader


def serial_response(*payloads: str) -> bytes:
  """Frame one simulated reply, including multiple measurement messages when needed."""
  return bytes((ACK, ACK)) + b"".join(
    bytes((ENQ,)) + build_frame(payload) + bytes((EOT,)) for payload in payloads
  )


def preflight_transcript() -> bytes:
  """Supply notebook identity and status replies for the synthetic filter slides."""
  identity = (
    "+\rAnthos Fluoro\rV1.2 b13 11.12.2014\r141\r1191\r"
    "Excitation Filter Slider Nr: 7\rEmission Filter Slider Nr: 9\r"
    "Barcode: not installed\rTemperatureControl: installed\rDevice Code: 57855\r"
    "Plate Control: TRUE\rPIC Firmware: + V1.0\rCPLD Version: 6"
  )
  return b"".join(
    serial_response(payload)
    for payload in (identity, identity, "+ READY", "+ 0 0 0 0 0 0 0 0 0 0 0", "+ 7 9", "+ 25.6")
  )


def absorbance_transcript(base: int) -> bytes:
  """Supply preparation, one row of milli-OD data, and temperature without ejection."""
  preparation = b"".join(
    serial_response(payload) for payload in ("+ 7 9", "+", "+ 0 1", "+", "+ -67")
  )
  row = "+ " + " ".join(str(base + column) for column in range(12))
  return preparation + serial_response(row, "+ INFO") + serial_response("+ 25.6")


class SerialSimulation:
  """Supply per-baud transcripts beneath the production driver and protocol parser."""

  def __init__(self, reader: FilterMaxF5, responses: dict[int, bytes]):
    """Record the production transport's settings without opening a serial device."""
    self.reader = reader
    self.responses = responses
    self.incoming = bytearray()
    self.events: list[tuple[str, int]] = []
    self.writes: list[tuple[int, bytes]] = []
    self.is_open = False

  async def open(self) -> None:
    """Require close before reopening and preserve the notebook's serial settings."""
    assert not self.is_open
    serial = self.reader.io
    assert (serial.bytesize, serial.parity, serial.stopbits) == (7, "E", 1)
    assert not serial.rtscts and not serial.dsrdtr and not serial.xonxoff
    self.events.append(("open", serial.baudrate))
    self.incoming = bytearray(self.responses[serial.baudrate])
    self.is_open = True

  async def close(self) -> None:
    """Record closure without changing the reader's selected baud rate."""
    self.events.append(("close", self.reader.io.baudrate))
    self.is_open = False

  async def write(self, data: bytes) -> None:
    """Capture the exact outgoing bytes and their transmission rate."""
    assert self.is_open
    self.writes.append((self.reader.io.baudrate, data))

  async def read(self, num_bytes: int = 1) -> bytes:
    """Return recorded bytes or expire a silent handshake without a wall-clock delay."""
    assert self.is_open
    if not self.incoming:
      raise FilterMaxTimeoutError("Simulated serial timeout")
    data = bytes(self.incoming[:num_bytes])
    del self.incoming[:num_bytes]
    return data

  @property
  def payloads(self) -> list[str]:
    """Decode transmitted command frames for independent workflow assertions."""
    return [data[2:-5].decode("ascii") for _, data in self.writes if data.startswith(b"\x02")]


@contextlib.contextmanager
def simulate_serial(reader: FilterMaxF5, responses: dict[int, bytes]) -> Iterator[SerialSimulation]:
  """Replace every serial I/O entry point while retaining production protocol code."""
  simulation = SerialSimulation(reader, responses)
  with (
    patch.object(reader.io, "setup", side_effect=simulation.open),
    patch.object(reader.io, "stop", side_effect=simulation.close),
    patch.object(reader.io, "write", side_effect=simulation.write),
    patch.object(reader.io, "read", side_effect=simulation.read),
  ):
    yield simulation


class TestSerialIntegration(unittest.IsolatedAsyncioTestCase):
  """Exercise the complete script and driver together without hardware or a display."""

  async def test_connection_check_reports_actual_baud_without_motion(self) -> None:
    """Either supported baud must complete the notebook's exact query sequence."""
    for baudrate in (38400, 9600):
      with self.subTest(baudrate=baudrate):
        reader = FilterMaxF5("offline", script.DEFAULT_PORT)
        with (
          simulate_serial(reader, {38400: b"", baudrate: preflight_transcript()}) as serial,
          contextlib.redirect_stdout(io.StringIO()) as output,
        ):
          await script.run_connection_check(reader)
        self.assertEqual(serial.payloads, ["?", "?", "STAT", "SGS", "CS", "TG"])
        self.assertIn(f"Connected serial settings: {baudrate} baud, 7E1", output.getvalue())
        self.assertIn("serial 1191 is READY", output.getvalue())
        self.assertEqual(serial.events[-1], ("close", baudrate))
        self.assertFalse(serial.is_open)
        self.assertFalse(serial.incoming)
        if baudrate == 9600:
          self.assertEqual(serial.events[:3], [("open", 38400), ("close", 38400), ("open", 9600)])
          self.assertEqual([data for rate, data in serial.writes if rate == 38400], [bytes((ENQ,))])

  async def test_a260_read_maps_all_values_and_ejects(self) -> None:
    """The real driver scans row A once at 260 nm, then ejects the plate."""
    for baudrate in (38400, 9600):
      with self.subTest(baudrate=baudrate):
        catalog = catalog_fixture()
        reader = FilterMaxF5("offline", script.DEFAULT_PORT, catalog)
        incoming = (
          preflight_transcript()
          + serial_response("+ 6512")
          + serial_response("+ -67")
          + absorbance_transcript(1000)
          + serial_response("+ 6512")
        )
        with (
          simulate_serial(reader, {38400: b"", baudrate: incoming}) as serial,
          contextlib.redirect_stdout(io.StringIO()) as output,
        ):

          def confirm() -> bool:
            """Confirm only after opening the tray, before loading or scanning."""
            self.assertEqual(serial.payloads, ["?", "?", "STAT", "SGS", "CS", "TG", "E P"])
            return True

          await script.run_scan(reader, catalog, confirm)
        commands = serial.payloads
        self.assertEqual(
          [command for command in commands if command.startswith("ABS ")],
          ["ABS 0 1 260 1 1 1 1 0 0 0 0 0 4 1 1 0 O e INFO"],
        )
        plate_command = (
          "PLATE Temp 12770 8570 1427 1069 1405 1118 12 8 902 900 640 640 0 300 0 1027 1 1"
        )
        self.assertEqual(commands.count(plate_command), 1)
        self.assertEqual(commands.count("SHIFT"), 1)
        self.assertEqual(commands.count("E P"), 2)
        self.assertEqual(commands[-1], "E P")
        self.assertNotIn("STOP", commands)
        table = output.getvalue().split("Raw absorbance (unblanked).\n")[1].splitlines()[2:14]
        self.assertEqual(len(table), 12)
        for column, row in enumerate(table):
          values = [cell.strip() for cell in row.split("|")]
          self.assertEqual(values, [f"A{column + 1}", f"{(1000 + column) / 1000:.6f}"])
        self.assertFalse(serial.is_open)
        self.assertFalse(serial.incoming)

  async def test_silent_reader_reports_both_rates_and_initial_stage(self) -> None:
    """The CLI's session wrapper retains useful driver diagnostics on failed setup."""
    reader = FilterMaxF5("offline", script.DEFAULT_PORT)
    with simulate_serial(reader, {38400: b"", 9600: b""}) as serial:
      with self.assertRaises(RuntimeError) as caught:
        await script.run_connection_check(reader)
    for expected in (script.DEFAULT_PORT, "38400, 9600", "initial acknowledgement"):
      self.assertIn(expected, str(caught.exception))
    self.assertEqual(serial.payloads, [])
    self.assertFalse(serial.is_open)


class TestFilterChecks(unittest.TestCase):
  """Filter availability must be tied to the installed slide and technique."""

  def setUp(self) -> None:
    """Prepare a complete synthetic catalog."""
    self.catalog = catalog_fixture()
    self.slides = InstalledFilterSlides(excitation_slide_id=7, emission_slide_id=9)

  def test_resolves_only_installed_excitation_absorbance_filters(self) -> None:
    """An uninstalled duplicate must not cause ambiguity."""
    selections = script.resolve_filters(self.catalog, self.slides)
    self.assertEqual(tuple(selections), (260,))
    self.assertEqual([selection.slide_id for selection in selections.values()], [7])
    self.assertEqual([selection.position for selection in selections.values()], [1])

  def test_inventory_distinguishes_installed_and_uninstalled_sets(self) -> None:
    """Inventory includes both bays, wavelengths, techniques, and stored sets."""
    table = script.filter_inventory(self.catalog, self.slides)
    self.assertIn("excitation=7, emission=9", table)
    self.assertIn("absorbance", table)
    stored = next(line for line in table.splitlines() if "Uninstalled set" in line)
    self.assertTrue(stored.rstrip().endswith("no"))
    installed = next(line for line in table.splitlines() if "Synthetic UV" in line)
    self.assertTrue(installed.rstrip().endswith("yes"))

  def test_reports_missing_a260(self) -> None:
    """Other wavelengths cannot substitute for the installed A260 filter."""
    catalog = FilterSlideCatalog(
      (replace(self.catalog.filters[0], wavelength=450), *self.catalog.filters[1:])
    )
    with self.assertRaises(FilterMaxCatalogError) as raised:
      script.resolve_filters(catalog, self.slides)
    self.assertIn("A260", str(raised.exception))

  def test_wrong_technique_or_uninstalled_filter_does_not_count(self) -> None:
    """A wavelength in a different set or detection mode is unavailable."""
    for replacement in (
      replace(self.catalog.filters[0], techniques=frozenset({"fluorescence"}), apply_to=2),
      replace(self.catalog.filters[0], slide_id=8),
    ):
      with self.subTest(replacement=replacement):
        catalog = FilterSlideCatalog((replacement, *self.catalog.filters[1:]))
        with self.assertRaisesRegex(FilterMaxCatalogError, "A260"):
          script.resolve_filters(catalog, self.slides)

  def test_unknown_installed_slide_is_rejected(self) -> None:
    """Unknown IDs in either bay require an accurate catalog."""
    for slides in (InstalledFilterSlides(99, 9), InstalledFilterSlides(7, 99)):
      with self.subTest(slides=slides):
        with self.assertRaisesRegex(FilterMaxCatalogError, "Unknown installed"):
          script.resolve_filters(self.catalog, slides)

  def test_duplicate_installed_wavelength_is_rejected(self) -> None:
    """Do not silently select between two matching slots."""
    catalog = FilterSlideCatalog(
      (*self.catalog.filters, replace(self.catalog.filters[0], position=4))
    )
    with self.assertRaisesRegex(FilterMaxCatalogError, "ambiguous"):
      script.resolve_filters(catalog, self.slides)


class TestResults(unittest.TestCase):
  """Verify the A260 table, well mapping, and unavailable readings."""

  def test_mapping_and_a260_only_columns(self) -> None:
    """Round only for display, and keep all twelve wells in order."""
    results = {wl: result_fixture(wl) for wl in script.WAVELENGTHS}
    results[260].data[0][0] = 0.1234564
    table = script.results_table(results)
    self.assertEqual(
      [cell.strip() for cell in table.splitlines()[0].split("|")], ["Well", "A260 (OD)"]
    )
    rows = [line.split("|") for line in table.splitlines()[2:]]
    self.assertEqual([row[0].strip() for row in rows], list(script.WELLS))
    self.assertEqual(rows[0][1].strip(), "0.123456")
    self.assertTrue(all(len(row) == 2 for row in rows))
    self.assertEqual(rows[-1][1].strip(), "0.210000")

  def test_unavailable_values(self) -> None:
    """Missing and non-finite readings show N/A without fabricating values."""
    for value in (None, math.nan, math.inf, -math.inf):
      with self.subTest(value=value):
        a260 = result_fixture(260)
        a260.data[0][0] = value
        row = script.results_table({260: a260}).splitlines()[2].split("|")
        self.assertEqual(row[-1].strip(), "N/A")

  def test_zero_is_a_valid_reading(self) -> None:
    """A zero OD should display zero, not unavailable."""
    table = script.results_table({260: result_fixture(260, 0)})
    self.assertEqual(table.splitlines()[2].split("|")[-1].strip(), "0.000000")


class TestWorkflow(unittest.IsolatedAsyncioTestCase):
  """Exercise workflow ordering without opening a serial connection."""

  def setUp(self) -> None:
    """Mock all hardware methods and capture terminal output."""
    self.reader = reader_fixture()
    self.reader.read_absorbance.side_effect = [[result_fixture(wl)] for wl in script.WAVELENGTHS]
    self.catalog = catalog_fixture()
    self.output = io.StringIO()
    self.capture = contextlib.redirect_stdout(self.output)
    self.capture.__enter__()
    self.addCleanup(self.capture.__exit__, None, None, None)

  async def test_identity_timeout_identifies_the_port_without_motion(self) -> None:
    """Report connection failure before trying to inspect slides or move the tray."""
    self.reader.setup.side_effect = FilterMaxTimeoutError("No protocol byte")
    with self.assertRaisesRegex(RuntimeError, "identity on /dev/ttyUSB0"):
      await script.run_scan(self.reader, self.catalog, lambda: True)
    self.reader.get_filter_slides.assert_not_awaited()
    self.reader.move_plate_tray_out.assert_not_awaited()
    self.reader.read_absorbance.assert_not_awaited()

  async def test_slide_timeout_is_distinguished_from_identity_failure(self) -> None:
    """A connected reader may time out at a later step; release its connection."""
    self.reader.get_filter_slides.side_effect = FilterMaxTimeoutError("No protocol byte")
    with self.assertRaisesRegex(RuntimeError, "querying filter slides on /dev/ttyUSB0"):
      await script.run_scan(self.reader, self.catalog, lambda: True)
    self.reader.move_plate_tray_out.assert_not_awaited()
    self.reader.stop.assert_awaited_once()

  async def test_unexpected_serial_blocks_scan_before_motion(self) -> None:
    """Measurements use the same expected instrument identity as connection checks."""
    self.reader.get_instrument_info.return_value = replace(
      self.reader.get_instrument_info.return_value, serial_number="9999"
    )
    confirm = MagicMock(return_value=True)
    with self.assertRaisesRegex(RuntimeError, "Expected F5 serial 1191, received 9999"):
      await script.run_scan(self.reader, self.catalog, confirm)
    confirm.assert_not_called()
    self.reader.get_status.assert_not_awaited()
    self.reader.move_plate_tray_out.assert_not_awaited()
    self.reader.read_absorbance.assert_not_awaited()
    self.reader.stop.assert_awaited_once()

  async def test_okay_closes_tray_before_single_a260_read(self) -> None:
    """No read or retraction is allowed until the confirmation returns true."""

    def confirm() -> bool:
      self.reader.move_plate_tray_out.assert_awaited_once()
      self.reader.move_plate_tray_in.assert_not_awaited()
      self.reader.read_absorbance.assert_not_awaited()
      return True

    await script.run_scan(self.reader, self.catalog, confirm)
    calls = self.reader.read_absorbance.await_args_list
    self.assertEqual([call.kwargs["wavelength"] for call in calls], [260])
    self.assertEqual([call.kwargs["eject_after_read"] for call in calls], [True])
    for call in calls:
      self.assertEqual(call.args[0].orientation, "portrait")
      self.assertEqual(call.kwargs["wells"], tuple(f"A{i}" for i in range(1, 13)))
      self.assertEqual(call.kwargs["slide_id"], 7)
    method_names = [call[0] for call in self.reader.mock_calls]
    self.assertLess(method_names.index("move_plate_tray_in"), method_names.index("read_absorbance"))
    self.reader.stop.assert_awaited_once()
    self.assertIn(script.READY_MESSAGE, self.output.getvalue())
    self.assertIn("A260 (OD)", self.output.getvalue())
    self.assertNotIn("A280", self.output.getvalue())
    self.assertNotIn("A230", self.output.getvalue())
    self.assertNotIn("PARTIAL", self.output.getvalue())

  async def test_cancel_does_not_close_tray_or_read(self) -> None:
    """Cancel and window-close share the same false confirmation path."""
    await script.run_scan(self.reader, self.catalog, lambda: False)
    self.reader.move_plate_tray_out.assert_awaited_once()
    self.reader.move_plate_tray_in.assert_not_awaited()
    self.reader.read_absorbance.assert_not_awaited()
    self.reader.stop.assert_awaited_once()

  async def test_missing_filters_prevent_motion_and_confirmation(self) -> None:
    """A failed preflight must not present a ready-to-scan prompt."""
    confirm = MagicMock(return_value=True)
    self.catalog = FilterSlideCatalog(
      (replace(self.catalog.filters[0], wavelength=450), *self.catalog.filters[1:])
    )
    with self.assertRaises(FilterMaxCatalogError):
      await script.run_scan(self.reader, self.catalog, confirm)
    confirm.assert_not_called()
    self.reader.move_plate_tray_out.assert_not_awaited()
    self.reader.move_plate_tray_in.assert_not_awaited()
    self.reader.read_absorbance.assert_not_awaited()
    self.reader.stop.assert_awaited_once()
    self.assertNotIn(script.READY_MESSAGE, self.output.getvalue())

  async def test_failed_scan_reports_unavailable_readings_and_disconnects(self) -> None:
    """A failed endpoint reports unavailable A260 readings and releases the port."""
    self.reader.read_absorbance.side_effect = RuntimeError("Read failed")
    with self.assertRaisesRegex(RuntimeError, "Read failed"):
      await script.run_scan(self.reader, self.catalog, lambda: True)
    self.reader.read_absorbance.assert_awaited_once()
    self.reader.stop.assert_awaited_once()
    output = self.output.getvalue()
    self.assertIn("PARTIAL RESULTS", output)
    self.assertNotIn("0.100000", output)
    self.assertIn("N/A", output)

  async def test_cancellation_releases_connection_and_reports_partial_table(self) -> None:
    """KeyboardInterrupt's task cancellation still closes the connection."""
    self.reader.read_absorbance.side_effect = asyncio.CancelledError()
    with self.assertRaises(asyncio.CancelledError):
      await script.run_scan(self.reader, self.catalog, lambda: True)
    self.reader.stop.assert_awaited_once()
    self.assertIn("PARTIAL RESULTS", self.output.getvalue())

  async def test_malformed_result_is_rejected(self) -> None:
    """Validate endpoint count, geometry, wavelength, and subtraction mode."""
    invalid = result_fixture(260)
    for measurements in (
      [],
      [invalid, invalid],
      [replace(invalid, data=[[1.0]])],
      [replace(invalid, wavelengths=(450,))],
      [replace(invalid, reference_subtracted=True)],
    ):
      with self.subTest(measurements=measurements):
        self.reader.read_absorbance.reset_mock()
        self.reader.read_absorbance.side_effect = [measurements]
        with self.assertRaises(RuntimeError):
          await script.run_scan(self.reader, self.catalog, lambda: True)
        self.reader.read_absorbance.assert_awaited_once()


class TestConnectionCheck(unittest.IsolatedAsyncioTestCase):
  """The diagnostic mode uses only read-only queries and always releases the port."""

  def setUp(self) -> None:
    """Capture query reports from a mocked reader."""
    self.reader = reader_fixture()
    self.output = io.StringIO()
    self.capture = contextlib.redirect_stdout(self.output)
    self.capture.__enter__()
    self.addCleanup(self.capture.__exit__, None, None, None)

  async def test_queries_match_the_notebook_without_motion(self) -> None:
    """The exact call list excludes tray, filter motion, and measurement methods."""
    await script.run_connection_check(self.reader)
    self.assertEqual(
      [call[0] for call in self.reader.mock_calls],
      [
        "setup",
        "get_instrument_info",
        "get_status",
        "get_filter_slides",
        "get_temperature",
        "stop",
      ],
    )
    self.assertIn("serial 1191 is READY", self.output.getvalue())
    self.assertIn("25.6", self.output.getvalue())

  async def test_unexpected_instrument_serial_fails_and_disconnects(self) -> None:
    """The notebook's expected instrument serial is checked before further queries."""
    self.reader.get_instrument_info.return_value = replace(
      self.reader.get_instrument_info.return_value, serial_number="9999"
    )
    with self.assertRaisesRegex(RuntimeError, "Expected F5 serial 1191, received 9999"):
      await script.run_connection_check(self.reader)
    self.reader.get_status.assert_not_awaited()
    self.reader.stop.assert_awaited_once()

  async def test_not_ready_fails_without_motion(self) -> None:
    """Connecting alone is insufficient for a successful preflight."""
    self.reader.get_status.return_value = FilterMaxStatus(ready=False, text="BUSY", state_bits=())
    with self.assertRaisesRegex(RuntimeError, "not ready: BUSY"):
      await script.run_connection_check(self.reader)
    self.reader.move_plate_tray_out.assert_not_awaited()
    self.reader.read_absorbance.assert_not_awaited()
    self.reader.stop.assert_awaited_once()

  async def test_cleanup_failure_preserves_query_error(self) -> None:
    """Report a failed shutdown separately without replacing the original failure."""
    original = ValueError("status query failed")
    self.reader.get_status.side_effect = original
    self.reader.stop.side_effect = RuntimeError("shutdown failed")
    with contextlib.redirect_stderr(io.StringIO()) as errors:
      with self.assertRaises(ValueError) as raised:
        await script.run_connection_check(self.reader)
    self.assertIs(raised.exception, original)
    self.assertIn("Reader shutdown also failed", errors.getvalue())

  async def test_cleanup_failure_preserves_cancellation(self) -> None:
    """Interrupted queries remain interrupted when closing also fails."""
    original = asyncio.CancelledError()
    self.reader.get_status.side_effect = original
    self.reader.stop.side_effect = RuntimeError("shutdown failed")
    with contextlib.redirect_stderr(io.StringIO()):
      with self.assertRaises(asyncio.CancelledError) as raised:
        await script.run_connection_check(self.reader)
    self.assertIs(raised.exception, original)

  async def test_cleanup_failure_after_success_is_not_hidden(self) -> None:
    """A standalone disconnect error must still produce a failure."""
    self.reader.stop.side_effect = RuntimeError("shutdown failed")
    with self.assertRaisesRegex(RuntimeError, "shutdown failed"):
      await script.run_connection_check(self.reader)


class TestGuiAndEntryPoint(unittest.TestCase):
  """Validate desktop failures and button handling without showing a window."""

  def test_default_adapter_and_explicit_port_override(self) -> None:
    """Use the specific FTDI adapter by default while retaining --port overrides."""
    for extra_args, expected_port in (
      ([], "/dev/serial/by-id/usb-FTDI_US232R_FT4QLERM-if00-port0"),
      (["--port", "/dev/ttyUSB0"], "/dev/ttyUSB0"),
    ):
      with (
        self.subTest(extra_args=extra_args),
        patch.object(script.FilterSlideCatalog, "from_softmax_xml", return_value=catalog_fixture()),
        patch.object(script.FilterMaxF5, "__init__", return_value=None) as construct,
        patch.object(script, "report_environment"),
        patch.object(script, "LoadingDialog"),
        patch.object(script, "run_scan", new_callable=AsyncMock) as scan,
      ):
        self.assertEqual(script.main(["--catalog", "synthetic.xml", *extra_args]), 0)
        self.assertEqual(construct.call_args.kwargs["port"], expected_port)
        scan.assert_awaited_once()

  def test_connection_mode_needs_no_catalog_or_tkinter(self) -> None:
    """Exercise the CLI using a mocked reader with no catalog or graphical display."""
    reader = reader_fixture()
    with (
      patch.dict(sys.modules, {"tkinter": None}),
      patch.object(script, "FilterMaxF5", return_value=reader),
      patch.object(script, "report_environment"),
      patch.object(script.FilterSlideCatalog, "from_softmax_xml") as catalog,
      patch.object(script, "LoadingDialog") as dialog,
      contextlib.redirect_stdout(io.StringIO()) as output,
    ):
      self.assertEqual(script.main(["--check-connection"]), 0)
      catalog.assert_not_called()
      dialog.assert_not_called()
      reader.move_plate_tray_out.assert_not_awaited()
      reader.read_absorbance.assert_not_awaited()
      reader.stop.assert_awaited_once()
      self.assertIn("serial 1191 is READY", output.getvalue())

  def test_scan_mode_requires_catalog_before_hardware_construction(self) -> None:
    """Only connection mode may omit the catalog argument."""
    with (
      patch.object(script, "FilterMaxF5") as construct,
      contextlib.redirect_stderr(io.StringIO()),
    ):
      with self.assertRaises(SystemExit) as raised:
        script.main([])
      self.assertEqual(raised.exception.code, 2)
      construct.assert_not_called()

  def test_debug_enables_serial_logging_and_returns_failure_on_timeout(self) -> None:
    """Diagnostics expose the failing stage and traceback without opening a real port."""
    reader = reader_fixture()
    reader.setup.side_effect = FilterMaxTimeoutError("initial acknowledgement after ENQ")
    with (
      patch.object(script, "FilterMaxF5", return_value=reader),
      patch.object(script, "report_environment"),
      patch.object(script.pylabrobot, "verbose") as verbose,
      contextlib.redirect_stdout(io.StringIO()),
      contextlib.redirect_stderr(io.StringIO()) as errors,
    ):
      self.assertEqual(script.main(["--check-connection", "--debug"]), 1)
      verbose.assert_called_once_with(True, level=script.LOG_LEVEL_IO)
      self.assertIn("initial acknowledgement after ENQ", errors.getvalue())
      self.assertIn("Traceback", errors.getvalue())
      reader.get_instrument_info.assert_not_awaited()

  def test_environment_report_uses_the_configured_transport(self) -> None:
    """Reporting the configured settings must not require opening a port."""
    reader = FilterMaxF5("offline", script.DEFAULT_PORT)
    with contextlib.redirect_stdout(io.StringIO()) as output:
      script.report_environment(reader)
    for expected in ("Python:", "PyLabRobot:", "filtermax_f5.py", "FT4QLERM", "38400 baud, 7E1"):
      self.assertIn(expected, output.getvalue())

  def test_no_display_fails_before_reader_construction(self) -> None:
    """Unavailable Tk must not leave an open tray or serial connection."""
    tkinter = MagicMock(TclError=RuntimeError)

    with (
      patch.dict(sys.modules, {"tkinter": tkinter}),
      patch.object(script.FilterSlideCatalog, "from_softmax_xml", return_value=catalog_fixture()),
      patch.object(script.FilterMaxF5, "__init__", return_value=None) as construct,
      patch.object(tkinter, "Tk", side_effect=tkinter.TclError("no display")),
      contextlib.redirect_stderr(io.StringIO()) as errors,
    ):
      self.assertEqual(script.main(["--catalog", "synthetic.xml"]), 1)
      construct.assert_not_called()
      self.assertIn("Cannot open the plate-loading GUI", errors.getvalue())

  def test_missing_tkinter_fails_before_reader_construction(self) -> None:
    """Scanning requires Tkinter even when the interpreter can run connection checks."""
    with (
      patch.dict(sys.modules, {"tkinter": None}),
      patch.object(script.FilterSlideCatalog, "from_softmax_xml", return_value=catalog_fixture()),
      patch.object(script.FilterMaxF5, "__init__", return_value=None) as construct,
      contextlib.redirect_stderr(io.StringIO()) as errors,
    ):
      self.assertEqual(script.main(["--catalog", "synthetic.xml"]), 1)
      construct.assert_not_called()
      self.assertIn("Tkinter is required", errors.getvalue())

  def test_missing_catalog_fails_before_dialog_or_reader(self) -> None:
    """Never silently fall back to a test catalog."""
    with (
      patch.object(
        script.FilterSlideCatalog, "from_softmax_xml", side_effect=FileNotFoundError("XML")
      ),
      patch.object(script, "LoadingDialog") as dialog,
      patch.object(script.FilterMaxF5, "__init__", return_value=None) as construct,
      contextlib.redirect_stderr(io.StringIO()),
    ):
      self.assertEqual(script.main(["--catalog", "missing.xml"]), 1)
      dialog.assert_not_called()
      construct.assert_not_called()

  def test_gui_buttons_and_close_are_explicit(self) -> None:
    """Only Okay confirms; Cancel and window-close return false."""
    tkinter = MagicMock(TclError=RuntimeError)

    with (
      patch.dict(sys.modules, {"tkinter": tkinter}),
      patch.object(tkinter, "Tk") as root,
      patch.object(tkinter, "Canvas"),
      patch.object(tkinter, "Frame"),
      patch.object(tkinter, "Label"),
      patch.object(tkinter, "Button") as button,
    ):
      dialog = script.LoadingDialog()
      actions = {call.kwargs["text"]: call.kwargs["command"] for call in button.call_args_list}
      self.assertFalse(dialog.confirmed)
      root.return_value.mainloop.side_effect = actions["Okay"]
      self.assertTrue(dialog.show())
      root.return_value.mainloop.side_effect = actions["Cancel"]
      self.assertFalse(dialog.show())
      close_callback = root.return_value.protocol.call_args.args[1]
      root.return_value.mainloop.side_effect = close_callback
      self.assertFalse(dialog.show())
      dialog.close()
      root.return_value.destroy.assert_called_once()


if __name__ == "__main__":
  unittest.main()
