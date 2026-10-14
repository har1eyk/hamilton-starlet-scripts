#!/usr/bin/env python3
"""Read raw A260 for A1–A12 on a FilterMax F5.

Run with the PyLabRobot checkout on add-filtermax-f5 (including eject_after_read):
  python helper-scripts/filtermax-A260-readings.py --catalog /path/to/filters.xml

Check the connection without a catalog, GUI, tray movement, or measurement:
  python helper-scripts/filtermax-A260-readings.py --check-connection --debug

Select this checkout explicitly when using another interpreter (from the project root):
  PYTHONPATH=/home/hamilton-robot/Documents/Hamilton-Starlet/pylabrobot /usr/bin/python3 \\
    helper-scripts/filtermax-A260-readings.py --check-connection --debug

Setup tries 38400 baud, then 9600 only if the initial handshake is silent.
The connected rate is retained for all queries and measurements in the session.

The XML must describe the instrument's actual filters. Close SoftMax Pro before
running. Use a UV-compatible plate matching the built-in Costar 96-well geometry.
Results are unblanked OD values, without reference or path-length correction.
"""

from __future__ import annotations

import argparse
import asyncio
import inspect
import math
import sys
import traceback
from contextlib import asynccontextmanager
from pathlib import Path
from typing import AsyncIterator, Callable, Mapping, Optional, Sequence

import pylabrobot
from pylabrobot.io.validation_utils import LOG_LEVEL_IO
from pylabrobot.molecular_devices.filtermax import (
  AbsorbanceResult,
  FilterMaxCatalogError,
  FilterMaxF5,
  FilterMaxTimeoutError,
  FilterSelection,
  FilterSlideCatalog,
  InstalledFilterSlides,
  PlateGeometry,
)

DEFAULT_PORT = "/dev/serial/by-id/usb-FTDI_US232R_FT4QLERM-if00-port0"
EXPECTED_SERIAL = "1191"
WAVELENGTHS = (260,)
WELLS = tuple(f"A{column}" for column in range(1, 13))
READY_MESSAGE = "Ready to scan wells for A260 readings in a 96 well plate in portrait orientation."


def format_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
  """Align a small terminal table without an additional dependency."""
  widths = [max(len(row[i]) for row in [headers, *rows]) for i in range(len(headers))]

  def line(row: Sequence[str]) -> str:
    """Pad one table row."""
    return " | ".join(value.ljust(width) for value, width in zip(row, widths))

  return "\n".join([line(headers), "-+-".join("-" * width for width in widths), *map(line, rows)])


def filter_inventory(catalog: FilterSlideCatalog, slides: InstalledFilterSlides) -> str:
  """List catalog definitions and distinguish installed slides from stored sets."""
  installed = {"excitation": slides.excitation_slide_id, "emission": slides.emission_slide_id}
  rows = [
    (
      item.kind,
      str(item.slide_id),
      item.slide_name,
      str(item.position),
      str(item.wavelength),
      str(item.bandwidth),
      ", ".join(sorted(item.techniques)) or "none",
      "yes" if installed[item.kind] == item.slide_id else "no",
    )
    for item in sorted(catalog.filters, key=lambda item: (item.kind, item.slide_id, item.position))
  ]
  return (
    f"Installed slide IDs: excitation={slides.excitation_slide_id}, "
    f"emission={slides.emission_slide_id}\n"
    + format_table(
      ("Kind", "Slide ID", "Set", "Slot", "nm", "Bandwidth nm", "Techniques", "Installed"), rows
    )
  )


def resolve_filters(
  catalog: FilterSlideCatalog, slides: InstalledFilterSlides
) -> dict[int, FilterSelection]:
  """Require an unambiguous A260 absorbance filter on the installed excitation slide."""
  errors = []
  if not catalog.filters_for_slide("excitation", slides.excitation_slide_id):
    errors.append(f"Unknown installed excitation slide ID {slides.excitation_slide_id}.")
  if not catalog.filters_for_slide("emission", slides.emission_slide_id):
    errors.append(f"Unknown installed emission slide ID {slides.emission_slide_id}.")
  selections = {}
  for wavelength in WAVELENGTHS:
    try:
      selections[wavelength] = catalog.resolve(
        kind="excitation",
        technique="absorbance",
        wavelength=wavelength,
        slide_id=slides.excitation_slide_id,
      )
    except FilterMaxCatalogError as error:
      errors.append(f"A{wavelength}: {error}")
  if errors:
    raise FilterMaxCatalogError("Filter checks failed:\n" + "\n".join(errors))
  return selections


def finite_value(value: Optional[float]) -> Optional[float]:
  """Treat missing and non-finite instrument values as unavailable."""
  return value if value is not None and math.isfinite(value) else None


def results_table(results: Mapping[int, AbsorbanceResult]) -> str:
  """Format the twelve row A absorbance values at 260 nm."""
  rows = []
  for column, well in enumerate(WELLS):
    value = finite_value(results[260].data[0][column]) if 260 in results else None
    rows.append([well, f"{value:.6f}" if value is not None else "N/A"])
  return format_table(("Well", "A260 (OD)"), rows)


class LoadingDialog:
  """Show a portrait loading diagram and require an explicit Okay."""

  def __init__(self) -> None:
    """Prepare the GUI before opening the serial port or moving the tray."""
    try:
      import tkinter as tk
    except ImportError as error:
      raise RuntimeError("Tkinter is required for the plate-loading dialog.") from error
    try:
      self.root = tk.Tk()
    except tk.TclError as error:
      raise RuntimeError(f"Cannot open the plate-loading GUI: {error}") from error
    self.root.withdraw()
    self.root.title("FilterMax F5 — load a 96-well plate")
    self.root.resizable(False, False)
    self.confirmed = False
    self.root.protocol("WM_DELETE_WINDOW", lambda: self._finish(False))
    tk.Label(self.root, text=READY_MESSAGE, wraplength=460, padx=15, pady=12).pack()
    canvas = tk.Canvas(self.root, width=460, height=525, background="white", highlightthickness=0)
    canvas.pack(padx=15)
    canvas.create_rectangle(45, 10, 375, 65, fill="#e2e8f0", outline="#475569")
    canvas.create_text(210, 36, text="FilterMax F5 — instrument / back", font=("TkDefaultFont", 11))
    canvas.create_text(210, 84, text="Top view of the open tray")
    canvas.create_rectangle(45, 106, 375, 475, fill="#f1f5f9", outline="#475569", width=2)
    # Manufacturer's portrait convention: A1 upper-right, H1 upper-left.
    # https://moldevkb.blob.core.windows.net/kb01/documents/br/softmaxpro/
    # smp-7.3/softmax-pro-user-guide_5312581a.pdf (Plate Type Settings).
    for index, row in enumerate("HGFEDCBA"):
      x = 76 + index * 38
      canvas.create_text(x, 122, text=row, fill="#0369a1" if row == "A" else "#334155")
      for column in range(1, 13):
        y = 149 + (column - 1) * 27
        canvas.create_oval(
          x - 10,
          y - 10,
          x + 10,
          y + 10,
          fill="#38bdf8" if row == "A" else "white",
          outline="#475569",
        )
        if index == 0:
          canvas.create_text(32, y, text=str(column))
    for label, y in (("A1", 149), ("A12", 446)):
      canvas.create_line(398, y, 356, y, arrow=tk.LAST, fill="#0369a1", width=2)
      canvas.create_text(419, y, text=label, fill="#0369a1")
    canvas.create_text(210, 494, text="You / front of reader")
    tk.Label(
      self.root,
      text="A1 is nearest the instrument, at the upper-right.\n"
      "Row A is the highlighted right-hand column.\n"
      "Load the matching UV-compatible plate and clear the tray path.\n"
      "Okay closes the tray and starts the A260 reading.",
      padx=15,
      pady=10,
    ).pack()
    buttons = tk.Frame(self.root)
    buttons.pack(pady=(0, 15))
    tk.Button(buttons, text="Cancel", width=12, command=lambda: self._finish(False)).pack(
      side=tk.LEFT, padx=8
    )
    tk.Button(buttons, text="Okay", width=12, command=lambda: self._finish(True)).pack(
      side=tk.LEFT, padx=8
    )

  def _finish(self, confirmed: bool) -> None:
    """Record the button choice and leave the modal event loop."""
    self.confirmed = confirmed
    self.root.withdraw()
    self.root.quit()

  def show(self) -> bool:
    """Wait for Okay, Cancel, or window close on the main thread."""
    self.root.deiconify()
    self.root.mainloop()
    return self.confirmed

  def close(self) -> None:
    """Release the hidden window after scanning or cancellation."""
    self.root.destroy()


def report_environment(reader: FilterMaxF5) -> None:
  """Report the interpreter, imported driver, and configured serial transport."""
  io = reader.io
  print(f"Python: {sys.version.split()[0]} ({sys.executable})")
  print(f"PyLabRobot: {pylabrobot.__file__}")
  print(f"F5 driver: {inspect.getfile(FilterMaxF5)}")
  print(f"Serial port: {reader.port} -> {Path(reader.port).resolve()}")
  print(
    f"Initial serial settings: {io.baudrate} baud, {io.bytesize}{io.parity}{io.stopbits}; "
    f"read timeout={io.timeout}s; write timeout={io.write_timeout}s; "
    f"RTS/CTS={io.rtscts}; DSR/DTR={io.dsrdtr}; XON/XOFF={io.xonxoff}",
    flush=True,
  )
  print("Baud detection: 38400, then 9600 only if the initial handshake is silent.", flush=True)


@asynccontextmanager
async def reader_session(reader: FilterMaxF5) -> AsyncIterator[None]:
  """Connect and disconnect while preserving the primary error on cleanup failure."""
  print(f"Connecting to FilterMax F5 on {reader.port}…", flush=True)
  try:
    await reader.setup()
  except FilterMaxTimeoutError as error:
    raise RuntimeError(
      f"Timed out querying FilterMax F5 identity on {reader.port}. " f"Protocol detail: {error}"
    ) from error
  try:
    io = reader.io
    print(
      f"Connected serial settings: {io.baudrate} baud, {io.bytesize}{io.parity}{io.stopbits}",
      flush=True,
    )
    yield
  except BaseException:
    try:
      await reader.stop()
    except BaseException as stop_error:
      print(f"Reader shutdown also failed: {stop_error!r}", file=sys.stderr, flush=True)
    raise
  else:
    await reader.stop()
    print("Serial connection closed.", flush=True)


async def inspect_reader(
  reader: FilterMaxF5, expected_serial: Optional[str] = None
) -> InstalledFilterSlides:
  """Follow the notebook's read-only preflight and require a ready instrument."""
  stage = "instrument identity"
  try:
    info = await reader.get_instrument_info()
    print(
      f"Instrument: {info.model}; serial={info.serial_number}; "
      f"device code={info.device_code}; firmware={info.firmware}",
      flush=True,
    )
    if expected_serial is not None and info.serial_number != expected_serial:
      raise RuntimeError(f"Expected F5 serial {expected_serial}, received {info.serial_number}.")
    stage = "reader status"
    status = await reader.get_status()
    print(f"Reader status: {status.text}", flush=True)
    stage = "filter slides"
    slides = await reader.get_filter_slides()
    print(
      f"Installed slides: excitation={slides.excitation_slide_id}, "
      f"emission={slides.emission_slide_id}",
      flush=True,
    )
    stage = "temperature"
    temperature = await reader.get_temperature()
    print(f"Reader temperature: {temperature:g} °C", flush=True)
  except FilterMaxTimeoutError as error:
    raise RuntimeError(
      f"Timed out querying {stage} on {reader.port}. Protocol detail: {error}"
    ) from error
  if not status.ready:
    raise RuntimeError(f"FilterMax is not ready: {status.text}")
  return slides


async def run_connection_check(reader: FilterMaxF5) -> None:
  """Check the expected F5's identity and status without motion or measurement."""
  async with reader_session(reader):
    await inspect_reader(reader, expected_serial=EXPECTED_SERIAL)
    print(f"Connection verified: FilterMax F5 serial {EXPECTED_SERIAL} is READY.", flush=True)


async def run_scan(
  reader: FilterMaxF5,
  catalog: FilterSlideCatalog,
  confirm_loading: Callable[[], bool],
) -> None:
  """Check the A260 filter, load the plate, and collect one endpoint measurement."""
  results: dict[int, AbsorbanceResult] = {}
  async with reader_session(reader):
    slides = await inspect_reader(reader, expected_serial=EXPECTED_SERIAL)
    print(filter_inventory(catalog, slides), flush=True)
    selections = resolve_filters(catalog, slides)
    plate = PlateGeometry.costar_96_clear_portrait()
    print(READY_MESSAGE, flush=True)
    await reader.move_plate_tray_out()
    if not confirm_loading():
      print("Cancelled: no readings taken; tray left open.", flush=True)
      return
    await reader.move_plate_tray_in()
    try:
      for wavelength in WAVELENGTHS:
        print(f"Reading A1–A12 at {wavelength} nm…", flush=True)
        selection = selections[wavelength]
        measurements = await reader.read_absorbance(
          plate,
          wavelength=wavelength,
          wells=WELLS,
          slide_id=selection.slide_id,
          filter_position=selection.position,
          eject_after_read=True,
        )
        if len(measurements) != 1:
          raise RuntimeError(f"Expected one endpoint result at {wavelength} nm.")
        result = measurements[0]
        if len(result.data) != 8 or any(len(row) != 12 for row in result.data):
          raise RuntimeError(f"Expected an 8 × 12 result at {wavelength} nm.")
        if result.wavelengths != (wavelength,) or result.reference_subtracted:
          raise RuntimeError(f"Unexpected wavelength or reference subtraction at {wavelength} nm.")
        results[wavelength] = result
    except BaseException:
      try:
        print("\nPARTIAL RESULTS — scan incomplete.")
        print(results_table(results), flush=True)
      except Exception as report_error:
        print(f"Reporting partial results also failed: {report_error!r}", file=sys.stderr)
      raise
    print("\nRaw absorbance (unblanked).")
    print(results_table(results), flush=True)


def main(argv: Optional[Sequence[str]] = None) -> int:
  """Parse configuration and run the interactive workflow."""
  parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
  )
  parser.add_argument(
    "--port",
    default=DEFAULT_PORT,
    help="F5 serial port (default: FTDI US232R FT4QLERM stable device path)",
  )
  parser.add_argument(
    "--catalog", type=Path, help="Actual SoftMax XML export (required for scanning)"
  )
  parser.add_argument(
    "--check-connection",
    action="store_true",
    help="Query identity and status only; requires no catalog or GUI and never moves the tray",
  )
  parser.add_argument(
    "--debug", action="store_true", help="Show serial I/O logging and error tracebacks"
  )
  args = parser.parse_args(argv)
  if not args.check_connection and args.catalog is None:
    parser.error("--catalog is required unless --check-connection is used")
  dialog = None
  exit_code = 0
  try:
    if args.debug:
      pylabrobot.verbose(True, level=LOG_LEVEL_IO)
    if args.check_connection:
      reader = FilterMaxF5(name="filtermax-connection-check", port=args.port)
      report_environment(reader)
      asyncio.run(run_connection_check(reader))
    else:
      if "eject_after_read" not in inspect.signature(FilterMaxF5.read_absorbance).parameters:
        raise RuntimeError(
          "Use the add-filtermax-f5 checkout with the eject_after_read driver option."
        )
      catalog = FilterSlideCatalog.from_softmax_xml(args.catalog.expanduser())
      dialog = LoadingDialog()
      reader = FilterMaxF5(name="filtermax-uv", port=args.port, filter_slide_catalog=catalog)
      report_environment(reader)
      asyncio.run(run_scan(reader, catalog, dialog.show))
  except KeyboardInterrupt:
    print("\nInterrupted.", file=sys.stderr)
    exit_code = 130
  except Exception as error:
    print(f"Error: {error}", file=sys.stderr)
    if args.debug:
      traceback.print_exc()
    exit_code = 1
  finally:
    if dialog is not None:
      try:
        dialog.close()
      except Exception as close_error:
        print(f"Dialog cleanup failed: {close_error!r}", file=sys.stderr)
        if exit_code == 0:
          exit_code = 1
  return exit_code


if __name__ == "__main__":
  raise SystemExit(main())
