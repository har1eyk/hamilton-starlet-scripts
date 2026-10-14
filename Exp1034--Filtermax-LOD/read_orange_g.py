#!/usr/bin/env python3
"""Read the Exp1034 Orange G test plate with PLR's FilterMax F5 driver.

Default: offline plan only. Use --run-reader for an interactively approved run.
Reads all 96 wells once at 450 nm; no DNA concentration or DNA LOD is calculated.
Keep this script beside well_assignments.json. See README.md for loading and use.
"""

from __future__ import annotations

import argparse
import asyncio
import csv
import hashlib
import json
import math
import sys
from collections import Counter
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean
from typing import Any, Callable, Sequence
from uuid import uuid4

import pylabrobot
from pylabrobot.molecular_devices.filtermax import (
  AbsorbanceResult,
  FilterMaxF5,
  FilterSlideCatalog,
  PlateGeometry,
)

HERE = Path(__file__).resolve().parent
DEFAULT_CATALOG = HERE.parent / "helper-scripts" / "filtermax_slides.xml"
DEFAULT_PORT = "/dev/serial/by-id/usb-FTDI_US232R_FT4QLERM-if00-port0"
EXPECTED_SERIAL = "1191"
ROWS = "ABCDEFGH"
WELLS = tuple(f"{row}{column}" for row in ROWS for column in range(1, 13))
VOLUMES = (50, 75, 100, 200)
SOURCES = ("A31", "A30", "A29", "A28", "A27", "A26", "A25")
APPROVAL = "APPROVE ORANGE G F5 RUN"


def require(condition: bool, message: str) -> None:
  """Reject invalid inputs without relying on assertions."""
  if not condition:
    raise ValueError(message)


def utc_now() -> str:
  """Return a timezone-aware run timestamp."""
  return datetime.now(timezone.utc).isoformat()


def json_safe(value: Any) -> Any:
  """Preserve nonfinite instrument numbers as strings in otherwise strict JSON."""
  if isinstance(value, float) and not math.isfinite(value):
    return str(value)
  if isinstance(value, dict):
    return {key: json_safe(item) for key, item in value.items()}
  if isinstance(value, (list, tuple)):
    return [json_safe(item) for item in value]
  return value


def write_json(path: Path, value: Any) -> None:
  """Save a readable strict JSON record."""
  path.write_text(json.dumps(json_safe(value), indent=2, allow_nan=False) + "\n")


def write_csv(path: Path, rows: list[dict]) -> None:
  """Save named columns, leaving unavailable values empty."""
  with path.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)


def load_manifest() -> list[dict]:
  """Reuse well, tube, volume and replicate identities without DNA concentration labels."""
  records = []
  for volume, _dna_concentration, source, wells in json.loads(
    (HERE / "well_assignments.json").read_text()
  ):
    require(volume in VOLUMES and source in SOURCES, "Unknown Exp1034 volume or source.")
    for replicate, well in enumerate(wells, 1):
      records.append(
        {
          "well": well,
          "sample_id": f"{source}_{volume:03d}uL_R{replicate}",
          "source_tube": source,
          "planned_volume_uL": volume,
          "replicate": replicate,
          "role": "planned_buffer_blank" if source == "A31" else "orange_g_test",
        }
      )
  require(
    len(records) == 96 and {record["well"] for record in records} == set(WELLS),
    "The plate map must contain A1-H12 exactly once.",
  )
  counts = Counter((record["planned_volume_uL"], record["source_tube"]) for record in records)
  require(
    all(
      counts[volume, source] == (6 if source == "A31" else 3)
      for volume in VOLUMES
      for source in SOURCES
    ),
    "Expected six A31 blanks and three wells per dye source at each volume.",
  )
  return sorted(records, key=lambda record: WELLS.index(record["well"]))


def load_geometry(path: Path | None) -> PlateGeometry:
  """Use the existing F5 Costar portrait preset or an explicit PLR geometry JSON."""
  geometry = (
    PlateGeometry.costar_96_clear_portrait()
    if path is None
    else PlateGeometry(**json.loads(path.read_text()))
  )
  require((geometry.rows, geometry.columns) == (8, 12), "An 8 by 12 plate is required.")
  dimensions = (
    geometry.length,
    geometry.width,
    geometry.height,
    geometry.well_depth,
    geometry.left_column_offset,
    geometry.top_row_offset,
    geometry.column_spacing,
    geometry.row_spacing,
    geometry.well_size_x,
    geometry.well_size_y,
    geometry.absorbance_z,
  )
  require(all(math.isfinite(value) and value > 0 for value in dimensions), "Invalid geometry.")
  require(geometry.well_depth < geometry.height <= 22, "Invalid F5 plate height/depth.")
  require(geometry.absorbance_z <= geometry.height, "Read height exceeds plate height.")
  require(
    geometry.well_size_x <= geometry.column_spacing
    and geometry.well_size_y <= geometry.row_spacing
    and geometry.left_column_offset >= geometry.well_size_x / 2
    and geometry.top_row_offset >= geometry.well_size_y / 2
    and geometry.left_column_offset + 11 * geometry.column_spacing + geometry.well_size_x / 2
    <= geometry.length
    and geometry.top_row_offset + 7 * geometry.row_spacing + geometry.well_size_y / 2
    <= geometry.width,
    "Well positions do not fit within the plate.",
  )
  return geometry


def result_rows(
  result: AbsorbanceResult, records: list[dict], wavelength: int, blank_correct: bool
) -> tuple[list[dict], list[dict]]:
  """Map row-major OD to wells and optionally subtract six blanks at matching volume."""
  require(
    result.wavelengths == (wavelength,) and not result.reference_subtracted,
    "Unexpected wavelength or reference subtraction in the returned result.",
  )
  require(
    len(result.data) == 8 and all(len(row) == 12 for row in result.data),
    "Expected an 8 by 12 result matrix.",
  )
  rows = []
  for record in records:
    well = record["well"]
    value = result.data[ROWS.index(well[0])][int(well[1:]) - 1]
    finite = value is not None and math.isfinite(value)
    rows.append(
      {
        **record,
        "wavelength_nm": wavelength,
        "raw_od": value if finite else None,
        "data_status": "finite" if finite else "missing_or_nonfinite",
      }
    )
  blanks = []
  if blank_correct:
    for volume in VOLUMES:
      group = [
        row for row in rows if row["source_tube"] == "A31" and row["planned_volume_uL"] == volume
      ]
      values = [row["raw_od"] for row in group if row["raw_od"] is not None]
      background = mean(values) if len(values) == 6 else None
      blanks.append(
        {
          "planned_volume_uL": volume,
          "blank_wells": ";".join(row["well"] for row in group),
          "finite_blank_count": len(values),
          "blank_mean_od": background,
        }
      )
      for row in rows:
        if row["planned_volume_uL"] != volume:
          continue
        row["blank_mean_od"] = background
        row["blank_corrected_od"] = (
          row["raw_od"] - background
          if row["raw_od"] is not None and background is not None
          else None
        )
        row["correction_status"] = (
          "incomplete_blank_group" if background is None else row["data_status"]
        )
  return rows, blanks


async def acquire(
  args: argparse.Namespace,
  catalog: FilterSlideCatalog,
  geometry: PlateGeometry,
  records: list[dict],
  plan: dict,
  confirm: Callable[[str], str] = input,
) -> None:
  """Require this run's approval before constructing the reader; persist failures and raw data."""
  require(not args.output.exists(), "Output directory already exists; choose a fresh directory.")
  print("Close SoftMax Pro and any other connection to the F5. Transfer the plate manually.")
  require(
    confirm(f"Type {APPROVAL} to authorize reader I/O and tray motion for this run: ").strip()
    == APPROVAL,
    "Run not approved; no hardware connection opened.",
  )
  args.output.mkdir(parents=True, exist_ok=False)
  metadata = {
    **plan,
    "mode": "PHYSICAL_ORANGE_G_READ",
    "status": "connecting",
    "started_at": utc_now(),
  }
  reader = None
  primary_error: BaseException | None = None
  try:
    write_json(args.output / "reader_metadata.json", metadata)
    write_csv(args.output / "well_manifest.csv", records)
    reader = FilterMaxF5("Exp1034_OrangeG", args.port, filter_slide_catalog=catalog)
    await reader.setup()
    info = await reader.get_instrument_info()
    metadata["instrument"] = asdict(info)
    require(info.serial_number == args.expected_serial, "Unexpected F5 serial number.")
    status = await reader.get_status()
    require(status.ready, f"Reader is not ready: {status.text}")
    slides = await reader.get_filter_slides()
    metadata["installed_slides"] = asdict(slides)
    selected = catalog.resolve(
      kind="excitation",
      technique="absorbance",
      wavelength=args.wavelength,
      slide_id=slides.excitation_slide_id,
    )
    metadata["selected_filter"] = asdict(selected)
    metadata["status"] = "awaiting_plate"
    write_json(args.output / "reader_metadata.json", metadata)
    print(
      f"F5 serial {info.serial_number}: {args.wavelength} nm, "
      f"{selected.bandwidth} nm bandwidth, slide {selected.slide_id}, slot {selected.position}."
    )
    print(f"Geometry: {geometry.name}; orientation: {geometry.orientation}.")
    print("Confirm that this profile matches the physical Corning 3635 plate before loading.")
    print("Default Costar geometry is a starting profile, not a 3635-specific F5 calibration.")
    await reader.move_plate_tray_out()
    corner = "upper-right" if geometry.orientation == "portrait" else "upper-left"
    print(f"Load with A1 at the {corner}, nearest the instrument, viewed from the front.")
    print("Remove the lid/seal, inspect for bubbles, seat the plate, and clear the tray path.")
    if args.blank_correct:
      print("Blank correction requires all A31 wells to contain dye-free matching buffer.")
    require(
      confirm("Type LOADED to confirm the plate/profile and start the full-plate read: ").strip()
      == "LOADED",
      "Loading cancelled; no measurement requested.",
    )
    metadata.update(status="reading", read_requested_at=utc_now())
    write_json(args.output / "reader_metadata.json", metadata)
    results = await reader.read_absorbance(
      geometry,
      wavelength=args.wavelength,
      wells=WELLS,
      slide_id=selected.slide_id,
      filter_position=selected.position,
    )
    # The installed driver loads and ejects the tray; it has no eject_after_read parameter.
    write_json(args.output / "raw_backend_results.json", [asdict(result) for result in results])
    require(len(results) == 1, "Expected one endpoint result.")
    rows, blanks = result_rows(results[0], records, args.wavelength, args.blank_correct)
    write_csv(args.output / "absorbance.csv", rows)
    if blanks:
      write_csv(args.output / "blank_summary.csv", blanks)
    missing = sum(row["raw_od"] is None for row in rows)
    metadata.update(
      status="complete" if missing == 0 else "complete_with_missing_values",
      missing_well_count=missing,
      temperature=results[0].temperature,
      measurement_timestamp=results[0].timestamp,
    )
    print(f"Saved {96 - missing}/96 finite readings to {args.output / 'absorbance.csv'}")
    print("The reader has ejected the plate. Raw OD is not pathlength-corrected.")
  except BaseException as error:
    primary_error = error
    metadata.update(status="failed_or_cancelled", error=repr(error))
    raise
  finally:
    cleanup_error: BaseException | None = None
    if reader is not None:
      try:
        await reader.stop()
      except BaseException as error:
        cleanup_error = error
        metadata["shutdown_error"] = repr(error)
        print(f"Reader shutdown failed: {error!r}", file=sys.stderr)
    metadata["finished_at"] = utc_now()
    try:
      write_json(args.output / "reader_metadata.json", metadata)
    except Exception as error:
      if primary_error is None and cleanup_error is None:
        raise
      print(f"Final metadata save also failed: {error!r}", file=sys.stderr)
    if cleanup_error is not None and primary_error is None:
      raise cleanup_error


def main(argv: Sequence[str] | None = None) -> int:
  """Print an offline plan by default; acquire only on explicit interactive request."""
  parser = argparse.ArgumentParser(description=__doc__)
  mode = parser.add_mutually_exclusive_group()
  mode.add_argument(
    "--dry-run", action="store_true", help="Validate and display the offline plan (default)"
  )
  mode.add_argument("--run-reader", action="store_true", help="Request an approved physical read")
  parser.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
  parser.add_argument(
    "--port", default=DEFAULT_PORT, help="F5 serial device; use --port if the adapter changed"
  )
  parser.add_argument("--expected-serial", default=EXPECTED_SERIAL)
  parser.add_argument(
    "--wavelength", type=int, default=450, help="Visible absorbance filter, default 450 nm"
  )
  parser.add_argument("--geometry", type=Path, help="JSON object containing PlateGeometry fields")
  parser.add_argument(
    "--blank-correct",
    action="store_true",
    help="Confirm A31 is dye-free buffer; subtract volume-matched blanks",
  )
  parser.add_argument(
    "--output", type=Path, help="Fresh result directory; auto-named under runs/ by default"
  )
  args = parser.parse_args(argv)
  try:
    require(400 <= args.wavelength <= 700, "This Orange G script requires a visible wavelength.")
    require(bool(args.expected_serial.strip()), "Expected reader serial cannot be empty.")
    records = load_manifest()
    geometry = load_geometry(args.geometry)
    catalog = FilterSlideCatalog.from_softmax_xml(args.catalog)
    selected = catalog.resolve(
      kind="excitation", technique="absorbance", wavelength=args.wavelength
    )
    if args.output is None:
      stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
      args.output = HERE / "runs" / f"orange-g-f5-{stamp}-{uuid4().hex[:8]}"
    args.output = args.output.expanduser().resolve()
    plan = {
      "experiment": "Exp1034",
      "sample": "Orange G test plate; no DNA",
      "physical_plate_catalog": "3635",
      "wavelength_nm": args.wavelength,
      "catalog_filter": asdict(selected),
      "plate_geometry": asdict(geometry),
      "geometry_source": str(args.geometry)
      if args.geometry
      else "PLR Costar 96 clear portrait preset",
      "geometry_calibration": "Not independently verified for this physical plate by this script",
      "well_count": len(records),
      "port": args.port,
      "expected_serial": args.expected_serial,
      "blank_correction": "A31 mean at matching planned volume" if args.blank_correct else "none",
      "pathlength_correction": False,
      "reference_subtraction": False,
      "catalog_path": str(args.catalog.resolve()),
      "catalog_sha256": hashlib.sha256(args.catalog.read_bytes()).hexdigest(),
      "assignment_sha256": hashlib.sha256(
        (HERE / "well_assignments.json").read_bytes()
      ).hexdigest(),
      "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      "pylabrobot_file": pylabrobot.__file__,
      "python_executable": sys.executable,
      "output": str(args.output),
    }
    if not args.run_reader:
      print(json.dumps(plan, indent=2))
      print("OFFLINE ONLY: no reader constructed, port opened, tray moved, or data generated.")
      return 0
    if args.port.startswith("/dev/"):
      require(
        Path(args.port).exists(),
        f"Serial device not found: {args.port}. Set --port to the F5 adapter.",
      )
    asyncio.run(acquire(args, catalog, geometry, records, plan))
    return 0
  except (KeyboardInterrupt, EOFError):
    print("Cancelled.", file=sys.stderr)
    return 130
  except Exception as error:
    print(f"Error: {error}", file=sys.stderr)
    return 1


if __name__ == "__main__":
  raise SystemExit(main())
