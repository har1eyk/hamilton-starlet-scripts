#!/usr/bin/env python3
"""Read all 96 Exp1034 DNA wells at 260 nm on the FilterMax F5.

Default: offline preview only. --run-reader requests an interactively approved read.
Save raw A260 and separate volume-matched blank corrections, without converting OD
to DNA concentration. Keep this script beside well_assignments.json.
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
from statistics import mean, stdev
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
DEFAULT_PORT = (
  "/dev/serial/by-id/usb-Prolific_Technology_Inc._USB-Serial_Controller_BFADb2A8J14-if00-port0"
)
EXPECTED_SERIAL = "1191"
WAVELENGTH = 260
ROWS = "ABCDEFGH"
WELLS = tuple(f"{row}{column}" for row in ROWS for column in range(1, 13))
VOLUMES = (50, 75, 100, 200)
CONCENTRATIONS = {"A31": 0, "A30": 16, "A29": 8, "A28": 4, "A27": 2, "A26": 1, "A25": 0.5}
APPROVAL = "APPROVE DNA A260 F5 RUN"


def require(condition: bool, message: str) -> None:
  """Reject invalid inputs without relying on assertions."""
  if not condition:
    raise ValueError(message)


def utc_now() -> str:
  """Return a timezone-aware run timestamp."""
  return datetime.now(timezone.utc).isoformat()


def json_safe(value: Any) -> Any:
  """Preserve nonfinite instrument values as strings in strict JSON."""
  if isinstance(value, float) and not math.isfinite(value):
    return str(value)
  if isinstance(value, dict):
    return {key: json_safe(item) for key, item in value.items()}
  if isinstance(value, (list, tuple)):
    return [json_safe(item) for item in value]
  return value


def write_json(path: Path, value: Any) -> None:
  """Write a strict JSON record without losing nonfinite instrument values."""
  path.write_text(json.dumps(json_safe(value), indent=2, allow_nan=False) + "\n")


def write_csv(path: Path, rows: list[dict]) -> None:
  """Write named columns with unavailable values left empty."""
  with path.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)


def load_manifest() -> list[dict]:
  """Preserve Exp1034 assignments and label concentrations as planned, not measured."""
  records = []
  for volume, concentration, source, wells in json.loads(
    (HERE / "well_assignments.json").read_text()
  ):
    require(volume in VOLUMES and source in CONCENTRATIONS, "Unknown Exp1034 volume or source.")
    require(concentration == CONCENTRATIONS[source], "Unexpected planned DNA concentration.")
    for replicate, well in enumerate(wells, 1):
      records.append(
        {
          "well": well,
          "sample_id": f"{source}_{volume:03d}uL_R{replicate}",
          "source_tube": source,
          "planned_volume_uL": volume,
          "planned_dna_concentration_ng_uL": concentration,
          "replicate": replicate,
          "replicate_type": "dispensing",
          "role": "buffer_blank" if source == "A31" else "dna_sample",
        }
      )
  require(
    len(records) == 96 and {record["well"] for record in records} == set(WELLS),
    "The plate map must contain A1-H12 exactly once.",
  )
  require(len({record["sample_id"] for record in records}) == 96, "Duplicate sample IDs.")
  counts = Counter((record["planned_volume_uL"], record["source_tube"]) for record in records)
  require(
    all(
      counts[volume, source] == (6 if source == "A31" else 3)
      for volume in VOLUMES
      for source in CONCENTRATIONS
    ),
    "Expected six A31 blanks and three wells per DNA source at each volume.",
  )
  return sorted(records, key=lambda record: WELLS.index(record["well"]))


def load_geometry(path: Path | None) -> PlateGeometry:
  """Reuse the Orange G F5 portrait profile or an explicit PLR geometry JSON."""
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
  result: AbsorbanceResult, records: list[dict]
) -> tuple[list[dict], list[dict], list[dict]]:
  """Export untouched A260 and separate corrections using six blanks at each volume."""
  require(
    result.wavelengths == (WAVELENGTH,) and not result.reference_subtracted,
    "Expected raw 260 nm absorbance without reference subtraction.",
  )
  require(
    len(result.data) == 8 and all(len(row) == 12 for row in result.data),
    "Expected an 8 by 12 result matrix.",
  )
  timestamp = datetime.fromtimestamp(result.timestamp, timezone.utc).isoformat()
  raw = []
  for record in records:
    well = record["well"]
    value = result.data[ROWS.index(well[0])][int(well[1:]) - 1]
    finite = value is not None and math.isfinite(value)
    raw.append(
      {
        **record,
        "wavelength_nm": WAVELENGTH,
        "a260": value if finite else None,
        "valid": "",
        "exclusion_reason": "",
        "read_id": "endpoint_1",
        "timestamp": timestamp,
        "data_status": "finite" if finite else "missing_or_nonfinite",
      }
    )
  blanks = []
  for volume in VOLUMES:
    group = [r for r in raw if r["source_tube"] == "A31" and r["planned_volume_uL"] == volume]
    values = [r["a260"] for r in group if r["a260"] is not None]
    complete = len(group) == len(values) == 6
    blanks.append(
      {
        "planned_volume_uL": volume,
        "blank_wells": ";".join(r["well"] for r in group),
        "finite_blank_count": len(values),
        "blank_mean_a260": mean(values) if complete else None,
        "blank_sample_sd_a260": stdev(values) if complete else None,
        "correction_status": "complete" if complete else "incomplete_blank_group",
      }
    )
  backgrounds = {b["planned_volume_uL"]: b["blank_mean_a260"] for b in blanks}
  corrected = []
  for record, reading in zip(records, raw):
    background = backgrounds[record["planned_volume_uL"]]
    value = reading["a260"]
    corrected.append(
      {
        **record,
        "read_id": reading["read_id"],
        "timestamp": timestamp,
        "raw_a260": value,
        "blank_mean_a260": background,
        "blank_corrected_a260": value - background
        if value is not None and background is not None
        else None,
        "data_status": reading["data_status"],
        "correction_status": "incomplete_blank_group"
        if background is None
        else reading["data_status"],
      }
    )
  return raw, corrected, blanks


def build_plan(
  args: argparse.Namespace,
  catalog: FilterSlideCatalog,
  geometry: PlateGeometry,
  records: list[dict],
) -> dict:
  """Describe the fixed A260 acquisition and raw-data provenance without opening hardware."""
  selected = catalog.resolve(kind="excitation", technique="absorbance", wavelength=WAVELENGTH)
  return {
    "experiment": "Exp1034",
    "sample": "DNA dilution series and matching buffer blanks",
    "physical_plate_catalog": "3635",
    "plate_uv_transparent": True,
    "wavelength": WAVELENGTH,
    "catalog_filter": asdict(selected),
    "plate_geometry": asdict(geometry),
    "geometry_source": str(args.geometry.resolve())
    if args.geometry
    else "PLR Costar 96 clear portrait profile, reused from the successful Orange G read",
    "geometry_calibration": "Not independently calibrated for Corning 3635 by this script",
    "well_count": len(records),
    "port": args.port,
    "expected_serial": args.expected_serial,
    "blank_correction": "Separate CSV: subtract six A31 blanks' mean at matching planned volume",
    "raw_csv": "raw_a260.csv",
    "raw_csv_blank_subtracted": False,
    "pathlength_correction_applied": False,
    "reference_subtracted": False,
    "catalog_path": str(args.catalog.resolve()),
    "catalog_sha256": hashlib.sha256(args.catalog.read_bytes()).hexdigest(),
    "assignment_sha256": hashlib.sha256((HERE / "well_assignments.json").read_bytes()).hexdigest(),
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "pylabrobot_file": pylabrobot.__file__,
    "python_executable": sys.executable,
    "output": str(args.output),
  }


async def acquire(
  args: argparse.Namespace,
  catalog: FilterSlideCatalog,
  geometry: PlateGeometry,
  records: list[dict],
  plan: dict,
  confirm: Callable[[str], str] = input,
) -> None:
  """Read once after approval, saving raw data before validation and preserving failures."""
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
    "mode": "PHYSICAL_DNA_A260_READ",
    "status": "connecting",
    "started_at": utc_now(),
  }
  reader = None
  primary_error: BaseException | None = None
  try:
    write_json(args.output / "reader_metadata.json", metadata)
    write_csv(args.output / "well_manifest.csv", records)
    reader = FilterMaxF5("Exp1034_DNA_A260", args.port, filter_slide_catalog=catalog)
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
      wavelength=WAVELENGTH,
      slide_id=slides.excitation_slide_id,
    )
    metadata.update(selected_filter=asdict(selected), status="awaiting_plate")
    write_json(args.output / "reader_metadata.json", metadata)
    print(
      f"F5 serial {info.serial_number}: {WAVELENGTH} nm, {selected.bandwidth} nm bandwidth, "
      f"slide {selected.slide_id}, slot {selected.position}."
    )
    print(f"Corning 3635 UV-transparent plate; {geometry.orientation} orientation.")
    print(f"Reader profile: {geometry.name}; absorbance height: {geometry.absorbance_z:g} mm.")
    await reader.move_plate_tray_out()
    corner = "upper-right" if geometry.orientation == "portrait" else "upper-left"
    print(f"Load with A1 at the {corner}, nearest the instrument, viewed from the front.")
    print("Remove the lid/seal, inspect for bubbles, seat the plate, and clear the tray path.")
    print("A31 blank wells must contain DNA-free matching buffer at their planned volumes.")
    require(
      confirm(
        "Type LOADED to confirm the plate/profile and start the full-plate A260 read: "
      ).strip()
      == "LOADED",
      "Loading cancelled; no measurement requested.",
    )
    metadata.update(status="reading", read_requested_at=utc_now())
    write_json(args.output / "reader_metadata.json", metadata)
    results = await reader.read_absorbance(
      geometry,
      wavelength=WAVELENGTH,
      wells=WELLS,
      slide_id=selected.slide_id,
      filter_position=selected.position,
    )
    # The installed driver loads and ejects the tray for this endpoint operation.
    write_json(args.output / "raw_backend_results.json", [asdict(result) for result in results])
    require(len(results) == 1, "Expected one endpoint result.")
    raw, corrected, blanks = result_rows(results[0], records)
    write_csv(args.output / "raw_a260.csv", raw)
    write_csv(args.output / "blank_corrected_a260.csv", corrected)
    write_csv(args.output / "blank_summary.csv", blanks)
    missing = sum(row["a260"] is None for row in raw)
    metadata.update(
      status="complete" if missing == 0 else "complete_with_missing_values",
      validity_review="pending",
      missing_well_count=missing,
      incomplete_blank_volumes=[
        b["planned_volume_uL"] for b in blanks if b["blank_mean_a260"] is None
      ],
      temperature=results[0].temperature,
      measurement_timestamp=results[0].timestamp,
    )
    print(f"Saved {96 - missing}/96 finite readings to {args.output / 'raw_a260.csv'}")
    print("Volume-matched blank corrections are saved separately; negative values are retained.")
    print("The reader has ejected the plate. Raw A260 is not pathlength-corrected.")
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
  """Preview offline by default; physical acquisition requires an explicit run request."""
  parser = argparse.ArgumentParser(description=__doc__)
  mode = parser.add_mutually_exclusive_group()
  mode.add_argument("--dry-run", action="store_true", help="Display the offline plan (default)")
  mode.add_argument("--run-reader", action="store_true", help="Request an approved physical read")
  parser.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
  parser.add_argument("--port", default=DEFAULT_PORT, help="F5 serial device; no adapter guessing")
  parser.add_argument("--expected-serial", default=EXPECTED_SERIAL)
  parser.add_argument("--geometry", type=Path, help="JSON object containing PlateGeometry fields")
  parser.add_argument("--output", type=Path, help="Fresh result directory; auto-named under runs/")
  args = parser.parse_args(argv)
  try:
    require(bool(args.expected_serial.strip()), "Expected reader serial cannot be empty.")
    require(bool(args.port.strip()), "Reader port cannot be empty.")
    records = load_manifest()
    geometry = load_geometry(args.geometry)
    catalog = FilterSlideCatalog.from_softmax_xml(args.catalog)
    if args.output is None:
      stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
      args.output = HERE / "runs" / f"dna-a260-f5-{stamp}-{uuid4().hex[:8]}"
    args.output = args.output.expanduser().resolve()
    plan = build_plan(args, catalog, geometry, records)
    if not args.run_reader:
      print(json.dumps(plan, indent=2))
      print("OFFLINE ONLY: no reader constructed, port opened, tray moved, or data generated.")
      return 0
    if args.port.startswith("/dev/"):
      require(
        Path(args.port).exists(), f"Serial device not found: {args.port}. Set --port explicitly."
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
