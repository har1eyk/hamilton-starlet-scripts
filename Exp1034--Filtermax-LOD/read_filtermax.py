#!/usr/bin/env python3
"""Optional, explicitly approved F5 A260 acquisition after MANUAL plate transfer.

Uses the installed FilterMaxF5 API. --check-profile is entirely offline.
There is no default port or assumed physical plate geometry.
"""

from __future__ import annotations

import argparse
import asyncio
import hashlib
import json
import math
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

from pylabrobot.molecular_devices.filtermax import FilterMaxF5, FilterSlideCatalog, PlateGeometry

from exp1034 import environment, manifest, require, write_csv, write_json


def load_profile(path: Path) -> tuple[dict, PlateGeometry]:
  """Reject unverified plate profiles before constructing a device or opening I/O."""
  profile = json.loads(path.read_text())
  require(profile["physical_plate_catalog"] == "3635", "This experiment requires Corning 3635.")
  require(
    profile["geometry_verified"] is True and bool(profile["geometry_evidence"].strip()),
    "A measured/verified 3635 reader profile is required; built-in Costar geometry is not evidence.",
  )
  require(bool(profile["expected_serial"].strip()), "Reader serial must be explicit.")
  geometry = PlateGeometry(**profile["plate_geometry"])
  require((geometry.rows, geometry.columns) == (8, 12), "An 8 by 12 plate is required.")
  numbers = [
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
  ]
  require(all(math.isfinite(v) and v > 0 for v in numbers), "Invalid plate dimensions.")
  require(
    geometry.height <= 22 and geometry.well_depth < geometry.height,
    "Invalid F5 plate height/depth.",
  )
  require(
    geometry.left_column_offset + 11 * geometry.column_spacing + geometry.well_size_x / 2
    < geometry.length,
    "Columns exceed plate dimensions.",
  )
  require(
    geometry.top_row_offset + 7 * geometry.row_spacing + geometry.well_size_y / 2 < geometry.width,
    "Rows exceed plate dimensions.",
  )
  return profile, geometry


async def acquire(
  profile: dict,
  geometry: PlateGeometry,
  catalog: FilterSlideCatalog,
  port: str,
  output: Path,
  preparation_record: Path,
) -> None:
  """Check actual filters, prompt for manual loading, then preserve the raw full plate."""
  preparation = json.loads(preparation_record.read_text())
  require(
    preparation["experiment"] == "Exp1034" and preparation["mode"] == "PHYSICAL_PREPARATION",
    "Supply the actual preparation record, not the dry-run metadata.",
  )
  require(
    preparation["a30_measured_concentration"] == 16 and preparation["a30_total_volume"] == 3000,
    "Physical preparation must document the required A30 starting solution.",
  )
  completed = datetime.fromisoformat(preparation["plate_completed_at"].replace("Z", "+00:00"))
  require(completed.tzinfo is not None, "Preparation timestamp must include timezone.")
  require(completed <= datetime.now(timezone.utc), "Plate completion cannot be in the future.")
  output.mkdir(parents=True, exist_ok=False)
  reader = FilterMaxF5("Exp1034_F5", port=port, filter_slide_catalog=catalog)
  metadata = {
    "experiment": "Exp1034",
    "mode": "PHYSICAL_READER_RUN",
    "wavelength": 260,
    "pathlength_correction_applied": False,
    "reference_subtracted": False,
    "environment": environment(),
    "reader_profile": profile,
    "port": port,
    "preparation_record": preparation,
    "status": "connecting",
    "started_at": datetime.now(timezone.utc).isoformat(),
  }
  write_json(output / "reader_metadata.json", metadata)
  try:
    await reader.setup()
    info = await reader.get_instrument_info()
    require(info.serial_number == profile["expected_serial"], "Unexpected reader serial.")
    status = await reader.get_status()
    require(status.ready, f"Reader is not ready: {status.text}")
    slides = await reader.get_filter_slides()
    selected = catalog.resolve(
      kind="excitation", technique="absorbance", wavelength=260, slide_id=slides.excitation_slide_id
    )
    metadata.update(
      {
        "instrument": asdict(info),
        "installed_slides": asdict(slides),
        "selected_filter": asdict(selected),
      }
    )
    write_json(output / "reader_metadata.json", metadata)
    await reader.move_plate_tray_out()
    print(f"Manually load the verified 3635 plate in {geometry.orientation} orientation.")
    print(
      "Confirm A1 against the reader profile, remove ordinary lid/seal, and inspect all wells for bubbles."
    )
    require(
      input("Type LOADED to close the tray and read A260: ").strip() == "LOADED",
      "Loading cancelled.",
    )
    metadata["read_requested_at"] = datetime.now(timezone.utc).isoformat()
    metadata["seconds_since_plate_completion"] = (
      datetime.now(timezone.utc) - completed
    ).total_seconds()
    results = await reader.read_absorbance(
      geometry,
      wavelength=260,
      wells=[r.well for r in manifest()],
      slide_id=selected.slide_id,
      filter_position=selected.position,
    )
    # This installed API closes/ejects the tray itself; no obsolete eject_after_read kwarg.
    require(len(results) == 1, "Expected one endpoint reader result.")
    result = results[0]
    require(
      result.wavelengths == (260,) and result.reference_subtracted is False,
      "Reader returned an unexpected optical measurement.",
    )
    rows = []
    for record in manifest():
      value = result.data["ABCDEFGH".index(record.well[0])][int(record.well[1:]) - 1]
      available = value is not None and math.isfinite(value)
      rows.append(
        {
          "well": record.well,
          "a260": value if available else "",
          "valid": "",
          "exclusion_reason": "" if available else "missing/nonfinite reader result",
          "read_id": "endpoint_1",
          "timestamp": datetime.fromtimestamp(result.timestamp, timezone.utc).isoformat(),
        }
      )
    # Preserve nonfinite raw values losslessly as strings in the original payload.
    payload = asdict(result)
    payload["data"] = [
      [v if v is None or math.isfinite(v) else str(v) for v in row] for row in result.data
    ]
    write_json(output / "raw_backend_result.json", payload)
    write_csv(output / "raw_a260.csv", rows)
    write_csv(output / "well_manifest.csv", [asdict(r) for r in manifest()])
    metadata.update(
      {
        "status": "complete; validity review required",
        "finished_at": datetime.now(timezone.utc).isoformat(),
        "temperature": result.temperature,
      }
    )
  except BaseException as error:
    metadata.update({"status": "failed_or_cancelled", "error": repr(error)})
    raise
  finally:
    write_json(output / "reader_metadata.json", metadata)
    await reader.stop()


def main() -> None:
  """Require both a run flag and fresh typed approval before any reader I/O."""
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("--profile", type=Path, required=True)
  parser.add_argument("--catalog", type=Path, required=True)
  modes = parser.add_mutually_exclusive_group(required=True)
  modes.add_argument("--check-profile", action="store_true")
  modes.add_argument("--run-reader", action="store_true")
  parser.add_argument("--port")
  parser.add_argument("--preparation-record", type=Path)
  parser.add_argument("--output", type=Path)
  args = parser.parse_args()
  profile, geometry = load_profile(args.profile)
  catalog = FilterSlideCatalog.from_softmax_xml(args.catalog)
  require(
    any(f.wavelength == 260 and "absorbance" in f.techniques for f in catalog.filters),
    "Catalog contains no absorbance 260 nm filter; installed slide is checked at runtime.",
  )
  profile["catalog_sha256"] = hashlib.sha256(args.catalog.read_bytes()).hexdigest()
  if args.check_profile:
    print(
      "Offline profile syntax/catalog checks passed; physical setup and installed filters remain unverified."
    )
    return
  require(
    args.port is not None and args.output is not None and args.preparation_record is not None,
    "--run-reader requires --port, --output and --preparation-record.",
  )
  require(
    input("Type APPROVE F5 RUN to authorize reader I/O and tray motion for this run: ").strip()
    == "APPROVE F5 RUN",
    "Run not approved; no reader connection opened.",
  )
  asyncio.run(acquire(profile, geometry, catalog, args.port, args.output, args.preparation_record))


if __name__ == "__main__":
  main()
