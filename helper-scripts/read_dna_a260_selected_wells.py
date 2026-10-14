"""Selected-well FilterMax F5 acquisition and fixed-calibration DNA analysis."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable
from uuid import uuid4

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.figure import Figure

import pylabrobot
from pylabrobot.molecular_devices.filtermax import (
  AbsorbanceResult,
  FilterMaxF5,
  FilterSlideCatalog,
  PlateGeometry,
)

HERE = Path(__file__).resolve().parent
DEFAULT_CATALOG = HERE / "filtermax_slides.xml"
DEFAULT_PORT = (
  "/dev/serial/by-id/usb-Prolific_Technology_Inc._USB-Serial_Controller_BFADb2A8J14-if00-port0"
)
EXPECTED_SERIAL = "1191"
WAVELENGTH = 260
ROWS = "ABCDEFGH"
WELLS = tuple(f"{row}{column}" for row in ROWS for column in range(1, 13))
SLOPE = 0.00338
INTERCEPT = 0.0022
DILUTION_FACTOR = 10
PLOT_TITLE = "Concentration vs Absorbance on the FilterMax F5 for 100ul"


def require(condition: bool, message: str) -> None:
  """Reject invalid inputs without relying on assertions."""
  if not condition:
    raise ValueError(message)


def utc_now() -> str:
  """Return a timezone-aware timestamp."""
  return datetime.now(timezone.utc).isoformat()


def json_safe(value: Any) -> Any:
  """Preserve nonfinite backend values as strings in strict JSON."""
  if isinstance(value, float) and not math.isfinite(value):
    return str(value)
  if isinstance(value, dict):
    return {key: json_safe(item) for key, item in value.items()}
  if isinstance(value, (list, tuple)):
    return [json_safe(item) for item in value]
  return value


def write_json(path: Path, value: Any) -> None:
  """Save a strict JSON record, including lossless nonfinite values."""
  path.write_text(json.dumps(json_safe(value), indent=2, allow_nan=False) + "\n")


def parse_wells(text: str) -> tuple[str, ...]:
  """Validate comma-separated well names and preserve their entered order."""
  wells = tuple(well.strip().upper() for well in text.split(","))
  require(all(wells), "Enter one or more wells separated by commas, e.g. A8, A9.")
  invalid = [well for well in wells if well not in WELLS]
  require(not invalid, f"Invalid wells: {', '.join(invalid)}. Use A1 through H12.")
  require(len(set(wells)) == len(wells), "Duplicate wells; list each well once.")
  return wells


def load_geometry(path: Path | None) -> PlateGeometry:
  """Validate the portrait reader profile or an explicitly supplied geometry."""
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


@dataclass(frozen=True)
class ReadPlan:
  """Offline-validated settings and provenance for one selected-well acquisition."""

  wells: tuple[str, ...]
  port: str
  expected_serial: str
  geometry: PlateGeometry
  catalog: FilterSlideCatalog
  output: Path
  metadata: dict


def prepare_read(
  text: str,
  *,
  port: str = DEFAULT_PORT,
  expected_serial: str = EXPECTED_SERIAL,
  catalog_path: Path = DEFAULT_CATALOG,
  geometry_path: Path | None = None,
  results_root: Path = HERE / "runs",
  notebook_path: Path | None = None,
) -> ReadPlan:
  """Build a fresh run plan without constructing a reader or creating files."""
  wells = parse_wells(text)
  require(bool(port.strip()), "Reader port cannot be empty.")
  require(bool(expected_serial.strip()), "Expected reader serial cannot be empty.")
  geometry_path = None if geometry_path is None else Path(geometry_path).expanduser().resolve()
  geometry = load_geometry(geometry_path)
  catalog_path = Path(catalog_path).expanduser().resolve()
  catalog = FilterSlideCatalog.from_softmax_xml(catalog_path)
  selected = catalog.resolve(kind="excitation", technique="absorbance", wavelength=WAVELENGTH)
  stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
  output = (
    Path(results_root).expanduser().resolve() / f"dna-a260-selected-{stamp}-{uuid4().hex[:8]}"
  )
  row_indices = [ROWS.index(well[0]) for well in wells]
  metadata = {
    "selected_wells": list(wells),
    "well_count": len(wells),
    "hardware_rows": list(ROWS[min(row_indices) : max(row_indices) + 1]),
    "selection_behavior": "Hardware scans the inclusive row span; selected wells are retained.",
    "wavelength": WAVELENGTH,
    "port": port,
    "expected_serial": expected_serial,
    "plate_geometry": asdict(geometry),
    "geometry_source": str(geometry_path)
    if geometry_path
    else "PLR Costar 96 clear portrait profile",
    "catalog_filter": asdict(selected),
    "catalog_path": str(catalog_path),
    "catalog_sha256": hashlib.sha256(catalog_path.read_bytes()).hexdigest(),
    "module_path": str(Path(__file__).resolve()),
    "module_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "pylabrobot_file": pylabrobot.__file__,
    "python_executable": sys.executable,
    "output": str(output),
    "raw_csv": "raw_a260.csv",
    "raw_csv_blank_subtracted": False,
    "pathlength_correction_applied": False,
    "reference_subtracted": False,
  }
  if notebook_path is not None:
    notebook_path = Path(notebook_path).expanduser().resolve()
    metadata.update(
      notebook_path=str(notebook_path),
      notebook_sha256=hashlib.sha256(notebook_path.read_bytes()).hexdigest(),
    )
  return ReadPlan(wells, port, expected_serial, geometry, catalog, output, metadata)


def raw_result_rows(result: AbsorbanceResult, wells: tuple[str, ...]) -> list[dict]:
  """Map raw A260 to selected wells and mark missing or nonfinite values."""
  require(
    result.wavelengths == (WAVELENGTH,) and not result.reference_subtracted,
    "Expected raw 260 nm absorbance without reference subtraction.",
  )
  require(
    len(result.data) == 8 and all(len(row) == 12 for row in result.data),
    "Expected an 8 by 12 result matrix.",
  )
  timestamp = datetime.fromtimestamp(result.timestamp, timezone.utc).isoformat()
  records = []
  for well in wells:
    value = result.data[ROWS.index(well[0])][int(well[1:]) - 1]
    finite = value is not None and math.isfinite(value)
    records.append(
      {
        "well": well,
        "wavelength_nm": WAVELENGTH,
        "a260": value if finite else None,
        "data_status": "finite" if finite else "missing_or_nonfinite",
        "read_id": "endpoint_1",
        "timestamp": timestamp,
      }
    )
  return records


async def acquire_selected(plan: ReadPlan, confirm: Callable[[str], str] = input) -> Path:
  """Execute a requested read, requiring LOADED before starting measurement."""
  require(not plan.output.exists(), "Output directory already exists; prepare a fresh run.")
  if plan.port.startswith("/dev/"):
    require(Path(plan.port).exists(), f"Serial device not found: {plan.port}. Set PORT explicitly.")
  print(f"Read {len(plan.wells)} wells at {WAVELENGTH} nm: {', '.join(plan.wells)}")
  print("Hardware scans rows:", ", ".join(plan.metadata["hardware_rows"]))
  print("Results directory:", plan.output)
  print("Close SoftMax Pro and any other connection to the F5. Transfer the plate manually.")
  plan.output.mkdir(parents=True, exist_ok=False)
  metadata = {
    **plan.metadata,
    "mode": "PHYSICAL_DNA_A260_SELECTED_WELLS_READ",
    "status": "connecting",
    "started_at": utc_now(),
  }
  reader = None
  primary_error: BaseException | None = None
  try:
    write_json(plan.output / "reader_metadata.json", metadata)
    reader = FilterMaxF5("DNA_A260_selected", plan.port, filter_slide_catalog=plan.catalog)
    await reader.setup()
    info = await reader.get_instrument_info()
    metadata["instrument"] = asdict(info)
    require(info.serial_number == plan.expected_serial, "Unexpected F5 serial number.")
    status = await reader.get_status()
    require(status.ready, f"Reader is not ready: {status.text}")
    slides = await reader.get_filter_slides()
    metadata["installed_slides"] = asdict(slides)
    selected = plan.catalog.resolve(
      kind="excitation",
      technique="absorbance",
      wavelength=WAVELENGTH,
      slide_id=slides.excitation_slide_id,
    )
    metadata.update(selected_filter=asdict(selected), status="awaiting_plate")
    write_json(plan.output / "reader_metadata.json", metadata)
    print(
      f"F5 serial {info.serial_number}: {WAVELENGTH} nm, {selected.bandwidth} nm bandwidth, "
      f"slide {selected.slide_id}, slot {selected.position}."
    )
    print(f"Reader profile: {plan.geometry.name}; height: {plan.geometry.absorbance_z:g} mm.")
    await reader.move_plate_tray_out()
    corner = "upper-right" if plan.geometry.orientation == "portrait" else "upper-left"
    print(f"Load a matching UV-transparent plate with A1 at the {corner}, viewed from the front.")
    print(
      "A1 is nearest the instrument. Remove the lid/seal, inspect for bubbles, and seat the plate."
    )
    print("Clear the tray path before continuing.")
    require(
      confirm(
        f"Type LOADED to confirm the plate/profile and read {', '.join(plan.wells)} at A260: "
      ).strip()
      == "LOADED",
      "Loading cancelled; no measurement requested.",
    )
    metadata.update(status="reading", read_requested_at=utc_now())
    write_json(plan.output / "reader_metadata.json", metadata)
    results = await reader.read_absorbance(
      plan.geometry,
      wavelength=WAVELENGTH,
      wells=plan.wells,
      slide_id=selected.slide_id,
      filter_position=selected.position,
    )
    # The driver loads and ejects the tray for this endpoint operation.
    write_json(plan.output / "raw_backend_results.json", [asdict(result) for result in results])
    require(len(results) == 1, "Expected one endpoint result.")
    raw = raw_result_rows(results[0], plan.wells)
    with (plan.output / "raw_a260.csv").open("w", newline="") as stream:
      writer = csv.DictWriter(stream, fieldnames=list(raw[0]))
      writer.writeheader()
      writer.writerows(raw)
    missing = sum(row["a260"] is None for row in raw)
    metadata.update(
      status="complete" if missing == 0 else "complete_with_missing_values",
      missing_well_count=missing,
      temperature=results[0].temperature,
      measurement_timestamp=results[0].timestamp,
    )
    print(
      f"Saved {len(raw) - missing}/{len(raw)} finite readings to {plan.output / 'raw_a260.csv'}"
    )
    print("The reader has ejected the plate. Run the analysis cells to calculate concentrations.")
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
      write_json(plan.output / "reader_metadata.json", metadata)
    except Exception as error:
      if primary_error is None and cleanup_error is None:
        raise
      print(f"Final metadata save also failed: {error!r}", file=sys.stderr)
    if cleanup_error is not None and primary_error is None:
      raise cleanup_error
  return plan.output


def analyze_readings(raw: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
  """Infer the lowest finite blank and apply the fixed calibration and 10x dilution."""
  require({"well", "a260"}.issubset(raw.columns), "Raw data must contain well and a260 columns.")
  wells = parse_wells(",".join(raw["well"].astype(str)))
  values = pd.to_numeric(raw["a260"], errors="coerce").astype(float).tolist()
  finite_indices = [index for index, value in enumerate(values) if math.isfinite(value)]
  require(
    len(finite_indices) >= 2, "At least two finite readings are needed: a blank and a sample."
  )
  blank_index = min(finite_indices, key=lambda index: values[index])
  blank = values[blank_index]
  records = []
  for index, (well, value) in enumerate(zip(wells, values)):
    finite = math.isfinite(value)
    role = "blank" if index == blank_index else "sample" if finite else "unavailable"
    corrected = value - blank if finite else None
    concentration = (corrected - INTERCEPT) / SLOPE if role == "sample" else None
    records.append(
      {
        "well": well,
        "role": role,
        "raw_a260": value if finite else None,
        "blank_corrected_a260": corrected,
        "calculated_concentration_ng_uL": concentration,
        "reported_concentration_ng_uL": concentration * DILUTION_FACTOR
        if concentration is not None
        else None,
        "data_status": "finite" if finite else "missing_or_nonfinite",
      }
    )
  metadata = {
    "blank_well": wells[blank_index],
    "blank_a260": blank,
    "blank_selection": "Lowest finite raw A260; first entered well wins ties.",
    "slope": SLOPE,
    "intercept": INTERCEPT,
    "dilution_factor": DILUTION_FACTOR,
    "equation": "blank_corrected_a260 = 0.00338 * calculated_concentration_ng_uL + 0.0022",
    "reported_concentration": "calculated_concentration_ng_uL * 10",
    "additional_volume_correction_applied": False,
    "pathlength_correction_applied": False,
  }
  return pd.DataFrame(records), metadata


def save_analysis(raw_path: Path) -> tuple[pd.DataFrame, dict]:
  """Save calculated concentrations and calibration provenance beside a raw CSV."""
  raw_path = Path(raw_path).expanduser().resolve()
  raw_bytes = raw_path.read_bytes()
  table, metadata = analyze_readings(pd.read_csv(raw_path))
  metadata.update(
    raw_csv=str(raw_path),
    raw_csv_sha256=hashlib.sha256(raw_bytes).hexdigest(),
    analyzed_at=utc_now(),
  )
  table.to_csv(raw_path.parent / "concentrations.csv", index=False)
  write_json(raw_path.parent / "analysis_metadata.json", metadata)
  return table, metadata


def concentration_display_table(table: pd.DataFrame) -> pd.DataFrame:
  """Format the notebook table without rounding the underlying calculations."""
  columns = {
    "well": "Well",
    "role": "Role",
    "raw_a260": "Raw A260",
    "blank_corrected_a260": "Blank-corrected A260",
    "calculated_concentration_ng_uL": "Calculated concentration (ng/ul)",
    "reported_concentration_ng_uL": "Reported concentration ×10 (ng/ul)",
  }
  formatted = table[list(columns)].copy()
  for column, digits in (
    ("raw_a260", 5),
    ("blank_corrected_a260", 5),
    ("calculated_concentration_ng_uL", 2),
    ("reported_concentration_ng_uL", 1),
  ):
    formatted[column] = [
      "—" if pd.isna(value) else f"{value:.{digits}f}" for value in table[column]
    ]
  return formatted.rename(columns=columns).set_index("Well")


def plot_concentrations(table: pd.DataFrame) -> Figure:
  """Plot calibration concentrations with labels showing 10x reported concentrations."""
  samples = table.loc[table["role"] == "sample"].sort_values("calculated_concentration_ng_uL")
  require(not samples.empty, "No finite sample concentrations to plot.")
  x = samples["calculated_concentration_ng_uL"]
  y = samples["blank_corrected_a260"]
  lower, upper = min(0.0, float(x.min())), max(0.0, float(x.max()))
  padding = max(upper - lower, 1.0) * 0.18
  line_x = [lower - padding, upper + padding]
  line_y = [SLOPE * value + INTERCEPT for value in line_x]
  fig, ax = plt.subplots(figsize=(11, 6.5))
  fig.subplots_adjust(left=0.11, right=0.96, top=0.88, bottom=0.19)
  ax.plot(line_x, line_y, color="#526c88", linewidth=1.8, label="Supplied calibration")
  ax.scatter(
    x, y, color="#007f8b", s=65, edgecolor="white", linewidth=0.8, zorder=3, label="Samples"
  )
  for index, row in enumerate(samples.to_dict("records")):
    offset = (-12, 24) if index % 2 == 0 else (16, -34)
    ax.annotate(
      f"{row['well']}\n{row['reported_concentration_ng_uL']:.1f}ng/ul",
      (row["calculated_concentration_ng_uL"], row["blank_corrected_a260"]),
      xytext=offset,
      textcoords="offset points",
      ha="right" if index % 2 == 0 else "left",
      va="center",
      fontsize=10,
      arrowprops={"arrowstyle": "-", "color": "#7b8794", "linewidth": 0.8},
      bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "none", "alpha": 0.9},
    )
  ax.text(
    0.035,
    0.96,
    "y = 0.00338x + 0.0022\nSlope = 0.00338\ny-intercept = 0.0022",
    transform=ax.transAxes,
    va="top",
    fontsize=11,
    bbox={"boxstyle": "round,pad=0.5", "fc": "#f2f5f8", "ec": "#d9e1e8"},
  )
  ax.set_title(PLOT_TITLE, fontsize=14, pad=18)
  ax.set_xlabel("Concentration (ng/ul)", fontsize=11)
  ax.set_ylabel("Blank-corrected A260 (OD)", fontsize=11)
  ax.set_xlim(line_x)
  ax.margins(y=0.18)
  ax.grid(alpha=0.2)
  ax.set_axisbelow(True)
  ax.spines["top"].set_visible(False)
  ax.spines["right"].set_visible(False)
  ax.legend(loc="lower right", frameon=False)
  fig.text(
    0.5,
    0.055,
    "X-position: calibration concentration before dilution correction.\n"
    "Sample labels: reported concentration after the ×10 dilution factor.",
    ha="center",
    va="center",
    fontsize=10,
    color="#4b5563",
  )
  return fig
