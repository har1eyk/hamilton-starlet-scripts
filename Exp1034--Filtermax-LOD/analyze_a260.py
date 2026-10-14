#!/usr/bin/env python3
"""Analyze one Exp1034 reader repeat with explicitly reviewed calibration ranges."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from datetime import datetime, timezone
from pathlib import Path

from exp1034 import VOLUMES, manifest, require, write_csv, write_json


def load_readings(path: Path) -> dict[str, dict]:
  """Require all 96 original rows, one reader repeat, and explicit validity labels."""
  with path.open(newline="") as stream:
    rows = list(csv.DictReader(stream))
  require(
    len(rows) == 96 and {r["well"] for r in rows} == {r.well for r in manifest()},
    "Raw input must contain each of the 96 wells exactly once.",
  )
  require(
    len({r["read_id"] for r in rows}) == 1 and all(r["read_id"].strip() for r in rows),
    "Analyze one reader repeat at a time; read_id is required.",
  )
  result = {}
  for row in rows:
    require(row["valid"].lower() in ("true", "false"), "Each well needs valid=true or false.")
    valid = row["valid"].lower() == "true"
    require(
      valid or bool(row["exclusion_reason"].strip()), "Invalid points need an exclusion reason."
    )
    require(
      not valid or not row["exclusion_reason"].strip(),
      "A valid point cannot have an exclusion reason.",
    )
    timestamp = datetime.fromisoformat(row["timestamp"].replace("Z", "+00:00"))
    require(timestamp.tzinfo is not None, "Read timestamps must include a timezone.")
    value = float(row["a260"]) if row["a260"].strip() else None
    require(
      not valid or (value is not None and math.isfinite(value)), "Valid readings must be finite."
    )
    require(
      value is None or math.isfinite(value),
      "Use an empty value, not NaN, for unavailable readings.",
    )
    result[row["well"]] = {**row, "a260": value, "valid": valid}
  return result


def fit_line(xs: list[float], ys: list[float]) -> dict:
  """Fit ordinary least squares with an intercept, never a forced zero."""
  require(
    len(set(xs)) >= 3,
    "At least three distinct valid concentrations are required for a screening fit.",
  )
  xmean, ymean = statistics.mean(xs), statistics.mean(ys)
  sxx = sum((x - xmean) ** 2 for x in xs)
  slope = sum((x - xmean) * (y - ymean) for x, y in zip(xs, ys)) / sxx
  require(math.isfinite(slope) and slope > 0, "Calibration slope must be finite and positive.")
  intercept = ymean - slope * xmean
  residuals = [y - intercept - slope * x for x, y in zip(xs, ys)]
  sse = sum(r**2 for r in residuals)
  syy = sum((y - ymean) ** 2 for y in ys)
  return {
    "slope_A260_per_ng_uL": slope,
    "intercept_A260": intercept,
    "r_squared": 1 - sse / syy if syy else None,
    "residual_sd": math.sqrt(sse / (len(xs) - 2)),
    "n": len(xs),
  }


def analyze(readings: dict[str, dict], ranges: dict) -> tuple[list[dict], list[dict]]:
  """Compute separate blank statistics, calibration and preliminary LOD per volume."""
  require(set(ranges) == {str(v) for v in VOLUMES}, "Provide reviewed ranges for all four volumes.")
  results, corrected = [], []
  for volume in VOLUMES:
    selection = ranges[str(volume)]
    selected = selection["concentrations"]
    require(
      selection["reviewed"] is True and bool(selection["reason"].strip()),
      "Inspect the measured range and document its selection before fitting.",
    )
    require(
      len(selected) == len(set(selected)) and len(selected) >= 3,
      "Select at least three distinct calibration concentrations.",
    )
    require(set(selected) <= {0.5, 1, 2, 4, 8, 16}, "Unknown calibration concentration.")
    group = [r for r in manifest() if r.final_volume == volume]
    blanks = [r for r in group if r.final_concentration == 0]
    require(
      all(readings[r.well]["valid"] for r in blanks),
      f"Volume {volume} requires six valid blanks; repeat preparation if blanks are unusable.",
    )
    blank_values = [readings[r.well]["a260"] for r in blanks]
    mean, sd = statistics.mean(blank_values), statistics.stdev(blank_values)
    xs, ys = [], []
    for record in group:
      raw = readings[record.well]
      if record.final_concentration == 0:
        continue
      value = raw["a260"] - mean if raw["a260"] is not None else None
      included = raw["valid"] and record.final_concentration in selected
      if included:
        xs.append(record.final_concentration)
        ys.append(value)
      corrected.append(
        {
          "well": record.well,
          "volume": volume,
          "concentration": record.final_concentration,
          "replicate": record.replicate,
          "replicate_type": record.replicate_type,
          "raw_a260": raw["a260"],
          "blank_mean": mean,
          "blank_corrected_a260": value,
          "used_in_fit": included,
          "valid": raw["valid"],
          "exclusion_reason": raw["exclusion_reason"]
          or ("outside selected calibration range" if not included else ""),
        }
      )
    require(set(xs) == set(selected), "A selected concentration has no valid readings.")
    fitted = fit_line(xs, ys)
    lod = 3.3 * sd / fitted["slope_A260_per_ng_uL"]
    warnings = []
    if sd == 0:
      warnings.append(
        "Blank SD is zero; zero calculated LOD is not evidence of zero detection limit."
      )
    if lod < min(xs) or lod > max(xs):
      warnings.append(
        "Calculated LOD lies outside the fitted concentration range; confirm with fresh samples."
      )
    if len(xs) != 3 * len(selected):
      warnings.append("Invalid sample points leave an unbalanced calibration; inspect exclusions.")
    results.append(
      {
        "volume": volume,
        "blank_wells": [r.well for r in blanks],
        "blank_values": blank_values,
        "blank_mean": mean,
        "blank_sample_sd": sd,
        "calibration_selection": selection,
        **fitted,
        "preliminary_concentration_lod_ng_uL": lod,
        "preliminary_mass_lod_ng": lod * volume,
        "warnings": warnings,
      }
    )
  return results, corrected


def main() -> None:
  """Write analysis separately from untouched reader measurements."""
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("--input", type=Path, required=True)
  parser.add_argument("--ranges", type=Path, required=True)
  parser.add_argument("--reader-metadata", type=Path, required=True)
  parser.add_argument("--output", type=Path, required=True)
  args = parser.parse_args()
  readings = load_readings(args.input)
  ranges = json.loads(args.ranges.read_text())
  metadata = json.loads(args.reader_metadata.read_text())
  require(metadata["wavelength"] == 260, "Reader metadata must identify A260.")
  require(
    metadata["pathlength_correction_applied"] is False,
    "This screening analysis expects raw, uncorrected A260.",
  )
  require(
    metadata["reference_subtracted"] is False, "Reference-subtracted values are not raw A260."
  )
  results, corrected = analyze(readings, ranges)
  args.output.mkdir(parents=True, exist_ok=False)
  # Keep the exact original bytes, including invalid readings and comments.
  (args.output / "raw_input.csv").write_bytes(args.input.read_bytes())
  write_csv(args.output / "blank_corrected_samples.csv", corrected)
  write_json(
    args.output / "analysis.json",
    {
      "experiment": "Exp1034",
      "classification": "screening; not an established LOD or LOQ",
      "generated_at": datetime.now(timezone.utc).isoformat(),
      "raw_input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
      "reader_metadata": metadata,
      "results": results,
      "interpretation": "Replicates are dispenses from one chain. Reader repeats are not independent "
      "preparations. Confirm with fresh blanks and independently prepared low-concentration samples.",
    },
  )
  print(f"Wrote screening analysis to {args.output.resolve()}")


if __name__ == "__main__":
  main()
