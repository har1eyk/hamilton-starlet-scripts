"""Reproducible analysis of the Exp1034 Orange G plate; no instrument access."""

from __future__ import annotations

import csv
import hashlib
import html
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

VOLUMES = [50, 75, 100, 200]
SOURCES = ["A30", "A29", "A28", "A27", "A26", "A25"]
RELATIVE = {source: 2.0 ** -i for i, source in enumerate(SOURCES)}
LABELS = dict(zip(SOURCES, ["1", "1/2", "1/4", "1/8", "1/16", "1/32"]))
BOOTSTRAPS = 10_000
SEED = 1034


def require(condition: bool, message: str) -> None:
  """Reject incomplete or inconsistent inputs before analysis."""
  if not condition:
    raise ValueError(message)


def digest(path: Path) -> str:
  """Hash the exact source bytes, including line endings."""
  return hashlib.sha256(path.read_bytes()).hexdigest()


def read_plate(path: Path) -> tuple[dict, dict, dict]:
  """Read explicitly named raw-data and sample-ID matrices by well coordinates."""
  rows = list(csv.reader(path.read_text(encoding="utf-8-sig").splitlines()))
  require(any("Test run no.: 151" in row for row in rows), "Expected reader run 151.")
  sections: dict[str, dict] = {"raw": {}, "ids": {}}
  section = None
  seen_sections = []
  header_seen = False
  for row in rows:
    if not row:
      continue
    if row[0] in ("1. Raw Data (492)", "2. Sample IDs"):
      section = "raw" if row[0].startswith("1.") else "ids"
      require(section not in seen_sections, "Duplicate measurement section.")
      seen_sections.append(section)
      header_seen = False
      continue
    if section and row[0] == "" and len(row) >= 13:
      require(row[1:13] == [str(i) for i in range(1, 13)], "Unexpected column order.")
      header_seen = True
      continue
    if section and row[0] in list("ABCDEFGH"):
      require(header_seen and len(row) == 13, "Malformed plate matrix row.")
      for column, value in enumerate(row[1:], 1):
        well = f"{row[0]}{column}"
        require(well not in sections[section], f"Duplicate {section} well: {well}")
        sections[section][well] = float(value) if section == "raw" else value.strip()
  expected = {f"{r}{c}" for r in "ABCDEFGH" for c in range(1, 13)}
  require(all(set(values) == expected for values in sections.values()), "Missing plate wells.")
  require(all(math.isfinite(x) for x in sections["raw"].values()), "Nonfinite raw OD.")
  metadata = {
    "run": 151, "wavelength_nm": 492, "export_header": rows[:5],
    "path_length_correction": False, "correction_setting_source": "user confirmation",
    "concentration_basis": "user-confirmed twofold series relative to A30",
  }
  return sections["raw"], sections["ids"], metadata


def load_inputs(folder: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
  """Reconcile the reader export, fixed-width IDs, and original dispensing manifest."""
  raw_path = folder / "OrangeG_dye_test_raw-data_w_sampleIDs.CSV"
  raw, exported, metadata = read_plate(raw_path)
  fixed = {}
  for line in (folder / "Exp1034_OrangeG_Sample_IDs_fixed.txt").read_text().splitlines():
    if not line.strip():
      continue
    well = line[:4].strip()
    require(well not in fixed, f"Duplicate fixed-file well: {well}")
    fixed[well] = (line[4:11].strip(), line[11:].strip())
  require(set(fixed) == set(raw), "Fixed-file well coverage differs from reader export.")
  rows, audit = [], []
  seen = set()
  for volume, _dna_concentration, source, wells in json.loads(
    (folder / "well_assignments.json").read_text()
  ):
    require(source in SOURCES + ["A31"] and volume in VOLUMES, "Unknown group.")
    for replicate, well in enumerate(wells, 1):
      require(well not in seen, f"Duplicate manifest well: {well}")
      seen.add(well)
      expected_id = f"{source}_{volume:03d}uL_R{replicate}"
      content = f"X{'ABCDEFGH'.index(well[0]) * 12 + int(well[1:])}"
      audit.append({
        "well": well, "expected_content": content, "txt_content": fixed[well][0],
        "expected_id": expected_id, "txt_id": fixed[well][1],
        "reader_id": exported[well], "content_matches": content == fixed[well][0],
        "txt_id_matches": expected_id == fixed[well][1],
        "reader_id_matches": expected_id == exported[well],
      })
      rows.append({
        "well": well, "sample_id": expected_id, "source": source,
        "volume_uL": volume, "replicate": replicate, "is_blank": source == "A31",
        "relative_concentration": 0.0 if source == "A31" else RELATIVE[source],
        "raw_od_492": raw[well],
      })
  require(seen == set(raw), "Manifest well coverage differs from reader export.")
  audit_frame = pd.DataFrame(audit)
  for flag in ["content_matches", "txt_id_matches", "reader_id_matches"]:
    require(bool(audit_frame[flag].all()), f"ID audit failed: {flag}.")
  frame = pd.DataFrame(rows)
  require(frame.sample_id.nunique() == 96, "Nonunique sample IDs.")
  require(int(frame.is_blank.sum()) == 24, "Expected 24 blanks.")
  for volume in VOLUMES:
    subset = frame[frame.volume_uL == volume]
    require(len(subset) == 24, f"Expected 24 wells at {volume} uL.")
    for source in SOURCES + ["A31"]:
      group = subset[subset.source == source]
      count = 6 if source == "A31" else 3
      require(set(group.replicate) == set(range(1, count + 1)), "Replicate coverage mismatch.")
  metadata["source_sha256"] = {p.name: digest(p) for p in folder.iterdir() if p.is_file()}
  return frame, audit_frame, metadata


def correct_blanks(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
  """Subtract each volume's six-well blank mean without clipping negative results."""
  blank_rows = []
  for volume in VOLUMES:
    group = frame[frame.is_blank & (frame.volume_uL == volume)]
    values = group.raw_od_492.to_numpy()
    blank_rows.append({
      "volume_uL": volume, "n": len(group), "blank_wells": ";".join(group.well),
      "blank_mean_od": values.mean(), "blank_sd_od": values.std(ddof=1),
      "blank_sem_od": values.std(ddof=1) / np.sqrt(len(values)),
      "blank_min_od": values.min(), "blank_max_od": values.max(),
    })
  blanks = pd.DataFrame(blank_rows)
  corrected = frame.merge(blanks[["volume_uL", "blank_mean_od", "blank_sd_od"]],
                          on="volume_uL", validate="many_to_one")
  corrected["corrected_od"] = corrected.raw_od_492 - corrected.blank_mean_od
  corrected["used_in_fits"] = ~corrected.is_blank
  corrected["signal_flag"] = ""
  below = ~corrected.is_blank & (corrected.corrected_od < 0)
  near = ~corrected.is_blank & (corrected.corrected_od >= 0) & (
    corrected.corrected_od <= 3 * corrected.blank_sd_od
  )
  corrected.loc[below, "signal_flag"] = "below matched blank mean"
  corrected.loc[near, "signal_flag"] = "within 3 blank SD of zero (screening flag)"
  summaries = []
  for source in SOURCES:
    for volume in VOLUMES:
      group = corrected[(corrected.source == source) & (corrected.volume_uL == volume)]
      mean = group.corrected_od.mean()
      sd = group.corrected_od.std(ddof=1)
      worst_index = (group.corrected_od - mean).abs().idxmax()
      summaries.append({
        "source": source, "relative_concentration": RELATIVE[source], "volume_uL": volume,
        "n": len(group), "wells": ";".join(group.well), "raw_mean_od": group.raw_od_492.mean(),
        "blank_mean_od": group.blank_mean_od.iloc[0], "mean_corrected_od": mean,
        "sample_sd_od": sd, "cv_percent": 100 * sd / mean if mean > 0 else np.nan,
        "combined_mean_se_od": np.sqrt(sd ** 2 / 3 + group.blank_sd_od.iloc[0] ** 2 / 6),
        "largest_replicate_deviation_well": group.loc[worst_index, "well"],
        "largest_replicate_deviation_od": group.loc[worst_index, "corrected_od"] - mean,
      })
  return corrected, blanks, pd.DataFrame(summaries)


def bootstrap_samples(frame: pd.DataFrame, count: int = BOOTSTRAPS, seed: int = SEED
                      ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
  """Draw a shared blank mean per volume and sample replicates within each group."""
  rng = np.random.default_rng(seed)
  blank_draws = np.empty((count, len(VOLUMES)))
  sample_draws = np.empty((count, len(SOURCES), len(VOLUMES), 3))
  for vi, volume in enumerate(VOLUMES):
    values = frame[frame.is_blank & (frame.volume_uL == volume)].sort_values(
      "replicate").raw_od_492.to_numpy()
    blank_draws[:, vi] = values[rng.integers(0, len(values), size=(count, len(values)))].mean(1)
  for si, source in enumerate(SOURCES):
    for vi, volume in enumerate(VOLUMES):
      values = frame[(frame.source == source) & (frame.volume_uL == volume)].sort_values(
        "replicate").raw_od_492.to_numpy()
      sample_draws[:, si, vi, :] = values[rng.integers(0, 3, size=(count, 3))]
  corrected = sample_draws - blank_draws[:, None, :, None]
  return corrected, sample_draws, blank_draws


def regression(x: np.ndarray, y: np.ndarray) -> dict:
  """Fit free-intercept OLS; RMSE uses n observations, not n minus two."""
  centered = x - x.mean()
  slope = float(centered @ y / (centered @ centered))
  intercept = float(y.mean() - slope * x.mean())
  residual = y - (intercept + slope * x)
  sse = float(residual @ residual)
  sst = float((y - y.mean()) @ (y - y.mean()))
  return {"slope": slope, "intercept": intercept, "r_squared": 1 - sse / sst,
          "residual_rmse_od": float(np.sqrt(sse / len(x))), "n": len(x)}


def analyze_fits(corrected: pd.DataFrame, draws: np.ndarray
                 ) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
  """Fit both orientations and bootstrap slopes and group-mean departures."""
  fits, departures, detail = [], [], {}
  specs = [("volume", source) for source in SOURCES] + [("concentration", v) for v in VOLUMES]
  for kind, key in specs:
    if kind == "volume":
      si = SOURCES.index(key)
      levels = np.array(VOLUMES, dtype=float)
      groups = [corrected[(corrected.source == key) & (corrected.volume_uL == v)].sort_values(
        "replicate") for v in VOLUMES]
      boot = draws[:, si, :, :]
      label = f"{key} ({LABELS[key]} of A30)"
      units = "OD/uL"
    else:
      vi = VOLUMES.index(key)
      order = list(reversed(SOURCES))
      levels = np.array([RELATIVE[s] for s in order])
      groups = [corrected[(corrected.source == s) & (corrected.volume_uL == key)].sort_values(
        "replicate") for s in order]
      boot = draws[:, ::-1, vi, :]
      label = f"{key} uL"
      units = "OD/relative concentration"
    x = np.repeat(levels, 3)
    y = np.concatenate([g.corrected_od.to_numpy() for g in groups])
    fit = regression(x, y)
    by = boot.reshape(len(draws), -1)
    centered = x - x.mean()
    slopes = by @ centered / (centered @ centered)
    intercepts = by.mean(1) - slopes * x.mean()
    boot_departures = boot.mean(2) - (intercepts[:, None] + slopes[:, None] * levels)
    slope_ci = np.percentile(slopes, [2.5, 97.5])
    intercept_ci = np.percentile(intercepts, [2.5, 97.5])
    departure_ci = np.percentile(boot_departures, [2.5, 97.5], axis=0)
    departing = (departure_ci[0] > 0) | (departure_ci[1] < 0)
    weak = slope_ci[0] <= 0 <= slope_ci[1]
    verdict = ("response too weak to assess confidently" if weak else
               "exploratory evidence of nonlinearity" if departing.any() else
               "no clear departure detected")
    fit_id = f"{kind}_{key}"
    fits.append({"fit_id": fit_id, "x_variable": kind, "series": label, **fit,
                 "slope_units": units, "slope_ci_low": slope_ci[0], "slope_ci_high": slope_ci[1],
                 "intercept_ci_low": intercept_ci[0], "intercept_ci_high": intercept_ci[1],
                 "departing_groups": int(departing.sum()), "interpretation": verdict})
    for i, group in enumerate(groups):
      mean = group.corrected_od.mean()
      fitted = fit["intercept"] + fit["slope"] * levels[i]
      departures.append({
        "fit_id": fit_id, "x_variable": kind, "series": label, "x": levels[i],
        "source": group.source.iloc[0], "volume_uL": int(group.volume_uL.iloc[0]),
        "wells": ";".join(group.well), "mean_corrected_od": mean, "fitted_od": fitted,
        "mean_departure_od": mean - fitted, "departure_ci_low": departure_ci[0, i],
        "departure_ci_high": departure_ci[1, i], "ci_excludes_zero": bool(departing[i]),
      })
    detail[fit_id] = {"x": x, "y": y, "levels": levels, "groups": groups,
                      "bootstrap_y": by, "bootstrap_slopes": slopes,
                      "bootstrap_intercepts": intercepts, "bootstrap_departures": boot_departures}
  return pd.DataFrame(fits), pd.DataFrame(departures), detail


def verify_analysis(frame: pd.DataFrame, corrected: pd.DataFrame, blanks: pd.DataFrame,
                    draws: np.ndarray, raw_draws: np.ndarray, blank_draws: np.ndarray,
                    fits: pd.DataFrame, departures: pd.DataFrame, detail: dict,
                    inputs: Path, metadata: dict) -> dict:
  """Cross-check OLS and bootstrap calculations using independent SciPy regressions."""
  from scipy.stats import linregress

  max_regression_difference = 0.0
  max_bootstrap_difference = 0.0
  for fit in fits.to_dict("records"):
    values = detail[fit["fit_id"]]
    independent = linregress(values["x"], values["y"])
    for value, expected in [(fit["slope"], independent.slope),
                            (fit["intercept"], independent.intercept),
                            (fit["r_squared"], independent.rvalue ** 2)]:
      require(np.isclose(value, expected, atol=1e-12, rtol=1e-12), "Independent OLS mismatch.")
      max_regression_difference = max(max_regression_difference, abs(value - expected))
    residual = values["y"] - (independent.intercept + independent.slope * values["x"])
    require(np.isclose(fit["residual_rmse_od"], np.sqrt(np.mean(residual ** 2)), atol=1e-12),
            "Independent RMSE mismatch.")
    for index in range(32):
      check = linregress(values["x"], values["bootstrap_y"][index])
      difference = abs(check.slope - values["bootstrap_slopes"][index])
      max_bootstrap_difference = max(max_bootstrap_difference, difference)
      require(difference < 1e-12, "Bootstrap OLS mismatch.")
      independent_departure = values["bootstrap_y"][index].reshape(-1, 3).mean(1) - (
        check.intercept + check.slope * values["levels"])
      require(np.allclose(independent_departure, values["bootstrap_departures"][index],
                          atol=1e-12), "Bootstrap group-departure mismatch.")
  by_well = corrected.set_index("well")
  for record in frame.to_dict("records"):
    volume = record["volume_uL"]
    blank_values = frame[frame.is_blank & (frame.volume_uL == volume)].raw_od_492.tolist()
    independent_blank = math.fsum(blank_values) / 6
    expected = record["raw_od_492"] - independent_blank
    require(math.isclose(by_well.loc[record["well"], "corrected_od"], expected, abs_tol=1e-14),
            "Volume-matched subtraction mismatch.")
  require(np.allclose(raw_draws - draws, blank_draws[:, None, :, None], atol=1e-14),
          "A blank draw was not shared across all concentrations and replicates.")
  repeat = bootstrap_samples(frame)[0]
  require(np.array_equal(repeat, draws), "Seeded bootstrap is not reproducible.")
  require(len(corrected) == 96 and int(corrected.used_in_fits.sum()) == 72, "Lost wells.")
  require(by_well.loc["A8", "corrected_od"] < 0, "Negative result was clipped.")
  for name, expected in metadata["source_sha256"].items():
    require(digest(inputs / name) == expected, "Source copy changed during analysis.")
  return {
    "status": "passed", "well_count": 96, "sample_count": 72, "blank_count": 24,
    "independent_regression": "scipy.stats.linregress for all 10 fits",
    "maximum_ols_difference": max_regression_difference,
    "bootstrap_regression_checks": 32 * len(fits),
    "maximum_bootstrap_slope_difference": max_bootstrap_difference,
    "bootstrap_iterations": len(draws), "seed": SEED, "seed_reproducible": True,
    "shared_blank_draws_checked": True, "all_96_corrections_independently_checked": True,
    "negative_values_preserved": True, "no_wells_excluded": True,
    "source_copy_hashes_unchanged": True,
  }


def save_figures(fits: pd.DataFrame, departures: pd.DataFrame, detail: dict,
                 output: Path) -> list[Path]:
  """Export numeric-axis response and residual panels in PNG and SVG formats."""
  import matplotlib
  matplotlib.use("Agg")
  import matplotlib.pyplot as plt
  from matplotlib.lines import Line2D

  plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
    "axes.titlesize": 12, "axes.labelsize": 10, "axes.spines.top": False,
    "axes.spines.right": False, "axes.edgecolor": "#96a0aa", "axes.labelcolor": "#20364a",
    "text.color": "#20364a", "xtick.color": "#334b5c", "ytick.color": "#334b5c",
    "svg.fonttype": "none", "savefig.facecolor": "white"})
  output.mkdir(exist_ok=True)
  paths = []
  for kind, size in [("volume", (15, 9)), ("concentration", (12, 9))]:
    selected = fits[fits.x_variable == kind]
    columns = 3 if kind == "volume" else 2
    for residual_mode in [False, True]:
      fig, axes = plt.subplots(2, columns, figsize=size)
      for ax, record in zip(axes.flat, selected.to_dict("records")):
        data = detail[record["fit_id"]]
        levels = data["levels"]
        means = np.array([g.corrected_od.mean() for g in data["groups"]])
        deviations = departures[departures.fit_id == record["fit_id"]]
        line = record["intercept"] + record["slope"] * levels
        ax.axhline(0, color="#8896a3", linewidth=0.85, zorder=0)
        ax.grid(axis="y", color="#e8edf1", linewidth=0.7)
        if residual_mode:
          individual = data["y"] - (record["intercept"] + record["slope"] * data["x"])
          ax.scatter(data["x"], individual, s=20, color="#a7b5c1", alpha=0.7, zorder=2)
          for x, y, row in zip(levels, means - line, deviations.to_dict("records")):
            color = "#b85a16" if row["ci_excludes_zero"] else "#14688a"
            ax.vlines(x, row["departure_ci_low"], row["departure_ci_high"], color=color,
                      linewidth=2, zorder=3)
            ax.scatter([x], [y], s=34, color=color, zorder=4)
          ax.set_ylabel("Residual OD at 492 nm")
        else:
          x_line = np.linspace(levels.min(), levels.max(), 100)
          ax.plot(x_line, record["intercept"] + record["slope"] * x_line,
                  color="#254f68", linewidth=1.5, zorder=2)
          for rep, marker in enumerate(["o", "^", "v"]):
            ys = np.array([g.corrected_od.to_numpy()[rep] for g in data["groups"]])
            ax.scatter(levels, ys, marker=marker, s=29, color="#2389a5", alpha=0.75, zorder=3)
          sd = np.array([g.corrected_od.std(ddof=1) for g in data["groups"]])
          ax.errorbar(levels, means, yerr=sd, fmt="s", markersize=4.5, mfc="white",
                      mec="#20364a", ecolor="#20364a", capsize=4, elinewidth=1.4, zorder=4)
          ax.set_ylabel("Blank-subtracted OD at 492 nm")
          ax.text(0.04, 0.96,
                  f"R² = {record['r_squared']:.4f}\nRMSE = {record['residual_rmse_od']:.4f} OD",
                  transform=ax.transAxes, va="top", fontsize=10,
                  bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none", "pad": 2})
        ax.set_title(record["series"].replace("uL", "µL"), loc="left", fontweight="bold", pad=10)
        if kind == "volume":
          ax.set_xticks(VOLUMES)
          ax.set_xlim(40, 212)
          ax.set_xlabel("Intended volume (µL)")
        else:
          ax.set_xlim(-0.015, 1.04)
          ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
          ax.set_xlabel("Relative concentration (A30 = 1)")
        ax.margins(y=0.13)
      prefix = "Residuals" if residual_mode else "Blank-subtracted absorbance"
      suffix = "volume at each concentration" if kind == "volume" else "concentration at each volume"
      fig.suptitle(f"{prefix}: {suffix}", x=0.055, ha="left", fontsize=17, fontweight="bold", y=0.98)
      if residual_mode:
        caption = ("Points: mean departures; bars: pointwise 95% bootstrap intervals. Orange intervals exclude zero. "
                   "Gray dots: individual residuals.\n10,000 resamples; seed 1034; shared volume-matched blank resampling. "
                   "Exploratory intervals, without multiplicity adjustment.")
      else:
        handles = [Line2D([], [], marker="o", color="#2389a5", linestyle="", label="Individual wells (n = 3)"),
                   Line2D([], [], marker="s", mfc="white", color="#20364a", linestyle="", label="Mean ± sample SD"),
                   Line2D([], [], color="#254f68", label="OLS fit, free intercept")]
        fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.045), frameon=False)
        caption = ("Run 151 · 492 nm · Each well minus the mean of its six volume-matched blanks. "
                   "SD bars show sample-well variation; blank uncertainty enters the bootstrap.\n"
                   "Numerical linear axes; each panel has its own y-axis scale. All wells retained.")
      fig.text(0.055, 0.015, caption, fontsize=9, color="#4b5c69", va="bottom")
      fig.subplots_adjust(left=0.065, right=0.98, top=0.9, bottom=0.17 if not residual_mode else 0.14,
                          wspace=0.29, hspace=0.42)
      stem = f"{'residuals' if residual_mode else 'response'}_by_{kind}"
      for extension in ["png", "svg"]:
        path = output / f"{stem}.{extension}"
        fig.savefig(path, dpi=180)
        paths.append(path)
      plt.close(fig)
  return paths


def write_tables(output: Path, audit: pd.DataFrame, corrected: pd.DataFrame,
                 blanks: pd.DataFrame, groups: pd.DataFrame, fits: pd.DataFrame,
                 departures: pd.DataFrame) -> None:
  """Save unrounded numerical results with explicit units and stable identifiers."""
  output.mkdir(exist_ok=True)
  tables = {"sample_id_audit": audit, "well_results": corrected,
            "blank_statistics": blanks, "group_summary": groups,
            "linear_fits": fits, "group_departures": departures}
  for name, frame in tables.items():
    frame.to_csv(output / f"{name}.csv", index=False, float_format="%.12g")


def build_report(output: Path, blanks: pd.DataFrame, corrected: pd.DataFrame,
                 groups: pd.DataFrame, fits: pd.DataFrame, departures: pd.DataFrame,
                 metadata: dict, verification: dict) -> Path:
  """Write a standalone HTML report with embedded figures and numerical findings."""
  import base64

  def table(frame: pd.DataFrame) -> str:
    """Render an escaped, scrollable scientific table."""
    return '<div class="table">' + frame.to_html(index=False, border=0,
      float_format=lambda x: f"{x:.5g}", escape=True) + '</div>'

  def figure(stem: str, caption: str) -> str:
    """Embed a figure so the HTML works without external image files."""
    encoded = base64.b64encode((output / "figures" / f"{stem}.png").read_bytes()).decode("ascii")
    return f'<figure><img src="data:image/png;base64,{encoded}" alt="{html.escape(caption)}"><figcaption>{html.escape(caption)}</figcaption></figure>'

  volume_fits = fits[fits.x_variable == "volume"]
  concentration_fits = fits[fits.x_variable == "concentration"]
  volume_result = volume_fits.set_index("fit_id")
  weakest = groups[groups.source == "A25"].set_index("volume_uL")
  highest = groups[groups.source == "A30"].set_index("volume_uL")
  ratio = highest.loc[200, "mean_corrected_od"] / highest.loc[50, "mean_corrected_od"]
  display_fits = fits[["series", "x_variable", "r_squared", "slope", "slope_ci_low",
                       "slope_ci_high", "intercept", "residual_rmse_od", "interpretation"]].copy()
  key_departures = departures[departures.ci_excludes_zero][[
    "series", "x_variable", "source", "volume_uL", "wells", "mean_departure_od",
    "departure_ci_low", "departure_ci_high"]]
  quiet = corrected[~corrected.is_blank & (corrected.signal_flag != "")][[
    "well", "sample_id", "raw_od_492", "blank_mean_od", "corrected_od", "signal_flag"]]
  low_groups = groups[groups.source.isin(["A26", "A25"])][[
    "source", "volume_uL", "wells", "mean_corrected_od", "sample_sd_od", "cv_percent",
    "largest_replicate_deviation_well", "largest_replicate_deviation_od"]]
  concentration_departures = concentration_fits[concentration_fits.departing_groups > 0].series.tolist()
  content = f'''<!doctype html><html lang="en"><meta charset="utf-8">
<title>Exp1034 Orange G analysis — run 151</title>
<style>body{{font:16px/1.55 system-ui,sans-serif;color:#20364a;max-width:1250px;margin:38px auto;padding:0 24px}}
h1{{font-size:32px;line-height:1.2}}h2{{margin-top:38px;font-size:23px}}h3{{font-size:19px}}
.lede{{font-size:18px}}.meta,figcaption{{color:#526575;font-size:14px}}.table{{overflow-x:auto;margin:18px 0}}
table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{padding:9px 10px;text-align:right;border-bottom:1px solid #dce4e9;white-space:nowrap}}
th{{background:#edf3f6}}td:first-child,th:first-child{{text-align:left}}figure{{margin:25px 0}}img{{width:100%;height:auto}}
code{{background:#edf3f6;padding:2px 4px}}.key{{border-left:4px solid #177c99;padding:12px 20px;background:#f3f8fa}}
a{{color:#096b90}}@media print{{body{{max-width:none;margin:0}}figure{{break-inside:avoid}}}}</style>
<h1>Exp1034 Orange G plate analysis</h1>
<p class="meta">FLUOstar OPTIMA run 151 · 6 October 2026, 10:44:04 AM (export time) · 492 nm · path-length correction off, confirmed by operator</p>
<div class="key"><p class="lede"><strong>All 96 electronic sample IDs are in the intended wells.</strong>
The data show strong concentration and volume trends, with small but detectable departures in several fits.
The lowest concentration is dominated by variability at the smaller volumes.</p></div>
<h2>1. Sample IDs and blank correction</h2>
<p>The raw-data export, fixed-width sample-ID file and original dispensing manifest agree for all 96 wells,
including source, intended volume and replicate. The TXT Content numbers also match the expected X1–X96 mapping.
There are 72 samples (three per concentration/volume group) and 24 blanks (six per volume).
This establishes correct electronic mapping, not independent verification of the physical liquid in each well.</p>
<p>Each sample was corrected as <strong>raw OD − mean blank OD at that same volume</strong>.
All individual readings remain in the analysis, including negative corrected values.</p>
{table(blanks[["volume_uL", "n", "blank_wells", "blank_mean_od", "blank_sd_od", "blank_sem_od"]])}
<h2>2. Absorbance versus volume at each concentration</h2>
<p>A30–A25 are labelled relative to A30: 1, 1/2, 1/4, 1/8, 1/16 and 1/32. These are confirmed dilution ratios,
not the DNA concentrations from the original experiment design.</p>
<p>A30 has R² = {volume_result.loc['volume_A30', 'r_squared']:.4f}, but the 100 and 200 µL group departures
exclude zero in their pointwise bootstrap intervals. A29–A26 show no clear departure by this rule;
their R² values range from {volume_result.loc['volume_A26', 'r_squared']:.4f} to {volume_result.loc['volume_A29', 'r_squared']:.4f}.
A25 has R² = {volume_result.loc['volume_A25', 'r_squared']:.4f}; its lack of a detected departure does not make it a precise calibration.</p>
{figure('response_by_volume', 'Individual wells, group means ± sample SD and free-intercept OLS lines. Each panel uses its own y-axis scale.')}
<p><strong>A straight-line fit is not the same as proportional response.</strong> For A30, the corrected mean is
{highest.loc[50, 'mean_corrected_od']:.4f} OD at 50 µL and {highest.loc[200, 'mean_corrected_od']:.4f} OD at 200 µL:
a {ratio:.2f}-fold signal change for a fourfold intended-volume change. Its fitted intercept is
{volume_result.loc['volume_A30', 'intercept']:.4f} OD. These data alone cannot separate dispensing effects from
optical/path-length effects or establish actual delivered volumes.</p>
{figure('residuals_by_volume', 'Volume-fit group departures with pointwise 95% bootstrap intervals; orange intervals exclude zero.')}
<h2>3. Dilution curves at each fixed volume</h2>
<p>The concentration fits have R² values from {concentration_fits.r_squared.min():.4f} to {concentration_fits.r_squared.max():.4f}.
The exploratory departure rule flags {', '.join(concentration_departures)}; 75 µL shows no clear departure.
These are small departures that can coexist with high R² values.</p>
{figure('response_by_concentration', 'Absorbance versus relative concentration at 50, 75, 100 and 200 µL. Both axes are linear.')}
{figure('residuals_by_concentration', 'Concentration-fit residuals; close low-concentration points retain their true numerical x positions.')}
<h2>4. Fit estimates and groups driving departures</h2>
<p>OLS uses every sample well, with a freely estimated intercept. R² and RMSE are calculated across individual
wells (n = 12 per volume fit; n = 18 per concentration fit). RMSE = √(sum of squared residuals / n).
Slopes have units OD/µL for volume fits and OD per relative-concentration unit for concentration fits.</p>
{table(display_fits)}
<p>Groups below have pointwise 95% intervals excluding zero. These are group effects, not proof that one well failed.</p>
{table(key_departures)}
<h2>5. Low signals and replicate variation</h2>
<p>A25 at 50 µL averages {weakest.loc[50, 'mean_corrected_od']:.4f} OD above blank. A8 is
{corrected.set_index('well').loc['A8', 'corrected_od']:.4f} OD after correction; this is a measured value below the
mean blank, not evidence by itself of an empty well. A25's 100 µL mean
({weakest.loc[100, 'mean_corrected_od']:.4f}) is below its 75 µL mean
({weakest.loc[75, 'mean_corrected_od']:.4f}). Blank and sample variation are substantial relative to these signals.</p>
{table(quiet)}
<p>The “within 3 blank SD” flag is a descriptive screen, not an established detection limit.
The following table identifies the largest within-group replicate deviations for the two lowest concentrations.
CV becomes unstable when the corrected mean is close to zero; no well is excluded.</p>
{table(low_groups)}
<h2>6. Method, uncertainty and reproducibility</h2>
<p>10,000 stratified bootstrap resamples, seed 1034. Each replicate resamples the six blanks within each volume
and the three samples within each concentration/volume group. Each resampled blank mean is shared across all
corresponding samples. The lines are refitted and 2.5th/97.5th percentiles form slope and group-departure intervals.
Sample-SD error bars describe sample-well spread; blank uncertainty is included in the bootstrap intervals.</p>
<p>Pointwise intervals are exploratory and not adjusted for multiple comparisons. Three sample wells per group
and one dilution preparation limit generalization. “No clear departure detected” means no detected departure by
this rule, not proof of linearity. A slope interval spanning zero would be labelled too weak to assess;
none of the ten slope intervals does so here. Instrument saturation, physical well contents and absolute
dispensed volumes are not independently established.</p>
<p>Verification passed: all 96 subtractions independently reproduced; all 10 regressions matched SciPy
(maximum coefficient/R² difference {verification['maximum_ols_difference']:.2g}); 320 bootstrap regressions
cross-checked; deterministic resampling and shared blank draws verified. Original source hashes were preserved.</p>
<p><a href="OrangeG_analysis.ipynb">Executed analysis notebook</a> ·
<a href="tables/well_results.csv">All well results</a> · <a href="tables/linear_fits.csv">Fit results</a> ·
<a href="tables/blank_statistics.csv">Blank statistics</a> · <a href="tables/sample_id_audit.csv">Sample-ID audit</a> ·
<a href="tables/group_summary.csv">Group summary</a> · <a href="tables/group_departures.csv">Group departures</a></p>
<p class="meta">Sources: OrangeG_dye_test_raw-data_w_sampleIDs.CSV, Exp1034_OrangeG_Sample_IDs_fixed.txt and
well_assignments.json. Exact source copies are retained in inputs/; hashes and acquisition metadata are in analysis_metadata.json.</p>
</html>'''
  path = output / "OrangeG_analysis_report.html"
  path.write_text(content, encoding="utf-8")
  return path
