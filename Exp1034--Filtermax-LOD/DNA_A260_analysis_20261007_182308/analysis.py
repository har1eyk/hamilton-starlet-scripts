"""Reproduce the Exp1034 A260 calculations from immutable input copies; no hardware I/O."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
VOLUMES = np.array([50, 75, 100, 200])
STOCK_NG_UL = 70.0  # User correction on 2026-10-07; original exports retain 16.
CONCENTRATION_SCALE = STOCK_NG_UL / 16.0
CONCENTRATIONS = np.array([0.5, 1, 2, 4, 8, 16]) * CONCENTRATION_SCALE
ROWS = "ABCDEFGH"


def save_json(path: Path, data: dict) -> None:
  """Write readable JSON with native scalar values."""
  path.write_text(json.dumps(data, indent=2, allow_nan=False,
                             default=lambda value: value.item()) + "\n")


def ols(x: np.ndarray, y: np.ndarray) -> dict:
  """Fit a free-intercept line to individual observations, not a forced-zero line."""
  x, y = np.asarray(x, float), np.asarray(y, float)
  intercept, slope = np.linalg.lstsq(np.column_stack([np.ones(len(x)), x]), y, rcond=None)[0]
  residual = y - intercept - slope * x
  sse = float(residual @ residual)
  return {"n": len(x), "slope": slope, "intercept": intercept,
          "r_squared": 1 - sse / np.sum((y-y.mean())**2),
          "residual_rmse": np.sqrt(sse / len(x)),
          "residual_sd_n_minus_2": np.sqrt(sse / (len(x)-2))}


def bootstrap_line(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
  """Fit each bootstrap's equally replicated group means; equivalent to individual-well OLS."""
  centered = x - x.mean()
  slope = (y @ centered) / (centered @ centered)
  return slope, y.mean(axis=1) - slope*x.mean()


def interval(values: np.ndarray) -> tuple[float, float]:
  """Return the conditional within-plate percentile bootstrap 95% interval."""
  lo, hi = np.quantile(values, [0.025, 0.975])
  return float(lo), float(hi)


def analyze(root: Path = HERE, n_boot: int = 10_000, seed: int = 1034) -> dict:
  """Audit identities, subtract matched blanks, fit lines and export every original observation."""
  inp, out = root / "inputs", root / "tables"
  out.mkdir(exist_ok=True)
  raw = pd.read_csv(inp / "raw_a260.csv")
  manifest = pd.read_csv(inp / "well_manifest.csv")
  metadata = json.loads((inp / "reader_metadata.json").read_text())
  expected = []
  for volume, concentration, source, wells in json.loads((inp / "well_assignments.json").read_text()):
    for replicate, well in enumerate(wells, 1):
      expected.append({"well": well, "sample_id": f"{source}_{volume:03d}uL_R{replicate}",
                       "source_tube": source, "planned_volume_uL": volume,
                       "planned_dna_concentration_ng_uL": concentration, "replicate": replicate,
                       "replicate_type": "dispensing",
                       "role": "buffer_blank" if concentration == 0 else "dna_sample"})
  expected = pd.DataFrame(expected).set_index("well").sort_index()
  assert len(raw) == raw.well.nunique() == len(manifest) == manifest.well.nunique() == 96
  assert set(raw.well) == {f"{r}{c}" for r in ROWS for c in range(1, 13)}
  audit = pd.DataFrame(index=expected.index)
  for name, frame in (("raw", raw), ("run_manifest", manifest)):
    ordered = frame.set_index("well").loc[expected.index]
    for field in expected.columns:
      audit[f"{name}_{field}_matches"] = ordered[field] == expected[field]
  assert audit.all().all(), "Electronic identity mismatch; review before analysis."
  audit["all_labels_match"] = audit.all(axis=1)
  audit["derived_content_number"] = [f"X{ROWS.index(w[0])*12 + int(w[1:])}" for w in audit.index]
  audit["content_number_note"] = "Derived row-major number; raw A260 export has no Content field"
  audit.reset_index().to_csv(out / "sample_id_audit.csv", index=False)
  assert set(raw.wavelength_nm) == {260} and metadata["wavelength"] == 260
  assert metadata["status"] == "complete" and metadata["missing_well_count"] == 0
  assert metadata["reference_subtracted"] is False
  assert metadata["pathlength_correction_applied"] is False
  assert np.isfinite(raw.a260).all() and raw.read_id.nunique() == 1
  assert raw.valid.isna().all() and raw.exclusion_reason.isna().all()
  assert hashlib.sha256((inp / "well_assignments.json").read_bytes()).hexdigest() == metadata["assignment_sha256"]
  payload = json.loads((inp / "raw_backend_results.json").read_text())
  assert len(payload) == 1 and payload[0]["wavelengths"] == [260]
  assert payload[0]["reference_subtracted"] is False
  for record in raw.itertuples():
    np.testing.assert_allclose(record.a260, payload[0]["data"][ROWS.index(record.well[0])][int(record.well[1:])-1], atol=1e-12)
  d = raw.rename(columns={"planned_volume_uL": "volume_uL"}).copy()
  d["concentration_ng_uL"] = d.planned_dna_concentration_ng_uL * CONCENTRATION_SCALE
  d["concentration_basis"] = "User-reported A30 70 ng/uL; nominal twofold dilutions, not inferred from A260"
  audit["exported_concentration_ng_uL"] = expected.planned_dna_concentration_ng_uL
  audit["analysis_concentration_ng_uL"] = expected.planned_dna_concentration_ng_uL * CONCENTRATION_SCALE
  audit["concentration_metadata_corrected"] = expected.planned_dna_concentration_ng_uL > 0
  audit.reset_index().to_csv(out / "sample_id_audit.csv", index=False)
  d["mass_ng"] = d.volume_uL*d.concentration_ng_uL
  d["plate_row"] = d.well.str[0]
  d["plate_column"] = d.well.str[1:].astype(int)
  d["analysis_included"] = True
  d["review_flag"] = ""
  d.loc[d.well == "D9", "review_flag"] = f"Elevated within the 50 uL / {CONCENTRATIONS[0]:g} ng/uL triplicate; retained"
  rng = np.random.default_rng(seed)
  boot = np.empty((n_boot, 4, 6))
  boot_blank_means, boot_blank_sd = np.empty((n_boot, 4)), np.empty((n_boot, 4))
  blanks, groups = [], []
  for vi, volume in enumerate(VOLUMES):
    b = d[(d.volume_uL == volume) & (d.concentration_ng_uL == 0)]
    assert len(b) == 6
    bv = b.a260.to_numpy()
    bdraw = rng.choice(bv, (n_boot, 6), replace=True)
    boot_blank_means[:, vi] = bdraw.mean(axis=1)
    boot_blank_sd[:, vi] = bdraw.std(axis=1, ddof=1)
    mean, sd = bv.mean(), bv.std(ddof=1)
    blo, bhi = interval(boot_blank_means[:, vi])
    slo, shi = interval(boot_blank_sd[:, vi])
    blanks.append({"volume_uL": volume, "n": 6, "blank_mean": mean, "blank_sd": sd,
                   "mean_ci95_low": blo, "mean_ci95_high": bhi,
                   "sd_boot_ci95_low": slo, "sd_boot_ci95_high": shi,
                   "raw_threshold_3_3sd": mean+3.3*sd, "raw_threshold_10sd": mean+10*sd,
                   "blank_wells": ";".join(b.well)})
    mask = d.volume_uL == volume
    d.loc[mask, "blank_mean"] = mean
    d.loc[mask, "blank_sd"] = sd
    d.loc[mask, "corrected_a260"] = d.loc[mask, "a260"] - mean
    for ci, concentration in enumerate(CONCENTRATIONS):
      g = d[mask & (d.concentration_ng_uL == concentration)]
      assert len(g) == 3
      values = g.a260.to_numpy()
      # Each resampled blank mean is shared by all six concentrations at its volume.
      boot[:, vi, ci] = rng.choice(values, (n_boot, 3), replace=True).mean(axis=1) - boot_blank_means[:, vi]
      corrected = values - mean
      m, s = corrected.mean(), corrected.std(ddof=1)
      lo, hi = interval(boot[:, vi, ci])
      groups.append({"volume_uL": volume, "concentration_ng_uL": concentration,
                     "source_tube": g.source_tube.iloc[0], "mass_ng": volume*concentration,
                     "n": 3, "mean_raw": values.mean(), "mean_corrected": m,
                     "sample_sd": s, "cv_percent": 100*s/m if m > 0 else np.nan,
                     "mean_ci95_low": lo, "mean_ci95_high": hi,
                     "mean_in_blank_sd": m/sd, "replicates_above_3_3sd": int(np.sum(corrected > 3.3*sd)),
                     "replicates_above_10sd": int(np.sum(corrected > 10*sd)),
                     "strict_screen_10sd_and_cv10": bool(np.all(corrected > 10*sd) and 100*s/m <= 10),
                     "wells": ";".join(g.well)})
  b, g = pd.DataFrame(blanks), pd.DataFrame(groups)
  exported = pd.read_csv(inp / "blank_summary.csv").set_index("planned_volume_uL")
  np.testing.assert_allclose(exported.loc[VOLUMES, "blank_mean_a260"], b.blank_mean, atol=1e-12)
  np.testing.assert_allclose(exported.loc[VOLUMES, "blank_sample_sd_a260"], b.blank_sd, atol=1e-12)
  existing = pd.read_csv(inp / "blank_corrected_a260.csv").set_index("well")
  np.testing.assert_allclose(existing.loc[d.well, "blank_corrected_a260"], d.corrected_a260, atol=1e-12)
  samples = d[d.concentration_ng_uL > 0]
  fits, departures, limit_rows = [], [], []
  for vi, volume in enumerate(VOLUMES):
    for fit_kind, selected in (("concentration_full", np.ones(6, dtype=bool)),
                               ("concentration_low", CONCENTRATIONS <= CONCENTRATIONS[3])):
      x = CONCENTRATIONS[selected]
      subset = samples[(samples.volume_uL == volume) & samples.concentration_ng_uL.isin(x)]
      f = ols(subset.concentration_ng_uL.to_numpy(), subset.corrected_a260.to_numpy())
      sb, ib = bootstrap_line(x, boot[:, vi, selected])
      sl, sh = interval(sb)
      il, ih = interval(ib)
      f.update(fit_kind=fit_kind, volume_uL=volume, concentration_ng_uL=np.nan,
               slope_ci95_low=sl, slope_ci95_high=sh, intercept_ci95_low=il, intercept_ci95_high=ih,
               fitted_concentrations=";".join(f"{c:g}" for c in x))
      fits.append(f)
      if fit_kind == "concentration_full":
        delta_boot = boot[:, vi, :] - ib[:, None] - sb[:, None]*CONCENTRATIONS
        means = g[g.volume_uL == volume].mean_corrected.to_numpy()
        for ci, concentration in enumerate(CONCENTRATIONS):
          lo, hi = interval(delta_boot[:, ci])
          departures.append({"fit_kind": fit_kind, "volume_uL": volume, "concentration_ng_uL": concentration,
                             "departure": means[ci] - f["intercept"] - f["slope"]*concentration,
                             "ci95_low": lo, "ci95_high": hi, "excludes_zero": lo > 0 or hi < 0})
        mask = d.volume_uL == volume
        d.loc[mask, "full_fit_prediction"] = f["intercept"] + f["slope"]*d.loc[mask, "concentration_ng_uL"]
      valid = sb > 0
      for multiplier, name in ((3.3, "DL"), (10, "QL")):
        estimates = multiplier*boot_blank_sd[valid, vi]/sb[valid]
        lo, hi = interval(estimates)
        estimate = multiplier*b.iloc[vi].blank_sd/f["slope"]
        limit_rows.append({"volume_uL": volume, "fit_kind": fit_kind, "metric": name,
                           "estimate_ng_uL": estimate, "ci95_low": lo, "ci95_high": hi,
                           "mass_ng": estimate*volume, "bootstrap_positive_slope_fraction": valid.mean(),
                           "status": "Exploratory blank-noise estimate; no independent validation"})
  for ci, concentration in enumerate(CONCENTRATIONS):
    subset = samples[samples.concentration_ng_uL == concentration]
    f = ols(subset.volume_uL.to_numpy(), subset.corrected_a260.to_numpy())
    sb, ib = bootstrap_line(VOLUMES, boot[:, :, ci])
    sl, sh = interval(sb)
    il, ih = interval(ib)
    f.update(fit_kind="volume", volume_uL=np.nan, concentration_ng_uL=concentration,
             slope_ci95_low=sl, slope_ci95_high=sh, intercept_ci95_low=il, intercept_ci95_high=ih,
             fitted_concentrations="")
    fits.append(f)
    delta_boot = boot[:, :, ci] - ib[:, None] - sb[:, None]*VOLUMES
    for vi, volume in enumerate(VOLUMES):
      lo, hi = interval(delta_boot[:, vi])
      m = g[(g.volume_uL == volume) & (g.concentration_ng_uL == concentration)].mean_corrected.iloc[0]
      departures.append({"fit_kind": "volume", "volume_uL": volume, "concentration_ng_uL": concentration,
                         "departure": m - f["intercept"] - f["slope"]*volume,
                         "ci95_low": lo, "ci95_high": hi, "excludes_zero": lo > 0 or hi < 0})
  f, dep, limits = pd.DataFrame(fits), pd.DataFrame(departures), pd.DataFrame(limit_rows)
  d["model_residual"] = d.corrected_a260 - d.full_fit_prediction
  d.loc[d.concentration_ng_uL == 0, "model_residual"] = d.loc[d.concentration_ng_uL == 0, "corrected_a260"]
  d["within_group_residual"] = d.corrected_a260 - d.groupby(["volume_uL", "concentration_ng_uL"]).corrected_a260.transform("mean")
  d["signal_in_blank_sd"] = d.corrected_a260 / d.blank_sd
  decisions = []
  for volume in VOLUMES:
    t = g[g.volume_uL == volume]
    low = t[t.concentration_ng_uL == CONCENTRATIONS[0]].iloc[0]
    eligible = t[t.strict_screen_10sd_and_cv10]
    decisions.append({"volume_uL": volume, "lowest_tested_ng_uL": CONCENTRATIONS[0],
                      "lowest_tested_mass_ng": CONCENTRATIONS[0]*volume,
                      "cv_at_lowest_percent": low.cv_percent,
                      "signal_in_blank_sd_at_lowest": low.mean_in_blank_sd,
                      "all_3_above_blank_screen_at_lowest": low.replicates_above_3_3sd == 3,
                      "lowest_strict_screen_ng_uL": eligible.concentration_ng_uL.min(),
                      "strict_screen_definition": "All 3 wells > 10 blank SD and corrected CV <= 10%; exploratory, not validated LOQ"})
  sensitivity = []
  for omit in (None, "D9"):
    t = samples[samples.volume_uL == 50]
    if omit:
      t = t[t.well != omit]
    low = t[t.concentration_ng_uL == CONCENTRATIONS[0]].corrected_a260
    sensitivity.append({"analysis": "Primary: all wells" if omit is None else "Sensitivity only: omit D9",
                        "low_group_n": len(low), "low_mean": low.mean(), "low_sd": low.std(ddof=1),
                        "low_cv_percent": 100*low.std(ddof=1)/low.mean(),
                        **ols(t.concentration_ng_uL.to_numpy(), t.corrected_a260.to_numpy())})
  ratio_rows = []
  for vi, volume in enumerate(VOLUMES):
    for ci, concentration in enumerate(CONCENTRATIONS):
      numerator = g[(g.volume_uL == volume) & (g.concentration_ng_uL == concentration)].mean_corrected.iloc[0]/volume
      denominator = g[(g.volume_uL == 100) & (g.concentration_ng_uL == concentration)].mean_corrected.iloc[0]/100
      lo, hi = interval((boot[:, vi, ci]/volume)/(boot[:, 2, ci]/100))
      ratio_rows.append({"volume_uL": volume, "concentration_ng_uL": concentration,
                         "relative_volume_normalized_signal": numerator/denominator, "ci95_low": lo, "ci95_high": hi})
  for name, frame in (("well_results", d), ("blank_statistics", b), ("group_summary", g),
                      ("linear_fits", f), ("group_departures", dep), ("exploratory_limits", limits),
                      ("volume_decisions", pd.DataFrame(decisions)), ("d9_sensitivity", pd.DataFrame(sensitivity)),
                      ("volume_normalized_ratios", pd.DataFrame(ratio_rows))):
    frame.to_csv(out / f"{name}.csv", index=False,
                 float_format=None if name == "well_results" else "%.12g")
  np.savez_compressed(root / "bootstrap_results.npz", corrected_means=boot,
                      blank_means=boot_blank_means, blank_sd=boot_blank_sd)
  summary = {"run": "dna-a260-f5-20261007T182308Z-4a000d5d", "wavelength_nm": 260,
             "seed": seed, "bootstrap_resamples": n_boot, "wells": 96, "samples": 72, "blanks": 24,
             "excluded_wells": [], "all_electronic_labels_match": True,
             "user_reported_A30_ng_uL": STOCK_NG_UL,
             "original_export_A30_ng_uL": 16.0,
             "concentration_scale_correction": CONCENTRATION_SCALE,
             "nominal_concentrations_ng_uL": CONCENTRATIONS.tolist(),
             "concentration_note": "User corrected A30 to 70 ng/uL. Original concentration metadata is preserved and superseded for analysis; other tubes assume the same twofold dilution. No independent concentration verification.",
             "physical_contents_independently_verified": False, "pathlength_correction": False,
             "reader_validity_review": "pending in original export; all finite readings retained for exploratory analysis",
             "signal_min_raw": float(d.a260.min()), "signal_max_raw": float(d.a260.max()),
             "decisions": decisions, "blank_statistics": blanks,
             "bootstrap_scope": "Resample 3 wells within each condition and 6 blanks within volume; reuse each resampled blank mean. Does not cover plate/day, dilution-chain or concentration uncertainty.",
             "linearity_rule": "Free-intercept OLS; group departure 95% bootstrap intervals excluding zero are exploratory, unadjusted for multiple comparisons.",
             "limit_rule": "3.3 or 10 times sample blank SD divided by slope; full 2.1875-70 and fixed lower-four-level 2.1875-17.5 ng/uL ranges. Extrapolation below 2.1875 is untested.",
             "sources": [
               {"title": "Corning 3635 product specification", "url": "https://ecatalog.corning.com/life-sciences/b2b/SEA/en/Microplates/Assay-Microplates/96-Well-Microplates/Corning%C2%AE-96-well-UV-Transparent-Microplates/p/3635"},
               {"title": "ICH Q2(R2), sections 3.2.3.3-3.2.3.5", "url": "https://database.ich.org/sites/default/files/ICH_Q2%28R2%29_Guideline_2023_1130.pdf"},
               {"title": "Molecular Devices: DNA/RNA absorbance and path length", "url": "https://www.moleculardevices.com/en/assets/app-note/br/dna-and-rna-absorbance-measurements-using-spectramax-microplate-readers"}]}
  save_json(root / "analysis_metadata.json", summary)
  save_json(root / "validation_checks.json", {"complete_well_coverage": True, "manifest_matches": True,
            "raw_backend_mapping_matches": True, "assignment_hash_matches": True,
            "reader_blank_calculations_reproduced": True, "reader_corrected_values_reproduced": True,
            "source_files_modified": False, "no_wells_excluded": True})
  return summary


if __name__ == "__main__":
  result = analyze()
  print(json.dumps(result["decisions"], indent=2, default=lambda value: value.item()))
