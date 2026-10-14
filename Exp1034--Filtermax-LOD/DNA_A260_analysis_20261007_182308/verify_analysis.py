"""Independent numerical/provenance checks; reads files only, except validation_checks.json."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import linregress

HERE = Path(__file__).resolve().parent


def verify(root=HERE):
  checks = json.loads((root / 'validation_checks.json').read_text())
  original_checks = []
  for item in json.loads((root / 'source_provenance.json').read_text()):
    assert hashlib.sha256((root / item['copy']).read_bytes()).hexdigest() == item['sha256']
    source = Path(item['source'])
    if source.exists():
      assert hashlib.sha256(source.read_bytes()).hexdigest() == item['sha256']
      original_checks.append(True)
  d = pd.read_csv(root / 'tables/well_results.csv')
  raw = pd.read_csv(root / 'inputs/raw_a260.csv').set_index('well')
  assert d.well.nunique() == len(d) == 96
  assert d.analysis_included.all()
  # CSV parsers can differ by one floating-point rounding unit; byte-level raw
  # preservation is checked above, separately from this numerical comparison.
  np.testing.assert_allclose(d.a260, raw.loc[d.well, 'a260'], rtol=0, atol=5e-16)
  by_source = {'A30': 70, 'A29': 35, 'A28': 17.5, 'A27': 8.75,
               'A26': 4.375, 'A25': 2.1875, 'A31': 0}
  np.testing.assert_array_equal(d.concentration_ng_uL, d.source_tube.map(by_source))
  assert set(d[d.concentration_ng_uL == 0].volume_uL.value_counts()) == {6}
  for volume, t in d.groupby('volume_uL'):
    blanks = t[t.concentration_ng_uL == 0].a260
    np.testing.assert_allclose(t.corrected_a260, t.a260 - blanks.mean(), atol=1e-11, rtol=0)
  sample = d[d.concentration_ng_uL > 0]
  fits = pd.read_csv(root / 'tables/linear_fits.csv')
  for f in fits.itertuples():
    if f.fit_kind == 'volume':
      t = sample[sample.concentration_ng_uL == f.concentration_ng_uL]
      x = t.volume_uL
    else:
      t = sample[sample.volume_uL == f.volume_uL]
      if f.fit_kind == 'concentration_low':
        t = t[t.concentration_ng_uL <= 17.5]
      x = t.concentration_ng_uL
    r = linregress(x, t.corrected_a260)
    rmse = np.sqrt(np.mean((t.corrected_a260 - r.intercept - r.slope*x)**2))
    np.testing.assert_allclose([f.slope, f.intercept, f.r_squared, f.residual_rmse],
                              [r.slope, r.intercept, r.rvalue**2, rmse], rtol=1e-8, atol=1e-11)
  boot = np.load(root / 'bootstrap_results.npz')
  assert boot['corrected_means'].shape == (10000, 4, 6)
  assert boot['blank_means'].shape == (10000, 4)
  # Use explicit algebra to verify the limits separately from the fitter.
  blank = pd.read_csv(root / 'tables/blank_statistics.csv').set_index('volume_uL')
  for r in pd.read_csv(root / 'tables/exploratory_limits.csv').itertuples():
    slope = fits[(fits.volume_uL == r.volume_uL) & (fits.fit_kind == r.fit_kind)].slope.iloc[0]
    expected = (3.3 if r.metric == 'DL' else 10)*blank.loc[r.volume_uL, 'blank_sd']/slope
    np.testing.assert_allclose(r.estimate_ng_uL, expected, rtol=1e-8)
  checks.update(independent_scipy_regressions_verified=len(fits),
                preserved_input_hashes_verified=7, original_source_hashes_verified=len(original_checks),
                user_corrected_concentration_mapping_verified=True,
                negative_blank_corrected_readings_retained=int((d.corrected_a260 < 0).sum()),
                bootstrap_shape_verified=True, exploratory_limit_arithmetic_verified=True)
  (root / 'validation_checks.json').write_text(json.dumps(checks, indent=2) + '\n')
  return checks


if __name__ == '__main__':
  print(json.dumps(verify(), indent=2))
