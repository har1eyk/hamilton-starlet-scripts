Exp1034 DNA A260 analysis - run 20261007T182308Z-4a000d5d

START HERE
  output/pdf/DNA_A260_lab_notebook.pdf - printable report with 15 figures
  DNA_A260_lab_notebook.html - browsable report; keep figures/ beside it
  DNA_A260_analysis.ipynb - executed, reproducible analysis and figures
  figures/ - every plot as high-resolution PNG and editable SVG
  tables/ - all audited well-level data, summaries, fits and intervals

CONCENTRATION CORRECTION
  User-reported A30 is 70 ng/uL, replacing the original export's 16 ng/uL.
  A30/A29/A28/A27/A26/A25 = 70/35/17.5/8.75/4.375/2.1875 ng/uL.
  These downstream values assume the same twofold dilution chain.
  Original files under inputs/ are unchanged. well_results.csv preserves
  planned_dna_concentration_ng_uL from the export and adds concentration_ng_uL
  for the corrected analysis. Read concentration_correction.json for provenance.

REPRODUCE
  Install requirements.txt into an analysis Python environment. No hardware
  software is needed. From this folder:
    python analysis.py
    python plot_analysis.py
    python verify_analysis.py
    python build_report.py
  Or run all cells in DNA_A260_analysis.ipynb for calculations and figures.
  build_report.py uses DejaVu fonts at /usr/share/fonts/truetype/dejavu;
  adjust FONT for another operating system. Original absolute source paths
  need not exist after download: input-copy hashes are always verified.

METHODS
  Separate six-blank means and sample SDs at 50, 75, 100 and 200 uL.
  Corrected OD = raw OD - mean blank OD at that volume. No exclusions.
  OLS with free intercept, fitted to individual sample wells.
  Full concentration range: 2.1875-70 ng/uL; lower four: 2.1875-17.5.
  10,000 bootstrap resamples, seed 1034, sharing each sampled blank mean
  across every concentration at its volume. Intervals are within-plate only.
  Exploratory DL = 3.3 blank SD / slope; QL = 10 blank SD / slope.
  These formulas do not validate detection or quantitation below tested levels.

TABLES
  sample_id_audit: identity comparisons and original vs corrected concentrations.
  well_results: all 96 raw/corrected wells; user-corrected nominal concentration;
    planned volume; source; replicate; original validity fields; residuals.
  blank_statistics: six blanks per volume, mean, SD and bootstrap intervals.
  group_summary: 24 groups (n=3), mean, SD, CV, signal/blank SD and thresholds.
  linear_fits: slope/intercept/R2/RMSE and bootstrap slope/intercept intervals.
    residual_rmse = sqrt(SSE/n); residual_sd_n_minus_2 = sqrt(SSE/(n-2)).
  group_departures: mean observed minus fitted response and bootstrap intervals.
  exploratory_limits: noise-based DL/QL, full and lower-range slopes.
  volume_decisions: low-level repeatability and illustrative strict screen.
  d9_sensitivity: all-well primary analysis vs secondary D9 omission.
  volume_normalized_ratios: OD/volume relative to 100 uL at each concentration.

LIMITATIONS
  One plate, one read, one serial dilution chain. Three wells per sample group
  are dispensing replicates, not independent dilutions or days. Original
  validity review is pending. No path-length correction, independent stock
  validation or actual fill measurement is available. No purity or actual DNA
  concentration is inferred. D9 is flagged and retained. The sample-ID audit
  proves electronic agreement, not physical contents. Historical Orange G
  readings are not pooled with this DNA read. No hardware was operated.
