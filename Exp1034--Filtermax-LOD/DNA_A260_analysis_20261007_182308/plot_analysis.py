"""Export scientific PNG/SVG figures from the audited A260 analysis tables."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LogNorm, TwoSlopeNorm
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
V = [50, 75, 100, 200]
from analysis import CONCENTRATIONS, CONCENTRATION_SCALE
C = CONCENTRATIONS.tolist()
SCALE = CONCENTRATION_SCALE
COLORS = {50: "#D97732", 75: "#8760A6", 100: "#008579", 200: "#2875B6"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.titlesize": 12,
                     "axes.labelsize": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#E6EAED", "grid.linewidth": 0.6,
                     "axes.axisbelow": True, "svg.fonttype": "none", "savefig.facecolor": "white"})


def plate_array(d: pd.DataFrame, field: str) -> np.ndarray:
  """Map wells explicitly to the row-major plate image."""
  arr = np.full((8, 12), np.nan)
  for r in d.itertuples():
    arr["ABCDEFGH".index(r.well[0]), int(r.well[1:])-1] = d.loc[d.well == r.well, field].iloc[0]
  return arr


def axes_plate(ax: plt.Axes) -> None:
  """Display physical well coordinates rather than zero-based array indices."""
  ax.set_xticks(range(12), range(1, 13))
  ax.set_yticks(range(8), list("ABCDEFGH"))
  ax.set_xlabel("Plate column")
  ax.set_ylabel("Plate row")
  ax.grid(False)
  ax.set_xticks(np.arange(-.5, 12, 1), minor=True)
  ax.set_yticks(np.arange(-.5, 8, 1), minor=True)
  ax.grid(which="minor", color="white", linewidth=1)
  ax.tick_params(which="minor", bottom=False, left=False)


def make_figures(root: Path = HERE) -> list[dict]:
  """Render complementary figures without changing source data or fitting ranges."""
  out = root / "figures"
  out.mkdir(exist_ok=True)
  d = pd.read_csv(root / "tables/well_results.csv")
  s = d[d.concentration_ng_uL > 0]
  b = pd.read_csv(root / "tables/blank_statistics.csv").set_index("volume_uL")
  g = pd.read_csv(root / "tables/group_summary.csv")
  fits = pd.read_csv(root / "tables/linear_fits.csv")
  dep = pd.read_csv(root / "tables/group_departures.csv")
  limits = pd.read_csv(root / "tables/exploratory_limits.csv")
  ratios = pd.read_csv(root / "tables/volume_normalized_ratios.csv")
  sens = pd.read_csv(root / "tables/d9_sensitivity.csv")
  figures = []

  def finish(fig: plt.Figure, slug: str, title: str, caption: str, note: str) -> None:
    """Save the same readable scientific figure in raster and editable vector formats."""
    fig.suptitle(title, fontsize=17, fontweight="bold", x=.04, ha="left", y=.985)
    fig.text(.04, .015, note, fontsize=9, color="#485462")
    fig.get_layout_engine().set(rect=(0, .045, 1, .90))
    fig.savefig(out / f"{slug}.png", dpi=220)
    fig.savefig(out / f"{slug}.svg")
    figures.append({"slug": slug, "title": title, "caption": caption,
                    "png": f"figures/{slug}.png", "svg": f"figures/{slug}.svg"})
    plt.close(fig)

  def condition(ax: plt.Axes, subset: pd.DataFrame, summary: pd.DataFrame, xcol: str, color: str) -> None:
    """Show all individual wells and mean plus sample SD at their actual numerical x values."""
    ax.scatter(subset[xcol], subset.corrected_a260, s=25, alpha=.60, c=color, zorder=3)
    ax.errorbar(summary[xcol], summary.mean_corrected, yerr=summary.sample_sd,
                fmt="D", ms=5, color="#202C38", capsize=4, lw=1.1, zorder=4)
    ax.axhline(0, color="#87929A", lw=.8)

  # 01: one full-sized, annotated plate map.
  fig, ax = plt.subplots(figsize=(13.5, 7.7), layout="constrained")
  dd = d.copy(); dd["volume_index"] = dd.volume_uL.map({v: i for i, v in enumerate(V)})
  ax.imshow(plate_array(dd, "volume_index"), cmap=ListedColormap([COLORS[v] for v in V]), vmin=-.5, vmax=3.5, alpha=.20, aspect="auto")
  for r in d.itertuples():
    label = f"{r.source_tube} | R{r.replicate}\n{r.concentration_ng_uL:g} ng/µL\n{r.volume_uL:g} µL"
    ax.text(int(r.well[1:])-1, "ABCDEFGH".index(r.well[0]), label, ha="center", va="center", fontsize=8.3)
  axes_plate(ax)
  finish(fig, "01_plate_layout", "01  Plate layout and electronic identity audit",
         "All 96 raw sample IDs match both the run manifest and the experiment layout. There are 72 DNA wells and 24 A31 buffer blanks. Original concentration metadata uses 16 ng/µL for A30; this map uses the user-corrected 70 ng/µL stock and nominal twofold dilutions. Electronic agreement does not independently prove physical contents.",
         "Each cell: source tube, nominal concentration, planned volume and dispensing replicate. A31 = buffer blank.")

  fig, axes = plt.subplots(1, 2, figsize=(13, 6.4), layout="constrained")
  for ax, field, title in zip(axes, ["a260", "corrected_a260"], ["Raw A260", "Volume-matched blank subtraction"]):
    a = plate_array(d, field)
    im = ax.imshow(a, cmap="viridis", vmin=0, vmax=.54, aspect="auto")
    for row in range(8):
      for col in range(12):
        ax.text(col, row, f"{a[row,col]:.3f}", color="white" if a[row,col] < .30 else "#17202B", ha="center", va="center", fontsize=7.3)
    axes_plate(ax); ax.set_title(title)
    ax.add_patch(Rectangle((7.5, 2.5), 1, 1, fill=False, lw=2, edgecolor="#E75B35"))
    fig.colorbar(im, ax=ax, shrink=.8, label="Absorbance (OD)")
  finish(fig, "02_plate_absorbance", "02  Raw and corrected readings retain every well",
         "All 96 readings are finite. Raw A260 ranges from 0.05649 to 0.53148 OD. Each well is corrected with its own volume's blank mean. D9 is outlined for follow-up, not excluded. Neither heatmap measures actual fill volumes.",
         "Color scales are identical. Labels are rounded for display; CSV tables retain calculation precision.")

  fig, axes = plt.subplots(1, 2, figsize=(12.5, 6.5), layout="constrained")
  for i, v in enumerate(V):
    t = d[(d.volume_uL == v) & (d.concentration_ng_uL == 0)]
    axes[0].scatter(np.full(6, v), t.a260, color=COLORS[v], s=40, alpha=.7)
    axes[0].errorbar(v, b.loc[v,"blank_mean"], yerr=b.loc[v,"blank_sd"], fmt="D", color="black", capsize=5)
    axes[1].vlines(v, b.loc[v,"sd_boot_ci95_low"], b.loc[v,"sd_boot_ci95_high"], color=COLORS[v], lw=2)
    axes[1].scatter(v, b.loc[v,"blank_sd"], s=80, color=COLORS[v])
  for ax in axes:
    ax.set_xticks(V); ax.set_xlabel("Planned volume (µL)")
  axes[0].set_title("Six blanks per volume: mean ± sample SD"); axes[0].set_ylabel("Raw A260 (OD)")
  axes[1].set_title("Blank SD with conditional bootstrap 95% interval"); axes[1].set_ylabel("Between-well blank SD (OD)"); axes[1].set_ylim(bottom=0)
  finish(fig, "03_blank_performance", "03  The 50 µL DNA-run blanks are quiet",
         "Blank SDs are 0.000330, 0.001054, 0.001046 and 0.002242 OD at 50, 75, 100 and 200 µL. The high 50 µL background seen in the earlier Orange G run is absent here. This is a different preparation and wavelength; it does not identify the cause of the earlier problem.",
         "These are six wells on one plate, not six independent preparations or repeat instrument scans.")

  fig, axes = plt.subplots(2, 2, figsize=(12.5, 8), layout="constrained")
  for ax, v in zip(axes.flat, V):
    t, q = s[s.volume_uL == v], g[g.volume_uL == v]
    f = fits[(fits.fit_kind == "concentration_full") & (fits.volume_uL == v)].iloc[0]
    condition(ax, t, q, "concentration_ng_uL", COLORS[v])
    x = np.linspace(0, C[-1], 100); ax.plot(x, f.intercept+f.slope*x, color=COLORS[v], lw=1.5)
    ax.set(xlim=(0,72.5), xticks=[0,10,20,35,50,70], xlabel="Nominal DNA concentration (ng/µL)", ylabel="Blank-corrected A260 (OD)", title=f"{v} µL")
    ax.text(.04,.96,f"slope {f.slope:.5f}  |  intercept {f.intercept:.4f}\nR² {f.r_squared:.5f}  |  RMSE {f.residual_rmse:.4f}", transform=ax.transAxes, va="top", fontsize=9, bbox=dict(facecolor="white", alpha=.85, edgecolor="none"))
  finish(fig, "04_concentration_series", "04  Concentration-response curves at each volume",
         "Responses rise strongly with nominal concentration at every volume. Full-range R² values are high, but this alone does not establish strict linearity or low-end accuracy. Fit all 18 sample wells per panel with a free intercept; the residual analysis shows where a straight line misses group means.",
         "Dots = all wells; diamonds and bars = group mean ± sample SD (n=3). Lines = ordinary least squares, free intercept.")

  fig, axes = plt.subplots(2,2,figsize=(12.5,8),layout="constrained")
  for ax,v in zip(axes.flat,V):
    t=s[(s.volume_uL==v)&(s.concentration_ng_uL<=C[2])];q=g[(g.volume_uL==v)&(g.concentration_ng_uL<=C[2])]
    condition(ax,t,q,"concentration_ng_uL",COLORS[v])
    bb=d[(d.volume_uL==v)&(d.concentration_ng_uL==0)]
    ax.scatter(np.zeros(6),bb.corrected_a260,c="#7A8792",s=24)
    ax.axhline(3.3*b.loc[v,"blank_sd"],ls="--",color="#9A6224",label="3.3 × blank SD")
    ax.axhline(10*b.loc[v,"blank_sd"],ls=":",color="#6F587A",label="10 × blank SD")
    ax.set(title=f"{v} µL",xlabel="Nominal DNA concentration (ng/µL)",ylabel="Blank-corrected A260 (OD)",xticks=[0,C[0],C[1],C[2]],xlim=(-.1*SCALE,2.12*SCALE))
    ax.legend(loc="upper left",fontsize=8)
    if v==50:
      row=t[t.well=="D9"].iloc[0];ax.annotate("D9",(C[0],row.corrected_a260),xytext=(.8*SCALE,row.corrected_a260+.003),arrowprops=dict(arrowstyle="-",color="black"),fontsize=10)
  finish(fig,"05_low_concentration_zoom","05  Detecting a signal is different from precise measurement",
         "All three 2.1875 ng/µL wells at each volume exceed the matching 3.3-blank-SD screen. Their corrected CVs are 100.5%, 18.9%, 8.5% and 8.0% as volume increases. Thus a favorable blank-noise ratio at 50 µL does not rescue its inconsistent low-concentration samples. Thresholds are exploratory screens, not validated detection probabilities.",
         "Gray points = volume-matched blanks after subtraction. Thresholds describe single-well signal, not standard errors of triplicate means.")

  fig,axes=plt.subplots(2,3,figsize=(13,8),layout="constrained")
  for ax,c in zip(axes.flat,C):
    t=s[s.concentration_ng_uL==c];q=g[g.concentration_ng_uL==c]
    condition(ax,t,q,"volume_uL","#488697")
    f=fits[(fits.fit_kind=="volume")&(fits.concentration_ng_uL==c)].iloc[0]
    x=np.linspace(50,200,100);ax.plot(x,f.intercept+f.slope*x,color="#4B626C",lw=1)
    ax.set(title=f"{c:g} ng/µL  |  R² {f.r_squared:.3f}",xlabel="Planned volume (µL)",ylabel="Corrected A260 (OD)",xticks=V)
  finish(fig,"06_volume_series","06  Signal generally increases with well volume",
         "At fixed concentration, larger planned volumes generally give larger corrected absorbance. The lowest concentration is disrupted by the elevated D9 result at 50 µL. At 2.1875 ng/µL, moving from 100 to 200 µL nearly doubles mean signal, but blank variability also increases, so signal gain alone is not a detection-limit improvement.",
         "Numerical linear volume axes. Dots = wells; diamonds = mean ± sample SD; free-intercept fits use 12 wells per panel.")

  fig,ax=plt.subplots(figsize=(10.5,6.5),layout="constrained")
  values=g.pivot(index="concentration_ng_uL",columns="volume_uL",values="cv_percent").loc[C,V].to_numpy()
  im=ax.imshow(values,cmap="YlOrRd",vmin=0,vmax=20,aspect="auto")
  for ci,c in enumerate(C):
    for vi,v in enumerate(V):
      ax.text(vi,ci,f"{values[ci,vi]:.1f}%",ha="center",va="center",fontsize=13,color="white" if values[ci,vi]>13 else "#14232F")
  ax.set_xticks(range(4),[f"{v} µL" for v in V]);ax.set_yticks(range(6),[f"{c:g}" for c in C])
  ax.set(xlabel="Planned well volume",ylabel="Nominal DNA concentration (ng/µL)");ax.grid(False)
  fig.colorbar(im,ax=ax,label="CV (%) of blank-corrected triplicates; color capped at 20%")
  finish(fig,"07_replicate_precision","07  Precision favors 100 µL at the lowest tested concentration",
         "At 100 and 200 µL, all six concentrations have within-condition corrected CV below 10%. At 75 µL, 2.1875 and 4.375 ng/µL have CVs of 18.9% and 12.9%. At 50 µL, CV is below 10% from 8.75 to 70 ng/µL, but fails badly at 2.1875 ng/µL. These three-well CVs measure within-plate repeatability, not dilution accuracy or between-day precision.",
         "CV = 100 × sample SD / corrected mean. The common blank mean affects the denominator; three replicates give uncertain precision estimates.")

  fig,axes=plt.subplots(1,2,figsize=(12.8,6.6),layout="constrained")
  values=g.pivot(index="concentration_ng_uL",columns="volume_uL",values="mean_in_blank_sd").loc[C,V].to_numpy()
  im=axes[0].imshow(values,cmap="viridis",norm=LogNorm(vmin=3.3,vmax=400),aspect="auto")
  for ci,c in enumerate(C):
    for vi,v in enumerate(V):
      row=g[(g.volume_uL==v)&(g.concentration_ng_uL==c)].iloc[0]
      axes[0].text(vi,ci,f"{values[ci,vi]:.1f}",ha="center",va="center",color="white" if values[ci,vi]<100 else "black",fontsize=11)
      axes[1].text(vi,ci,f"{row.replicates_above_10sd}/3\nCV {row.cv_percent:.1f}%",ha="center",va="center",fontsize=10)
  axes[1].imshow(np.zeros((6,4)),cmap=ListedColormap(["#F0F4F6"]),aspect="auto")
  # Put text above the image explicitly.
  for t in axes[1].texts:t.set_zorder(3)
  for ax in axes:
    ax.set_xticks(range(4),V);ax.set_yticks(range(6),[f"{c:g}" for c in C]);ax.set_xlabel("Planned volume (µL)");ax.set_ylabel("Nominal concentration (ng/µL)");ax.grid(False)
  axes[0].set_title("Mean signal / matched blank SD");axes[1].set_title("Wells above 10 × blank SD, with CV")
  fig.colorbar(im,ax=axes[0],label="Signal / blank SD (log color scale)",shrink=.75)
  finish(fig,"08_signal_and_precision","08  Blank noise and sample precision must be considered together",
         "The 50 µL, 2.1875 ng/µL mean is 31.7 blank SD above zero, yet its CV is 100.5%. In contrast, 100 µL gives a smaller ratio of 8.24 but much tighter samples. The blank-noise ratio is not instrument signal-to-noise and is not a substitute for repeatable sample results or accuracy validation.",
         "Left: mean corrected A260 / six-blank sample SD. Right: individual wells above the illustrative 10-SD screen.")

  fig,axes=plt.subplots(1,3,figsize=(13.5,6.4),layout="constrained")
  for vi,v in enumerate(V):
    for kind,marker,offset in (("concentration_full","o",-2),("concentration_low","^",2)):
      f=fits[(fits.fit_kind==kind)&(fits.volume_uL==v)].iloc[0]
      axes[0].vlines(v+offset,f.slope_ci95_low,f.slope_ci95_high,color=COLORS[v],lw=1)
      axes[0].scatter(v+offset,f.slope,color=COLORS[v],marker=marker,s=45)
      for ax,metric in zip(axes[1:],["DL","QL"]):
        r=limits[(limits.volume_uL==v)&(limits.fit_kind==kind)&(limits.metric==metric)].iloc[0]
        ax.vlines(v+offset,r.ci95_low,r.ci95_high,color=COLORS[v],lw=1.4)
        ax.scatter(v+offset,r.estimate_ng_uL,c=COLORS[v],marker=marker,s=45)
  axes[0].set(title="Calibration slope",ylabel="OD per ng/µL")
  for ax,metric in zip(axes[1:],["DL: 3.3 × blank SD / slope","QL: 10 × blank SD / slope"]):
    ax.set(title=metric,ylabel="Estimated concentration (ng/µL)")
    ax.axhline(C[0],color="#8B5D39",ls="--",lw=1,label="Lowest tested: 2.1875")
    ax.legend(fontsize=8)
  for ax in axes:ax.set_xlabel("Planned volume (µL)");ax.set_xticks(V);ax.set_ylim(bottom=0)
  finish(fig,"09_exploratory_limits","09  Noise-based limits are hypotheses to validate",
         "Full-range blank-noise DL estimates are approximately 0.73, 1.44, 1.02 and 1.12 ng/µL; QL estimates are 2.22, 4.35, 3.09 and 3.40 ng/µL. The apparently favorable 50 µL DL conflicts with its low-end sample variability. The fixed 2.1875-17.5 ng/µL fits show range sensitivity. Values below 2.1875 ng/µL were not tested, and neither these intervals nor high R² establish a validated limit.",
         "Circles = 2.1875-70 ng/µL fit; triangles = fixed 2.1875-17.5 fit. Bars = conditional 95% bootstrap intervals (10,000 resamples).")

  fig,axes=plt.subplots(2,3,figsize=(13,8),layout="constrained")
  for ax,c in zip(axes.flat,C):
    q=ratios[ratios.concentration_ng_uL==c]
    for r in q.itertuples():
      ax.vlines(r.volume_uL,r.ci95_low,r.ci95_high,color=COLORS[r.volume_uL],lw=1.3)
      ax.scatter(r.volume_uL,r.relative_volume_normalized_signal,c=COLORS[r.volume_uL],s=50)
    ax.axhline(1,ls="--",c="#6A7880",lw=1)
    ax.set(title=f"{c:g} ng/µL",xticks=V,xlabel="Planned volume (µL)",ylabel="(OD / volume) relative to 100 µL")
  finish(fig,"10_volume_normalization","10  Volume normalization diagnoses departures from simple scaling",
         "Most volume-normalized responses are close to the 100 µL reference, particularly at moderate concentrations. The 50 µL, 2.1875 ng/µL point is elevated and uncertain because D9 raises its mean. Dividing by planned volume is a descriptive normalization, not a measured path-length correction or evidence of actual dispensed volume.",
         "Points = ratios of group means; bars = bootstrap 95% intervals. The 100 µL reference equals one by definition.")

  fig,axes=plt.subplots(2,2,figsize=(12.5,8),layout="constrained")
  for ax,v in zip(axes.flat,V):
    q=dep[(dep.fit_kind=="concentration_full")&(dep.volume_uL==v)]
    t=s[s.volume_uL==v]
    ax.scatter(t.concentration_ng_uL,t.model_residual,c=COLORS[v],s=22,alpha=.45)
    for r in q.itertuples():
      color="#B54F37" if r.excludes_zero else "#273A48"
      ax.vlines(r.concentration_ng_uL,r.ci95_low,r.ci95_high,color=color,lw=2)
      ax.scatter(r.concentration_ng_uL,r.departure,c=color,marker="D",s=40)
    ax.axhline(0,c="#68757D",lw=1)
    ax.set(title=f"{v} µL",xlabel="Nominal concentration (ng/µL)",ylabel="Observed minus fitted A260 (OD)",xticks=[C[0],C[2],C[3],C[4],C[5]])
  finish(fig,"11_concentration_residuals","11  High R² does not rule out systematic departures",
         "Red group means have 95% bootstrap departure intervals that exclude zero. Several conditions depart detectably from their fitted lines despite high R². These are exploratory, unadjusted comparisons based on one plate. They can reflect optical response, dilution error, dispensing differences or isolated wells; they do not isolate a mechanism.",
         "Dots = individual residuals; diamonds and vertical lines = group departures with bootstrap 95% intervals. No points were removed.")

  fig,axes=plt.subplots(1,2,figsize=(13,6.4),layout="constrained")
  for ax,field,title in zip(axes,["model_residual","within_group_residual"],["Residual from concentration fit","Residual from condition mean"]):
    a=plate_array(d,field);lim=max(abs(a.min()),abs(a.max()))
    im=ax.imshow(a,cmap="RdBu_r",norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim),aspect="auto")
    axes_plate(ax);ax.set_title(title)
    for _,r in d.reindex(d[field].abs().sort_values(ascending=False).index).head(4).iterrows():
      ax.text(r.plate_column-1,"ABCDEFGH".index(r.plate_row),r.well,ha="center",va="center",fontsize=9,color="white",fontweight="bold")
    fig.colorbar(im,ax=ax,label="Residual (OD)",shrink=.8)
  finish(fig,"12_spatial_residuals","12  Spatial checks locate influential wells without assigning causes",
         "D9 is conspicuous within its own concentration-volume condition. Spatial residuals also identify the strongest discrepancies elsewhere, but each location has only one observation and variance changes with signal. A single map cannot distinguish a positional optical effect from pipetting, dilution or a bubble. All wells remain included.",
         "Left: blanks centered on their own volume mean. Right: every well centered on its concentration-volume group mean.")

  fig,axes=plt.subplots(1,2,figsize=(12.5,6.5),layout="constrained")
  low=s[(s.volume_uL==50)&(s.concentration_ng_uL==C[0])].sort_values("well")
  for i,r in enumerate(low.itertuples()):
    axes[0].scatter(i,r.corrected_a260,c="#BA5035" if r.well=="D9" else COLORS[50],s=85)
    axes[0].text(i,r.corrected_a260+.0007,f"{r.corrected_a260:.4f}",ha="center",fontsize=10)
  baseline=g[(g.volume_uL==100)&(g.concentration_ng_uL==C[0])].mean_corrected.iloc[0]/2
  axes[0].axhline(baseline,color=COLORS[100],ls="--",label="Half of 100 µL mean (scaling reference)")
  axes[0].set(xticks=range(3),xticklabels=low.well.tolist(),ylabel="Corrected A260 (OD)",title="50 µL at 2.1875 ng/µL",ylim=(0,.027))
  axes[0].legend(loc="upper left",fontsize=8)
  t=s[(s.volume_uL==50)&(s.concentration_ng_uL<=C[3])];q=g[(g.volume_uL==50)&(g.concentration_ng_uL<=C[3])]
  condition(axes[1],t,q,"concentration_ng_uL",COLORS[50])
  x=np.linspace(0,C[3],100)
  for r,ls,color in zip(sens.itertuples(),["-","--"],["#C16E36","#557688"]):
    axes[1].plot(x,r.intercept+r.slope*x,ls=ls,c=color,label=r.analysis)
  axes[1].set(title="Influence on the full-range fitted line (zoom)",xlabel="Nominal concentration (ng/µL)",ylabel="Corrected A260 (OD)",xticks=[0,C[0],C[1],C[2],C[3]])
  axes[1].legend(fontsize=8)
  finish(fig,"13_d9_sensitivity","13  D9 drives the largest low-volume inconsistency",
         "At 50 µL and 2.1875 ng/µL, A8, C10 and D9 read 0.005148, 0.003658 and 0.022558 OD after subtraction. The three-well CV is 100.5%. Omitting D9 only as a sensitivity calculation leaves two wells with CV 23.9%, still uncertain. D9 has no documented physical exclusion reason and is retained in every primary result.",
         "The secondary omission is diagnostic only. No raw value or original validity field was changed.")

  fig,axes=plt.subplots(1,2,figsize=(12.5,6.5),layout="constrained")
  for ax,maximum,title in zip(axes,[3300*SCALE,420*SCALE],["All tested DNA masses","Low-mass detail"]):
    for v in V:
      q=g[(g.volume_uL==v)&(g.mass_ng<=maximum)]
      t=s[(s.volume_uL==v)&(s.mass_ng<=maximum)]
      ax.scatter(t.mass_ng,t.corrected_a260,c=COLORS[v],s=20,alpha=.45)
      ax.errorbar(q.mass_ng,q.mean_corrected,yerr=q.sample_sd,fmt="o",ms=4,c=COLORS[v],capsize=3,label=f"{v} µL")
    ax.set(xlim=(0,maximum),xlabel="Nominal DNA mass per well (ng)",ylabel="Corrected A260 (OD)",title=title)
    ax.legend(fontsize=9)
  finish(fig,"14_dna_mass","14  Minimum concentration and minimum sample mass are different goals",
         "At the lowest tested concentration of 2.1875 ng/µL, the wells contain nominally 109.375, 164.063, 218.75 and 437.5 ng at increasing volumes. The 100 µL condition resolves that level with 8.5% within-plate CV using 218.75 ng. Reducing volume saves DNA but reduces optical signal and can expose dispensing or meniscus variability. Actual DNA mass was not independently measured.",
         "Mass = nominal concentration × planned volume. Points = wells; bars = sample SD, n=3. Axes are numerical and linear.")
  fig, axes = plt.subplots(2, 3, figsize=(13, 8), layout="constrained")
  for ax, c in zip(axes.flat, C):
    q = dep[(dep.fit_kind == "volume") & (dep.concentration_ng_uL == c)]
    f = fits[(fits.fit_kind == "volume") & (fits.concentration_ng_uL == c)].iloc[0]
    t = s[s.concentration_ng_uL == c]
    ax.scatter(t.volume_uL, t.corrected_a260 - f.intercept - f.slope*t.volume_uL,
               color="#6D8790", alpha=.5, s=22)
    for r in q.itertuples():
      color = "#B54F37" if r.excludes_zero else "#273A48"
      ax.vlines(r.volume_uL, r.ci95_low, r.ci95_high, color=color, lw=2)
      ax.scatter(r.volume_uL, r.departure, color=color, marker="D", s=40)
    ax.axhline(0, color="#68757D", lw=1)
    ax.set(title=f"{c:g} ng/µL", xlabel="Planned volume (µL)",
           ylabel="Observed minus fitted A260 (OD)", xticks=V)
  finish(fig, "15_volume_residuals", "15  Volume-response residuals expose departures from a line",
         "The 35 ng/µL series has an exploratory departure at 200 µL; the 70 ng/µL series has departures at 50, 100 and 200 µL. Lower-concentration series show no clear group departure by this conditional interval rule. The lowest concentration remains imprecise because of D9; absence of a clear departure is not proof of linearity.",
         "Diamonds and bars = group departures and conditional 95% bootstrap intervals. Red intervals exclude zero; no multiplicity adjustment.")
  (root / "figure_index.json").write_text(json.dumps(figures,indent=2)+"\n")
  return figures


if __name__ == "__main__":
  print(f"Exported {len(make_figures())} PNG/SVG figure pairs.")
