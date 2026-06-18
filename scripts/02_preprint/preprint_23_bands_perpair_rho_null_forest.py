#!/usr/bin/env python3
"""Per-pair cophenetic trace, shown the way the GATE actually tests it.

The joint-density renders (preprint_17/18) fail on the cophenetic substrate
because ``Δ^coph`` is ~90% tied: a 2D density piles the non-mover atom on the
rank-centre (= the diagonal) and fakes a trace in every band. Spearman ρ — the
§5.3 statistic of record — is tie-corrected, so it ignores that centre and is
carried by the pairs that actually move. This figure renders ρ directly, so the
ties never enter.

Per band (2×3, δ θ α / β γ_l γ_h), one row per patient:

  • light grey span  = the patient's matched-strength surrogate null, p5–p95
  • tick at centre   = surrogate median ρ (null centre)
  • coloured stem    = obs ρ − null-median (the paired difference the Wilcoxon
                       gate ranks); GREEN = trace direction (obs > null median),
                       RED = anti
  • dot at obs ρ     = filled if it clears its OWN null (obs > p95, or < p5 for
                       a strong anti); open if it sits inside the null band

A cohort trace = most stems green and most dots filled-and-right-of-grey
(β, α). A split cohort = green stems at the top but a red cluster at the bottom
(γ_l: Pat_05/02/06 strongly pro, Pat_14/15/07/13 anti → fails the gate). A
non-trace = stems straddling the null both ways (δ, θ, γ_h).

Row ordering is selectable (``--sort``):
  • ``rho``     — patients sorted by observed ρ within each panel (default; the
                  classic forest, best at top). Rows are anonymous.
  • ``patient`` — patients in fixed canonical order (Pat_02 … Pat_15, Pat_02 at
                  top) identical across all six panels, so the SAME patient sits
                  at the SAME row in every band — read down a column of panels to
                  follow one patient across bands. Patient labels on the left.
The sort type is appended to the output filename.

The band symbol is coloured by the locked gate verdict (green = clears the
matched-strength gate at p<0.05, grey = does not). Gate p, coherence counts and
the per-patient table are printed to stdout (kept out of the figure per the
no-in-axes-text rule); they belong in the caption.

Input (locked, instant — no recompute): the §5.3 split-baseline matched-strength
per-patient table
``data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv``
(cophenetic substrate). Output:
``data/preprint/figures/all_bands/patients_forest_plot/
fig_bands_perpair_rho_null_forest_coph_sort-{rho,patient}.pdf``
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import wilcoxon

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

# Single typography knob — every text size below is set RELATIVE to this (via
# matplotlib's named sizes), so bumping BASE_FONTSIZE rescales the whole
# figure's text at once. No per-element fontsize is hardcoded.
BASE_FONTSIZE = 16
plt.rcParams.update({
    "font.size": BASE_FONTSIZE,
    "axes.labelsize": "large",      # axis labels
    "xtick.labelsize": "medium",    # ρ-axis numbers
    "ytick.labelsize": "medium",    # patient row labels (sort=patient)
    "legend.fontsize": "large",
})

PERPAT_CSV = (ROOT / "data" / "audit" / "matched_strength_surrogate_split_baseline"
             / "per_patient_per_band.csv")
BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
# canonical cohort order (n=10), used when --sort patient (Pat_02 at top)
CANON = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
         "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]

C_PRO = "#1a7c3e"      # obs ρ above null median — trace direction
C_ANTI = "#c0392b"     # obs ρ below null median — anti
C_NULL = "#b9b9b9"     # matched-strength null span (p5–p95)
C_NULLMID = "#7a7a7a"  # null median tick


def gate_p(d: pd.DataFrame) -> float:
    """Locked §5.3 gate: paired Wilcoxon of (obs ρ − own surrogate median), >."""
    diff = (d.obs_rho - d.surr_p50).values
    return float(wilcoxon(diff, alternative="greater").pvalue)


def plot_band(ax, d: pd.DataFrame, xlim, p_gate: float, sort: str,
              show_ylabels: bool) -> None:
    if sort == "patient":
        order = [p for p in CANON if p in set(d.patient)]      # top → bottom
        d = d.set_index("patient").loc[order].reset_index()
    else:                                                      # sort == "rho"
        order = None
        d = d.sort_values("obs_rho").reset_index(drop=True)

    for i, row in d.iterrows():
        # matched-strength null band (p5–p95) + median tick
        ax.plot([row.surr_p5, row.surr_p95], [i, i], color=C_NULL, lw=5.5,
                solid_capstyle="round", alpha=0.6, zorder=1)
        ax.plot([row.surr_p50], [i], marker="|", ms=8, color=C_NULLMID,
                mew=1.3, zorder=2)
        # paired-difference stem from null centre to observed ρ
        pro = row.obs_rho > row.surr_p50
        col = C_PRO if pro else C_ANTI
        ax.plot([row.surr_p50, row.obs_rho], [i, i], color=col, lw=1.8, zorder=3)
        beyond = (row.obs_rho > row.surr_p95) or (row.obs_rho < row.surr_p5)
        ax.plot([row.obs_rho], [i], marker="o", ms=7.5, zorder=4,
                mfc=(col if beyond else "white"), mec=col, mew=1.5)

    ax.axvline(0.0, color="0.45", ls="--", lw=0.8, zorder=0)
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.8, len(d) - 0.2)
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["bottom"].set_color("0.6")

    if sort == "patient":
        ax.invert_yaxis()                       # Pat_02 (i=0) at the top
        ax.spines["left"].set_color("0.6")
        if show_ylabels:
            ax.set_yticks(range(len(d)))
            ax.set_yticklabels([p.replace("Pat_", "Pat ") for p in order],
                               size="medium")
        else:
            ax.set_yticks([])
    else:
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)

    # band symbol, coloured by the locked gate verdict
    band = d.band.iloc[0]
    tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")
    ax.text(0.97, 0.07, tex, transform=ax.transAxes, ha="right", va="bottom",
            size="xx-large", fontweight="bold",
            color=(C_PRO if p_gate < 0.05 else "0.45"))


def main() -> Path:
    ap = argparse.ArgumentParser(
        description="Per-pair cophenetic trace forest (band-grouped).")
    ap.add_argument("--sort", choices=["rho", "patient"], default="rho",
                    help="row order within each band panel: 'rho' (default, "
                         "forest sorted by observed ρ) or 'patient' (fixed "
                         "canonical order, same patient at the same row across "
                         "all bands).")
    args = ap.parse_args()

    df = pd.read_csv(PERPAT_CSV)
    lo = float(df.obs_rho.min()); hi = float(df.obs_rho.max())
    pad = 0.06 * (hi - lo)
    xlim = (lo - pad, hi + pad)

    left = 0.085 if args.sort == "patient" else 0.04
    fig = plt.figure(figsize=(15.5, 9.5))
    outer = fig.add_gridspec(2, 3, wspace=0.10, hspace=0.20,
                             left=left, right=0.985, top=0.93, bottom=0.20)
    print(f"Per-pair cophenetic trace — locked gate (paired Wilcoxon, "
          f"one-sided); sort={args.sort}:")
    print(f"  {'band':11s} {'gate_p':>8s} {'n_pro':>6s} {'n_anti':>6s} "
          f"{'n_clear_p95':>11s} {'verdict':>8s}")
    for idx, band in enumerate(BAND_ORDER):
        d = df[df.band == band]
        rr, cc = idx // 3, idx % 3
        ax = fig.add_subplot(outer[rr, cc])
        p = gate_p(d)
        plot_band(ax, d, xlim, p, args.sort, show_ylabels=(cc == 0))
        if rr == 1:
            ax.set_xlabel(r"$\rho^{\mathrm{coph}}$", color="0.20")
        n_pro = int((d.obs_rho > d.surr_p50).sum())
        n_anti = int((d.obs_rho < d.surr_p50).sum())
        n_clear = int((d.obs_rho > d.surr_p95).sum())
        verdict = "TRACE" if p < 0.05 else "no"
        print(f"  {band:11s} {p:8.4f} {n_pro:6d} {n_anti:6d} {n_clear:11d} "
              f"{verdict:>8s}")

    handles = [
        Line2D([0], [0], color=C_NULL, lw=9, solid_capstyle="round",
               alpha=0.6, label="matched-strength null (p5–p95)"),
        Line2D([0], [0], marker="o", color=C_PRO, lw=2.6, mfc=C_PRO, mec=C_PRO,
               ms=13, label="trace, clears null (filled)"),
        Line2D([0], [0], marker="o", color=C_PRO, lw=2.6, mfc="white", mec=C_PRO,
               ms=13, label="trace, within null (open)"),
        Line2D([0], [0], marker="o", color=C_ANTI, lw=2.6, mfc=C_ANTI,
               mec=C_ANTI, ms=13, label="anti"),
    ]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, 0.005), ncol=2, frameon=False,
               handlelength=2.4, handletextpad=0.7, columnspacing=2.5)

    out_dir = (ROOT / "data" / "preprint" / "figures" / "all_bands"
               / "patients_forest_plot")
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"fig_bands_perpair_rho_null_forest_coph_sort-{args.sort}.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
