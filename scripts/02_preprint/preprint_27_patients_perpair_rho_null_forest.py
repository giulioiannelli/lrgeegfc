#!/usr/bin/env python3
"""Per-pair cophenetic trace, grouped BY PATIENT across bands.

The patient-grouped twin of ``preprint_23`` (which is grouped by band, patients
as rows). Here each panel is one patient and the six bands are the forest rows,
so you read a single patient's trace across all bands at a glance — and, by
scanning the same row across panels, how coherent each band is across the
cohort. Same locked gate data, same encoding as preprint_23.

Per patient (2×5 grid, canonical cohort order), one row per band
(δ θ α β γ_l γ_h, δ at top):

  • light grey span  = the patient's matched-strength surrogate null, p5–p95
  • tick at centre   = surrogate median ρ (null centre)
  • coloured stem    = obs ρ − null-median (the paired difference the Wilcoxon
                       gate ranks); GREEN = trace direction (obs > null median),
                       RED = anti
  • dot at obs ρ     = filled if it clears its OWN null (obs > p95, or < p5);
                       open if it sits inside the null band

The band labels (left column) are coloured by the locked cohort gate verdict
(green = clears the matched-strength gate at p<0.05 → α, β; grey = does not).
A cohort trace band reads as a green stem in (almost) every patient panel;
a split band (γ_l) reads as green stems in most panels but a red cluster in a
few; a non-trace band straddles the null both ways across panels.

Gate p, coherence counts and the per-patient×band table are printed to stdout
(kept out of the figure per the no-in-axes-text rule); they belong in the
caption.

Input (locked, instant — no recompute): the §5.3 split-baseline matched-strength
per-patient table
``data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv``
(cophenetic substrate). Output:
``data/preprint/figures/all_bands/patients_forest_plot/
fig_patients_perpair_rho_null_forest_coph.pdf``
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
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
    "ytick.labelsize": "medium",
    "legend.fontsize": "x-small",   # 5 long entries on one row → keep compact
})

PERPAT_CSV = (ROOT / "data" / "audit" / "matched_strength_surrogate_split_baseline"
             / "per_patient_per_band.csv")
# δ at top → γ_h at bottom (axis is inverted after plotting)
BAND_ORDER = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]

C_PRO = "#1a7c3e"      # obs ρ above null median — trace direction
C_ANTI = "#c0392b"     # obs ρ below null median — anti
C_NULL = "#b9b9b9"     # matched-strength null span (p5–p95)
C_NULLMID = "#7a7a7a"  # null median tick


def gate_p(d: pd.DataFrame) -> float:
    """Locked §5.3 gate: paired Wilcoxon of (obs ρ − own surrogate median), >."""
    diff = (d.obs_rho - d.surr_p50).values
    return float(wilcoxon(diff, alternative="greater").pvalue)


def plot_patient(ax, d: pd.DataFrame, xlim) -> None:
    """One patient: a forest row per band, fixed band order (δ top → γ_h bottom)."""
    d = d.set_index("band")
    for j, band in enumerate(BAND_ORDER):
        row = d.loc[band]
        ax.plot([row.surr_p5, row.surr_p95], [j, j], color=C_NULL, lw=5.5,
                solid_capstyle="round", alpha=0.6, zorder=1)
        ax.plot([row.surr_p50], [j], marker="|", ms=8, color=C_NULLMID,
                mew=1.3, zorder=2)
        pro = row.obs_rho > row.surr_p50
        col = C_PRO if pro else C_ANTI
        ax.plot([row.surr_p50, row.obs_rho], [j, j], color=col, lw=1.8, zorder=3)
        beyond = (row.obs_rho > row.surr_p95) or (row.obs_rho < row.surr_p5)
        ax.plot([row.obs_rho], [j], marker="o", ms=7.5, zorder=4,
                mfc=(col if beyond else "white"), mec=col, mew=1.5)

    ax.axvline(0.0, color="0.45", ls="--", lw=0.8, zorder=0)
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.6, len(BAND_ORDER) - 0.4)
    ax.invert_yaxis()                       # δ (j=0) at the top
    ax.set_yticks(range(len(BAND_ORDER)))
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines["left"].set_color("0.6")
    ax.spines["bottom"].set_color("0.6")


def main() -> Path:
    df = pd.read_csv(PERPAT_CSV)
    patients = list(dict.fromkeys(df.patient))           # canonical cohort order
    lo = float(min(df.obs_rho.min(), df.surr_p5.min()))
    hi = float(max(df.obs_rho.max(), df.surr_p95.max()))
    pad = 0.06 * (hi - lo)
    xlim = (lo - pad, hi + pad)

    # locked cohort gate verdict per band → colours the band labels
    verdict = {b: gate_p(df[df.band == b]) for b in BAND_ORDER}
    band_label_col = {b: (C_PRO if verdict[b] < 0.05 else "0.45") for b in BAND_ORDER}
    band_tex = [BRAIN_BAND_TEX_DICT.get(b, rf"${b}$") for b in BAND_ORDER]

    nrow, ncol = 2, 5
    fig = plt.figure(figsize=(16.5, 8.0))
    outer = fig.add_gridspec(nrow, ncol, wspace=0.12, hspace=0.22,
                             left=0.045, right=0.99, top=0.95, bottom=0.13)

    print("Per-pair cophenetic trace — per patient across bands "
          "(locked gate verdict colours band labels):")
    print(f"  cohort gate p: " + "  ".join(
        f"{b}={verdict[b]:.4f}{'*' if verdict[b] < 0.05 else ''}"
        for b in BAND_ORDER))
    print(f"  {'patient':9s} " + " ".join(f"{b[:5]:>6s}" for b in BAND_ORDER)
          + "   n_pro n_clear")

    for idx, pat in enumerate(patients):
        rr, cc = idx // ncol, idx % ncol
        ax = fig.add_subplot(outer[rr, cc])
        d = df[df.patient == pat]
        plot_patient(ax, d, xlim)

        if cc == 0:                          # band labels only on the left column
            ax.set_yticklabels(band_tex, size="large")
            for tick, b in zip(ax.get_yticklabels(), BAND_ORDER):
                tick.set_color(band_label_col[b])
        else:
            ax.set_yticklabels([])
        if rr == nrow - 1:
            ax.set_xlabel(r"$\rho^{\mathrm{coph}}$", color="0.20")

        # patient tag (per-panel identifier, analogue of the band symbol in p23)
        ax.text(0.97, 0.04, pat.replace("Pat_", "Pat "), transform=ax.transAxes,
                ha="right", va="bottom", size="large", fontweight="bold",
                color="0.18")

        ds = d.set_index("band")
        n_pro = int((ds.obs_rho > ds.surr_p50).sum())
        n_clear = int((ds.obs_rho > ds.surr_p95).sum())
        print(f"  {pat:9s} " + " ".join(
            f"{ds.loc[b].obs_rho:+6.2f}" for b in BAND_ORDER)
            + f"   {n_pro:5d} {n_clear:7d}")

    handles = [
        Line2D([0], [0], color=C_NULL, lw=5.5, solid_capstyle="round",
               alpha=0.6, label="matched-strength null (p5–p95)"),
        Line2D([0], [0], marker="o", color=C_PRO, lw=1.8, mfc=C_PRO, mec=C_PRO,
               ms=7.5, label="trace direction, clears own null (filled)"),
        Line2D([0], [0], marker="o", color=C_PRO, lw=1.8, mfc="white", mec=C_PRO,
               ms=7.5, label="trace direction, within null (open)"),
        Line2D([0], [0], marker="o", color=C_ANTI, lw=1.8, mfc="white",
               mec=C_ANTI, ms=7.5, label="anti direction"),
        Line2D([0], [0], marker="s", color="none", mfc=C_PRO, mec="none",
               ms=11, label="green band label = clears cohort gate (p<0.05)"),
    ]
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, 0.005), ncol=5, frameon=False)

    out_dir = (ROOT / "data" / "preprint" / "figures" / "all_bands"
               / "patients_forest_plot")
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_patients_perpair_rho_null_forest_coph.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
