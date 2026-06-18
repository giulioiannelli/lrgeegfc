#!/usr/bin/env python3
"""HEADLINE figure: the β per-pair cophenetic trace, as a clean joint density.

Why a smooth pooled density is honest HERE (where it was not for the all-bands
comparison): the cross-band magnitude trap was caused by bands with one or two
extreme-ρ patients (γ_l: Pat_05 0.85; γ_h: Pat_05 0.92) dominating a pooled
density. β has NO such outlier — its eight positive patients are all moderate
(ρ 0.21–0.51) — so pooling β's movers is representative of the cohort, not
driven by anyone. The smooth density is therefore a faithful picture of the
β trace, and it is catchy.

Three elements, one figure:
  • main panel — smooth KDE of the pooled β movers on within-patient rank
    coordinates u = rank Δ^coph_task, v = rank Δ^coph_rest (non-mover atom
    dropped, survivors re-ranked → tie-free). A trace is the green diagonal
    ridge (u≈v): pairs reorganized by the task stay reorganized into rest.
  • marginals — conditional split: P(u) for the upper vs lower half of v
    (top), P(v) for the upper vs lower half of u (right). A positive trace
    pushes the green (high-other-rank) histogram to the right of the red.
  • bottom strip — the statistics: each patient's ρ^coph_split as a dot against
    its own matched-strength null (p5–p95 grey bar). 8/10 right of zero, 7/10
    clear their own null → the density is a cohort-consistent effect, not a
    pooling artifact. (Gate: paired Wilcoxon p=0.005; ρ̄=0.22.)

Inputs (locked, instant): per-patient gate table
``matched_strength_surrogate_split_baseline/per_patient_per_band.csv`` for the
strip; cached cophenetic per-pair splits for the density.
Output: data/preprint/figures/all_bands/fig_beta_perpair_headline_density_coph.pdf
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
sys.path.insert(0, str(Path(__file__).resolve().parent))
p18 = importlib.import_module("preprint_18_bands_joint_density_rawfc")
use_lrg_style()

BAND = "beta"
COHORT = p18.COHORT
GRID = 200
KDE_BW = 0.16
SEED = 20260605
PERPAT_CSV = (ROOT / "data" / "audit" / "matched_strength_surrogate_split_baseline"
             / "per_patient_per_band.csv")

# Sequential density ramp: transparent-white (low) → saturated green (ridge),
# on-theme (green = trace) and clean on a transparent background.
DENS_CMAP = LinearSegmentedColormap.from_list(
    "trace_density",
    ["#ffffff", "#d6ece0", "#9fd4b4", "#5fac7e", "#1a7c3e", "#0a3a1d"],
)
C_PRO, C_ANTI, C_NULL = "#1a7c3e", "#c0392b", "#b9b9b9"


def movers_uv(dt, dr, rng):
    dt = np.round(dt, 12); dr = np.round(dr, 12)
    vt, ct = np.unique(dt, return_counts=True)
    vr, cr = np.unique(dr, return_counts=True)
    k = (dt != vt[ct.argmax()]) & (dr != vr[cr.argmax()])
    return (p18.rank_random_tiebreak(dt[k], rng),
            p18.rank_random_tiebreak(dr[k], rng))


def main() -> Path:
    rng = np.random.default_rng(SEED)
    pairs = p18.coph_obs_pairs(BAND)
    U, V = [], []
    for dt, dr in pairs.values():
        u, v = movers_uv(dt, dr, rng)
        U.append(u); V.append(v)
    U = np.concatenate(U); V = np.concatenate(V)

    xy = np.linspace(0.0, 1.0, GRID)
    xx, yy = np.meshgrid(xy, xy)
    dens = gaussian_kde(np.vstack([U, V]), bw_method=KDE_BW)(
        np.vstack([xx.ravel(), yy.ravel()])).reshape(GRID, GRID)   # [v, u]

    fig = plt.figure(figsize=(8.6, 10.4))
    gs = fig.add_gridspec(3, 2, width_ratios=[4.2, 1.0],
                          height_ratios=[1.0, 4.2, 1.25],
                          wspace=0.04, hspace=0.06,
                          left=0.11, right=0.97, top=0.97, bottom=0.08)
    ax = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax)
    ax_lab = fig.add_subplot(gs[0, 1])
    ax_strip = fig.add_subplot(gs[2, :])

    # ---- main density ----
    levels = np.linspace(dens.max() * 0.06, dens.max(), 11)
    ax.contourf(xx, yy, dens, levels=levels, cmap=DENS_CMAP, extend="max")
    ax.contour(xx, yy, dens, levels=levels[::2], colors="white", linewidths=0.5,
               alpha=0.6)
    ax.plot([0, 1], [0, 1], color="0.25", lw=1.1, ls="--", alpha=0.8, zorder=4)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_box_aspect(1.0)
    ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
    ax.set_xticklabels(["0", "½", "1"]); ax.set_yticklabels(["0", "½", "1"])
    ax.set_xlabel(r"$u$ = within-patient rank of "
                  r"$\Delta^{\mathrm{coph}}_{\mathrm{task}}(i,j)$", color="0.18")
    ax.set_ylabel(r"$v$ = within-patient rank of "
                  r"$\Delta^{\mathrm{coph}}_{\mathrm{rest}}(i,j)$", color="0.18")
    for s in ax.spines.values():
        s.set_edgecolor("0.55"); s.set_linewidth(0.9)

    # ---- conditional-split marginals ----
    edges = np.linspace(0, 1, 31)
    def split_hist(a, b):                       # P(a) for b in upper vs lower half
        hi, _ = np.histogram(a[b > 0.5], bins=edges, density=True)
        lo, _ = np.histogram(a[b <= 0.5], bins=edges, density=True)
        n = max(hi.max(), lo.max(), 1e-9)
        return lo / n, hi / n
    tl, th = split_hist(U, V)
    ax_top.stairs(tl, edges, fill=True, facecolor=C_ANTI, alpha=0.40,
                  edgecolor="#7f1d1d", lw=0.9)
    ax_top.stairs(th, edges, fill=True, facecolor=C_PRO, alpha=0.45,
                  edgecolor="#0a3a1d", lw=0.9)
    rl, rh = split_hist(V, U)
    ax_right.stairs(rl, edges, fill=True, orientation="horizontal",
                    facecolor=C_ANTI, alpha=0.40, edgecolor="#7f1d1d", lw=0.9)
    ax_right.stairs(rh, edges, fill=True, orientation="horizontal",
                    facecolor=C_PRO, alpha=0.45, edgecolor="#0a3a1d", lw=0.9)
    for a in (ax_top, ax_right):
        a.set_xticks([]); a.set_yticks([])
        for s in a.spines.values():
            s.set_visible(False)
    ax_top.set_ylim(0, 1.05); ax_right.set_xlim(0, 1.05)

    ax_lab.axis("off")
    ax_lab.text(0.5, 0.5, BRAIN_BAND_TEX_DICT.get(BAND, BAND),
                transform=ax_lab.transAxes, ha="center", va="center",
                fontsize=34, fontweight="bold", color="0.18")

    # ---- bottom strip: per-patient ρ vs own matched-strength null ----
    df = pd.read_csv(PERPAT_CSV)
    d = df[df.band == BAND].sort_values("obs_rho").reset_index(drop=True)
    for i, row in d.iterrows():
        ax_strip.plot([row.surr_p5, row.surr_p95], [i, i], color=C_NULL, lw=4.5,
                      solid_capstyle="round", alpha=0.6, zorder=1)
        ax_strip.plot([row.surr_p50], [i], marker="|", ms=7, color="0.5",
                      mew=1.2, zorder=2)
        pro = row.obs_rho > row.surr_p50
        col = C_PRO if pro else C_ANTI
        ax_strip.plot([row.surr_p50, row.obs_rho], [i, i], color=col, lw=1.6,
                      zorder=3)
        beyond = (row.obs_rho > row.surr_p95) or (row.obs_rho < row.surr_p5)
        ax_strip.plot([row.obs_rho], [i], marker="o", ms=7, zorder=4,
                      mfc=(col if beyond else "white"), mec=col, mew=1.4)
    ax_strip.axvline(0.0, color="0.45", ls="--", lw=0.8)
    ax_strip.set_ylim(-0.8, len(d) - 0.2)
    ax_strip.set_yticks([])
    ax_strip.set_xlabel(r"per-patient $\rho^{\mathrm{coph}}_{\mathrm{split}}$ "
                        r"vs own matched-strength null  (8/10 $>$ null median, "
                        r"7/10 clear it)", color="0.18")
    ax_strip.spines[["left", "right", "top"]].set_visible(False)
    ax_strip.spines["bottom"].set_color("0.6")

    out_dir = ROOT / "data" / "preprint" / "figures" / "all_bands"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "fig_beta_perpair_headline_density_coph.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  pooled movers: n={U.size}  density max={dens.max():.3f}")
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    main()
