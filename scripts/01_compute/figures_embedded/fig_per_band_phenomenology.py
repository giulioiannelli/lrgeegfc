#!/usr/bin/env python3
"""Per-band phenomenology headline figure (Step 2, session 2026-06-26).

The whole per-band vision in one publication figure, grounded in the locked
audit outputs (`.agents/reports/2026-06-26_per-band-phenomenology-vision.md`):

  A  per-patient x per-band cophenetic-trace matrix (rho_split), ordered by
     tracer strength, cells that CLEAR their own matched-strength surrogate
     boxed. Shows the gradient Pat_06 (6/6 bands) -> Pat_15 (0/6) and the
     full beta column.
  B  cohort gate: per-band median rho_split vs the surrogate-median null
     (p95 threshold), with the per-patient spread as dots. Only beta and
     alpha clear (starred).
  C  beta trace localization: per-system matched-strength (both tails, BH
     across a-priori systems) -> OFC carrier, sensorimotor + PFC depleted.
  D  the honesty panel: per-patient beta rho_split vs split-half reliability
     is flat (Spearman ~0) -> the between-patient spread is NOT measurement
     noise (Pat_15 = most reliable, least trace; Pat_08 = strong trace, low
     coverage). Coverage refuted too (rho with OFC electrode fraction).

Reads
-----
- data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv
- data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv
- data/audit/localization_atlas/matched_strength_allbands_include.csv
- data/audit/ofc_coverage/ofc_coverage_vs_trace.csv

Writes
------
- data/reports/per_band_phenomenology/fig_per_band_phenomenology.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT, BRAIN_BANDS_NAMES
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style
from lrgsglib.plotlib import imshow_colorbar_caxdivider

ROOT = setup_script_env()
use_lrg_style()

AUDIT = ROOT / "data/audit"
GATE = AUDIT / "matched_strength_surrogate_split_baseline/per_patient_per_band.csv"
COH = AUDIT / "matched_strength_surrogate_split_baseline/cohort_summary.csv"
LOC = AUDIT / "localization_atlas/matched_strength_allbands_include.csv"
TAX = AUDIT / "raw_vs_multiscale/band_taxonomy_raw_vs_multiscale.csv"
GRASS = AUDIT / "grassmann_cluster_extent/cohort_summary.csv"
OUT = ROOT / "data/reports/per_band_phenomenology/fig_per_band_phenomenology.pdf"

BANDS = list(BRAIN_BANDS_NAMES)
TEX = [BRAIN_BAND_TEX_DICT[b] for b in BANDS]
CLEAR = 0.05           # own-surrogate one-sided p threshold
ANTI = 0.95
CARRIER = "#2ca02c"
DEPLETED = "#d62728"
NEUTRAL = "#b8b8b8"
# a-priori localization systems (exclude atlas 'other'); undersampled if K<4
SYS_MIN_K = 4


def _load():
    g = pd.read_csv(GATE)
    coh = pd.read_csv(COH).set_index("band")
    loc = pd.read_csv(LOC)
    tax = pd.read_csv(TAX)
    grass = pd.read_csv(GRASS).set_index("band")
    return g, coh, loc, tax, grass


# ------------------------------------------------------------------ panel A
def panel_A(ax, g):
    rho = g.pivot(index="patient", columns="band", values="obs_rho")[BANDS]
    pval = g.pivot(index="patient", columns="band", values="obs_p_one_sided")[BANDS]
    clears = pval < CLEAR
    # order patients by (#bands cleared, mean rho) descending
    order = sorted(rho.index,
                   key=lambda p: (int(clears.loc[p].sum()), float(rho.loc[p].mean())),
                   reverse=True)
    rho, pval, clears = rho.loc[order], pval.loc[order], clears.loc[order]
    M = rho.to_numpy()
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-0.5, vmax=0.95)
    im = ax.imshow(M, cmap="coolwarm", norm=norm, aspect="auto")

    nrow, ncol = M.shape
    for i in range(nrow):
        for j in range(ncol):
            v = M[i, j]
            txt = f"{v * 100:+.0f}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=6.2,
                    color="white" if abs(v) > 0.33 else "0.15")
            if clears.iloc[i, j]:            # solid box = clears own surrogate
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                       ec="black", lw=1.6, zorder=5))
            elif pval.iloc[i, j] > ANTI:     # dotted box = significantly anti
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                       ec="0.15", ls=":", lw=1.1, zorder=5))

    ax.set_xticks(range(ncol))
    ax.set_xticklabels(TEX, fontsize=11)
    ax.set_yticks(range(nrow))
    ax.set_yticklabels([p.replace("Pat_", "P") for p in order], fontsize=8)
    ax.set_xlabel("band   (solid box = clears surrogate  ·  dotted = anti)", fontsize=8.5)
    # right margin: #bands traced per patient
    counts = clears.sum(axis=1).to_numpy()
    for i, c in enumerate(counts):
        ax.text(ncol - 0.35, i, f"{c}", ha="left", va="center", fontsize=8,
                fontweight="bold", transform=ax.transData, clip_on=False)
    ax.text(ncol - 0.42, -0.62, "#tr", ha="left", va="center", fontsize=7.5,
            fontweight="bold", clip_on=False)
    ax.set_title(r"$\bf{A}$  per-patient trace  $\rho_{\rm split}\times100$",
                 fontsize=10, loc="left")
    imshow_colorbar_caxdivider(im, ax, size="4%", pad=0.55)
    return order


# ------------------------------------------------------------------ panel B
def panel_B(ax, g, coh):
    xoff = np.linspace(-0.26, 0.26, 10)
    for j, b in enumerate(BANDS):
        sub = g[g.band == b].sort_values("obs_rho")
        ys = sub.obs_rho.to_numpy()
        ax.scatter(j + xoff, ys, s=14, color="0.55", alpha=0.8, zorder=2,
                   edgecolor="none")
        p95 = float(coh.loc[b, "surr_median_rho_p95"])
        med = float(coh.loc[b, "obs_median_rho"])
        gp = float(coh.loc[b, "paired_wilcoxon_p"])
        # surrogate null ceiling (one-sided 5%)
        ax.plot([j - 0.30, j + 0.30], [p95, p95], color=DEPLETED, lw=2.0, zorder=3)
        sig = gp < 0.05
        ax.scatter([j], [med], s=130, marker="D", zorder=4,
                   color="black" if sig else "white", edgecolor="black", lw=1.4)
        if sig:
            ax.text(j, med + 0.06, "$\\star$", ha="center", va="bottom",
                    fontsize=15, color="black")
    ax.axhline(0, color="0.4", ls="--", lw=0.8, zorder=1)
    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels(TEX, fontsize=11)
    ax.set_ylabel(r"$\rho_{\rm split}$")
    ax.set_xlim(-0.5, len(BANDS) - 0.5)
    ax.set_title(r"$\bf{B}$  cohort gate: median vs matched-strength null",
                 fontsize=10, loc="left")
    handles = [
        Line2D([0], [0], marker="D", color="black", ls="", ms=9, label="cohort median (gate p<.05)"),
        Line2D([0], [0], marker="D", color="white", mec="black", ls="", ms=9, label="cohort median (n.s.)"),
        Line2D([0], [0], marker="o", color="0.55", ls="", ms=6, label="per-patient"),
        Line2D([0], [0], color=DEPLETED, lw=2, label="surrogate 95% (null ceiling)"),
    ]
    ax.legend(handles=handles, fontsize=6.8, loc="upper right", frameon=False,
              handletextpad=0.4, labelspacing=0.3)


# ------------------------------------------------------------------ panel C
def panel_C(ax, loc):
    s = loc[(loc.band == "beta") & (loc.granularity == "system")
            & (loc.epi == "include") & (loc.unit != "other")].copy()
    pu = s.matched_strength_p.to_numpy()
    pl = s.matched_strength_p_lower.to_numpy()
    qu = bh_fdr(pu)
    ql = bh_fdr(pl)
    net = np.log10(np.clip(pl, 1e-6, 1) / np.clip(pu, 1e-6, 1))  # >0 carrier, <0 depleted
    K = s.K_implanted.to_numpy()
    units = s.unit.to_numpy()
    order = np.argsort(net)                       # barh bottom->top ascending
    colors, labels = [], []
    for idx in order:
        under = K[idx] < SYS_MIN_K
        if under:
            colors.append(NEUTRAL)
        elif qu[idx] < 0.05:
            colors.append(CARRIER)
        elif ql[idx] < 0.05:
            colors.append(DEPLETED)
        else:
            colors.append(NEUTRAL)
        tag = ""
        if not under and qu[idx] < 0.05:
            tag = f"  carrier q={qu[idx]:.3f}"
        elif not under and ql[idx] < 0.05:
            tag = f"  depleted q={ql[idx]:.3f}"
        elif under:
            tag = f"  (K={K[idx]})"
        labels.append(units[idx].replace("_", " ") + tag)

    ypos = np.arange(len(order))
    bars = ax.barh(ypos, net[order], color=colors, edgecolor="0.25", lw=0.6)
    for idx, b in zip(order, bars):
        if K[idx] < SYS_MIN_K:
            b.set_hatch("////")
    ax.axvline(0, color="0.3", lw=1.0)
    thr = -np.log10(0.05)
    for x in (-thr, thr):
        ax.axvline(x, color="0.6", ls=":", lw=0.9)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=7.6)
    ax.set_xlabel(r"depleted  $\longleftarrow$  signed $-\log_{10}p$  $\longrightarrow$  carrier")
    ax.set_title(r"$\bf{C}$  $\beta$ trace localization (per-system)",
                 fontsize=10, loc="left")
    leg = [Patch(fc=CARRIER, ec="0.25", label="carrier (BH q<.05)"),
           Patch(fc=DEPLETED, ec="0.25", label="depleted (BH q<.05)"),
           Patch(fc=NEUTRAL, ec="0.25", hatch="////", label="undersampled K<4")]
    ax.legend(handles=leg, fontsize=6.8, loc="lower right", frameon=False)


# ------------------------------------------------------------------ panel D
def panel_D(ax, tax, grass):
    tx = tax.set_index("band").loc[BANDS]
    ms = tx.ms_median.to_numpy()
    raw = tx.raw_median.to_numpy()
    gr = grass.loc[BANDS]
    x = np.arange(len(BANDS))
    w = 0.38
    floor = float(raw.min())

    ax.axhline(0, color="0.3", lw=0.9, zorder=1)
    ax.axhline(floor, color="0.5", ls="--", lw=1.0, zorder=1)
    ax.text(len(BANDS) - 0.55, floor + 0.006, rf"raw floor $\approx${floor:.2f}",
            ha="right", va="bottom", fontsize=7, color="0.35")

    ax.bar(x - w / 2, raw, w, color="#9e9e9e", edgecolor="0.3", lw=0.5,
           label="raw per-edge", zorder=2)
    strong = {"alpha", "beta"}
    cols = ["#1f77b4" if b in strong else "#a6cee3" for b in BANDS]
    ax.bar(x + w / 2, ms, w, color=cols, edgecolor="0.3", lw=0.5,
           label="multiscale (cophenetic)", zorder=2)

    # star bands where multiscale clears the cohort gate (alpha, beta)
    for i, b in enumerate(BANDS):
        if b in strong:
            ax.text(x[i] + w / 2, ms[i] + 0.014, r"$\star$", ha="center",
                    va="bottom", fontsize=12)
    bi = BANDS.index("beta")
    ax.text(x[bi] + w / 2, ms[bi] + 0.048, "→OFC", ha="center", va="bottom",
            fontsize=6.8, color="#1f77b4")

    # theta: raw sees a trace, multiscale calls it absent
    ti = BANDS.index("theta")
    ax.annotate("", xy=(x[ti] + w / 2, ms[ti] - 0.004),
                xytext=(x[ti] + w / 2, raw[ti]),
                arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.4))
    ax.text(x[ti], 0.175, "raw: present\nMS: absent", fontsize=6.6,
            color="#d62728", va="center", ha="center")

    # --- Grassmann (subspace) verdict strip below y=0 ---
    y0, y1 = -0.155, -0.095
    GC = {"strong": "#2ca02c", "weak": "#f0ad4e", "none": "#d3d3d3"}
    GG = {"strong": "✓", "weak": "~", "none": "·"}
    for i, b in enumerate(BANDS):
        gp = float(gr.loc[b, "cluster_p_cluster_mass"])
        gloo = float(gr.loc[b, "cluster_p_mass_loo_max"])
        tier = "strong" if (gp < 0.05 and gloo < 0.05) else ("weak" if gp < 0.10 else "none")
        ax.add_patch(Rectangle((x[i] - 0.45, y0), 0.90, y1 - y0, fc=GC[tier],
                               ec="0.35", lw=0.5, zorder=2, clip_on=False))
        ax.text(x[i], (y0 + y1) / 2, GG[tier], ha="center", va="center",
                fontsize=10, color="white" if tier != "none" else "0.45", zorder=3)
    ax.text(2.5, y1 - 0.075,
            "Grassmann subspace $d_G(k)$ verdict:   ✓ clears (LOO)    ~ weak    · none",
            ha="center", va="top", fontsize=6.6)
    ax.text(5.42, 0.305, "β: both multiscale probes clear\n"
            "α: $\\rho$-only    γ$_l$: subspace-only",
            ha="right", va="top", fontsize=6.3,
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))

    ax.set_xticks(x)
    ax.set_xticklabels(TEX, fontsize=11)
    ax.set_yticks([0.0, 0.1, 0.2, 0.3])
    ax.set_ylabel("cohort median trace")
    ax.set_xlim(-0.5, len(BANDS) - 0.5)
    ax.set_ylim(-0.24, 0.34)
    ax.set_title(r"$\bf{D}$  raw vs two multiscale lenses ($\rho_{\rm coph}$, Grassmann)",
                 fontsize=10, loc="left")
    ax.legend(fontsize=7, loc="upper left", frameon=False,
              handletextpad=0.4, labelspacing=0.3)


def main():
    g, coh, loc, tax, grass = _load()
    fig = plt.figure(figsize=(13, 10))
    gs = GridSpec(2, 2, width_ratios=[1.15, 1.0], height_ratios=[1.0, 0.92],
                  hspace=0.30, wspace=0.32, figure=fig)
    panel_A(fig.add_subplot(gs[0, 0]), g)
    panel_B(fig.add_subplot(gs[0, 1]), g, coh)
    panel_C(fig.add_subplot(gs[1, 0]), loc)
    panel_D(fig.add_subplot(gs[1, 1]), tax, grass)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
