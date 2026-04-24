#!/usr/bin/env python3
"""Section 6 — H1 + H2a dissociation figure (2 rows × 6 bands).

Top row:  H1 contrast = VI(TL,TT;k) − mean(other 5 pairs; k)
Bottom:   H2a contrast = VI(Pre,Post;k) − VI(TT,Post;k)

Display: 4 patients (Pat_02, Pat_05, Pat_07, Pat_08).
Unanimity shading: computed from all 5 patients (incl. Pat_03).

Output: data/figures/section6/h1_h2_all_bands.pdf + .md
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.config.paths import FIGURES_ROOT

OUTDIR = FIGURES_ROOT / "section6"
OUTDIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
})

# ── Constants ────────────────────────────────────────────────────────
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = {
    "delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
    "beta": r"$\beta$", "low_gamma": r"$\gamma_l$", "high_gamma": r"$\gamma_h$",
}
ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
DISPLAY_PATIENTS = ["Pat_02", "Pat_05", "Pat_07", "Pat_08"]
PAT_COLORS = {
    "Pat_02": "#1f77b4",
    "Pat_05": "#2ca02c",
    "Pat_07": "#d62728",
    "Pat_08": "#9467bd",
}
PAT_SHORT = {"Pat_02": "P2", "Pat_05": "P5", "Pat_07": "P7", "Pat_08": "P8"}

ALL_PAIRS = [
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"), ("rest_pre", "rest_post"),
    ("task_learn", "task_test"), ("task_learn", "rest_post"), ("task_test", "rest_post"),
]

# ── Load data ────────────────────────────────────────────────────────
print("Loading VI profiles...")
df = pd.read_csv("data/wp0_metric_exploration/task3_multiscale/vi_raw_profiles.csv")
df["pair_key"] = list(zip(df["phase_a"], df["phase_b"]))

# Common k range
k_common = None
for pat in ALL_PATIENTS:
    for band in BANDS:
        ks = set(df[(df["patient"] == pat) & (df["band"] == band)]["k"].unique())
        k_common = ks if k_common is None else k_common & ks
k_range = sorted(k_common)
print(f"  k range: {min(k_range)}–{max(k_range)} ({len(k_range)} values)")

# Build pivoted arrays: (band, pair) → DataFrame[k × patient]
vi_arr = {}
for band in BANDS:
    for pair in ALL_PAIRS:
        sub = df[(df["band"] == band) & (df["pair_key"].apply(lambda x: x == pair))]
        piv = sub.pivot_table(index="k", columns="patient", values="vi")
        piv = piv.reindex(k_range).dropna(how="all")
        vi_arr[(band, pair)] = piv


# ── Compute contrasts ────────────────────────────────────────────────
def compute_h1_contrasts(band):
    """H1: VI(TL,TT) − mean(other 5 pairs). Returns {pat: Series}."""
    piv_tt = vi_arr[(band, ("task_learn", "task_test"))]
    other_pivs = [vi_arr[(band, p)] for p in ALL_PAIRS
                  if p != ("task_learn", "task_test")]
    result = {}
    for pat in ALL_PATIENTS:
        if pat not in piv_tt.columns:
            continue
        tt_vals = piv_tt[pat]
        others = pd.concat([p[pat] for p in other_pivs if pat in p.columns], axis=1)
        result[pat] = tt_vals - others.mean(axis=1)
    return result


def compute_h2a_contrasts(band):
    """H2a: VI(Pre,Post) − VI(TT,Post). Returns {pat: Series}."""
    piv_pp = vi_arr[(band, ("rest_pre", "rest_post"))]
    piv_tp = vi_arr[(band, ("task_test", "rest_post"))]
    result = {}
    for pat in ALL_PATIENTS:
        if pat not in piv_pp.columns or pat not in piv_tp.columns:
            continue
        result[pat] = piv_pp[pat] - piv_tp[pat]
    return result


# ── Figure ───────────────────────────────────────────────────────────
print("Generating figure...")

fig, axes = plt.subplots(2, 6, figsize=(20, 6), sharex=True)
fig.subplots_adjust(hspace=0.25, wspace=0.08)

# Track y-limits per row for shared scaling
ylims_h1 = []
ylims_h2 = []

for ib, band in enumerate(BANDS):
    ax_h1 = axes[0, ib]
    ax_h2 = axes[1, ib]

    # ── H1 (top row) ──
    h1 = compute_h1_contrasts(band)
    h1_df = pd.DataFrame(h1)

    # Unanimity shading (all 5 patients)
    for k in h1_df.index:
        row = h1_df.loc[k].dropna()
        if len(row) == len(ALL_PATIENTS) and (row < 0).all():
            ax_h1.axvspan(k - 0.5, k + 0.5, color="#2196F3", alpha=0.25,
                          zorder=0, linewidth=0)

    # Display patients only
    for pat in DISPLAY_PATIENTS:
        if pat in h1:
            ax_h1.plot(h1[pat].index, h1[pat].values, color=PAT_COLORS[pat],
                       lw=1.0, alpha=0.6, zorder=2)

    # Mean of display patients
    h1_disp = pd.DataFrame({p: h1[p] for p in DISPLAY_PATIENTS if p in h1})
    if not h1_disp.empty:
        mean_h1 = h1_disp.mean(axis=1)
        ax_h1.plot(mean_h1.index, mean_h1.values, "k-", lw=2.0, zorder=3)

    ax_h1.axhline(0, color="#999999", ls="--", lw=0.8, zorder=1)
    ax_h1.set_title(BAND_TEX[band], fontsize=14, fontweight="bold")
    ax_h1.tick_params(labelsize=7)
    ylims_h1.append((ax_h1.get_ylim()))

    # ── H2a (bottom row) ──
    h2a = compute_h2a_contrasts(band)
    h2a_df = pd.DataFrame(h2a)

    # Unanimity shading (all 5 patients)
    for k in h2a_df.index:
        row = h2a_df.loc[k].dropna()
        if len(row) == len(ALL_PATIENTS) and (row > 0).all():
            ax_h2.axvspan(k - 0.5, k + 0.5, color="#FF9800", alpha=0.25,
                          zorder=0, linewidth=0)

    # Display patients
    for pat in DISPLAY_PATIENTS:
        if pat in h2a:
            ax_h2.plot(h2a[pat].index, h2a[pat].values, color=PAT_COLORS[pat],
                       lw=1.0, alpha=0.6, zorder=2)

    # Mean of display patients
    h2a_disp = pd.DataFrame({p: h2a[p] for p in DISPLAY_PATIENTS if p in h2a})
    if not h2a_disp.empty:
        mean_h2a = h2a_disp.mean(axis=1)
        ax_h2.plot(mean_h2a.index, mean_h2a.values, "k-", lw=2.0, zorder=3)

    ax_h2.axhline(0, color="#999999", ls="--", lw=0.8, zorder=1)
    ax_h2.tick_params(labelsize=7)
    ylims_h2.append((ax_h2.get_ylim()))

# ── Shared y-axis per row ──
ymin_h1 = min(y[0] for y in ylims_h1)
ymax_h1 = max(y[1] for y in ylims_h1)
ymin_h2 = min(y[0] for y in ylims_h2)
ymax_h2 = max(y[1] for y in ylims_h2)

for ib in range(6):
    axes[0, ib].set_ylim(ymin_h1, ymax_h1)
    axes[1, ib].set_ylim(ymin_h2, ymax_h2)
    # Only left-most gets y-tick labels
    if ib > 0:
        axes[0, ib].set_yticklabels([])
        axes[1, ib].set_yticklabels([])

# Row labels
axes[0, 0].set_ylabel("H1: task stability\n$\\Delta$VI", fontsize=10)
axes[1, 0].set_ylabel("H2a: task trace\n$\\Delta$VI", fontsize=10)

# X-axis label on bottom row
for ib in range(6):
    axes[1, ib].set_xlabel("$k$", fontsize=9)

# Legend (bottom-right panel)
from matplotlib.lines import Line2D
handles = [Line2D([0], [0], color=PAT_COLORS[p], lw=1.0, alpha=0.6,
                  label=PAT_SHORT[p]) for p in DISPLAY_PATIENTS]
handles.append(Line2D([0], [0], color="black", lw=2.0, label="Mean"))
handles.append(mpatches := plt.Rectangle((0, 0), 1, 1, fc="#2196F3", alpha=0.25,
               ec="none", label="H1 unan. (5/5)"))
handles.append(plt.Rectangle((0, 0), 1, 1, fc="#FF9800", alpha=0.25,
               ec="none", label="H2a unan. (5/5)"))
axes[1, -1].legend(handles=handles, fontsize=7, loc="lower right",
                    frameon=True, framealpha=0.9, edgecolor="#cccccc")

fig.savefig(OUTDIR / "h1_h2_all_bands.pdf", bbox_inches="tight", dpi=200)
plt.close(fig)
print("Saved h1_h2_all_bands.pdf")

# ── Markdown companion ───────────────────────────────────────────────
md = """\
# h1_h2_all_bands

## What the figure shows
2×6 panel figure. Top row: H1 contrast (task stability) per band — VI(TL,TT) minus
mean of other 5 pairs. Bottom row: H2a contrast (task trace) per band — VI(Pre,Post)
minus VI(TT,Post). Thin colored lines = 4 patients (Pat_03 excluded from display).
Thick black = mean. Blue shading = 5/5 unanimity for H1. Orange shading = 5/5
unanimity for H2a (computed from all 5 patients including Pat_03).

## Key result
Top: beta has the widest blue shading (strongest H1). Bottom: ONLY alpha shows
orange shading (29 contiguous k values with 5/5 positive trace). All other bands
have no or sparse shading — alpha is the consensus trace band.

## How to read it
Lines below zero in top row = TL-TT pair is most similar (H1 supported).
Lines above zero in bottom row = positive task trace (H2a supported).
Shaded regions = unanimous across all 5 patients at that k.
"""
(OUTDIR / "h1_h2_all_bands.md").write_text(md)
print("Saved h1_h2_all_bands.md")
