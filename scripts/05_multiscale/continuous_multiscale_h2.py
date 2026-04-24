#!/usr/bin/env python3
"""Continuous multiscale H2 task-trace profiles.

Instead of cutting dendrograms at fixed k, we sweep the threshold height h
continuously through the dendrogram and compute VI(h) between phase pairs.

For a dendrogram with merge heights h_1 < h_2 < ... < h_{N-1}, cutting at
height h gives a partition: all merges below h are kept, above h are cut.
This is scipy.cluster.hierarchy.fcluster(Z, t=h, criterion='distance').

The number of clusters k(h) varies continuously from N (h=0) to 1 (h=max).
VI(h) between two dendrograms measures how different their partitions are
at threshold h.

H2a(h) = VI(Pre,Post; h) - VI(TT,Post; h)
  Positive = at height h, rest_post looks more like task_test than like rest_pre.

This gives a truly continuous profile over the full hierarchy.

Output:
  data/figures/metric_exploration/continuous_h2/
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster

import argparse
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.config.paths import FIGURES_ROOT, lrg_cache_for

_p = argparse.ArgumentParser(add_help=False)
_p.add_argument("--fc-method", default="imcoh_abs",
                choices=["msc", "corr", "imcoh_abs", "imcoh_sq"])
_FC_ARGS, _ = _p.parse_known_args()
FC_METHOD = _FC_ARGS.fc_method
CACHE_ROOT = lrg_cache_for(FC_METHOD)

OUTDIR = FIGURES_ROOT / "metric_exploration" / "continuous_h2" / FC_METHOD
OUTDIR.mkdir(parents=True, exist_ok=True)

from lrg_eegfc.config.const import PATIENTS_4PHASE
PATIENTS = PATIENTS_4PHASE
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_LABELS = {
    "delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
    "beta": r"$\beta$", "low_gamma": r"$\gamma_l$", "high_gamma": r"$\gamma_h$",
}
BAND_COLORS = {
    "delta": "#1f77b4", "theta": "#ff7f0e", "alpha": "#2ca02c",
    "beta": "#d62728", "low_gamma": "#9467bd", "high_gamma": "#8c564b",
}
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]

# Height grid: 200 points from 0.01 to 0.99
H_GRID = np.linspace(0.01, 0.99, 200)


from lrg_eegfc.utils.metrics import compute_vi  # canonical implementation


def get_partition_at_height(Z: np.ndarray, h: float) -> np.ndarray:
    """Get flat partition by cutting dendrogram at height h."""
    return fcluster(Z, t=h, criterion="distance")


def get_k_at_height(Z: np.ndarray, h: float) -> int:
    """Number of clusters at height h."""
    return len(np.unique(get_partition_at_height(Z, h)))


# ── Load all linkage matrices ─────────────────────────────────────────
print("Loading LRG results...")
linkages = {}
for pat in PATIENTS:
    for band in BANDS:
        for phase in PHASES:
            try:
                res = load_lrg_result(pat, phase, band, fc_method=FC_METHOD,
                                      cache_root=CACHE_ROOT)
                linkages[(pat, band, phase)] = res.linkage_matrix
            except Exception as e:
                print(f"  SKIP {pat}/{band}/{phase}: {e}")

print(f"Loaded {len(linkages)} linkage matrices")

# ── Compute VI(h) for key pairs ───────────────────────────────────────
print("Computing VI(h) curves...")

KEY_PAIRS = [
    ("rest_pre", "rest_post", "Pre↔Post"),
    ("task_test", "rest_post", "TT↔Post"),
    ("rest_pre", "task_test", "Pre↔TT"),
    ("task_learn", "task_test", "TL↔TT"),
]

# Store: vi_curves[pat][band][pair_label] = array of VI values at H_GRID
vi_curves = {}
for pat in PATIENTS:
    vi_curves[pat] = {}
    for band in BANDS:
        vi_curves[pat][band] = {}
        for p1, p2, label in KEY_PAIRS:
            key1 = (pat, band, p1)
            key2 = (pat, band, p2)
            if key1 not in linkages or key2 not in linkages:
                continue
            Z1 = linkages[key1]
            Z2 = linkages[key2]
            vis = np.zeros(len(H_GRID))
            for ih, h in enumerate(H_GRID):
                lab1 = get_partition_at_height(Z1, h)
                lab2 = get_partition_at_height(Z2, h)
                vis[ih] = compute_vi(lab1, lab2)
            vi_curves[pat][band][label] = vis
    print(f"  {pat} done")

# Also store k(h) for reference (from any dendrogram — they all have ~117 nodes)
# Use Pat_02/alpha/rest_pre as reference
ref_Z = linkages[("Pat_02", "alpha", "rest_pre")]
k_at_h = np.array([get_k_at_height(ref_Z, h) for h in H_GRID])


# ── Compute H2 contrasts ─────────────────────────────────────────────
# H2a(h) = VI(Pre,Post; h) - VI(TT,Post; h)
# H2b(h) = VI(Pre,TT; h) - VI(TT,Post; h)
h2a = {}  # h2a[pat][band] = array
h2b = {}
for pat in PATIENTS:
    h2a[pat] = {}
    h2b[pat] = {}
    for band in BANDS:
        c = vi_curves[pat][band]
        if "Pre↔Post" in c and "TT↔Post" in c:
            h2a[pat][band] = c["Pre↔Post"] - c["TT↔Post"]
        if "Pre↔TT" in c and "TT↔Post" in c:
            h2b[pat][band] = c["Pre↔TT"] - c["TT↔Post"]


# ── Figure 1: H2a continuous profiles per band ───────────────────────
fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True, sharey=True)
axes_flat = axes.ravel()

for ib, band in enumerate(BANDS):
    ax = axes_flat[ib]
    color = BAND_COLORS[band]

    # Individual patients
    all_curves = []
    for pat in PATIENTS:
        if band in h2a[pat]:
            curve = h2a[pat][band]
            all_curves.append(curve)
            ax.plot(H_GRID, curve, color=color, alpha=0.2, linewidth=0.8)

    if all_curves:
        all_curves = np.array(all_curves)
        mean_c = np.mean(all_curves, axis=0)
        ax.plot(H_GRID, mean_c, color=color, linewidth=2.5, zorder=5)

        # Shade unanimous regions
        all_pos = np.all(all_curves > 0, axis=0)
        all_neg = np.all(all_curves < 0, axis=0)
        ax.fill_between(H_GRID, mean_c.min() - 0.1, mean_c.max() + 0.1,
                         where=all_pos, color="#4CAF50", alpha=0.15, zorder=0)
        ax.fill_between(H_GRID, mean_c.min() - 0.1, mean_c.max() + 0.1,
                         where=all_neg, color="#F44336", alpha=0.15, zorder=0)

    ax.axhline(0, color="k", linewidth=0.8, linestyle="--")
    ax.set_title(BAND_LABELS[band], fontsize=13, fontweight="bold")
    if ib >= 3:
        ax.set_xlabel("Threshold height h")
    if ib % 3 == 0:
        ax.set_ylabel("H2a contrast")
    ax.grid(alpha=0.3)

    # Add secondary x-axis showing approximate k
    if ib == 2:
        ax2 = ax.twiny()
        k_ticks = [2, 5, 10, 20, 50, 100]
        h_for_k = []
        for kt in k_ticks:
            idx = np.argmin(np.abs(k_at_h - kt))
            h_for_k.append(H_GRID[idx])
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(h_for_k)
        ax2.set_xticklabels([f"k≈{kt}" for kt in k_ticks], fontsize=7)

from matplotlib.patches import Patch
from matplotlib.lines import Line2D
handles = [
    Line2D([0], [0], color="gray", linewidth=2.5, label="Mean (4 patients)"),
    Line2D([0], [0], color="gray", linewidth=0.8, alpha=0.4, label="Individual patients"),
    Patch(facecolor="#4CAF50", alpha=0.3, label="Unanimous: task persists"),
    Patch(facecolor="#F44336", alpha=0.3, label="Unanimous: brain recovers"),
]
fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9,
           bbox_to_anchor=(0.5, -0.02))
fig.suptitle("H2a: VI(Pre,Post; h) − VI(TT,Post; h)\n"
             "Continuous profile over dendrogram height — "
             "positive = rest_post retains task structure",
             fontsize=12, y=1.02)
fig.tight_layout()
fig.savefig(OUTDIR / "h2a_continuous_per_band.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved h2a_continuous_per_band.pdf")


# ── Figure 2: Money plot — all bands on one axis, mean only ──────────
fig, ax = plt.subplots(figsize=(12, 6))

for band in BANDS:
    all_curves = []
    for pat in PATIENTS:
        if band in h2a[pat]:
            all_curves.append(h2a[pat][band])
    if all_curves:
        mean_c = np.mean(all_curves, axis=0)
        ax.plot(H_GRID, mean_c, color=BAND_COLORS[band], linewidth=2.5,
                label=BAND_LABELS[band], zorder=5)
        # Thin envelope
        lo = np.min(all_curves, axis=0)
        hi = np.max(all_curves, axis=0)
        ax.fill_between(H_GRID, lo, hi, color=BAND_COLORS[band], alpha=0.08)

ax.axhline(0, color="k", linewidth=1, linestyle="--")
ax.set_xlabel("Threshold height h (low = fine communities, high = coarse)", fontsize=11)
ax.set_ylabel("H2a: VI(Pre,Post) − VI(TT,Post)", fontsize=11)
ax.set_title("Task trace across the full hierarchy\n"
             "Positive = task persists in rest_post | Negative = brain recovers\n"
             "Envelopes = min-max across 4 patients",
             fontsize=12)

# Secondary axis with approximate k
ax2 = ax.twiny()
k_ticks = [2, 3, 5, 10, 20, 50, 100]
h_for_k = []
for kt in k_ticks:
    idx = np.argmin(np.abs(k_at_h - kt))
    h_for_k.append(H_GRID[idx])
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(h_for_k)
ax2.set_xticklabels([f"k≈{kt}" for kt in k_ticks], fontsize=8)
ax2.set_xlabel("Approximate number of communities", fontsize=9)

ax.legend(fontsize=10, loc="best")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(OUTDIR / "h2a_all_bands_continuous.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved h2a_all_bands_continuous.pdf")


# ── Figure 3: Unanimity map — continuous version ─────────────────────
# For each (band, h), count patients with positive H2a
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ic, (contrast_dict, title, fname_suffix) in enumerate([
    (h2a, "H2a: VI(Pre,Post) − VI(TT,Post)", "h2a"),
    (h2b, "H2b: VI(Pre,TT) − VI(TT,Post)", "h2b"),
]):
    ax = axes[ic]
    mat = np.zeros((len(BANDS), len(H_GRID)))
    for ib, band in enumerate(BANDS):
        curves = []
        for pat in PATIENTS:
            if band in contrast_dict[pat]:
                curves.append(contrast_dict[pat][band])
        if curves:
            curves = np.array(curves)
            n_pos = np.sum(curves > 0, axis=0)  # 0-4
            mat[ib, :] = n_pos

    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("unanimity",
        ["#D32F2F", "#EF9A9A", "white", "#90CAF9", "#1565C0"], N=256)

    im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=4, aspect="auto",
                   extent=[H_GRID[0], H_GRID[-1], len(BANDS) - 0.5, -0.5])
    cbar = plt.colorbar(im, ax=ax, shrink=0.8, ticks=[0, 1, 2, 3, 4])
    cbar.set_ticklabels(["0/4\nrecovers", "1/4", "2/4", "3/4", "4/4\npersists"])

    ax.set_yticks(range(len(BANDS)))
    ax.set_yticklabels([BAND_LABELS[b] for b in BANDS])
    ax.set_xlabel("Threshold height h")
    ax.set_title(title, fontsize=10)

    # Add k reference lines
    for kt in [3, 5, 10, 20]:
        idx = np.argmin(np.abs(k_at_h - kt))
        ax.axvline(H_GRID[idx], color="k", linewidth=0.5, alpha=0.3, linestyle=":")
        ax.text(H_GRID[idx], -0.45, f"k≈{kt}", fontsize=7, ha="center")

fig.suptitle("Continuous unanimity map: how many patients show task persistence\n"
             "Blue = all persist | Red = all recover | White = split",
             fontsize=12, y=1.05)
fig.tight_layout()
fig.savefig(OUTDIR / "h2_continuous_unanimity.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved h2_continuous_unanimity.pdf")


# ── Figure 4: VI(h) curves for all pairs — per band ──────────────────
# Show the raw VI curves (not contrasts) so we can see the full picture
fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True, sharey=True)
axes_flat = axes.ravel()

PAIR_COLORS_MAP = {
    "TL↔TT": "#2196F3",
    "Pre↔Post": "#4CAF50",
    "Pre↔TT": "#9E9E9E",
    "TT↔Post": "#E91E63",
}

for ib, band in enumerate(BANDS):
    ax = axes_flat[ib]
    for label, pcolor in PAIR_COLORS_MAP.items():
        all_c = []
        for pat in PATIENTS:
            if label in vi_curves[pat][band]:
                all_c.append(vi_curves[pat][band][label])
        if all_c:
            mean_c = np.mean(all_c, axis=0)
            ax.plot(H_GRID, mean_c, color=pcolor, linewidth=2, label=label)
            lo = np.min(all_c, axis=0)
            hi = np.max(all_c, axis=0)
            ax.fill_between(H_GRID, lo, hi, color=pcolor, alpha=0.1)

    ax.set_title(BAND_LABELS[band], fontsize=13, fontweight="bold")
    if ib >= 3:
        ax.set_xlabel("Height h")
    if ib % 3 == 0:
        ax.set_ylabel("VI(h)")
    ax.grid(alpha=0.3)
    if ib == 0:
        ax.legend(fontsize=7, loc="upper left")

fig.suptitle("Raw VI(h) curves for key phase pairs (mean ± range across 4 patients)\n"
             "Lower = more similar at that hierarchy level",
             fontsize=12, y=1.02)
fig.tight_layout()
fig.savefig(OUTDIR / "vi_curves_all_pairs.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved vi_curves_all_pairs.pdf")


# ── Print where unanimous regions are ─────────────────────────────────
print("\n" + "=" * 70)
print("CONTINUOUS UNANIMITY REGIONS")
print("=" * 70)
for contrast_name, contrast_dict in [("H2a", h2a), ("H2b", h2b)]:
    print(f"\n{contrast_name}:")
    for band in BANDS:
        curves = []
        for pat in PATIENTS:
            if band in contrast_dict[pat]:
                curves.append(contrast_dict[pat][band])
        if not curves:
            continue
        curves = np.array(curves)
        all_pos = np.all(curves > 0, axis=0)
        all_neg = np.all(curves < 0, axis=0)

        # Find contiguous unanimous regions
        pos_regions = []
        neg_regions = []
        for mask, regions in [(all_pos, pos_regions), (all_neg, neg_regions)]:
            in_region = False
            start = 0
            for i in range(len(H_GRID)):
                if mask[i] and not in_region:
                    start = i
                    in_region = True
                elif not mask[i] and in_region:
                    regions.append((H_GRID[start], H_GRID[i - 1],
                                    k_at_h[start], k_at_h[i - 1]))
                    in_region = False
            if in_region:
                regions.append((H_GRID[start], H_GRID[-1],
                                k_at_h[start], k_at_h[-1]))

        line = f"  {BAND_LABELS[band]:>12s}: "
        parts = []
        for h_lo, h_hi, k_hi, k_lo in pos_regions:
            parts.append(f"PERSIST h=[{h_lo:.2f},{h_hi:.2f}] (k≈{k_lo}-{k_hi})")
        for h_lo, h_hi, k_hi, k_lo in neg_regions:
            parts.append(f"RECOVER h=[{h_lo:.2f},{h_hi:.2f}] (k≈{k_lo}-{k_hi})")
        if parts:
            print(line + " | ".join(parts))
        else:
            print(line + "no unanimous region")

print("\nDone!")
