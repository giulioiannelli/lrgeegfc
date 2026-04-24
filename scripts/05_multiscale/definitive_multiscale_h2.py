#!/usr/bin/env python3
"""Definitive multiscale H2 visualization.

The key figure: a heatmap where x = hierarchy height (continuous), y = band,
and color = signed unanimity.  This shows AT A GLANCE where in the hierarchy
each band persists vs recovers, with no arbitrary scale choices.

Also produces:
  - Per-band panels with patient-level traces and unanimous shading
  - Summary table of contiguous unanimous regions

Output: data/figures/metric_exploration/continuous_h2/
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
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
BAND_TEX = {
    "delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
    "beta": r"$\beta$", "low_gamma": r"$\gamma_l$", "high_gamma": r"$\gamma_h$",
}
BAND_COLORS = {
    "delta": "#1f77b4", "theta": "#ff7f0e", "alpha": "#2ca02c",
    "beta": "#d62728", "low_gamma": "#9467bd", "high_gamma": "#8c564b",
}

# Dense height grid
H_GRID = np.linspace(0.005, 0.995, 300)


from lrg_eegfc.utils.metrics import compute_vi  # canonical implementation


# ── Load linkage matrices ─────────────────────────────────────────────
print("Loading LRG results...")
linkages = {}
for pat in PATIENTS:
    for phase in ["rest_pre", "task_learn", "task_test", "rest_post"]:
        for band in BANDS:
            try:
                res = load_lrg_result(pat, phase, band, fc_method=FC_METHOD,
                                      cache_root=CACHE_ROOT)
                linkages[(pat, band, phase)] = res.linkage_matrix
            except Exception:
                pass
print(f"Loaded {len(linkages)} linkages")

# k(h) reference
ref_Z = linkages[("Pat_02", "alpha", "rest_pre")]
k_at_h = np.array([len(np.unique(fcluster(ref_Z, h, criterion="distance")))
                    for h in H_GRID])

# ── Compute VI(h) and H2 contrasts ───────────────────────────────────
print("Computing VI(h) curves (this takes a minute)...")

# h2a[pat][band] = VI(Pre,Post;h) - VI(TT,Post;h)
# h2b[pat][band] = VI(Pre,TT;h) - VI(TT,Post;h)
# Also store raw VI for the overlay figure
h2a_curves = {p: {} for p in PATIENTS}
h2b_curves = {p: {} for p in PATIENTS}
vi_raw = {p: {b: {} for b in BANDS} for p in PATIENTS}

PAIRS_TO_COMPUTE = [
    ("rest_pre", "rest_post"),
    ("task_test", "rest_post"),
    ("rest_pre", "task_test"),
    ("task_learn", "task_test"),
]

for pat in PATIENTS:
    for band in BANDS:
        vi_dict = {}
        for p1, p2 in PAIRS_TO_COMPUTE:
            k1 = (pat, band, p1)
            k2 = (pat, band, p2)
            if k1 not in linkages or k2 not in linkages:
                continue
            Z1, Z2 = linkages[k1], linkages[k2]
            vis = np.zeros(len(H_GRID))
            for ih, h in enumerate(H_GRID):
                vis[ih] = compute_vi(
                    fcluster(Z1, h, criterion="distance"),
                    fcluster(Z2, h, criterion="distance"),
                )
            vi_dict[(p1, p2)] = vis
            vi_raw[pat][band][f"{p1}↔{p2}"] = vis

        if ("rest_pre", "rest_post") in vi_dict and ("task_test", "rest_post") in vi_dict:
            h2a_curves[pat][band] = vi_dict[("rest_pre", "rest_post")] - vi_dict[("task_test", "rest_post")]
        if ("rest_pre", "task_test") in vi_dict and ("task_test", "rest_post") in vi_dict:
            h2b_curves[pat][band] = vi_dict[("rest_pre", "task_test")] - vi_dict[("task_test", "rest_post")]
    print(f"  {pat} done")


def find_regions(mask, min_width=0.02):
    """Find contiguous True regions wider than min_width."""
    regions = []
    in_r = False
    start = 0
    for i in range(len(H_GRID)):
        if mask[i] and not in_r:
            start = i
            in_r = True
        elif not mask[i] and in_r:
            if H_GRID[i - 1] - H_GRID[start] >= min_width:
                regions.append((start, i - 1))
            in_r = False
    if in_r and H_GRID[-1] - H_GRID[start] >= min_width:
        regions.append((start, len(H_GRID) - 1))
    return regions


# ── FIGURE 1: The definitive heatmap ─────────────────────────────────
# x = height h, y = band, color = signed unanimity (continuous)
# Value = mean sign * agreement: mean(sign(H2a_i)) gives [-1, +1]
# +1 = all patients positive (persist), -1 = all negative (recover)

fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=True)

cmap = LinearSegmentedColormap.from_list(
    "persist_recover",
    ["#B71C1C", "#EF9A9A", "#FAFAFA", "#90CAF9", "#0D47A1"],
    N=256,
)

for ic, (contrast_curves, title) in enumerate([
    (h2a_curves, "H2a: does rest_post resemble task_test more than rest_pre?"),
    (h2b_curves, "H2b: does rest_post resemble task_test more than rest_pre did?"),
]):
    ax = axes[ic]

    # Build signed-unanimity matrix
    mat = np.full((len(BANDS), len(H_GRID)), np.nan)
    for ib, band in enumerate(BANDS):
        curves = []
        for pat in PATIENTS:
            if band in contrast_curves[pat]:
                curves.append(contrast_curves[pat][band])
        if curves:
            curves = np.array(curves)
            # Signed mean: +1 if all positive, -1 if all negative
            signed = np.mean(np.sign(curves), axis=0)
            mat[ib, :] = signed

    im = ax.imshow(mat, cmap=cmap, vmin=-1, vmax=1, aspect="auto",
                   extent=[H_GRID[0], H_GRID[-1], len(BANDS) - 0.5, -0.5],
                   interpolation="nearest")

    # Overlay contour lines for unanimous regions (|signed| = 1)
    for ib, band in enumerate(BANDS):
        if np.isnan(mat[ib, 0]):
            continue
        # 4/4 unanimous
        unan = np.abs(mat[ib, :]) == 1.0
        for start, end in find_regions(unan, min_width=0.01):
            h_lo, h_hi = H_GRID[start], H_GRID[end]
            sign = mat[ib, start]
            ax.plot([h_lo, h_hi], [ib, ib], linewidth=4, color="k", alpha=0.5, solid_capstyle="round")

    ax.set_yticks(range(len(BANDS)))
    ax.set_yticklabels([BAND_TEX[b] for b in BANDS], fontsize=12)
    ax.set_title(title, fontsize=11, fontweight="bold")

    if ic == 1:
        ax.set_xlabel("Dendrogram height h (fine structure → coarse structure)", fontsize=11)

    cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
    cbar.set_ticklabels(["4/4\nrecover", "3/4", "split", "3/4", "4/4\npersist"])

    # k reference
    for kt in [3, 5, 10, 20, 50]:
        idx = np.argmin(np.abs(k_at_h - kt))
        ax.axvline(H_GRID[idx], color="k", linewidth=0.5, alpha=0.3, linestyle=":")
        ax.text(H_GRID[idx], -0.45, f"k≈{kt}", fontsize=7, ha="center")

fig.suptitle(
    "Multiscale task trace: where in the hierarchy does the task leave a mark?\n"
    "Blue = task persists in rest_post | Red = brain recovers | "
    "Black bars = 4/4 patients unanimous",
    fontsize=12, y=1.03,
)
fig.tight_layout()
fig.savefig(OUTDIR / "definitive_h2_heatmap.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved definitive_h2_heatmap.pdf")


# ── FIGURE 2: Per-band detail panels with patient traces ─────────────
fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True)
axes_flat = axes.ravel()

for ib, band in enumerate(BANDS):
    ax = axes_flat[ib]
    color = BAND_COLORS[band]

    curves = []
    for ip, pat in enumerate(PATIENTS):
        if band in h2a_curves[pat]:
            c = h2a_curves[pat][band]
            curves.append(c)
            ax.plot(H_GRID, c, color=color, alpha=0.25, linewidth=0.8,
                    label=pat if ib == 0 else None)

    if curves:
        curves_arr = np.array(curves)
        mean_c = np.mean(curves_arr, axis=0)
        ax.plot(H_GRID, mean_c, color=color, linewidth=2.5, zorder=5)

        # Shade: 4/4 unanimous
        all_pos = np.all(curves_arr > 0, axis=0)
        all_neg = np.all(curves_arr < 0, axis=0)
        ylim = max(0.5, np.abs(curves_arr).max() * 1.1)
        ax.fill_between(H_GRID, -ylim, ylim, where=all_pos,
                         color="#1565C0", alpha=0.12, zorder=0, label="4/4 persist" if ib == 0 else None)
        ax.fill_between(H_GRID, -ylim, ylim, where=all_neg,
                         color="#B71C1C", alpha=0.12, zorder=0, label="4/4 recover" if ib == 0 else None)

        # Shade: 3/4 (lighter)
        n_pos = np.sum(curves_arr > 0, axis=0)
        n_neg = np.sum(curves_arr < 0, axis=0)
        ax.fill_between(H_GRID, -ylim, ylim, where=(n_pos >= 3) & ~all_pos,
                         color="#1565C0", alpha=0.05, zorder=0, label="3/4 persist" if ib == 0 else None)
        ax.fill_between(H_GRID, -ylim, ylim, where=(n_neg >= 3) & ~all_neg,
                         color="#B71C1C", alpha=0.05, zorder=0, label="3/4 recover" if ib == 0 else None)

    ax.axhline(0, color="k", linewidth=0.8, linestyle="--")
    ax.set_title(BAND_TEX[band], fontsize=13, fontweight="bold")
    if ib >= 3:
        ax.set_xlabel("Height h")
    if ib % 3 == 0:
        ax.set_ylabel("H2a contrast")
    ax.grid(alpha=0.3)

    # k reference ticks on top
    if ib < 3:
        ax2 = ax.twiny()
        for kt in [3, 5, 10, 20, 50]:
            idx = np.argmin(np.abs(k_at_h - kt))
            ax2.axvline(H_GRID[idx], color="k", linewidth=0.3, alpha=0.2, linestyle=":")
        ax2.set_xlim(ax.get_xlim())
        k_ticks_show = [3, 5, 10, 20, 50]
        h_pos = [H_GRID[np.argmin(np.abs(k_at_h - kt))] for kt in k_ticks_show]
        ax2.set_xticks(h_pos)
        ax2.set_xticklabels([f"k≈{kt}" for kt in k_ticks_show], fontsize=6)

if axes_flat[0].get_legend_handles_labels()[1]:
    axes_flat[0].legend(fontsize=6, loc="upper right")

fig.suptitle(
    "H2a per band: VI(Pre,Post) − VI(TT,Post)\n"
    "Positive = task persists | Thick = mean, thin = individual patients\n"
    "Dark shading = 4/4 unanimous, light = 3/4",
    fontsize=11, y=1.04,
)
fig.tight_layout()
fig.savefig(OUTDIR / "h2a_per_band_detail.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved h2a_per_band_detail.pdf")


# ── FIGURE 3: Raw VI(h) for the 4 key pairs, per band ────────────────
fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True, sharey=True)
axes_flat = axes.ravel()

PAIR_STYLES = {
    "task_learn↔task_test": ("#2196F3", "-", "TL↔TT"),
    "rest_pre↔rest_post": ("#4CAF50", "-", "Pre↔Post"),
    "rest_pre↔task_test": ("#9E9E9E", "--", "Pre↔TT"),
    "task_test↔rest_post": ("#E91E63", "-", "TT↔Post"),
}

for ib, band in enumerate(BANDS):
    ax = axes_flat[ib]
    for pair_key, (pcolor, lstyle, plabel) in PAIR_STYLES.items():
        all_c = []
        for pat in PATIENTS:
            if pair_key in vi_raw[pat][band]:
                all_c.append(vi_raw[pat][band][pair_key])
        if all_c:
            mean_c = np.mean(all_c, axis=0)
            ax.plot(H_GRID, mean_c, color=pcolor, linewidth=2, linestyle=lstyle,
                    label=plabel if ib == 0 else None)
            lo, hi = np.min(all_c, axis=0), np.max(all_c, axis=0)
            ax.fill_between(H_GRID, lo, hi, color=pcolor, alpha=0.08)

    ax.set_title(BAND_TEX[band], fontsize=13, fontweight="bold")
    if ib >= 3:
        ax.set_xlabel("Height h")
    if ib % 3 == 0:
        ax.set_ylabel("VI(h)")
    ax.grid(alpha=0.3)

if axes_flat[0].get_legend_handles_labels()[1]:
    axes_flat[0].legend(fontsize=8, loc="upper left")

fig.suptitle(
    "Raw VI(h) for key phase pairs (mean ± range across 4 patients)\n"
    "Lower = more similar at that hierarchy level",
    fontsize=11, y=1.02,
)
fig.tight_layout()
fig.savefig(OUTDIR / "vi_raw_all_pairs.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved vi_raw_all_pairs.pdf")


# ── Summary table ─────────────────────────────────────────────────────
print("\n" + "=" * 80)
print("DEFINITIVE MULTISCALE H2 SUMMARY")
print("=" * 80)
print("\nContiguous regions with ≥3/4 patients agreeing (width ≥ 0.02):")
print(f"{'Band':>10s}  {'Type':>8s}  {'h range':>16s}  {'k range':>10s}  {'width':>6s}  {'unanimity':>10s}")
print("-" * 70)

for contrast_name, contrast_dict in [("H2a", h2a_curves), ("H2b", h2b_curves)]:
    for band in BANDS:
        curves = []
        for pat in PATIENTS:
            if band in contrast_dict[pat]:
                curves.append(contrast_dict[pat][band])
        if not curves:
            continue
        curves_arr = np.array(curves)
        n_pos = np.sum(curves_arr > 0, axis=0)
        n_neg = np.sum(curves_arr < 0, axis=0)

        for threshold, unan_label in [(4, "4/4"), (3, "3/4")]:
            for sign_label, mask in [
                ("PERSIST", n_pos >= threshold),
                ("RECOVER", n_neg >= threshold),
            ]:
                # Exclude already-covered regions for 3/4
                if threshold == 3:
                    if sign_label == "PERSIST":
                        mask = mask & ~(n_pos >= 4)
                    else:
                        mask = mask & ~(n_neg >= 4)
                regions = find_regions(mask, min_width=0.02)
                for start, end in regions:
                    h_lo, h_hi = H_GRID[start], H_GRID[end]
                    k_hi, k_lo = k_at_h[start], k_at_h[end]
                    width = h_hi - h_lo
                    print(
                        f"{BAND_TEX[band]:>10s}  {sign_label:>8s}  "
                        f"[{h_lo:.2f}, {h_hi:.2f}]  "
                        f"k=[{k_lo:>3d},{k_hi:>3d}]  "
                        f"{width:>5.2f}  "
                        f"{unan_label:>10s}  "
                        f"({contrast_name})"
                    )

print("\nDone!")
