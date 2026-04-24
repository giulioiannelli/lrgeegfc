#!/usr/bin/env python3
"""Continuous multiscale comparison of ALL phase pairs.

Uses cophenetic threshold height h as the continuous scale variable.
At each h, cutting the dendrogram gives a partition — VI(h) measures
how different two partitions are.

h is intrinsic to the ultrametric — it means the same thing for all trees.
h near 0 = fine communities (many clusters), h near 1 = coarse (few clusters).

Computes VI(h) for ALL 6 phase pairs × 6 bands × 4 patients.

Key figures:
  1. Per-band panels showing all 6 pairs as VI(h) curves
  2. Band × h heatmaps for each pair (similarity maps)
  3. Band × h heatmaps for key contrasts
  4. Fine vs coarse decomposition showing where effects live

Output: data/figures/metric_exploration/continuous_all_pairs/
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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

OUTDIR = FIGURES_ROOT / "metric_exploration" / "continuous_all_pairs" / FC_METHOD
OUTDIR.mkdir(parents=True, exist_ok=True)

from lrg_eegfc.config.const import PATIENTS_4PHASE
PATIENTS_4PH = PATIENTS_4PHASE
ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]

PATIENT_PHASES = {
    "Pat_02": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_03": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_05": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "rest_post"],
    "Pat_08": ["rest_pre", "task_learn", "task_test", "rest_post"],
}

BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = {
    "delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
    "beta": r"$\beta$", "low_gamma": r"$\gamma_l$", "high_gamma": r"$\gamma_h$",
}
BAND_COLORS = {
    "delta": "#1f77b4", "theta": "#ff7f0e", "alpha": "#2ca02c",
    "beta": "#d62728", "low_gamma": "#9467bd", "high_gamma": "#8c564b",
}

ALL_PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
ALL_PAIRS = [
    ("task_learn", "task_test"),   # within task
    ("rest_pre", "rest_post"),         # within rest
    ("rest_pre", "task_learn"),      # cross
    ("rest_pre", "task_test"),       # cross
    ("task_learn", "rest_post"),     # cross
    ("task_test", "rest_post"),      # cross (task trace)
]
# Use PATIENTS as alias for 4-phase patients in contrast functions (backward compat)
PATIENTS = PATIENTS_4PH
PAIR_LABELS = {
    ("task_learn", "task_test"): "TL↔TT",
    ("rest_pre", "rest_post"): "Pre↔Post",
    ("rest_pre", "task_learn"): "Pre↔TL",
    ("rest_pre", "task_test"): "Pre↔TT",
    ("task_learn", "rest_post"): "TL↔Post",
    ("task_test", "rest_post"): "TT↔Post",
}
PAIR_COLORS = {
    ("task_learn", "task_test"): "#1565C0",   # blue — within task
    ("rest_pre", "rest_post"): "#2E7D32",         # green — within rest
    ("rest_pre", "task_learn"): "#9E9E9E",      # grey
    ("rest_pre", "task_test"): "#BDBDBD",       # light grey
    ("task_learn", "rest_post"): "#E65100",     # orange
    ("task_test", "rest_post"): "#C62828",      # red — task trace
}
PAIR_LINESTYLES = {
    ("task_learn", "task_test"): "-",
    ("rest_pre", "rest_post"): "-",
    ("rest_pre", "task_learn"): "--",
    ("rest_pre", "task_test"): "--",
    ("task_learn", "rest_post"): "-.",
    ("task_test", "rest_post"): "-",
}

# Log-spaced height grid: the LRG diffusion distances produce a log-scaled
# hierarchy, so we sample logarithmically to give uniform resolution across scales
H_GRID = np.geomspace(0.003, 0.995, 400)


from lrg_eegfc.utils.metrics import compute_vi  # canonical implementation


# ── Load linkage matrices ─────────────────────────────────────────────
print("Loading LRG results...")
linkages = {}
for pat in ALL_PATIENTS:
    for phase in PATIENT_PHASES[pat]:
        for band in BANDS:
            try:
                res = load_lrg_result(pat, phase, band, fc_method=FC_METHOD,
                                      cache_root=CACHE_ROOT)
                linkages[(pat, band, phase)] = res.linkage_matrix
            except Exception:
                pass
print(f"Loaded {len(linkages)} linkages")

# ── Compute k(h) per (patient, band, phase) for reference ────────────
# Show the ACTUAL k at each h for each tree (not a single reference)
print("Computing k(h) profiles...")
k_profiles = {}
for pat in ALL_PATIENTS:
    for band in BANDS:
        for phase in PATIENT_PHASES[pat]:
            key = (pat, band, phase)
            if key not in linkages:
                continue
            Z = linkages[key]
            ks = np.array([len(np.unique(fcluster(Z, h, criterion="distance")))
                           for h in H_GRID])
            k_profiles[key] = ks

# Mean k(h) for annotation
all_k = np.array([k_profiles[k] for k in k_profiles])
mean_k_at_h = np.mean(all_k, axis=0)

# ── Compute VI(h) for all combinations ────────────────────────────────
print("Computing VI(h) for all pairs (this takes a few minutes)...")
# vi_h[pat][band][(p1,p2)] = array of VI values at each h
vi_h = {pat: {band: {} for band in BANDS} for pat in ALL_PATIENTS}

# Compute all available pairs for each patient
total = 0
for pat in ALL_PATIENTS:
    phases = PATIENT_PHASES[pat]
    from itertools import combinations
    for p1, p2 in combinations(phases, 2):
        total += len(BANDS)

done = 0
for pat in ALL_PATIENTS:
    phases = PATIENT_PHASES[pat]
    from itertools import combinations
    pat_pairs = list(combinations(phases, 2))
    for band in BANDS:
        for p1, p2 in pat_pairs:
            k1, k2 = (pat, band, p1), (pat, band, p2)
            if k1 not in linkages or k2 not in linkages:
                done += 1
                continue
            Z1, Z2 = linkages[k1], linkages[k2]
            vals = np.zeros(len(H_GRID))
            for ih, h in enumerate(H_GRID):
                l1 = fcluster(Z1, h, criterion="distance")
                l2 = fcluster(Z2, h, criterion="distance")
                vals[ih] = compute_vi(l1, l2)
            vi_h[pat][band][(p1, p2)] = vals
            done += 1
    print(f"  {pat} done ({done}/{total})")


# ── Normalize: VI(h) / log(mean_k(h)) to get NVI ∈ [0, ~2] ──────────
# This normalizes by the "capacity" at each scale: fine scale can have
# more entropy so raw VI is naturally larger
# Actually use VI / log(k_mean) where k_mean is the mean k at that h
# This way values are comparable across scales
norm_factor = np.log(np.maximum(mean_k_at_h, 1.001))  # log(k), avoid log(1)=0


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 1: Per-band, all 6 pairs as curves
# ═══════════════════════════════════════════════════════════════════════
print("\nFigure 1: per-band all-pairs curves...")

fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True)
axes_flat = axes.ravel()

for ib, band in enumerate(BANDS):
    ax = axes_flat[ib]

    for pair in ALL_PAIRS:
        curves = []
        for pat in ALL_PATIENTS:
            v = vi_h[pat][band].get(pair)
            if v is not None:
                curves.append(v)
        if not curves:
            continue
        mean_c = np.mean(curves, axis=0)
        # Normalize by log(k) at each scale
        mean_c_norm = mean_c / np.maximum(norm_factor, 0.01)

        n_pat = len(curves)
        lbl = f"{PAIR_LABELS[pair]} (n={n_pat})"
        ax.plot(H_GRID, mean_c_norm,
                color=PAIR_COLORS[pair],
                linestyle=PAIR_LINESTYLES[pair],
                linewidth=2.5 if pair in [("task_learn", "task_test"), ("rest_pre", "rest_post"),
                                           ("task_test", "rest_post")] else 1.5,
                label=lbl,
                alpha=0.9)

        # Light envelope
        lo = np.min(curves, axis=0) / np.maximum(norm_factor, 0.01)
        hi = np.max(curves, axis=0) / np.maximum(norm_factor, 0.01)
        ax.fill_between(H_GRID, lo, hi, color=PAIR_COLORS[pair], alpha=0.06)

    ax.set_title(BAND_TEX[band], fontsize=14, fontweight="bold")
    ax.set_xscale("log")
    if ib >= 3:
        ax.set_xlabel("Height h (fine → coarse)")
    if ib % 3 == 0:
        ax.set_ylabel("NVI = VI(h) / log(k)")
    ax.grid(alpha=0.3)
    if ib == 0:
        ax.legend(fontsize=7, loc="upper right")

    # k reference on top
    if ib < 3:
        ax2 = ax.twiny()
        ax2.set_xscale("log")
        ax2.set_xlim(ax.get_xlim())
        k_ticks = [100, 50, 20, 10, 5, 3, 2]
        h_pos = [H_GRID[np.argmin(np.abs(mean_k_at_h - kt))] for kt in k_ticks]
        ax2.set_xticks(h_pos)
        ax2.set_xticklabels([f"{kt}" for kt in k_ticks], fontsize=7)
        ax2.set_xlabel("≈ k (communities)", fontsize=8)

fig.suptitle(
    "All phase pairs across the full hierarchy — NVI = VI / log(k)\n"
    "Normalized so values are comparable across scales (higher = more different)\n"
    "Blue = TL↔TT  |  Green = Pre↔Post  |  Red = TT↔Post  |  Grey = cross",
    fontsize=12, y=1.04,
)
fig.tight_layout()
fig.savefig(OUTDIR / "all_pairs_per_band.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved all_pairs_per_band.pdf")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 2: Band × h heatmaps — one per pair
# Shows the mean NVI for each (band, h) — where are phases similar/different?
# ═══════════════════════════════════════════════════════════════════════
print("Figure 2: band × h heatmaps per pair...")

fig, axes = plt.subplots(2, 3, figsize=(18, 9))
axes_flat = axes.ravel()

cmap_vi = LinearSegmentedColormap.from_list(
    "vi", ["#1B5E20", "#A5D6A7", "#FFF9C4", "#EF9A9A", "#B71C1C"], N=256)

vmax_global = 0
for pair in ALL_PAIRS:
    for band in BANDS:
        curves = [vi_h[pat][band].get(pair) for pat in ALL_PATIENTS
                  if vi_h[pat][band].get(pair) is not None]
        if curves:
            mn = np.mean(curves, axis=0) / np.maximum(norm_factor, 0.01)
            vmax_global = max(vmax_global, np.max(mn))

for ip, pair in enumerate(ALL_PAIRS):
    ax = axes_flat[ip]
    mat = np.full((len(BANDS), len(H_GRID)), np.nan)
    n_pats_for_pair = 0
    for ib, band in enumerate(BANDS):
        curves = [vi_h[pat][band].get(pair) for pat in ALL_PATIENTS
                  if vi_h[pat][band].get(pair) is not None]
        if curves:
            mat[ib, :] = np.mean(curves, axis=0) / np.maximum(norm_factor, 0.01)
            n_pats_for_pair = max(n_pats_for_pair, len(curves))

    # Use pcolormesh for proper log-scale x-axis
    Y = np.arange(len(BANDS))
    im = ax.pcolormesh(H_GRID, Y, mat, cmap=cmap_vi, vmin=0, vmax=vmax_global * 0.8,
                        shading="nearest")
    ax.set_xscale("log")

    ax.set_yticks(range(len(BANDS)))
    ax.set_yticklabels([BAND_TEX[b] for b in BANDS], fontsize=10)
    ax.set_ylim(len(BANDS) - 0.5, -0.5)
    ax.set_xlabel("h (fine → coarse)" if ip >= 3 else "")
    lbl = PAIR_LABELS[pair]
    ptype = "within" if pair in [("task_learn", "task_test"), ("rest_pre", "rest_post")] else "cross"
    ax.set_title(f"{lbl} ({ptype}, n={n_pats_for_pair})", fontsize=11, fontweight="bold",
                 color=PAIR_COLORS[pair])

    # k annotations
    for kt in [50, 20, 10, 5, 3]:
        idx = np.argmin(np.abs(mean_k_at_h - kt))
        ax.axvline(H_GRID[idx], color="k", linewidth=0.3, alpha=0.3, linestyle=":")
        if ip < 3:
            ax.text(H_GRID[idx], -0.45, f"k≈{kt}", fontsize=6, ha="center")

plt.colorbar(im, ax=axes_flat[-1], shrink=0.7, label="NVI (higher = more different)")

fig.suptitle(
    "Phase pair distance across scales — NVI = VI(h)/log(k)\n"
    "Green = similar | Red = different",
    fontsize=13, y=1.02,
)
fig.tight_layout()
fig.savefig(OUTDIR / "heatmaps_per_pair.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved heatmaps_per_pair.pdf")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 3: Contrast heatmaps — the key comparisons
# ═══════════════════════════════════════════════════════════════════════
print("Figure 3: contrast heatmaps...")

CONTRASTS = [
    ("H1: TL↔TT is closest\n(n varies by patient)",
     lambda vi_b: _relative_gain(vi_b, ("task_learn", "task_test"))),
    ("H2a: TT↔Post < Pre↔Post\n(task persists in rest_post)",
     lambda vi_b: _relative_diff(vi_b, ("rest_pre", "rest_post"), ("task_test", "rest_post"))),
    ("H2a': TT↔Post < Pre↔TL\n(rest_post ≈ task, not rest_pre≈task)",
     lambda vi_b: _relative_diff(vi_b, ("rest_pre", "task_learn"), ("task_test", "rest_post"))),
    ("TL↔Post vs Pre↔Post\n(task trace via TL)",
     lambda vi_b: _relative_diff(vi_b, ("rest_pre", "rest_post"), ("task_learn", "rest_post"))),
    ("TL↔TT vs Pre↔Post\n(task stable vs rest stable)",
     lambda vi_b: _relative_diff(vi_b, ("rest_pre", "rest_post"), ("task_learn", "task_test"))),
    ("Within vs Cross gap\n(H3)",
     lambda vi_b: _within_cross_gap(vi_b)),
]


def _relative_diff(vi_band_data, pair_a, pair_b, patients=None):
    """(VI_a - VI_b) / mean(VI_a, VI_b). Positive = b is more similar."""
    if patients is None:
        patients = ALL_PATIENTS
    per_pat = []
    for pat in patients:
        if pat not in vi_band_data:
            continue
        va = vi_band_data[pat].get(pair_a)
        vb = vi_band_data[pat].get(pair_b)
        if va is not None and vb is not None:
            denom = (va + vb) / 2
            with np.errstate(divide="ignore", invalid="ignore"):
                rel = np.where(denom > 1e-10, (va - vb) / denom, 0.0)
            per_pat.append(rel)
    if not per_pat:
        return np.full(len(H_GRID), np.nan), np.full(len(H_GRID), np.nan), 0
    arr = np.array(per_pat)
    return np.mean(arr, axis=0), np.mean(np.sign(arr), axis=0), len(per_pat)


def _relative_gain(vi_band_data, target_pair, patients=None):
    """(mean_others - target) / mean_others. Positive = target is closest."""
    if patients is None:
        patients = ALL_PATIENTS
    per_pat = []
    for pat in patients:
        if pat not in vi_band_data:
            continue
        target = vi_band_data[pat].get(target_pair)
        if target is None:
            continue
        # Use all available pairs for this patient
        others = [vi_band_data[pat][p] for p in vi_band_data[pat]
                  if p != target_pair]
        if not others:
            continue
        mo = np.mean(others, axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            g = np.where(mo > 1e-10, (mo - target) / mo, 0.0)
        per_pat.append(g)
    if not per_pat:
        return np.full(len(H_GRID), np.nan), np.full(len(H_GRID), np.nan), 0
    arr = np.array(per_pat)
    return np.mean(arr, axis=0), np.mean(np.sign(arr), axis=0), len(per_pat)


def _within_cross_gap(vi_band_data, patients=None):
    """(mean_cross - mean_within) / mean_all. Positive = H3 holds."""
    if patients is None:
        patients = ALL_PATIENTS
    within = [("task_learn", "task_test"), ("rest_pre", "rest_post")]
    cross = [p for p in ALL_PAIRS if p not in within]
    per_pat = []
    for pat in patients:
        if pat not in vi_band_data:
            continue
        wv = [vi_band_data[pat][p] for p in within if p in vi_band_data[pat]]
        cv = [vi_band_data[pat][p] for p in cross if p in vi_band_data[pat]]
        if wv and cv:
            wm = np.mean(wv, axis=0)
            cm = np.mean(cv, axis=0)
            denom = (cm + wm) / 2
            with np.errstate(divide="ignore", invalid="ignore"):
                g = np.where(denom > 1e-10, (cm - wm) / denom, 0.0)
            per_pat.append(g)
    if not per_pat:
        return np.full(len(H_GRID), np.nan), np.full(len(H_GRID), np.nan), 0
    arr = np.array(per_pat)
    return np.mean(arr, axis=0), np.mean(np.sign(arr), axis=0), len(per_pat)


fig, axes = plt.subplots(2, 3, figsize=(18, 10))
axes_flat = axes.ravel()

cmap_contrast = LinearSegmentedColormap.from_list(
    "contrast", ["#B71C1C", "#EF9A9A", "#FAFAFA", "#90CAF9", "#0D47A1"], N=256)

for ic, (title, func) in enumerate(CONTRASTS):
    ax = axes_flat[ic]

    mat_val = np.full((len(BANDS), len(H_GRID)), np.nan)
    mat_unan = np.full((len(BANDS), len(H_GRID)), np.nan)

    n_pats_contrast = 0
    for ib, band in enumerate(BANDS):
        # Build per-band dict with ALL patients
        vi_band = {pat: vi_h[pat][band] for pat in ALL_PATIENTS}
        val, unan, n_p = func(vi_band)
        mat_val[ib, :] = val
        mat_unan[ib, :] = unan
        n_pats_contrast = max(n_pats_contrast, n_p)

    # Use signed unanimity for color (stronger signal)
    vmax_c = 1.0
    Y_centers = np.arange(len(BANDS))
    im = ax.pcolormesh(H_GRID, Y_centers, mat_unan, cmap=cmap_contrast,
                        vmin=-vmax_c, vmax=vmax_c, shading="nearest")
    ax.set_xscale("log")
    ax.set_ylim(len(BANDS) - 0.5, -0.5)

    # Overlay: thick border where ALL patients unanimous (|unan| == 1)
    for ib in range(len(BANDS)):
        unanimous = np.abs(mat_unan[ib, :]) >= 0.99  # == 1.0 with float tolerance
        if np.any(unanimous):
            # Find contiguous regions
            starts = np.where(np.diff(np.concatenate(([0], unanimous.astype(int)))) == 1)[0]
            ends = np.where(np.diff(np.concatenate((unanimous.astype(int), [0]))) == -1)[0]
            for s, e in zip(starts, ends):
                if H_GRID[e] - H_GRID[s] > 0.005:  # skip tiny fragments
                    rect = plt.Rectangle(
                        (H_GRID[s], ib - 0.45), H_GRID[e] - H_GRID[s], 0.9,
                        linewidth=2, edgecolor="black", facecolor="none", zorder=10)
                    ax.add_patch(rect)

    ax.set_yticks(range(len(BANDS)))
    ax.set_yticklabels([BAND_TEX[b] for b in BANDS], fontsize=10)
    ax.set_xlabel("h (fine → coarse)" if ic >= 3 else "")
    ax.set_title(title, fontsize=9, fontweight="bold")

    # k annotations
    for kt in [50, 20, 10, 5, 3]:
        idx = np.argmin(np.abs(mean_k_at_h - kt))
        ax.axvline(H_GRID[idx], color="k", linewidth=0.3, alpha=0.4, linestyle=":")
        if ic < 3:
            ax.text(H_GRID[idx], -0.45, f"k≈{kt}", fontsize=6, ha="center")

cbar = fig.colorbar(im, ax=axes_flat, shrink=0.5, pad=0.02, aspect=30,
                     location="bottom")
cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
cbar.set_ticklabels(["All negative", "Majority neg", "Split", "Majority pos", "All positive"])
cbar.set_label("Signed unanimity across patients (black box = all agree)", fontsize=10)

fig.suptitle(
    "Multiscale contrasts — all hypotheses, all scales\n"
    "Blue = positive (hypothesis holds) | Red = negative | "
    "Black box = all 4 patients agree",
    fontsize=12, y=1.02,
)
fig.tight_layout()
fig.savefig(OUTDIR / "contrast_heatmaps.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved contrast_heatmaps.pdf")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 4: k(h) profile — what does "fine" mean?
# ═══════════════════════════════════════════════════════════════════════
print("Figure 4: k(h) profile...")

fig, ax = plt.subplots(figsize=(10, 4))
# Show range across all trees
all_k_vals = np.array([k_profiles[k] for k in k_profiles])
mean_k = np.mean(all_k_vals, axis=0)
lo_k = np.min(all_k_vals, axis=0)
hi_k = np.max(all_k_vals, axis=0)

ax.plot(H_GRID, mean_k, "k-", linewidth=2, label="Mean across all trees")
ax.fill_between(H_GRID, lo_k, hi_k, color="gray", alpha=0.2, label="Range")
ax.set_xscale("log")
ax.set_xlabel("Cophenetic height h (log scale)", fontsize=11)
ax.set_ylabel("Number of clusters k", fontsize=11)
ax.set_title("How many communities at each scale?\nThis maps h → k for interpreting the heatmaps",
             fontsize=12)
ax.legend()
ax.set_yscale("log")
ax.grid(alpha=0.3)
ax.set_ylim(1, 120)

# Mark key h values
for kt in [2, 3, 5, 10, 20, 50, 100]:
    idx = np.argmin(np.abs(mean_k - kt))
    ax.axhline(kt, color="gray", linewidth=0.5, alpha=0.3)
    ax.text(H_GRID[-1] + 0.01, kt, f"k={kt}", fontsize=8, va="center")

fig.tight_layout()
fig.savefig(OUTDIR / "k_at_h_profile.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved k_at_h_profile.pdf")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 5: Pat_06 and Pat_07 — individual profiles
# ═══════════════════════════════════════════════════════════════════════
print("Figure 5: Pat_06 and Pat_07 individual profiles...")

for pat, phases in [("Pat_06", ["rest_pre", "rest_post"]),
                     ("Pat_07", ["rest_pre", "task_learn", "rest_post"])]:
    from itertools import combinations as combs
    pat_pairs = list(combs(phases, 2))
    n_pairs = len(pat_pairs)

    fig, axes = plt.subplots(2, 3, figsize=(18, 9), sharex=True)
    axes_flat = axes.ravel()

    for ib, band in enumerate(BANDS):
        ax = axes_flat[ib]
        for p1, p2 in pat_pairs:
            v = vi_h[pat][band].get((p1, p2))
            if v is None:
                continue
            v_norm = v / np.maximum(norm_factor, 0.01)
            short_p1 = {"rest_pre": "Pre", "task_learn": "TL", "task_test": "TT", "rest_post": "Post"}
            lbl = f"{short_p1[p1]}↔{short_p1[p2]}"
            # Color: use matching pair color if available
            pair_key = (p1, p2)
            color = PAIR_COLORS.get(pair_key, "#333333")
            ax.plot(H_GRID, v_norm, color=color, linewidth=2.5, label=lbl)

        # Also overlay the 4-patient mean for Pre↔Post as reference
        ref_curves = [vi_h[p][band].get(("rest_pre", "rest_post"))
                      for p in PATIENTS_4PH
                      if ("rest_pre", "rest_post") in vi_h[p][band]]
        if ref_curves:
            ref_mean = np.mean(ref_curves, axis=0) / np.maximum(norm_factor, 0.01)
            ax.plot(H_GRID, ref_mean, color="#2E7D32", linewidth=1, linestyle=":",
                    alpha=0.5, label="Pre↔Post (4-pat mean)")

        ax.set_title(BAND_TEX[band], fontsize=14, fontweight="bold")
        ax.set_xscale("log")
        if ib >= 3:
            ax.set_xlabel("Height h (log scale)")
        if ib % 3 == 0:
            ax.set_ylabel("NVI = VI(h) / log(k)")
        ax.grid(alpha=0.3)
        if ib == 0:
            ax.legend(fontsize=8, loc="upper right")

    phase_str = ", ".join(phases)
    fig.suptitle(
        f"{pat} — Available pairs ({phase_str})\n"
        f"NVI(h) for {n_pairs} pair(s) | Dotted = 4-patient Pre↔Post reference",
        fontsize=12, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(OUTDIR / f"{pat}_profile.pdf", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Saved {pat}_profile.pdf")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE 6: Pre↔Post comparison across ALL 6 patients
# ═══════════════════════════════════════════════════════════════════════
print("Figure 6: Pre↔Post all 6 patients...")

fig, axes = plt.subplots(2, 3, figsize=(18, 9), sharex=True)
axes_flat = axes.ravel()

PAT_LINE_COLORS = {
    "Pat_02": "#1f77b4", "Pat_03": "#ff7f0e", "Pat_05": "#2ca02c",
    "Pat_06": "#d62728", "Pat_07": "#9467bd", "Pat_08": "#8c564b",
}

for ib, band in enumerate(BANDS):
    ax = axes_flat[ib]
    for pat in ALL_PATIENTS:
        v = vi_h[pat][band].get(("rest_pre", "rest_post"))
        if v is None:
            continue
        v_norm = v / np.maximum(norm_factor, 0.01)
        ax.plot(H_GRID, v_norm, color=PAT_LINE_COLORS[pat], linewidth=2,
                label=pat, alpha=0.8)

    ax.set_title(BAND_TEX[band], fontsize=14, fontweight="bold")
    ax.set_xscale("log")
    if ib >= 3:
        ax.set_xlabel("Height h (log scale)")
    if ib % 3 == 0:
        ax.set_ylabel("NVI Pre↔Post")
    ax.grid(alpha=0.3)
    if ib == 0:
        ax.legend(fontsize=8)

fig.suptitle(
    "Pre↔Post distance across all 6 patients\n"
    "NVI(h) = VI / log(k) — higher = more reorganization between rest phases",
    fontsize=12, y=1.02,
)
fig.tight_layout()
fig.savefig(OUTDIR / "pre_post_all_patients.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved pre_post_all_patients.pdf")


# ═══════════════════════════════════════════════════════════════════════
# Print summary
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("KEY FINDINGS FROM ALL-PAIRS CONTINUOUS ANALYSIS")
print("=" * 80)

print("\nH2a unanimous regions — VI(Pre,Post) > VI(TT,Post):")
for band in BANDS:
    vi_band = {pat: vi_h[pat][band] for pat in ALL_PATIENTS}
    val, unan, n_p = _relative_diff(vi_band, ("rest_pre", "rest_post"), ("task_test", "rest_post"))
    # Find h ranges where unan == 1 (persist) or -1 (recover)
    persist_h = H_GRID[unan >= 0.99]
    recover_h = H_GRID[unan <= -0.99]

    parts = []
    if len(persist_h) > 0:
        # Find contiguous ranges
        groups = np.split(persist_h, np.where(np.diff(persist_h) > 0.01)[0] + 1)
        for g in groups:
            if len(g) > 1:
                k_lo = mean_k_at_h[np.argmin(np.abs(H_GRID - g[-1]))]
                k_hi = mean_k_at_h[np.argmin(np.abs(H_GRID - g[0]))]
                parts.append(f"PERSIST h=[{g[0]:.2f},{g[-1]:.2f}] k≈[{k_lo:.0f},{k_hi:.0f}]")
    if len(recover_h) > 0:
        groups = np.split(recover_h, np.where(np.diff(recover_h) > 0.01)[0] + 1)
        for g in groups:
            if len(g) > 1:
                k_lo = mean_k_at_h[np.argmin(np.abs(H_GRID - g[-1]))]
                k_hi = mean_k_at_h[np.argmin(np.abs(H_GRID - g[0]))]
                parts.append(f"RECOVER h=[{g[0]:.2f},{g[-1]:.2f}] k≈[{k_lo:.0f},{k_hi:.0f}]")
    print(f"  {BAND_TEX[band]:>12s}: {' | '.join(parts) if parts else 'none'}")

print("\nH1 unanimous regions — TL↔TT relative gain > 0:")
for band in BANDS:
    vi_band = {pat: vi_h[pat][band] for pat in ALL_PATIENTS}
    val, unan, n_p = _relative_gain(vi_band, ("task_learn", "task_test"))
    pos_h = H_GRID[unan >= 0.99]
    parts = []
    if len(pos_h) > 0:
        groups = np.split(pos_h, np.where(np.diff(pos_h) > 0.01)[0] + 1)
        for g in groups:
            if len(g) > 1:
                k_lo = mean_k_at_h[np.argmin(np.abs(H_GRID - g[-1]))]
                k_hi = mean_k_at_h[np.argmin(np.abs(H_GRID - g[0]))]
                parts.append(f"h=[{g[0]:.2f},{g[-1]:.2f}] k≈[{k_lo:.0f},{k_hi:.0f}]")
    print(f"  {BAND_TEX[band]:>12s}: {' | '.join(parts) if parts else 'none'}")

print("\nDone!")
