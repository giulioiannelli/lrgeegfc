#!/usr/bin/env python3
"""Band × Phase reorganization analysis using logCosine of D(τ).

Two variability axes:
  1. BAND: which frequency bands reorganize vs stay stable?
  2. PHASE PAIR: between which specific phases does reorganization occur?

Three complementary analyses:
  A) Z-score sign consistency (strict + relaxed)
  B) Rank consistency: is the ordering of pairs consistent across patients?
  C) Meaningful contrasts: within-type vs cross-type, task vs rest stability

All results must be patient-independent (consistent across all patients).
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from scipy.spatial.distance import squareform
from scipy.stats import rankdata

from lrg_eegfc.config import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.config.paths import LRG_CACHE, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result

# ---------------------------------------------------------------------------
PATIENT_PHASES = {
    "Pat_02": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_03": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_05": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "rest_post"],
    "Pat_08": ["rest_pre", "task_learn", "task_test", "rest_post"],
}
PATIENTS_ALL = list(PATIENT_PHASES.keys())
PATIENTS_4PH = PATIENTS_4PHASE
PATIENTS_TASK = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
BANDS = BRAIN_BANDS_NAMES
OUT_DIR = FIGURES_ROOT / "metric_exploration"

# All possible phase pairs (ordered for 4-phase patients)
ALL_PAIRS_4PH = [
    ("rest_pre", "rest_post"),
    ("task_learn", "task_test"),
    ("rest_pre", "task_learn"),
    ("rest_pre", "task_test"),
    ("task_learn", "rest_post"),
    ("task_test", "rest_post"),
]

# Pairs available for Pat_07 too
COMMON_PAIRS = [
    ("rest_pre", "rest_post"),
    ("rest_pre", "task_learn"),
    ("task_learn", "rest_post"),
]

PAIR_SHORT = {
    ("rest_pre", "rest_post"): "rest_pre↔rest_post",
    ("task_learn", "task_test"): "tLearn↔tTest",
    ("rest_pre", "task_learn"): "rest_pre↔tLearn",
    ("rest_pre", "task_test"): "rest_pre↔tTest",
    ("task_learn", "rest_post"): "tLearn↔rest_post",
    ("task_test", "rest_post"): "tTest↔rest_post",
}

PAIR_TYPE = {
    ("rest_pre", "rest_post"): "REST↔REST",
    ("task_learn", "task_test"): "TASK↔TASK",
    ("rest_pre", "task_learn"): "REST↔TASK",
    ("rest_pre", "task_test"): "REST↔TASK",
    ("task_learn", "rest_post"): "TASK↔REST",
    ("task_test", "rest_post"): "TASK↔REST",
}

# Category assignment for each pair
PAIR_CATEGORY = {
    ("rest_pre", "rest_post"): "within",
    ("task_learn", "task_test"): "within",
    ("rest_pre", "task_learn"): "cross",
    ("rest_pre", "task_test"): "cross",
    ("task_learn", "rest_post"): "cross",
    ("task_test", "rest_post"): "cross",
}


def band_label(b):
    return BRAIN_BAND_TEX_DICT.get(b, b)


def log_cosine(U1, U2):
    lU1 = np.log(np.clip(U1, 1e-15, None))
    lU2 = np.log(np.clip(U2, 1e-15, None))
    return np.dot(lU1, lU2) / max(np.linalg.norm(lU1) * np.linalg.norm(lU2), 1e-12)


# ===================================================================
# Data loading
# ===================================================================
def load_ultrametrics():
    data = {}
    for pat, phases in PATIENT_PHASES.items():
        for band in BANDS:
            for ph in phases:
                lrg = load_lrg_result(pat, ph, band, "msc", cache_root=LRG_CACHE)
                if lrg is None:
                    continue
                U = lrg.ultrametric_matrix
                if U.ndim == 2:
                    U = squareform(U)
                data[(pat, ph, band)] = U
    return data


def compute_all_sims(data):
    sims = {}
    for pat, phases in PATIENT_PHASES.items():
        for band in BANDS:
            for p1, p2 in combinations(phases, 2):
                k1, k2 = (pat, p1, band), (pat, p2, band)
                if k1 in data and k2 in data:
                    sims[(pat, band, p1, p2)] = log_cosine(data[k1], data[k2])
    return sims


# ===================================================================
# Analysis A: Z-score consistency (strict + relaxed)
# ===================================================================
def zscore_analysis(sims, patients, pairs):
    """Z-score within each (patient, band), then check sign consistency."""
    zscores = {}
    for pat in patients:
        for band in BANDS:
            vals = {}
            for p1, p2 in pairs:
                key = (pat, band, p1, p2)
                if key in sims:
                    vals[(p1, p2)] = sims[key]
            if len(vals) < 2:
                continue
            arr = np.array(list(vals.values()))
            mu, sigma = arr.mean(), arr.std()
            if sigma < 1e-10:
                for (p1, p2) in vals:
                    zscores[(pat, band, p1, p2)] = 0.0
            else:
                for (p1, p2), v in vals.items():
                    zscores[(pat, band, p1, p2)] = (v - mu) / sigma

    results = {}
    for band in BANDS:
        for pair in pairs:
            p1, p2 = pair
            pat_zs = []
            for pat in patients:
                key = (pat, band, p1, p2)
                if key in zscores:
                    pat_zs.append((pat, zscores[key]))
            if not pat_zs:
                continue
            z_vals = [z for _, z in pat_zs]
            n_pos = sum(z > 0 for z in z_vals)
            n_neg = sum(z < 0 for z in z_vals)
            N = len(z_vals)
            results[(band, pair)] = {
                "mean_z": np.mean(z_vals),
                "std_z": np.std(z_vals),
                "n_pos": n_pos,
                "n_neg": n_neg,
                "N": N,
                "unanimous": n_pos == N or n_neg == N,
                "relaxed": max(n_pos, n_neg) >= N - 1,  # allow 1 dissenter
                "direction": "similar" if n_pos > n_neg else "different",
                "agreement": max(n_pos, n_neg) / N,
                "patient_zs": pat_zs,
            }
    return results


# ===================================================================
# Analysis B: Rank consistency
# ===================================================================
def rank_analysis(sims, patients, pairs):
    """For each (patient, band), rank pairs from most similar (1) to least (N).
    Then check if rankings are consistent across patients."""
    all_ranks = {}  # (pat, band) → {pair: rank}
    for pat in patients:
        for band in BANDS:
            vals = {}
            for pair in pairs:
                p1, p2 = pair
                key = (pat, band, p1, p2)
                if key in sims:
                    vals[pair] = sims[key]
            if len(vals) < 2:
                continue
            # Rank: higher similarity = rank 1 (descending)
            sorted_pairs = sorted(vals.items(), key=lambda x: -x[1])
            ranks = {pair: rank + 1 for rank, (pair, _) in enumerate(sorted_pairs)}
            all_ranks[(pat, band)] = ranks

    results = {}
    for band in BANDS:
        for pair in pairs:
            ranks_for_pair = []
            for pat in patients:
                key = (pat, band)
                if key in all_ranks and pair in all_ranks[key]:
                    ranks_for_pair.append((pat, all_ranks[key][pair]))
            if not ranks_for_pair:
                continue
            rank_vals = [r for _, r in ranks_for_pair]
            results[(band, pair)] = {
                "mean_rank": np.mean(rank_vals),
                "std_rank": np.std(rank_vals),
                "min_rank": min(rank_vals),
                "max_rank": max(rank_vals),
                "patient_ranks": ranks_for_pair,
                "N": len(rank_vals),
                "always_top2": all(r <= 2 for r in rank_vals),
                "always_bottom2": all(r >= len(pairs) - 1 for r in rank_vals),
            }
    return results


# ===================================================================
# Analysis C: Meaningful contrasts
# ===================================================================
def contrast_analysis(sims, patients):
    """Compute scientifically meaningful contrasts for each (patient, band).

    Contrasts (4-phase patients only):
      C1: within_vs_cross = mean(within_pairs) - mean(cross_pairs)
      C2: task_stability = sim(tLearn,tTest) - sim(rest_pre,rest_post)
      C3: task_vs_cross = sim(tLearn,tTest) - mean(cross_pairs)
      C4: rest_vs_cross = sim(rest_pre,rest_post) - mean(cross_pairs)
    """
    within_pairs = [p for p in ALL_PAIRS_4PH if PAIR_CATEGORY[p] == "within"]
    cross_pairs = [p for p in ALL_PAIRS_4PH if PAIR_CATEGORY[p] == "cross"]

    contrasts_per_pat = {}  # (pat, band, contrast_name) → value
    for pat in patients:
        for band in BANDS:
            within_vals = []
            cross_vals = []
            pair_vals = {}
            for pair in ALL_PAIRS_4PH:
                p1, p2 = pair
                key = (pat, band, p1, p2)
                if key not in sims:
                    break
                pair_vals[pair] = sims[key]
                if PAIR_CATEGORY[pair] == "within":
                    within_vals.append(sims[key])
                else:
                    cross_vals.append(sims[key])
            else:
                # All pairs present
                mean_within = np.mean(within_vals)
                mean_cross = np.mean(cross_vals)
                tt = pair_vals[("task_learn", "task_test")]
                rr = pair_vals[("rest_pre", "rest_post")]

                contrasts_per_pat[(pat, band, "within_vs_cross")] = mean_within - mean_cross
                contrasts_per_pat[(pat, band, "task_stability")] = tt - mean_cross
                contrasts_per_pat[(pat, band, "rest_stability")] = rr - mean_cross
                contrasts_per_pat[(pat, band, "task_minus_rest")] = tt - rr

    # Aggregate across patients
    contrast_names = ["within_vs_cross", "task_stability", "rest_stability", "task_minus_rest"]
    results = {}
    for band in BANDS:
        for cname in contrast_names:
            vals = []
            pat_vals = []
            for pat in patients:
                key = (pat, band, cname)
                if key in contrasts_per_pat:
                    vals.append(contrasts_per_pat[key])
                    pat_vals.append((pat, contrasts_per_pat[key]))
            if not vals:
                continue
            arr = np.array(vals)
            n_pos = sum(arr > 0)
            N = len(arr)
            results[(band, cname)] = {
                "mean": arr.mean(),
                "std": arr.std(),
                "n_pos": n_pos,
                "N": N,
                "unanimous": n_pos == N or n_pos == 0,
                "relaxed": max(n_pos, N - n_pos) >= N - 1,
                "direction": "positive" if n_pos > N - n_pos else "negative",
                "agreement": max(n_pos, N - n_pos) / N,
                "patient_vals": pat_vals,
            }
    return results


# ===================================================================
# Band reorganization strength
# ===================================================================
def band_reorg_strength(sims, patients, pairs):
    strengths = {}
    for pat in patients:
        for band in BANDS:
            vals = []
            for p1, p2 in pairs:
                if (pat, band, p1, p2) in sims:
                    vals.append(sims[(pat, band, p1, p2)])
            if len(vals) >= 2:
                strengths[(pat, band)] = np.std(vals)
    return strengths


# ===================================================================
# FIGURE 1: Main heatmap (z-score, with relaxed consistency)
# ===================================================================
def fig_main_heatmap(zscore_results, pairs, title_suffix=""):
    n_bands = len(BANDS)
    n_pairs = len(pairs)

    fig, ax = plt.subplots(figsize=(max(10, n_pairs * 1.8), 7))

    mat = np.full((n_bands, n_pairs), np.nan)
    for i, band in enumerate(BANDS):
        for j, pair in enumerate(pairs):
            info = zscore_results.get((band, pair))
            if info:
                mat[i, j] = info["mean_z"]

    vmax = max(1.5, np.nanmax(np.abs(mat)))
    im = ax.imshow(mat, cmap='RdBu', vmin=-vmax, vmax=vmax, aspect='auto')

    for i in range(n_bands):
        for j in range(n_pairs):
            info = zscore_results.get((BANDS[i], pairs[j]))
            if info is None:
                ax.text(j, i, "N/A", ha='center', va='center', fontsize=7, color='gray')
                continue

            z = info["mean_z"]
            N = info["N"]

            # Markers for consistency level
            if info["unanimous"]:
                marker = "★"
                border_color = 'gold'
                border_width = 3.0
            elif info["relaxed"]:
                marker = "◆"
                border_color = 'orange'
                border_width = 2.0
            else:
                marker = ""
                border_color = None
                border_width = 0

            txt = f"{marker}\nz={z:+.2f}\n{info['n_pos']}+/{info['n_neg']}−"
            color = 'white' if abs(z) > vmax * 0.4 else 'black'
            weight = 'bold' if info["unanimous"] else 'normal'
            ax.text(j, i, txt, ha='center', va='center', fontsize=7,
                    color=color, fontweight=weight)

            if border_color:
                rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                      linewidth=border_width, edgecolor=border_color,
                                      facecolor='none')
                ax.add_patch(rect)

    ax.set_xticks(range(n_pairs))
    pair_labels = [f"{PAIR_SHORT[p]}\n({PAIR_TYPE[p]})" for p in pairs]
    ax.set_xticklabels(pair_labels, fontsize=8)
    ax.set_yticks(range(n_bands))
    ax.set_yticklabels([band_label(b) for b in BANDS], fontsize=10)

    plt.colorbar(im, ax=ax, shrink=0.8, label="z-score (+ = more similar, − = less similar)")
    ax.set_title(f"Band × Phase-pair: z-score heatmap{title_suffix}\n"
                 "★ gold = unanimous  ◆ orange = N−1 agree  Blue=stable  Red=reorganized",
                 fontsize=11)
    fig.tight_layout()
    return fig


# ===================================================================
# FIGURE 2: Rank consistency heatmap
# ===================================================================
def fig_rank_heatmap(rank_results, pairs, n_patients, title_suffix=""):
    n_bands = len(BANDS)
    n_pairs = len(pairs)
    N = n_pairs

    fig, axes = plt.subplots(1, 2, figsize=(max(12, n_pairs * 2.4), 7),
                              gridspec_kw={"width_ratios": [3, 1]})
    ax_rank = axes[0]
    ax_std = axes[1]

    mat_rank = np.full((n_bands, n_pairs), np.nan)
    mat_std = np.full((n_bands, n_pairs), np.nan)
    for i, band in enumerate(BANDS):
        for j, pair in enumerate(pairs):
            info = rank_results.get((band, pair))
            if info:
                mat_rank[i, j] = info["mean_rank"]
                mat_std[i, j] = info["std_rank"]

    # Mean rank heatmap (1 = most similar, N = least)
    im = ax_rank.imshow(mat_rank, cmap='RdYlGn_r', vmin=1, vmax=N, aspect='auto')
    for i in range(n_bands):
        for j in range(n_pairs):
            info = rank_results.get((BANDS[i], pairs[j]))
            if info is None:
                ax_rank.text(j, i, "N/A", ha='center', va='center', fontsize=7, color='gray')
                continue
            mr = info["mean_rank"]
            sr = info["std_rank"]
            ranks_str = ",".join(str(r) for _, r in info["patient_ranks"])
            txt = f"{mr:.1f}±{sr:.1f}\n[{ranks_str}]"
            color = 'white' if mr < 2 or mr > N - 1 else 'black'
            weight = 'bold' if sr < 0.5 else 'normal'  # very consistent
            ax_rank.text(j, i, txt, ha='center', va='center', fontsize=6,
                         color=color, fontweight=weight)

            # Border: low std = consistent
            if sr < 0.6:
                rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                      linewidth=2.5, edgecolor='gold', facecolor='none')
                ax_rank.add_patch(rect)
            elif sr < 1.0:
                rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                      linewidth=1.5, edgecolor='orange', facecolor='none')
                ax_rank.add_patch(rect)

    ax_rank.set_xticks(range(n_pairs))
    pair_labels = [f"{PAIR_SHORT[p]}\n({PAIR_TYPE[p]})" for p in pairs]
    ax_rank.set_xticklabels(pair_labels, fontsize=8)
    ax_rank.set_yticks(range(n_bands))
    ax_rank.set_yticklabels([band_label(b) for b in BANDS], fontsize=10)
    plt.colorbar(im, ax=ax_rank, shrink=0.8, label="Mean rank (1=most similar)")
    ax_rank.set_title(f"Rank consistency{title_suffix}\n"
                      "Gold border = std<0.6, Orange = std<1.0",
                      fontsize=10)

    # Right panel: rank spread per band (mean of std across pairs)
    band_mean_std = [np.nanmean(mat_std[i, :]) for i in range(n_bands)]
    colors_bar = plt.cm.YlOrRd(Normalize(vmin=0, vmax=max(band_mean_std) * 1.2)(band_mean_std))
    ax_std.barh(range(n_bands), band_mean_std, color=colors_bar, edgecolor='black', linewidth=0.5)
    ax_std.set_yticks(range(n_bands))
    ax_std.set_yticklabels([])
    ax_std.set_xlabel("Mean rank std\n(lower = more consistent)")
    ax_std.set_title("Rank\nspread", fontsize=9)
    ax_std.invert_yaxis()

    fig.suptitle(f"Pair ranking consistency across patients{title_suffix}", fontsize=12, y=1.01)
    fig.tight_layout()
    return fig


# ===================================================================
# FIGURE 3: Contrast analysis
# ===================================================================
def fig_contrasts(contrast_results, patients, title_suffix=""):
    contrast_names = ["within_vs_cross", "task_stability", "rest_stability", "task_minus_rest"]
    contrast_labels = [
        "Within − Cross\n(category effect)",
        "Task↔Task − Cross\n(task stability)",
        "Rest↔Rest − Cross\n(rest stability)",
        "Task↔Task − Rest↔Rest\n(task > rest?)",
    ]

    fig, axes = plt.subplots(1, 4, figsize=(20, 6), sharey=True)

    for c_idx, (cname, clabel) in enumerate(zip(contrast_names, contrast_labels)):
        ax = axes[c_idx]

        for b_idx, band in enumerate(BANDS):
            info = contrast_results.get((band, cname))
            if info is None:
                continue

            # Plot individual patients
            for pat, val in info["patient_vals"]:
                ax.scatter(val, b_idx, marker='o', s=30, alpha=0.6,
                           color='#1976D2' if val > 0 else '#D32F2F',
                           edgecolors='black', linewidths=0.3)

            # Plot mean
            m = info["mean"]
            ax.scatter(m, b_idx, marker='D', s=80, color='black', zorder=10)

            # Mark consistency
            if info["unanimous"]:
                ax.annotate("★", (m, b_idx), textcoords="offset points",
                            xytext=(12, 0), fontsize=12, color='gold', fontweight='bold',
                            va='center')
            elif info["relaxed"]:
                ax.annotate("◆", (m, b_idx), textcoords="offset points",
                            xytext=(12, 0), fontsize=10, color='orange',
                            va='center')

        ax.axvline(0, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
        ax.set_title(clabel, fontsize=9, fontweight='bold')
        ax.set_xlabel("Contrast value (logCosine D(τ))")
        if c_idx == 0:
            ax.set_yticks(range(len(BANDS)))
            ax.set_yticklabels([band_label(b) for b in BANDS], fontsize=10)

    fig.suptitle(f"Reorganization contrasts{title_suffix}\n"
                 "Dots = patients, Diamond = mean. ★ = unanimous, ◆ = N−1 agree",
                 fontsize=12, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    return fig


# ===================================================================
# FIGURE 4: Band reorganization strength
# ===================================================================
def fig_band_strength(strengths, patients, title_suffix=""):
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(BANDS))
    width = 0.12

    for i, pat in enumerate(patients):
        vals = [strengths.get((pat, band), 0) for band in BANDS]
        offset = (i - len(patients) / 2 + 0.5) * width
        ax.bar(x + offset, vals, width * 0.9,
               label=pat.replace("Pat_0", "P"), alpha=0.8)

    means = [np.mean([strengths.get((pat, band), 0) for pat in patients])
             for band in BANDS]
    ax.plot(x, means, 'k-o', linewidth=2, markersize=6, label='Mean', zorder=10)

    ax.set_xticks(x)
    ax.set_xticklabels([band_label(b) for b in BANDS], fontsize=10)
    ax.set_ylabel("Reorganization strength\n(std of pairwise similarities)")
    ax.set_title(f"Per-band reorganization strength{title_suffix}\n"
                 "(high = strong differentiation between phase pairs)",
                 fontsize=11)
    ax.legend(fontsize=8, ncol=len(patients))
    fig.tight_layout()
    return fig


# ===================================================================
# FIGURE 5: Raw similarity per-pair profiles
# ===================================================================
def fig_pair_profiles(sims, patients, pairs, zscore_results, title_suffix=""):
    n_pairs = len(pairs)
    ncols = min(3, n_pairs)
    nrows = (n_pairs + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows), sharey=True)
    if nrows * ncols == 1:
        axes_flat = [axes]
    else:
        axes_flat = np.array(axes).ravel()

    for p_idx, pair in enumerate(pairs):
        ax = axes_flat[p_idx]
        p1, p2 = pair

        for pat in patients:
            vals = [sims.get((pat, band, p1, p2), np.nan) for band in BANDS]
            ax.plot(range(len(BANDS)), vals, 'o-', alpha=0.5, markersize=4,
                    label=pat.replace("Pat_0", "P"))

        # Highlight consistent/relaxed bands
        for b_idx, band in enumerate(BANDS):
            info = zscore_results.get((band, pair))
            if info:
                if info["unanimous"]:
                    c = '#1976D2' if info["direction"] == "similar" else '#D32F2F'
                    ax.axvspan(b_idx - 0.35, b_idx + 0.35, alpha=0.2, color=c)
                elif info["relaxed"]:
                    c = '#1976D2' if info["direction"] == "similar" else '#D32F2F'
                    ax.axvspan(b_idx - 0.35, b_idx + 0.35, alpha=0.08, color=c)

        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels([band_label(b) for b in BANDS], fontsize=8, rotation=30)
        ax.set_title(f"{PAIR_SHORT[pair]} ({PAIR_TYPE[pair]})", fontsize=10)
        if p_idx % ncols == 0:
            ax.set_ylabel("logCosine D(τ)")
        ax.legend(fontsize=6, ncol=2)

    for idx in range(len(pairs), len(axes_flat)):
        axes_flat[idx].set_visible(False)

    fig.suptitle(f"Per-pair similarity profiles across bands{title_suffix}\n"
                 "(shaded = consistent across patients: blue=stable, red=reorganized)",
                 fontsize=12, fontweight='bold')
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    return fig


# ===================================================================
# Print comprehensive summary
# ===================================================================
def print_summary(zscore_results, rank_results, contrast_results,
                  strengths, patients, pairs, label):
    print(f"\n{'='*80}")
    print(f"SUMMARY: {label}")
    print(f"{'='*80}")

    # --- Z-score consistency ---
    n_unan = sum(1 for info in zscore_results.values() if info.get("unanimous"))
    n_relax = sum(1 for info in zscore_results.values() if info.get("relaxed"))
    n_total = len(zscore_results)
    print(f"\nZ-score consistency: {n_unan} unanimous, {n_relax} relaxed (N−1), out of {n_total}")

    print(f"\n{'Band':<14}", end="")
    for pair in pairs:
        print(f"  {PAIR_SHORT[pair]:>14}", end="")
    print(f"  {'Reorg STD':>10}")

    for band in BANDS:
        mean_s = np.mean([strengths.get((pat, band), 0) for pat in patients])
        print(f"{band_label(band):<14}", end="")
        for pair in pairs:
            info = zscore_results.get((band, pair))
            if info is None:
                print(f"  {'N/A':>14}", end="")
            else:
                z = info["mean_z"]
                if info["unanimous"]:
                    print(f"  {'★ z=' + f'{z:+.2f}':>14}", end="")
                elif info["relaxed"]:
                    print(f"  {'◆ z=' + f'{z:+.2f}':>14}", end="")
                else:
                    frac = f"{info['n_pos']}+/{info['n_neg']}−"
                    print(f"  {frac:>14}", end="")
        print(f"  {mean_s:>10.4f}")

    # --- Rank consistency ---
    if rank_results:
        print(f"\n--- Rank consistency ---")
        print(f"{'Band':<14}", end="")
        for pair in pairs:
            print(f"  {PAIR_SHORT[pair]:>14}", end="")
        print()
        for band in BANDS:
            print(f"{band_label(band):<14}", end="")
            for pair in pairs:
                info = rank_results.get((band, pair))
                if info is None:
                    print(f"  {'N/A':>14}", end="")
                else:
                    mr = info["mean_rank"]
                    sr = info["std_rank"]
                    flag = "★" if sr < 0.6 else ("◆" if sr < 1.0 else " ")
                    print(f"  {flag}{mr:.1f}±{sr:.1f}".rjust(14), end="")
            print()

    # --- Contrasts ---
    if contrast_results:
        print(f"\n--- Contrasts ---")
        contrast_names = ["within_vs_cross", "task_stability", "rest_stability", "task_minus_rest"]
        for cname in contrast_names:
            print(f"\n  {cname}:")
            for band in BANDS:
                info = contrast_results.get((band, cname))
                if info is None:
                    continue
                flag = "★" if info["unanimous"] else ("◆" if info["relaxed"] else " ")
                sign = "+" if info["mean"] > 0 else "−"
                pat_str = " | ".join(f"{p}:{v:+.4f}" for p, v in info["patient_vals"])
                print(f"    {flag} {band_label(band):<14} mean={info['mean']:+.5f} "
                      f"({info['n_pos']}/{info['N']} pos)  [{pat_str}]")

    # --- Final interpretation ---
    print(f"\n{'='*80}")
    print("CONSISTENT FINDINGS:")
    print(f"{'='*80}")

    # Z-score unanimous
    for band in BANDS:
        for pair in pairs:
            info = zscore_results.get((band, pair))
            if info and info["unanimous"]:
                arrow = "STABLE" if info["direction"] == "similar" else "REORGANIZED"
                print(f"  ★ {band_label(band):>14} × {PAIR_SHORT[pair]:<20} → {arrow} "
                      f"(z={info['mean_z']:+.2f}, N={info['N']})")

    # Z-score relaxed (not unanimous)
    relax_only = []
    for band in BANDS:
        for pair in pairs:
            info = zscore_results.get((band, pair))
            if info and info["relaxed"] and not info["unanimous"]:
                relax_only.append((band, pair, info))
    if relax_only:
        print(f"\n  --- Relaxed (N−1 agree) ---")
        for band, pair, info in relax_only:
            arrow = "STABLE" if info["direction"] == "similar" else "REORGANIZED"
            print(f"  ◆ {band_label(band):>14} × {PAIR_SHORT[pair]:<20} → {arrow} "
                  f"(z={info['mean_z']:+.2f}, {max(info['n_pos'],info['n_neg'])}/{info['N']})")

    # Contrasts unanimous
    if contrast_results:
        print(f"\n  --- Contrasts ---")
        contrast_labels = {
            "within_vs_cross": "Within > Cross",
            "task_stability": "Task stable",
            "rest_stability": "Rest stable",
            "task_minus_rest": "Task > Rest",
        }
        for cname in ["within_vs_cross", "task_stability", "rest_stability", "task_minus_rest"]:
            for band in BANDS:
                info = contrast_results.get((band, cname))
                if info and info["unanimous"]:
                    sign = "YES" if info["mean"] > 0 else "NO"
                    print(f"  ★ {band_label(band):>14} {contrast_labels[cname]:<20} → {sign} "
                          f"(Δ={info['mean']:+.5f}, N={info['N']})")


# ===================================================================
# Main
# ===================================================================
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading ultrametric data...")
    data = load_ultrametrics()
    print(f"  Loaded {len(data)} entries")

    print("Computing pairwise logCosine D(τ)...")
    sims = compute_all_sims(data)
    print(f"  Computed {len(sims)} pairs")

    # ================================================================
    # Analysis 1: 4-phase patients, all 6 pairs
    # ================================================================
    print("\n" + "=" * 80)
    print("ANALYSIS 1: 4-phase patients (Pat_02,03,05,08), 6 pairs")
    print("=" * 80)

    zs_4 = zscore_analysis(sims, PATIENTS_4PH, ALL_PAIRS_4PH)
    rk_4 = rank_analysis(sims, PATIENTS_4PH, ALL_PAIRS_4PH)
    ct_4 = contrast_analysis(sims, PATIENTS_4PH)
    st_4 = band_reorg_strength(sims, PATIENTS_4PH, ALL_PAIRS_4PH)

    print_summary(zs_4, rk_4, ct_4, st_4, PATIENTS_4PH, ALL_PAIRS_4PH, "4-phase patients")

    fig1 = fig_main_heatmap(zs_4, ALL_PAIRS_4PH, " (4-phase patients)")
    fig1.savefig(OUT_DIR / "band_phase_zscore_4ph.pdf", bbox_inches='tight')
    plt.close(fig1)

    fig2 = fig_rank_heatmap(rk_4, ALL_PAIRS_4PH, len(PATIENTS_4PH), " (4-phase patients)")
    fig2.savefig(OUT_DIR / "band_phase_ranks_4ph.pdf", bbox_inches='tight')
    plt.close(fig2)

    fig3 = fig_contrasts(ct_4, PATIENTS_4PH, " (4-phase patients)")
    fig3.savefig(OUT_DIR / "band_contrasts_4ph.pdf", bbox_inches='tight')
    plt.close(fig3)

    fig4 = fig_band_strength(st_4, PATIENTS_4PH, " (4-phase patients)")
    fig4.savefig(OUT_DIR / "band_reorg_strength_4ph.pdf", bbox_inches='tight')
    plt.close(fig4)

    fig5 = fig_pair_profiles(sims, PATIENTS_4PH, ALL_PAIRS_4PH, zs_4, " (4-phase patients)")
    fig5.savefig(OUT_DIR / "pair_profiles_4ph.pdf", bbox_inches='tight')
    plt.close(fig5)

    print("  Saved: band_phase_zscore_4ph.pdf, band_phase_ranks_4ph.pdf,")
    print("         band_contrasts_4ph.pdf, band_reorg_strength_4ph.pdf, pair_profiles_4ph.pdf")

    # ================================================================
    # Analysis 2: All task patients, 3 common pairs
    # ================================================================
    print("\n" + "=" * 80)
    print("ANALYSIS 2: All task patients (02,03,05,07,08), 3 common pairs")
    print("=" * 80)

    zs_a = zscore_analysis(sims, PATIENTS_TASK, COMMON_PAIRS)
    rk_a = rank_analysis(sims, PATIENTS_TASK, COMMON_PAIRS)
    st_a = band_reorg_strength(sims, PATIENTS_TASK, COMMON_PAIRS)

    print_summary(zs_a, rk_a, {}, st_a, PATIENTS_TASK, COMMON_PAIRS,
                  "All task patients (3 common pairs)")

    fig6 = fig_main_heatmap(zs_a, COMMON_PAIRS, " (all task patients)")
    fig6.savefig(OUT_DIR / "band_phase_zscore_all.pdf", bbox_inches='tight')
    plt.close(fig6)

    fig7 = fig_rank_heatmap(rk_a, COMMON_PAIRS, len(PATIENTS_TASK), " (all task patients)")
    fig7.savefig(OUT_DIR / "band_phase_ranks_all.pdf", bbox_inches='tight')
    plt.close(fig7)

    fig8 = fig_band_strength(st_a, PATIENTS_TASK, " (all task patients)")
    fig8.savefig(OUT_DIR / "band_reorg_strength_all.pdf", bbox_inches='tight')
    plt.close(fig8)

    fig9 = fig_pair_profiles(sims, PATIENTS_TASK, COMMON_PAIRS, zs_a, " (all task patients)")
    fig9.savefig(OUT_DIR / "pair_profiles_all.pdf", bbox_inches='tight')
    plt.close(fig9)

    print("  Saved: band_phase_zscore_all.pdf, band_phase_ranks_all.pdf,")
    print("         band_reorg_strength_all.pdf, pair_profiles_all.pdf")

    # ================================================================
    # Pat_06 control + raw data dump
    # ================================================================
    print("\n\n=== Pat_06 (rest only control) ===")
    for band in BANDS:
        key = ("Pat_06", band, "rest_pre", "rest_post")
        if key in sims:
            print(f"  {band_label(band):<14} rest_pre↔rest_post logCosine = {sims[key]:.6f}")

    # Dump raw data to CSV
    rows = []
    for (pat, band, p1, p2), val in sorted(sims.items()):
        rows.append({
            "patient": pat, "band": band, "phase1": p1, "phase2": p2,
            "logCosine": val, "pair_type": PAIR_CATEGORY.get((p1, p2), "?"),
        })
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "logcosine_dtau_results.csv", index=False)
    print(f"\n  Saved logcosine_dtau_results.csv ({len(df)} rows)")

    print("\nDone!")


if __name__ == "__main__":
    main()
