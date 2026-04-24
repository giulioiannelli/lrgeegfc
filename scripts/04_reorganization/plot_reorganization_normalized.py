#!/usr/bin/env python3
"""Reorganization analysis with per-patient normalization.

For each patient, all pairwise distances (6 bands × 6 pairs = 36 values) are
min-max normalized to [0, 1].  This removes inter-patient scale differences
while preserving within-patient band contrasts.

Then computes RI from normalized distances and produces summary figures.
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cophenet
from scipy.stats import pearsonr

from lrg_eegfc.config import BRAIN_BANDS_NAMES, PHASE_LABELS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.config.paths import LRG_CACHE, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result

PATIENTS = PATIENTS_4PHASE
BANDS = BRAIN_BANDS_NAMES
PHASES = list(PHASE_LABELS)
OUTPUT_DIR = FIGURES_ROOT / "metric_exploration"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

BAND_TEX = [BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS]
ALL_PAIRS = list(combinations(PHASES, 2))

REST_PAIR = ("rest_pre", "rest_post")
TASK_PAIR = ("task_learn", "task_test")
WITHIN_PAIRS = [REST_PAIR, TASK_PAIR]
CROSS_PAIRS = [
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"),
    ("task_learn", "rest_post"), ("task_test", "rest_post"),
]

PAIR_LABELS = {
    ("rest_pre", "rest_post"): "rest_pre↔rest_post",
    ("rest_pre", "task_learn"): "rest_pre↔tLearn",
    ("rest_pre", "task_test"): "rest_pre↔tTest",
    ("task_learn", "task_test"): "tLearn↔tTest",
    ("task_learn", "rest_post"): "tLearn↔rest_post",
    ("task_test", "rest_post"): "tTest↔rest_post",
}

PAIR_CATEGORY = {}
for p in ALL_PAIRS:
    if set(p) == set(REST_PAIR):
        PAIR_CATEGORY[p] = "within-rest"
    elif set(p) == set(TASK_PAIR):
        PAIR_CATEGORY[p] = "within-task"
    else:
        PAIR_CATEGORY[p] = "cross"

CAT_COLORS = {
    "within-rest": "#2166ac",
    "within-task": "#4393c3",
    "cross": "#b2182b",
}


def log_cophenetic_distance(Z1, Z2):
    c1 = cophenet(Z1)
    c2 = cophenet(Z2)
    log_c1 = np.log(np.clip(c1, 1e-12, None))
    log_c2 = np.log(np.clip(c2, 1e-12, None))
    corr, _ = pearsonr(log_c1, log_c2)
    return 1.0 - corr


def load_and_compute():
    """Load all LRG, compute raw distances, then min-max normalize per patient."""
    # Load
    lrg_data = {}
    for pat in PATIENTS:
        for band in BANDS:
            for phase in PHASES:
                lrg = load_lrg_result(pat, phase, band, "msc", cache_root=LRG_CACHE)
                if lrg is not None:
                    lrg_data[(pat, phase, band)] = lrg.linkage_matrix

    # Raw distances
    raw = {}
    for pat in PATIENTS:
        for band in BANDS:
            for pair in ALL_PAIRS:
                ki = (pat, pair[0], band)
                kj = (pat, pair[1], band)
                if ki in lrg_data and kj in lrg_data:
                    raw[(pat, band, pair)] = log_cophenetic_distance(
                        lrg_data[ki], lrg_data[kj]
                    )

    # Min-max normalize per patient (across all bands and pairs)
    normed = {}
    for pat in PATIENTS:
        vals = [raw[(pat, b, p)] for b in BANDS for p in ALL_PAIRS
                if (pat, b, p) in raw]
        vmin, vmax = min(vals), max(vals)
        span = vmax - vmin if vmax > vmin else 1e-12
        for band in BANDS:
            for pair in ALL_PAIRS:
                key = (pat, band, pair)
                if key in raw:
                    normed[key] = (raw[key] - vmin) / span

    return raw, normed


def compute_ri(dists, pat, band):
    """RI = D_cross / D_within from distance dict."""
    d_within_vals = [dists.get((pat, band, p), np.nan) for p in WITHIN_PAIRS]
    d_cross_vals = [dists.get((pat, band, p), np.nan) for p in CROSS_PAIRS]
    d_within = np.nanmean(d_within_vals)
    d_cross = np.nanmean(d_cross_vals)
    if d_within == 0 or np.isnan(d_within):
        return np.nan
    return d_cross / d_within


def build_dataframe(raw, normed):
    rows = []
    for pat in PATIENTS:
        for band in BANDS:
            ri_raw = compute_ri(raw, pat, band)
            ri_norm = compute_ri(normed, pat, band)

            d_rest_raw = raw.get((pat, band, REST_PAIR), np.nan)
            d_task_raw = raw.get((pat, band, TASK_PAIR), np.nan)
            d_rest_norm = normed.get((pat, band, REST_PAIR), np.nan)
            d_task_norm = normed.get((pat, band, TASK_PAIR), np.nan)

            d_cross_raw = np.nanmean([raw.get((pat, band, p), np.nan)
                                      for p in CROSS_PAIRS])
            d_cross_norm = np.nanmean([normed.get((pat, band, p), np.nan)
                                       for p in CROSS_PAIRS])

            rows.append({
                "patient": pat, "band": band,
                "RI_raw": ri_raw, "RI_norm": ri_norm,
                "D_rest_raw": d_rest_raw, "D_task_raw": d_task_raw,
                "D_cross_raw": d_cross_raw,
                "D_rest_norm": d_rest_norm, "D_task_norm": d_task_norm,
                "D_cross_norm": d_cross_norm,
            })
    return pd.DataFrame(rows)


# ===================================================================
# Figures
# ===================================================================
def plot_main_summary(df, normed):
    fig = plt.figure(figsize=(24, 18))
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.4, wspace=0.35)

    pat_colors = {"Pat_02": "#66c2a5", "Pat_03": "#fc8d62",
                  "Pat_05": "#8da0cb", "Pat_08": "#e78ac3"}
    x = np.arange(len(BANDS))

    # -----------------------------------------------------------
    # (A) RI raw vs normalized comparison
    # -----------------------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    for pat in PATIENTS:
        pdf = df[df["patient"] == pat]
        ris = [pdf[pdf["band"] == b]["RI_raw"].values[0] for b in BANDS]
        ax.scatter(x + {"Pat_02": -0.15, "Pat_03": -0.05,
                        "Pat_05": 0.05, "Pat_08": 0.15}[pat],
                   ris, s=70, color=pat_colors[pat], alpha=0.5,
                   edgecolors="gray", linewidths=0.3)
    for j, band in enumerate(BANDS):
        vals = df[df["band"] == band]["RI_raw"].values
        ax.plot([j - 0.2, j + 0.2], [np.mean(vals)] * 2, color="black", lw=2)
    ax.axhline(1.0, color="gray", ls="--", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(BAND_TEX, fontsize=11)
    ax.set_ylabel("RI (raw)", fontsize=11)
    ax.set_title("(A) RI before normalization", fontsize=12, fontweight="bold")
    ax.grid(alpha=0.2, axis="y")

    ax2 = fig.add_subplot(gs[0, 1])
    for pat in PATIENTS:
        pdf = df[df["patient"] == pat]
        ris = [pdf[pdf["band"] == b]["RI_norm"].values[0] for b in BANDS]
        ax2.scatter(x + {"Pat_02": -0.15, "Pat_03": -0.05,
                         "Pat_05": 0.05, "Pat_08": 0.15}[pat],
                    ris, s=70, color=pat_colors[pat], label=pat,
                    edgecolors="black", linewidths=0.5)
    for j, band in enumerate(BANDS):
        vals = df[df["band"] == band]["RI_norm"].values
        m = np.mean(vals)
        s = np.std(vals)
        ax2.errorbar(j, m, yerr=s, fmt="s", color="black", ms=7,
                     capsize=5, capthick=2, lw=2, zorder=6)
    ax2.axhline(1.0, color="gray", ls="--", lw=1)
    ax2.set_xticks(x)
    ax2.set_xticklabels(BAND_TEX, fontsize=11)
    ax2.set_ylabel("RI (normalized)", fontsize=11)
    ax2.set_title("(B) RI after per-patient normalization", fontsize=12, fontweight="bold")
    ax2.legend(fontsize=8, loc="upper left")
    ax2.grid(alpha=0.2, axis="y")

    # -----------------------------------------------------------
    # (C) Normalized within vs cross bars
    # -----------------------------------------------------------
    ax3 = fig.add_subplot(gs[0, 2])
    width = 0.35
    w_means = [df[df["band"] == b]["D_cross_norm"].mean() for b in BANDS]
    c_means = [df[df["band"] == b]["D_cross_norm"].mean() for b in BANDS]
    wi_means = [np.mean([df[df["band"] == b]["D_rest_norm"].mean(),
                         df[df["band"] == b]["D_task_norm"].mean()])
                for b in BANDS]
    cr_means = [df[df["band"] == b]["D_cross_norm"].mean() for b in BANDS]
    wi_stds = [np.std([
        *df[df["band"] == b]["D_rest_norm"].values,
        *df[df["band"] == b]["D_task_norm"].values]) for b in BANDS]
    cr_stds = [df[df["band"] == b]["D_cross_norm"].std() for b in BANDS]

    ax3.bar(x - width / 2, wi_means, width, yerr=wi_stds,
            label="Within-cond", color="#2166ac", alpha=0.8, capsize=4)
    ax3.bar(x + width / 2, cr_means, width, yerr=cr_stds,
            label="Cross-cond", color="#b2182b", alpha=0.8, capsize=4)
    ax3.set_xticks(x)
    ax3.set_xticklabels(BAND_TEX, fontsize=11)
    ax3.set_ylabel("Normalized distance", fontsize=11)
    ax3.set_title("(C) Normalized distances", fontsize=12, fontweight="bold")
    ax3.legend(fontsize=9)
    ax3.grid(alpha=0.2, axis="y")

    # -----------------------------------------------------------
    # (D) Per-patient normalized bar charts (all 6 pairs)
    # -----------------------------------------------------------
    for idx, pat in enumerate(PATIENTS):
        ax = fig.add_subplot(gs[1 + idx // 2, idx % 2 + (1 if idx >= 2 else 0)])
        if idx < 2:
            ax = fig.add_subplot(gs[1, idx])
        else:
            ax = fig.add_subplot(gs[2, idx - 2])

        n_pairs = len(ALL_PAIRS)
        bar_width = 0.12

        for k, pair in enumerate(ALL_PAIRS):
            cat = PAIR_CATEGORY[pair]
            color = CAT_COLORS[cat]
            vals = [normed.get((pat, band, pair), np.nan) for band in BANDS]
            offset = (k - n_pairs / 2 + 0.5) * bar_width
            ax.bar(x + offset, vals, bar_width, color=color, alpha=0.85,
                   edgecolor="white", linewidth=0.5,
                   label=PAIR_LABELS[pair] if idx == 0 else None)

        ax.set_xticks(x)
        ax.set_xticklabels(BAND_TEX, fontsize=11)
        ax.set_ylabel("Norm. distance", fontsize=10)
        ax.set_title(f"({chr(68 + idx)}) {pat} — all pairs (normalized)",
                     fontsize=12, fontweight="bold")
        ax.grid(alpha=0.2, axis="y")
        ax.set_ylim(0, 1.05)

    # Legend from first patient panel
    handles = []
    labels_list = []
    for pair in ALL_PAIRS:
        cat = PAIR_CATEGORY[pair]
        h = plt.Rectangle((0, 0), 1, 1, color=CAT_COLORS[cat], alpha=0.85)
        handles.append(h)
        labels_list.append(PAIR_LABELS[pair])
    fig.legend(handles, labels_list, loc="lower center", ncol=6,
               fontsize=10, bbox_to_anchor=(0.5, -0.01))

    fig.suptitle(
        "Cross-Phase Structural Reorganization — per-patient normalized (1-LogCoph)\n"
        "Each patient's distances min-max scaled to [0, 1] before computing RI",
        fontsize=15, fontweight="bold", y=1.01,
    )
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    path = OUTPUT_DIR / "reorganization_normalized_summary.pdf"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


def plot_normalized_heatmaps(normed):
    """Per-patient 4x4 heatmaps on normalized scale."""
    for pat in PATIENTS:
        fig, axes = plt.subplots(2, 3, figsize=(20, 12))
        axes_flat = axes.ravel()

        for j, band in enumerate(BANDS):
            ax = axes_flat[j]
            matrix = np.zeros((4, 4))
            for i1, p1 in enumerate(PHASES):
                for i2, p2 in enumerate(PHASES):
                    if i1 < i2:
                        pair = (p1, p2)
                        val = normed.get((pat, band, pair), np.nan)
                        matrix[i1, i2] = val
                        matrix[i2, i1] = val

            mask = np.eye(4, dtype=bool)
            display = np.where(mask, np.nan, matrix)
            im = ax.imshow(display, cmap="YlOrRd", vmin=0, vmax=1, aspect="equal")
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.06, shrink=0.85)

            for i1 in range(4):
                for i2 in range(4):
                    if i1 != i2:
                        val = matrix[i1, i2]
                        color = "white" if val > 0.55 else "black"
                        ax.text(i2, i1, f"{val:.2f}", ha="center", va="center",
                                fontsize=11, fontweight="bold", color=color)

            ax.set_xticks(range(4))
            ax.set_yticks(range(4))
            ax.set_xticklabels(PHASES, rotation=45, ha="right", fontsize=9)
            ax.set_yticklabels(PHASES, fontsize=9)
            ax.set_title(BRAIN_BAND_TEX_DICT.get(band, band),
                         fontsize=13, fontweight="bold")

        fig.suptitle(
            f"{pat} — Normalized phase distances (1-LogCoph, [0-1] per patient)\n"
            "All heatmaps share [0, 1] scale — directly comparable across bands",
            fontsize=14, fontweight="bold", y=1.01,
        )
        plt.tight_layout()
        path = OUTPUT_DIR / f"heatmaps_normalized_{pat}.pdf"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        fig.savefig(path.with_suffix(".png"), dpi=100, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {path}")


def main():
    print("Loading and computing distances...")
    raw, normed = load_and_compute()
    df = build_dataframe(raw, normed)

    csv_path = OUTPUT_DIR / "reorganization_normalized_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved: {csv_path}")

    print("\nGenerating figures...")
    plot_main_summary(df, normed)
    plot_normalized_heatmaps(normed)

    # Summary table
    print("\n" + "=" * 70)
    print("NORMALIZED RI — summary")
    print("=" * 70)
    pivot = df.pivot(index="patient", columns="band", values="RI_norm")
    pivot = pivot[BANDS]
    print("\nRI per patient/band:")
    print(pivot.round(3).to_string())
    print(f"\n{'Band':<12} {'Mean':>7} {'Std':>7} {'CV':>7} {'Min':>7} {'Max':>7} {'All>1':>6}")
    print("-" * 58)
    for band in BANDS:
        vals = pivot[band].values
        m, s = np.mean(vals), np.std(vals)
        cv = s / m if m > 0 else 0
        print(f"{band:<12} {m:>7.3f} {s:>7.3f} {cv:>7.3f} {min(vals):>7.3f} "
              f"{max(vals):>7.3f} {str(np.all(vals > 1)):>6}")


if __name__ == "__main__":
    main()
