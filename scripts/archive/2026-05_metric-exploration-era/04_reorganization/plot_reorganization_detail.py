#!/usr/bin/env python3
"""Detailed reorganization analysis for the best metric (1-LogCoph).

Produces:
  1. Summary figure: RI by band (all patients overlaid) + group bar chart
  2. Per-patient detail: 6 pairwise distances for all bands, phase distance
     heatmaps, RI profile
  3. Grand summary: what is universal vs patient-specific

Data source: data/figures/metric_exploration/metric_exploration_results.csv
Plus recomputes full 4x4 distance matrices for heatmap panels.
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
from scipy.spatial.distance import squareform
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

# Phase pair categories
REST_PAIR = ("rest_pre", "rest_post")
TASK_PAIR = ("task_learn", "task_test")
CROSS_PAIRS = [
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"),
    ("task_learn", "rest_post"), ("task_test", "rest_post"),
]
ALL_PAIRS = list(combinations(PHASES, 2))

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
    """1 - Pearson(log(coph1), log(coph2))."""
    c1 = cophenet(Z1)
    c2 = cophenet(Z2)
    log_c1 = np.log(np.clip(c1, 1e-12, None))
    log_c2 = np.log(np.clip(c2, 1e-12, None))
    corr, _ = pearsonr(log_c1, log_c2)
    return 1.0 - corr


def load_all_lrg():
    data = {}
    for pat in PATIENTS:
        for band in BANDS:
            for phase in PHASES:
                lrg = load_lrg_result(pat, phase, band, "msc", cache_root=LRG_CACHE)
                if lrg is not None:
                    data[(pat, phase, band)] = lrg.linkage_matrix
    return data


def compute_full_distances(lrg_data):
    """Compute all pairwise distances → dict[(pat, band, pair)] = distance."""
    dists = {}
    for pat in PATIENTS:
        for band in BANDS:
            for pi, pj in ALL_PAIRS:
                ki = (pat, pi, band)
                kj = (pat, pj, band)
                if ki in lrg_data and kj in lrg_data:
                    d = log_cophenetic_distance(lrg_data[ki], lrg_data[kj])
                    dists[(pat, band, (pi, pj))] = d
    return dists


# ===================================================================
# Figure 1: Grand summary — RI by band with patient variability
# ===================================================================
def plot_grand_summary(csv_path):
    df = pd.read_csv(csv_path)
    mdf = df[df["metric"] == "1-LogCoph"].copy()

    fig, axes = plt.subplots(1, 3, figsize=(22, 7), width_ratios=[1.2, 1, 1])

    # --- Panel A: RI dot plot by band ---
    ax = axes[0]
    x = np.arange(len(BANDS))
    pat_colors = {"Pat_02": "#66c2a5", "Pat_03": "#fc8d62",
                  "Pat_05": "#8da0cb", "Pat_08": "#e78ac3"}
    offsets = {"Pat_02": -0.15, "Pat_03": -0.05, "Pat_05": 0.05, "Pat_08": 0.15}

    for pat in PATIENTS:
        pdf = mdf[mdf["patient"] == pat]
        ris = [pdf[pdf["band"] == b]["RI"].values[0] for b in BANDS]
        ax.scatter(x + offsets[pat], ris, s=100, color=pat_colors[pat],
                   label=pat, zorder=5, edgecolors="black", linewidths=0.5)

    # Band means + error bars
    for j, band in enumerate(BANDS):
        vals = mdf[mdf["band"] == band]["RI"].values
        mean_ri = np.mean(vals)
        std_ri = np.std(vals)
        ax.errorbar(j, mean_ri, yerr=std_ri, fmt="s", color="black",
                    markersize=8, capsize=6, capthick=2, linewidth=2, zorder=6)

    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1.5, alpha=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(BAND_TEX, fontsize=13)
    ax.set_ylabel("Reorganization Index (RI)", fontsize=12)
    ax.set_title("(A) RI by frequency band\n(all patients)", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(alpha=0.2, axis="y")
    ax.set_ylim(0.7, None)

    # --- Panel B: D_within vs D_cross grouped bars ---
    ax = axes[1]
    width = 0.35
    d_within_means = []
    d_cross_means = []
    d_within_stds = []
    d_cross_stds = []
    for band in BANDS:
        bdf = mdf[mdf["band"] == band]
        d_within_means.append(bdf["D_within"].mean())
        d_cross_means.append(bdf["D_cross"].mean())
        d_within_stds.append(bdf["D_within"].std())
        d_cross_stds.append(bdf["D_cross"].std())

    ax.bar(x - width / 2, d_within_means, width, yerr=d_within_stds,
           label="Within-condition", color="#2166ac", alpha=0.8, capsize=4)
    ax.bar(x + width / 2, d_cross_means, width, yerr=d_cross_stds,
           label="Cross-condition", color="#b2182b", alpha=0.8, capsize=4)
    ax.set_xticks(x)
    ax.set_xticklabels(BAND_TEX, fontsize=13)
    ax.set_ylabel("Mean distance (1 - LogCoph corr)", fontsize=11)
    ax.set_title("(B) Within vs cross-condition\ndistances", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(alpha=0.2, axis="y")

    # --- Panel C: D_rest vs D_task decomposition ---
    ax = axes[2]
    d_rest_means = [mdf[mdf["band"] == b]["D_rest"].mean() for b in BANDS]
    d_task_means = [mdf[mdf["band"] == b]["D_task"].mean() for b in BANDS]
    d_rest_stds = [mdf[mdf["band"] == b]["D_rest"].std() for b in BANDS]
    d_task_stds = [mdf[mdf["band"] == b]["D_task"].std() for b in BANDS]

    ax.bar(x - width / 2, d_rest_means, width, yerr=d_rest_stds,
           label="D(rest_pre↔rest_post)", color="#2166ac", alpha=0.8, capsize=4)
    ax.bar(x + width / 2, d_task_means, width, yerr=d_task_stds,
           label="D(task_learn↔task_test)", color="#4393c3", alpha=0.8, capsize=4)
    ax.set_xticks(x)
    ax.set_xticklabels(BAND_TEX, fontsize=13)
    ax.set_ylabel("Distance", fontsize=11)
    ax.set_title("(C) Within-condition decomposition\nrest vs task stability",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(alpha=0.2, axis="y")

    fig.suptitle(
        "Cross-Phase Structural Reorganization — 1-LogCoph metric\n"
        "RI = D(rest↔task) / D(within-condition)  |  RI > 1 → reorganization",
        fontsize=15, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    path = OUTPUT_DIR / "reorganization_grand_summary.pdf"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


# ===================================================================
# Figure 2: Per-patient detail — all 6 pairwise distances by band
# ===================================================================
def plot_per_patient_detail(dists):
    fig, axes = plt.subplots(2, 2, figsize=(22, 16))

    for idx, pat in enumerate(PATIENTS):
        ax = axes[idx // 2, idx % 2]

        n_bands = len(BANDS)
        n_pairs = len(ALL_PAIRS)
        x = np.arange(n_bands)
        bar_width = 0.12

        for k, pair in enumerate(ALL_PAIRS):
            cat = PAIR_CATEGORY[pair]
            color = CAT_COLORS[cat]
            label = PAIR_LABELS[pair]
            vals = [dists.get((pat, band, pair), np.nan) for band in BANDS]
            offset = (k - n_pairs / 2 + 0.5) * bar_width
            bars = ax.bar(x + offset, vals, bar_width, color=color, alpha=0.85,
                          edgecolor="white", linewidth=0.5)
            # Add label only on first band iteration
            if idx == 0:
                bars.set_label(label)

        ax.set_xticks(x)
        ax.set_xticklabels(BAND_TEX, fontsize=12)
        ax.set_ylabel("Distance (1 - LogCoph corr)", fontsize=11)
        ax.set_title(f"{pat}", fontsize=14, fontweight="bold")
        ax.grid(alpha=0.2, axis="y")
        ax.set_ylim(0, None)

    # Add legend from first axes
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, fontsize=10,
               bbox_to_anchor=(0.5, -0.02))

    fig.suptitle(
        "Per-patient pairwise phase distances — 1-LogCoph\n"
        "Blue = within-condition (rest↔rest, task↔task)  |  Red = cross-condition (rest↔task)",
        fontsize=15, fontweight="bold", y=1.01,
    )
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    path = OUTPUT_DIR / "reorganization_per_patient_detail.pdf"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


# ===================================================================
# Figure 3: Per-patient 4x4 heatmaps for each band
# ===================================================================
def plot_patient_heatmaps(dists):
    for pat in PATIENTS:
        fig, axes = plt.subplots(2, 3, figsize=(20, 12))
        axes = axes.ravel()

        for j, band in enumerate(BANDS):
            ax = axes[j]
            matrix = np.zeros((4, 4))
            for i1, p1 in enumerate(PHASES):
                for i2, p2 in enumerate(PHASES):
                    if i1 < i2:
                        pair = (p1, p2)
                        val = dists.get((pat, band, pair), np.nan)
                        matrix[i1, i2] = val
                        matrix[i2, i1] = val

            mask = np.eye(4, dtype=bool)
            display = np.where(mask, np.nan, matrix)
            im = ax.imshow(display, cmap="YlOrRd", vmin=0,
                           vmax=max(0.4, np.nanmax(display)), aspect="equal")
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.06, shrink=0.85)

            for i1 in range(4):
                for i2 in range(4):
                    if i1 != i2:
                        val = matrix[i1, i2]
                        vmax = max(0.4, np.nanmax(display))
                        color = "white" if val / vmax > 0.55 else "black"
                        ax.text(i2, i1, f"{val:.3f}", ha="center", va="center",
                                fontsize=10, fontweight="bold", color=color)

            ax.set_xticks(range(4))
            ax.set_yticks(range(4))
            ax.set_xticklabels(PHASES, rotation=45, ha="right", fontsize=9)
            ax.set_yticklabels(PHASES, fontsize=9)
            ax.set_title(f"{BRAIN_BAND_TEX_DICT.get(band, band)}",
                         fontsize=13, fontweight="bold")

        fig.suptitle(
            f"{pat} — Phase distance matrices (1-LogCoph)\n"
            "Low values (yellow) = similar structure  |  "
            "High values (red) = reorganized structure",
            fontsize=15, fontweight="bold", y=1.01,
        )
        plt.tight_layout()
        path = OUTPUT_DIR / f"heatmaps_{pat}.pdf"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        fig.savefig(path.with_suffix(".png"), dpi=100, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {path}")


# ===================================================================
# Figure 4: Universal vs variable — coefficient of variation by band
# ===================================================================
def plot_universality(dists):
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    # Panel A: Mean distance per pair type across patients (by band)
    ax = axes[0]
    pair_types = {"within-rest": [], "within-task": [], "cross": []}
    for band in BANDS:
        for ptype in pair_types:
            vals = []
            for pat in PATIENTS:
                for pair in ALL_PAIRS:
                    if PAIR_CATEGORY[pair] == ptype:
                        v = dists.get((pat, band, pair), np.nan)
                        if not np.isnan(v):
                            vals.append(v)
            pair_types[ptype].append(np.mean(vals))

    x = np.arange(len(BANDS))
    width = 0.25
    for i, (ptype, vals) in enumerate(pair_types.items()):
        ax.bar(x + (i - 1) * width, vals, width, color=CAT_COLORS[ptype],
               label=ptype, alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(BAND_TEX, fontsize=13)
    ax.set_ylabel("Mean distance", fontsize=11)
    ax.set_title("(A) Distances by pair category\n(averaged across patients)",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)
    ax.grid(alpha=0.2, axis="y")

    # Panel B: CV of RI across patients per band
    csv_path = OUTPUT_DIR / "metric_exploration_results.csv"
    df = pd.read_csv(csv_path)
    mdf = df[df["metric"] == "1-LogCoph"]

    ax = axes[1]
    means = []
    stds = []
    cvs = []
    for band in BANDS:
        bdf = mdf[mdf["band"] == band]
        m = bdf["RI"].mean()
        s = bdf["RI"].std()
        means.append(m)
        stds.append(s)
        cvs.append(s / m if m > 0 else 0)

    colors_bar = ["#d73027" if m > 1.5 else "#fc8d59" if m > 1.2
                  else "#fee090" for m in means]

    bars = ax.bar(x, means, yerr=stds, color=colors_bar, alpha=0.85,
                  capsize=6, edgecolor="black", linewidth=0.5)
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1.5)

    # Annotate with CV
    for j, (m, cv) in enumerate(zip(means, cvs)):
        ax.text(j, m + stds[j] + 0.05, f"CV={cv:.2f}",
                ha="center", va="bottom", fontsize=9, fontstyle="italic")

    ax.set_xticks(x)
    ax.set_xticklabels(BAND_TEX, fontsize=13)
    ax.set_ylabel("Reorganization Index (RI)", fontsize=11)
    ax.set_title("(B) RI with cross-patient variability\n"
                 "(error bars = std, CV annotated)",
                 fontsize=13, fontweight="bold")
    ax.grid(alpha=0.2, axis="y")
    ax.set_ylim(0, None)

    fig.suptitle(
        "Universality of reorganization across patients — 1-LogCoph",
        fontsize=15, fontweight="bold", y=1.01,
    )
    plt.tight_layout()
    path = OUTPUT_DIR / "reorganization_universality.pdf"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=100, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


def main():
    csv_path = OUTPUT_DIR / "metric_exploration_results.csv"

    print("Loading LRG data for full pairwise distances...")
    lrg_data = load_all_lrg()
    dists = compute_full_distances(lrg_data)
    print(f"  Computed {len(dists)} pairwise distances")

    print("\nGenerating figures...")
    plot_grand_summary(csv_path)
    plot_per_patient_detail(dists)
    plot_patient_heatmaps(dists)
    plot_universality(dists)

    # Print summary table
    df = pd.read_csv(csv_path)
    mdf = df[df["metric"] == "1-LogCoph"]
    print("\n" + "=" * 70)
    print("SUMMARY TABLE: 1-LogCoph RI by patient and band")
    print("=" * 70)
    pivot = mdf.pivot(index="patient", columns="band", values="RI")
    pivot = pivot[BANDS]
    print(pivot.round(3).to_string())
    print("\nMeans:")
    print(pivot.mean().round(3).to_string())
    print("\nStd:")
    print(pivot.std().round(3).to_string())
    print("\nCV:")
    print((pivot.std() / pivot.mean()).round(3).to_string())


if __name__ == "__main__":
    main()
