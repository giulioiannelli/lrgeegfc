#!/usr/bin/env python3
"""Final synthesis figure: Band × Phase reorganization of LRG ultrametric D(τ).

Creates a single comprehensive figure summarizing:
  - Which bands reorganize (frequency gradient)
  - Between which phases (task stability, rest variability)
  - Which results are consistent across ALL patients

Uses logCosine of D(τ) on MSC-derived LRG dendrograms.
4-phase patients: Pat_02, Pat_03, Pat_05, Pat_08.
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize, LinearSegmentedColormap
from scipy.spatial.distance import squareform

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
PATIENTS_4PH = PATIENTS_4PHASE
BANDS = BRAIN_BANDS_NAMES
OUT_DIR = FIGURES_ROOT / "metric_exploration"

ALL_PAIRS = [
    ("rest_pre", "rest_post"),
    ("task_learn", "task_test"),
    ("rest_pre", "task_learn"),
    ("rest_pre", "task_test"),
    ("task_learn", "rest_post"),
    ("task_test", "rest_post"),
]

PAIR_SHORT = {
    ("rest_pre", "rest_post"): "rest_pre↔rest_post",
    ("task_learn", "task_test"): "tLearn↔tTest",
    ("rest_pre", "task_learn"): "rest_pre↔tLearn",
    ("rest_pre", "task_test"): "rest_pre↔tTest",
    ("task_learn", "rest_post"): "tLearn↔rest_post",
    ("task_test", "rest_post"): "tTest↔rest_post",
}

PAIR_CATEGORY = {
    ("rest_pre", "rest_post"): "within",
    ("task_learn", "task_test"): "within",
    ("rest_pre", "task_learn"): "cross",
    ("rest_pre", "task_test"): "cross",
    ("task_learn", "rest_post"): "cross",
    ("task_test", "rest_post"): "cross",
}


def bl(b):
    return BRAIN_BAND_TEX_DICT.get(b, b)


def log_cosine(U1, U2):
    lU1 = np.log(np.clip(U1, 1e-15, None))
    lU2 = np.log(np.clip(U2, 1e-15, None))
    return np.dot(lU1, lU2) / max(np.linalg.norm(lU1) * np.linalg.norm(lU2), 1e-12)


def load_all():
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


def compute_sims(data):
    sims = {}
    for pat, phases in PATIENT_PHASES.items():
        for band in BANDS:
            for p1, p2 in combinations(phases, 2):
                k1, k2 = (pat, p1, band), (pat, p2, band)
                if k1 in data and k2 in data:
                    sims[(pat, band, p1, p2)] = log_cosine(data[k1], data[k2])
    return sims


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    data = load_all()
    sims = compute_sims(data)
    print(f"  {len(data)} ultrametrics, {len(sims)} pairwise similarities")

    patients = PATIENTS_4PH

    # ================================================================
    # Compute contrasts per (patient, band)
    # ================================================================
    within_pairs = [p for p in ALL_PAIRS if PAIR_CATEGORY[p] == "within"]
    cross_pairs = [p for p in ALL_PAIRS if PAIR_CATEGORY[p] == "cross"]

    # Store raw values and contrasts
    raw = {}  # (band, pair) → list of patient values
    contrasts = {}  # (band, contrast_name) → list of patient values
    strengths = {}  # (pat, band) → std

    for pat in patients:
        for band in BANDS:
            pair_vals = {}
            for pair in ALL_PAIRS:
                p1, p2 = pair
                key = (pat, band, p1, p2)
                if key in sims:
                    pair_vals[pair] = sims[key]
                    raw.setdefault((band, pair), []).append(sims[key])

            if len(pair_vals) == 6:
                within_v = np.mean([pair_vals[p] for p in within_pairs])
                cross_v = np.mean([pair_vals[p] for p in cross_pairs])
                tt = pair_vals[("task_learn", "task_test")]
                rr = pair_vals[("rest_pre", "rest_post")]

                contrasts.setdefault((band, "within_minus_cross"), []).append(within_v - cross_v)
                contrasts.setdefault((band, "task_stability"), []).append(tt - cross_v)
                contrasts.setdefault((band, "rest_stability"), []).append(rr - cross_v)
                contrasts.setdefault((band, "task_minus_rest"), []).append(tt - rr)

                strengths[(pat, band)] = np.std(list(pair_vals.values()))

    # ================================================================
    # SYNTHESIS FIGURE
    # ================================================================
    fig = plt.figure(figsize=(22, 14))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.35,
                           width_ratios=[2.5, 1.5, 1.5])

    # ---- Panel A: Raw similarity per pair (band × pair heatmap) ----
    ax_a = fig.add_subplot(gs[0, 0])

    mat_mean = np.full((len(BANDS), len(ALL_PAIRS)), np.nan)
    mat_std = np.full((len(BANDS), len(ALL_PAIRS)), np.nan)
    for i, band in enumerate(BANDS):
        for j, pair in enumerate(ALL_PAIRS):
            vals = raw.get((band, pair))
            if vals:
                mat_mean[i, j] = np.mean(vals)
                mat_std[i, j] = np.std(vals)

    im_a = ax_a.imshow(mat_mean, cmap='YlOrRd', aspect='auto',
                        vmin=np.nanmin(mat_mean) - 0.01,
                        vmax=np.nanmax(mat_mean) + 0.01)

    for i in range(len(BANDS)):
        for j in range(len(ALL_PAIRS)):
            m = mat_mean[i, j]
            s = mat_std[i, j]
            if np.isnan(m):
                continue
            # Check unanimity of z-score direction
            vals = raw.get((BANDS[i], ALL_PAIRS[j]), [])
            if vals:
                arr = np.array(vals)
                band_mean = np.mean([np.mean(raw.get((BANDS[i], p), [np.nan]))
                                     for p in ALL_PAIRS if raw.get((BANDS[i], p))])
                above = sum(v > band_mean for v in vals)
                below = len(vals) - above
                if above == len(vals) or below == len(vals):
                    marker = "★"
                elif max(above, below) >= len(vals) - 1:
                    marker = "◆"
                else:
                    marker = ""
            else:
                marker = ""
            color = 'white' if m > np.nanmean(mat_mean) + 0.01 else 'black'
            ax_a.text(j, i, f"{m:.4f}\n±{s:.4f}\n{marker}", ha='center', va='center',
                      fontsize=6.5, color=color, fontweight='bold' if marker == "★" else 'normal')

    ax_a.set_xticks(range(len(ALL_PAIRS)))
    xlabels = []
    for p in ALL_PAIRS:
        cat = PAIR_CATEGORY[p]
        color_tag = "W" if cat == "within" else "C"
        xlabels.append(f"{PAIR_SHORT[p]}\n[{color_tag}]")
    ax_a.set_xticklabels(xlabels, fontsize=7)
    ax_a.set_yticks(range(len(BANDS)))
    ax_a.set_yticklabels([bl(b) for b in BANDS], fontsize=10)
    plt.colorbar(im_a, ax=ax_a, shrink=0.7, label="logCosine D(τ)")
    ax_a.set_title("A) Mean pairwise similarity (4 patients)\n"
                    "★=unanimous direction, ◆=N−1 agree", fontsize=10)

    # ---- Panel B: Contrasts as 4 horizontal strip charts ----
    gs_b = gs[0, 1:].subgridspec(4, 1, hspace=0.5)
    contrast_names = ["within_minus_cross", "task_stability", "rest_stability", "task_minus_rest"]
    contrast_labels = [
        "Within − Cross (category effect)",
        "Task↔Task − Cross (task stability)",
        "Rest↔Rest − Cross (rest stability)",
        "Task↔Task − Rest↔Rest",
    ]
    colors_contrast = ['#2196F3', '#4CAF50', '#FF9800', '#9C27B0']

    for c_idx, (cname, clabel) in enumerate(zip(contrast_names, contrast_labels)):
        ax_b = fig.add_subplot(gs_b[c_idx])
        band_ys = np.arange(len(BANDS))

        for b_idx, band in enumerate(BANDS):
            vals = contrasts.get((band, cname), [])
            if not vals:
                continue
            arr = np.array(vals)
            m = arr.mean()
            n_pos = sum(arr > 0)
            N = len(arr)
            unanimous = n_pos == N or n_pos == 0
            relaxed = max(n_pos, N - n_pos) >= N - 1

            # Individual patients
            for v in vals:
                ax_b.scatter(v * 1000, b_idx, s=20, alpha=0.5,
                             color=colors_contrast[c_idx], edgecolors='none')
            # Mean diamond
            ax_b.scatter(m * 1000, b_idx, s=70, marker='D',
                         color=colors_contrast[c_idx], edgecolors='black',
                         linewidths=0.5, zorder=10)
            if unanimous:
                ax_b.annotate("★", (m * 1000, b_idx), textcoords="offset points",
                              xytext=(10, 0), fontsize=11, color='gold',
                              fontweight='bold', va='center')
            elif relaxed:
                ax_b.annotate("◆", (m * 1000, b_idx), textcoords="offset points",
                              xytext=(10, 0), fontsize=9, color='orange', va='center')

        ax_b.axvline(0, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
        ax_b.set_yticks(band_ys)
        ax_b.set_yticklabels([bl(b) for b in BANDS], fontsize=8)
        ax_b.set_title(f"B{c_idx+1}) {clabel}", fontsize=8, fontweight='bold',
                        pad=2)
        ax_b.invert_yaxis()
        if c_idx == 3:
            ax_b.set_xlabel("Δ × 10³ (logCosine D(τ))", fontsize=8)
        else:
            ax_b.set_xticklabels([])

    # ---- Panel C: Reorganization strength (frequency gradient) ----
    ax_c = fig.add_subplot(gs[1, 0])

    mean_str = [np.mean([strengths.get((pat, band), 0) for pat in patients])
                for band in BANDS]
    pat_str = {pat: [strengths.get((pat, band), 0) for band in BANDS]
               for pat in patients}

    x = np.arange(len(BANDS))
    width = 0.15
    for i, pat in enumerate(patients):
        offset = (i - len(patients) / 2 + 0.5) * width
        ax_c.bar(x + offset, pat_str[pat], width * 0.9,
                 label=pat.replace("Pat_0", "P"), alpha=0.7)
    ax_c.plot(x, mean_str, 'k-o', linewidth=2.5, markersize=7,
              label='Mean', zorder=10)

    ax_c.set_xticks(x)
    ax_c.set_xticklabels([bl(b) for b in BANDS], fontsize=10)
    ax_c.set_ylabel("Reorganization strength\n(std of pairwise logCosine)")
    ax_c.set_title("C) Frequency gradient of reorganization\n"
                    "(higher = more phase-dependent hierarchy)", fontsize=10)
    ax_c.legend(fontsize=7, ncol=5)

    # ---- Panel D: Pat_06 control ----
    ax_d = fig.add_subplot(gs[1, 1])

    pat06_vals = [sims.get(("Pat_06", band, "rest_pre", "rest_post"), np.nan) for band in BANDS]
    mean_4ph_rr = [np.mean(raw.get((band, ("rest_pre", "rest_post")), [np.nan])) for band in BANDS]

    ax_d.barh(range(len(BANDS)), pat06_vals, height=0.35, alpha=0.7,
              color='#E91E63', label='Pat_06 (rest only)')
    ax_d.barh([y + 0.35 for y in range(len(BANDS))], mean_4ph_rr, height=0.35,
              alpha=0.7, color='#2196F3', label='4-phase mean')
    ax_d.set_yticks([y + 0.175 for y in range(len(BANDS))])
    ax_d.set_yticklabels([bl(b) for b in BANDS], fontsize=10)
    ax_d.set_xlabel("logCosine D(τ) [rest_pre↔rest_post]")
    ax_d.set_title("D) Rest stability control\n"
                    "Pat_06 (epileptic, no task)", fontsize=10)
    ax_d.legend(fontsize=7, loc='lower right')
    ax_d.invert_yaxis()

    # ---- Panel E: Summary table ----
    ax_e = fig.add_subplot(gs[1, 2])
    ax_e.axis('off')

    summary_text = (
        "CONSISTENT FINDINGS\n"
        "(unanimous across 4 patients)\n"
        "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        "\n"
        "▶ Within > Cross (category effect):\n"
        "   θ, β, γ_l, γ_h\n"
        "\n"
        "▶ Task hierarchy stable\n"
        "  (tLearn↔tTest > cross avg):\n"
        "   α, β, γ_l, γ_h\n"
        "\n"
        "▶ Task > Rest stability:\n"
        "   α, β\n"
        "\n"
        "▶ Rest NOT consistently stable\n"
        "  in ANY band\n"
        "\n"
        "▶ Frequency gradient:\n"
        "   δ > θ > α > β > γ_l > γ_h\n"
        "   (low freq = more reorganization)\n"
        "\n"
        "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
        "Metric: logCosine of D(τ)\n"
        "FC: MSC (nperseg=4096)\n"
        "Hierarchy: LRG ultrametric"
    )
    ax_e.text(0.05, 0.95, summary_text, transform=ax_e.transAxes,
              fontsize=9, fontfamily='monospace', verticalalignment='top',
              bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow',
                        edgecolor='gray', alpha=0.9))

    fig.suptitle("Ultrametric reorganization: Band × Phase analysis\n"
                 "logCosine D(τ) — LRG on MSC matrices — 4-phase patients",
                 fontsize=14, fontweight='bold', y=0.98)

    fig.savefig(OUT_DIR / "band_phase_synthesis.pdf", bbox_inches='tight', dpi=150)
    plt.close(fig)
    print("Saved band_phase_synthesis.pdf")

    # ================================================================
    # Also include Pat_06 and Pat_07 in a separate 6-patient overview
    # ================================================================
    fig2, axes2 = plt.subplots(2, 3, figsize=(18, 10))
    fig2.suptitle("Per-patient logCosine D(τ) similarity profiles\n"
                  "All 6 patients — logCosine of LRG ultrametric on MSC",
                  fontsize=13, fontweight='bold')

    all_patients = list(PATIENT_PHASES.keys())
    for p_idx, pat in enumerate(all_patients):
        ax = axes2[p_idx // 3, p_idx % 3]
        phases = PATIENT_PHASES[pat]
        pairs = list(combinations(phases, 2))

        for pair in pairs:
            p1, p2 = pair
            vals = [sims.get((pat, band, p1, p2), np.nan) for band in BANDS]
            cat = PAIR_CATEGORY.get(pair, "?")
            ls = '-' if cat == "within" else '--'
            lw = 2 if cat == "within" else 1
            ax.plot(range(len(BANDS)), vals, f'o{ls}', markersize=4,
                    linewidth=lw, label=PAIR_SHORT.get(pair, f"{p1}↔{p2}"),
                    alpha=0.8)

        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels([bl(b) for b in BANDS], fontsize=8, rotation=30)
        ax.set_ylabel("logCosine D(τ)")
        n_phases = len(phases)
        n_pairs = len(pairs)
        ax.set_title(f"{pat} ({n_phases} phases, {n_pairs} pairs)", fontsize=10)
        ax.legend(fontsize=6, ncol=2, loc='lower right')
        ax.set_ylim(bottom=min(0.55, ax.get_ylim()[0]))

    fig2.tight_layout(rect=[0, 0, 1, 0.93])
    fig2.savefig(OUT_DIR / "all_patients_profiles.pdf", bbox_inches='tight', dpi=150)
    plt.close(fig2)
    print("Saved all_patients_profiles.pdf")

    print("\nDone!")


if __name__ == "__main__":
    main()
