#!/usr/bin/env python3
"""Final synthesis figure for the mean_ARI reorganization metric.

mean_ARI = mean Adjusted Rand Index across flat clusterings at k=2..20.
This metric achieves 6/6 bands unanimous for within > cross, making it
the most consistent measure of hierarchical reorganization.

Produces a single comprehensive multi-panel figure showing:
  A) Raw mean_ARI per pair and band (all patients overlaid)
  B) Within−Cross gap per patient × band (the core result)
  C) Band reorganization profile: which pairs are stable vs reorganized
  D) Per-patient phase×phase heatmaps (one per band, 4-phase patients)
  E) Pat_06 and Pat_07 control panels
  F) ARI across scales k: where does the signal come from?
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize

from lrg_eegfc.config import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.config.paths import FIGURES_ROOT

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
OUT_DIR = FIGURES_ROOT / "metric_exploration" / "mean_ari_final"

ALL_PAIRS = [
    ("rest_pre", "rest_post"), ("task_learn", "task_test"),
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"),
    ("task_learn", "rest_post"), ("task_test", "rest_post"),
]
PAIR_SHORT = {
    ("rest_pre", "rest_post"): "rest_pre↔rest_post",
    ("task_learn", "task_test"): "tLearn↔tTest",
    ("rest_pre", "task_learn"): "rest_pre↔tLearn",
    ("rest_pre", "task_test"): "rest_pre↔tTest",
    ("task_learn", "rest_post"): "tLearn↔rest_post",
    ("task_test", "rest_post"): "tTest↔rest_post",
}
PAIR_CAT = {
    ("rest_pre", "rest_post"): "within",
    ("task_learn", "task_test"): "within",
    ("rest_pre", "task_learn"): "cross",
    ("rest_pre", "task_test"): "cross",
    ("task_learn", "rest_post"): "cross",
    ("task_test", "rest_post"): "cross",
}
# Colors for pair types
WITHIN_COLOR = '#1565C0'
CROSS_COLOR = '#C62828'
TASK_TASK_COLOR = '#2E7D32'
REST_REST_COLOR = '#6A1B9A'

PAIR_COLORS = {
    ("rest_pre", "rest_post"): REST_REST_COLOR,
    ("task_learn", "task_test"): TASK_TASK_COLOR,
    ("rest_pre", "task_learn"): '#E65100',
    ("rest_pre", "task_test"): '#BF360C',
    ("task_learn", "rest_post"): '#AD1457',
    ("task_test", "rest_post"): '#880E4F',
}

def bl(b):
    return BRAIN_BAND_TEX_DICT.get(b, b)


def load_data():
    df = pd.read_csv(FIGURES_ROOT / "metric_exploration" / "partition_multiscale" / "results.csv")
    return df


def get_val(df, pat, band, p1, p2, col="mean_ARI"):
    row = df[(df["patient"] == pat) & (df["band"] == band) &
             (df["phase1"] == p1) & (df["phase2"] == p2)]
    if len(row) == 0:
        return np.nan
    return row[col].values[0]


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_data()
    print(f"Loaded {len(df)} rows")

    # ==========================================================
    # FIGURE 1: Core result — within>cross gap
    # ==========================================================
    fig1 = plt.figure(figsize=(20, 14))
    gs = gridspec.GridSpec(2, 3, figure=fig1, hspace=0.38, wspace=0.35,
                            height_ratios=[1.2, 1])

    # --- Panel A: Raw mean_ARI per pair across bands (mean ± patients) ---
    ax_a = fig1.add_subplot(gs[0, 0:2])

    x = np.arange(len(BANDS))
    for pair in ALL_PAIRS:
        p1, p2 = pair
        cat = PAIR_CAT[pair]
        color = PAIR_COLORS[pair]
        ls = '-' if cat == "within" else '--'
        lw = 2.5 if cat == "within" else 1.3
        ms = 8 if cat == "within" else 5
        marker = 's' if cat == "within" else 'o'

        # Per-patient thin lines
        for pat in PATIENTS_4PH:
            vals = [get_val(df, pat, b, p1, p2) for b in BANDS]
            ax_a.plot(x, vals, marker, color=color, alpha=0.2, markersize=3,
                      linestyle=ls, linewidth=0.5)

        # Mean thick line
        means = [np.nanmean([get_val(df, pat, b, p1, p2) for pat in PATIENTS_4PH])
                 for b in BANDS]
        stds = [np.nanstd([get_val(df, pat, b, p1, p2) for pat in PATIENTS_4PH])
                for b in BANDS]
        ax_a.errorbar(x, means, yerr=stds, fmt=f'{marker}{ls}', color=color,
                       linewidth=lw, markersize=ms, capsize=3,
                       label=f"{PAIR_SHORT[pair]} [{cat[0].upper()}]")

    ax_a.set_xticks(x)
    ax_a.set_xticklabels([bl(b) for b in BANDS], fontsize=11)
    ax_a.set_ylabel("mean ARI", fontsize=12)
    ax_a.set_title("A) Mean ARI per phase pair across frequency bands\n"
                    "(solid = within-type, dashed = cross-type, error bars = ±1 std)",
                    fontsize=11)
    ax_a.legend(fontsize=8, ncol=3, loc='lower right')
    ax_a.set_ylim(bottom=0.35)

    # --- Panel B: Within−Cross gap (THE core result) ---
    ax_b = fig1.add_subplot(gs[0, 2])

    gap_mat = np.full((len(BANDS), len(PATIENTS_4PH)), np.nan)
    for b_idx, band in enumerate(BANDS):
        for p_idx, pat in enumerate(PATIENTS_4PH):
            within_vals = [get_val(df, pat, band, p1, p2) for p1, p2 in ALL_PAIRS
                           if PAIR_CAT[(p1, p2)] == "within"]
            cross_vals = [get_val(df, pat, band, p1, p2) for p1, p2 in ALL_PAIRS
                          if PAIR_CAT[(p1, p2)] == "cross"]
            within_vals = [v for v in within_vals if not np.isnan(v)]
            cross_vals = [v for v in cross_vals if not np.isnan(v)]
            if within_vals and cross_vals:
                gap_mat[b_idx, p_idx] = np.mean(within_vals) - np.mean(cross_vals)

    vmax = max(0.05, np.nanmax(np.abs(gap_mat)) * 1.1)
    im_b = ax_b.imshow(gap_mat, cmap='Greens', vmin=0, vmax=vmax, aspect='auto')
    for i in range(len(BANDS)):
        for j in range(len(PATIENTS_4PH)):
            v = gap_mat[i, j]
            if not np.isnan(v):
                color = 'white' if v > vmax * 0.6 else 'black'
                ax_b.text(j, i, f"{v:+.3f}", ha='center', va='center',
                          fontsize=9, fontweight='bold', color=color)
    ax_b.set_xticks(range(len(PATIENTS_4PH)))
    ax_b.set_xticklabels([p.replace("Pat_0", "P") for p in PATIENTS_4PH], fontsize=10)
    ax_b.set_yticks(range(len(BANDS)))
    ax_b.set_yticklabels([bl(b) for b in BANDS], fontsize=11)
    plt.colorbar(im_b, ax=ax_b, shrink=0.7, label="Within − Cross gap")
    ax_b.set_title("B) Within−Cross gap\n(ALL positive = 6/6 bands ★)",
                    fontsize=11, fontweight='bold')

    # --- Panel C: Z-scored pair profile per band ---
    ax_c = fig1.add_subplot(gs[1, 0:2])

    # Z-score within each (patient, band), then average across patients
    z_means = np.full((len(BANDS), len(ALL_PAIRS)), np.nan)
    z_agreement = np.full((len(BANDS), len(ALL_PAIRS)), np.nan)
    for b_idx, band in enumerate(BANDS):
        for pr_idx, pair in enumerate(ALL_PAIRS):
            p1, p2 = pair
            zs = []
            for pat in PATIENTS_4PH:
                all_vals = [get_val(df, pat, band, pp1, pp2) for pp1, pp2 in ALL_PAIRS]
                all_vals = [v for v in all_vals if not np.isnan(v)]
                if not all_vals:
                    continue
                mu, sigma = np.mean(all_vals), np.std(all_vals)
                v = get_val(df, pat, band, p1, p2)
                if sigma > 1e-10 and not np.isnan(v):
                    zs.append((v - mu) / sigma)
            if zs:
                z_means[b_idx, pr_idx] = np.mean(zs)
                n_pos = sum(z > 0 for z in zs)
                z_agreement[b_idx, pr_idx] = max(n_pos, len(zs) - n_pos) / len(zs)

    vmax_z = max(1.5, np.nanmax(np.abs(z_means)) * 1.1)
    im_c = ax_c.imshow(z_means, cmap='RdBu', vmin=-vmax_z, vmax=vmax_z, aspect='auto')

    for i in range(len(BANDS)):
        for j in range(len(ALL_PAIRS)):
            z = z_means[i, j]
            ag = z_agreement[i, j]
            if np.isnan(z):
                continue
            unanimous = ag == 1.0
            relaxed = ag >= 0.75

            if unanimous:
                marker = "★"
                fw = 'bold'
            elif relaxed:
                marker = "◆"
                fw = 'normal'
            else:
                marker = ""
                fw = 'normal'

            txt = f"{marker}\n{z:+.2f}"
            color = 'white' if abs(z) > vmax_z * 0.45 else 'black'
            ax_c.text(j, i, txt, ha='center', va='center', fontsize=8,
                      color=color, fontweight=fw)

            if unanimous:
                rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                      linewidth=2.5, edgecolor='gold', facecolor='none')
                ax_c.add_patch(rect)

    ax_c.set_xticks(range(len(ALL_PAIRS)))
    xlabels = [f"{PAIR_SHORT[p]}\n[{'W' if PAIR_CAT[p]=='within' else 'C'}]" for p in ALL_PAIRS]
    ax_c.set_xticklabels(xlabels, fontsize=8)
    ax_c.set_yticks(range(len(BANDS)))
    ax_c.set_yticklabels([bl(b) for b in BANDS], fontsize=11)
    plt.colorbar(im_c, ax=ax_c, shrink=0.7, label="z-score (+ = above avg = stable)")
    ax_c.set_title("C) Z-scored mean ARI: which pairs are stable vs reorganized\n"
                    "★ gold = unanimous across 4 patients, ◆ = 3/4 agree",
                    fontsize=10)

    # --- Panel D: Reorganization strength per band ---
    ax_d = fig1.add_subplot(gs[1, 2])

    strengths = {}
    for pat in PATIENTS_4PH:
        for band in BANDS:
            vals = [get_val(df, pat, band, p1, p2) for p1, p2 in ALL_PAIRS]
            vals = [v for v in vals if not np.isnan(v)]
            if len(vals) >= 2:
                strengths[(pat, band)] = np.std(vals)

    band_x = np.arange(len(BANDS))
    width = 0.15
    for i, pat in enumerate(PATIENTS_4PH):
        vals = [strengths.get((pat, b), 0) for b in BANDS]
        offset = (i - len(PATIENTS_4PH) / 2 + 0.5) * width
        ax_d.bar(band_x + offset, vals, width * 0.9, alpha=0.7,
                 label=pat.replace("Pat_0", "P"))

    means = [np.mean([strengths.get((pat, b), 0) for pat in PATIENTS_4PH]) for b in BANDS]
    ax_d.plot(band_x, means, 'k-D', linewidth=2, markersize=6, label='Mean', zorder=10)
    ax_d.set_xticks(band_x)
    ax_d.set_xticklabels([bl(b) for b in BANDS], fontsize=9, rotation=30)
    ax_d.set_ylabel("std(pairwise ARI)")
    ax_d.set_title("D) Reorganization strength\n(std across pairs per band)", fontsize=10)
    ax_d.legend(fontsize=7, ncol=3)

    fig1.suptitle("Hierarchical reorganization analysis — mean ARI metric\n"
                  "LRG dendrograms on MSC matrices, flat clustering at k=2..20",
                  fontsize=14, fontweight='bold', y=0.99)
    fig1.savefig(OUT_DIR / "mean_ari_core_result.pdf", bbox_inches='tight', dpi=150)
    plt.close(fig1)
    print("Saved mean_ari_core_result.pdf")

    # ==========================================================
    # FIGURE 2: Per-band phase×phase heatmaps (4-phase patients)
    # ==========================================================
    fig2, axes2 = plt.subplots(2, 3, figsize=(18, 11))
    fig2.suptitle("Per-band phase×phase mean ARI (averaged over 4 patients)\n"
                  "Rows/cols = experimental phases, color = mean ARI similarity",
                  fontsize=13, fontweight='bold')

    phases_4 = ["rest_pre", "task_learn", "task_test", "rest_post"]
    phase_labels = ["rest_pre", "tLearn", "tTest", "rest_post"]

    for b_idx, band in enumerate(BANDS):
        ax = axes2[b_idx // 3, b_idx % 3]
        mat = np.eye(4)  # diagonal = 1
        for i, pi in enumerate(phases_4):
            for j, pj in enumerate(phases_4):
                if i >= j:
                    continue
                vals = [get_val(df, pat, band, pi, pj) for pat in PATIENTS_4PH]
                vals = [v for v in vals if not np.isnan(v)]
                if vals:
                    m = np.mean(vals)
                    mat[i, j] = m
                    mat[j, i] = m

        im = ax.imshow(mat, cmap='YlOrRd', vmin=0.4, vmax=1.0, aspect='equal')
        for i in range(4):
            for j in range(4):
                color = 'white' if mat[i, j] > 0.75 else 'black'
                ax.text(j, i, f"{mat[i,j]:.3f}", ha='center', va='center',
                        fontsize=9, fontweight='bold', color=color)
        ax.set_xticks(range(4))
        ax.set_xticklabels(phase_labels, fontsize=9, rotation=30)
        ax.set_yticks(range(4))
        ax.set_yticklabels(phase_labels, fontsize=9)
        ax.set_title(f"{bl(band)}", fontsize=12, fontweight='bold')

    fig2.tight_layout(rect=[0, 0, 0.92, 0.93])
    cbar_ax = fig2.add_axes([0.93, 0.15, 0.02, 0.7])
    fig2.colorbar(im, cax=cbar_ax, label="mean ARI")
    fig2.savefig(OUT_DIR / "phase_phase_heatmaps.pdf", bbox_inches='tight', dpi=150)
    plt.close(fig2)
    print("Saved phase_phase_heatmaps.pdf")

    # ==========================================================
    # FIGURE 3: All 6 patients — per-patient profiles
    # ==========================================================
    fig3, axes3 = plt.subplots(2, 3, figsize=(18, 10))
    fig3.suptitle("Mean ARI per phase pair — all 6 patients\n"
                  "(solid = within-type, dashed = cross-type)",
                  fontsize=13, fontweight='bold')

    for p_idx, (pat, phases) in enumerate(PATIENT_PHASES.items()):
        ax = axes3[p_idx // 3, p_idx % 3]
        pairs = list(combinations(phases, 2))

        for pair in pairs:
            p1, p2 = pair
            cat = PAIR_CAT.get(pair, "?")
            color = PAIR_COLORS.get(pair, 'gray')
            ls = '-' if cat == "within" else '--'
            lw = 2.5 if cat == "within" else 1.2
            ms = 7 if cat == "within" else 4
            marker = 's' if cat == "within" else 'o'

            vals = [get_val(df, pat, b, p1, p2) for b in BANDS]
            ax.plot(range(len(BANDS)), vals, f'{marker}{ls}', color=color,
                    linewidth=lw, markersize=ms,
                    label=PAIR_SHORT.get(pair, f"{p1}↔{p2}"), alpha=0.8)

        ax.set_xticks(range(len(BANDS)))
        ax.set_xticklabels([bl(b) for b in BANDS], fontsize=8, rotation=30)
        ax.set_ylabel("mean ARI")
        ax.set_title(f"{pat} ({len(phases)} phases)", fontsize=11)
        ax.legend(fontsize=6, ncol=2, loc='lower right')
        ax.set_ylim(0.3, 1.02)

    fig3.tight_layout(rect=[0, 0, 1, 0.93])
    fig3.savefig(OUT_DIR / "all_patients_profiles.pdf", bbox_inches='tight', dpi=150)
    plt.close(fig3)
    print("Saved all_patients_profiles.pdf")

    # ==========================================================
    # FIGURE 4: ARI across scales — where does the signal come from?
    # ==========================================================
    k_values = [2, 3, 4, 5, 6, 8, 10, 15, 20]
    k_cols = [f"ARI_k{k}" for k in k_values]

    fig4, axes4 = plt.subplots(2, 3, figsize=(18, 10))
    fig4.suptitle("ARI at each scale k: where does the within>cross signal come from?\n"
                  "(green = within, red = cross, averaged over 4 patients)",
                  fontsize=13, fontweight='bold')

    for b_idx, band in enumerate(BANDS):
        ax = axes4[b_idx // 3, b_idx % 3]

        within_means, cross_means = [], []
        within_stds, cross_stds = [], []
        for kcol in k_cols:
            wv = []
            cv = []
            for pat in PATIENTS_4PH:
                for pair in ALL_PAIRS:
                    p1, p2 = pair
                    v = get_val(df, pat, band, p1, p2, col=kcol)
                    if not np.isnan(v):
                        if PAIR_CAT[pair] == "within":
                            wv.append(v)
                        else:
                            cv.append(v)
            within_means.append(np.mean(wv) if wv else np.nan)
            within_stds.append(np.std(wv) if wv else np.nan)
            cross_means.append(np.mean(cv) if cv else np.nan)
            cross_stds.append(np.std(cv) if cv else np.nan)

        ax.fill_between(range(len(k_values)),
                         np.array(within_means) - np.array(within_stds),
                         np.array(within_means) + np.array(within_stds),
                         alpha=0.15, color='green')
        ax.fill_between(range(len(k_values)),
                         np.array(cross_means) - np.array(cross_stds),
                         np.array(cross_means) + np.array(cross_stds),
                         alpha=0.15, color='red')
        ax.plot(range(len(k_values)), within_means, 's-', color='green',
                linewidth=2, markersize=6, label='Within')
        ax.plot(range(len(k_values)), cross_means, 'o--', color='red',
                linewidth=2, markersize=5, label='Cross')

        # Gap
        gaps = [w - c for w, c in zip(within_means, cross_means)]
        ax2 = ax.twinx()
        ax2.bar(range(len(k_values)), gaps, alpha=0.2, color='blue', width=0.6)
        ax2.set_ylabel("Gap", fontsize=8, color='blue')
        ax2.tick_params(axis='y', labelcolor='blue', labelsize=7)

        ax.set_xticks(range(len(k_values)))
        ax.set_xticklabels([str(k) for k in k_values], fontsize=8)
        ax.set_xlabel("k (number of clusters)")
        ax.set_ylabel("ARI")
        ax.set_title(bl(band), fontsize=12, fontweight='bold')
        ax.legend(fontsize=7, loc='upper right')

    fig4.tight_layout(rect=[0, 0, 1, 0.93])
    fig4.savefig(OUT_DIR / "ari_across_scales.pdf", bbox_inches='tight', dpi=150)
    plt.close(fig4)
    print("Saved ari_across_scales.pdf")

    # ==========================================================
    # FIGURE 5: Specific contrasts — task stability, rest stability
    # ==========================================================
    fig5, axes5 = plt.subplots(1, 4, figsize=(20, 5))

    contrast_defs = [
        ("Within − Cross", lambda pv: np.mean([pv[p] for p in pv if PAIR_CAT[p] == "within"]) -
                                       np.mean([pv[p] for p in pv if PAIR_CAT[p] == "cross"])),
        ("Task↔Task − Cross\n(task stability)", lambda pv: pv.get(("task_learn", "task_test"), np.nan) -
                                       np.mean([pv[p] for p in pv if PAIR_CAT[p] == "cross"])),
        ("Rest↔Rest − Cross\n(rest stability)", lambda pv: pv.get(("rest_pre", "rest_post"), np.nan) -
                                       np.mean([pv[p] for p in pv if PAIR_CAT[p] == "cross"])),
        ("Task↔Task − Rest↔Rest", lambda pv: pv.get(("task_learn", "task_test"), np.nan) -
                                               pv.get(("rest_pre", "rest_post"), np.nan)),
    ]

    pat_markers = {'Pat_02': 'o', 'Pat_03': 's', 'Pat_05': 'D', 'Pat_08': '^'}

    for c_idx, (cname, cfunc) in enumerate(contrast_defs):
        ax = axes5[c_idx]
        for b_idx, band in enumerate(BANDS):
            vals = []
            for pat in PATIENTS_4PH:
                pv = {}
                for pair in ALL_PAIRS:
                    p1, p2 = pair
                    v = get_val(df, pat, band, p1, p2)
                    if not np.isnan(v):
                        pv[pair] = v
                if len(pv) == 6:
                    try:
                        c = cfunc(pv)
                        vals.append(c)
                        ax.scatter(c, b_idx, marker=pat_markers[pat], s=40,
                                   alpha=0.6, color='#1565C0', edgecolors='black',
                                   linewidths=0.3)
                    except Exception:
                        pass

            if vals:
                m = np.mean(vals)
                unanimous = all(v > 0 for v in vals) or all(v < 0 for v in vals)
                ax.scatter(m, b_idx, marker='D', s=100, color='black', zorder=10)
                if unanimous:
                    ax.annotate("★", (m, b_idx), textcoords="offset points",
                                xytext=(12, 0), fontsize=14, color='gold',
                                fontweight='bold', va='center')

        ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)
        ax.set_yticks(range(len(BANDS)))
        ax.set_yticklabels([bl(b) for b in BANDS] if c_idx == 0 else [], fontsize=10)
        ax.invert_yaxis()
        ax.set_xlabel("Contrast (mean ARI)")
        ax.set_title(cname, fontsize=10, fontweight='bold')

    fig5.suptitle("Reorganization contrasts — mean ARI\n"
                  "Dots = individual patients, Diamond = mean, ★ = unanimous sign",
                  fontsize=12, fontweight='bold')
    fig5.tight_layout(rect=[0, 0, 1, 0.88])
    fig5.savefig(OUT_DIR / "contrasts.pdf", bbox_inches='tight', dpi=150)
    plt.close(fig5)
    print("Saved contrasts.pdf")

    # ==========================================================
    # Print final summary
    # ==========================================================
    print("\n" + "=" * 70)
    print("MEAN ARI — FINAL RESULTS SUMMARY")
    print("=" * 70)

    print("\n--- Within−Cross gap (unanimous ALL 6 bands) ---")
    for band in BANDS:
        gaps = []
        for pat in PATIENTS_4PH:
            wv = [get_val(df, pat, band, p1, p2) for p1, p2 in ALL_PAIRS
                  if PAIR_CAT[(p1, p2)] == "within"]
            cv = [get_val(df, pat, band, p1, p2) for p1, p2 in ALL_PAIRS
                  if PAIR_CAT[(p1, p2)] == "cross"]
            wv = [v for v in wv if not np.isnan(v)]
            cv = [v for v in cv if not np.isnan(v)]
            if wv and cv:
                gaps.append(np.mean(wv) - np.mean(cv))
        all_pos = all(g > 0 for g in gaps)
        print(f"  {bl(band):<14} gaps: {' '.join(f'{g:+.4f}' for g in gaps)}  "
              f"{'★ ALL+' if all_pos else 'MIXED'}")

    print("\n--- Mean ARI per pair (averaged over 4 patients) ---")
    print(f"  {'Pair':<20} {'Type':<8}", end="")
    for band in BANDS:
        print(f"  {bl(band):>8}", end="")
    print()
    for pair in ALL_PAIRS:
        p1, p2 = pair
        cat = PAIR_CAT[pair]
        print(f"  {PAIR_SHORT[pair]:<20} {cat:<8}", end="")
        for band in BANDS:
            vals = [get_val(df, pat, band, p1, p2) for pat in PATIENTS_4PH]
            vals = [v for v in vals if not np.isnan(v)]
            m = np.mean(vals) if vals else np.nan
            print(f"  {m:>8.4f}", end="")
        print()

    print("\n--- Task stability: tLearn↔tTest is rank 1? ---")
    for band in BANDS:
        ranks = []
        for pat in PATIENTS_4PH:
            pair_vals = []
            for pair in ALL_PAIRS:
                v = get_val(df, pat, band, pair[0], pair[1])
                pair_vals.append((pair, v))
            pair_vals.sort(key=lambda x: -x[1])
            rank = next(i + 1 for i, (p, _) in enumerate(pair_vals)
                        if p == ("task_learn", "task_test"))
            ranks.append(rank)
        all_rank1 = all(r == 1 for r in ranks)
        print(f"  {bl(band):<14} ranks: {ranks}  "
              f"{'★ ALL RANK 1' if all_rank1 else f'best={min(ranks)}, worst={max(ranks)}'}")

    print("\n--- Reorganization strength (std of pairwise ARI) ---")
    for band in BANDS:
        stds = []
        for pat in PATIENTS_4PH:
            vals = [get_val(df, pat, band, p1, p2) for p1, p2 in ALL_PAIRS]
            vals = [v for v in vals if not np.isnan(v)]
            if len(vals) >= 2:
                stds.append(np.std(vals))
        print(f"  {bl(band):<14} mean_std={np.mean(stds):.4f}  "
              f"per_pat: {' '.join(f'{s:.4f}' for s in stds)}")

    print(f"\nAll figures saved to: {OUT_DIR}")
    print("Done!")


if __name__ == "__main__":
    main()
