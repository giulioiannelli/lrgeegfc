#!/usr/bin/env python3
"""Clear reorganization visualization: within-group vs cross-group ARI curves.

For each (patient, band), computes:
  - Mean ARI(h) across within-condition pairs (rest_pre↔rest_post, task_learn↔task_test)
  - Mean ARI(h) across cross-condition pairs (4 rest↔task pairs)

If reorganization exists, within-condition ARI > cross-condition ARI, i.e.
the hierarchical structure is more preserved within rest/task than across them.

Main figure: 6 panels (one per band), each with mean±std across patients.
Summary: bar chart of the reorganization gap per band.

All patients included according to their available phases.
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from scipy.cluster.hierarchy import fcluster
from sklearn.metrics import adjusted_rand_score

from lrg_eegfc.config import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import LRG_CACHE as _LRG_CACHE, FIGURES_ROOT
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
PATIENTS = list(PATIENT_PHASES.keys())
BANDS = BRAIN_BANDS_NAMES
BAND_TEX = [BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS]

USE_CREMA = "--crema" in sys.argv

from lrg_eegfc.config.paths import LRG_CREMA_CACHE
LRG_CACHE = LRG_CREMA_CACHE if USE_CREMA else _LRG_CACHE
OUTPUT_DIR = FIGURES_ROOT / "crema_metric_exploration" if USE_CREMA else FIGURES_ROOT / "metric_exploration"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

REST = {"rest_pre", "rest_post"}
TASK = {"task_learn", "task_test"}

N_THRESHOLDS = 200
THRESHOLDS = np.linspace(0.01, 1.0, N_THRESHOLDS)

PAT_COLORS = {
    "Pat_02": "#66c2a5", "Pat_03": "#fc8d62", "Pat_05": "#8da0cb",
    "Pat_06": "#a6d854", "Pat_07": "#ffd92f", "Pat_08": "#e78ac3",
}


# ===================================================================
# Compute ARI curve between two linkage matrices
# ===================================================================
def ari_curve(Z1, Z2):
    """Return ARI(h) for h in THRESHOLDS."""
    mh1, mh2 = Z1[:, 2].max(), Z2[:, 2].max()
    ari = np.empty(N_THRESHOLDS)
    for i, t in enumerate(THRESHOLDS):
        l1 = fcluster(Z1, t=t * mh1, criterion="distance")
        l2 = fcluster(Z2, t=t * mh2, criterion="distance")
        ari[i] = adjusted_rand_score(l1, l2)
    return ari


def classify_pair(p1, p2):
    s = {p1, p2}
    if s <= REST:
        return "within"
    elif s <= TASK:
        return "within"
    else:
        return "cross"


# ===================================================================
# Compute everything
# ===================================================================
def compute_all():
    """For each (patient, band), compute mean within/cross ARI curves.

    Returns dict: (patient, band) → {"within": ndarray, "cross": ndarray,
                                       "within_all": list, "cross_all": list}
    Only includes patients that have at least 1 within pair AND 1 cross pair.
    """
    results = {}

    for pat, phases in PATIENT_PHASES.items():
        # Load LRG
        lrgs = {}
        for band in BANDS:
            for ph in phases:
                lrg = load_lrg_result(pat, ph, band, "msc", cache_root=LRG_CACHE)
                if lrg is not None:
                    lrgs[(band, ph)] = lrg

        for band in BANDS:
            available = [ph for ph in phases if (band, ph) in lrgs]
            pairs = list(combinations(available, 2))

            within_curves = []
            cross_curves = []
            for p1, p2 in pairs:
                Z1 = lrgs[(band, p1)].linkage_matrix
                Z2 = lrgs[(band, p2)].linkage_matrix
                ari = ari_curve(Z1, Z2)
                cat = classify_pair(p1, p2)
                if cat == "within":
                    within_curves.append(ari)
                else:
                    cross_curves.append(ari)

            if within_curves and cross_curves:
                results[(pat, band)] = {
                    "within": np.mean(within_curves, axis=0),
                    "cross": np.mean(cross_curves, axis=0),
                    "within_all": within_curves,
                    "cross_all": cross_curves,
                    "n_within": len(within_curves),
                    "n_cross": len(cross_curves),
                }

    return results


# ===================================================================
# Figure 1: ARI curves — within vs cross (one panel per band)
# ===================================================================
def plot_ari_curves(results):
    """6 panels: mean±std ARI curves for within (blue) and cross (red)."""
    fig, axes = plt.subplots(2, 3, figsize=(20, 12), sharex=True, sharey=True)
    axes_flat = axes.ravel()

    for bi, band in enumerate(BANDS):
        ax = axes_flat[bi]

        # Collect per-patient within and cross curves
        within_all = []
        cross_all = []
        pats_used = []
        for pat in PATIENTS:
            key = (pat, band)
            if key in results:
                within_all.append(results[key]["within"])
                cross_all.append(results[key]["cross"])
                pats_used.append(pat)

        if not within_all:
            ax.set_title(BAND_TEX[bi], fontsize=13, fontweight="bold")
            continue

        W = np.array(within_all)
        C = np.array(cross_all)
        w_mean, w_std = W.mean(axis=0), W.std(axis=0)
        c_mean, c_std = C.mean(axis=0), C.std(axis=0)

        # Plot mean curves
        ax.plot(THRESHOLDS, w_mean, color="#2166ac", lw=2.5, label="Within-condition")
        ax.fill_between(THRESHOLDS, w_mean - w_std, w_mean + w_std,
                        color="#2166ac", alpha=0.15)
        ax.plot(THRESHOLDS, c_mean, color="#b2182b", lw=2.5, label="Cross-condition")
        ax.fill_between(THRESHOLDS, c_mean - c_std, c_mean + c_std,
                        color="#b2182b", alpha=0.15)

        # Shade the gap where within > cross (reorganization signal)
        ax.fill_between(THRESHOLDS, c_mean, w_mean,
                        where=w_mean > c_mean,
                        color="#fddbc7", alpha=0.4)

        # Compute gap summary
        gap = w_mean - c_mean
        mask_h50 = THRESHOLDS >= 0.5
        gap_h50 = gap[mask_h50].mean()

        ax.axhline(0, color="gray", lw=0.5, alpha=0.5)
        ax.axvline(0.5, color="gray", ls=":", lw=1, alpha=0.5)

        ax.set_title(f"{BAND_TEX[bi]}   (gap h≥0.5: {gap_h50:+.3f})",
                     fontsize=13, fontweight="bold")
        ax.set_xlim(0, 1.02)
        ax.set_ylim(-0.1, 1.05)
        ax.grid(alpha=0.2)

        if bi == 0:
            ax.legend(fontsize=10, loc="lower right")

        # Add per-patient thin lines for transparency
        for i, pat in enumerate(pats_used):
            ax.plot(THRESHOLDS, W[i], color="#2166ac", lw=0.5, alpha=0.3)
            ax.plot(THRESHOLDS, C[i], color="#b2182b", lw=0.5, alpha=0.3)

    # Shared labels
    for ax in axes[1, :]:
        ax.set_xlabel("Normalized threshold (h)", fontsize=11)
    for ax in axes[:, 0]:
        ax.set_ylabel("ARI", fontsize=11)

    fig.suptitle(
        "Multiscale ARI: Within-condition vs Cross-condition\n"
        "Blue = mean ARI for rest_pre↔rest_post & task_learn↔task_test  |  "
        "Red = mean ARI for rest↔task pairs\n"
        f"({len([p for p in PATIENTS if any((p,b) in results for b in BANDS)])} patients, "
        f"shaded = ±1 std across patients)",
        fontsize=13, fontweight="bold", y=1.03,
    )
    plt.tight_layout()
    path = OUTPUT_DIR / "reorganization_ari_curves.pdf"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


# ===================================================================
# Figure 2: Summary bar chart — reorganization gap per band
# ===================================================================
def plot_gap_summary(results):
    """Bar chart: reorganization gap = mean(ARI_within - ARI_cross) per band.

    Shows gap at full range and at h>=0.5, plus per-patient dots.
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    x = np.arange(len(BANDS))

    for ax_idx, (h_lo, range_label) in enumerate([(0.01, "Full range (h ∈ [0, 1])"),
                                                    (0.5, "Macro range (h ≥ 0.5)")]):
        ax = axes[ax_idx]
        mask = THRESHOLDS >= h_lo

        band_gaps = []
        band_gaps_per_pat = []
        for band in BANDS:
            pat_gaps = []
            for pat in PATIENTS:
                key = (pat, band)
                if key in results:
                    gap_curve = results[key]["within"] - results[key]["cross"]
                    pat_gaps.append(gap_curve[mask].mean())
            band_gaps.append(np.mean(pat_gaps) if pat_gaps else 0)
            band_gaps_per_pat.append(pat_gaps)

        # Bars
        colors = ["#b2182b" if g < 0 else "#2166ac" for g in band_gaps]
        ax.bar(x, band_gaps, color=colors, alpha=0.7, edgecolor="black", linewidth=0.5)

        # Per-patient dots
        for j, (band, pat_gaps) in enumerate(zip(BANDS, band_gaps_per_pat)):
            pats_used = [p for p in PATIENTS if (p, band) in results]
            for k, (pat, g) in enumerate(zip(pats_used, pat_gaps)):
                ax.scatter(j + (k - len(pats_used)/2 + 0.5) * 0.08, g,
                           s=50, color=PAT_COLORS[pat], edgecolors="black",
                           linewidths=0.5, zorder=5)

        ax.axhline(0, color="black", lw=1)
        ax.set_xticks(x)
        ax.set_xticklabels(BAND_TEX, fontsize=12)
        ax.set_ylabel("Reorganization gap\n(ARI_within − ARI_cross)", fontsize=11)
        ax.set_title(range_label, fontsize=13, fontweight="bold")
        ax.grid(alpha=0.2, axis="y")

        # Annotate bars
        for j, g in enumerate(band_gaps):
            ax.text(j, g + 0.01 * (1 if g >= 0 else -1),
                    f"{g:+.3f}", ha="center", va="bottom" if g >= 0 else "top",
                    fontsize=9, fontweight="bold")

    # Legend
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=PAT_COLORS[p],
                       markeredgecolor='black', markersize=8, label=p)
               for p in PATIENTS if any((p, b) in results for b in BANDS)]
    axes[1].legend(handles=handles, fontsize=8, loc="upper right")

    fig.suptitle(
        "Reorganization Gap per Band\n"
        "Positive = within-condition more similar than cross-condition = REORGANIZATION",
        fontsize=14, fontweight="bold", y=1.04,
    )
    plt.tight_layout()
    path = OUTPUT_DIR / "reorganization_gap_summary.pdf"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


# ===================================================================
# Figure 3: Per-patient gap curves
# ===================================================================
def plot_per_patient_gaps(results):
    """For each band, show gap(h) = ARI_within(h) - ARI_cross(h) per patient."""
    fig, axes = plt.subplots(2, 3, figsize=(20, 12), sharex=True, sharey=True)
    axes_flat = axes.ravel()

    for bi, band in enumerate(BANDS):
        ax = axes_flat[bi]

        for pat in PATIENTS:
            key = (pat, band)
            if key not in results:
                continue
            gap = results[key]["within"] - results[key]["cross"]
            ax.plot(THRESHOLDS, gap, color=PAT_COLORS[pat], lw=1.5,
                    label=pat, alpha=0.8)

        ax.axhline(0, color="black", lw=1, alpha=0.7)
        ax.axvline(0.5, color="gray", ls=":", lw=1, alpha=0.5)
        ax.fill_between([0, 1.02], 0, 1, alpha=0.03, color="blue")
        ax.fill_between([0, 1.02], -1, 0, alpha=0.03, color="red")

        ax.set_title(BAND_TEX[bi], fontsize=13, fontweight="bold")
        ax.set_xlim(0, 1.02)
        ax.set_ylim(-0.6, 0.8)
        ax.grid(alpha=0.2)

        if bi == 0:
            ax.legend(fontsize=9, loc="upper left")

    for ax in axes[1, :]:
        ax.set_xlabel("Normalized threshold (h)", fontsize=11)
    for ax in axes[:, 0]:
        ax.set_ylabel("Gap (within − cross)", fontsize=11)

    fig.suptitle(
        "Per-Patient Reorganization Gap  (ARI_within − ARI_cross)\n"
        "Above 0 (blue zone) = reorganization  |  Below 0 (red zone) = no reorganization",
        fontsize=14, fontweight="bold", y=1.03,
    )
    plt.tight_layout()
    path = OUTPUT_DIR / "reorganization_gap_per_patient.pdf"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


# ===================================================================
# Figure 4: Individual pair comparison — 3 key pairs per band
# ===================================================================
def plot_key_pairs(results):
    """For each band: bar chart of ARI at key pairs, showing per-patient values.

    Three pairs: rest_pre↔rest_post (rest), task_learn↔task_test (task), rest↔task mean (cross).
    """
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    axes_flat = axes.ravel()
    x = np.arange(3)
    pair_labels = ["rest_pre↔rest_post\n(rest stability)", "tLearn↔tTest\n(task stability)",
                   "Cross mean\n(rest↔task)"]

    for bi, band in enumerate(BANDS):
        ax = axes_flat[bi]

        for pi, pat in enumerate(PATIENTS):
            key = (pat, band)
            if key not in results:
                continue

            r = results[key]
            # Compute mean ARI for each category
            rest_ari = None
            task_ari = None
            for curve_list, pair_type in [(r["within_all"], "within"),
                                          (r["cross_all"], "cross")]:
                pass  # handled below

            # within_all contains: rest_pre↔rest_post and task_learn↔task_test (if available)
            # For 4-phase patients, within_all has 2 curves. Index 0 = rest_pre↔rest_post.
            # But the ordering depends on combinations() which goes in PHASES order.
            # Since PHASES = [rest_pre, task_learn, task_test, rest_post]:
            #   combinations gives: (rest_pre,task_learn), (rest_pre,task_test), (rest_pre,rest_post),
            #                       (task_learn,task_test), (task_learn,rest_post), (task_test,rest_post)
            # Within: (rest_pre,rest_post) and (task_learn,task_test) => indices in within_all

            # More robust: recompute for specific pairs
            phases = PATIENT_PHASES[pat]
            lrgs = {}
            for ph in phases:
                lrg = load_lrg_result(pat, ph, band, "msc", cache_root=LRG_CACHE)
                if lrg is not None:
                    lrgs[ph] = lrg

            # Rest pair
            if "rest_pre" in lrgs and "rest_post" in lrgs:
                rest_mean_ari = ari_curve(lrgs["rest_pre"].linkage_matrix,
                                          lrgs["rest_post"].linkage_matrix).mean()
            else:
                rest_mean_ari = np.nan

            # Task pair
            if "task_learn" in lrgs and "task_test" in lrgs:
                task_mean_ari = ari_curve(lrgs["task_learn"].linkage_matrix,
                                          lrgs["task_test"].linkage_matrix).mean()
            else:
                task_mean_ari = np.nan

            cross_mean_ari = np.mean(r["cross"])  # mean across all h

            vals = [rest_mean_ari, task_mean_ari, cross_mean_ari]
            offset = (pi - len(PATIENTS)/2 + 0.5) * 0.12
            ax.bar(x + offset, vals, 0.12, color=PAT_COLORS[pat],
                   edgecolor="black", linewidth=0.3, alpha=0.85,
                   label=pat if bi == 0 else None)

        ax.set_xticks(x)
        ax.set_xticklabels(pair_labels, fontsize=9)
        ax.set_ylabel("Mean ARI", fontsize=10)
        ax.set_title(BAND_TEX[bi], fontsize=13, fontweight="bold")
        ax.set_ylim(0, 1.05)
        ax.grid(alpha=0.2, axis="y")

    axes_flat[0].legend(fontsize=8, loc="upper right")

    fig.suptitle(
        "Mean ARI per Phase Pair — within vs cross\n"
        "Higher = more similar structure | Rest/Task bars should be higher than Cross for reorganization",
        fontsize=14, fontweight="bold", y=1.03,
    )
    plt.tight_layout()
    path = OUTPUT_DIR / "reorganization_key_pairs.pdf"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")


# ===================================================================
# Main
# ===================================================================
def main():
    print("Computing ARI curves for all patients/bands...")
    results = compute_all()
    print(f"  Computed for {len(results)} (patient, band) combinations")

    pats_in = sorted(set(p for p, _ in results.keys()))
    print(f"  Patients with data: {pats_in}")

    print("\n--- Figure 1: ARI curves within vs cross ---")
    plot_ari_curves(results)

    print("\n--- Figure 2: Reorganization gap summary ---")
    plot_gap_summary(results)

    print("\n--- Figure 3: Per-patient gap curves ---")
    plot_per_patient_gaps(results)

    print("\n--- Figure 4: Key pair comparison ---")
    plot_key_pairs(results)

    # Print summary table
    print("\n" + "=" * 70)
    print("REORGANIZATION GAP SUMMARY")
    print("=" * 70)
    print("Gap = mean(ARI_within) - mean(ARI_cross)")
    print("Positive = within-condition more similar = REORGANIZATION")
    print()

    mask_full = THRESHOLDS >= 0.01
    mask_h50 = THRESHOLDS >= 0.5

    print(f"{'Band':<12} {'Gap(full)':>10} {'Gap(h≥.5)':>10} "
          f"{'Within(full)':>12} {'Cross(full)':>12} {'N_patients':>10}")
    print("-" * 68)

    for band in BANDS:
        gaps_full = []
        gaps_h50 = []
        w_full = []
        c_full = []
        for pat in PATIENTS:
            key = (pat, band)
            if key in results:
                w = results[key]["within"]
                c = results[key]["cross"]
                gaps_full.append((w - c)[mask_full].mean())
                gaps_h50.append((w - c)[mask_h50].mean())
                w_full.append(w[mask_full].mean())
                c_full.append(c[mask_full].mean())

        if gaps_full:
            print(f"{band:<12} {np.mean(gaps_full):>+10.3f} {np.mean(gaps_h50):>+10.3f} "
                  f"{np.mean(w_full):>12.3f} {np.mean(c_full):>12.3f} "
                  f"{len(gaps_full):>10d}")

    # Per-patient detail
    print("\n\nPER-PATIENT GAPS (full range):")
    print(f"{'Band':<12}", end="")
    for pat in PATIENTS:
        if any((pat, b) in results for b in BANDS):
            print(f"  {pat:>8}", end="")
    print()
    for band in BANDS:
        print(f"{band:<12}", end="")
        for pat in PATIENTS:
            key = (pat, band)
            if key in results:
                g = (results[key]["within"] - results[key]["cross"]).mean()
                print(f"  {g:>+8.3f}", end="")
            elif any((pat, b) in results for b in BANDS):
                print(f"  {'N/A':>8}", end="")
        print()


if __name__ == "__main__":
    main()
