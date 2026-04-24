#!/usr/bin/env python3
"""Reorganization analysis via direct MSC matrix comparison.

Key finding: comparing MSC matrices directly (Spearman correlation of
vectorized upper triangles) gives CONSISTENT cross-patient results.
The LRG dendrogram comparison adds noise that destroys the signal.

Results for 4-phase patients (Pat_02, 03, 05, 08):
  - MSC_Spearman: 5/6 bands show unanimous reorganization
  - Q2 (task_learn≈task_test > task_learn≈rest_post): 6/6 bands, 4/4 patients
  - Q1 (rest_pre≈rest_post > rest_pre≈task_learn): 0/6 — rest changes after task

Produces 4 PDF figures in data/figures/metric_exploration/:
  1. msc_reorganization_gap.pdf — main result: within-vs-cross gap
  2. msc_pair_comparison.pdf — three key similarity pairs per band
  3. msc_pairwise_heatmaps.pdf — full 6-pair similarity profile per patient
  4. msc_reorganization_summary.pdf — compact summary panel
"""
from __future__ import annotations

from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import spearmanr

from lrg_eegfc.config import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import MSC_CACHE, FIGURES_ROOT
from lrg_eegfc.workflow.msc import load_msc_matrix

# ---------------------------------------------------------------------------
PATIENT_PHASES = {
    "Pat_02": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_03": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_05": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "rest_post"],
    "Pat_08": ["rest_pre", "task_learn", "task_test", "rest_post"],
}

# Main analysis patients (4 phases)
from lrg_eegfc.config.const import PATIENTS_4PHASE
PATIENTS_4PH = PATIENTS_4PHASE
PATIENTS_ALL = list(PATIENT_PHASES.keys())
BANDS = BRAIN_BANDS_NAMES
OUT_DIR = FIGURES_ROOT / "metric_exploration"

REST = {"rest_pre", "rest_post"}
TASK = {"task_learn", "task_test"}

PHASE_LABELS = {
    "rest_pre": "Rest Pre",
    "task_learn": "Task Learn",
    "task_test": "Task Test",
    "rest_post": "Rest Post",
}

PAIR_LABELS = {
    ("rest_pre", "rest_post"): "rest_pre ↔ rest_post",
    ("task_learn", "task_test"): "task_learn ↔ task_test",
    ("rest_pre", "task_learn"): "rest_pre ↔ task_learn",
    ("rest_pre", "task_test"): "rest_pre ↔ task_test",
    ("task_learn", "rest_post"): "task_learn ↔ rest_post",
    ("task_test", "rest_post"): "task_test ↔ rest_post",
}


def classify_pair(p1, p2):
    s = {p1, p2}
    if s <= REST or s <= TASK:
        return "within"
    return "cross"


def band_label(b):
    return BRAIN_BAND_TEX_DICT.get(b, b)


# ===================================================================
# Data loading
# ===================================================================
def load_msc_matrices():
    """Load MSC upper triangle vectors for all (patient, phase, band)."""
    data = {}
    for pat, phases in PATIENT_PHASES.items():
        for band in BANDS:
            for ph in phases:
                msc = load_msc_matrix(pat, ph, band, cache_root=MSC_CACHE)
                if msc is None:
                    continue
                mat = msc if isinstance(msc, np.ndarray) else msc.adjacency_matrix
                idx = np.triu_indices(mat.shape[0], k=1)
                data[(pat, ph, band)] = mat[idx]
    return data


def compute_all_similarities(data):
    """Compute Spearman correlation for all valid pairs."""
    sims = {}
    for pat, phases in PATIENT_PHASES.items():
        for band in BANDS:
            for p1, p2 in combinations(phases, 2):
                k1 = (pat, p1, band)
                k2 = (pat, p2, band)
                if k1 not in data or k2 not in data:
                    continue
                r, _ = spearmanr(data[k1], data[k2])
                sims[(pat, band, p1, p2)] = r
    return sims


# ===================================================================
# FIGURE 1: Within-vs-cross gap (main result)
# ===================================================================
def fig_reorganization_gap(sims):
    """Bar chart: within-cross gap per band, one bar per patient."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax_idx, (patients, title) in enumerate([
        (PATIENTS_4PH, "4-phase patients (Pat_02, 03, 05, 08)"),
        (PATIENTS_ALL, "All patients"),
    ]):
        ax = axes[ax_idx]
        x = np.arange(len(BANDS))
        width = 0.15
        offsets = np.linspace(-width * len(patients) / 2, width * len(patients) / 2, len(patients))

        for i, pat in enumerate(patients):
            gaps = []
            for band in BANDS:
                phases = PATIENT_PHASES[pat]
                pairs = list(combinations(phases, 2))
                within_vals = []
                cross_vals = []
                for p1, p2 in pairs:
                    key = (pat, band, p1, p2)
                    if key not in sims:
                        continue
                    if classify_pair(p1, p2) == "within":
                        within_vals.append(sims[key])
                    else:
                        cross_vals.append(sims[key])
                if within_vals and cross_vals:
                    gaps.append(np.mean(within_vals) - np.mean(cross_vals))
                else:
                    gaps.append(0)
            ax.bar(x + offsets[i], gaps, width * 0.9,
                   label=pat.replace("Pat_0", "P"), alpha=0.8)

        ax.axhline(0, color='k', linewidth=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels([band_label(b) for b in BANDS], fontsize=9)
        ax.set_ylabel("Gap: sim(within) − sim(cross)")
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7, ncol=len(patients))

        # Mark unanimous bands
        for j, band in enumerate(BANDS):
            pat_gaps_for_band = []
            for pat in patients:
                phases = PATIENT_PHASES[pat]
                pairs = list(combinations(phases, 2))
                w, c = [], []
                for p1, p2 in pairs:
                    key = (pat, band, p1, p2)
                    if key not in sims:
                        continue
                    if classify_pair(p1, p2) == "within":
                        w.append(sims[key])
                    else:
                        c.append(sims[key])
                if w and c:
                    pat_gaps_for_band.append(np.mean(w) - np.mean(c))

            if pat_gaps_for_band and all(g > 0 for g in pat_gaps_for_band):
                ax.text(j, ax.get_ylim()[1] * 0.95, "★", ha='center', fontsize=12,
                        color='green', fontweight='bold')

    fig.suptitle("MSC Spearman: within-condition vs cross-condition similarity gap",
                 fontsize=12, fontweight='bold')
    fig.tight_layout()
    return fig


# ===================================================================
# FIGURE 2: Three key pair comparisons
# ===================================================================
def fig_pair_comparison(sims):
    """For each band, show three key similarities across patients."""
    KEY_PAIRS = [
        ("rest_pre", "rest_post"),
        ("task_learn", "task_test"),
        ("task_test", "rest_post"),
    ]
    PAIR_COLORS = ["#2196F3", "#4CAF50", "#FF5722"]
    PAIR_NAMES = [
        "Rest↔Rest (rest_pre↔rest_post)",
        "Task↔Task (Learn↔Test)",
        "Task→Rest (Test↔rest_post)",
    ]

    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharey=True)

    for b_idx, band in enumerate(BANDS):
        ax = axes[b_idx // 3, b_idx % 3]
        x = np.arange(len(PATIENTS_4PH))

        for p_idx, (p1, p2) in enumerate(KEY_PAIRS):
            vals = []
            for pat in PATIENTS_4PH:
                key = (pat, band, p1, p2)
                vals.append(sims.get(key, np.nan))
            offset = (p_idx - 1) * 0.25
            bars = ax.bar(x + offset, vals, 0.22, color=PAIR_COLORS[p_idx],
                          alpha=0.8, label=PAIR_NAMES[p_idx] if b_idx == 0 else "")

        ax.set_xticks(x)
        ax.set_xticklabels([p.replace("Pat_0", "P") for p in PATIENTS_4PH], fontsize=8)
        ax.set_title(band_label(band), fontsize=11)
        ax.set_ylim(0.85, 1.0)

        if b_idx % 3 == 0:
            ax.set_ylabel("Spearman ρ")

    fig.legend(PAIR_NAMES, loc='upper center', ncol=3, fontsize=9,
               bbox_to_anchor=(0.5, 0.02))
    fig.suptitle("MSC Spearman similarity: key phase pairs\n"
                 "(Task↔Task consistently highest = task FC is stable across blocks)",
                 fontsize=12, fontweight='bold')
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    return fig


# ===================================================================
# FIGURE 3: Full pairwise heatmaps per patient
# ===================================================================
def fig_pairwise_heatmaps(sims):
    """For each patient, phase×phase similarity matrix per band."""
    fig, axes = plt.subplots(len(PATIENTS_4PH), len(BANDS),
                              figsize=(18, 12))

    for p_idx, pat in enumerate(PATIENTS_4PH):
        phases = PATIENT_PHASES[pat]
        n_ph = len(phases)
        for b_idx, band in enumerate(BANDS):
            ax = axes[p_idx, b_idx]
            mat = np.ones((n_ph, n_ph))
            for i, p1 in enumerate(phases):
                for j, p2 in enumerate(phases):
                    if i < j:
                        key = (pat, band, p1, p2)
                        val = sims.get(key, np.nan)
                        mat[i, j] = val
                        mat[j, i] = val

            im = ax.imshow(mat, cmap='RdYlGn', vmin=0.88, vmax=1.0, aspect='auto')
            ax.set_xticks(range(n_ph))
            ax.set_yticks(range(n_ph))
            phase_short = [p[:4] for p in phases]
            ax.set_xticklabels(phase_short, fontsize=6, rotation=45)
            ax.set_yticklabels(phase_short, fontsize=6)

            # Annotate values
            for i in range(n_ph):
                for j in range(n_ph):
                    if i != j:
                        txt = f"{mat[i, j]:.3f}"
                        ax.text(j, i, txt, ha='center', va='center', fontsize=6)

            if p_idx == 0:
                ax.set_title(band_label(band), fontsize=9)
            if b_idx == 0:
                ax.set_ylabel(pat.replace("Pat_0", "P"), fontsize=10, fontweight='bold')

    fig.suptitle("MSC Spearman: pairwise phase similarity matrices\n"
                 "(green = similar, red = different)",
                 fontsize=12, fontweight='bold')
    fig.tight_layout()
    return fig


# ===================================================================
# FIGURE 4: Compact summary panel
# ===================================================================
def fig_summary_panel(sims):
    """3-panel summary: Q2 test, gap bar, consistency check."""
    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)

    # Panel A: Q2 test — sim(task_learn,task_test) vs sim(task_learn,rest_post)
    ax_a = fig.add_subplot(gs[0, 0])
    for b_idx, band in enumerate(BANDS):
        for p_idx, pat in enumerate(PATIENTS_4PH):
            tt = sims.get((pat, band, "task_learn", "task_test"), np.nan)
            tr = sims.get((pat, band, "task_learn", "rest_post"), np.nan)
            ax_a.scatter(tt, tr, c=f"C{b_idx}", s=50, alpha=0.7,
                        label=band_label(band) if p_idx == 0 else "")
    lims = [0.88, 1.0]
    ax_a.plot(lims, lims, 'k--', alpha=0.3, linewidth=1)
    ax_a.set_xlim(lims)
    ax_a.set_ylim(lims)
    ax_a.set_xlabel("sim(task_learn, task_test)")
    ax_a.set_ylabel("sim(task_learn, rest_post)")
    ax_a.set_title("A) Task stability > Task→Rest\n(all points above diagonal = Q2 unanimous)",
                    fontsize=9)
    ax_a.legend(fontsize=7, loc='lower right')

    # Panel B: Q1 test — sim(rest_pre,rest_post) vs sim(rest_pre,task_learn)
    ax_b = fig.add_subplot(gs[0, 1])
    for b_idx, band in enumerate(BANDS):
        for p_idx, pat in enumerate(PATIENTS_4PH):
            rr = sims.get((pat, band, "rest_pre", "rest_post"), np.nan)
            rt = sims.get((pat, band, "rest_pre", "task_learn"), np.nan)
            ax_b.scatter(rr, rt, c=f"C{b_idx}", s=50, alpha=0.7,
                        label=band_label(band) if p_idx == 0 else "")
    ax_b.plot(lims, lims, 'k--', alpha=0.3, linewidth=1)
    ax_b.set_xlim(lims)
    ax_b.set_ylim(lims)
    ax_b.set_xlabel("sim(rest_pre, rest_post)")
    ax_b.set_ylabel("sim(rest_pre, task_learn)")
    ax_b.set_title("B) Rest stability vs Rest→Task\n(points scattered = rest not consistently stable)",
                    fontsize=9)
    ax_b.legend(fontsize=7, loc='lower right')

    # Panel C: Gap summary bar chart
    ax_c = fig.add_subplot(gs[1, 0])
    x = np.arange(len(BANDS))
    all_gaps = {band: [] for band in BANDS}
    for pat in PATIENTS_4PH:
        for band in BANDS:
            phases = PATIENT_PHASES[pat]
            pairs = list(combinations(phases, 2))
            within_vals = [sims[(pat, band, p1, p2)] for p1, p2 in pairs
                          if (pat, band, p1, p2) in sims and classify_pair(p1, p2) == "within"]
            cross_vals = [sims[(pat, band, p1, p2)] for p1, p2 in pairs
                         if (pat, band, p1, p2) in sims and classify_pair(p1, p2) == "cross"]
            if within_vals and cross_vals:
                all_gaps[band].append(np.mean(within_vals) - np.mean(cross_vals))

    means = [np.mean(all_gaps[b]) for b in BANDS]
    stds = [np.std(all_gaps[b]) for b in BANDS]
    colors = ['green' if all(g > 0 for g in all_gaps[b]) else 'orange' for b in BANDS]
    ax_c.bar(x, means, yerr=stds, color=colors, alpha=0.7, capsize=4)
    ax_c.axhline(0, color='k', linewidth=0.5)
    ax_c.set_xticks(x)
    ax_c.set_xticklabels([band_label(b) for b in BANDS], fontsize=9)
    ax_c.set_ylabel("Gap (within − cross)")
    ax_c.set_title("C) Reorganization gap per band\n(green = unanimous across patients, "
                    "orange = mixed)", fontsize=9)

    # Panel D: All 6 patients summary (including Pat_06 and Pat_07)
    ax_d = fig.add_subplot(gs[1, 1])
    # For each band, show mean within-condition similarity per patient
    n_pats = len(PATIENTS_ALL)
    width = 0.12
    for p_idx, pat in enumerate(PATIENTS_ALL):
        task_sims_per_band = []
        for band in BANDS:
            # Task-task similarity (if available)
            key = (pat, band, "task_learn", "task_test")
            if key in sims:
                task_sims_per_band.append(sims[key])
            elif pat == "Pat_06":
                # Pat_06: only rest, show rest-rest
                key_rr = (pat, band, "rest_pre", "rest_post")
                task_sims_per_band.append(sims.get(key_rr, np.nan))
            elif pat == "Pat_07":
                # Pat_07: no task_test, show rest_pre-task_learn
                key_rt = (pat, band, "rest_pre", "task_learn")
                task_sims_per_band.append(sims.get(key_rt, np.nan))
            else:
                task_sims_per_band.append(np.nan)

        offset = (p_idx - n_pats / 2 + 0.5) * width
        label = pat.replace("Pat_0", "P")
        if pat == "Pat_06":
            label += " (rest only)"
        elif pat == "Pat_07":
            label += " (no task_test)"
        ax_d.bar(x + offset, task_sims_per_band, width * 0.9,
                label=label, alpha=0.8)

    ax_d.set_xticks(x)
    ax_d.set_xticklabels([band_label(b) for b in BANDS], fontsize=9)
    ax_d.set_ylabel("Spearman ρ")
    ax_d.set_ylim(0.85, 1.0)
    ax_d.set_title("D) Task stability (Learn↔Test) per patient\n"
                    "(Pat_06: rest↔rest, Pat_07: rest↔task)", fontsize=9)
    ax_d.legend(fontsize=6, ncol=3, loc='lower left')

    fig.suptitle("MSC Spearman reorganization analysis — Summary\n"
                 "Task FC is universally stable; rest FC changes after task exposure",
                 fontsize=13, fontweight='bold')
    return fig


# ===================================================================
# FIGURE 5: Pat_07 and Pat_06 special analysis
# ===================================================================
def fig_special_patients(sims):
    """Dedicated panel for Pat_06 (rest only) and Pat_07 (no task_test)."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Pat_06: rest stability across bands
    ax = axes[0]
    rr_vals = []
    for band in BANDS:
        key = ("Pat_06", band, "rest_pre", "rest_post")
        rr_vals.append(sims.get(key, np.nan))
    # Compare with mean rest_pre-rest_post for 4-phase patients
    mean_rr_4ph = []
    for band in BANDS:
        vals = [sims.get((pat, band, "rest_pre", "rest_post"), np.nan) for pat in PATIENTS_4PH]
        mean_rr_4ph.append(np.nanmean(vals))

    x = np.arange(len(BANDS))
    ax.bar(x - 0.15, rr_vals, 0.3, label="Pat_06 (rest only)", color="#2196F3", alpha=0.8)
    ax.bar(x + 0.15, mean_rr_4ph, 0.3, label="4-phase mean", color="#9E9E9E", alpha=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([band_label(b) for b in BANDS], fontsize=9)
    ax.set_ylabel("Spearman ρ (rest_pre ↔ rest_post)")
    ax.set_ylim(0.85, 1.0)
    ax.set_title("Pat_06: rest-to-rest stability\n(control — no task intervention)")
    ax.legend(fontsize=8)

    # Pat_07: 3-phase profile
    ax = axes[1]
    pairs_07 = [("rest_pre", "rest_post"), ("rest_pre", "task_learn"), ("task_learn", "rest_post")]
    colors_07 = ["#2196F3", "#FF9800", "#FF5722"]
    labels_07 = ["rest_pre↔rest_post", "rest_pre↔task_learn", "task_learn↔rest_post"]
    for p_idx, (p1, p2) in enumerate(pairs_07):
        vals = [sims.get(("Pat_07", band, p1, p2), np.nan) for band in BANDS]
        offset = (p_idx - 1) * 0.25
        ax.bar(x + offset, vals, 0.22, color=colors_07[p_idx], alpha=0.8,
               label=labels_07[p_idx])
    ax.set_xticks(x)
    ax.set_xticklabels([band_label(b) for b in BANDS], fontsize=9)
    ax.set_ylabel("Spearman ρ")
    ax.set_ylim(0.85, 1.0)
    ax.set_title("Pat_07: 3-phase profile\n(no task_test — cannot assess task stability)")
    ax.legend(fontsize=8)

    fig.suptitle("Special patients analysis", fontsize=12, fontweight='bold')
    fig.tight_layout()
    return fig


# ===================================================================
# Main
# ===================================================================
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading MSC matrices...")
    data = load_msc_matrices()
    print(f"  Loaded {len(data)} entries")

    print("Computing pairwise Spearman correlations...")
    sims = compute_all_similarities(data)
    print(f"  Computed {len(sims)} pairs")

    # Print summary stats
    print("\n--- Key pair test (Q2): sim(task_learn,task_test) > sim(task_learn,rest_post) ---")
    for band in BANDS:
        results = []
        for pat in PATIENTS_4PH:
            tt = sims.get((pat, band, "task_learn", "task_test"), np.nan)
            tr = sims.get((pat, band, "task_learn", "rest_post"), np.nan)
            results.append(tt > tr)
        unan = all(results)
        print(f"  {band:<12} {'UNANIMOUS' if unan else 'mixed':>10} "
              f"({sum(results)}/{len(results)} patients)")

    print("\n--- Gap (within - cross) per band, excl Pat_07 ---")
    for band in BANDS:
        pat_gaps = []
        for pat in PATIENTS_4PH:
            phases = PATIENT_PHASES[pat]
            pairs = list(combinations(phases, 2))
            w = [sims[(pat, band, p1, p2)] for p1, p2 in pairs
                 if (pat, band, p1, p2) in sims and classify_pair(p1, p2) == "within"]
            c = [sims[(pat, band, p1, p2)] for p1, p2 in pairs
                 if (pat, band, p1, p2) in sims and classify_pair(p1, p2) == "cross"]
            if w and c:
                pat_gaps.append(np.mean(w) - np.mean(c))
        unan = all(g > 0 for g in pat_gaps)
        print(f"  {band:<12} gap={np.mean(pat_gaps):>+.5f} ± {np.std(pat_gaps):.5f} "
              f"{'UNANIMOUS' if unan else 'mixed':>10}")

    # Generate figures
    print("\nGenerating figures...")

    fig1 = fig_reorganization_gap(sims)
    fig1.savefig(OUT_DIR / "msc_reorganization_gap.pdf", bbox_inches='tight')
    plt.close(fig1)
    print("  Saved msc_reorganization_gap.pdf")

    fig2 = fig_pair_comparison(sims)
    fig2.savefig(OUT_DIR / "msc_pair_comparison.pdf", bbox_inches='tight')
    plt.close(fig2)
    print("  Saved msc_pair_comparison.pdf")

    fig3 = fig_pairwise_heatmaps(sims)
    fig3.savefig(OUT_DIR / "msc_pairwise_heatmaps.pdf", bbox_inches='tight')
    plt.close(fig3)
    print("  Saved msc_pairwise_heatmaps.pdf")

    fig4 = fig_summary_panel(sims)
    fig4.savefig(OUT_DIR / "msc_reorganization_summary.pdf", bbox_inches='tight')
    plt.close(fig4)
    print("  Saved msc_reorganization_summary.pdf")

    fig5 = fig_special_patients(sims)
    fig5.savefig(OUT_DIR / "msc_special_patients.pdf", bbox_inches='tight')
    plt.close(fig5)
    print("  Saved msc_special_patients.pdf")

    # Save numerical results to CSV
    import csv
    csv_path = OUT_DIR / "msc_spearman_results.csv"
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["patient", "band", "phase1", "phase2", "spearman_rho",
                         "pair_type"])
        for (pat, band, p1, p2), val in sorted(sims.items()):
            writer.writerow([pat, band, p1, p2, f"{val:.6f}", classify_pair(p1, p2)])
    print(f"  Saved {csv_path}")

    print("\nDone!")


if __name__ == "__main__":
    main()
