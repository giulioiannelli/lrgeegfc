#!/usr/bin/env python3
"""Hierarchical Bridge Nodes Analysis.

For each node, track community membership as k (number of clusters) sweeps
from 2 to 20 in the LRG dendrogram.  Nodes that change community across
consecutive k values (after Hungarian label alignment) are "bridge" nodes.

Bridge score = fraction of scale transitions where the node switches community.

Outputs
-------
data/figures/multiscale_investigation/bridge_nodes/
    fig1_bridge_score_histogram.png
    fig2_bridge_score_scatter.png
    fig3_cross_phase_corr_heatmap.png
    fig4_cross_patient_corr.png
    fig5_jaccard_similarity.png
    fig6_consistent_bridges.png
    bridge_scores.csv
    cross_phase_correlations.csv
    cross_patient_correlations.csv
    FINDINGS.md
"""

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy.cluster.hierarchy import fcluster
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr
from pathlib import Path
from itertools import combinations

from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, PHASE_LABELS, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE
from lrg_eegfc.config.paths import FIGURES_ROOT

# ── Configuration ──────────────────────────────────────────────────────
FC_METHOD = "msc"
ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]
FOUR_PHASE_PATIENTS = PATIENTS_4PHASE
ALL_PHASES = list(PHASE_LABELS)
ALL_BANDS = BRAIN_BANDS_NAMES
K_RANGE = range(2, 21)  # sweep k from 2 to 20
TOP_FRAC = 0.20  # top 20% bridge nodes
OUT_DIR = FIGURES_ROOT / "multiscale_investigation" / "bridge_nodes"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Phase availability per patient
PATIENT_PHASES = {
    p: ALL_PHASES if p in FOUR_PHASE_PATIENTS else ["task_learn", "task_test"]
    for p in ALL_PATIENTS
}

# ── Helper functions ───────────────────────────────────────────────────

def align_labels(labels_a, labels_b):
    """Align labels_b to labels_a via Hungarian algorithm on contingency."""
    unique_a = np.unique(labels_a)
    unique_b = np.unique(labels_b)
    n = max(len(unique_a), len(unique_b))
    cost = np.zeros((n, n), dtype=int)
    map_a = {v: i for i, v in enumerate(unique_a)}
    map_b = {v: i for i, v in enumerate(unique_b)}
    for la, lb in zip(labels_a, labels_b):
        cost[map_a[la], map_b[lb]] += 1
    # Maximize overlap = minimize negative overlap
    row_ind, col_ind = linear_sum_assignment(-cost)
    # Build mapping: unique_b[col_ind[j]] -> unique_a[row_ind[j]]
    mapping = {}
    for r, c in zip(row_ind, col_ind):
        if c < len(unique_b):
            if r < len(unique_a):
                mapping[unique_b[c]] = unique_a[r]
            else:
                mapping[unique_b[c]] = unique_b[c] + 1000  # unmapped
    aligned = np.array([mapping.get(lb, lb) for lb in labels_b])
    return aligned


def compute_bridge_scores(linkage_matrix, k_range=K_RANGE):
    """Compute per-node bridge scores from a linkage matrix.

    Returns array of shape (n_nodes,) with bridge score in [0, 1].
    """
    k_values = list(k_range)
    # Get labels at each k
    all_labels = {}
    for k in k_values:
        all_labels[k] = fcluster(linkage_matrix, k, criterion="maxclust")
    n_nodes = len(all_labels[k_values[0]])
    n_transitions = len(k_values) - 1
    change_count = np.zeros(n_nodes, dtype=int)

    for i in range(n_transitions):
        k_curr = k_values[i]
        k_next = k_values[i + 1]
        labels_curr = all_labels[k_curr]
        labels_next = all_labels[k_next]
        # Align labels_next to labels_curr
        labels_next_aligned = align_labels(labels_curr, labels_next)
        changed = labels_curr != labels_next_aligned
        change_count += changed.astype(int)

    bridge_scores = change_count / n_transitions
    return bridge_scores


def jaccard(set_a, set_b):
    """Jaccard similarity between two sets."""
    if len(set_a) == 0 and len(set_b) == 0:
        return 1.0
    inter = len(set_a & set_b)
    union = len(set_a | set_b)
    return inter / union if union > 0 else 0.0


def top_bridge_set(scores, frac=TOP_FRAC):
    """Return set of indices in top frac of bridge scores."""
    threshold = np.percentile(scores, 100 * (1 - frac))
    return set(np.where(scores >= threshold)[0])


# ── Step 1: Compute bridge scores for all patient/band/phase ──────────
print("=" * 70)
print("STEP 1: Computing bridge scores")
print("=" * 70)

# Store results: bridge_data[patient][band][phase] = bridge_scores array
bridge_data = {}
records = []  # for CSV

for patient in ALL_PATIENTS:
    bridge_data[patient] = {}
    phases = PATIENT_PHASES[patient]
    for band in ALL_BANDS:
        bridge_data[patient][band] = {}
        for phase in phases:
            result = load_lrg_result(patient, phase, band, FC_METHOD)
            if result is None:
                print(f"  SKIP {patient} {phase} {band}: no LRG cache")
                continue
            scores = compute_bridge_scores(result.linkage_matrix)
            bridge_data[patient][band][phase] = scores
            n_nodes = len(scores)
            for node_idx, s in enumerate(scores):
                records.append({
                    "patient": patient,
                    "band": band,
                    "phase": phase,
                    "node": node_idx,
                    "bridge_score": s,
                    "n_nodes": n_nodes,
                })
            print(f"  {patient} {phase} {band}: {n_nodes} nodes, "
                  f"mean={scores.mean():.3f}, max={scores.max():.3f}")

df = pd.DataFrame(records)
df.to_csv(OUT_DIR / "bridge_scores.csv", index=False)
print(f"\nSaved bridge_scores.csv ({len(df)} rows)")

# ── Step 2: Cross-phase Spearman correlations ─────────────────────────
print("\n" + "=" * 70)
print("STEP 2: Cross-phase correlations of bridge scores")
print("=" * 70)

cross_phase_records = []
# For heatmap: per band, average across patients
cross_phase_by_band = {band: {} for band in ALL_BANDS}

for band in ALL_BANDS:
    for patient in FOUR_PHASE_PATIENTS:
        bd = bridge_data[patient][band]
        phases_avail = [p for p in ALL_PHASES if p in bd]
        for p1, p2 in combinations(phases_avail, 2):
            s1, s2 = bd[p1], bd[p2]
            if len(s1) != len(s2):
                # Different giant components -- skip
                continue
            rho, pval = spearmanr(s1, s2)
            cross_phase_records.append({
                "patient": patient, "band": band,
                "phase_1": p1, "phase_2": p2,
                "spearman_rho": rho, "p_value": pval,
                "n_nodes": len(s1),
            })

df_xphase = pd.DataFrame(cross_phase_records)
df_xphase.to_csv(OUT_DIR / "cross_phase_correlations.csv", index=False)
print(f"Saved cross_phase_correlations.csv ({len(df_xphase)} rows)")

# Average cross-phase correlations by band, for heatmap
phase_pairs = list(combinations(ALL_PHASES, 2))
heatmap_data = {}
for band in ALL_BANDS:
    sub = df_xphase[df_xphase["band"] == band]
    mat = np.full((4, 4), np.nan)
    for i in range(4):
        mat[i, i] = 1.0
    for i, j in combinations(range(4), 2):
        p1, p2 = ALL_PHASES[i], ALL_PHASES[j]
        vals = sub[(sub["phase_1"] == p1) & (sub["phase_2"] == p2)]["spearman_rho"]
        if len(vals) > 0:
            m = vals.mean()
            mat[i, j] = m
            mat[j, i] = m
    heatmap_data[band] = mat

# ── Step 3: Cross-patient correlations ────────────────────────────────
print("\n" + "=" * 70)
print("STEP 3: Cross-patient correlations of bridge scores")
print("=" * 70)

cross_patient_records = []
for band in ALL_BANDS:
    for phase in ALL_PHASES:
        # Collect patients that have this phase/band with same n_nodes
        available = []
        for patient in ALL_PATIENTS:
            if phase in bridge_data.get(patient, {}).get(band, {}):
                available.append(patient)
        for p1, p2 in combinations(available, 2):
            s1 = bridge_data[p1][band][phase]
            s2 = bridge_data[p2][band][phase]
            if len(s1) != len(s2):
                continue
            rho, pval = spearmanr(s1, s2)
            cross_patient_records.append({
                "patient_1": p1, "patient_2": p2,
                "band": band, "phase": phase,
                "spearman_rho": rho, "p_value": pval,
                "n_nodes": len(s1),
            })

df_xpat = pd.DataFrame(cross_patient_records)
df_xpat.to_csv(OUT_DIR / "cross_patient_correlations.csv", index=False)
print(f"Saved cross_patient_correlations.csv ({len(df_xpat)} rows)")

# ── Step 4: Jaccard similarity of top bridge node sets ────────────────
print("\n" + "=" * 70)
print("STEP 4: Jaccard similarity of top bridge node sets")
print("=" * 70)

jaccard_records = []
for band in ALL_BANDS:
    for patient in FOUR_PHASE_PATIENTS:
        bd = bridge_data[patient][band]
        phases_avail = [p for p in ALL_PHASES if p in bd]
        bridge_sets = {p: top_bridge_set(bd[p]) for p in phases_avail if len(bd[p]) == len(bd[phases_avail[0]])}
        for p1, p2 in combinations(bridge_sets.keys(), 2):
            j = jaccard(bridge_sets[p1], bridge_sets[p2])
            jaccard_records.append({
                "patient": patient, "band": band,
                "phase_1": p1, "phase_2": p2,
                "jaccard": j,
            })

df_jaccard = pd.DataFrame(jaccard_records)

# Classify pairs as within-condition vs cross-condition
rest_phases = {"rest_pre", "rest_post"}
task_phases = {"task_learn", "task_test"}

def pair_type(p1, p2):
    s = {p1, p2}
    if s <= rest_phases:
        return "within-rest"
    elif s <= task_phases:
        return "within-task"
    else:
        return "cross-condition"

if len(df_jaccard) > 0:
    df_jaccard["pair_type"] = df_jaccard.apply(lambda r: pair_type(r["phase_1"], r["phase_2"]), axis=1)

# ── FIGURES ────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("GENERATING FIGURES")
print("=" * 70)

# Color scheme
PHASE_COLORS = {
    "rest_pre": "#2196F3",
    "task_learn": "#FF9800",
    "task_test": "#F44336",
    "rest_post": "#4CAF50",
}
BAND_COLORS = plt.cm.viridis(np.linspace(0.1, 0.9, len(ALL_BANDS)))

# ── Fig 1: Bridge score histograms ────────────────────────────────────
print("  Fig 1: Bridge score histograms (Pat_02, alpha)")
fig, ax = plt.subplots(figsize=(8, 5))
band_demo = "alpha"
patient_demo = "Pat_02"
bd = bridge_data[patient_demo][band_demo]
bins = np.linspace(0, 1, 25)
for phase in ALL_PHASES:
    if phase in bd:
        ax.hist(bd[phase], bins=bins, alpha=0.45, label=phase,
                color=PHASE_COLORS[phase], edgecolor="white", linewidth=0.5)
ax.set_xlabel("Bridge Score", fontsize=13)
ax.set_ylabel("Number of Nodes", fontsize=13)
ax.set_title(f"Bridge Score Distribution — {patient_demo}, {BRAIN_BAND_TEX_DICT[band_demo]}",
             fontsize=14)
ax.legend(fontsize=11)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(OUT_DIR / "fig1_bridge_score_histogram.png", dpi=200)
plt.close(fig)

# ── Fig 2: Scatter plot rest_pre vs task_learn, rest_pre vs rest_post ──────────
print("  Fig 2: Bridge score scatter")
fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
pairs_to_plot = [("rest_pre", "task_learn"), ("rest_pre", "rest_post")]
for ax, (pa, pb) in zip(axes, pairs_to_plot):
    bd = bridge_data[patient_demo][band_demo]
    if pa in bd and pb in bd and len(bd[pa]) == len(bd[pb]):
        sa, sb = bd[pa], bd[pb]
        # Color by quartile of average
        avg = (sa + sb) / 2
        quartile = np.digitize(avg, np.percentile(avg, [25, 50, 75]))
        colors = plt.cm.RdYlBu_r(quartile / 3.0)
        ax.scatter(sa, sb, c=colors, s=20, alpha=0.7, edgecolors="none")
        rho, pval = spearmanr(sa, sb)
        ax.set_xlabel(f"Bridge Score ({pa})", fontsize=12)
        ax.set_ylabel(f"Bridge Score ({pb})", fontsize=12)
        ax.set_title(f"{pa} vs {pb}\n" + r"Spearman $\rho$" + f" = {rho:.3f}, p = {pval:.2e}",
                     fontsize=12)
        ax.plot([0, 1], [0, 1], "k--", alpha=0.3, linewidth=0.8)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_aspect("equal")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
fig.suptitle(f"Bridge Score Comparison — {patient_demo}, {BRAIN_BAND_TEX_DICT[band_demo]}",
             fontsize=14, y=1.02)
fig.tight_layout()
fig.savefig(OUT_DIR / "fig2_bridge_score_scatter.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# ── Fig 3: Cross-phase Spearman heatmap per band (avg across patients) ─
print("  Fig 3: Cross-phase correlation heatmap")
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes_flat = axes.flat
phase_short = ["rest_pre", "tLearn", "tTest", "rest_post"]
for idx, band in enumerate(ALL_BANDS):
    ax = axes_flat[idx]
    mat = heatmap_data[band]
    im = ax.imshow(mat, vmin=-0.2, vmax=1.0, cmap="RdBu_r", aspect="equal")
    ax.set_xticks(range(4))
    ax.set_xticklabels(phase_short, fontsize=9, rotation=45, ha="right")
    ax.set_yticks(range(4))
    ax.set_yticklabels(phase_short, fontsize=9)
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=13)
    for i in range(4):
        for j in range(4):
            if not np.isnan(mat[i, j]):
                ax.text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center",
                        fontsize=8, color="white" if mat[i, j] > 0.6 or mat[i, j] < -0.1 else "black")
fig.suptitle("Cross-Phase Spearman Correlation of Bridge Scores\n(averaged across 4-phase patients)",
             fontsize=14)
fig.colorbar(im, ax=list(axes.flat), shrink=0.6, label=r"Spearman $\rho$")
fig.tight_layout(rect=[0, 0, 0.92, 0.93])
fig.savefig(OUT_DIR / "fig3_cross_phase_corr_heatmap.png", dpi=200, bbox_inches="tight")
plt.close(fig)

# ── Fig 4: Cross-patient correlation ─────────────────────────────────
print("  Fig 4: Cross-patient correlation")
if len(df_xpat) > 0:
    fig, ax = plt.subplots(figsize=(10, 5))
    # Group by band, show distribution of rho across patient pairs/phases
    band_positions = np.arange(len(ALL_BANDS))
    bp_data = []
    for band in ALL_BANDS:
        vals = df_xpat[df_xpat["band"] == band]["spearman_rho"].values
        bp_data.append(vals)
    parts = ax.violinplot(bp_data, positions=band_positions, showmeans=True, showmedians=True)
    for pc, c in zip(parts["bodies"], BAND_COLORS):
        pc.set_facecolor(c)
        pc.set_alpha(0.6)
    ax.set_xticks(band_positions)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in ALL_BANDS], fontsize=12)
    ax.set_ylabel(r"Spearman $\rho$", fontsize=13)
    ax.set_title("Cross-Patient Correlation of Bridge Scores\n(all phase/patient pairs with matching node count)",
                 fontsize=13)
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig4_cross_patient_corr.png", dpi=200)
    plt.close(fig)
else:
    print("    No cross-patient pairs with matching node count -- skipping Fig 4")

# ── Fig 5: Jaccard similarity ────────────────────────────────────────
print("  Fig 5: Jaccard similarity")
if len(df_jaccard) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # 5a: By band, colored by pair type
    ax = axes[0]
    pair_types = ["within-rest", "within-task", "cross-condition"]
    type_colors = {"within-rest": "#2196F3", "within-task": "#FF9800", "cross-condition": "#9C27B0"}
    x_pos = np.arange(len(ALL_BANDS))
    width = 0.25
    for i, pt in enumerate(pair_types):
        means = []
        sems = []
        for band in ALL_BANDS:
            vals = df_jaccard[(df_jaccard["band"] == band) & (df_jaccard["pair_type"] == pt)]["jaccard"]
            means.append(vals.mean() if len(vals) > 0 else 0)
            sems.append(vals.std() / np.sqrt(len(vals)) if len(vals) > 1 else 0)
        ax.bar(x_pos + (i - 1) * width, means, width, yerr=sems,
               label=pt, color=type_colors[pt], alpha=0.8, capsize=3)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in ALL_BANDS], fontsize=11)
    ax.set_ylabel("Jaccard Similarity", fontsize=12)
    ax.set_title("Top 20% Bridge Nodes Overlap\nby Condition Type", fontsize=12)
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # 5b: Overall distribution by pair type
    ax = axes[1]
    for pt in pair_types:
        vals = df_jaccard[df_jaccard["pair_type"] == pt]["jaccard"].values
        if len(vals) > 0:
            ax.hist(vals, bins=15, alpha=0.5, label=f"{pt} (n={len(vals)})",
                    color=type_colors[pt], edgecolor="white")
    ax.set_xlabel("Jaccard Similarity", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("Distribution of Jaccard Similarities\nacross All Bands and Patients", fontsize=12)
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig5_jaccard_similarity.png", dpi=200)
    plt.close(fig)
else:
    print("    No Jaccard data -- skipping Fig 5")

# ── Fig 6: Fraction of consistent bridges across all phases ──────────
print("  Fig 6: Consistent bridge nodes")
fig, ax = plt.subplots(figsize=(10, 6))
consistent_data = []
for patient in FOUR_PHASE_PATIENTS:
    for band in ALL_BANDS:
        bd = bridge_data[patient][band]
        phases_avail = [p for p in ALL_PHASES if p in bd]
        if len(phases_avail) < 4:
            continue
        # Check all have same n_nodes
        n_nodes_list = [len(bd[p]) for p in phases_avail]
        if len(set(n_nodes_list)) != 1:
            continue
        n_nodes = n_nodes_list[0]
        # For each node: is it a top bridge in ALL phases?
        bridge_sets = [top_bridge_set(bd[p]) for p in phases_avail]
        consistent = bridge_sets[0]
        for bs in bridge_sets[1:]:
            consistent = consistent & bs
        frac = len(consistent) / n_nodes
        consistent_data.append({
            "patient": patient, "band": band,
            "n_consistent": len(consistent), "n_nodes": n_nodes,
            "fraction": frac,
        })

if consistent_data:
    df_cons = pd.DataFrame(consistent_data)
    # Grouped bar: patients on x, bands as groups
    x = np.arange(len(FOUR_PHASE_PATIENTS))
    width = 0.12
    for i, band in enumerate(ALL_BANDS):
        vals = []
        for patient in FOUR_PHASE_PATIENTS:
            sub = df_cons[(df_cons["patient"] == patient) & (df_cons["band"] == band)]
            vals.append(sub["fraction"].values[0] if len(sub) > 0 else 0)
        ax.bar(x + (i - 2.5) * width, vals, width, label=BRAIN_BAND_TEX_DICT[band],
               color=BAND_COLORS[i], alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(FOUR_PHASE_PATIENTS, fontsize=11)
    ax.set_ylabel("Fraction of Nodes", fontsize=13)
    ax.set_title("Fraction of Nodes that are Consistently Bridges\nacross ALL 4 Phases (top 20%)",
                 fontsize=13)
    ax.legend(fontsize=9, ncol=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.tight_layout()
fig.savefig(OUT_DIR / "fig6_consistent_bridges.png", dpi=200)
plt.close(fig)

# ── Summary statistics for FINDINGS.md ────────────────────────────────
print("\n" + "=" * 70)
print("SUMMARY STATISTICS")
print("=" * 70)

findings_lines = ["# Hierarchical Bridge Nodes Analysis - Findings\n"]
findings_lines.append(f"**Date:** Generated by bridge_nodes_analysis.py\n")
findings_lines.append(f"**FC method:** {FC_METHOD}\n")
findings_lines.append(f"**K range:** {K_RANGE.start} to {K_RANGE.stop - 1} ({len(K_RANGE) - 1} transitions)\n")
findings_lines.append(f"**Bridge threshold:** top {TOP_FRAC*100:.0f}%\n\n")

# Overall bridge score stats
findings_lines.append("## 1. Bridge Score Distribution\n")
overall_mean = df["bridge_score"].mean()
overall_std = df["bridge_score"].std()
overall_median = df["bridge_score"].median()
findings_lines.append(f"- Overall mean bridge score: {overall_mean:.4f} (std={overall_std:.4f}, median={overall_median:.4f})\n")

by_band = df.groupby("band")["bridge_score"].agg(["mean", "std"]).reindex(ALL_BANDS)
findings_lines.append("- By band:\n")
for band in ALL_BANDS:
    findings_lines.append(f"  - {band}: mean={by_band.loc[band, 'mean']:.4f}, std={by_band.loc[band, 'std']:.4f}\n")

by_phase = df.groupby("phase")["bridge_score"].agg(["mean", "std"]).reindex(ALL_PHASES)
findings_lines.append("- By phase:\n")
for phase in ALL_PHASES:
    if phase in by_phase.index and not np.isnan(by_phase.loc[phase, "mean"]):
        findings_lines.append(f"  - {phase}: mean={by_phase.loc[phase, 'mean']:.4f}, std={by_phase.loc[phase, 'std']:.4f}\n")

# Band with most bridge nodes
band_bridge_frac = {}
for band in ALL_BANDS:
    sub = df[df["band"] == band]
    threshold_val = np.percentile(sub["bridge_score"], 100 * (1 - TOP_FRAC))
    band_bridge_frac[band] = (sub["bridge_score"] >= threshold_val).mean()
most_bridges_band = max(band_bridge_frac, key=band_bridge_frac.get)
findings_lines.append(f"\n- Band with highest mean bridge score: **{most_bridges_band}** "
                      f"(mean={by_band.loc[most_bridges_band, 'mean']:.4f})\n")

# Cross-phase correlations
findings_lines.append("\n## 2. Cross-Phase Stability of Bridge Scores\n")
if len(df_xphase) > 0:
    mean_rho = df_xphase["spearman_rho"].mean()
    median_rho = df_xphase["spearman_rho"].median()
    findings_lines.append(f"- Mean cross-phase Spearman rho: {mean_rho:.4f} (median={median_rho:.4f})\n")

    # By condition pair type
    for pa, pb in phase_pairs:
        sub = df_xphase[(df_xphase["phase_1"] == pa) & (df_xphase["phase_2"] == pb)]
        if len(sub) > 0:
            findings_lines.append(f"  - {pa} vs {pb}: rho={sub['spearman_rho'].mean():.4f} "
                                  f"(n={len(sub)} patient-band combos)\n")

    # Interpretation
    if mean_rho > 0.5:
        findings_lines.append("\n**Interpretation:** Bridge scores are HIGHLY correlated across phases, "
                              "suggesting bridge nodes are structurally determined.\n")
    elif mean_rho > 0.2:
        findings_lines.append("\n**Interpretation:** Bridge scores show MODERATE cross-phase correlation, "
                              "suggesting partial structural determination with state-dependent modulation.\n")
    else:
        findings_lines.append("\n**Interpretation:** Bridge scores show WEAK cross-phase correlation, "
                              "suggesting bridge identity is primarily state-dependent.\n")

# Cross-patient correlations
findings_lines.append("\n## 3. Cross-Patient Correlation\n")
if len(df_xpat) > 0:
    mean_xpat = df_xpat["spearman_rho"].mean()
    median_xpat = df_xpat["spearman_rho"].median()
    findings_lines.append(f"- Mean cross-patient Spearman rho: {mean_xpat:.4f} (median={median_xpat:.4f})\n")
    findings_lines.append(f"- Number of valid patient pairs (matching node count): {len(df_xpat)}\n")

    if mean_xpat > 0.3:
        findings_lines.append("\n**Interpretation:** Same nodes tend to be bridges across patients, "
                              "suggesting structural origin (anatomical rather than functional).\n")
    elif mean_xpat > 0.1:
        findings_lines.append("\n**Interpretation:** Weak cross-patient similarity suggests "
                              "bridge identity is mostly patient-specific.\n")
    else:
        findings_lines.append("\n**Interpretation:** Bridge identity is NOT shared across patients.\n")
else:
    findings_lines.append("- No valid cross-patient pairs (node counts differ). "
                          "Bridge identity is patient-specific by construction (different electrode placements).\n")

# Jaccard analysis
findings_lines.append("\n## 4. Jaccard Overlap of Top Bridge Sets\n")
if len(df_jaccard) > 0:
    for pt in ["within-rest", "within-task", "cross-condition"]:
        vals = df_jaccard[df_jaccard["pair_type"] == pt]["jaccard"]
        if len(vals) > 0:
            findings_lines.append(f"- {pt}: mean Jaccard = {vals.mean():.4f} (std={vals.std():.4f}, n={len(vals)})\n")

    within_vals = df_jaccard[df_jaccard["pair_type"].str.startswith("within")]["jaccard"]
    cross_vals = df_jaccard[df_jaccard["pair_type"] == "cross-condition"]["jaccard"]
    if len(within_vals) > 0 and len(cross_vals) > 0:
        if within_vals.mean() > cross_vals.mean() + 0.05:
            findings_lines.append("\n**Interpretation:** Within-condition bridge overlap is HIGHER than cross-condition, "
                                  "suggesting task-specific reorganization at the node level.\n")
        else:
            findings_lines.append("\n**Interpretation:** Bridge overlap is SIMILAR within and across conditions, "
                                  "suggesting stable bridge identity regardless of condition.\n")

# Consistent bridges
findings_lines.append("\n## 5. Consistently Bridge Nodes (across all 4 phases)\n")
if consistent_data:
    df_cons = pd.DataFrame(consistent_data)
    mean_frac = df_cons["fraction"].mean()
    findings_lines.append(f"- Mean fraction of consistent bridges: {mean_frac:.4f} "
                          f"({mean_frac*100:.1f}% of nodes)\n")
    for band in ALL_BANDS:
        sub = df_cons[df_cons["band"] == band]
        if len(sub) > 0:
            findings_lines.append(f"  - {band}: {sub['fraction'].mean():.4f}\n")

findings_lines.append("\n## 6. Key Takeaways\n")
findings_lines.append("- Bridge scores quantify how frequently a node changes community "
                      "membership as hierarchical resolution varies.\n")
findings_lines.append("- Cross-phase and cross-patient correlations reveal whether "
                      "bridge identity is a structural property or state-dependent.\n")
findings_lines.append("- Jaccard overlap distinguishes task-specific from stable bridge nodes.\n")

findings_text = "".join(findings_lines)
(OUT_DIR / "FINDINGS.md").write_text(findings_text)
print(findings_text)

print("\n" + "=" * 70)
print(f"ALL OUTPUTS SAVED TO: {OUT_DIR}")
print("=" * 70)
