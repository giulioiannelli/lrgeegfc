#!/usr/bin/env python3
"""Dendrogram Merge Sequence Analysis.

Compares the ORDER in which sub-communities merge across phases to reveal
shared hierarchical structure independent of metric distances.
"""

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
from pathlib import Path
from itertools import combinations
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr, pearsonr
from scipy.cluster.hierarchy import fcluster
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE
from lrg_eegfc.config.paths import FIGURES_ROOT

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
FC_METHOD = "msc"
ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]
FOUR_PHASE_PATIENTS = PATIENTS_4PHASE
BANDS = list(BRAIN_BANDS.keys())
PHASES = list(PHASE_LABELS)  # rest_pre, task_learn, task_test, rest_post

WITHIN_PAIRS = [("rest_pre", "rest_post"), ("task_learn", "task_test")]
CROSS_PAIRS = [
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"),
    ("task_learn", "rest_post"), ("task_test", "rest_post"),
]

OUTPUT_DIR = FIGURES_ROOT / "multiscale_investigation" / "merge_sequence"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Plot styling
plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "figure.dpi": 150,
    "savefig.dpi": 200,
})


# ---------------------------------------------------------------------------
# Core algorithm: compute merge ranks for all node pairs
# ---------------------------------------------------------------------------
def compute_merge_ranks(linkage_matrix):
    """For each pair of original nodes, find the merge step (rank) at which
    they first end up in the same cluster.

    Parameters
    ----------
    linkage_matrix : np.ndarray, shape (n-1, 4)
        Standard scipy linkage matrix.

    Returns
    -------
    merge_rank_matrix : np.ndarray, shape (n, n)
        Symmetric matrix where entry (i,j) = merge step at which nodes
        i and j first share a cluster. Diagonal = 0.
    """
    n = linkage_matrix.shape[0] + 1  # number of original nodes

    # Track which original nodes belong to each cluster
    # Clusters 0..n-1 are original nodes; cluster n+i is formed at step i
    cluster_members = {i: {i} for i in range(n)}

    merge_rank = np.zeros((n, n), dtype=np.float64)

    for step, row in enumerate(linkage_matrix):
        c1, c2 = int(row[0]), int(row[1])
        new_cluster_id = n + step

        members_1 = cluster_members[c1]
        members_2 = cluster_members[c2]

        # All pairs across the two merging clusters get this merge step
        for a in members_1:
            for b in members_2:
                merge_rank[a, b] = step
                merge_rank[b, a] = step

        # Create new cluster as union
        cluster_members[new_cluster_id] = members_1 | members_2

        # Clean up (optional, saves memory)
        del cluster_members[c1]
        del cluster_members[c2]

    return merge_rank


def upper_triangle_vector(matrix):
    """Extract upper triangle (no diagonal) as a flat vector."""
    n = matrix.shape[0]
    idx = np.triu_indices(n, k=1)
    return matrix[idx]


# ---------------------------------------------------------------------------
# Load all data
# ---------------------------------------------------------------------------
print("Loading LRG results...")
results = {}
skipped = []

for patient in ALL_PATIENTS:
    for band in BANDS:
        for phase in PHASES:
            # Pat_06 and Pat_07 only have task_learn, task_test
            if patient in ("Pat_06", "Pat_07") and phase in ("rest_pre", "rest_post"):
                continue
            res = load_lrg_result(patient, phase, band, FC_METHOD)
            if res is not None:
                results[(patient, band, phase)] = res
            else:
                skipped.append((patient, band, phase))

print(f"Loaded {len(results)} LRG results, skipped {len(skipped)}")
if skipped:
    print(f"  Skipped: {skipped[:5]}{'...' if len(skipped) > 5 else ''}")


# ---------------------------------------------------------------------------
# Compute merge ranks and ultrametric vectors for all results
# ---------------------------------------------------------------------------
print("Computing merge ranks...")
merge_rank_vectors = {}
ultrametric_vectors = {}

for key, res in results.items():
    Z = res.linkage_matrix
    n_nodes = res.n_nodes

    # Merge rank matrix
    mr = compute_merge_ranks(Z)
    merge_rank_vectors[key] = upper_triangle_vector(mr)

    # Ultrametric (cophenetic) distance vector
    ultra_sq = squareform(res.ultrametric_matrix)
    ultrametric_vectors[key] = upper_triangle_vector(ultra_sq)

print(f"Computed merge ranks for {len(merge_rank_vectors)} conditions")


# ---------------------------------------------------------------------------
# Compute pairwise correlations between phases
# ---------------------------------------------------------------------------
print("Computing pairwise phase correlations...")

records = []

for patient in FOUR_PHASE_PATIENTS:
    for band in BANDS:
        for p1, p2 in list(combinations(PHASES, 2)):
            key1 = (patient, band, p1)
            key2 = (patient, band, p2)

            if key1 not in merge_rank_vectors or key2 not in merge_rank_vectors:
                continue

            # Check compatible sizes (same giant component)
            v1_mr = merge_rank_vectors[key1]
            v2_mr = merge_rank_vectors[key2]
            v1_um = ultrametric_vectors[key1]
            v2_um = ultrametric_vectors[key2]

            if len(v1_mr) != len(v2_mr):
                # Different giant component sizes - skip
                continue

            # Merge-rank Spearman correlation (topological - order only)
            rho_mr, p_mr = spearmanr(v1_mr, v2_mr)

            # Ultrametric-distance Pearson correlation (metric - sensitive to actual distances)
            rho_um, p_um = pearsonr(v1_um, v2_um)

            pair_label = f"{p1}-{p2}"
            pair_type = "within" if (p1, p2) in WITHIN_PAIRS else "cross"

            records.append({
                "patient": patient,
                "band": band,
                "phase_1": p1,
                "phase_2": p2,
                "pair_label": pair_label,
                "pair_type": pair_type,
                "merge_rank_rho": rho_mr,
                "merge_rank_p": p_mr,
                "ultrametric_rho": rho_um,
                "ultrametric_p": p_um,
                "n_pairs": len(v1_mr),
            })

df = pd.DataFrame(records)
print(f"Computed {len(df)} pairwise comparisons")
print(f"  Within-condition: {(df['pair_type'] == 'within').sum()}")
print(f"  Cross-condition:  {(df['pair_type'] == 'cross').sum()}")

# Save CSV
csv_path = OUTPUT_DIR / "merge_sequence_correlations.csv"
df.to_csv(csv_path, index=False)
print(f"Saved CSV to {csv_path}")

# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("SUMMARY STATISTICS")
print("=" * 70)

for ptype in ["within", "cross"]:
    sub = df[df["pair_type"] == ptype]
    print(f"\n{ptype.upper()} condition (n={len(sub)}):")
    print(f"  Merge-rank rho:   {sub['merge_rank_rho'].mean():.4f} +/- {sub['merge_rank_rho'].std():.4f}")
    print(f"  Ultrametric rho:  {sub['ultrametric_rho'].mean():.4f} +/- {sub['ultrametric_rho'].std():.4f}")

print(f"\nMerge-rank within vs cross delta: "
      f"{df[df['pair_type']=='within']['merge_rank_rho'].mean() - df[df['pair_type']=='cross']['merge_rank_rho'].mean():.4f}")

# Per-band summary
print("\nPer-band merge-rank correlation (mean +/- std):")
for band in BANDS:
    sub = df[df["band"] == band]
    within = sub[sub["pair_type"] == "within"]["merge_rank_rho"]
    cross = sub[sub["pair_type"] == "cross"]["merge_rank_rho"]
    tex = BRAIN_BAND_TEX_DICT[band]
    print(f"  {band:12s}  within={within.mean():.4f}+/-{within.std():.4f}  "
          f"cross={cross.mean():.4f}+/-{cross.std():.4f}  "
          f"delta={within.mean()-cross.mean():.4f}")


# ===========================================================================
# FIGURE 1: Phase x Phase correlation matrices, one per band
# ===========================================================================
print("\nGenerating Figure 1: Phase correlation matrices...")

fig, axes = plt.subplots(2, 3, figsize=(14, 9))
axes = axes.ravel()

for idx, band in enumerate(BANDS):
    ax = axes[idx]
    n_phases = len(PHASES)
    corr_mat = np.full((n_phases, n_phases), np.nan)
    np.fill_diagonal(corr_mat, 1.0)

    for i, p1 in enumerate(PHASES):
        for j, p2 in enumerate(PHASES):
            if i >= j:
                continue
            sub = df[(df["band"] == band) &
                     (df["phase_1"] == p1) & (df["phase_2"] == p2)]
            if len(sub) > 0:
                val = sub["merge_rank_rho"].mean()
                corr_mat[i, j] = val
                corr_mat[j, i] = val

    im = ax.imshow(corr_mat, cmap="RdYlBu_r", vmin=0.0, vmax=1.0,
                   aspect="equal")

    # Annotate cells
    for i in range(n_phases):
        for j in range(n_phases):
            if not np.isnan(corr_mat[i, j]):
                val = corr_mat[i, j]
                color = "white" if val > 0.7 or val < 0.3 else "black"
                ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                        fontsize=9, color=color, fontweight="bold")

    phase_short = ["rest_pre", "tLearn", "tTest", "rest_post"]
    ax.set_xticks(range(n_phases))
    ax.set_xticklabels(phase_short, rotation=45, ha="right")
    ax.set_yticks(range(n_phases))
    ax.set_yticklabels(phase_short)
    ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]} ({band})", fontsize=12)

fig.suptitle("Merge-Rank Spearman Correlation Between Phases\n"
             "(averaged across 4-phase patients, MSC method)",
             fontsize=14, fontweight="bold", y=1.02)
fig.colorbar(im, ax=axes.tolist(), shrink=0.6, label="Spearman rho")
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig1_phase_corr_matrices.png")
plt.close(fig)
print("  Saved fig1_phase_corr_matrices.png")


# ===========================================================================
# FIGURE 2: Boxplot within vs cross condition
# ===========================================================================
print("Generating Figure 2: Within vs cross boxplot...")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: all bands pooled
ax = axes[0]
within_vals = df[df["pair_type"] == "within"]["merge_rank_rho"].values
cross_vals = df[df["pair_type"] == "cross"]["merge_rank_rho"].values

bp = ax.boxplot([within_vals, cross_vals], labels=["Within\n(rest-rest, task-task)",
                                                     "Cross\n(rest-task)"],
                patch_artist=True, widths=0.5)
bp["boxes"][0].set_facecolor("#4ECDC4")
bp["boxes"][1].set_facecolor("#FF6B6B")

# Overlay individual points
for i, (vals, x_pos) in enumerate(zip([within_vals, cross_vals], [1, 2])):
    jitter = np.random.default_rng(42).uniform(-0.12, 0.12, size=len(vals))
    ax.scatter(np.full(len(vals), x_pos) + jitter, vals,
               alpha=0.4, s=20, c="gray", zorder=3)

ax.set_ylabel("Merge-Rank Spearman Correlation")
ax.set_title("A) All bands pooled", fontweight="bold")
ax.set_ylim(-0.05, 1.05)

# Panel B: per band
ax = axes[1]
band_positions = []
band_labels = []
colors_within = "#4ECDC4"
colors_cross = "#FF6B6B"

x = 0
for band in BANDS:
    sub = df[df["band"] == band]
    w = sub[sub["pair_type"] == "within"]["merge_rank_rho"].values
    c = sub[sub["pair_type"] == "cross"]["merge_rank_rho"].values

    bp_w = ax.boxplot([w], positions=[x], widths=0.35, patch_artist=True,
                      manage_ticks=False)
    bp_w["boxes"][0].set_facecolor(colors_within)
    bp_c = ax.boxplot([c], positions=[x + 0.4], widths=0.35, patch_artist=True,
                      manage_ticks=False)
    bp_c["boxes"][0].set_facecolor(colors_cross)

    band_positions.append(x + 0.2)
    band_labels.append(BRAIN_BAND_TEX_DICT[band])
    x += 1.2

ax.set_xticks(band_positions)
ax.set_xticklabels(band_labels)
ax.set_ylabel("Merge-Rank Spearman Correlation")
ax.set_title("B) Per band (teal=within, red=cross)", fontweight="bold")
ax.set_ylim(-0.05, 1.05)

fig.suptitle("Within- vs Cross-Condition Merge Sequence Similarity",
             fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig2_within_vs_cross_boxplot.png")
plt.close(fig)
print("  Saved fig2_within_vs_cross_boxplot.png")


# ===========================================================================
# FIGURE 3: Merge-RANK vs Ultrametric-distance correlation scatter
# ===========================================================================
print("Generating Figure 3: Rank vs ultrametric scatter...")

fig, ax = plt.subplots(figsize=(7, 7))

within_mask = df["pair_type"] == "within"
cross_mask = df["pair_type"] == "cross"

ax.scatter(df.loc[within_mask, "ultrametric_rho"],
           df.loc[within_mask, "merge_rank_rho"],
           c="#4ECDC4", s=50, alpha=0.7, edgecolors="black", linewidths=0.5,
           label="Within-condition", zorder=3)
ax.scatter(df.loc[cross_mask, "ultrametric_rho"],
           df.loc[cross_mask, "merge_rank_rho"],
           c="#FF6B6B", s=50, alpha=0.7, edgecolors="black", linewidths=0.5,
           label="Cross-condition", zorder=3)

# Identity line
lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1])]
ax.plot(lims, lims, "k--", alpha=0.5, linewidth=1, label="y = x")
ax.set_xlim(lims)
ax.set_ylim(lims)

ax.set_xlabel("Ultrametric Distance Pearson Correlation")
ax.set_ylabel("Merge-Rank Spearman Correlation")
ax.set_title("Topological (Spearman rank) vs Metric (Pearson distance) Preservation\n"
             "Points above diagonal: merge ORDER more preserved than distances",
             fontweight="bold", fontsize=11)
ax.legend(loc="lower right")
ax.set_aspect("equal")

# Count above/below diagonal
above = (df["merge_rank_rho"] > df["ultrametric_rho"]).sum()
below = (df["merge_rank_rho"] < df["ultrametric_rho"]).sum()
equal = (df["merge_rank_rho"] == df["ultrametric_rho"]).sum()
ax.text(0.05, 0.95, f"Above diagonal: {above}/{len(df)}\nBelow diagonal: {below}/{len(df)}",
        transform=ax.transAxes, fontsize=10, va="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.8))

fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig3_rank_vs_ultrametric_scatter.png")
plt.close(fig)
print("  Saved fig3_rank_vs_ultrametric_scatter.png")


# ===========================================================================
# FIGURE 4: Per-band effect size (within - cross) for merge-rank
# ===========================================================================
print("Generating Figure 4: Per-band effect size...")

fig, ax = plt.subplots(figsize=(9, 5))

band_deltas_mr = []
band_deltas_um = []
band_within_mr = []
band_cross_mr = []

for band in BANDS:
    sub = df[df["band"] == band]
    w_mr = sub[sub["pair_type"] == "within"]["merge_rank_rho"].mean()
    c_mr = sub[sub["pair_type"] == "cross"]["merge_rank_rho"].mean()
    w_um = sub[sub["pair_type"] == "within"]["ultrametric_rho"].mean()
    c_um = sub[sub["pair_type"] == "cross"]["ultrametric_rho"].mean()
    band_deltas_mr.append(w_mr - c_mr)
    band_deltas_um.append(w_um - c_um)
    band_within_mr.append(w_mr)
    band_cross_mr.append(c_mr)

x = np.arange(len(BANDS))
width = 0.35

bars1 = ax.bar(x - width / 2, band_deltas_mr, width, label="Merge-Rank delta",
               color="#4ECDC4", edgecolor="black", linewidth=0.5)
bars2 = ax.bar(x + width / 2, band_deltas_um, width, label="Ultrametric delta",
               color="#FF6B6B", edgecolor="black", linewidth=0.5)

ax.set_xticks(x)
ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS])
ax.set_ylabel("Within - Cross Correlation (effect size)")
ax.set_title("Per-Band Effect Size: Within vs Cross Condition\n"
             "(positive = within-condition pairs more similar)",
             fontweight="bold")
ax.axhline(0, color="black", linewidth=0.8)
ax.legend()

# Annotate values
for bar in bars1:
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, h + 0.002, f"{h:.3f}",
            ha="center", va="bottom", fontsize=8)
for bar in bars2:
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, h + 0.002, f"{h:.3f}",
            ha="center", va="bottom", fontsize=8)

fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig4_per_band_effect_size.png")
plt.close(fig)
print("  Saved fig4_per_band_effect_size.png")


# ===========================================================================
# FIGURE 5: Per-patient profiles
# ===========================================================================
print("Generating Figure 5: Per-patient profiles...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.ravel()

for idx, patient in enumerate(FOUR_PHASE_PATIENTS):
    ax = axes[idx]
    sub = df[df["patient"] == patient]

    # Get unique phase pairs
    pair_labels = []
    within_mr = []
    cross_mr = []

    for band in BANDS:
        band_sub = sub[sub["band"] == band]
        w = band_sub[band_sub["pair_type"] == "within"]["merge_rank_rho"].values
        c = band_sub[band_sub["pair_type"] == "cross"]["merge_rank_rho"].values
        within_mr.append(w.mean() if len(w) > 0 else np.nan)
        cross_mr.append(c.mean() if len(c) > 0 else np.nan)

    x = np.arange(len(BANDS))
    ax.plot(x, within_mr, "o-", color="#4ECDC4", linewidth=2, markersize=8,
            label="Within", zorder=3)
    ax.plot(x, cross_mr, "s-", color="#FF6B6B", linewidth=2, markersize=8,
            label="Cross", zorder=3)
    ax.fill_between(x, within_mr, cross_mr, alpha=0.15, color="gray")

    ax.set_xticks(x)
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS])
    ax.set_ylabel("Merge-Rank Spearman rho")
    ax.set_title(f"{patient}", fontweight="bold")
    ax.legend(loc="lower left")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(alpha=0.3)

fig.suptitle("Per-Patient Merge Sequence Preservation Across Bands\n"
             "(Within vs Cross condition, MSC method)",
             fontsize=14, fontweight="bold")
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig5_per_patient_profiles.png")
plt.close(fig)
print("  Saved fig5_per_patient_profiles.png")


# ===========================================================================
# FINDINGS.md
# ===========================================================================
print("\nWriting FINDINGS.md...")

# Compute key stats for findings
within_mr_mean = df[df["pair_type"] == "within"]["merge_rank_rho"].mean()
within_mr_std = df[df["pair_type"] == "within"]["merge_rank_rho"].std()
cross_mr_mean = df[df["pair_type"] == "cross"]["merge_rank_rho"].mean()
cross_mr_std = df[df["pair_type"] == "cross"]["merge_rank_rho"].std()

within_um_mean = df[df["pair_type"] == "within"]["ultrametric_rho"].mean()
cross_um_mean = df[df["pair_type"] == "cross"]["ultrametric_rho"].mean()

above_diag = (df["merge_rank_rho"] > df["ultrametric_rho"]).sum()
below_diag = (df["merge_rank_rho"] < df["ultrametric_rho"]).sum()
total = len(df)

# Best and worst bands
band_effects = {}
for band in BANDS:
    sub = df[df["band"] == band]
    w = sub[sub["pair_type"] == "within"]["merge_rank_rho"].mean()
    c = sub[sub["pair_type"] == "cross"]["merge_rank_rho"].mean()
    band_effects[band] = w - c

best_band = max(band_effects, key=band_effects.get)
worst_band = min(band_effects, key=band_effects.get)

# Per-pair breakdown
pair_summary = df.groupby("pair_label").agg(
    merge_rank_rho_mean=("merge_rank_rho", "mean"),
    merge_rank_rho_std=("merge_rank_rho", "std"),
    ultrametric_rho_mean=("ultrametric_rho", "mean"),
    ultrametric_rho_std=("ultrametric_rho", "std"),
    count=("merge_rank_rho", "count"),
).round(4)

findings = f"""# Dendrogram Merge Sequence Analysis - Findings

## Method
For each LRG result, the linkage matrix encodes the order in which sub-communities
merge as you traverse the hierarchy from leaves to root. For each pair of original
nodes (i, j), we compute the "merge rank" = the step at which they first share a
cluster. We then compare these merge-rank vectors across phases using Spearman
correlation.

- **FC method**: MSC (magnitude-squared coherence)
- **Patients**: {', '.join(FOUR_PHASE_PATIENTS)} (4-phase)
- **Bands**: {', '.join(BANDS)}
- **Total comparisons**: {total}

## Key Results

### 1. Within vs Cross-Condition Merge Sequence Correlation
- **Within-condition** (rest-rest, task-task): rho = {within_mr_mean:.4f} +/- {within_mr_std:.4f}
- **Cross-condition** (rest-task): rho = {cross_mr_mean:.4f} +/- {cross_mr_std:.4f}
- **Effect size** (within - cross): {within_mr_mean - cross_mr_mean:.4f}

### 2. Topological vs Metric Preservation
- Merge-RANK correlation (topological): within={within_mr_mean:.4f}, cross={cross_mr_mean:.4f}
- Ultrametric-distance Pearson correlation (metric): within={within_um_mean:.4f}, cross={cross_um_mean:.4f}
- Points above diagonal (rank > ultrametric): {above_diag}/{total} ({100*above_diag/total:.1f}%)
- Points below diagonal (rank < ultrametric): {below_diag}/{total} ({100*below_diag/total:.1f}%)

### 3. Per-Band Effect Sizes (within - cross delta)
"""

for band in BANDS:
    findings += f"- **{band}**: {band_effects[band]:.4f}\n"

findings += f"""
- **Largest effect**: {best_band} ({band_effects[best_band]:.4f})
- **Smallest effect**: {worst_band} ({band_effects[worst_band]:.4f})

### 4. Phase-Pair Breakdown
| Pair | Merge-Rank rho | Ultrametric rho | N |
|------|---------------|-----------------|---|
"""

for pair_label, row in pair_summary.iterrows():
    findings += (f"| {pair_label} | {row['merge_rank_rho_mean']:.4f} +/- "
                 f"{row['merge_rank_rho_std']:.4f} | {row['ultrametric_rho_mean']:.4f} +/- "
                 f"{row['ultrametric_rho_std']:.4f} | {int(row['count'])} |\n")

findings += f"""
### 5. Per-Patient Consistency
"""

for patient in FOUR_PHASE_PATIENTS:
    sub = df[df["patient"] == patient]
    w = sub[sub["pair_type"] == "within"]["merge_rank_rho"].mean()
    c = sub[sub["pair_type"] == "cross"]["merge_rank_rho"].mean()
    findings += f"- **{patient}**: within={w:.4f}, cross={c:.4f}, delta={w-c:.4f}\n"

findings += """
## Figures
1. `fig1_phase_corr_matrices.png` - Phase x phase merge-rank correlation, one subplot per band
2. `fig2_within_vs_cross_boxplot.png` - Boxplot of within vs cross condition correlations
3. `fig3_rank_vs_ultrametric_scatter.png` - Merge-rank vs ultrametric correlation scatter
4. `fig4_per_band_effect_size.png` - Per-band effect sizes for merge-rank and ultrametric
5. `fig5_per_patient_profiles.png` - Per-patient profiles across bands

## Interpretation
- High merge-rank correlation means the same node pairs tend to merge early/late across
  both phases, indicating shared hierarchical community structure.
- If merge-rank correlation exceeds ultrametric correlation, the TOPOLOGY (merge order) is
  more preserved than the exact distances, suggesting robust hierarchical scaffolding
  that is metrically modulated by cognitive state.
"""

findings_path = OUTPUT_DIR / "FINDINGS.md"
with open(findings_path, "w") as f:
    f.write(findings)
print(f"Saved FINDINGS.md to {findings_path}")

print("\nDone!")
