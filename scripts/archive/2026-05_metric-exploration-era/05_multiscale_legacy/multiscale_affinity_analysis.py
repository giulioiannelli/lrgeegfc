#!/usr/bin/env python3
"""Multiscale Co-Classification Probability Matrices.

For each pair of nodes (i,j), compute P(same cluster) averaged across ALL
dendrogram cut levels. This gives a "multiscale affinity matrix" that encodes
the full tree structure into a single N×N matrix.

Compare these affinity matrices across phases and contrast with raw FC-based
comparison.
"""

import csv
from pathlib import Path
from itertools import combinations

import numpy as np
from scipy.cluster.hierarchy import fcluster
from scipy.stats import spearmanr, pearsonr
from scipy.spatial.distance import squareform

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.gridspec as gridspec

from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE
from lrg_eegfc.config.paths import FIGURES_ROOT, LRG_CACHE

OUTPUT_DIR = FIGURES_ROOT / "multiscale_investigation" / "multiscale_affinity"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]
FOUR_PHASE_PATIENTS = PATIENTS_4PHASE
ALL_PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
BANDS = list(BRAIN_BANDS.keys())
FC_METHOD = "msc"
N_THRESHOLDS = 100  # number of evenly-spaced cut thresholds

WITHIN_PAIRS = [("rest_pre", "rest_post"), ("task_learn", "task_test")]
CROSS_PAIRS = [
    ("rest_pre", "task_learn"),
    ("rest_pre", "task_test"),
    ("task_learn", "rest_post"),
    ("task_test", "rest_post"),
]

# ---------------------------------------------------------------------------
# Core algorithm
# ---------------------------------------------------------------------------

def compute_affinity_matrix(linkage_matrix: np.ndarray, n_nodes: int,
                            n_thresholds: int = N_THRESHOLDS) -> np.ndarray:
    """Compute multiscale co-classification probability matrix.

    Sweep through evenly-spaced distance thresholds across the full linkage
    range. At each threshold, compute the flat partition and the co-classification
    matrix C_k[i,j] = 1 if i,j are in the same cluster.  Average all C_k.

    Returns
    -------
    affinity : np.ndarray, shape (n_nodes, n_nodes)
        Symmetric matrix with values in [0, 1].
    """
    distances = linkage_matrix[:, 2]
    d_min, d_max = distances.min(), distances.max()
    # Thresholds from just above d_min to d_max (inclusive)
    thresholds = np.linspace(d_min + 1e-12, d_max, n_thresholds)

    # Accumulate co-classification counts
    affinity = np.zeros((n_nodes, n_nodes), dtype=np.float64)

    for t in thresholds:
        labels = fcluster(linkage_matrix, t=t, criterion="distance")
        # Build co-classification: same label → 1
        # Vectorized: outer comparison
        co = (labels[:, None] == labels[None, :]).astype(np.float64)
        affinity += co

    affinity /= n_thresholds
    return affinity


def upper_tri_vec(mat: np.ndarray) -> np.ndarray:
    """Extract upper triangle (excluding diagonal) as a flat vector."""
    idx = np.triu_indices_from(mat, k=1)
    return mat[idx]


# ---------------------------------------------------------------------------
# Load / compute all affinity matrices and FC matrices
# ---------------------------------------------------------------------------

print("=" * 70)
print("Computing multiscale affinity matrices")
print("=" * 70)

# Store results: affinity_matrices[patient][band][phase] = ndarray
affinity_matrices = {}
fc_matrices = {}

# Track available phases per patient
patient_phases = {}

for pat in ALL_PATIENTS:
    affinity_matrices[pat] = {}
    fc_matrices[pat] = {}

    # Determine available phases
    avail = []
    for phase in ALL_PHASES:
        lrg_path = LRG_CACHE / pat / f"alpha_{phase}_lrg_{FC_METHOD}.npz"
        if lrg_path.exists():
            avail.append(phase)
    patient_phases[pat] = avail
    print(f"\n{pat}: phases = {avail}")

    for band in BANDS:
        affinity_matrices[pat][band] = {}
        fc_matrices[pat][band] = {}

        for phase in avail:
            # Load LRG result
            lrg_path = LRG_CACHE / pat / f"{band}_{phase}_lrg_{FC_METHOD}.npz"
            if not lrg_path.exists():
                print(f"  SKIP {pat}/{band}/{phase} - no LRG cache")
                continue

            data = np.load(str(lrg_path), allow_pickle=True)
            linkage_matrix = data["linkage_matrix"]
            n_nodes = int(data["n_nodes"])

            aff = compute_affinity_matrix(linkage_matrix, n_nodes)
            affinity_matrices[pat][band][phase] = aff

            # Load raw MSC FC matrix
            msc = load_msc_matrix(pat, phase, band)
            if msc is not None:
                if hasattr(msc, "adjacency_matrix"):
                    fc_matrices[pat][band][phase] = msc.adjacency_matrix
                else:
                    fc_matrices[pat][band][phase] = np.array(msc)

        print(f"  {band}: {len(affinity_matrices[pat][band])} phases done")


# ---------------------------------------------------------------------------
# Pairwise comparisons (4-phase patients only)
# ---------------------------------------------------------------------------

print("\n" + "=" * 70)
print("Computing pairwise correlations")
print("=" * 70)

records = []  # list of dicts for CSV

for pat in FOUR_PHASE_PATIENTS:
    for band in BANDS:
        aff_dict = affinity_matrices[pat][band]
        fc_dict = fc_matrices[pat][band]

        phases_avail = sorted(set(aff_dict.keys()) & set(fc_dict.keys()))
        if len(phases_avail) < 2:
            continue

        for p1, p2 in combinations(ALL_PHASES, 2):
            if p1 not in phases_avail or p2 not in phases_avail:
                continue

            pair_label = f"{p1}-{p2}"
            pair_tuple = (p1, p2)

            if pair_tuple in WITHIN_PAIRS or (p2, p1) in WITHIN_PAIRS:
                condition = "within"
            else:
                condition = "cross"

            # Affinity correlation
            v1 = upper_tri_vec(aff_dict[p1])
            v2 = upper_tri_vec(aff_dict[p2])
            aff_pearson, _ = pearsonr(v1, v2)
            aff_spearman, _ = spearmanr(v1, v2)

            # FC correlation
            f1 = upper_tri_vec(fc_dict[p1])
            f2 = upper_tri_vec(fc_dict[p2])
            fc_pearson, _ = pearsonr(f1, f2)
            fc_spearman, _ = spearmanr(f1, f2)

            # Frobenius distance (normalized)
            n = aff_dict[p1].shape[0]
            n_pairs = n * (n - 1) / 2
            aff_frob = np.linalg.norm(v1 - v2) / np.sqrt(n_pairs)
            fc_frob = np.linalg.norm(f1 - f2) / np.sqrt(n_pairs)

            records.append({
                "patient": pat,
                "band": band,
                "phase_pair": pair_label,
                "condition": condition,
                "aff_pearson": aff_pearson,
                "aff_spearman": aff_spearman,
                "fc_pearson": fc_pearson,
                "fc_spearman": fc_spearman,
                "aff_frob": aff_frob,
                "fc_frob": fc_frob,
            })

print(f"  Total pairwise records: {len(records)}")


# ---------------------------------------------------------------------------
# Save CSV
# ---------------------------------------------------------------------------

csv_path = OUTPUT_DIR / "pairwise_results.csv"
fieldnames = list(records[0].keys())
with open(csv_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(records)
print(f"  CSV saved to {csv_path}")


# ===========================================================================
# FIGURE 1: Pat_02 alpha affinity matrices across 4 phases
# ===========================================================================

print("\n" + "=" * 70)
print("Figure 1: Pat_02 alpha affinity matrices")
print("=" * 70)

fig, axes = plt.subplots(1, 4, figsize=(20, 5))
for i, phase in enumerate(ALL_PHASES):
    aff = affinity_matrices["Pat_02"]["alpha"][phase]
    im = axes[i].imshow(aff, cmap="viridis", vmin=0, vmax=1, aspect="equal")
    axes[i].set_title(phase, fontsize=14, fontweight="bold")
    axes[i].set_xlabel("Node")
    if i == 0:
        axes[i].set_ylabel("Node")

fig.suptitle("Multiscale Co-Classification Affinity — Pat_02, alpha band",
             fontsize=16, fontweight="bold", y=1.02)
cbar = fig.colorbar(im, ax=axes, shrink=0.8, label="P(same cluster)")
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig1_affinity_heatmaps_Pat02_alpha.png",
            dpi=200, bbox_inches="tight")
plt.close(fig)
print("  Saved fig1")


# ===========================================================================
# FIGURE 2: Scatter — affinity-based vs FC-based correlation, colored by condition
# ===========================================================================

print("\nFigure 2: Affinity vs FC correlation scatter")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

within_aff_p = [r["aff_pearson"] for r in records if r["condition"] == "within"]
within_fc_p = [r["fc_pearson"] for r in records if r["condition"] == "within"]
cross_aff_p = [r["aff_pearson"] for r in records if r["condition"] == "cross"]
cross_fc_p = [r["fc_pearson"] for r in records if r["condition"] == "cross"]

ax = axes[0]
ax.scatter(cross_fc_p, cross_aff_p, c="tab:red", alpha=0.5, s=40, label="Cross-condition", zorder=2)
ax.scatter(within_fc_p, within_aff_p, c="tab:blue", alpha=0.7, s=50, label="Within-condition", zorder=3)
ax.set_xlabel("FC-based Pearson r", fontsize=12)
ax.set_ylabel("Affinity-based Pearson r", fontsize=12)
ax.set_title("Pearson correlation", fontsize=13, fontweight="bold")
ax.legend(fontsize=10)
ax.axline((0, 0), slope=1, color="gray", ls="--", alpha=0.5, zorder=1)
ax.set_aspect("equal")
lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]), max(ax.get_xlim()[1], ax.get_ylim()[1])]
ax.set_xlim(lims)
ax.set_ylim(lims)

within_aff_s = [r["aff_spearman"] for r in records if r["condition"] == "within"]
within_fc_s = [r["fc_spearman"] for r in records if r["condition"] == "within"]
cross_aff_s = [r["aff_spearman"] for r in records if r["condition"] == "cross"]
cross_fc_s = [r["fc_spearman"] for r in records if r["condition"] == "cross"]

ax = axes[1]
ax.scatter(cross_fc_s, cross_aff_s, c="tab:red", alpha=0.5, s=40, label="Cross-condition", zorder=2)
ax.scatter(within_fc_s, within_aff_s, c="tab:blue", alpha=0.7, s=50, label="Within-condition", zorder=3)
ax.set_xlabel("FC-based Spearman rho", fontsize=12)
ax.set_ylabel("Affinity-based Spearman rho", fontsize=12)
ax.set_title("Spearman correlation", fontsize=13, fontweight="bold")
ax.legend(fontsize=10)
ax.axline((0, 0), slope=1, color="gray", ls="--", alpha=0.5, zorder=1)
ax.set_aspect("equal")
lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]), max(ax.get_xlim()[1], ax.get_ylim()[1])]
ax.set_xlim(lims)
ax.set_ylim(lims)

fig.suptitle("Pairwise phase similarity: Affinity vs FC representation",
             fontsize=15, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig2_affinity_vs_fc_scatter.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("  Saved fig2")


# ===========================================================================
# FIGURE 3: Bar chart — within vs cross affinity correlation, per band
# ===========================================================================

print("\nFigure 3: Within vs cross condition, per band")

fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=False)

for ax_idx, metric in enumerate(["aff_pearson", "fc_pearson"]):
    ax = axes[ax_idx]
    within_means, cross_means = [], []
    within_sems, cross_sems = [], []

    for band in BANDS:
        w_vals = [r[metric] for r in records if r["band"] == band and r["condition"] == "within"]
        c_vals = [r[metric] for r in records if r["band"] == band and r["condition"] == "cross"]
        within_means.append(np.mean(w_vals) if w_vals else 0)
        cross_means.append(np.mean(c_vals) if c_vals else 0)
        within_sems.append(np.std(w_vals) / np.sqrt(len(w_vals)) if len(w_vals) > 1 else 0)
        cross_sems.append(np.std(c_vals) / np.sqrt(len(c_vals)) if len(c_vals) > 1 else 0)

    x = np.arange(len(BANDS))
    w = 0.35
    band_labels = [BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS]

    bars1 = ax.bar(x - w/2, within_means, w, yerr=within_sems, label="Within-condition",
                   color="tab:blue", alpha=0.8, capsize=3)
    bars2 = ax.bar(x + w/2, cross_means, w, yerr=cross_sems, label="Cross-condition",
                   color="tab:red", alpha=0.8, capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels(band_labels, fontsize=10)
    ax.set_ylabel("Pearson r", fontsize=12)
    title = "Affinity-based" if "aff" in metric else "FC-based"
    ax.set_title(f"{title} similarity", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)
    ax.axhline(0, color="gray", ls="-", alpha=0.3)

fig.suptitle("Within vs cross-condition similarity by frequency band",
             fontsize=15, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig3_within_vs_cross_per_band.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("  Saved fig3")


# ===========================================================================
# FIGURE 4: Effect size comparison — affinity vs FC
# ===========================================================================

print("\nFigure 4: Effect size comparison")

fig, ax = plt.subplots(figsize=(10, 6))

aff_effects, fc_effects = [], []
for band in BANDS:
    w_aff = [r["aff_pearson"] for r in records if r["band"] == band and r["condition"] == "within"]
    c_aff = [r["aff_pearson"] for r in records if r["band"] == band and r["condition"] == "cross"]
    w_fc = [r["fc_pearson"] for r in records if r["band"] == band and r["condition"] == "within"]
    c_fc = [r["fc_pearson"] for r in records if r["band"] == band and r["condition"] == "cross"]

    aff_effect = np.mean(w_aff) - np.mean(c_aff) if w_aff and c_aff else 0
    fc_effect = np.mean(w_fc) - np.mean(c_fc) if w_fc and c_fc else 0
    aff_effects.append(aff_effect)
    fc_effects.append(fc_effect)

x = np.arange(len(BANDS))
w = 0.35
band_labels = [BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS]

ax.bar(x - w/2, aff_effects, w, label="Affinity-based", color="tab:purple", alpha=0.8)
ax.bar(x + w/2, fc_effects, w, label="FC-based", color="tab:orange", alpha=0.8)
ax.set_xticks(x)
ax.set_xticklabels(band_labels, fontsize=11)
ax.set_ylabel("Effect size (within - cross mean r)", fontsize=12)
ax.set_title("Condition discrimination: Affinity vs FC representation",
             fontsize=14, fontweight="bold")
ax.legend(fontsize=11)
ax.axhline(0, color="gray", ls="-", alpha=0.3)

fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig4_effect_size_comparison.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("  Saved fig4")


# ===========================================================================
# FIGURE 5: Grand summary heatmap — ratio of effect sizes per patient x band
# ===========================================================================

print("\nFigure 5: Grand summary heatmap")

ratio_matrix = np.full((len(FOUR_PHASE_PATIENTS), len(BANDS)), np.nan)

for i, pat in enumerate(FOUR_PHASE_PATIENTS):
    for j, band in enumerate(BANDS):
        w_aff = [r["aff_pearson"] for r in records
                 if r["patient"] == pat and r["band"] == band and r["condition"] == "within"]
        c_aff = [r["aff_pearson"] for r in records
                 if r["patient"] == pat and r["band"] == band and r["condition"] == "cross"]
        w_fc = [r["fc_pearson"] for r in records
                if r["patient"] == pat and r["band"] == band and r["condition"] == "within"]
        c_fc = [r["fc_pearson"] for r in records
                if r["patient"] == pat and r["band"] == band and r["condition"] == "cross"]

        aff_effect = np.mean(w_aff) - np.mean(c_aff) if w_aff and c_aff else 0
        fc_effect = np.mean(w_fc) - np.mean(c_fc) if w_fc and c_fc else 0

        if abs(fc_effect) > 1e-6:
            ratio_matrix[i, j] = aff_effect / fc_effect
        else:
            ratio_matrix[i, j] = np.nan

fig, ax = plt.subplots(figsize=(10, 5))
band_labels = [BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS]

# Use diverging colormap centered at 1
vmax = max(abs(np.nanmax(ratio_matrix)), abs(np.nanmin(ratio_matrix)), 3)
im = ax.imshow(ratio_matrix, cmap="RdBu_r", aspect="auto",
               norm=matplotlib.colors.TwoSlopeNorm(vcenter=1.0, vmin=min(0, np.nanmin(ratio_matrix)),
                                                    vmax=max(2, np.nanmax(ratio_matrix))))

ax.set_xticks(range(len(BANDS)))
ax.set_xticklabels(band_labels, fontsize=11)
ax.set_yticks(range(len(FOUR_PHASE_PATIENTS)))
ax.set_yticklabels(FOUR_PHASE_PATIENTS, fontsize=11)
ax.set_xlabel("Frequency band", fontsize=12)
ax.set_ylabel("Patient", fontsize=12)

# Annotate cells
for i in range(ratio_matrix.shape[0]):
    for j in range(ratio_matrix.shape[1]):
        val = ratio_matrix[i, j]
        if not np.isnan(val):
            color = "white" if abs(val - 1) > 0.8 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=10,
                    fontweight="bold", color=color)

cbar = fig.colorbar(im, ax=ax, label="Affinity / FC effect ratio")
ax.set_title("Multiscale advantage ratio (>1 = affinity better at discriminating conditions)",
             fontsize=13, fontweight="bold")

fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig5_advantage_ratio_heatmap.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("  Saved fig5")


# ===========================================================================
# FINDINGS.md
# ===========================================================================

print("\n" + "=" * 70)
print("Generating FINDINGS.md")
print("=" * 70)

# Compute summary statistics
within_aff = [r["aff_pearson"] for r in records if r["condition"] == "within"]
cross_aff = [r["aff_pearson"] for r in records if r["condition"] == "cross"]
within_fc = [r["fc_pearson"] for r in records if r["condition"] == "within"]
cross_fc = [r["fc_pearson"] for r in records if r["condition"] == "cross"]

aff_separation = np.mean(within_aff) - np.mean(cross_aff)
fc_separation = np.mean(within_fc) - np.mean(cross_fc)

# Per-band summary
band_summaries = []
for band in BANDS:
    w_a = [r["aff_pearson"] for r in records if r["band"] == band and r["condition"] == "within"]
    c_a = [r["aff_pearson"] for r in records if r["band"] == band and r["condition"] == "cross"]
    w_f = [r["fc_pearson"] for r in records if r["band"] == band and r["condition"] == "within"]
    c_f = [r["fc_pearson"] for r in records if r["band"] == band and r["condition"] == "cross"]
    aff_e = np.mean(w_a) - np.mean(c_a)
    fc_e = np.mean(w_f) - np.mean(c_f)
    ratio = aff_e / fc_e if abs(fc_e) > 1e-6 else float("nan")
    band_summaries.append((band, np.mean(w_a), np.mean(c_a), aff_e, np.mean(w_f), np.mean(c_f), fc_e, ratio))

# Count where affinity wins
n_wins = sum(1 for _, _, _, ae, _, _, fe, _ in band_summaries if ae > fe)
n_ratio_gt1 = np.sum(ratio_matrix[~np.isnan(ratio_matrix)] > 1)
n_ratio_total = np.sum(~np.isnan(ratio_matrix))

findings = f"""# Multiscale Co-Classification Affinity Analysis — Findings

## Method
For each LRG dendrogram, we sweep {N_THRESHOLDS} evenly-spaced distance thresholds
and compute the co-classification matrix at each level (C[i,j] = 1 if nodes i,j
are in the same cluster). The **affinity matrix** A[i,j] = mean(C_k[i,j]) across
all k thresholds, giving a value in [0,1] representing the probability that two
nodes are grouped together across scales.

## Key Question
Does the multiscale affinity matrix capture condition-dependent structure better
than the raw MSC functional connectivity matrix?

## Results

### Overall separation (within-condition minus cross-condition Pearson r)
- **Affinity-based**: {aff_separation:.4f}
- **FC-based**: {fc_separation:.4f}
- **Ratio (affinity/FC)**: {aff_separation/fc_separation:.2f} (>1 means affinity better)

### Within vs cross-condition mean Pearson r
- Affinity within: {np.mean(within_aff):.4f} +/- {np.std(within_aff):.4f}
- Affinity cross:  {np.mean(cross_aff):.4f} +/- {np.std(cross_aff):.4f}
- FC within:       {np.mean(within_fc):.4f} +/- {np.std(within_fc):.4f}
- FC cross:        {np.mean(cross_fc):.4f} +/- {np.std(cross_fc):.4f}

### Per-band effect sizes

| Band | Aff within | Aff cross | Aff effect | FC within | FC cross | FC effect | Ratio |
|------|-----------|-----------|------------|----------|----------|-----------|-------|
"""

for band, aw, ac, ae, fw, fc_, fe, ratio in band_summaries:
    tex = BRAIN_BAND_TEX_DICT.get(band, band)
    findings += f"| {tex} | {aw:.4f} | {ac:.4f} | {ae:.4f} | {fw:.4f} | {fc_:.4f} | {fe:.4f} | {ratio:.2f} |\n"

findings += f"""
### Patient-level advantage ratio
- Cells where affinity outperforms FC (ratio > 1): {int(n_ratio_gt1)}/{int(n_ratio_total)}
- Bands where affinity has larger overall effect: {n_wins}/{len(BANDS)}

## Interpretation
"""

if aff_separation > fc_separation:
    findings += """The multiscale affinity representation provides **better** separation between
within-condition and cross-condition phase pairs than raw FC, indicating that
hierarchical grouping structure is more condition-specific than edge weights alone.
"""
else:
    findings += """The raw FC representation provides better separation between within-condition
and cross-condition phase pairs, suggesting that edge weights already capture
condition-dependent information that is partially lost in the hierarchical
co-classification averaging.
"""

findings += f"""
## Files
- `pairwise_results.csv` — Full pairwise comparison data
- `fig1_affinity_heatmaps_Pat02_alpha.png` — Example affinity matrices
- `fig2_affinity_vs_fc_scatter.png` — Scatter: affinity vs FC correlation
- `fig3_within_vs_cross_per_band.png` — Within/cross bars per band
- `fig4_effect_size_comparison.png` — Effect size: affinity vs FC
- `fig5_advantage_ratio_heatmap.png` — Per-patient per-band advantage ratio

Generated: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M')}
"""

findings_path = OUTPUT_DIR / "FINDINGS.md"
with open(findings_path, "w") as f:
    f.write(findings)
print(f"  Saved {findings_path}")

print("\n" + "=" * 70)
print("DONE — all outputs in", OUTPUT_DIR)
print("=" * 70)
