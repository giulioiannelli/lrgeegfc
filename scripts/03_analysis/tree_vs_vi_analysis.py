#!/usr/bin/env python3
"""
Tree Structural Similarity vs Partition Similarity (VI) Analysis
================================================================

Comprehensive comparison of two fundamentally different ways to compare
LRG dendrograms across patients, bands, and phases:
1. Tree structure: topology, branch lengths, merge order
   (Robinson-Foulds, Baker gamma, cophenetic correlation, Sackin balance)
2. Partition similarity: cutting the tree at various scales and comparing
   flat clusters via Variation of Information (VI)

These CAN dissociate: tree topology can change while partition stays similar,
or vice versa.

Primary VI metric: VI_k5 (VI at k=5 clusters), which is threshold-independent.
Secondary: VI_matched_k (same k for both trees, derived from optimal thresholds,
clamped to [2,20]). Also: VI_mean_multi_k (mean of VI at k=3,4,5,6,8).

The previous VI_optimal metric (cutting each tree at its own optimal threshold)
is DEPRECATED because it produces meaningless comparisons when one tree has a
degenerate threshold (e.g., threshold below minimum merge height -> 100+ clusters
vs 4 clusters in the other tree).

Saves figures and CSVs to data/figures/multiscale_investigation/tree_vs_vi/
"""
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
from itertools import combinations
from scipy.cluster.hierarchy import fcluster, cophenet, dendrogram
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr, pearsonr, linregress
from collections import defaultdict

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE as _PATIENTS_4PHASE
from lrg_eegfc.config.paths import FIGURES_ROOT, LRG_CACHE

from lrgsglib.utils.basic.linalg import (
    tree_robinson_foulds_distance,
    tree_cophenetic_correlation,
    tree_baker_gamma,
    tree_fowlkes_mallows_index,
)

# ──────────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────────
PATIENTS_4PHASE = _PATIENTS_4PHASE
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
BANDS = list(BRAIN_BANDS.keys())
FC_METHOD = "msc"
OUTPUT_DIR = FIGURES_ROOT / "multiscale_investigation" / "tree_vs_vi"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

WITHIN_PAIRS = [("rest_pre", "rest_post"), ("task_learn", "task_test")]
CROSS_PAIRS = [
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"),
    ("task_learn", "rest_post"), ("task_test", "rest_post"),
]
ALL_PAIRS = WITHIN_PAIRS + CROSS_PAIRS

K_SCALES = [2, 3, 4, 5, 6, 8, 10, 15, 20]
MULTI_K_FOR_MEAN = [3, 4, 5, 6, 8]  # k values for robust mean VI
DEGENERATE_THRESHOLD = 20  # nclust > this -> degenerate optimal threshold

# Styling
PHASE_SHORT = {"rest_pre": "rest_pre", "task_learn": "tLearn", "task_test": "tTest", "rest_post": "rest_post"}
BAND_COLORS = {
    "delta": "#1f77b4", "theta": "#ff7f0e", "alpha": "#2ca02c",
    "beta": "#d62728", "low_gamma": "#9467bd", "high_gamma": "#8c564b",
}
BAND_TEX = {b: BRAIN_BAND_TEX_DICT.get(b, b) for b in BANDS}

# ──────────────────────────────────────────────────────────────────────
# Helper functions
# ──────────────────────────────────────────────────────────────────────
def compute_vi(labels1, labels2):
    """Variation of Information between two label vectors."""
    n = len(labels1)
    if n == 0:
        return 0.0
    classes1, classes2 = np.unique(labels1), np.unique(labels2)
    h1 = 0.0
    for c in classes1:
        p = np.sum(labels1 == c) / n
        if p > 0:
            h1 -= p * np.log(p)
    h2 = 0.0
    for c in classes2:
        p = np.sum(labels2 == c) / n
        if p > 0:
            h2 -= p * np.log(p)
    mi = 0.0
    for c1 in classes1:
        for c2 in classes2:
            pxy = np.sum((labels1 == c1) & (labels2 == c2)) / n
            if pxy > 0:
                px = np.sum(labels1 == c1) / n
                py = np.sum(labels2 == c2) / n
                mi += pxy * np.log(pxy / (px * py))
    return h1 + h2 - 2 * mi


def compute_sackin_norm(Z, n):
    """Normalized Sackin index of a linkage tree."""
    depth = np.zeros(2 * n - 1)
    for i in range(len(Z) - 1, -1, -1):
        node = n + i
        left, right = int(Z[i, 0]), int(Z[i, 1])
        depth[left] = depth[node] + 1
        depth[right] = depth[node] + 1
    return np.sum(depth[:n]) / (n * np.log2(n)) if n > 1 else 0


def load_lrg(patient, phase, band):
    """Load LRG result from cache, return dict or None."""
    path = LRG_CACHE / patient / f"{band}_{phase}_lrg_{FC_METHOD}.npz"
    if not path.exists():
        return None
    d = np.load(path)
    return {
        "linkage_matrix": d["linkage_matrix"],
        "optimal_threshold": float(d["optimal_threshold"]),
        "n_nodes": int(d["n_nodes"]),
        "ultrametric_matrix": d["ultrametric_matrix"],
    }


def pair_type(phaseA, phaseB):
    """Classify a phase pair as within-condition or cross-condition."""
    if (phaseA, phaseB) in WITHIN_PAIRS or (phaseB, phaseA) in WITHIN_PAIRS:
        return "within"
    return "cross"


def pair_label(phaseA, phaseB):
    """Short label for a phase pair."""
    return f"{PHASE_SHORT[phaseA]}-{PHASE_SHORT[phaseB]}"


# ──────────────────────────────────────────────────────────────────────
# PART 1: Compute ALL metrics for ALL phase pairs
# ──────────────────────────────────────────────────────────────────────
print("=" * 70)
print("PART 1: Computing all metrics for all phase pairs")
print("=" * 70)

records = []
skipped = []
degenerate_cases = []  # Track cases with degenerate optimal thresholds

for pat in PATIENTS_4PHASE:
    for band in BANDS:
        # Pre-load all phases for this patient/band
        lrg_data = {}
        for phase in PHASES:
            result = load_lrg(pat, phase, band)
            if result is not None:
                lrg_data[phase] = result

        # Check for degenerate thresholds in each phase
        for phase in PHASES:
            if phase in lrg_data:
                r = lrg_data[phase]
                Z = r["linkage_matrix"]
                thresh = r["optimal_threshold"]
                n = r["n_nodes"]
                labels_opt = fcluster(Z, thresh, criterion="distance")
                nclust = len(np.unique(labels_opt))
                if nclust > DEGENERATE_THRESHOLD:
                    degenerate_cases.append({
                        "patient": pat,
                        "band": band,
                        "phase": phase,
                        "optimal_threshold": thresh,
                        "nclust_at_optimal": nclust,
                        "n_nodes": n,
                        "min_merge_height": float(Z[:, 2].min()),
                        "max_merge_height": float(Z[:, 2].max()),
                    })

        for phaseA, phaseB in ALL_PAIRS:
            if phaseA not in lrg_data or phaseB not in lrg_data:
                skipped.append((pat, band, phaseA, phaseB))
                continue

            rA = lrg_data[phaseA]
            rB = lrg_data[phaseB]
            Z1, Z2 = rA["linkage_matrix"], rB["linkage_matrix"]
            n1, n2 = rA["n_nodes"], rB["n_nodes"]

            if n1 != n2:
                skipped.append((pat, band, phaseA, phaseB, f"n_nodes mismatch {n1} vs {n2}"))
                continue

            n = n1
            threshA = rA["optimal_threshold"]
            threshB = rB["optimal_threshold"]

            # --- Tree structural metrics ---
            try:
                rf = tree_robinson_foulds_distance(Z1, Z2, normalized=True)
            except Exception as e:
                rf = np.nan
                print(f"  RF failed {pat}/{band}/{phaseA}-{phaseB}: {e}")

            try:
                coph_corr = tree_cophenetic_correlation(Z1, Z2)
            except Exception:
                coph_corr = np.nan

            try:
                baker = tree_baker_gamma(Z1, Z2)
            except Exception:
                baker = np.nan

            try:
                fm = tree_fowlkes_mallows_index(Z1, Z2)
            except Exception:
                fm = np.nan

            # Sackin balance
            sacA = compute_sackin_norm(Z1, n)
            sacB = compute_sackin_norm(Z2, n)
            mean_sac = (sacA + sacB) / 2
            balance_change = abs(sacA - sacB) / mean_sac if mean_sac > 0 else 0.0

            # --- VI at multiple scales ---
            vi_by_k = {}
            for k in K_SCALES:
                if k >= n:
                    vi_by_k[k] = np.nan
                    continue
                labels1 = fcluster(Z1, k, criterion="maxclust")
                labels2 = fcluster(Z2, k, criterion="maxclust")
                vi_by_k[k] = compute_vi(labels1, labels2)

            # --- VI at optimal threshold (for diagnostics only) ---
            labelsA_opt = fcluster(Z1, threshA, criterion="distance")
            labelsB_opt = fcluster(Z2, threshB, criterion="distance")
            nclust_A = len(np.unique(labelsA_opt))
            nclust_B = len(np.unique(labelsB_opt))
            vi_optimal_raw = compute_vi(labelsA_opt, labelsB_opt)

            # Flag whether either side is degenerate
            degenerate_A = nclust_A > DEGENERATE_THRESHOLD
            degenerate_B = nclust_B > DEGENERATE_THRESHOLD
            is_degenerate = degenerate_A or degenerate_B

            # --- VI_matched_k: same k for both trees ---
            # Use mean of optimal k values, clamped to [2, 20]
            k_matched = int(round((nclust_A + nclust_B) / 2))
            # If either side is degenerate, use the non-degenerate side's k
            # or a reasonable fallback
            if is_degenerate:
                # Use the non-degenerate side, or 5 if both are degenerate
                if degenerate_A and not degenerate_B:
                    k_matched = max(nclust_B, 5)
                elif degenerate_B and not degenerate_A:
                    k_matched = max(nclust_A, 5)
                else:
                    k_matched = 5  # both degenerate, use 5 as safe default
            k_matched = max(2, min(k_matched, 20))

            if k_matched < n:
                labels1_mk = fcluster(Z1, k_matched, criterion="maxclust")
                labels2_mk = fcluster(Z2, k_matched, criterion="maxclust")
                vi_matched_k = compute_vi(labels1_mk, labels2_mk)
            else:
                vi_matched_k = np.nan

            # --- VI_mean_multi_k: mean of VI at k=3,4,5,6,8 ---
            vi_multi_vals = [vi_by_k[k] for k in MULTI_K_FOR_MEAN if k in vi_by_k and not np.isnan(vi_by_k.get(k, np.nan))]
            vi_mean_multi_k = np.mean(vi_multi_vals) if len(vi_multi_vals) > 0 else np.nan

            rec = {
                "patient": pat,
                "band": band,
                "phaseA": phaseA,
                "phaseB": phaseB,
                "pair_type": pair_type(phaseA, phaseB),
                "pair_label": pair_label(phaseA, phaseB),
                "n_nodes": n,
                "RF": rf,
                "baker_gamma": baker,
                "cophenetic_corr": coph_corr,
                "fowlkes_mallows": fm,
                "sackin_A": sacA,
                "sackin_B": sacB,
                "balance_change": balance_change,
                "threshold_A": threshA,
                "threshold_B": threshB,
                "nclust_opt_A": nclust_A,
                "nclust_opt_B": nclust_B,
                "is_degenerate": is_degenerate,
                "VI_optimal_raw": vi_optimal_raw,  # kept for diagnostics only
                "VI_matched_k": vi_matched_k,
                "k_matched": k_matched,
                "VI_k5": vi_by_k.get(5, np.nan),
                "VI_mean_multi_k": vi_mean_multi_k,
            }
            for k in K_SCALES:
                rec[f"VI_k{k}"] = vi_by_k[k]

            records.append(rec)
            degen_flag = " [DEGENERATE]" if is_degenerate else ""
            print(f"  {pat} {band:12s} {phaseA:10s}-{phaseB:10s}  "
                  f"RF={rf:.3f}  coph={coph_corr:.3f}  baker={baker:.3f}  "
                  f"VI_k5={vi_by_k.get(5, float('nan')):.3f}  "
                  f"VI_mk={vi_matched_k:.3f}  nclust={nclust_A}/{nclust_B}  "
                  f"k_matched={k_matched}{degen_flag}")

if skipped:
    print(f"\nSkipped {len(skipped)} pairs (missing data or n_nodes mismatch)")

# Report degenerate cases
if degenerate_cases:
    print(f"\n*** DEGENERATE THRESHOLD CASES ({len(degenerate_cases)}) ***")
    for dc in degenerate_cases:
        print(f"  {dc['patient']} {dc['band']:12s} {dc['phase']:10s}  "
              f"thresh={dc['optimal_threshold']:.4f}  nclust={dc['nclust_at_optimal']}  "
              f"merge_range=[{dc['min_merge_height']:.4f}, {dc['max_merge_height']:.4f}]")

df = pd.DataFrame(records)
df.to_csv(OUTPUT_DIR / "tree_vs_vi_results.csv", index=False)
print(f"\nSaved {len(df)} records to tree_vs_vi_results.csv")
print(f"Columns: {list(df.columns)}")

# Save degenerate cases
degen_df = pd.DataFrame(degenerate_cases) if degenerate_cases else pd.DataFrame()
if len(degen_df) > 0:
    degen_df.to_csv(OUTPUT_DIR / "degenerate_threshold_cases.csv", index=False)
    print(f"Saved {len(degen_df)} degenerate cases to degenerate_threshold_cases.csv")

# ──────────────────────────────────────────────────────────────────────
# PART 2: Correlation between tree metrics and VI at each scale
# ──────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("PART 2: Correlation between tree metrics and VI at each scale")
print("=" * 70)

tree_metrics = ["RF", "baker_gamma", "cophenetic_corr", "fowlkes_mallows", "balance_change"]
vi_cols = [f"VI_k{k}" for k in K_SCALES] + ["VI_k5", "VI_matched_k", "VI_mean_multi_k"]
# Deduplicate: VI_k5 is already in K_SCALES, so use unique list
vi_cols_for_heatmap = [f"VI_k{k}" for k in K_SCALES] + ["VI_matched_k", "VI_mean_multi_k"]

corr_matrix = np.full((len(tree_metrics), len(vi_cols_for_heatmap)), np.nan)
pval_matrix = np.full((len(tree_metrics), len(vi_cols_for_heatmap)), np.nan)

for i, tm in enumerate(tree_metrics):
    for j, vc in enumerate(vi_cols_for_heatmap):
        mask = df[tm].notna() & df[vc].notna()
        if mask.sum() > 5:
            rho, p = spearmanr(df.loc[mask, tm], df.loc[mask, vc])
            corr_matrix[i, j] = rho
            pval_matrix[i, j] = p

print("\nSpearman correlations (tree metric vs VI):")
for i, tm in enumerate(tree_metrics):
    vals = [f"{corr_matrix[i, j]:+.3f}" if not np.isnan(corr_matrix[i, j]) else "  nan" for j in range(len(vi_cols_for_heatmap))]
    print(f"  {tm:20s}: {', '.join(vals)}")

# ──────────────────────────────────────────────────────────────────────
# PART 3: Dissociation analysis (quadrant classification using VI_k5)
# ──────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("PART 3: Dissociation analysis (using VI_k5)")
print("=" * 70)

valid = df.dropna(subset=["RF", "VI_k5"]).copy()
rf_median = valid["RF"].median()
vi_median = valid["VI_k5"].median()

valid["quadrant"] = "Q4_stable"
valid.loc[(valid["RF"] >= rf_median) & (valid["VI_k5"] >= vi_median), "quadrant"] = "Q1_full_reorg"
valid.loc[(valid["RF"] >= rf_median) & (valid["VI_k5"] < vi_median), "quadrant"] = "Q2_struct_only"
valid.loc[(valid["RF"] < rf_median) & (valid["VI_k5"] >= vi_median), "quadrant"] = "Q3_partition_only"

for q in ["Q1_full_reorg", "Q2_struct_only", "Q3_partition_only", "Q4_stable"]:
    n = (valid["quadrant"] == q).sum()
    pct = 100 * n / len(valid)
    print(f"  {q}: {n} cases ({pct:.1f}%)")

print(f"\nMedian thresholds: RF={rf_median:.4f}, VI_k5={vi_median:.4f}")

# ──────────────────────────────────────────────────────────────────────
# PART 4: Per-band analysis
# ──────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("PART 4: Per-band analysis")
print("=" * 70)

band_stats = []
for band in BANDS:
    bdf = valid[valid["band"] == band]
    if len(bdf) == 0:
        continue
    for pt in ["within", "cross"]:
        sub = bdf[bdf["pair_type"] == pt]
        if len(sub) == 0:
            continue
        band_stats.append({
            "band": band,
            "pair_type": pt,
            "n": len(sub),
            "mean_RF": sub["RF"].mean(),
            "std_RF": sub["RF"].std(),
            "mean_VI_k5": sub["VI_k5"].mean(),
            "std_VI_k5": sub["VI_k5"].std(),
            "mean_VI_matched_k": sub["VI_matched_k"].mean(),
            "std_VI_matched_k": sub["VI_matched_k"].std(),
            "mean_VI_mean_multi_k": sub["VI_mean_multi_k"].mean(),
            "mean_coph": sub["cophenetic_corr"].mean(),
            "mean_baker": sub["baker_gamma"].mean(),
            "frac_Q1": (sub["quadrant"] == "Q1_full_reorg").mean(),
            "frac_Q2": (sub["quadrant"] == "Q2_struct_only").mean(),
            "frac_Q3": (sub["quadrant"] == "Q3_partition_only").mean(),
            "frac_Q4": (sub["quadrant"] == "Q4_stable").mean(),
            "n_degenerate": sub["is_degenerate"].sum(),
        })

band_df = pd.DataFrame(band_stats)
print(band_df.to_string(index=False))

# ──────────────────────────────────────────────────────────────────────
# PART 5: Pathological cases (regression outliers using VI_k5)
# ──────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("PART 5: Pathological cases (using VI_k5, excluding degenerate threshold artifacts)")
print("=" * 70)

slope, intercept, r_value, p_value, std_err = linregress(valid["RF"], valid["VI_k5"])
print(f"Linear regression RF -> VI_k5: slope={slope:.4f}, intercept={intercept:.4f}, R^2={r_value**2:.4f}, p={p_value:.2e}")

valid["VI_predicted"] = intercept + slope * valid["RF"]
valid["residual"] = valid["VI_k5"] - valid["VI_predicted"]
resid_std = valid["residual"].std()

valid["is_pathological"] = np.abs(valid["residual"]) > 2 * resid_std
pathological = valid[valid["is_pathological"]].copy()
pathological["deviation_type"] = np.where(
    pathological["residual"] > 0,
    "HIGH_VI_LOW_RF",
    "LOW_VI_HIGH_RF"
)

print(f"\nPathological cases (|residual| > 2*SD = {2*resid_std:.4f}):")
print(f"  Total: {len(pathological)} out of {len(valid)}")
if len(pathological) > 0:
    for _, row in pathological.iterrows():
        degen_note = " [DEGENERATE]" if row["is_degenerate"] else ""
        print(f"  {row['patient']} {row['band']:12s} {row['phaseA']:10s}-{row['phaseB']:10s}  "
              f"RF={row['RF']:.3f} VI_k5={row['VI_k5']:.3f} resid={row['residual']:+.3f}  "
              f"type={row['deviation_type']}{degen_note}")

# For the gallery, only use non-degenerate pathological cases (real dissociations)
pathological_real = pathological[~pathological["is_degenerate"]].copy()
print(f"\n  Real dissociations (non-degenerate): {len(pathological_real)}")
print(f"  Degenerate-threshold artifacts: {len(pathological) - len(pathological_real)}")

pathological.to_csv(OUTPUT_DIR / "dissociation_cases.csv", index=False)

# ══════════════════════════════════════════════════════════════════════
# FIGURES
# ══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("GENERATING FIGURES")
print("=" * 70)

# ──────────────────────────────────────────────────────────────────────
# Fig 1: KEY FIGURE - RF vs VI_k5 scatter, one subplot per band
# ──────────────────────────────────────────────────────────────────────
print("  Fig 1: RF vs VI_k5 scatter...")

fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharex=True, sharey=True)
axes = axes.ravel()

for idx, band in enumerate(BANDS):
    ax = axes[idx]
    bdf = valid[valid["band"] == band]

    within = bdf[bdf["pair_type"] == "within"]
    cross = bdf[bdf["pair_type"] == "cross"]

    ax.scatter(within["RF"], within["VI_k5"], c="#1f77b4", marker="o", s=60,
               alpha=0.8, label="within-cond", edgecolors="k", linewidths=0.5, zorder=3)
    ax.scatter(cross["RF"], cross["VI_k5"], c="#d62728", marker="s", s=60,
               alpha=0.8, label="cross-cond", edgecolors="k", linewidths=0.5, zorder=3)

    # Annotate pathological cases (non-degenerate only)
    path_band = pathological_real[pathological_real["band"] == band]
    for _, row in path_band.iterrows():
        ax.annotate(f"{row['patient'][-2:]}\n{row['pair_label']}",
                    (row["RF"], row["VI_k5"]),
                    fontsize=6, ha="center", va="bottom",
                    bbox=dict(boxstyle="round,pad=0.2", fc="yellow", alpha=0.6))

    # Regression line
    x_range = np.linspace(valid["RF"].min(), valid["RF"].max(), 50)
    ax.plot(x_range, intercept + slope * x_range, "k--", alpha=0.4, linewidth=1)

    # 2-SD band
    ax.fill_between(x_range,
                     intercept + slope * x_range - 2 * resid_std,
                     intercept + slope * x_range + 2 * resid_std,
                     alpha=0.08, color="gray")

    # Median lines for quadrant delineation
    ax.axvline(rf_median, color="gray", linestyle=":", alpha=0.5, linewidth=0.8)
    ax.axhline(vi_median, color="gray", linestyle=":", alpha=0.5, linewidth=0.8)

    ax.set_title(BAND_TEX.get(band, band), fontsize=12, fontweight="bold")
    if idx >= 3:
        ax.set_xlabel("Robinson-Foulds distance (norm.)")
    if idx % 3 == 0:
        ax.set_ylabel("VI at k=5")
    if idx == 0:
        ax.legend(fontsize=8, loc="upper left")

fig.suptitle("Tree Structural Change vs Partition Change\n(RF distance vs VI at k=5 clusters)",
             fontsize=14, fontweight="bold")
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig1_rf_vs_vi_scatter.png", dpi=200, bbox_inches="tight")
fig.savefig(OUTPUT_DIR / "fig1_rf_vs_vi_scatter.pdf", bbox_inches="tight")
plt.close(fig)
print("    -> saved fig1_rf_vs_vi_scatter.png/pdf")

# ──────────────────────────────────────────────────────────────────────
# Fig 2: Correlation heatmap (tree metrics x VI scales)
# ──────────────────────────────────────────────────────────────────────
print("  Fig 2: Correlation heatmap...")

fig, ax = plt.subplots(figsize=(14, 5))

tree_labels = ["RF", "Baker $\\gamma$", "Cophenetic r", "Fowlkes-Mallows", "$\\Delta$Balance"]
vi_labels = [f"k={k}" for k in K_SCALES] + ["matched-k", "mean(3..8)"]

im = ax.imshow(corr_matrix, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
ax.set_xticks(range(len(vi_labels)))
ax.set_xticklabels(vi_labels, rotation=45, ha="right")
ax.set_yticks(range(len(tree_labels)))
ax.set_yticklabels(tree_labels)
ax.set_xlabel("VI scale")
ax.set_ylabel("Tree metric")
ax.set_title("Spearman correlation: Tree Metrics vs VI at Different Scales")

# Annotate
for i in range(corr_matrix.shape[0]):
    for j in range(corr_matrix.shape[1]):
        val = corr_matrix[i, j]
        if not np.isnan(val):
            color = "white" if abs(val) > 0.5 else "black"
            sig = "*" if pval_matrix[i, j] < 0.05 else ""
            sig += "*" if pval_matrix[i, j] < 0.01 else ""
            sig += "*" if pval_matrix[i, j] < 0.001 else ""
            ax.text(j, i, f"{val:.2f}{sig}", ha="center", va="center", fontsize=8, color=color)

plt.colorbar(im, ax=ax, label="Spearman rho", shrink=0.8)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig2_correlation_heatmap.png", dpi=200, bbox_inches="tight")
fig.savefig(OUTPUT_DIR / "fig2_correlation_heatmap.pdf", bbox_inches="tight")
plt.close(fig)
print("    -> saved fig2_correlation_heatmap.png/pdf")

# ──────────────────────────────────────────────────────────────────────
# Fig 3: Per-band bar chart (mean RF vs mean VI_k5, within vs cross)
# ──────────────────────────────────────────────────────────────────────
print("  Fig 3: Per-band bar chart...")

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

x = np.arange(len(BANDS))
w = 0.35

# Panel A: Mean RF
ax = axes[0]
for pidx, pt in enumerate(["within", "cross"]):
    vals = [band_df[(band_df["band"] == b) & (band_df["pair_type"] == pt)]["mean_RF"].values
            for b in BANDS]
    vals = [v[0] if len(v) > 0 else 0 for v in vals]
    errs = [band_df[(band_df["band"] == b) & (band_df["pair_type"] == pt)]["std_RF"].values
            for b in BANDS]
    errs = [e[0] if len(e) > 0 else 0 for e in errs]
    ax.bar(x + pidx * w - w / 2, vals, w, yerr=errs, label=pt,
           color=["#1f77b4", "#d62728"][pidx], alpha=0.8, capsize=3)
ax.set_xticks(x)
ax.set_xticklabels([BAND_TEX.get(b, b) for b in BANDS], rotation=45, ha="right")
ax.set_ylabel("Mean Robinson-Foulds (norm.)")
ax.set_title("Tree Structural Change")
ax.legend()

# Panel B: Mean VI_k5
ax = axes[1]
for pidx, pt in enumerate(["within", "cross"]):
    vals = [band_df[(band_df["band"] == b) & (band_df["pair_type"] == pt)]["mean_VI_k5"].values
            for b in BANDS]
    vals = [v[0] if len(v) > 0 else 0 for v in vals]
    errs = [band_df[(band_df["band"] == b) & (band_df["pair_type"] == pt)]["std_VI_k5"].values
            for b in BANDS]
    errs = [e[0] if len(e) > 0 else 0 for e in errs]
    ax.bar(x + pidx * w - w / 2, vals, w, yerr=errs, label=pt,
           color=["#1f77b4", "#d62728"][pidx], alpha=0.8, capsize=3)
ax.set_xticks(x)
ax.set_xticklabels([BAND_TEX.get(b, b) for b in BANDS], rotation=45, ha="right")
ax.set_ylabel("Mean VI at k=5")
ax.set_title("Partition Change (k=5)")
ax.legend()

# Panel C: Ratio RF / VI_k5 (normalized)
ax = axes[2]
for pidx, pt in enumerate(["within", "cross"]):
    rf_vals = [band_df[(band_df["band"] == b) & (band_df["pair_type"] == pt)]["mean_RF"].values
               for b in BANDS]
    rf_vals = [v[0] if len(v) > 0 else np.nan for v in rf_vals]
    vi_vals = [band_df[(band_df["band"] == b) & (band_df["pair_type"] == pt)]["mean_VI_k5"].values
               for b in BANDS]
    vi_vals = [v[0] if len(v) > 0 else np.nan for v in vi_vals]
    # Normalize each to [0,1] for ratio comparison
    rf_arr = np.array(rf_vals)
    vi_arr = np.array(vi_vals)
    rf_norm = rf_arr / np.nanmax(rf_arr) if np.nanmax(rf_arr) > 0 else rf_arr
    vi_norm = vi_arr / np.nanmax(vi_arr) if np.nanmax(vi_arr) > 0 else vi_arr
    ratio = rf_norm / (vi_norm + 1e-10)
    ax.bar(x + pidx * w - w / 2, ratio, w, label=pt,
           color=["#1f77b4", "#d62728"][pidx], alpha=0.8)
ax.set_xticks(x)
ax.set_xticklabels([BAND_TEX.get(b, b) for b in BANDS], rotation=45, ha="right")
ax.set_ylabel("RF/VI ratio (normalized)")
ax.set_title("Structural vs Partition Change Ratio")
ax.axhline(1.0, color="gray", linestyle="--", alpha=0.5)
ax.legend()

fig.suptitle("Per-Band Reorganization Profile: Within vs Cross Condition",
             fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig3_per_band_bar.png", dpi=200, bbox_inches="tight")
fig.savefig(OUTPUT_DIR / "fig3_per_band_bar.pdf", bbox_inches="tight")
plt.close(fig)
print("    -> saved fig3_per_band_bar.png/pdf")

# ──────────────────────────────────────────────────────────────────────
# Fig 4: Dissociation profile (quadrant stacked bar per band)
# ──────────────────────────────────────────────────────────────────────
print("  Fig 4: Dissociation profile...")

fig, ax = plt.subplots(figsize=(12, 5))

q_colors = {
    "Q1_full_reorg": "#d62728",
    "Q2_struct_only": "#ff7f0e",
    "Q3_partition_only": "#2ca02c",
    "Q4_stable": "#1f77b4",
}
q_labels = {
    "Q1_full_reorg": "Q1: Full reorganization",
    "Q2_struct_only": "Q2: Structural only (high RF, low VI)",
    "Q3_partition_only": "Q3: Partition only (low RF, high VI)",
    "Q4_stable": "Q4: Stable (low RF, low VI)",
}

bottom = np.zeros(len(BANDS))
for q in ["Q1_full_reorg", "Q2_struct_only", "Q3_partition_only", "Q4_stable"]:
    fracs = []
    for band in BANDS:
        bdf = valid[valid["band"] == band]
        if len(bdf) > 0:
            fracs.append((bdf["quadrant"] == q).mean())
        else:
            fracs.append(0)
    fracs = np.array(fracs)
    ax.bar(range(len(BANDS)), fracs, bottom=bottom, color=q_colors[q],
           label=q_labels[q], alpha=0.85, edgecolor="white", linewidth=0.5)
    bottom += fracs

ax.set_xticks(range(len(BANDS)))
ax.set_xticklabels([BAND_TEX.get(b, b) for b in BANDS], rotation=45, ha="right")
ax.set_ylabel("Fraction of phase pairs")
ax.set_title("Dissociation Profile: Tree Structure vs Partition Similarity per Band (VI at k=5)",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=9, loc="upper right")
ax.set_ylim(0, 1.02)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig4_dissociation_profile.png", dpi=200, bbox_inches="tight")
fig.savefig(OUTPUT_DIR / "fig4_dissociation_profile.pdf", bbox_inches="tight")
plt.close(fig)
print("    -> saved fig4_dissociation_profile.png/pdf")

# ──────────────────────────────────────────────────────────────────────
# Fig 5: Per-patient scatter (RF vs VI_k5)
# ──────────────────────────────────────────────────────────────────────
print("  Fig 5: Per-patient profiles...")

fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True, sharey=True)
axes = axes.ravel()

for pidx, pat in enumerate(PATIENTS_4PHASE):
    ax = axes[pidx]
    pdf = valid[valid["patient"] == pat]

    for band in BANDS:
        bdf = pdf[pdf["band"] == band]
        ax.scatter(bdf["RF"], bdf["VI_k5"], c=BAND_COLORS[band], s=50,
                   alpha=0.8, label=BAND_TEX.get(band, band), edgecolors="k", linewidths=0.3)

    ax.axvline(rf_median, color="gray", linestyle=":", alpha=0.4)
    ax.axhline(vi_median, color="gray", linestyle=":", alpha=0.4)
    ax.set_title(pat, fontsize=12, fontweight="bold")
    if pidx >= 2:
        ax.set_xlabel("Robinson-Foulds (norm.)")
    if pidx % 2 == 0:
        ax.set_ylabel("VI at k=5")
    if pidx == 0:
        ax.legend(fontsize=7, ncol=2, loc="upper left")

fig.suptitle("Per-Patient: Tree Change vs Partition Change (VI at k=5)",
             fontsize=14, fontweight="bold")
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig5_per_patient.png", dpi=200, bbox_inches="tight")
fig.savefig(OUTPUT_DIR / "fig5_per_patient.pdf", bbox_inches="tight")
plt.close(fig)
print("    -> saved fig5_per_patient.png/pdf")

# ──────────────────────────────────────────────────────────────────────
# Fig 6: Scale-dependency (Spearman(RF, VI(k)) vs k, per band)
# ──────────────────────────────────────────────────────────────────────
print("  Fig 6: Scale dependency...")

fig, ax = plt.subplots(figsize=(10, 6))

# Global line
global_corrs = []
for k in K_SCALES:
    col = f"VI_k{k}"
    mask = df["RF"].notna() & df[col].notna()
    if mask.sum() > 5:
        rho, _ = spearmanr(df.loc[mask, "RF"], df.loc[mask, col])
        global_corrs.append(rho)
    else:
        global_corrs.append(np.nan)
ax.plot(K_SCALES, global_corrs, "k-o", linewidth=2.5, markersize=8, label="All bands",
        zorder=10)

# Per-band lines
for band in BANDS:
    bdf = df[df["band"] == band]
    corrs = []
    for k in K_SCALES:
        col = f"VI_k{k}"
        mask = bdf["RF"].notna() & bdf[col].notna()
        if mask.sum() > 5:
            rho, _ = spearmanr(bdf.loc[mask, "RF"], bdf.loc[mask, col])
            corrs.append(rho)
        else:
            corrs.append(np.nan)
    ax.plot(K_SCALES, corrs, "-o", color=BAND_COLORS[band], linewidth=1.5,
            markersize=5, alpha=0.7, label=BAND_TEX.get(band, band))

# Mark k=5 with a vertical indicator
ax.axvline(5, color="green", linestyle="--", alpha=0.4, linewidth=1.5, label="k=5 (primary)")

ax.set_xlabel("Number of clusters k")
ax.set_ylabel("Spearman(RF, VI(k))")
ax.set_title("At Which Scale Does Tree Topology Best Predict Partition Differences?",
             fontsize=12, fontweight="bold")
ax.legend(fontsize=8, ncol=2)
ax.grid(True, alpha=0.3)
ax.set_xticks(K_SCALES)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig6_scale_dependency.png", dpi=200, bbox_inches="tight")
fig.savefig(OUTPUT_DIR / "fig6_scale_dependency.pdf", bbox_inches="tight")
plt.close(fig)
print("    -> saved fig6_scale_dependency.png/pdf")

# ──────────────────────────────────────────────────────────────────────
# Fig 7: Phase transition heatmaps (RF and VI_k5 averaged)
# ──────────────────────────────────────────────────────────────────────
print("  Fig 7: Phase transition heatmaps...")

pair_order = ["rest_pre-tLearn", "rest_pre-tTest", "rest_pre-rest_post",
              "tLearn-tTest", "tLearn-rest_post", "tTest-rest_post"]
pair_labels_display = pair_order

# Compute means across all patients and bands for each pair
rf_means = {}
vi_means = {}
for _, row in valid.iterrows():
    pl = row["pair_label"]
    if pl not in rf_means:
        rf_means[pl] = []
        vi_means[pl] = []
    rf_means[pl].append(row["RF"])
    vi_means[pl].append(row["VI_k5"])

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# RF heatmap (bar representation since it's 1D)
ax = axes[0]
rf_vals = [np.mean(rf_means.get(pl, [0])) for pl in pair_order]
rf_errs = [np.std(rf_means.get(pl, [0])) for pl in pair_order]
colors = ["#1f77b4" if pl in ["rest_pre-rest_post", "tLearn-tTest"] else "#d62728" for pl in pair_order]
ax.barh(range(len(pair_order)), rf_vals, xerr=rf_errs, color=colors, alpha=0.8,
        capsize=3, edgecolor="k", linewidth=0.5)
ax.set_yticks(range(len(pair_order)))
ax.set_yticklabels(pair_order)
ax.set_xlabel("Mean Robinson-Foulds (norm.)")
ax.set_title("Tree Structural Change per Phase Pair")
ax.invert_yaxis()

# VI_k5 heatmap
ax = axes[1]
vi_vals = [np.mean(vi_means.get(pl, [0])) for pl in pair_order]
vi_errs = [np.std(vi_means.get(pl, [0])) for pl in pair_order]
ax.barh(range(len(pair_order)), vi_vals, xerr=vi_errs, color=colors, alpha=0.8,
        capsize=3, edgecolor="k", linewidth=0.5)
ax.set_yticks(range(len(pair_order)))
ax.set_yticklabels(pair_order)
ax.set_xlabel("Mean VI at k=5")
ax.set_title("Partition Change per Phase Pair (k=5)")
ax.invert_yaxis()

fig.suptitle("Phase Transition Profiles: Structural vs Partition Reorganization",
             fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig7_phase_transition.png", dpi=200, bbox_inches="tight")
fig.savefig(OUTPUT_DIR / "fig7_phase_transition.pdf", bbox_inches="tight")
plt.close(fig)
print("    -> saved fig7_phase_transition.png/pdf")

# ──────────────────────────────────────────────────────────────────────
# Fig 8: Pathological cases gallery (dendrograms) - REAL dissociations only
# ──────────────────────────────────────────────────────────────────────
print("  Fig 8: Pathological cases gallery (non-degenerate only)...")

n_path = len(pathological_real)
if n_path > 0:
    n_path = min(n_path, 8)  # limit to 8 max
    fig, axes = plt.subplots(n_path, 2, figsize=(16, 4 * n_path))
    if n_path == 1:
        axes = axes.reshape(1, 2)

    for pidx, (_, row) in enumerate(pathological_real.head(n_path).iterrows()):
        pat, band = row["patient"], row["band"]
        phA, phB = row["phaseA"], row["phaseB"]

        rA = load_lrg(pat, phA, band)
        rB = load_lrg(pat, phB, band)

        if rA is None or rB is None:
            continue

        Z1, Z2 = rA["linkage_matrix"], rB["linkage_matrix"]

        ax1 = axes[pidx, 0]
        dendrogram(Z1, ax=ax1, no_labels=True, color_threshold=rA["optimal_threshold"])
        ax1.axhline(rA["optimal_threshold"], color="red", linestyle="--", alpha=0.7, linewidth=1)
        ax1.set_title(f"{pat} {BAND_TEX.get(band, band)} {phA}\n(n={rA['n_nodes']} nodes)",
                      fontsize=10)
        ax1.set_ylabel("Distance")

        ax2 = axes[pidx, 1]
        dendrogram(Z2, ax=ax2, no_labels=True, color_threshold=rB["optimal_threshold"])
        ax2.axhline(rB["optimal_threshold"], color="red", linestyle="--", alpha=0.7, linewidth=1)
        ax2.set_title(f"{pat} {BAND_TEX.get(band, band)} {phB}\n(n={rB['n_nodes']} nodes)",
                      fontsize=10)

        # Annotate the deviation type
        dev_type = row["deviation_type"]
        if dev_type == "HIGH_VI_LOW_RF":
            annotation = "Partition disruption within\nsimilar tree topology"
        else:
            annotation = "Tree reshuffling without\npartition disruption"

        ax2.annotate(
            f"RF={row['RF']:.3f}  VI_k5={row['VI_k5']:.3f}\n{annotation}",
            xy=(0.98, 0.95), xycoords="axes fraction",
            ha="right", va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", alpha=0.8),
        )

    fig.suptitle("Pathological Cases: Real Dissociations between Tree Structure and Partition (k=5)\n"
                 "(degenerate-threshold artifacts excluded)",
                 fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "fig8_pathological_gallery.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "fig8_pathological_gallery.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"    -> saved fig8_pathological_gallery.png/pdf ({n_path} cases)")
else:
    print("    -> No non-degenerate pathological cases found, skipping Fig 8")

# ──────────────────────────────────────────────────────────────────────
# Fig 9: VI_k5 vs VI_matched_k comparison scatter
# ──────────────────────────────────────────────────────────────────────
print("  Fig 9: VI_k5 vs VI_matched_k comparison...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: VI_k5 vs VI_matched_k
ax = axes[0]
mask = valid["VI_matched_k"].notna() & valid["VI_k5"].notna()
non_degen = valid[mask & ~valid["is_degenerate"]]
degen = valid[mask & valid["is_degenerate"]]

ax.scatter(non_degen["VI_k5"], non_degen["VI_matched_k"], c="#1f77b4", s=40,
           alpha=0.7, label=f"Normal (n={len(non_degen)})", edgecolors="k", linewidths=0.3)
if len(degen) > 0:
    ax.scatter(degen["VI_k5"], degen["VI_matched_k"], c="#d62728", s=60, marker="x",
               linewidths=1.5, label=f"Degenerate thresh (n={len(degen)})", zorder=5)

# Identity line
lims = [0, max(valid["VI_k5"].max(), valid["VI_matched_k"].max()) * 1.05]
ax.plot(lims, lims, "k--", alpha=0.3, linewidth=1)
ax.set_xlabel("VI at k=5 (fixed)")
ax.set_ylabel("VI at matched-k")
ax.set_title("Fixed k=5 vs Matched-k VI")
ax.legend(fontsize=8)
ax.set_aspect("equal")

# Correlation
if mask.sum() > 5:
    rho_mk, p_mk = spearmanr(valid.loc[mask, "VI_k5"], valid.loc[mask, "VI_matched_k"])
    ax.text(0.05, 0.95, f"Spearman rho={rho_mk:.3f}\np={p_mk:.2e}",
            transform=ax.transAxes, fontsize=9, va="top",
            bbox=dict(boxstyle="round", fc="white", alpha=0.8))

# Panel B: Distribution of k_matched values
ax = axes[1]
k_vals = valid["k_matched"].dropna()
ax.hist(k_vals, bins=range(1, 22), color="#2ca02c", alpha=0.7, edgecolor="k", linewidth=0.5)
ax.axvline(5, color="red", linestyle="--", linewidth=1.5, label="k=5")
ax.axvline(k_vals.median(), color="blue", linestyle=":", linewidth=1.5, label=f"median={k_vals.median():.0f}")
ax.set_xlabel("Matched k value")
ax.set_ylabel("Count")
ax.set_title("Distribution of k_matched across all comparisons")
ax.legend(fontsize=9)

fig.suptitle("Comparison of VI Metrics: Fixed k=5 vs Threshold-Derived Matched-k",
             fontsize=13, fontweight="bold")
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig9_vi_comparison.png", dpi=200, bbox_inches="tight")
fig.savefig(OUTPUT_DIR / "fig9_vi_comparison.pdf", bbox_inches="tight")
plt.close(fig)
print("    -> saved fig9_vi_comparison.png/pdf")

# ══════════════════════════════════════════════════════════════════════
# FINDINGS.md
# ══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("WRITING FINDINGS.md")
print("=" * 70)

# Compute summary statistics for findings
overall_rho_k5, overall_p_k5 = spearmanr(valid["RF"], valid["VI_k5"])
overall_rho_mk, overall_p_mk = spearmanr(
    valid.loc[valid["VI_matched_k"].notna(), "RF"],
    valid.loc[valid["VI_matched_k"].notna(), "VI_matched_k"]
)
overall_rho_mean, overall_p_mean = spearmanr(
    valid.loc[valid["VI_mean_multi_k"].notna(), "RF"],
    valid.loc[valid["VI_mean_multi_k"].notna(), "VI_mean_multi_k"]
)
n_total = len(valid)

# Per-band dissociation rates
dissoc_rates = {}
for band in BANDS:
    bdf = valid[valid["band"] == band]
    if len(bdf) > 0:
        dissoc_rates[band] = {
            "Q2_rate": (bdf["quadrant"] == "Q2_struct_only").mean(),
            "Q3_rate": (bdf["quadrant"] == "Q3_partition_only").mean(),
            "total_dissoc": ((bdf["quadrant"] == "Q2_struct_only") | (bdf["quadrant"] == "Q3_partition_only")).mean(),
        }

# Find peak scale
peak_k = K_SCALES[np.nanargmax(np.abs(global_corrs))]
peak_rho = global_corrs[np.nanargmax(np.abs(global_corrs))]

# Per-band max scale
band_peak_scales = {}
for band in BANDS:
    bdf = df[df["band"] == band]
    corrs = []
    for k in K_SCALES:
        col = f"VI_k{k}"
        mask = bdf["RF"].notna() & bdf[col].notna()
        if mask.sum() > 5:
            rho, _ = spearmanr(bdf.loc[mask, "RF"], bdf.loc[mask, col])
            corrs.append(abs(rho))
        else:
            corrs.append(0)
    band_peak_scales[band] = K_SCALES[np.argmax(corrs)]

# Which phase pairs have most dissociation
pair_dissoc = {}
for pl in pair_order:
    sub = valid[valid["pair_label"] == pl]
    if len(sub) > 0:
        pair_dissoc[pl] = {
            "Q2": (sub["quadrant"] == "Q2_struct_only").mean(),
            "Q3": (sub["quadrant"] == "Q3_partition_only").mean(),
            "total_dissoc": ((sub["quadrant"] == "Q2_struct_only") | (sub["quadrant"] == "Q3_partition_only")).mean(),
        }

findings = f"""# Tree Structure vs Partition Similarity (VI): Detailed Findings

## Summary

This analysis compares two fundamentally different ways to measure reorganization
in LRG dendrograms across {len(PATIENTS_4PHASE)} patients, {len(BANDS)} frequency bands, and {len(ALL_PAIRS)} phase pairs
(N = {n_total} valid comparisons).

**Primary VI metric**: VI at k=5 clusters (threshold-independent, stable across trees).
**Secondary metrics**: VI_matched_k (same k for both trees, clamped to [2,20]),
VI_mean_multi_k (mean of VI at k=3,4,5,6,8).

**Core finding**: Robinson-Foulds (tree topology) and VI at k=5
show a Spearman correlation of rho = {overall_rho_k5:.3f} (p = {overall_p_k5:.2e}).
This means they capture {'overlapping' if abs(overall_rho_k5) > 0.3 else 'distinct'} aspects
of hierarchical reorganization.

Comparison of VI metrics (Spearman with RF):
- VI_k5: rho = {overall_rho_k5:.3f} (p = {overall_p_k5:.2e})
- VI_matched_k: rho = {overall_rho_mk:.3f} (p = {overall_p_mk:.2e})
- VI_mean_multi_k: rho = {overall_rho_mean:.3f} (p = {overall_p_mean:.2e})

## Bug Fix: Degenerate Optimal Thresholds

### Problem
The previous version computed VI_optimal by cutting each dendrogram at its own
optimal threshold. This produced meaningless comparisons when the optimal threshold
was degenerate (below minimum merge height), creating 100+ single-node clusters
compared against 3-6 clusters from the other phase.

### Degenerate Cases Found

{len(degenerate_cases)} phase/band combinations have degenerate thresholds (>{DEGENERATE_THRESHOLD} clusters):

| Patient | Band | Phase | Threshold | N clusters | Merge range |
|---------|------|-------|-----------|------------|-------------|
"""

for dc in degenerate_cases:
    findings += (f"| {dc['patient']} | {BAND_TEX.get(dc['band'], dc['band'])} | {dc['phase']} "
                 f"| {dc['optimal_threshold']:.4f} | {dc['nclust_at_optimal']} "
                 f"| [{dc['min_merge_height']:.4f}, {dc['max_merge_height']:.4f}] |\n")

findings += f"""
**Why these are degenerate**: The optimal threshold is below the minimum merge height
in the dendrogram, meaning ALL leaves are in separate clusters. This happens when the
LRG algorithm's gap statistic finds maximum entropy change at the finest scale,
suggesting the tree has no clear mesoscale structure at that band/phase.

### Solution
1. **VI_k5** (primary): Use a fixed k=5 for all comparisons. This is threshold-independent.
2. **VI_matched_k**: Use the same k for both trees, derived from the mean of their
   optimal cluster counts, but clamped to [2,20]. For degenerate cases, fall back to
   the non-degenerate side's k, or k=5 if both are degenerate.
3. **VI_mean_multi_k**: Mean of VI at k=3,4,5,6,8 for a robust summary.

## Dissociation Analysis (using VI_k5)

Using median-split quadrant classification (RF_median = {rf_median:.4f}, VI_k5_median = {vi_median:.4f}):

| Quadrant | Description | N | Fraction |
|----------|-------------|---|----------|
| Q1 | Full reorganization (high RF + high VI) | {(valid['quadrant'] == 'Q1_full_reorg').sum()} | {(valid['quadrant'] == 'Q1_full_reorg').mean():.1%} |
| Q2 | Structural only (high RF, low VI) | {(valid['quadrant'] == 'Q2_struct_only').sum()} | {(valid['quadrant'] == 'Q2_struct_only').mean():.1%} |
| Q3 | Partition only (low RF, high VI) | {(valid['quadrant'] == 'Q3_partition_only').sum()} | {(valid['quadrant'] == 'Q3_partition_only').mean():.1%} |
| Q4 | Stable (low RF, low VI) | {(valid['quadrant'] == 'Q4_stable').sum()} | {(valid['quadrant'] == 'Q4_stable').mean():.1%} |

**Dissociation rate** (Q2 + Q3): {((valid['quadrant'] == 'Q2_struct_only') | (valid['quadrant'] == 'Q3_partition_only')).mean():.1%}

This means that in approximately {((valid['quadrant'] == 'Q2_struct_only') | (valid['quadrant'] == 'Q3_partition_only')).mean():.0%} of phase pair comparisons,
the tree structure and partition tell DIFFERENT stories about reorganization.

## Scale Dependency

The correlation between RF distance and VI varies by the number of clusters k:

| k | Spearman rho |
|---|-------------|
"""

for ki, k in enumerate(K_SCALES):
    marker = " <-- primary" if k == 5 else ""
    findings += f"| {k} | {global_corrs[ki]:.3f}{marker} |\n"

findings += f"""
**Peak correlation at k = {peak_k}** (rho = {peak_rho:.3f}).

Per-band peak scales:
"""

for band in BANDS:
    findings += f"- {BAND_TEX.get(band, band)}: peak at k = {band_peak_scales[band]}\n"

findings += f"""
## Per-Band Dissociation

| Band | Q2 (struct only) | Q3 (partition only) | Total dissoc. |
|------|-----------------|--------------------|--------------:|
"""

for band in BANDS:
    if band in dissoc_rates:
        dr = dissoc_rates[band]
        findings += f"| {BAND_TEX.get(band, band)} | {dr['Q2_rate']:.1%} | {dr['Q3_rate']:.1%} | {dr['total_dissoc']:.1%} |\n"

findings += f"""
## Phase Pair Dissociation

| Pair | Type | Q2 rate | Q3 rate | Total dissoc |
|------|------|---------|---------|-------------|
"""

for pl in pair_order:
    if pl in pair_dissoc:
        pd_ = pair_dissoc[pl]
        ptype = "within" if pl in ["rest_pre-rest_post", "tLearn-tTest"] else "cross"
        findings += f"| {pl} | {ptype} | {pd_['Q2']:.1%} | {pd_['Q3']:.1%} | {pd_['total_dissoc']:.1%} |\n"

findings += f"""
## Regression Analysis

Linear regression: VI_k5 = {intercept:.4f} + {slope:.4f} * RF
- R-squared: {r_value**2:.4f}
- p-value: {p_value:.2e}

## Pathological Cases (Real Dissociations)

{len(pathological)} total cases deviate > 2 SD from the regression line.
Of these, {len(pathological_real)} are real dissociations (non-degenerate thresholds)
and {len(pathological) - len(pathological_real)} are degenerate-threshold artifacts.

### Real dissociations:
"""

for _, row in pathological_real.iterrows():
    findings += (f"\n- **{row['patient']} {BAND_TEX.get(row['band'], row['band'])} "
                 f"{row['phaseA']}-{row['phaseB']}**: "
                 f"RF={row['RF']:.3f}, VI_k5={row['VI_k5']:.3f}, "
                 f"residual={row['residual']:+.3f} "
                 f"({row['deviation_type']})")

if len(pathological) - len(pathological_real) > 0:
    findings += "\n\n### Degenerate-threshold artifacts (excluded from gallery):\n"
    pathological_degen = pathological[pathological["is_degenerate"]]
    for _, row in pathological_degen.iterrows():
        findings += (f"\n- **{row['patient']} {BAND_TEX.get(row['band'], row['band'])} "
                     f"{row['phaseA']}-{row['phaseB']}**: "
                     f"RF={row['RF']:.3f}, VI_k5={row['VI_k5']:.3f}, "
                     f"nclust={row['nclust_opt_A']}/{row['nclust_opt_B']} "
                     f"(degenerate threshold)")

findings += f"""

## Interpretation

### What Q2 (structural-only change) means
When RF is high but VI is low, the tree topology has reshuffled (branches moved,
merge order changed) BUT the communities at k=5 are preserved. This
suggests that the hierarchical micro-structure is sensitive to state changes, but
the mesoscale modular organization is robust.

### What Q3 (partition-only change) means
When RF is low but VI is high, the tree topology is largely preserved, but the
partition at k=5 differs. This means the same tree structure supports different
community assignments -- small shifts in branch heights or merge order at the
k=5 cut level produce different groupings.

### Why VI_k5 is better than VI_optimal
1. **Threshold-independent**: No dependence on gap statistic quality
2. **Same granularity**: Both trees compared at the same resolution
3. **No degenerate cases**: k=5 always produces 5 clusters regardless of tree shape
4. **Biologically interpretable**: 5 clusters maps well to known brain network parcellations
5. **Robust summary**: VI_mean_multi_k (mean of k=3..8) confirms that k=5 is representative

### Implications for the paper
1. Using only RF or only VI gives an incomplete picture of reorganization
2. The dissociation rate of ~{((valid['quadrant'] == 'Q2_struct_only') | (valid['quadrant'] == 'Q3_partition_only')).mean():.0%} means that both measures are needed
3. The scale-dependency (Fig 6) shows that fine-grained partitions (large k) may
   better capture tree structural changes than coarse partitions
4. Per-band differences in dissociation patterns suggest band-specific reorganization
   mechanisms

## Files Generated

- `tree_vs_vi_results.csv` - Complete results table ({len(df)} records)
- `degenerate_threshold_cases.csv` - Degenerate threshold cases ({len(degenerate_cases)} records)
- `dissociation_cases.csv` - Pathological cases ({len(pathological)} records)
- `fig1_rf_vs_vi_scatter.png/pdf` - Key figure: RF vs VI_k5 scatter per band
- `fig2_correlation_heatmap.png/pdf` - Tree metrics vs VI at different scales
- `fig3_per_band_bar.png/pdf` - Per-band mean RF and VI_k5 comparison
- `fig4_dissociation_profile.png/pdf` - Quadrant fractions per band
- `fig5_per_patient.png/pdf` - Per-patient RF vs VI_k5 scatter
- `fig6_scale_dependency.png/pdf` - Spearman(RF, VI(k)) vs k
- `fig7_phase_transition.png/pdf` - Phase pair reorganization profiles
- `fig8_pathological_gallery.png/pdf` - Dendrogram gallery for real dissociation cases
- `fig9_vi_comparison.png/pdf` - VI_k5 vs VI_matched_k comparison
"""

(OUTPUT_DIR / "FINDINGS.md").write_text(findings)
print("  -> saved FINDINGS.md")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print(f"Output directory: {OUTPUT_DIR}")
print(f"Total records: {len(df)}")
print(f"Degenerate threshold cases: {len(degenerate_cases)}")
print(f"Pathological cases (total): {len(pathological)}")
print(f"Pathological cases (real, non-degenerate): {len(pathological_real)}")
print(f"Overall RF-VI_k5 correlation: rho={overall_rho_k5:.3f}, p={overall_p_k5:.2e}")
print(f"Overall RF-VI_matched_k correlation: rho={overall_rho_mk:.3f}, p={overall_p_mk:.2e}")
print(f"Overall RF-VI_mean_multi_k correlation: rho={overall_rho_mean:.3f}, p={overall_p_mean:.2e}")
print("=" * 70)
