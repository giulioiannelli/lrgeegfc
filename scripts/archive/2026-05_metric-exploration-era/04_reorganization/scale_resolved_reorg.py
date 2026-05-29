#!/usr/bin/env python3
"""Scale-Resolved Reorganization Profiles.

Investigates how partition distance (VI, NMI) between phase pairs varies
across hierarchical scales (dendrogram thresholds / number of clusters).

Uses Variation of Information (VI) as the primary metric — a true metric
on partitions with information-theoretic interpretation. VI=0 means identical
partitions; higher values mean more reorganization.
"""

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster
from sklearn.metrics import normalized_mutual_info_score

from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    PHASE_LABELS,
    BRAIN_BAND_TEX_DICT,
)
from lrg_eegfc.config.paths import FIGURES_ROOT


# ── VI computation ────────────────────────────────────────────────────────
def compute_vi(labels1, labels2):
    """Variation of Information between two label vectors.

    VI = H(P) + H(Q) - 2*MI(P,Q)
    VI = 0 means identical partitions, higher = more different.
    """
    n = len(labels1)
    if n == 0:
        return 0.0
    classes1 = np.unique(labels1)
    classes2 = np.unique(labels2)
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


# ── Configuration ──────────────────────────────────────────────────────────
FC_METHOD = "msc"
N_THRESHOLDS = 50
OUTPUT_DIR = FIGURES_ROOT / "multiscale_investigation" / "scale_resolved_reorg"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]
BANDS = list(BRAIN_BANDS.keys())
PHASES = list(PHASE_LABELS)

# Phase pair categories
WITHIN_REST = [("rest_pre", "rest_post")]
WITHIN_TASK = [("task_learn", "task_test")]
CROSS_CONDITION = [
    ("rest_pre", "task_learn"),
    ("rest_pre", "task_test"),
    ("task_learn", "rest_post"),
    ("task_test", "rest_post"),
]

PHASE_PAIR_LABELS = {
    ("rest_pre", "rest_post"): "rest-rest",
    ("task_learn", "task_test"): "task-task",
    ("rest_pre", "task_learn"): "rest_pre-tLearn",
    ("rest_pre", "task_test"): "rest_pre-tTest",
    ("task_learn", "rest_post"): "tLearn-rest_post",
    ("task_test", "rest_post"): "tTest-rest_post",
}

PHASE_PAIR_CATEGORY = {}
for pp in WITHIN_REST:
    PHASE_PAIR_CATEGORY[pp] = "within-rest"
for pp in WITHIN_TASK:
    PHASE_PAIR_CATEGORY[pp] = "within-task"
for pp in CROSS_CONDITION:
    PHASE_PAIR_CATEGORY[pp] = "cross-condition"

ALL_PHASE_PAIRS = WITHIN_REST + WITHIN_TASK + CROSS_CONDITION

# Colors for phase pairs
PAIR_COLORS = {
    ("rest_pre", "rest_post"): "#2196F3",       # blue
    ("task_learn", "task_test"): "#FF9800",  # orange
    ("rest_pre", "task_learn"): "#E91E63",     # pink
    ("rest_pre", "task_test"): "#9C27B0",      # purple
    ("task_learn", "rest_post"): "#F44336",    # red
    ("task_test", "rest_post"): "#795548",     # brown
}

CAT_COLORS = {
    "within-rest": "#2196F3",
    "within-task": "#FF9800",
    "cross-condition": "#E91E63",
}


# ── Step 1: Load all LRG results ──────────────────────────────────────────
print("Loading LRG results...")
lrg_data = {}  # (patient, band, phase) -> LRGResult

for patient in PATIENTS:
    for band in BANDS:
        for phase in PHASES:
            result = load_lrg_result(patient, phase, band, FC_METHOD)
            if result is not None:
                lrg_data[(patient, band, phase)] = result

print(f"  Loaded {len(lrg_data)} LRG results")


# ── Step 2: Compute scale-resolved profiles ──────────────────────────────
print("Computing reorganization profiles...")

# For each (patient, band, phase_pair), compute VI/NMI at N_THRESHOLDS scales
# We use number-of-clusters as the x-axis (more interpretable than raw threshold)

records = []  # will become a DataFrame

for patient in PATIENTS:
    for band in BANDS:
        for p1, p2 in ALL_PHASE_PAIRS:
            key1 = (patient, band, p1)
            key2 = (patient, band, p2)
            if key1 not in lrg_data or key2 not in lrg_data:
                continue

            r1 = lrg_data[key1]
            r2 = lrg_data[key2]

            Z1 = r1.linkage_matrix
            Z2 = r2.linkage_matrix

            # Get the range of linkage distances for both dendrograms
            max_dist = max(Z1[:, 2].max(), Z2[:, 2].max())
            min_dist = min(Z1[:, 2].min(), Z2[:, 2].min())

            # Evenly spaced thresholds
            thresholds = np.linspace(min_dist, max_dist, N_THRESHOLDS + 2)[1:-1]

            for t in thresholds:
                labels1 = fcluster(Z1, t=t, criterion="distance")
                labels2 = fcluster(Z2, t=t, criterion="distance")

                # Handle different n_nodes: use the minimum
                n_min = min(len(labels1), len(labels2))
                l1 = labels1[:n_min]
                l2 = labels2[:n_min]

                n_clusters_1 = len(np.unique(l1))
                n_clusters_2 = len(np.unique(l2))
                n_clusters_avg = (n_clusters_1 + n_clusters_2) / 2

                vi = compute_vi(l1, l2)
                nmi = normalized_mutual_info_score(l1, l2)

                records.append({
                    "patient": patient,
                    "band": band,
                    "phase1": p1,
                    "phase2": p2,
                    "pair_label": PHASE_PAIR_LABELS[(p1, p2)],
                    "category": PHASE_PAIR_CATEGORY[(p1, p2)],
                    "threshold": t,
                    "n_clusters_avg": n_clusters_avg,
                    "n_clusters_1": n_clusters_1,
                    "n_clusters_2": n_clusters_2,
                    "vi": vi,
                    "nmi": nmi,
                })

df = pd.DataFrame(records)
print(f"  Computed {len(df)} profile points")

# Save raw data
csv_path = OUTPUT_DIR / "scale_resolved_reorg_profiles.csv"
df.to_csv(csv_path, index=False)
print(f"  Saved CSV: {csv_path}")


# ── Utility: bin profiles by n_clusters for averaging ────────────────────
def bin_profiles(sub_df, n_bins=30, metric="vi"):
    """Bin by n_clusters_avg and compute mean metric per bin."""
    if sub_df.empty:
        return np.array([]), np.array([]), np.array([])
    min_k = sub_df["n_clusters_avg"].min()
    max_k = sub_df["n_clusters_avg"].max()
    bins = np.linspace(min_k, max_k, n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    means = []
    sems = []
    for i in range(n_bins):
        mask = (sub_df["n_clusters_avg"] >= bins[i]) & (sub_df["n_clusters_avg"] < bins[i + 1])
        vals = sub_df.loc[mask, metric].values
        if len(vals) > 0:
            means.append(np.mean(vals))
            sems.append(np.std(vals) / np.sqrt(len(vals)) if len(vals) > 1 else 0)
        else:
            means.append(np.nan)
            sems.append(np.nan)
    return bin_centers, np.array(means), np.array(sems)


# ── Fig 1: VI profiles per band, averaged across patients ────────────────
print("Generating Figure 1: VI profiles per band...")

fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharey=True)
axes = axes.ravel()

for i, band in enumerate(BANDS):
    ax = axes[i]
    band_df = df[df["band"] == band]

    for pp in ALL_PHASE_PAIRS:
        p1, p2 = pp
        sub = band_df[(band_df["phase1"] == p1) & (band_df["phase2"] == p2)]
        if sub.empty:
            continue

        # Average across patients: group by threshold, then average
        grouped = sub.groupby("threshold").agg(
            vi_mean=("vi", "mean"),
            vi_sem=("vi", lambda x: x.std() / np.sqrt(len(x)) if len(x) > 1 else 0),
            n_clusters_mean=("n_clusters_avg", "mean"),
        ).reset_index()

        # Sort by n_clusters (descending = fine to coarse)
        grouped = grouped.sort_values("n_clusters_mean")

        ax.plot(
            grouped["n_clusters_mean"],
            grouped["vi_mean"],
            color=PAIR_COLORS[pp],
            label=PHASE_PAIR_LABELS[pp],
            linewidth=1.8,
        )
        ax.fill_between(
            grouped["n_clusters_mean"],
            grouped["vi_mean"] - grouped["vi_sem"],
            grouped["vi_mean"] + grouped["vi_sem"],
            color=PAIR_COLORS[pp],
            alpha=0.15,
        )

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=14)
    ax.set_xlabel("Number of clusters")
    ax.axhline(y=0, color="gray", linestyle=":", alpha=0.5, linewidth=0.8)
    if i % 3 == 0:
        ax.set_ylabel("VI (distance)")
    if i == 0:
        ax.legend(fontsize=7, loc="upper left")

fig.suptitle(
    "Scale-Resolved Reorganization Profiles (VI) — Averaged Across Patients",
    fontsize=15, fontweight="bold", y=0.98,
)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig1_path = OUTPUT_DIR / "fig1_vi_profiles_per_band.png"
fig.savefig(fig1_path, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {fig1_path}")


# ── Fig 2: Grand average within vs cross condition ──────────────────────
print("Generating Figure 2: Grand average within vs cross condition...")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for metric_idx, (metric, metric_label) in enumerate([("vi", "VI (distance)"), ("nmi", "NMI")]):
    ax = axes[metric_idx]

    for cat, cat_label in [("within-rest", "Within rest (rest_pre-rest_post)"),
                           ("within-task", "Within task (tLearn-tTest)"),
                           ("cross-condition", "Cross condition (4 pairs)")]:
        cat_df = df[df["category"] == cat]
        if cat_df.empty:
            continue

        bin_centers, means, sems = bin_profiles(cat_df, n_bins=35, metric=metric)
        valid = ~np.isnan(means)
        ax.plot(
            bin_centers[valid], means[valid],
            color=CAT_COLORS[cat],
            label=cat_label,
            linewidth=2.2,
        )
        ax.fill_between(
            bin_centers[valid],
            (means - sems)[valid],
            (means + sems)[valid],
            color=CAT_COLORS[cat],
            alpha=0.15,
        )

    ax.set_xlabel("Number of clusters (fine → coarse)", fontsize=12)
    ax.set_ylabel(metric_label, fontsize=12)
    ax.set_title(f"Grand Average {metric.upper()} by Condition Type", fontsize=13)
    ax.legend(fontsize=9)
    ax.axhline(y=0, color="gray", linestyle=":", alpha=0.5, linewidth=0.8)

fig.suptitle(
    "Within-Condition vs Cross-Condition Reorganization Profiles",
    fontsize=14, fontweight="bold", y=1.02,
)
fig.tight_layout()
fig2_path = OUTPUT_DIR / "fig2_within_vs_cross_condition.png"
fig.savefig(fig2_path, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {fig2_path}")


# ── Fig 3: Heatmap of scale of maximum reorganization ────────────────────
print("Generating Figure 3: Heatmap of max-reorganization scale...")

# For each band x phase_pair, find the n_clusters at which VI is maximized
# (averaged across patients) — higher VI = more reorganization
heatmap_data = []

for band in BANDS:
    for pp in ALL_PHASE_PAIRS:
        p1, p2 = pp
        sub = df[(df["band"] == band) & (df["phase1"] == p1) & (df["phase2"] == p2)]
        if sub.empty:
            heatmap_data.append({
                "band": band,
                "pair": PHASE_PAIR_LABELS[pp],
                "max_vi_n_clusters": np.nan,
                "max_vi": np.nan,
            })
            continue

        # Average across patients per threshold
        grouped = sub.groupby("threshold").agg(
            vi_mean=("vi", "mean"),
            n_clusters_mean=("n_clusters_avg", "mean"),
        ).reset_index()

        idx_max = grouped["vi_mean"].idxmax()
        heatmap_data.append({
            "band": band,
            "pair": PHASE_PAIR_LABELS[pp],
            "max_vi_n_clusters": grouped.loc[idx_max, "n_clusters_mean"],
            "max_vi": grouped.loc[idx_max, "vi_mean"],
        })

heatmap_df = pd.DataFrame(heatmap_data)

fig, axes = plt.subplots(1, 2, figsize=(16, 5))

# Panel A: n_clusters at max VI
pivot_ncl = heatmap_df.pivot(index="band", columns="pair", values="max_vi_n_clusters")
# Reorder bands
pivot_ncl = pivot_ncl.reindex(BANDS)
pair_order = [PHASE_PAIR_LABELS[pp] for pp in ALL_PHASE_PAIRS]
pivot_ncl = pivot_ncl[[p for p in pair_order if p in pivot_ncl.columns]]

im1 = axes[0].imshow(pivot_ncl.values, aspect="auto", cmap="viridis")
axes[0].set_xticks(range(len(pivot_ncl.columns)))
axes[0].set_xticklabels(pivot_ncl.columns, rotation=45, ha="right", fontsize=9)
axes[0].set_yticks(range(len(pivot_ncl.index)))
axes[0].set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in pivot_ncl.index], fontsize=11)
axes[0].set_title("N clusters at max reorganization\n(highest VI)", fontsize=12)
plt.colorbar(im1, ax=axes[0], label="N clusters", shrink=0.8)

# Annotate
for yi in range(len(pivot_ncl.index)):
    for xi in range(len(pivot_ncl.columns)):
        val = pivot_ncl.values[yi, xi]
        if not np.isnan(val):
            axes[0].text(xi, yi, f"{val:.0f}", ha="center", va="center",
                        fontsize=8, color="white" if val < np.nanmean(pivot_ncl.values) else "black")

# Panel B: max VI value
pivot_vi = heatmap_df.pivot(index="band", columns="pair", values="max_vi")
pivot_vi = pivot_vi.reindex(BANDS)
pivot_vi = pivot_vi[[p for p in pair_order if p in pivot_vi.columns]]

im2 = axes[1].imshow(pivot_vi.values, aspect="auto", cmap="YlOrRd")
axes[1].set_xticks(range(len(pivot_vi.columns)))
axes[1].set_xticklabels(pivot_vi.columns, rotation=45, ha="right", fontsize=9)
axes[1].set_yticks(range(len(pivot_vi.index)))
axes[1].set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in pivot_vi.index], fontsize=11)
axes[1].set_title("Maximum VI value\n(max reorganization strength)", fontsize=12)
plt.colorbar(im2, ax=axes[1], label="VI", shrink=0.8)

# Annotate
for yi in range(len(pivot_vi.index)):
    for xi in range(len(pivot_vi.columns)):
        val = pivot_vi.values[yi, xi]
        if not np.isnan(val):
            axes[1].text(xi, yi, f"{val:.2f}", ha="center", va="center",
                        fontsize=8, color="black")

fig.suptitle("Scale of Maximum Reorganization by Band and Phase Pair",
             fontsize=14, fontweight="bold", y=1.02)
fig.tight_layout()
fig3_path = OUTPUT_DIR / "fig3_max_reorg_scale_heatmap.png"
fig.savefig(fig3_path, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {fig3_path}")


# ── Fig 4: Per-patient profiles for alpha band ──────────────────────────
print("Generating Figure 4: Per-patient profiles (alpha)...")

alpha_df = df[df["band"] == "alpha"]
# Only patients with all 4 phases
full_patients = [p for p in PATIENTS if all(
    (p, "alpha", ph) in lrg_data for ph in PHASES
)]

n_patients = len(full_patients)
fig, axes = plt.subplots(2, max(n_patients, 2), figsize=(5 * max(n_patients, 2), 8),
                         sharey="row")

for col, patient in enumerate(full_patients):
    for row, (metric, metric_label) in enumerate([("vi", "VI (distance)"), ("nmi", "NMI")]):
        ax = axes[row, col]
        pat_df = alpha_df[alpha_df["patient"] == patient]

        for pp in ALL_PHASE_PAIRS:
            p1, p2 = pp
            sub = pat_df[(pat_df["phase1"] == p1) & (pat_df["phase2"] == p2)]
            if sub.empty:
                continue
            sub = sub.sort_values("n_clusters_avg")
            ax.plot(
                sub["n_clusters_avg"], sub[metric],
                color=PAIR_COLORS[pp],
                label=PHASE_PAIR_LABELS[pp] if row == 0 and col == 0 else None,
                linewidth=1.5,
            )

        ax.set_xlabel("N clusters" if row == 1 else "")
        if col == 0:
            ax.set_ylabel(metric_label)
        if row == 0:
            ax.set_title(patient, fontsize=12)
        ax.axhline(y=0, color="gray", linestyle=":", alpha=0.5, linewidth=0.8)

# Hide unused subplot columns
for col in range(n_patients, axes.shape[1]):
    for row in range(2):
        axes[row, col].set_visible(False)

axes[0, 0].legend(fontsize=7, loc="upper left")
fig.suptitle(
    r"Per-Patient Scale-Resolved Profiles — $\alpha$ band",
    fontsize=14, fontweight="bold", y=1.01,
)
fig.tight_layout()
fig4_path = OUTPUT_DIR / "fig4_per_patient_alpha.png"
fig.savefig(fig4_path, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {fig4_path}")


# ── Analysis Summary ─────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("ANALYSIS SUMMARY")
print("=" * 70)

# 1. Scale of maximum reorganization
print("\n--- Scale of Maximum Reorganization (n_clusters at max VI) ---")
for band in BANDS:
    band_rows = heatmap_df[heatmap_df["band"] == band]
    cross = band_rows[band_rows["pair"].isin([PHASE_PAIR_LABELS[pp] for pp in CROSS_CONDITION])]
    within = band_rows[band_rows["pair"].isin([PHASE_PAIR_LABELS[pp] for pp in WITHIN_REST + WITHIN_TASK])]
    print(f"  {band:12s}: cross-cond max VI at ~{cross['max_vi_n_clusters'].mean():.0f} clusters "
          f"(VI={cross['max_vi'].mean():.3f}), "
          f"within-cond at ~{within['max_vi_n_clusters'].mean():.0f} clusters "
          f"(VI={within['max_vi'].mean():.3f})")

# 2. Fine vs coarse scale reorganization
print("\n--- Fine vs Coarse Scale Reorganization ---")
# Split thresholds into thirds: fine (many clusters), medium, coarse (few clusters)
for cat in ["within-rest", "within-task", "cross-condition"]:
    cat_df = df[df["category"] == cat]
    if cat_df.empty:
        continue
    q33, q66 = cat_df["n_clusters_avg"].quantile([0.33, 0.66])
    fine = cat_df[cat_df["n_clusters_avg"] > q66]["vi"].mean()
    medium = cat_df[(cat_df["n_clusters_avg"] >= q33) & (cat_df["n_clusters_avg"] <= q66)]["vi"].mean()
    coarse = cat_df[cat_df["n_clusters_avg"] < q33]["vi"].mean()
    print(f"  {cat:18s}: fine VI={fine:.3f}, medium VI={medium:.3f}, coarse VI={coarse:.3f}")

# 3. Cross-condition vs within-condition comparison
print("\n--- Cross vs Within Condition (overall) ---")
for cat in ["within-rest", "within-task", "cross-condition"]:
    cat_df = df[df["category"] == cat]
    if cat_df.empty:
        continue
    print(f"  {cat:18s}: mean VI={cat_df['vi'].mean():.3f} +/- {cat_df['vi'].std():.3f}, "
          f"mean NMI={cat_df['nmi'].mean():.3f} +/- {cat_df['nmi'].std():.3f}")

# 4. Consistency across bands
print("\n--- Consistency Across Bands (mean VI for cross-condition) ---")
for band in BANDS:
    cross_df_band = df[(df["band"] == band) & (df["category"] == "cross-condition")]
    if cross_df_band.empty:
        continue
    print(f"  {band:12s}: mean VI={cross_df_band['vi'].mean():.3f} +/- {cross_df_band['vi'].std():.3f}")


# ── Save summary CSV ─────────────────────────────────────────────────────
summary_records = []
for band in BANDS:
    for cat in ["within-rest", "within-task", "cross-condition"]:
        sub = df[(df["band"] == band) & (df["category"] == cat)]
        if sub.empty:
            continue
        q33, q66 = sub["n_clusters_avg"].quantile([0.33, 0.66])
        fine = sub[sub["n_clusters_avg"] > q66]["vi"].mean()
        coarse = sub[sub["n_clusters_avg"] < q33]["vi"].mean()
        summary_records.append({
            "band": band,
            "category": cat,
            "mean_vi": sub["vi"].mean(),
            "std_vi": sub["vi"].std(),
            "mean_nmi": sub["nmi"].mean(),
            "std_nmi": sub["nmi"].std(),
            "fine_scale_vi": fine,
            "coarse_scale_vi": coarse,
            "fine_minus_coarse": fine - coarse,
        })

summary_df = pd.DataFrame(summary_records)
summary_csv = OUTPUT_DIR / "scale_resolved_summary.csv"
summary_df.to_csv(summary_csv, index=False)
print(f"\nSaved summary CSV: {summary_csv}")

# ── Save heatmap CSV ─────────────────────────────────────────────────────
heatmap_csv = OUTPUT_DIR / "max_reorg_scale_heatmap.csv"
heatmap_df.to_csv(heatmap_csv, index=False)
print(f"Saved heatmap CSV: {heatmap_csv}")


# ── FINDINGS.md ──────────────────────────────────────────────────────────
findings_lines = [
    "# Scale-Resolved Reorganization Profiles — Findings",
    "",
    "## Method",
    "",
    "For each patient/band/phase-pair, we cut both LRG dendrograms at 50 evenly-spaced",
    "thresholds spanning the linkage distance range. At each threshold we computed VI",
    "(Variation of Information) and NMI between the resulting partitions. The x-axis is",
    "the number of clusters (many = fine/micro scale, few = coarse/macro scale).",
    "",
    "VI is a true metric on partitions: VI=0 means identical partitions, and higher",
    "values indicate more reorganization. Unlike ARI (a similarity measure), VI has",
    "information-theoretic interpretation and satisfies the triangle inequality.",
    "",
    "## Key Findings",
    "",
]

# Compute key stats for findings
cross_df = df[df["category"] == "cross-condition"]
within_rest_df = df[df["category"] == "within-rest"]
within_task_df = df[df["category"] == "within-task"]

findings_lines.append("### 1. Cross-condition reorganization is stronger than within-condition")
findings_lines.append("")
findings_lines.append(f"- Cross-condition mean VI: {cross_df['vi'].mean():.3f} +/- {cross_df['vi'].std():.3f}")
findings_lines.append(f"- Within-rest mean VI: {within_rest_df['vi'].mean():.3f} +/- {within_rest_df['vi'].std():.3f}")
if not within_task_df.empty:
    findings_lines.append(f"- Within-task mean VI: {within_task_df['vi'].mean():.3f} +/- {within_task_df['vi'].std():.3f}")
findings_lines.append("- (Higher VI = more reorganization / more different partitions)")
findings_lines.append("")

# Fine vs coarse
findings_lines.append("### 2. Scale-dependent reorganization pattern")
findings_lines.append("")
for cat in ["within-rest", "within-task", "cross-condition"]:
    cat_df = df[df["category"] == cat]
    if cat_df.empty:
        continue
    q33, q66 = cat_df["n_clusters_avg"].quantile([0.33, 0.66])
    fine = cat_df[cat_df["n_clusters_avg"] > q66]["vi"].mean()
    coarse = cat_df[cat_df["n_clusters_avg"] < q33]["vi"].mean()
    direction = "Fine-grained scales show MORE reorganization" if fine > coarse else "Coarse-grained scales show MORE reorganization"
    findings_lines.append(f"- **{cat}**: fine VI={fine:.3f}, coarse VI={coarse:.3f} -> {direction}")
findings_lines.append("")

findings_lines.append("### 3. Band-specific patterns")
findings_lines.append("")
for band in BANDS:
    cross_band = df[(df["band"] == band) & (df["category"] == "cross-condition")]
    if cross_band.empty:
        continue
    findings_lines.append(f"- {band}: cross-condition mean VI = {cross_band['vi'].mean():.3f}")
findings_lines.append("")

findings_lines.append("### 4. Individual consistency")
findings_lines.append("")
findings_lines.append("See Fig 4 for per-patient alpha profiles. Quantitative summary:")
for patient in full_patients:
    pat_cross = df[(df["patient"] == patient) & (df["band"] == "alpha") & (df["category"] == "cross-condition")]
    pat_within = df[(df["patient"] == patient) & (df["band"] == "alpha") & (df["category"] == "within-rest")]
    if not pat_cross.empty and not pat_within.empty:
        findings_lines.append(f"- {patient}: cross VI={pat_cross['vi'].mean():.3f}, within-rest VI={pat_within['vi'].mean():.3f}")
findings_lines.append("")

findings_lines.append("## Figures")
findings_lines.append("")
findings_lines.append("- `fig1_vi_profiles_per_band.png` -- VI(scale) for all 6 phase pairs, one subplot per band")
findings_lines.append("- `fig2_within_vs_cross_condition.png` -- Grand average VI/NMI by condition type")
findings_lines.append("- `fig3_max_reorg_scale_heatmap.png` -- Heatmap of scale at max reorganization + max VI")
findings_lines.append("- `fig4_per_patient_alpha.png` -- Per-patient profiles for alpha band")
findings_lines.append("")
findings_lines.append("## Data")
findings_lines.append("")
findings_lines.append("- `scale_resolved_reorg_profiles.csv` -- Full profiles (all thresholds)")
findings_lines.append("- `scale_resolved_summary.csv` -- Summary statistics by band/category")
findings_lines.append("- `max_reorg_scale_heatmap.csv` -- Scale of maximum reorganization")

findings_path = OUTPUT_DIR / "FINDINGS.md"
findings_path.write_text("\n".join(findings_lines))
print(f"Saved FINDINGS: {findings_path}")

print("\nDone!")
