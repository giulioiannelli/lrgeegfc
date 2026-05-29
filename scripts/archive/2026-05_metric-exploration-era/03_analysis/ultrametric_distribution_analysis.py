#!/usr/bin/env python3
"""Ultrametric Distance Distribution Fingerprints Analysis.

Investigates the statistical properties of ultrametric distance distributions
from LRG analysis to characterize multiscale brain network organization.
"""

import warnings

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import squareform
from itertools import combinations

from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE
from lrg_eegfc.config.paths import FIGURES_ROOT

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ---------- Configuration ----------
FC_METHOD = "msc"
ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]
FULL_PHASE_PATIENTS = PATIENTS_4PHASE  # have all 4 phases
BANDS = list(BRAIN_BANDS.keys())
PHASES = list(PHASE_LABELS)
OUTPUT_DIR = str(FIGURES_ROOT / "multiscale_investigation" / "ultrametric_distributions")
FIXED_BINS = 100  # fixed bin count for histograms / entropy

PHASE_COLORS = {
    "rest_pre": "#2196F3",
    "task_learn": "#FF5722",
    "task_test": "#FF9800",
    "rest_post": "#4CAF50",
}
PHASE_TYPE = {
    "rest_pre": "rest",
    "task_learn": "task",
    "task_test": "task",
    "rest_post": "rest",
}

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ---------- Helper functions ----------
def gini_coefficient(values):
    """Compute Gini coefficient of a distribution."""
    v = np.sort(np.abs(values))
    n = len(v)
    if n == 0 or v.sum() == 0:
        return 0.0
    index = np.arange(1, n + 1)
    return (2.0 * np.sum(index * v) / (n * np.sum(v))) - (n + 1) / n


def shannon_entropy_hist(values, n_bins=FIXED_BINS):
    """Shannon entropy of histogram-discretized distribution."""
    counts, _ = np.histogram(values, bins=n_bins)
    probs = counts / counts.sum()
    probs = probs[probs > 0]
    return -np.sum(probs * np.log2(probs))


def count_distinct_levels(linkage_matrix):
    """Count distinct merge heights in linkage."""
    heights = linkage_matrix[:, 2]
    return len(np.unique(np.round(heights, decimals=10)))


# ---------- Step 1: Load all LRG results and compute statistics ----------
print("=" * 70)
print("ULTRAMETRIC DISTANCE DISTRIBUTION FINGERPRINTS")
print("=" * 70)

records = []
all_data = {}  # (patient, band, phase) -> ultrametric_matrix

for patient in ALL_PATIENTS:
    for band in BANDS:
        for phase in PHASES:
            # Skip phases that don't exist for Pat_06, Pat_07
            if patient in ("Pat_06", "Pat_07") and phase in ("rest_pre", "rest_post"):
                continue

            result = load_lrg_result(patient, phase, band, FC_METHOD)
            if result is None:
                print(f"  MISSING: {patient} {phase} {band}")
                continue

            um = result.ultrametric_matrix  # condensed form
            key = (patient, band, phase)
            all_data[key] = um

            # Distribution moments
            mean_val = np.mean(um)
            std_val = np.std(um)
            skew_val = stats.skew(um)
            kurt_val = stats.kurtosis(um)

            # Shannon entropy
            h_val = shannon_entropy_hist(um)

            # Distinct hierarchical levels
            n_levels = count_distinct_levels(result.linkage_matrix)

            # Gini coefficient
            gini_val = gini_coefficient(um)

            # Percentiles
            p10, p25, p50, p75, p90 = np.percentile(um, [10, 25, 50, 75, 90])

            records.append({
                "patient": patient,
                "band": band,
                "phase": phase,
                "phase_type": PHASE_TYPE[phase],
                "n_nodes": result.n_nodes,
                "n_distances": len(um),
                "mean": mean_val,
                "std": std_val,
                "skewness": skew_val,
                "kurtosis": kurt_val,
                "entropy": h_val,
                "n_levels": n_levels,
                "gini": gini_val,
                "p10": p10,
                "p25": p25,
                "median": p50,
                "p75": p75,
                "p90": p90,
                "optimal_threshold": result.optimal_threshold,
            })

df = pd.DataFrame(records)
print(f"\nLoaded {len(df)} condition entries from {df['patient'].nunique()} patients")
print(f"Bands: {df['band'].nunique()}, Phases per patient: varies")

# Save full stats CSV
csv_path = os.path.join(OUTPUT_DIR, "ultrametric_distribution_stats.csv")
df.to_csv(csv_path, index=False, float_format="%.6f")
print(f"Saved stats to {csv_path}")


# ---------- Step 2: KS distances between phase pairs ----------
print("\nComputing KS distances between phase pairs...")

ks_records = []
for patient in ALL_PATIENTS:
    for band in BANDS:
        available_phases = [p for p in PHASES if (patient, band, p) in all_data]
        for p1, p2 in combinations(available_phases, 2):
            um1 = all_data[(patient, band, p1)]
            um2 = all_data[(patient, band, p2)]
            ks_stat, ks_pval = stats.ks_2samp(um1, um2)

            # Classify comparison type
            t1, t2 = PHASE_TYPE[p1], PHASE_TYPE[p2]
            if t1 == t2:
                comp_type = f"within-{t1}"
            else:
                comp_type = "cross-condition"

            ks_records.append({
                "patient": patient,
                "band": band,
                "phase_1": p1,
                "phase_2": p2,
                "comparison_type": comp_type,
                "ks_statistic": ks_stat,
                "ks_pvalue": ks_pval,
            })

df_ks = pd.DataFrame(ks_records)
ks_csv_path = os.path.join(OUTPUT_DIR, "ks_distances.csv")
df_ks.to_csv(ks_csv_path, index=False, float_format="%.6f")
print(f"Saved KS distances to {ks_csv_path}")


# ========================================================================
# FIGURE 1: Overlaid histograms per band for Pat_02
# ========================================================================
print("\n--- Figure 1: Overlaid histograms (Pat_02) ---")
patient_rep = "Pat_02"

fig, axes = plt.subplots(2, 3, figsize=(15, 9))
axes = axes.ravel()

for idx, band in enumerate(BANDS):
    ax = axes[idx]
    # Determine global range for this band
    all_vals = []
    for phase in PHASES:
        key = (patient_rep, band, phase)
        if key in all_data:
            all_vals.append(all_data[key])
    if not all_vals:
        continue
    vmin = min(v.min() for v in all_vals)
    vmax = max(v.max() for v in all_vals)
    bins = np.linspace(vmin, vmax, FIXED_BINS + 1)

    for phase in PHASES:
        key = (patient_rep, band, phase)
        if key not in all_data:
            continue
        um = all_data[key]
        ax.hist(um, bins=bins, alpha=0.45, label=phase,
                color=PHASE_COLORS[phase], density=True, edgecolor="none")

    ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]} band", fontsize=13)
    ax.set_xlabel("Ultrametric distance", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    if idx == 0:
        ax.legend(fontsize=9, framealpha=0.8)

fig.suptitle(f"Ultrametric Distance Distributions by Phase ({patient_rep})",
             fontsize=15, fontweight="bold", y=0.98)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig1_path = os.path.join(OUTPUT_DIR, "fig1_histograms_pat02.png")
fig.savefig(fig1_path, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {fig1_path}")


# ========================================================================
# FIGURE 2: Distribution moment profiles across phases (patient-averaged)
# ========================================================================
print("--- Figure 2: Moment profiles across phases ---")

# Use only patients with all 4 phases
df_full = df[df["patient"].isin(FULL_PHASE_PATIENTS)].copy()

metrics = ["mean", "skewness", "kurtosis"]
metric_labels = ["Mean", "Skewness", "Excess Kurtosis"]

fig, axes = plt.subplots(len(metrics), len(BANDS), figsize=(18, 10),
                         sharey="row")

for row, (metric, mlabel) in enumerate(zip(metrics, metric_labels)):
    for col, band in enumerate(BANDS):
        ax = axes[row, col]
        sub = df_full[df_full["band"] == band]
        grouped = sub.groupby("phase")[metric]
        means = grouped.mean().reindex(PHASES)
        sems = grouped.sem().reindex(PHASES)

        x = np.arange(len(PHASES))
        colors = [PHASE_COLORS[p] for p in PHASES]
        ax.bar(x, means, yerr=sems, color=colors, alpha=0.8,
               capsize=4, edgecolor="gray", linewidth=0.5)
        ax.set_xticks(x)
        if row == len(metrics) - 1:
            ax.set_xticklabels(PHASES, rotation=35, ha="right", fontsize=8)
        else:
            ax.set_xticklabels([])

        if row == 0:
            ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=12)
        if col == 0:
            ax.set_ylabel(mlabel, fontsize=11)

        ax.grid(axis="y", alpha=0.3)

fig.suptitle("Ultrametric Distribution Moments Across Phases\n(Mean +/- SEM, N=4 patients)",
             fontsize=14, fontweight="bold", y=1.01)
fig.tight_layout()
fig2_path = os.path.join(OUTPUT_DIR, "fig2_moment_profiles.png")
fig.savefig(fig2_path, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {fig2_path}")


# ========================================================================
# FIGURE 3: KS distance matrix (averaged across patients and bands)
# ========================================================================
print("--- Figure 3: KS distance heatmap ---")

# Use only full-phase patients
df_ks_full = df_ks[df_ks["patient"].isin(FULL_PHASE_PATIENTS)].copy()

# Build average KS matrix
ks_matrix = np.zeros((len(PHASES), len(PHASES)))
for i, p1 in enumerate(PHASES):
    for j, p2 in enumerate(PHASES):
        if i == j:
            continue
        mask = ((df_ks_full["phase_1"] == p1) & (df_ks_full["phase_2"] == p2)) | \
               ((df_ks_full["phase_1"] == p2) & (df_ks_full["phase_2"] == p1))
        vals = df_ks_full.loc[mask, "ks_statistic"]
        if len(vals) > 0:
            ks_matrix[i, j] = vals.mean()

fig, ax = plt.subplots(figsize=(6.5, 5.5))
im = ax.imshow(ks_matrix, cmap="YlOrRd", vmin=0)
ax.set_xticks(range(len(PHASES)))
ax.set_xticklabels(PHASES, fontsize=11)
ax.set_yticks(range(len(PHASES)))
ax.set_yticklabels(PHASES, fontsize=11)
cbar = fig.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label("Mean KS Statistic", fontsize=11)

# Annotate cells
for i in range(len(PHASES)):
    for j in range(len(PHASES)):
        val = ks_matrix[i, j]
        color = "white" if val > ks_matrix.max() * 0.6 else "black"
        ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                fontsize=11, fontweight="bold", color=color)

# Add comparison-type boxes
# rest-rest block: (0,3) and (3,0)
# task-task block: (1,2) and (2,1)
from matplotlib.patches import Rectangle
# within-rest
for (r, c) in [(0, 3), (3, 0)]:
    ax.add_patch(Rectangle((c - 0.5, r - 0.5), 1, 1, fill=False,
                            edgecolor="#2196F3", linewidth=2.5, linestyle="--"))
# within-task
for (r, c) in [(1, 2), (2, 1)]:
    ax.add_patch(Rectangle((c - 0.5, r - 0.5), 1, 1, fill=False,
                            edgecolor="#FF5722", linewidth=2.5, linestyle="--"))

ax.set_title("KS Distance Between Phase Pairs\n(Averaged over patients & bands)",
             fontsize=13, fontweight="bold")
fig.tight_layout()
fig3_path = os.path.join(OUTPUT_DIR, "fig3_ks_distance_matrix.png")
fig.savefig(fig3_path, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {fig3_path}")


# ========================================================================
# FIGURE 4: Skewness vs Kurtosis scatter, colored by phase type
# ========================================================================
print("--- Figure 4: Skewness vs Kurtosis scatter ---")

fig, ax = plt.subplots(figsize=(8, 6.5))

type_colors = {"rest": "#2196F3", "task": "#FF9800"}
markers_phase = {"rest_pre": "o", "task_learn": "s", "task_test": "D", "rest_post": "^"}

for phase in PHASES:
    sub = df[df["phase"] == phase]
    pt = PHASE_TYPE[phase]
    ax.scatter(sub["skewness"], sub["kurtosis"],
               c=type_colors[pt], marker=markers_phase[phase],
               s=60, alpha=0.7, edgecolors="gray", linewidth=0.4,
               label=f"{phase} ({pt})")

ax.set_xlabel("Skewness", fontsize=12)
ax.set_ylabel("Excess Kurtosis", fontsize=12)
ax.set_title("Ultrametric Distribution Shape: Skewness vs Kurtosis\nAll patients and bands",
             fontsize=13, fontweight="bold")
ax.legend(fontsize=10, framealpha=0.8)
ax.grid(alpha=0.3)
fig.tight_layout()
fig4_path = os.path.join(OUTPUT_DIR, "fig4_skewness_vs_kurtosis.png")
fig.savefig(fig4_path, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {fig4_path}")


# ========================================================================
# FIGURE 5: Number of distinct hierarchical levels per condition
# ========================================================================
print("--- Figure 5: Distinct hierarchical levels ---")

fig, axes = plt.subplots(1, len(BANDS), figsize=(18, 5), sharey=True)

for col, band in enumerate(BANDS):
    ax = axes[col]
    sub = df_full[df_full["band"] == band]

    # Individual patient lines
    for pat in FULL_PHASE_PATIENTS:
        pat_sub = sub[sub["patient"] == pat].set_index("phase").reindex(PHASES)
        ax.plot(range(len(PHASES)), pat_sub["n_levels"], "o-",
                alpha=0.4, markersize=5, color="gray")

    # Mean line
    grouped = sub.groupby("phase")["n_levels"]
    means = grouped.mean().reindex(PHASES)
    sems = grouped.sem().reindex(PHASES)
    x = np.arange(len(PHASES))
    ax.errorbar(x, means, yerr=sems, fmt="s-", color="black",
                markersize=8, capsize=5, linewidth=2, label="Mean +/- SEM")

    ax.set_xticks(x)
    ax.set_xticklabels(PHASES, rotation=35, ha="right", fontsize=9)
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=12)
    if col == 0:
        ax.set_ylabel("Distinct Hierarchical Levels", fontsize=11)
    ax.grid(axis="y", alpha=0.3)

fig.suptitle("Number of Distinct Hierarchical Levels Across Phases\n(Gray: individual patients, Black: mean)",
             fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout()
fig5_path = os.path.join(OUTPUT_DIR, "fig5_hierarchical_levels.png")
fig.savefig(fig5_path, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {fig5_path}")


# ========================================================================
# Summary statistics and findings
# ========================================================================
print("\n" + "=" * 70)
print("SUMMARY ANALYSIS")
print("=" * 70)

findings = []

# --- Q1: Do moments change systematically between rest and task? ---
print("\n--- Q1: Rest vs Task distribution moments ---")
rest_df = df_full[df_full["phase_type"] == "rest"]
task_df = df_full[df_full["phase_type"] == "task"]

for metric in ["mean", "std", "skewness", "kurtosis", "entropy", "gini"]:
    rest_vals = rest_df[metric].values
    task_vals = task_df[metric].values
    t_stat, t_pval = stats.ttest_ind(rest_vals, task_vals)
    u_stat, u_pval = stats.mannwhitneyu(rest_vals, task_vals, alternative="two-sided")
    effect_d = (rest_vals.mean() - task_vals.mean()) / np.sqrt(
        (rest_vals.var() + task_vals.var()) / 2)
    print(f"  {metric:12s}: rest={rest_vals.mean():.4f}+/-{rest_vals.std():.4f}  "
          f"task={task_vals.mean():.4f}+/-{task_vals.std():.4f}  "
          f"t={t_stat:+.3f} p={t_pval:.4f}  Cohen's d={effect_d:+.3f}")
    findings.append(f"  {metric}: rest={rest_vals.mean():.4f}+/-{rest_vals.std():.4f}, "
                    f"task={task_vals.mean():.4f}+/-{task_vals.std():.4f}, "
                    f"t={t_stat:+.3f}, p={t_pval:.4f}, d={effect_d:+.3f}")

# Band-specific rest vs task
print("\n--- Q1b: Rest vs Task moments by band ---")
band_findings = []
for band in BANDS:
    sub_rest = df_full[(df_full["phase_type"] == "rest") & (df_full["band"] == band)]
    sub_task = df_full[(df_full["phase_type"] == "task") & (df_full["band"] == band)]
    for metric in ["mean", "skewness", "kurtosis"]:
        rv = sub_rest[metric].values
        tv = sub_task[metric].values
        if len(rv) < 2 or len(tv) < 2:
            continue
        t_stat, t_pval = stats.ttest_ind(rv, tv)
        if t_pval < 0.1:
            msg = f"  {band:12s} {metric:12s}: t={t_stat:+.3f} p={t_pval:.4f}"
            print(msg)
            band_findings.append(msg)

# --- Q2: Stability of hierarchical levels ---
print("\n--- Q2: Hierarchical levels stability ---")
for band in BANDS:
    sub = df_full[df_full["band"] == band]
    cv_per_patient = []
    for pat in FULL_PHASE_PATIENTS:
        pat_levels = sub[sub["patient"] == pat]["n_levels"].values
        if len(pat_levels) > 1 and pat_levels.mean() > 0:
            cv_per_patient.append(pat_levels.std() / pat_levels.mean())
    mean_cv = np.mean(cv_per_patient) if cv_per_patient else float("nan")
    print(f"  {band:12s}: mean CV across patients = {mean_cv:.4f} "
          f"(mean levels = {sub['n_levels'].mean():.1f})")

# --- Q3: Shape consistency within conditions ---
print("\n--- Q3: Within-condition consistency (skewness/kurtosis) ---")
for phase in PHASES:
    sub = df_full[df_full["phase"] == phase]
    skew_cv = sub["skewness"].std() / abs(sub["skewness"].mean()) if abs(sub["skewness"].mean()) > 0.01 else float("inf")
    kurt_cv = sub["kurtosis"].std() / abs(sub["kurtosis"].mean()) if abs(sub["kurtosis"].mean()) > 0.01 else float("inf")
    print(f"  {phase:12s}: skew CV={skew_cv:.3f}, kurt CV={kurt_cv:.3f} "
          f"(skew mean={sub['skewness'].mean():.3f}, kurt mean={sub['kurtosis'].mean():.3f})")

# --- Q4: KS within vs cross condition ---
print("\n--- Q4: KS distance: within-condition vs cross-condition ---")
for comp_type in ["within-rest", "within-task", "cross-condition"]:
    sub = df_ks_full[df_ks_full["comparison_type"] == comp_type]
    print(f"  {comp_type:20s}: KS mean={sub['ks_statistic'].mean():.4f} +/- {sub['ks_statistic'].std():.4f} "
          f"(N={len(sub)}, sig fraction={( sub['ks_pvalue'] < 0.05).mean():.2f})")

# Compare within vs cross
within_ks = df_ks_full[df_ks_full["comparison_type"].str.startswith("within")]["ks_statistic"]
cross_ks = df_ks_full[df_ks_full["comparison_type"] == "cross-condition"]["ks_statistic"]
u_stat, u_pval = stats.mannwhitneyu(within_ks, cross_ks, alternative="two-sided")
print(f"\n  Within vs Cross KS: U={u_stat:.1f}, p={u_pval:.4f}")
print(f"    Within mean={within_ks.mean():.4f}, Cross mean={cross_ks.mean():.4f}")


# ========================================================================
# Write FINDINGS.md
# ========================================================================
findings_path = os.path.join(OUTPUT_DIR, "FINDINGS.md")
with open(findings_path, "w") as f:
    f.write("# Ultrametric Distance Distribution Fingerprints\n\n")
    f.write("## Overview\n\n")
    f.write(f"- **FC method**: {FC_METHOD}\n")
    f.write(f"- **Patients (full phases)**: {', '.join(FULL_PHASE_PATIENTS)}\n")
    f.write(f"- **All patients**: {', '.join(ALL_PATIENTS)}\n")
    f.write(f"- **Bands**: {', '.join(BANDS)}\n")
    f.write(f"- **Phases**: {', '.join(PHASES)}\n")
    f.write(f"- **Total conditions analyzed**: {len(df)}\n\n")

    f.write("## Q1: Do distribution moments change systematically between rest and task?\n\n")
    rest_df2 = df_full[df_full["phase_type"] == "rest"]
    task_df2 = df_full[df_full["phase_type"] == "task"]
    f.write("| Metric | Rest (mean +/- std) | Task (mean +/- std) | t-stat | p-value | Cohen's d |\n")
    f.write("|--------|-------------------|-------------------|--------|---------|----------|\n")
    for metric in ["mean", "std", "skewness", "kurtosis", "entropy", "gini"]:
        rv = rest_df2[metric].values
        tv = task_df2[metric].values
        t_stat, t_pval = stats.ttest_ind(rv, tv)
        d = (rv.mean() - tv.mean()) / np.sqrt((rv.var() + tv.var()) / 2)
        sig = " *" if t_pval < 0.05 else ""
        f.write(f"| {metric} | {rv.mean():.4f} +/- {rv.std():.4f} | "
                f"{tv.mean():.4f} +/- {tv.std():.4f} | {t_stat:+.3f} | "
                f"{t_pval:.4f}{sig} | {d:+.3f} |\n")

    f.write("\n## Q2: Stability of distinct hierarchical levels\n\n")
    f.write("CV (coefficient of variation) of the number of distinct hierarchical levels\n")
    f.write("across phases within each patient, per band:\n\n")
    f.write("| Band | Mean CV | Mean Levels |\n")
    f.write("|------|---------|-------------|\n")
    for band in BANDS:
        sub = df_full[df_full["band"] == band]
        cv_per_patient = []
        for pat in FULL_PHASE_PATIENTS:
            pat_levels = sub[sub["patient"] == pat]["n_levels"].values
            if len(pat_levels) > 1 and pat_levels.mean() > 0:
                cv_per_patient.append(pat_levels.std() / pat_levels.mean())
        mean_cv = np.mean(cv_per_patient) if cv_per_patient else float("nan")
        f.write(f"| {band} | {mean_cv:.4f} | {sub['n_levels'].mean():.1f} |\n")

    f.write("\nLow CV indicates the number of hierarchical levels is relatively stable ")
    f.write("across experimental conditions within patients, supporting its use as a fingerprint.\n")

    f.write("\n## Q3: Within-condition shape consistency\n\n")
    f.write("CV of skewness and kurtosis across patients within each phase:\n\n")
    f.write("| Phase | Skew Mean | Skew CV | Kurt Mean | Kurt CV |\n")
    f.write("|-------|-----------|---------|-----------|--------|\n")
    for phase in PHASES:
        sub = df_full[df_full["phase"] == phase]
        skew_cv = sub["skewness"].std() / abs(sub["skewness"].mean()) if abs(sub["skewness"].mean()) > 0.01 else float("inf")
        kurt_cv = sub["kurtosis"].std() / abs(sub["kurtosis"].mean()) if abs(sub["kurtosis"].mean()) > 0.01 else float("inf")
        f.write(f"| {phase} | {sub['skewness'].mean():.3f} | {skew_cv:.3f} | "
                f"{sub['kurtosis'].mean():.3f} | {kurt_cv:.3f} |\n")

    f.write("\n## Q4: KS distance within vs cross condition\n\n")
    f.write("| Comparison | Mean KS | Std KS | N | Fraction sig (p<0.05) |\n")
    f.write("|------------|---------|--------|---|----------------------|\n")
    for comp_type in ["within-rest", "within-task", "cross-condition"]:
        sub = df_ks_full[df_ks_full["comparison_type"] == comp_type]
        f.write(f"| {comp_type} | {sub['ks_statistic'].mean():.4f} | "
                f"{sub['ks_statistic'].std():.4f} | {len(sub)} | "
                f"{(sub['ks_pvalue'] < 0.05).mean():.2f} |\n")

    u_stat, u_pval = stats.mannwhitneyu(within_ks, cross_ks, alternative="two-sided")
    f.write(f"\nMann-Whitney U test (within vs cross): U={u_stat:.1f}, p={u_pval:.4f}\n")
    if u_pval < 0.05:
        f.write("**Cross-condition comparisons show significantly different KS distances** ")
        f.write("than within-condition comparisons, suggesting ultrametric distributions ")
        f.write("capture condition-dependent reorganization.\n")
    else:
        f.write("No significant difference between within- and cross-condition KS distances.\n")

    f.write("\n## Figures\n\n")
    f.write("1. `fig1_histograms_pat02.png` - Overlaid ultrametric distance histograms per band (Pat_02)\n")
    f.write("2. `fig2_moment_profiles.png` - Distribution moments across phases (patient-averaged)\n")
    f.write("3. `fig3_ks_distance_matrix.png` - KS distance heatmap between phase pairs\n")
    f.write("4. `fig4_skewness_vs_kurtosis.png` - Skewness vs kurtosis scatter by condition\n")
    f.write("5. `fig5_hierarchical_levels.png` - Distinct hierarchical levels per condition\n")

print(f"\nSaved findings to {findings_path}")
print("\nDone!")
