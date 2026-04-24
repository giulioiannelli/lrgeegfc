#!/usr/bin/env python3
"""Dendrogram Tree Balance Analysis.

Computes Sackin and Colless indices from LRG dendrograms to characterize
hierarchical tree shape across patients, phases, and frequency bands.

Output: data/figures/multiscale_investigation/tree_balance/
"""

import warnings
warnings.filterwarnings("ignore")

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
from scipy import stats

from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.config.const import (
    BRAIN_BANDS, PHASE_LABELS, BRAIN_BAND_TEX_DICT, BRAIN_BANDS_NAMES,
)
from lrg_eegfc.config.paths import FIGURES_ROOT

# ── Configuration ──────────────────────────────────────────────────────
FC_METHOD = "msc"
PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]
ALL_PHASES = list(PHASE_LABELS)  # rest_pre, task_learn, task_test, rest_post
BANDS = list(BRAIN_BANDS.keys())

# Pat_06 and Pat_07 only have task phases
PATIENT_PHASES = {
    p: ALL_PHASES if p not in ("Pat_06", "Pat_07")
    else ["task_learn", "task_test"]
    for p in PATIENTS
}

OUTPUT_DIR = FIGURES_ROOT / "multiscale_investigation" / "tree_balance"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Phase colors
PHASE_COLORS = {
    "rest_pre": "#2196F3",
    "task_learn": "#FF9800",
    "task_test": "#F44336",
    "rest_post": "#4CAF50",
}
PHASE_SHORT = {
    "rest_pre": "Rest-Pre",
    "task_learn": "Task-Learn",
    "task_test": "Task-Test",
    "rest_post": "Rest-Post",
}

# ── Tree balance computation ───────────────────────────────────────────
def compute_tree_balance(Z, n):
    """Compute Sackin and Colless indices from scipy linkage matrix.

    Parameters
    ----------
    Z : ndarray, shape (n-1, 4)
        Scipy linkage matrix.
    n : int
        Number of leaves (original observations).

    Returns
    -------
    dict with sackin, colless, sackin_norm, colless_norm, max_depth, mean_depth
    """
    n_internal = len(Z)  # = n - 1

    # Subtree sizes
    subtree_size = np.ones(2 * n - 1)
    for i in range(n_internal):
        left, right = int(Z[i, 0]), int(Z[i, 1])
        subtree_size[n + i] = subtree_size[left] + subtree_size[right]

    # Depths: process from root downward
    depth = np.zeros(2 * n - 1)
    root = 2 * n - 2
    depth[root] = 0
    for i in range(n_internal - 1, -1, -1):
        node = n + i
        left, right = int(Z[i, 0]), int(Z[i, 1])
        depth[left] = depth[node] + 1
        depth[right] = depth[node] + 1

    leaf_depths = depth[:n]
    sackin = float(np.sum(leaf_depths))
    max_depth = float(np.max(leaf_depths))
    mean_depth = float(np.mean(leaf_depths))

    # Colless index
    colless = 0.0
    for i in range(n_internal):
        left, right = int(Z[i, 0]), int(Z[i, 1])
        colless += abs(subtree_size[left] - subtree_size[right])

    # Normalizations
    # Sackin: divide by n*log2(n) (expected for balanced tree ~ n*log2(n))
    sackin_norm = sackin / (n * np.log2(n)) if n > 1 else 0.0
    # Colless: divide by n*(n-1)/2 (maximum possible)
    colless_norm = colless / (n * (n - 1) / 2) if n > 1 else 0.0

    # Expected Sackin for a random binary tree ~ 2*n*ln(n)
    sackin_random_expected = 2 * n * np.log(n) if n > 1 else 1.0
    sackin_ratio = sackin / sackin_random_expected

    return {
        "sackin": sackin,
        "colless": colless,
        "sackin_norm": sackin_norm,
        "colless_norm": colless_norm,
        "sackin_ratio": sackin_ratio,
        "max_depth": max_depth,
        "mean_depth": mean_depth,
    }


# ── Data collection ────────────────────────────────────────────────────
print("Loading LRG results and computing tree balance metrics...")
rows = []
missing = []

for patient in PATIENTS:
    for band in BANDS:
        for phase in PATIENT_PHASES[patient]:
            result = load_lrg_result(patient, phase, band, FC_METHOD)
            if result is None:
                missing.append((patient, band, phase))
                continue

            Z = result.linkage_matrix
            n = result.n_nodes
            metrics = compute_tree_balance(Z, n)

            is_rest = phase in ("rest_pre", "rest_post")
            rows.append({
                "patient": patient,
                "band": band,
                "phase": phase,
                "condition": "rest" if is_rest else "task",
                "n_nodes": n,
                **metrics,
            })

df = pd.DataFrame(rows)
print(f"Collected {len(df)} conditions, {len(missing)} missing")
if missing:
    print(f"  Missing: {missing[:10]}{'...' if len(missing) > 10 else ''}")

# Save CSV
csv_path = OUTPUT_DIR / "tree_balance_results.csv"
df.to_csv(csv_path, index=False)
print(f"Saved results to {csv_path}")

# Print summary stats
print("\n── Summary Statistics ──")
for metric in ["sackin_norm", "colless_norm"]:
    print(f"\n{metric}:")
    summary = df.groupby("phase")[metric].agg(["mean", "std", "count"])
    print(summary.to_string())

# ── Statistical tests ──────────────────────────────────────────────────
print("\n── Statistical Tests: Rest vs Task ──")
rest_df = df[df["condition"] == "rest"]
task_df = df[df["condition"] == "task"]

for metric in ["sackin_norm", "colless_norm", "sackin_ratio"]:
    rest_vals = rest_df[metric].values
    task_vals = task_df[metric].values
    u_stat, u_p = stats.mannwhitneyu(rest_vals, task_vals, alternative="two-sided")
    t_stat, t_p = stats.ttest_ind(rest_vals, task_vals)
    print(f"\n{metric}:")
    print(f"  Rest: {rest_vals.mean():.4f} +/- {rest_vals.std():.4f} (n={len(rest_vals)})")
    print(f"  Task: {task_vals.mean():.4f} +/- {task_vals.std():.4f} (n={len(task_vals)})")
    print(f"  Mann-Whitney U={u_stat:.1f}, p={u_p:.4f}")
    print(f"  t-test t={t_stat:.3f}, p={t_p:.4f}")

# Band-specific tests
print("\n── Band-specific Rest vs Task ──")
for band in BANDS:
    band_df = df[df["band"] == band]
    r = band_df[band_df["condition"] == "rest"]["sackin_norm"].values
    t = band_df[band_df["condition"] == "task"]["sackin_norm"].values
    if len(r) >= 3 and len(t) >= 3:
        u_stat, u_p = stats.mannwhitneyu(r, t, alternative="two-sided")
        print(f"  {band:12s}: rest={r.mean():.3f}+/-{r.std():.3f}, task={t.mean():.3f}+/-{t.std():.3f}, U={u_stat:.0f}, p={u_p:.3f}")

# ── Figure 1: Sackin_norm across phases per band ───────────────────────
print("\nGenerating figures...")

fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharey=True)
axes = axes.ravel()

for idx, band in enumerate(BANDS):
    ax = axes[idx]
    band_df = df[df["band"] == band]

    phase_means = []
    phase_stds = []
    for phase in ALL_PHASES:
        vals = band_df[band_df["phase"] == phase]["sackin_norm"]
        phase_means.append(vals.mean() if len(vals) > 0 else np.nan)
        phase_stds.append(vals.std() if len(vals) > 0 else np.nan)

    x = np.arange(len(ALL_PHASES))
    colors = [PHASE_COLORS[p] for p in ALL_PHASES]
    bars = ax.bar(x, phase_means, yerr=phase_stds, capsize=4,
                  color=colors, alpha=0.7, edgecolor="black", linewidth=0.5)

    # Overlay patient dots
    for pi, phase in enumerate(ALL_PHASES):
        vals = band_df[band_df["phase"] == phase]["sackin_norm"]
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(vals))
        ax.scatter(pi + jitter, vals, color="black", s=15, alpha=0.6, zorder=5)

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels([PHASE_SHORT[p] for p in ALL_PHASES], rotation=30, ha="right", fontsize=8)
    if idx % 3 == 0:
        ax.set_ylabel("Sackin (normalized)", fontsize=10)
    ax.grid(axis="y", alpha=0.3)

fig.suptitle("Dendrogram Tree Balance: Sackin Index by Phase and Band", fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])
path1 = OUTPUT_DIR / "fig1_sackin_by_phase_band.png"
fig.savefig(path1, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {path1}")

# ── Figure 2: Colless_norm across phases per band ──────────────────────
fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharey=True)
axes = axes.ravel()

for idx, band in enumerate(BANDS):
    ax = axes[idx]
    band_df = df[df["band"] == band]

    phase_means = []
    phase_stds = []
    for phase in ALL_PHASES:
        vals = band_df[band_df["phase"] == phase]["colless_norm"]
        phase_means.append(vals.mean() if len(vals) > 0 else np.nan)
        phase_stds.append(vals.std() if len(vals) > 0 else np.nan)

    x = np.arange(len(ALL_PHASES))
    colors = [PHASE_COLORS[p] for p in ALL_PHASES]
    ax.bar(x, phase_means, yerr=phase_stds, capsize=4,
           color=colors, alpha=0.7, edgecolor="black", linewidth=0.5)

    for pi, phase in enumerate(ALL_PHASES):
        vals = band_df[band_df["phase"] == phase]["colless_norm"]
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(vals))
        ax.scatter(pi + jitter, vals, color="black", s=15, alpha=0.6, zorder=5)

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels([PHASE_SHORT[p] for p in ALL_PHASES], rotation=30, ha="right", fontsize=8)
    if idx % 3 == 0:
        ax.set_ylabel("Colless (normalized)", fontsize=10)
    ax.grid(axis="y", alpha=0.3)

fig.suptitle("Dendrogram Tree Balance: Colless Index by Phase and Band", fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])
path2 = OUTPUT_DIR / "fig2_colless_by_phase_band.png"
fig.savefig(path2, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {path2}")

# ── Figure 3: Rest vs Task boxplot ────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(14, 5))

for ax, metric, label in zip(
    axes,
    ["sackin_norm", "colless_norm", "sackin_ratio"],
    ["Sackin (normalized)", "Colless (normalized)", "Sackin / E[random]"],
):
    rest_vals = df[df["condition"] == "rest"][metric].values
    task_vals = df[df["condition"] == "task"][metric].values

    bp = ax.boxplot(
        [rest_vals, task_vals],
        labels=["Rest", "Task"],
        patch_artist=True,
        widths=0.5,
    )
    bp["boxes"][0].set_facecolor("#2196F3")
    bp["boxes"][1].set_facecolor("#FF9800")
    for box in bp["boxes"]:
        box.set_alpha(0.6)

    # Overlay individual points
    for i, vals in enumerate([rest_vals, task_vals]):
        jitter = np.random.default_rng(42).uniform(-0.1, 0.1, len(vals))
        ax.scatter(np.full(len(vals), i + 1) + jitter, vals,
                   color="black", s=10, alpha=0.4, zorder=5)

    u_stat, u_p = stats.mannwhitneyu(rest_vals, task_vals, alternative="two-sided")
    sig = "***" if u_p < 0.001 else "**" if u_p < 0.01 else "*" if u_p < 0.05 else "n.s."
    ax.set_title(f"{label}\np={u_p:.4f} ({sig})", fontsize=11)
    ax.set_ylabel(label)
    ax.grid(axis="y", alpha=0.3)

fig.suptitle("Rest vs Task: Tree Balance Comparison", fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.93])
path3 = OUTPUT_DIR / "fig3_rest_vs_task_boxplot.png"
fig.savefig(path3, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {path3}")

# ── Figure 4: Band-specific balance profiles ──────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

for ax, metric, label in zip(
    axes,
    ["sackin_norm", "colless_norm"],
    ["Sackin (normalized)", "Colless (normalized)"],
):
    for phase in ALL_PHASES:
        means = []
        sems = []
        for band in BANDS:
            vals = df[(df["band"] == band) & (df["phase"] == phase)][metric]
            means.append(vals.mean() if len(vals) > 0 else np.nan)
            sems.append(vals.sem() if len(vals) > 1 else 0)
        ax.errorbar(
            range(len(BANDS)), means, yerr=sems,
            marker="o", label=PHASE_SHORT[phase],
            color=PHASE_COLORS[phase], linewidth=2, markersize=6, capsize=3,
        )

    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=11)
    ax.set_ylabel(label, fontsize=11)
    ax.set_xlabel("Frequency Band", fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

fig.suptitle("Tree Balance Across Frequency Bands", fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.93])
path4 = OUTPUT_DIR / "fig4_band_profiles.png"
fig.savefig(path4, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {path4}")

# ── Figure 5: Per-patient profiles ────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(16, 9))
axes = axes.ravel()

for idx, patient in enumerate(PATIENTS):
    ax = axes[idx]
    pat_df = df[df["patient"] == patient]

    for phase in PATIENT_PHASES[patient]:
        means = []
        for band in BANDS:
            vals = pat_df[(pat_df["band"] == band) & (pat_df["phase"] == phase)]["sackin_norm"]
            means.append(vals.values[0] if len(vals) > 0 else np.nan)
        ax.plot(
            range(len(BANDS)), means,
            marker="o", label=PHASE_SHORT[phase],
            color=PHASE_COLORS[phase], linewidth=2, markersize=5,
        )

    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=9)
    ax.set_title(patient, fontsize=12, fontweight="bold")
    ax.set_ylabel("Sackin (norm)", fontsize=9)
    ax.legend(fontsize=7, loc="best")
    ax.grid(alpha=0.3)

fig.suptitle("Per-Patient Sackin Index Profiles", fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])
path5 = OUTPUT_DIR / "fig5_patient_profiles.png"
fig.savefig(path5, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {path5}")

# ── Figure 6: Sackin vs Colless scatter ────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 7))

for phase in ALL_PHASES:
    sub = df[df["phase"] == phase]
    ax.scatter(
        sub["sackin_norm"], sub["colless_norm"],
        c=PHASE_COLORS[phase], label=PHASE_SHORT[phase],
        s=50, alpha=0.7, edgecolors="black", linewidth=0.3,
    )

# Correlation
r, p = stats.pearsonr(df["sackin_norm"], df["colless_norm"])
ax.set_xlabel("Sackin Index (normalized)", fontsize=12)
ax.set_ylabel("Colless Index (normalized)", fontsize=12)
ax.set_title(f"Sackin vs Colless (r={r:.3f}, p={p:.2e})", fontsize=13, fontweight="bold")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# Fit line
slope, intercept = np.polyfit(df["sackin_norm"], df["colless_norm"], 1)
x_range = np.linspace(df["sackin_norm"].min(), df["sackin_norm"].max(), 100)
ax.plot(x_range, slope * x_range + intercept, "k--", alpha=0.5, linewidth=1)

fig.tight_layout()
path6 = OUTPUT_DIR / "fig6_sackin_vs_colless.png"
fig.savefig(path6, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {path6}")

# ── FINDINGS.md ────────────────────────────────────────────────────────
# Compute key findings programmatically
rest_sackin = df[df["condition"] == "rest"]["sackin_norm"]
task_sackin = df[df["condition"] == "task"]["sackin_norm"]
u_sackin, p_sackin = stats.mannwhitneyu(rest_sackin, task_sackin, alternative="two-sided")

rest_colless = df[df["condition"] == "rest"]["colless_norm"]
task_colless = df[df["condition"] == "task"]["colless_norm"]
u_colless, p_colless = stats.mannwhitneyu(rest_colless, task_colless, alternative="two-sided")

r_sc, p_sc = stats.pearsonr(df["sackin_norm"], df["colless_norm"])

# Band with highest/lowest balance
band_sackin = df.groupby("band")["sackin_norm"].mean()
most_balanced_band = band_sackin.idxmin()
least_balanced_band = band_sackin.idxmax()

# Phase means
phase_sackin = df.groupby("phase")["sackin_norm"].mean()

# Patient consistency: ICC-like measure (coefficient of variation across patients)
patient_means = df.groupby("patient")["sackin_norm"].mean()
cv_patients = patient_means.std() / patient_means.mean()

# Band-by-condition interaction
band_condition = df.groupby(["band", "condition"])["sackin_norm"].mean().unstack()
band_deltas = band_condition["task"] - band_condition["rest"]

findings = f"""# Dendrogram Tree Balance Analysis - Findings

## Overview
Computed Sackin and Colless tree-balance indices for {len(df)} LRG dendrograms
across {len(PATIENTS)} patients, {len(BANDS)} frequency bands, and 4 phases.
FC method: {FC_METHOD}.

## Key Metrics
- **Sackin index**: Sum of leaf depths. Higher = more unbalanced (caterpillar-like).
- **Colless index**: Sum of |left_size - right_size| at each internal node.
- Both normalized for tree size to allow cross-condition comparison.

## Results

### 1. Rest vs Task (Global)
| Metric | Rest (mean +/- SD) | Task (mean +/- SD) | Mann-Whitney p |
|--------|-------------------|-------------------|----------------|
| Sackin (norm) | {rest_sackin.mean():.4f} +/- {rest_sackin.std():.4f} | {task_sackin.mean():.4f} +/- {task_sackin.std():.4f} | {p_sackin:.4f} |
| Colless (norm) | {rest_colless.mean():.4f} +/- {rest_colless.std():.4f} | {task_colless.mean():.4f} +/- {task_colless.std():.4f} | {p_colless:.4f} |

**Interpretation**: {"Task phases show MORE balanced trees (lower Sackin/Colless) than rest" if task_sackin.mean() < rest_sackin.mean() else "Task phases show LESS balanced trees (higher Sackin/Colless) than rest" if task_sackin.mean() > rest_sackin.mean() else "No clear difference between rest and task"}.
{f"This is statistically significant (p < 0.05)." if p_sackin < 0.05 else "This difference is not statistically significant."}

### 2. Phase-Specific Means (Sackin normalized)
| Phase | Mean Sackin (norm) |
|-------|-------------------|
{chr(10).join(f"| {PHASE_SHORT[p]} | {phase_sackin[p]:.4f} |" for p in ALL_PHASES if p in phase_sackin)}

### 3. Frequency Band Patterns
- Most balanced band: **{most_balanced_band}** (Sackin_norm = {band_sackin[most_balanced_band]:.4f})
- Least balanced band: **{least_balanced_band}** (Sackin_norm = {band_sackin[least_balanced_band]:.4f})

Band-specific rest-to-task shifts (Sackin_norm, task - rest):
{chr(10).join(f"- {band}: {band_deltas[band]:+.4f}" for band in BANDS if band in band_deltas.index)}

### 4. Sackin-Colless Correlation
Pearson r = {r_sc:.3f}, p = {p_sc:.2e}
These two indices are {"strongly" if abs(r_sc) > 0.7 else "moderately" if abs(r_sc) > 0.4 else "weakly"} correlated,
{"confirming they capture similar aspects of tree imbalance." if abs(r_sc) > 0.7 else "suggesting they capture partially different aspects of tree geometry."}

### 5. Patient Consistency
Coefficient of variation across patient means: {cv_patients:.3f}
{"Low CV indicates consistent tree balance across patients." if cv_patients < 0.15 else "Moderate CV suggests some patient-level variability." if cv_patients < 0.3 else "High CV indicates substantial between-patient differences."}

Patient mean Sackin (norm):
{chr(10).join(f"- {p}: {patient_means[p]:.4f}" for p in PATIENTS if p in patient_means.index)}

## Interpretation
The tree balance indices capture the **geometry** of the multiscale hierarchy from LRG.
A more balanced tree means the brain's functional modules split into roughly equal-sized
sub-modules at each scale -- a distributed, egalitarian organization. An unbalanced tree
indicates a hub-dominated topology where one large cluster successively sheds small
peripheral groups.

{"**The shift toward more balanced trees during task suggests that cognitive engagement promotes a more distributed, parallel hierarchical organization.**" if task_sackin.mean() < rest_sackin.mean() and p_sackin < 0.05 else "**The shift toward less balanced trees during task suggests that cognitive engagement promotes a more centralized, hub-dominated hierarchical organization.**" if task_sackin.mean() > rest_sackin.mean() and p_sackin < 0.05 else "**No significant global shift in tree balance between rest and task was detected, though band-specific patterns may exist.**"}

## Files
- `tree_balance_results.csv` -- raw data
- `fig1_sackin_by_phase_band.png` -- Sackin index per phase/band
- `fig2_colless_by_phase_band.png` -- Colless index per phase/band
- `fig3_rest_vs_task_boxplot.png` -- Rest vs Task comparison
- `fig4_band_profiles.png` -- Balance across frequency bands by phase
- `fig5_patient_profiles.png` -- Per-patient profiles
- `fig6_sackin_vs_colless.png` -- Sackin-Colless correlation scatter
"""

findings_path = OUTPUT_DIR / "FINDINGS.md"
findings_path.write_text(findings)
print(f"\nSaved {findings_path}")
print("\nDone.")
