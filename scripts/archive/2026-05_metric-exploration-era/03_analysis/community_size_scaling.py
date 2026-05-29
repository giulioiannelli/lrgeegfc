#!/usr/bin/env python3
"""Community Size Scaling Analysis.

At each hierarchical level k, the dendrogram produces communities of various
sizes.  The SIZE DISTRIBUTION of communities encodes whether the network has
core-periphery, distributed, or modular structure.

This script tracks how the size distribution evolves across scales and compares
across conditions (rest vs task phases).
"""

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import csv
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.cluster.hierarchy import fcluster

from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    PHASE_LABELS,
    PATIENTS_4PHASE as _PATIENTS_4PHASE,
)
from lrg_eegfc.config.paths import FIGURES_ROOT

# ── Configuration ────────────────────────────────────────────────────────────
FC_METHOD = "msc"
PATIENTS_4PHASE = _PATIENTS_4PHASE
PATIENTS_TASK_ONLY = ["Pat_06", "Pat_07"]
ALL_PATIENTS = PATIENTS_4PHASE + PATIENTS_TASK_ONLY
BANDS = list(BRAIN_BANDS.keys())
PHASES = list(PHASE_LABELS)  # rest_pre, task_learn, task_test, rest_post
K_RANGE = np.arange(2, 31)  # 2..30
OUTPUT_DIR = FIGURES_ROOT / "multiscale_investigation" / "size_scaling"

# Phase colours (consistent with project style)
PHASE_COLORS = {
    "rest_pre": "#1f77b4",
    "task_learn": "#ff7f0e",
    "task_test": "#2ca02c",
    "rest_post": "#d62728",
}
PHASE_LABELS_TEX = {
    "rest_pre": "Rest Pre",
    "task_learn": "Task Learn",
    "task_test": "Task Test",
    "rest_post": "Rest Post",
}

REST_PHASES = ["rest_pre", "rest_post"]
TASK_PHASES = ["task_learn", "task_test"]


# ── Metric functions ─────────────────────────────────────────────────────────
def gini(sizes):
    """Gini coefficient of community sizes."""
    sizes = np.sort(np.asarray(sizes, dtype=float))
    n = len(sizes)
    if n == 0 or np.sum(sizes) == 0:
        return 0.0
    index = np.arange(1, n + 1)
    return float((2 * np.sum(index * sizes) - (n + 1) * np.sum(sizes)) / (n * np.sum(sizes)))


def size_ratio(sizes):
    """Ratio of largest to second-largest community size."""
    s = sorted(sizes, reverse=True)
    if len(s) < 2 or s[1] == 0:
        return np.inf
    return s[0] / s[1]


def cv(sizes):
    """Coefficient of variation of community sizes."""
    sizes = np.asarray(sizes, dtype=float)
    if len(sizes) == 0 or np.mean(sizes) == 0:
        return 0.0
    return float(np.std(sizes) / np.mean(sizes))


def n_meaningful(sizes, n_nodes, threshold=0.05):
    """Number of communities with size > threshold * n_nodes."""
    cutoff = threshold * n_nodes
    return int(np.sum(np.array(sizes) > cutoff))


def compute_scaling_profile(linkage, n_nodes, k_range):
    """Compute all metrics across k values for one LRG result."""
    gini_vals = np.full(len(k_range), np.nan)
    ratio_vals = np.full(len(k_range), np.nan)
    cv_vals = np.full(len(k_range), np.nan)
    n_mean_vals = np.full(len(k_range), np.nan)

    for i, k in enumerate(k_range):
        if k > n_nodes:
            break
        labels = fcluster(linkage, k, criterion="maxclust")
        _, counts = np.unique(labels, return_counts=True)
        sizes = sorted(counts, reverse=True)

        gini_vals[i] = gini(sizes)
        ratio_vals[i] = size_ratio(sizes)
        cv_vals[i] = cv(sizes)
        n_mean_vals[i] = n_meaningful(sizes, n_nodes)

    return {
        "gini": gini_vals,
        "size_ratio": ratio_vals,
        "cv": cv_vals,
        "n_meaningful": n_mean_vals,
    }


def get_optimal_k(linkage, optimal_threshold):
    """Get the number of clusters at the LRG optimal threshold."""
    labels = fcluster(linkage, optimal_threshold, criterion="distance")
    return len(np.unique(labels))


# ── Data loading ─────────────────────────────────────────────────────────────
print("=" * 70)
print("Community Size Scaling Analysis")
print("=" * 70)

# Determine valid phases per patient
def get_phases(patient):
    if patient in PATIENTS_TASK_ONLY:
        return ["task_learn", "task_test"]
    return list(PHASES)

# Load all LRG results and compute profiles
results = {}  # (patient, band, phase) -> profile dict
metadata = {}  # (patient, band, phase) -> {n_nodes, optimal_threshold, optimal_k}
n_loaded = 0
n_missing = 0

for patient in ALL_PATIENTS:
    for band in BANDS:
        for phase in get_phases(patient):
            lrg = load_lrg_result(patient, phase, band, FC_METHOD)
            if lrg is None:
                n_missing += 1
                continue
            n_loaded += 1
            profile = compute_scaling_profile(lrg.linkage_matrix, lrg.n_nodes, K_RANGE)
            results[(patient, band, phase)] = profile

            opt_k = get_optimal_k(lrg.linkage_matrix, lrg.optimal_threshold)
            metadata[(patient, band, phase)] = {
                "n_nodes": lrg.n_nodes,
                "optimal_threshold": lrg.optimal_threshold,
                "optimal_k": opt_k,
            }

print(f"\nLoaded {n_loaded} LRG results, {n_missing} missing")

# ── Output directory ─────────────────────────────────────────────────────────
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: Gini(k) profiles for all 4 phases, per band, Pat_02
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[Fig 1] Gini(k) profiles — Pat_02, all bands")

fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharex=True, sharey=True)
axes = axes.ravel()

for idx, band in enumerate(BANDS):
    ax = axes[idx]
    for phase in PHASES:
        key = ("Pat_02", band, phase)
        if key not in results:
            continue
        ax.plot(K_RANGE, results[key]["gini"], label=PHASE_LABELS_TEX[phase],
                color=PHASE_COLORS[phase], linewidth=1.8, alpha=0.85)
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=14)
    ax.set_xlabel("k (number of communities)")
    ax.set_ylabel("Gini coefficient")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)
    if idx == 0:
        ax.legend(fontsize=9, loc="upper right")

fig.suptitle("Community Size Inequality (Gini) Across Scales — Pat_02", fontsize=16, y=1.01)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig1_gini_profiles_pat02.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTPUT_DIR / 'fig1_gini_profiles_pat02.png'}")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: Size ratio(k) profiles — Pat_02
# ═══════════════════════════════════════════════════════════════════════════════
print("[Fig 2] Size ratio(k) profiles — Pat_02, all bands")

fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharex=True)
axes = axes.ravel()

for idx, band in enumerate(BANDS):
    ax = axes[idx]
    for phase in PHASES:
        key = ("Pat_02", band, phase)
        if key not in results:
            continue
        vals = results[key]["size_ratio"].copy()
        vals = np.clip(vals, 0, 20)  # clip for visualization
        ax.plot(K_RANGE, vals, label=PHASE_LABELS_TEX[phase],
                color=PHASE_COLORS[phase], linewidth=1.8, alpha=0.85)
    ax.axhline(y=1.0, color="gray", linestyle="--", alpha=0.5, label="Equal split")
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=14)
    ax.set_xlabel("k (number of communities)")
    ax.set_ylabel("Size ratio (largest / 2nd largest)")
    ax.set_ylim(0, 15)
    ax.grid(True, alpha=0.3)
    if idx == 0:
        ax.legend(fontsize=8, loc="upper right")

fig.suptitle("Community Size Ratio Across Scales — Pat_02", fontsize=16, y=1.01)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig2_size_ratio_pat02.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTPUT_DIR / 'fig2_size_ratio_pat02.png'}")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 3: Grand average Gini(k) — rest vs task
# ═══════════════════════════════════════════════════════════════════════════════
print("[Fig 3] Grand average Gini(k) — rest vs task")

# Collect per-condition profiles (only 4-phase patients)
rest_profiles = []
task_profiles = []
for patient in PATIENTS_4PHASE:
    for band in BANDS:
        for phase in REST_PHASES:
            key = (patient, band, phase)
            if key in results:
                rest_profiles.append(results[key]["gini"])
        for phase in TASK_PHASES:
            key = (patient, band, phase)
            if key in results:
                task_profiles.append(results[key]["gini"])

rest_arr = np.array(rest_profiles)
task_arr = np.array(task_profiles)

fig, ax = plt.subplots(1, 1, figsize=(10, 6))

# Mean and SEM
with warnings.catch_warnings():
    warnings.simplefilter("ignore", RuntimeWarning)
    rest_mean = np.nanmean(rest_arr, axis=0)
    rest_sem = np.nanstd(rest_arr, axis=0) / np.sqrt(np.sum(~np.isnan(rest_arr), axis=0))
    task_mean = np.nanmean(task_arr, axis=0)
    task_sem = np.nanstd(task_arr, axis=0) / np.sqrt(np.sum(~np.isnan(task_arr), axis=0))

ax.plot(K_RANGE, rest_mean, color="#1f77b4", linewidth=2.5, label="Rest (rest_pre + rest_post)")
ax.fill_between(K_RANGE, rest_mean - rest_sem, rest_mean + rest_sem, color="#1f77b4", alpha=0.2)
ax.plot(K_RANGE, task_mean, color="#ff7f0e", linewidth=2.5, label="Task (task_learn + task_test)")
ax.fill_between(K_RANGE, task_mean - task_sem, task_mean + task_sem, color="#ff7f0e", alpha=0.2)

ax.set_xlabel("k (number of communities)", fontsize=13)
ax.set_ylabel("Gini coefficient", fontsize=13)
ax.set_title("Grand Average: Community Size Inequality — Rest vs Task", fontsize=15)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
ax.set_ylim(-0.05, 1.05)

fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig3_grand_avg_gini_rest_vs_task.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTPUT_DIR / 'fig3_grand_avg_gini_rest_vs_task.png'}")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 4: N meaningful communities at LRG optimal threshold — heatmap
# ═══════════════════════════════════════════════════════════════════════════════
print("[Fig 4] N meaningful communities at optimal threshold — heatmap")

# For 4-phase patients, build a (bands x phases) heatmap averaged across patients
n_mean_at_opt = np.full((len(BANDS), len(PHASES)), np.nan)

for bi, band in enumerate(BANDS):
    for pi, phase in enumerate(PHASES):
        vals = []
        for patient in PATIENTS_4PHASE:
            key = (patient, band, phase)
            if key in metadata and key in results:
                opt_k = metadata[key]["optimal_k"]
                n_nodes = metadata[key]["n_nodes"]
                # Get community sizes at optimal k
                lrg = load_lrg_result(patient, phase, band, FC_METHOD)
                if lrg is not None:
                    labels = fcluster(lrg.linkage_matrix, opt_k, criterion="maxclust")
                    _, counts = np.unique(labels, return_counts=True)
                    vals.append(n_meaningful(counts, n_nodes))
        if vals:
            n_mean_at_opt[bi, pi] = np.mean(vals)

fig, ax = plt.subplots(1, 1, figsize=(8, 6))
im = ax.imshow(n_mean_at_opt, aspect="auto", cmap="YlOrRd", interpolation="nearest")
ax.set_xticks(range(len(PHASES)))
ax.set_xticklabels([PHASE_LABELS_TEX[p] for p in PHASES], fontsize=11)
ax.set_yticks(range(len(BANDS)))
ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=11)

# Annotate cells
for bi in range(len(BANDS)):
    for pi in range(len(PHASES)):
        val = n_mean_at_opt[bi, pi]
        if not np.isnan(val):
            ax.text(pi, bi, f"{val:.1f}", ha="center", va="center", fontsize=11,
                    color="white" if val > np.nanmax(n_mean_at_opt) * 0.6 else "black")

cb = fig.colorbar(im, ax=ax, shrink=0.8)
cb.set_label("Mean # meaningful communities", fontsize=11)
ax.set_title("Meaningful Communities at LRG Optimal Threshold\n(averaged over 4-phase patients)", fontsize=14)

fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig4_n_meaningful_heatmap.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTPUT_DIR / 'fig4_n_meaningful_heatmap.png'}")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 5: Per-patient Gini profiles — alpha band
# ═══════════════════════════════════════════════════════════════════════════════
print("[Fig 5] Per-patient Gini profiles — alpha band")

fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey=True)
axes = axes.ravel()

for idx, patient in enumerate(ALL_PATIENTS):
    ax = axes[idx]
    phases_for_pat = get_phases(patient)
    for phase in phases_for_pat:
        key = (patient, "alpha", phase)
        if key not in results:
            continue
        ax.plot(K_RANGE, results[key]["gini"], label=PHASE_LABELS_TEX[phase],
                color=PHASE_COLORS[phase], linewidth=1.8, alpha=0.85)
    ax.set_title(patient, fontsize=13)
    ax.set_xlabel("k (number of communities)")
    ax.set_ylabel("Gini coefficient")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)
    if idx == 0:
        ax.legend(fontsize=9, loc="upper right")

fig.suptitle(r"Community Size Inequality (Gini) — $\alpha$ band, all patients", fontsize=16, y=1.01)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig5_gini_per_patient_alpha.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTPUT_DIR / 'fig5_gini_per_patient_alpha.png'}")


# ═══════════════════════════════════════════════════════════════════════════════
# FIGURE 6: CV(k) profiles — Pat_02
# ═══════════════════════════════════════════════════════════════════════════════
print("[Fig 6] CV(k) profiles — Pat_02, all bands")

fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharex=True, sharey=True)
axes = axes.ravel()

for idx, band in enumerate(BANDS):
    ax = axes[idx]
    for phase in PHASES:
        key = ("Pat_02", band, phase)
        if key not in results:
            continue
        ax.plot(K_RANGE, results[key]["cv"], label=PHASE_LABELS_TEX[phase],
                color=PHASE_COLORS[phase], linewidth=1.8, alpha=0.85)
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=14)
    ax.set_xlabel("k (number of communities)")
    ax.set_ylabel("Coefficient of Variation")
    ax.grid(True, alpha=0.3)
    if idx == 0:
        ax.legend(fontsize=9, loc="upper right")

fig.suptitle("Community Size CV Across Scales — Pat_02", fontsize=16, y=1.01)
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig6_cv_profiles_pat02.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTPUT_DIR / 'fig6_cv_profiles_pat02.png'}")


# ═══════════════════════════════════════════════════════════════════════════════
# CSV RESULTS
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[CSV] Writing results...")

csv_path = OUTPUT_DIR / "community_size_scaling_results.csv"
with open(csv_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "patient", "band", "phase", "k",
        "gini", "size_ratio", "cv", "n_meaningful",
        "n_nodes", "optimal_threshold", "optimal_k",
    ])
    for (patient, band, phase), profile in sorted(results.items()):
        meta = metadata[(patient, band, phase)]
        for i, k in enumerate(K_RANGE):
            writer.writerow([
                patient, band, phase, k,
                f"{profile['gini'][i]:.4f}" if not np.isnan(profile['gini'][i]) else "",
                f"{profile['size_ratio'][i]:.4f}" if not np.isnan(profile['size_ratio'][i]) else "",
                f"{profile['cv'][i]:.4f}" if not np.isnan(profile['cv'][i]) else "",
                f"{profile['n_meaningful'][i]:.0f}" if not np.isnan(profile['n_meaningful'][i]) else "",
                meta["n_nodes"],
                f"{meta['optimal_threshold']:.4f}",
                meta["optimal_k"],
            ])

print(f"  Saved: {csv_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# FINDINGS
# ═══════════════════════════════════════════════════════════════════════════════
print("\n[Analysis] Computing key findings...")

findings_lines = []
findings_lines.append("# Community Size Scaling Analysis — Findings\n")
findings_lines.append(f"**FC method**: {FC_METHOD}")
findings_lines.append(f"**Patients (4-phase)**: {', '.join(PATIENTS_4PHASE)}")
findings_lines.append(f"**Patients (task-only)**: {', '.join(PATIENTS_TASK_ONLY)}")
findings_lines.append(f"**k range**: {K_RANGE[0]}–{K_RANGE[-1]}")
findings_lines.append(f"**Total LRG results loaded**: {n_loaded}\n")

# 1. Rest vs task Gini difference
findings_lines.append("## 1. Rest vs Task: Gini Coefficient\n")
for k_idx, k_val in enumerate(K_RANGE):
    if k_val in [2, 5, 10, 15, 20, 30]:
        r_val = rest_mean[k_idx] if not np.isnan(rest_mean[k_idx]) else None
        t_val = task_mean[k_idx] if not np.isnan(task_mean[k_idx]) else None
        if r_val is not None and t_val is not None:
            diff = r_val - t_val
            direction = "rest > task (more centralized at rest)" if diff > 0 else "task > rest (more centralized during task)"
            findings_lines.append(f"- k={k_val}: Rest={r_val:.3f}, Task={t_val:.3f}, diff={diff:+.3f} ({direction})")

# 2. Equalization scale: at which k does size_ratio drop to ~1?
findings_lines.append("\n## 2. Equalization Scale (k where size ratio first drops below 2)\n")
for patient in PATIENTS_4PHASE:
    for band in ["alpha", "beta"]:
        for phase in PHASES:
            key = (patient, band, phase)
            if key not in results:
                continue
            ratios = results[key]["size_ratio"]
            eq_k = None
            for i, k_val in enumerate(K_RANGE):
                if not np.isnan(ratios[i]) and ratios[i] < 2.0:
                    eq_k = k_val
                    break
            if eq_k is not None:
                findings_lines.append(f"- {patient} {band} {PHASE_LABELS_TEX[phase]}: k={eq_k}")
            else:
                findings_lines.append(f"- {patient} {band} {PHASE_LABELS_TEX[phase]}: never (ratio stays >= 2)")

# 3. Meaningful communities at optimal threshold
findings_lines.append("\n## 3. Meaningful Communities at LRG Optimal Threshold\n")
findings_lines.append("| Band | Rest Pre | Task Learn | Task Test | Rest Post |")
findings_lines.append("|------|----------|------------|-----------|-----------|")
for bi, band in enumerate(BANDS):
    row = f"| {BRAIN_BAND_TEX_DICT[band]} |"
    for pi, phase in enumerate(PHASES):
        val = n_mean_at_opt[bi, pi]
        if np.isnan(val):
            row += " — |"
        else:
            row += f" {val:.1f} |"
    findings_lines.append(row)

# 4. Band-specific inequality patterns
findings_lines.append("\n## 4. Band-Specific Inequality at k=5 (Grand Average)\n")
for band in BANDS:
    vals = []
    for patient in PATIENTS_4PHASE:
        for phase in PHASES:
            key = (patient, band, phase)
            if key in results:
                k5_idx = 3  # k=5 is index 3 (K_RANGE starts at 2)
                g = results[key]["gini"][k5_idx]
                if not np.isnan(g):
                    vals.append(g)
    if vals:
        findings_lines.append(f"- {band}: Gini(k=5) = {np.mean(vals):.3f} +/- {np.std(vals):.3f} (n={len(vals)})")

# 5. Patient consistency
findings_lines.append("\n## 5. Cross-Patient Consistency (alpha band, k=5)\n")
for phase in PHASES:
    vals = []
    for patient in PATIENTS_4PHASE:
        key = (patient, "alpha", phase)
        if key in results:
            g = results[key]["gini"][3]
            if not np.isnan(g):
                vals.append(g)
    if vals:
        findings_lines.append(
            f"- {PHASE_LABELS_TEX[phase]}: Gini = {np.mean(vals):.3f} +/- {np.std(vals):.3f} "
            f"(range {np.min(vals):.3f}–{np.max(vals):.3f}, n={len(vals)})"
        )

# 6. Optimal k distribution
findings_lines.append("\n## 6. Optimal k from LRG Threshold\n")
for band in BANDS:
    opt_ks = []
    for patient in PATIENTS_4PHASE:
        for phase in PHASES:
            key = (patient, band, phase)
            if key in metadata:
                opt_ks.append(metadata[key]["optimal_k"])
    if opt_ks:
        findings_lines.append(
            f"- {band}: optimal k = {np.mean(opt_ks):.1f} +/- {np.std(opt_ks):.1f} "
            f"(range {np.min(opt_ks)}–{np.max(opt_ks)})"
        )

findings_text = "\n".join(findings_lines) + "\n"
findings_path = OUTPUT_DIR / "FINDINGS.md"
with open(findings_path, "w") as f:
    f.write(findings_text)
print(f"  Saved: {findings_path}")

print("\n" + "=" * 70)
print("Analysis complete. All outputs saved to:")
print(f"  {OUTPUT_DIR}")
print("=" * 70)
