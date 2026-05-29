#!/usr/bin/env python3
"""Cophenetic correlation analysis: ultrametric fit of FC networks.

Measures how well the LRG-derived ultrametric (hierarchical) representation
preserves the original FC distance structure. A high cophenetic correlation
means the network is inherently hierarchical.

Outputs
-------
- 5 figures in data/figures/multiscale_investigation/cophenetic_fit/
- CSV results table
- FINDINGS.md summary
"""

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr, pearsonr

from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.config.const import (
    BRAIN_BANDS, BRAIN_BANDS_NAMES, PHASE_LABELS, BRAIN_BAND_TEX_DICT,
    PATIENTS_4PHASE,
)
from lrg_eegfc.config.paths import FIGURES_ROOT

import networkx as nx

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
OUTPUT_DIR = FIGURES_ROOT / "multiscale_investigation" / "cophenetic_fit"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]
BANDS = list(BRAIN_BANDS_NAMES)
PHASES = list(PHASE_LABELS)  # rest_pre, task_learn, task_test, rest_post
FC_METHOD = "msc"

# Patients with all 4 phases available
FULL_PATIENTS = PATIENTS_4PHASE

# Phase availability per patient (based on cache)
PATIENT_PHASES = {
    "Pat_02": list(PHASE_LABELS),
    "Pat_03": list(PHASE_LABELS),
    "Pat_05": list(PHASE_LABELS),
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "rest_post"],
    "Pat_08": list(PHASE_LABELS),
}

# Plotting style
PHASE_COLORS = {
    "rest_pre": "#2196F3",
    "task_learn": "#FF9800",
    "task_test": "#E91E63",
    "rest_post": "#4CAF50",
}
PHASE_LABELS_PRETTY = {
    "rest_pre": "Rest Pre",
    "task_learn": "Task Learn",
    "task_test": "Task Test",
    "rest_post": "Rest Post",
}
BAND_COLORS = {
    "delta": "#1f77b4",
    "theta": "#ff7f0e",
    "alpha": "#2ca02c",
    "beta": "#d62728",
    "low_gamma": "#9467bd",
    "high_gamma": "#8c564b",
}

# ---------------------------------------------------------------------------
# Computation
# ---------------------------------------------------------------------------
print("=" * 70)
print("COPHENETIC CORRELATION ANALYSIS")
print("=" * 70)

rows = []
skipped = []

for patient in ALL_PATIENTS:
    for band in BANDS:
        for phase in PATIENT_PHASES[patient]:
            # Load LRG result
            lrg = load_lrg_result(patient, phase, band, fc_method=FC_METHOD)
            if lrg is None:
                skipped.append((patient, band, phase, "no LRG cache"))
                continue

            # Load MSC matrix (dense, nperseg=4096)
            msc_mat = load_msc_matrix(patient, phase, band)
            if msc_mat is None:
                skipped.append((patient, band, phase, "no MSC cache"))
                continue

            # --- Extract giant component from MSC to match LRG node set ---
            G = nx.from_numpy_array(msc_mat)
            giant_nodes = sorted(max(nx.connected_components(G), key=len))
            n_giant = len(giant_nodes)

            # Verify consistency with LRG result
            if n_giant != lrg.n_nodes:
                skipped.append((patient, band, phase,
                                f"node mismatch: giant={n_giant} lrg={lrg.n_nodes}"))
                continue

            # Sub-matrix for giant component
            msc_giant = msc_mat[np.ix_(giant_nodes, giant_nodes)]

            # --- FC distances: 1 - MSC weight ---
            fc_dist = 1.0 - msc_giant
            np.fill_diagonal(fc_dist, 0.0)
            # Ensure symmetry and non-negative
            fc_dist = np.maximum(fc_dist, 0.0)
            fc_dist = (fc_dist + fc_dist.T) / 2.0
            fc_dist_condensed = squareform(fc_dist, checks=False)

            # --- Cophenetic distances from linkage ---
            # When Y is passed, cophenet returns (c, d); without Y just d
            coph_result = cophenet(lrg.linkage_matrix, fc_dist_condensed)
            coph_corr_direct = coph_result[0]  # cophenetic correlation
            coph_dists = coph_result[1]         # cophenetic distances

            # Verify same length
            if len(coph_dists) != len(fc_dist_condensed):
                skipped.append((patient, band, phase,
                                f"condensed length mismatch: coph={len(coph_dists)} fc={len(fc_dist_condensed)}"))
                continue

            # --- Cophenetic correlation (Pearson) ---
            r_coph_pearson, p_coph_pearson = pearsonr(fc_dist_condensed, coph_dists)

            # --- Spearman rank correlation between FC weights and ultrametric distances ---
            # ultrametric_matrix is in condensed form
            um_condensed = lrg.ultrametric_matrix
            fc_weights_condensed = squareform(msc_giant, checks=False)
            r_spearman, p_spearman = spearmanr(fc_weights_condensed, um_condensed)

            # --- Simple Pearson correlation between FC weights and ultrametric ---
            r_pearson_fc_um, p_pearson_fc_um = pearsonr(fc_weights_condensed, um_condensed)

            rows.append({
                "patient": patient,
                "band": band,
                "phase": phase,
                "n_nodes": n_giant,
                "n_total_nodes": msc_mat.shape[0],
                "cophenetic_r": r_coph_pearson,
                "cophenetic_p": p_coph_pearson,
                "spearman_fc_um": r_spearman,
                "spearman_p": p_spearman,
                "pearson_fc_um": r_pearson_fc_um,
                "pearson_fc_um_p": p_pearson_fc_um,
                "mean_fc_weight": float(np.mean(fc_weights_condensed)),
                "mean_um_dist": float(np.mean(um_condensed)),
                "mean_coph_dist": float(np.mean(coph_dists)),
            })

            print(f"  {patient} {band:12s} {phase:12s} | "
                  f"coph_r={r_coph_pearson:+.4f}  "
                  f"spearman={r_spearman:+.4f}  "
                  f"n={n_giant}")

if skipped:
    print(f"\nSkipped {len(skipped)} conditions:")
    for s in skipped:
        print(f"  {s}")

# ---------------------------------------------------------------------------
# Build DataFrame
# ---------------------------------------------------------------------------
df = pd.DataFrame(rows)
csv_path = OUTPUT_DIR / "cophenetic_results.csv"
df.to_csv(csv_path, index=False)
print(f"\nSaved results to {csv_path}")
print(f"Total conditions: {len(df)}")
print(f"\n--- Summary statistics ---")
print(df.groupby("band")["cophenetic_r"].describe().round(4))

# ---------------------------------------------------------------------------
# Figure 1: Cophenetic correlation across phases per band (bar chart)
#   averaged over FULL_PATIENTS only (all 4 phases)
# ---------------------------------------------------------------------------
print("\n--- Figure 1: Mean cophenetic correlation per band x phase ---")
df_full = df[df["patient"].isin(FULL_PATIENTS)]

fig1, ax1 = plt.subplots(figsize=(12, 5))

x_pos = np.arange(len(PHASES))
width = 0.12
offsets = np.linspace(-0.3, 0.3, len(BANDS))

for i, band in enumerate(BANDS):
    means = []
    sems = []
    for phase in PHASES:
        vals = df_full[(df_full["band"] == band) & (df_full["phase"] == phase)]["cophenetic_r"]
        means.append(vals.mean() if len(vals) > 0 else np.nan)
        sems.append(vals.std() / np.sqrt(len(vals)) if len(vals) > 1 else 0)
    ax1.bar(x_pos + offsets[i], means, width, yerr=sems,
            label=BRAIN_BAND_TEX_DICT[band], color=BAND_COLORS[band],
            edgecolor="white", linewidth=0.5, capsize=3)

ax1.set_xticks(x_pos)
ax1.set_xticklabels([PHASE_LABELS_PRETTY[p] for p in PHASES])
ax1.set_ylabel("Cophenetic Correlation (r)")
ax1.set_title("Ultrametric Fit of FC Networks Across Phases and Bands\n"
              f"(mean +/- SEM, N={len(FULL_PATIENTS)} patients with all phases)")
ax1.legend(ncol=len(BANDS), loc="upper center", bbox_to_anchor=(0.5, -0.12),
           frameon=False)
ax1.axhline(0, color="gray", linewidth=0.5, linestyle="--")
ax1.set_ylim(bottom=min(0, df_full["cophenetic_r"].min() - 0.05))
fig1.tight_layout()
fig1.savefig(OUTPUT_DIR / "fig1_cophenetic_by_phase_band.png", dpi=200,
             bbox_inches="tight")
plt.close(fig1)
print("  Saved fig1_cophenetic_by_phase_band.png")

# ---------------------------------------------------------------------------
# Figure 2: Per-patient cophenetic correlation profiles
#   One panel per band, one line per patient, x=phase
# ---------------------------------------------------------------------------
print("--- Figure 2: Per-patient profiles ---")
fig2, axes2 = plt.subplots(2, 3, figsize=(15, 9), sharey=True)
axes2 = axes2.flatten()

patient_markers = {"Pat_02": "o", "Pat_03": "s", "Pat_05": "D",
                   "Pat_06": "^", "Pat_07": "v", "Pat_08": "p"}

for idx, band in enumerate(BANDS):
    ax = axes2[idx]
    for patient in ALL_PATIENTS:
        phases_avail = PATIENT_PHASES[patient]
        vals = []
        phase_idx = []
        for pi, phase in enumerate(PHASES):
            if phase in phases_avail:
                row = df[(df["patient"] == patient) &
                         (df["band"] == band) &
                         (df["phase"] == phase)]
                if len(row) > 0:
                    vals.append(row["cophenetic_r"].values[0])
                    phase_idx.append(pi)
        if vals:
            ax.plot(phase_idx, vals, marker=patient_markers[patient],
                    label=patient, linewidth=1.5, markersize=6, alpha=0.8)

    ax.set_xticks(range(len(PHASES)))
    ax.set_xticklabels([PHASE_LABELS_PRETTY[p] for p in PHASES],
                       fontsize=8, rotation=20)
    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=13)
    ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
    if idx % 3 == 0:
        ax.set_ylabel("Cophenetic Correlation (r)")

axes2[0].legend(fontsize=7, loc="best")
fig2.suptitle("Cophenetic Correlation per Patient Across Phases", fontsize=14, y=1.01)
fig2.tight_layout()
fig2.savefig(OUTPUT_DIR / "fig2_cophenetic_per_patient.png", dpi=200,
             bbox_inches="tight")
plt.close(fig2)
print("  Saved fig2_cophenetic_per_patient.png")

# ---------------------------------------------------------------------------
# Figure 3: Scatter plot - cophenetic correlation vs simple FC matrix
#   correlation across phases
#   For each patient x band: compare the cophenetic_r with simple Pearson
#   correlation of flattened FC matrices between rest_pre and each other phase
# ---------------------------------------------------------------------------
print("--- Figure 3: Cophenetic vs FC matrix correlation ---")

scatter_rows = []
for patient in FULL_PATIENTS:
    for band in BANDS:
        # Get rest_pre MSC as reference
        msc_pre = load_msc_matrix(patient, "rest_pre", band)
        if msc_pre is None:
            continue
        pre_flat = squareform(msc_pre, checks=False)

        for phase in ["task_learn", "task_test", "rest_post"]:
            msc_other = load_msc_matrix(patient, phase, band)
            if msc_other is None:
                continue
            other_flat = squareform(msc_other, checks=False)

            # Simple FC matrix correlation (rest_pre vs phase)
            r_fc, _ = pearsonr(pre_flat, other_flat)

            # Get cophenetic correlation for this phase
            row_this = df[(df["patient"] == patient) &
                          (df["band"] == band) &
                          (df["phase"] == phase)]
            row_pre = df[(df["patient"] == patient) &
                         (df["band"] == band) &
                         (df["phase"] == "rest_pre")]
            if len(row_this) == 0 or len(row_pre) == 0:
                continue

            coph_this = row_this["cophenetic_r"].values[0]
            coph_pre = row_pre["cophenetic_r"].values[0]

            scatter_rows.append({
                "patient": patient,
                "band": band,
                "phase": phase,
                "fc_matrix_corr": r_fc,
                "cophenetic_r_phase": coph_this,
                "cophenetic_r_pre": coph_pre,
                "delta_cophenetic": coph_this - coph_pre,
            })

df_scatter = pd.DataFrame(scatter_rows)

fig3, ax3 = plt.subplots(figsize=(8, 7))
for band in BANDS:
    subset = df_scatter[df_scatter["band"] == band]
    ax3.scatter(subset["fc_matrix_corr"], subset["cophenetic_r_phase"],
                label=BRAIN_BAND_TEX_DICT[band], color=BAND_COLORS[band],
                s=60, alpha=0.7, edgecolors="white", linewidth=0.5)

# Add identity line
lims = [min(ax3.get_xlim()[0], ax3.get_ylim()[0]),
        max(ax3.get_xlim()[1], ax3.get_ylim()[1])]
ax3.plot(lims, lims, "k--", alpha=0.3, label="identity")

ax3.set_xlabel("FC Matrix Correlation (rest_pre vs Phase)")
ax3.set_ylabel("Cophenetic Correlation (r) in Phase")
ax3.set_title("FC Matrix Similarity vs Ultrametric Fit\n"
              "(each point = patient x band x phase)")
ax3.legend(loc="best", frameon=True, framealpha=0.9)
fig3.tight_layout()
fig3.savefig(OUTPUT_DIR / "fig3_cophenetic_vs_fc_corr.png", dpi=200,
             bbox_inches="tight")
plt.close(fig3)
print("  Saved fig3_cophenetic_vs_fc_corr.png")

# ---------------------------------------------------------------------------
# Figure 4: Heatmap of cophenetic correlation (patients x bands)
#   averaged over all available phases per patient
# ---------------------------------------------------------------------------
print("--- Figure 4: Heatmap patients x bands ---")
pivot = df.groupby(["patient", "band"])["cophenetic_r"].mean().unstack("band")
pivot = pivot.reindex(columns=BANDS)
pivot = pivot.reindex(ALL_PATIENTS)

fig4, ax4 = plt.subplots(figsize=(10, 5))
im = ax4.imshow(pivot.values, aspect="auto", cmap="RdYlGn", vmin=-0.2, vmax=0.8)
ax4.set_xticks(range(len(BANDS)))
ax4.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=12)
ax4.set_yticks(range(len(ALL_PATIENTS)))
ax4.set_yticklabels(ALL_PATIENTS)
ax4.set_title("Cophenetic Correlation by Patient and Band\n(averaged over phases)")

# Annotate cells
for i in range(len(ALL_PATIENTS)):
    for j in range(len(BANDS)):
        val = pivot.values[i, j]
        if np.isnan(val):
            continue
        color = "white" if val < 0.3 else "black"
        ax4.text(j, i, f"{val:.3f}", ha="center", va="center",
                 fontsize=9, color=color)

plt.colorbar(im, ax=ax4, label="Cophenetic Correlation (r)", shrink=0.8)
fig4.tight_layout()
fig4.savefig(OUTPUT_DIR / "fig4_cophenetic_heatmap.png", dpi=200,
             bbox_inches="tight")
plt.close(fig4)
print("  Saved fig4_cophenetic_heatmap.png")

# ---------------------------------------------------------------------------
# Figure 5: Delta(cophenetic) between rest and task phases
#   For each patient x band: delta = mean(task phases) - mean(rest phases)
# ---------------------------------------------------------------------------
print("--- Figure 5: Delta cophenetic rest vs task ---")

delta_rows = []
for patient in FULL_PATIENTS:
    for band in BANDS:
        rest_vals = df[(df["patient"] == patient) &
                       (df["band"] == band) &
                       (df["phase"].isin(["rest_pre", "rest_post"]))]["cophenetic_r"]
        task_vals = df[(df["patient"] == patient) &
                       (df["band"] == band) &
                       (df["phase"].isin(["task_learn", "task_test"]))]["cophenetic_r"]
        if len(rest_vals) > 0 and len(task_vals) > 0:
            delta_rows.append({
                "patient": patient,
                "band": band,
                "mean_rest": rest_vals.mean(),
                "mean_task": task_vals.mean(),
                "delta_task_rest": task_vals.mean() - rest_vals.mean(),
            })

df_delta = pd.DataFrame(delta_rows)

fig5, (ax5a, ax5b) = plt.subplots(1, 2, figsize=(14, 5))

# Panel A: grouped bar chart (mean delta per band)
band_deltas_mean = []
band_deltas_sem = []
for band in BANDS:
    vals = df_delta[df_delta["band"] == band]["delta_task_rest"]
    band_deltas_mean.append(vals.mean() if len(vals) > 0 else 0)
    band_deltas_sem.append(vals.std() / np.sqrt(len(vals)) if len(vals) > 1 else 0)

bars = ax5a.bar(range(len(BANDS)), band_deltas_mean, yerr=band_deltas_sem,
                color=[BAND_COLORS[b] for b in BANDS], edgecolor="white",
                linewidth=0.5, capsize=4)
ax5a.set_xticks(range(len(BANDS)))
ax5a.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=11)
ax5a.axhline(0, color="gray", linewidth=0.8, linestyle="--")
ax5a.set_ylabel(r"$\Delta$ Cophenetic (Task $-$ Rest)")
ax5a.set_title(f"Change in Ultrametric Fit: Task vs Rest\n"
               f"(mean +/- SEM, N={len(FULL_PATIENTS)})")

# Annotate bars
for bar, val in zip(bars, band_deltas_mean):
    ax5a.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
              f"{val:+.4f}", ha="center",
              va="bottom" if val >= 0 else "top", fontsize=8)

# Panel B: per-patient dots
for i, band in enumerate(BANDS):
    subset = df_delta[df_delta["band"] == band]
    jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(subset))
    ax5b.scatter(i + jitter, subset["delta_task_rest"],
                 color=BAND_COLORS[band], s=50, alpha=0.7,
                 edgecolors="white", linewidth=0.5, zorder=3)

ax5b.set_xticks(range(len(BANDS)))
ax5b.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=11)
ax5b.axhline(0, color="gray", linewidth=0.8, linestyle="--")
ax5b.set_ylabel(r"$\Delta$ Cophenetic (Task $-$ Rest)")
ax5b.set_title("Per-Patient Changes")

fig5.suptitle("Hierarchy Changes Between Rest and Task", fontsize=14, y=1.02)
fig5.tight_layout()
fig5.savefig(OUTPUT_DIR / "fig5_delta_cophenetic_rest_task.png", dpi=200,
             bbox_inches="tight")
plt.close(fig5)
print("  Saved fig5_delta_cophenetic_rest_task.png")

# ---------------------------------------------------------------------------
# FINDINGS
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("KEY FINDINGS")
print("=" * 70)

findings = []

# 1. Overall cophenetic correlation level
overall_mean = df["cophenetic_r"].mean()
overall_std = df["cophenetic_r"].std()
findings.append(f"Overall cophenetic correlation: {overall_mean:.4f} +/- {overall_std:.4f}")
print(f"\n1. {findings[-1]}")

# 2. Band-specific hierarchy
band_means = df.groupby("band")["cophenetic_r"].mean().reindex(BANDS)
most_hier = band_means.idxmax()
least_hier = band_means.idxmin()
findings.append(f"Most hierarchical band: {most_hier} (r={band_means[most_hier]:.4f})")
findings.append(f"Least hierarchical band: {least_hier} (r={band_means[least_hier]:.4f})")
print(f"2. {findings[-2]}")
print(f"   {findings[-1]}")
print(f"   Band ranking: {band_means.sort_values(ascending=False).to_dict()}")

# 3. Rest vs Task
if len(df_delta) > 0:
    for band in BANDS:
        vals = df_delta[df_delta["band"] == band]["delta_task_rest"]
        if len(vals) > 0:
            findings.append(f"  {band}: delta(task-rest) = {vals.mean():+.4f} +/- {vals.std():.4f}")
    print(f"\n3. Rest vs Task cophenetic changes:")
    for f in findings[-len(BANDS):]:
        print(f"   {f}")

# 4. Phase-specific patterns
print(f"\n4. Phase-specific cophenetic means (full patients):")
phase_means = df_full.groupby("phase")["cophenetic_r"].mean().reindex(PHASES)
for phase in PHASES:
    print(f"   {phase}: {phase_means[phase]:.4f}")

# 5. Spearman correlation (FC weights vs ultrametric distances)
spearman_mean = df["spearman_fc_um"].mean()
findings.append(f"Mean Spearman(FC weights, ultrametric distances): {spearman_mean:.4f}")
print(f"\n5. {findings[-1]}")
print(f"   (negative = higher FC weight <-> shorter ultrametric distance, as expected)")

# 6. Patient consistency
print(f"\n6. Per-patient mean cophenetic correlation:")
pat_means = df.groupby("patient")["cophenetic_r"].mean()
for pat in ALL_PATIENTS:
    if pat in pat_means.index:
        print(f"   {pat}: {pat_means[pat]:.4f}")

# 7. FC correlation vs cophenetic divergence
if len(df_scatter) > 0:
    r_div, p_div = pearsonr(df_scatter["fc_matrix_corr"],
                            df_scatter["cophenetic_r_phase"])
    findings.append(f"Correlation between FC matrix similarity and cophenetic fit: "
                    f"r={r_div:.4f}, p={p_div:.2e}")
    print(f"\n7. {findings[-1]}")

# ---------------------------------------------------------------------------
# Write FINDINGS.md
# ---------------------------------------------------------------------------
findings_md = OUTPUT_DIR / "FINDINGS.md"
with open(findings_md, "w") as f:
    f.write("# Cophenetic Correlation Analysis: Ultrametric Fit of FC Networks\n\n")
    f.write("## Overview\n\n")
    f.write("The cophenetic correlation measures how well the LRG-derived ultrametric\n")
    f.write("(hierarchical) representation preserves the original FC distance structure.\n")
    f.write("A high cophenetic correlation means the network is inherently hierarchical.\n\n")
    f.write(f"- **Patients analyzed**: {', '.join(ALL_PATIENTS)}\n")
    f.write(f"- **Full-phase patients** (all 4 phases): {', '.join(FULL_PATIENTS)}\n")
    f.write(f"- **Bands**: {', '.join(BANDS)}\n")
    f.write(f"- **Phases**: {', '.join(PHASES)}\n")
    f.write(f"- **FC method**: MSC (dense, nperseg=4096)\n")
    f.write(f"- **Total conditions**: {len(df)}\n\n")

    f.write("## Key Findings\n\n")
    for i, finding in enumerate(findings, 1):
        f.write(f"{i}. {finding}\n")

    f.write("\n## Band-Specific Hierarchy\n\n")
    f.write("| Band | Mean Cophenetic r | Std |\n")
    f.write("|------|-------------------|-----|\n")
    for band in BANDS:
        bdf = df[df["band"] == band]
        f.write(f"| {band} | {bdf['cophenetic_r'].mean():.4f} | {bdf['cophenetic_r'].std():.4f} |\n")

    f.write("\n## Phase-Specific Patterns (Full Patients)\n\n")
    f.write("| Phase | Mean Cophenetic r |\n")
    f.write("|-------|-------------------|\n")
    for phase in PHASES:
        val = phase_means[phase]
        f.write(f"| {PHASE_LABELS_PRETTY[phase]} | {val:.4f} |\n")

    f.write("\n## Rest vs Task Changes\n\n")
    f.write("| Band | Delta (Task - Rest) | Direction |\n")
    f.write("|------|--------------------|-----------|\n")
    for band in BANDS:
        vals = df_delta[df_delta["band"] == band]["delta_task_rest"]
        if len(vals) > 0:
            m = vals.mean()
            direction = "more hierarchical" if m > 0 else "less hierarchical"
            f.write(f"| {band} | {m:+.4f} | {direction} |\n")

    f.write("\n## Patient Consistency\n\n")
    f.write("| Patient | Mean Cophenetic r | Phases Available |\n")
    f.write("|---------|-------------------|------------------|\n")
    for pat in ALL_PATIENTS:
        if pat in pat_means.index:
            f.write(f"| {pat} | {pat_means[pat]:.4f} | {len(PATIENT_PHASES[pat])} |\n")

    f.write("\n## Figures\n\n")
    f.write("1. `fig1_cophenetic_by_phase_band.png` - Bar chart: cophenetic correlation across phases for each band\n")
    f.write("2. `fig2_cophenetic_per_patient.png` - Line plots: per-patient cophenetic profiles\n")
    f.write("3. `fig3_cophenetic_vs_fc_corr.png` - Scatter: cophenetic r vs FC matrix correlation\n")
    f.write("4. `fig4_cophenetic_heatmap.png` - Heatmap: patients x bands (phase-averaged)\n")
    f.write("5. `fig5_delta_cophenetic_rest_task.png` - Delta cophenetic between rest and task\n")

print(f"\nSaved {findings_md}")
print(f"\nDone! All outputs in: {OUTPUT_DIR}")
