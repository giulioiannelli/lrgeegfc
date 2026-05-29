#!/usr/bin/env python3
"""Partition Entropy Profile Analysis.

At each dendrogram level k (number of clusters), compute:
- H(k): partition entropy = -sum_i (n_i/N) * log(n_i/N)
- dH/dk: entropy gradient (where biggest splits happen)
- largest_frac(k): fraction of nodes in the biggest community
- n_eff(k) = exp(H(k)): effective number of communities

Generates 7 figures + CSV + FINDINGS.md
"""

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.cluster.hierarchy import fcluster
from scipy.stats import pearsonr
from pathlib import Path

from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE
from lrg_eegfc.config.paths import FIGURES_ROOT

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
FC_METHOD = "msc"
ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]
FOUR_PHASE_PATIENTS = PATIENTS_4PHASE
TWO_PHASE_PATIENTS = ["Pat_06", "Pat_07"]
BANDS = list(BRAIN_BANDS.keys())
PHASES = list(PHASE_LABELS)  # rest_pre, task_learn, task_test, rest_post
K_RANGE = np.arange(2, 41)  # k = 2..40
OUTPUT_DIR = FIGURES_ROOT / "multiscale_investigation" / "partition_entropy"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Phase colours (consistent across project)
PHASE_COLORS = {
    "rest_pre": "#1f77b4",
    "task_learn": "#ff7f0e",
    "task_test": "#2ca02c",
    "rest_post": "#d62728",
}
PHASE_LABELS_PRETTY = {
    "rest_pre": "Rest Pre",
    "task_learn": "Task Learn",
    "task_test": "Task Test",
    "rest_post": "Rest Post",
}

# Condition groupings for analysis
REST_PHASES = ["rest_pre", "rest_post"]
TASK_PHASES = ["task_learn", "task_test"]


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------
def compute_partition_entropy(linkage_matrix, n_nodes, k_range):
    """Compute partition entropy profile H(k) for a range of k values."""
    H = np.zeros(len(k_range))
    largest_frac = np.zeros(len(k_range))
    n_communities = np.zeros(len(k_range), dtype=int)

    for i, k in enumerate(k_range):
        labels = fcluster(linkage_matrix, t=k, criterion="maxclust")
        _, counts = np.unique(labels, return_counts=True)
        N = len(labels)
        probs = counts / N
        H[i] = -np.sum(probs * np.log(probs))
        largest_frac[i] = np.max(counts) / N
        n_communities[i] = len(counts)

    # Effective number of communities
    n_eff = np.exp(H)
    # Entropy gradient
    dH_dk = np.gradient(H, k_range)

    return {
        "H": H,
        "dH_dk": dH_dk,
        "largest_frac": largest_frac,
        "n_eff": n_eff,
        "n_communities": n_communities,
    }


# ---------------------------------------------------------------------------
# Load all data
# ---------------------------------------------------------------------------
print("Loading LRG results and computing partition entropy profiles...")
results = {}  # (patient, band, phase) -> dict
rows = []

for patient in ALL_PATIENTS:
    phases_for_patient = PHASES if patient in FOUR_PHASE_PATIENTS else ["task_learn", "task_test"]
    for band in BANDS:
        for phase in phases_for_patient:
            lrg = load_lrg_result(patient, phase, band, FC_METHOD)
            if lrg is None:
                print(f"  SKIP {patient} {phase} {band} -- no cached LRG result")
                continue
            prof = compute_partition_entropy(lrg.linkage_matrix, lrg.n_nodes, K_RANGE)
            results[(patient, band, phase)] = prof

            # Find scale of max dH/dk
            idx_max_dh = np.argmax(prof["dH_dk"])
            k_max_dh = K_RANGE[idx_max_dh]
            # Fragmentation scale: first k where largest_frac < 0.5
            below_half = np.where(prof["largest_frac"] < 0.5)[0]
            k_frag = int(K_RANGE[below_half[0]]) if len(below_half) > 0 else -1

            rows.append({
                "patient": patient,
                "band": band,
                "phase": phase,
                "n_nodes": lrg.n_nodes,
                "k_max_dH": int(k_max_dh),
                "max_dH_value": float(prof["dH_dk"][idx_max_dh]),
                "k_fragmentation": k_frag,
                "H_at_k10": float(prof["H"][K_RANGE == 10][0]) if 10 in K_RANGE else np.nan,
                "H_at_k20": float(prof["H"][K_RANGE == 20][0]) if 20 in K_RANGE else np.nan,
                "n_eff_at_k10": float(prof["n_eff"][K_RANGE == 10][0]) if 10 in K_RANGE else np.nan,
                "n_eff_at_k20": float(prof["n_eff"][K_RANGE == 20][0]) if 20 in K_RANGE else np.nan,
            })

df = pd.DataFrame(rows)
df.to_csv(OUTPUT_DIR / "partition_entropy_profiles.csv", index=False)
print(f"  Loaded {len(results)} profiles, saved CSV.")

# ---------------------------------------------------------------------------
# Figure 1: H(k) profiles for Pat_02, one subplot per band
# ---------------------------------------------------------------------------
print("Generating Figure 1: H(k) profiles for Pat_02...")
fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharey=True, sharex=True)
axes = axes.ravel()
for i, band in enumerate(BANDS):
    ax = axes[i]
    for phase in PHASES:
        key = ("Pat_02", band, phase)
        if key not in results:
            continue
        prof = results[key]
        ax.plot(K_RANGE, prof["H"], color=PHASE_COLORS[phase],
                label=PHASE_LABELS_PRETTY[phase], linewidth=1.8)
    ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]} ({band})", fontsize=13)
    ax.set_xlabel("k (number of clusters)")
    ax.set_ylabel("H(k)")
    ax.grid(True, alpha=0.3)
    if i == 0:
        ax.legend(fontsize=9)
fig.suptitle("Partition Entropy H(k) -- Pat_02", fontsize=15, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUTPUT_DIR / "fig1_Hk_Pat02.png", dpi=180, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 2: dH/dk profiles for Pat_02
# ---------------------------------------------------------------------------
print("Generating Figure 2: dH/dk profiles for Pat_02...")
fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharey=True, sharex=True)
axes = axes.ravel()
for i, band in enumerate(BANDS):
    ax = axes[i]
    for phase in PHASES:
        key = ("Pat_02", band, phase)
        if key not in results:
            continue
        prof = results[key]
        ax.plot(K_RANGE, prof["dH_dk"], color=PHASE_COLORS[phase],
                label=PHASE_LABELS_PRETTY[phase], linewidth=1.8)
    ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]} ({band})", fontsize=13)
    ax.set_xlabel("k (number of clusters)")
    ax.set_ylabel("dH/dk")
    ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
    ax.grid(True, alpha=0.3)
    if i == 0:
        ax.legend(fontsize=9)
fig.suptitle("Entropy Gradient dH/dk -- Pat_02", fontsize=15, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUTPUT_DIR / "fig2_dHdk_Pat02.png", dpi=180, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 3: Largest community fraction, averaged across 4-phase patients
# ---------------------------------------------------------------------------
print("Generating Figure 3: Largest community fraction (patient average)...")
fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharey=True, sharex=True)
axes = axes.ravel()
for i, band in enumerate(BANDS):
    ax = axes[i]
    for phase in PHASES:
        fracs = []
        for patient in FOUR_PHASE_PATIENTS:
            key = (patient, band, phase)
            if key in results:
                fracs.append(results[key]["largest_frac"])
        if not fracs:
            continue
        mean_frac = np.mean(fracs, axis=0)
        sem_frac = np.std(fracs, axis=0) / np.sqrt(len(fracs))
        ax.plot(K_RANGE, mean_frac, color=PHASE_COLORS[phase],
                label=PHASE_LABELS_PRETTY[phase], linewidth=1.8)
        ax.fill_between(K_RANGE, mean_frac - sem_frac, mean_frac + sem_frac,
                        color=PHASE_COLORS[phase], alpha=0.15)
    ax.axhline(y=0.5, color="gray", linestyle="--", alpha=0.5, label="50%")
    ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]} ({band})", fontsize=13)
    ax.set_xlabel("k (number of clusters)")
    ax.set_ylabel("Largest community fraction")
    ax.grid(True, alpha=0.3)
    if i == 0:
        ax.legend(fontsize=8)
fig.suptitle("Largest Community Fraction vs k (mean over 4-phase patients)",
             fontsize=15, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUTPUT_DIR / "fig3_largest_frac_avg.png", dpi=180, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 4: Heatmap of k at max dH/dk across phases/bands
# ---------------------------------------------------------------------------
print("Generating Figure 4: k(max dH/dk) heatmap...")
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# 4a: Average across 4-phase patients
k_max_matrix = np.full((len(BANDS), len(PHASES)), np.nan)
for ib, band in enumerate(BANDS):
    for ip, phase in enumerate(PHASES):
        vals = []
        for patient in FOUR_PHASE_PATIENTS:
            key = (patient, band, phase)
            if key in results:
                prof = results[key]
                vals.append(K_RANGE[np.argmax(prof["dH_dk"])])
        if vals:
            k_max_matrix[ib, ip] = np.mean(vals)

ax = axes[0]
im = ax.imshow(k_max_matrix, aspect="auto", cmap="YlOrRd")
ax.set_xticks(range(len(PHASES)))
ax.set_xticklabels([PHASE_LABELS_PRETTY[p] for p in PHASES], rotation=30, ha="right")
ax.set_yticks(range(len(BANDS)))
ax.set_yticklabels([f"{BRAIN_BAND_TEX_DICT[b]} ({b})" for b in BANDS])
for ib in range(len(BANDS)):
    for ip in range(len(PHASES)):
        if not np.isnan(k_max_matrix[ib, ip]):
            ax.text(ip, ib, f"{k_max_matrix[ib, ip]:.1f}",
                    ha="center", va="center", fontsize=10, fontweight="bold")
ax.set_title("Scale of max dH/dk\n(mean over 4-phase patients)", fontsize=13)
plt.colorbar(im, ax=ax, label="k at max dH/dk")

# 4b: max dH/dk VALUE
dh_max_matrix = np.full((len(BANDS), len(PHASES)), np.nan)
for ib, band in enumerate(BANDS):
    for ip, phase in enumerate(PHASES):
        vals = []
        for patient in FOUR_PHASE_PATIENTS:
            key = (patient, band, phase)
            if key in results:
                prof = results[key]
                vals.append(np.max(prof["dH_dk"]))
        if vals:
            dh_max_matrix[ib, ip] = np.mean(vals)

ax = axes[1]
im2 = ax.imshow(dh_max_matrix, aspect="auto", cmap="viridis")
ax.set_xticks(range(len(PHASES)))
ax.set_xticklabels([PHASE_LABELS_PRETTY[p] for p in PHASES], rotation=30, ha="right")
ax.set_yticks(range(len(BANDS)))
ax.set_yticklabels([f"{BRAIN_BAND_TEX_DICT[b]} ({b})" for b in BANDS])
for ib in range(len(BANDS)):
    for ip in range(len(PHASES)):
        if not np.isnan(dh_max_matrix[ib, ip]):
            ax.text(ip, ib, f"{dh_max_matrix[ib, ip]:.3f}",
                    ha="center", va="center", fontsize=9, fontweight="bold",
                    color="white" if dh_max_matrix[ib, ip] > np.nanmedian(dh_max_matrix) else "black")
plt.colorbar(im2, ax=ax, label="max dH/dk value")
ax.set_title("Peak entropy gradient magnitude\n(mean over 4-phase patients)", fontsize=13)

fig.suptitle("Most Informative Split: Scale and Magnitude", fontsize=15, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(OUTPUT_DIR / "fig4_max_dHdk_heatmap.png", dpi=180, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 5: n_eff(k), rest vs task grand average
# ---------------------------------------------------------------------------
print("Generating Figure 5: Effective number of communities (rest vs task)...")
fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharey=False, sharex=True)
axes = axes.ravel()
for i, band in enumerate(BANDS):
    ax = axes[i]
    for condition_name, condition_phases, color, ls in [
        ("Rest", REST_PHASES, "#1f77b4", "-"),
        ("Task", TASK_PHASES, "#ff7f0e", "--"),
    ]:
        all_neff = []
        for patient in FOUR_PHASE_PATIENTS:
            for phase in condition_phases:
                key = (patient, band, phase)
                if key in results:
                    all_neff.append(results[key]["n_eff"])
        if not all_neff:
            continue
        mean_neff = np.mean(all_neff, axis=0)
        sem_neff = np.std(all_neff, axis=0) / np.sqrt(len(all_neff))
        ax.plot(K_RANGE, mean_neff, color=color, linestyle=ls,
                label=condition_name, linewidth=2)
        ax.fill_between(K_RANGE, mean_neff - sem_neff, mean_neff + sem_neff,
                        color=color, alpha=0.15)
    # Reference line: n_eff = k (perfect equipartition)
    ax.plot(K_RANGE, K_RANGE, color="gray", linestyle=":", alpha=0.5, label="n_eff = k")
    ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]} ({band})", fontsize=13)
    ax.set_xlabel("k (number of clusters)")
    ax.set_ylabel("n_eff = exp(H)")
    ax.grid(True, alpha=0.3)
    if i == 0:
        ax.legend(fontsize=9)
fig.suptitle("Effective Number of Communities: Rest vs Task (grand average)",
             fontsize=15, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUTPUT_DIR / "fig5_neff_rest_vs_task.png", dpi=180, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 6: Per-patient H(k) for alpha band
# ---------------------------------------------------------------------------
print("Generating Figure 6: Per-patient H(k) for alpha band...")
fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharey=True, sharex=True)
axes = axes.ravel()
patient_colors = plt.cm.Set2(np.linspace(0, 1, len(ALL_PATIENTS)))
for ip, phase in enumerate(PHASES):
    ax = axes[ip]
    for j, patient in enumerate(ALL_PATIENTS):
        key = (patient, "alpha", phase)
        if key not in results:
            continue
        ax.plot(K_RANGE, results[key]["H"], color=patient_colors[j],
                label=patient, linewidth=1.5)
    ax.set_title(f"{PHASE_LABELS_PRETTY[phase]}", fontsize=13)
    ax.set_xlabel("k")
    ax.set_ylabel("H(k)")
    ax.grid(True, alpha=0.3)
    if ip == 0:
        ax.legend(fontsize=8, ncol=2)

# Use remaining subplots for variance across patients
ax = axes[4]
for phase in PHASES:
    Hs = []
    for patient in FOUR_PHASE_PATIENTS:
        key = (patient, "alpha", phase)
        if key in results:
            Hs.append(results[key]["H"])
    if Hs:
        std_H = np.std(Hs, axis=0)
        ax.plot(K_RANGE, std_H, color=PHASE_COLORS[phase],
                label=PHASE_LABELS_PRETTY[phase], linewidth=1.5)
ax.set_title("Std dev of H(k) across patients", fontsize=13)
ax.set_xlabel("k")
ax.set_ylabel("std(H)")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8)

# Coefficient of variation
ax = axes[5]
for phase in PHASES:
    Hs = []
    for patient in FOUR_PHASE_PATIENTS:
        key = (patient, "alpha", phase)
        if key in results:
            Hs.append(results[key]["H"])
    if Hs:
        mean_H = np.mean(Hs, axis=0)
        std_H = np.std(Hs, axis=0)
        cv = std_H / (mean_H + 1e-12)
        ax.plot(K_RANGE, cv, color=PHASE_COLORS[phase],
                label=PHASE_LABELS_PRETTY[phase], linewidth=1.5)
ax.set_title("CV of H(k) across patients", fontsize=13)
ax.set_xlabel("k")
ax.set_ylabel("CV = std/mean")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8)

fig.suptitle(r"Partition Entropy H(k) -- $\alpha$ band, per patient",
             fontsize=15, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(OUTPUT_DIR / "fig6_Hk_alpha_per_patient.png", dpi=180, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Figure 7: Cross-phase Pearson correlation of H(k) profiles
# ---------------------------------------------------------------------------
print("Generating Figure 7: Cross-phase H(k) correlation...")
# For each band, compute within-rest, within-task, cross-condition correlations
# Average across 4-phase patients

# Panel A: Per-band correlation matrices (mean across patients)
fig, axes = plt.subplots(2, 3, figsize=(16, 11))
axes = axes.ravel()
for ib, band in enumerate(BANDS):
    ax = axes[ib]
    corr_mat = np.full((len(PHASES), len(PHASES)), np.nan)
    for i1, p1 in enumerate(PHASES):
        for i2, p2 in enumerate(PHASES):
            rs = []
            for patient in FOUR_PHASE_PATIENTS:
                k1 = (patient, band, p1)
                k2 = (patient, band, p2)
                if k1 in results and k2 in results:
                    r, _ = pearsonr(results[k1]["H"], results[k2]["H"])
                    rs.append(r)
            if rs:
                corr_mat[i1, i2] = np.mean(rs)
    im = ax.imshow(corr_mat, vmin=0.9, vmax=1.0, cmap="RdYlGn", aspect="equal")
    ax.set_xticks(range(len(PHASES)))
    ax.set_xticklabels([PHASE_LABELS_PRETTY[p] for p in PHASES], rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(PHASES)))
    ax.set_yticklabels([PHASE_LABELS_PRETTY[p] for p in PHASES], fontsize=9)
    for i1 in range(len(PHASES)):
        for i2 in range(len(PHASES)):
            if not np.isnan(corr_mat[i1, i2]):
                ax.text(i2, i1, f"{corr_mat[i1, i2]:.3f}",
                        ha="center", va="center", fontsize=9, fontweight="bold")
    ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]} ({band})", fontsize=13)
    plt.colorbar(im, ax=ax, shrink=0.8)

fig.suptitle("Cross-phase Pearson correlation of H(k) profiles\n(mean over 4-phase patients)",
             fontsize=15, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(OUTPUT_DIR / "fig7_cross_phase_Hk_correlation.png", dpi=180, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Analysis summary
# ---------------------------------------------------------------------------
print("\n" + "=" * 70)
print("PARTITION ENTROPY PROFILE ANALYSIS -- KEY FINDINGS")
print("=" * 70)

findings = []

# 1. Profile shape consistency
findings.append("## 1. Profile Shape Consistency Across Phases")
for band in BANDS:
    corrs_within_rest = []
    corrs_within_task = []
    corrs_cross = []
    for patient in FOUR_PHASE_PATIENTS:
        k_rsPre = (patient, band, "rest_pre")
        k_rsPost = (patient, band, "rest_post")
        k_tL = (patient, band, "task_learn")
        k_tT = (patient, band, "task_test")
        if all(k in results for k in [k_rsPre, k_rsPost, k_tL, k_tT]):
            r_rest, _ = pearsonr(results[k_rsPre]["H"], results[k_rsPost]["H"])
            r_task, _ = pearsonr(results[k_tL]["H"], results[k_tT]["H"])
            corrs_within_rest.append(r_rest)
            corrs_within_task.append(r_task)
            for kr in [k_rsPre, k_rsPost]:
                for kt in [k_tL, k_tT]:
                    r_cross, _ = pearsonr(results[kr]["H"], results[kt]["H"])
                    corrs_cross.append(r_cross)
    if corrs_within_rest:
        msg = (f"  {band:12s}: within-rest r={np.mean(corrs_within_rest):.4f}, "
               f"within-task r={np.mean(corrs_within_task):.4f}, "
               f"cross-condition r={np.mean(corrs_cross):.4f}")
        print(msg)
        findings.append(msg)

# 2. Scale of most informative split
findings.append("\n## 2. Scale of Most Informative Split (k at max dH/dk)")
for band in BANDS:
    for phase in PHASES:
        vals = df[(df["band"] == band) & (df["phase"] == phase) & df["patient"].isin(FOUR_PHASE_PATIENTS)]["k_max_dH"]
        if len(vals) > 0:
            msg = f"  {band:12s} {phase:12s}: k_max_dH = {vals.mean():.1f} +/- {vals.std():.1f}"
            print(msg)
            findings.append(msg)

# 3. Fragmentation scale
findings.append("\n## 3. Fragmentation Scale (k where largest community < 50%)")
for band in BANDS:
    for phase in PHASES:
        vals = df[(df["band"] == band) & (df["phase"] == phase) & df["patient"].isin(FOUR_PHASE_PATIENTS)]["k_fragmentation"]
        vals = vals[vals > 0]
        if len(vals) > 0:
            msg = f"  {band:12s} {phase:12s}: k_frag = {vals.mean():.1f} +/- {vals.std():.1f}"
            print(msg)
            findings.append(msg)

# 4. Effective communities at reference scales
findings.append("\n## 4. Effective Number of Communities at k=10 and k=20")
for band in BANDS:
    rest_10 = df[(df["band"] == band) & df["phase"].isin(REST_PHASES) & df["patient"].isin(FOUR_PHASE_PATIENTS)]["n_eff_at_k10"]
    task_10 = df[(df["band"] == band) & df["phase"].isin(TASK_PHASES) & df["patient"].isin(FOUR_PHASE_PATIENTS)]["n_eff_at_k10"]
    rest_20 = df[(df["band"] == band) & df["phase"].isin(REST_PHASES) & df["patient"].isin(FOUR_PHASE_PATIENTS)]["n_eff_at_k20"]
    task_20 = df[(df["band"] == band) & df["phase"].isin(TASK_PHASES) & df["patient"].isin(FOUR_PHASE_PATIENTS)]["n_eff_at_k20"]
    if len(rest_10) > 0:
        msg = (f"  {band:12s}: k=10 rest={rest_10.mean():.2f} task={task_10.mean():.2f} | "
               f"k=20 rest={rest_20.mean():.2f} task={task_20.mean():.2f}")
        print(msg)
        findings.append(msg)

# 5. Band-specific non-linearity
findings.append("\n## 5. Non-linearity of Entropy Profiles (deviation from linear)")
for band in BANDS:
    devs = []
    for patient in FOUR_PHASE_PATIENTS:
        for phase in PHASES:
            key = (patient, band, phase)
            if key in results:
                H = results[key]["H"]
                # Linear reference
                H_linear = np.linspace(H[0], H[-1], len(H))
                rmse = np.sqrt(np.mean((H - H_linear) ** 2))
                devs.append(rmse)
    if devs:
        msg = f"  {band:12s}: RMSE from linear = {np.mean(devs):.4f} +/- {np.std(devs):.4f}"
        print(msg)
        findings.append(msg)

# 6. Rest-task difference in entropy at k=10
findings.append("\n## 6. Rest vs Task Difference in H(k=10)")
for band in BANDS:
    rest_h10 = df[(df["band"] == band) & df["phase"].isin(REST_PHASES) & df["patient"].isin(FOUR_PHASE_PATIENTS)]["H_at_k10"]
    task_h10 = df[(df["band"] == band) & df["phase"].isin(TASK_PHASES) & df["patient"].isin(FOUR_PHASE_PATIENTS)]["H_at_k10"]
    if len(rest_h10) > 0 and len(task_h10) > 0:
        diff = task_h10.mean() - rest_h10.mean()
        msg = f"  {band:12s}: H_task - H_rest = {diff:+.4f}  (task {'higher' if diff > 0 else 'lower'})"
        print(msg)
        findings.append(msg)

# ---------------------------------------------------------------------------
# Write FINDINGS.md
# ---------------------------------------------------------------------------
findings_text = f"""# Partition Entropy Profile Analysis

**Date:** 2026-03-10
**FC method:** {FC_METHOD}
**Patients:** {', '.join(ALL_PATIENTS)} (4-phase: {', '.join(FOUR_PHASE_PATIENTS)})
**k range:** {K_RANGE[0]}--{K_RANGE[-1]}

## Method

For each dendrogram (LRG linkage), at each scale k (number of clusters):
- Partition entropy: H(k) = -sum_i (n_i/N) * log(n_i/N)
- Entropy gradient: dH/dk via finite differences
- Largest community fraction: max(n_i) / N
- Effective number of communities: n_eff(k) = exp(H(k))

""" + "\n".join(findings) + """

## Figures

1. `fig1_Hk_Pat02.png` -- H(k) profiles for all 4 phases, per band (Pat_02)
2. `fig2_dHdk_Pat02.png` -- dH/dk entropy gradient profiles (Pat_02)
3. `fig3_largest_frac_avg.png` -- Largest community fraction vs k (patient average)
4. `fig4_max_dHdk_heatmap.png` -- Scale and magnitude of most informative split (heatmap)
5. `fig5_neff_rest_vs_task.png` -- Effective number of communities, rest vs task
6. `fig6_Hk_alpha_per_patient.png` -- Per-patient H(k) for alpha band
7. `fig7_cross_phase_Hk_correlation.png` -- Cross-phase correlation of H(k) profiles
"""

with open(OUTPUT_DIR / "FINDINGS.md", "w") as f:
    f.write(findings_text)

print(f"\nAll outputs saved to {OUTPUT_DIR}/")
print("Done.")
