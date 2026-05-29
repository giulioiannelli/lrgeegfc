#!/usr/bin/env python3
"""Cross-Phase Multiscale Community Flow analysis.

For each pair of phases, at each hierarchical level k (number of clusters),
compute the contingency table between two partitions and measure mixing entropy,
NMI, and Variation of Information (VI) across scales. The shape of the
VI-vs-scale profile is a multiscale reorganization fingerprint.

VI is a true metric on the space of partitions with information-theoretic
interpretation: VI(P,Q) = H(P|Q) + H(Q|P).  VI=0 means identical partitions;
higher values indicate more different partitions.
"""

import warnings
from itertools import combinations
from pathlib import Path
from collections import defaultdict

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import fcluster
from sklearn.metrics import normalized_mutual_info_score

from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    PHASE_LABELS,
    PATIENTS_4PHASE,
)
from lrg_eegfc.config.paths import FIGURES_ROOT

warnings.filterwarnings("ignore", category=FutureWarning)

# ── configuration ──────────────────────────────────────────────────────
FC_METHOD = "msc"
K_VALUES = np.arange(2, 31)  # scales 2..30
OUTPUT_DIR = FIGURES_ROOT / "multiscale_investigation" / "community_flow"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

FOUR_PHASE_PATIENTS = PATIENTS_4PHASE
ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]
PHASES = list(PHASE_LABELS)  # rest_pre, task_learn, task_test, rest_post
BANDS = list(BRAIN_BANDS.keys())

# Phase-pair categories
WITHIN_PAIRS = [("rest_pre", "rest_post"), ("task_learn", "task_test")]
CROSS_PAIRS = [
    ("rest_pre", "task_learn"),
    ("rest_pre", "task_test"),
    ("task_learn", "rest_post"),
    ("task_test", "rest_post"),
]
ALL_PAIRS = WITHIN_PAIRS + CROSS_PAIRS

PHASE_TEX = {
    "rest_pre": "rest_pre",
    "task_learn": "task_learn",
    "task_test": "task_test",
    "rest_post": "rest_post",
}

# ── helper functions ───────────────────────────────────────────────────

def conditional_entropy(labels_a, labels_b):
    """Compute H(B|A) = how uncertain B is given A.

    H(B|A) = - sum_a P(a) sum_b P(b|a) log2 P(b|a)
    """
    n = len(labels_a)
    unique_a = np.unique(labels_a)
    h = 0.0
    for a in unique_a:
        mask = labels_a == a
        p_a = mask.sum() / n
        b_given_a = labels_b[mask]
        unique_b, counts_b = np.unique(b_given_a, return_counts=True)
        probs_b = counts_b / counts_b.sum()
        h -= p_a * np.sum(probs_b * np.log2(probs_b + 1e-30))
    return h


from lrg_eegfc.utils.metrics import compute_vi  # canonical implementation


def compute_scale_profiles(linkage_a, linkage_b, n_nodes, k_values):
    """Compute mixing entropy, NMI, VI profiles across scales k."""
    mixing_entropy = np.full(len(k_values), np.nan)
    nmi = np.full(len(k_values), np.nan)
    vi = np.full(len(k_values), np.nan)

    for i, k in enumerate(k_values):
        if k > n_nodes:
            continue
        labels_a = fcluster(linkage_a, t=k, criterion="maxclust")
        labels_b = fcluster(linkage_b, t=k, criterion="maxclust")

        mixing_entropy[i] = conditional_entropy(labels_a, labels_b)
        nmi[i] = normalized_mutual_info_score(labels_a, labels_b)
        vi[i] = compute_vi(labels_a, labels_b)

    return mixing_entropy, nmi, vi


# ── main computation ───────────────────────────────────────────────────

print("=" * 70)
print("Cross-Phase Multiscale Community Flow Analysis")
print("  (using Variation of Information as primary partition distance)")
print("=" * 70)

# Storage: results[patient][band][(phaseA, phaseB)] = {entropy, nmi, vi}
results = defaultdict(lambda: defaultdict(dict))
skipped = []

for patient in ALL_PATIENTS:
    for band in BANDS:
        # Load all available phases for this patient/band
        lrg_data = {}
        for phase in PHASES:
            r = load_lrg_result(patient, phase, band, FC_METHOD)
            if r is not None:
                lrg_data[phase] = r

        available_phases = list(lrg_data.keys())

        for pA, pB in ALL_PAIRS:
            if pA not in lrg_data or pB not in lrg_data:
                continue

            rA = lrg_data[pA]
            rB = lrg_data[pB]

            # Use minimum n_nodes between the two
            n_nodes = min(rA.n_nodes, rB.n_nodes)

            me, nmi_vals, vi_vals = compute_scale_profiles(
                rA.linkage_matrix, rB.linkage_matrix, n_nodes, K_VALUES
            )
            results[patient][band][(pA, pB)] = {
                "entropy": me,
                "nmi": nmi_vals,
                "vi": vi_vals,
                "n_nodes": n_nodes,
            }

        print(f"  {patient} / {band}: {len(available_phases)} phases, "
              f"{sum(1 for p in ALL_PAIRS if p in results[patient][band])} pairs")

print(f"\nComputation complete.")

# ── Figure 1: Mixing entropy profiles per band ────────────────────────
print("\nGenerating Figure 1: Mixing entropy profiles...")

fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey=True)
axes = axes.flatten()

pair_colors = plt.cm.tab10(np.linspace(0, 1, len(ALL_PAIRS)))

for idx, band in enumerate(BANDS):
    ax = axes[idx]

    for pair_idx, (pA, pB) in enumerate(ALL_PAIRS):
        # Average across 4-phase patients
        profiles = []
        for patient in FOUR_PHASE_PATIENTS:
            if (pA, pB) in results[patient][band]:
                profiles.append(results[patient][band][(pA, pB)]["entropy"])

        if not profiles:
            continue

        mean_profile = np.nanmean(profiles, axis=0)
        sem_profile = np.nanstd(profiles, axis=0) / np.sqrt(len(profiles))

        is_within = (pA, pB) in WITHIN_PAIRS
        ls = "-" if is_within else "--"
        lw = 2.5 if is_within else 1.5
        label = f"{PHASE_TEX[pA]}-{PHASE_TEX[pB]}"

        ax.plot(K_VALUES, mean_profile, ls, color=pair_colors[pair_idx],
                lw=lw, label=label)
        ax.fill_between(K_VALUES, mean_profile - sem_profile,
                        mean_profile + sem_profile,
                        alpha=0.15, color=pair_colors[pair_idx])

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=14)
    ax.set_xlabel("Number of clusters (k)")
    ax.set_ylabel("Mixing entropy H(B|A) [bits]")
    ax.grid(True, alpha=0.3)

axes[0].legend(fontsize=7, loc="upper left", ncol=2)
fig.suptitle("Multiscale Mixing Entropy Profiles\n(4-phase patients, solid=within, dashed=cross)",
             fontsize=15, fontweight="bold")
plt.tight_layout()
fig.savefig(OUTPUT_DIR / "fig1_mixing_entropy_profiles.png", dpi=200,
            bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTPUT_DIR / 'fig1_mixing_entropy_profiles.png'}")

# ── Figure 2: NMI profiles per band ───────────────────────────────────
print("Generating Figure 2: NMI profiles...")

fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey=True)
axes = axes.flatten()

for idx, band in enumerate(BANDS):
    ax = axes[idx]

    for pair_idx, (pA, pB) in enumerate(ALL_PAIRS):
        profiles = []
        for patient in FOUR_PHASE_PATIENTS:
            if (pA, pB) in results[patient][band]:
                profiles.append(results[patient][band][(pA, pB)]["nmi"])

        if not profiles:
            continue

        mean_profile = np.nanmean(profiles, axis=0)
        sem_profile = np.nanstd(profiles, axis=0) / np.sqrt(len(profiles))

        is_within = (pA, pB) in WITHIN_PAIRS
        ls = "-" if is_within else "--"
        lw = 2.5 if is_within else 1.5
        label = f"{PHASE_TEX[pA]}-{PHASE_TEX[pB]}"

        ax.plot(K_VALUES, mean_profile, ls, color=pair_colors[pair_idx],
                lw=lw, label=label)
        ax.fill_between(K_VALUES, mean_profile - sem_profile,
                        mean_profile + sem_profile,
                        alpha=0.15, color=pair_colors[pair_idx])

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=14)
    ax.set_xlabel("Number of clusters (k)")
    ax.set_ylabel("NMI")
    ax.grid(True, alpha=0.3)

axes[0].legend(fontsize=7, loc="upper right", ncol=2)
fig.suptitle("Multiscale NMI Profiles\n(4-phase patients, solid=within, dashed=cross)",
             fontsize=15, fontweight="bold")
plt.tight_layout()
fig.savefig(OUTPUT_DIR / "fig2_nmi_profiles.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTPUT_DIR / 'fig2_nmi_profiles.png'}")

# ── Figure 3: VI profiles per band ────────────────────────────────────
print("Generating Figure 3: VI profiles per band...")

fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True, sharey=True)
axes = axes.flatten()

for idx, band in enumerate(BANDS):
    ax = axes[idx]

    for pair_idx, (pA, pB) in enumerate(ALL_PAIRS):
        profiles = []
        for patient in FOUR_PHASE_PATIENTS:
            if (pA, pB) in results[patient][band]:
                profiles.append(results[patient][band][(pA, pB)]["vi"])

        if not profiles:
            continue

        mean_profile = np.nanmean(profiles, axis=0)
        sem_profile = np.nanstd(profiles, axis=0) / np.sqrt(len(profiles))

        is_within = (pA, pB) in WITHIN_PAIRS
        ls = "-" if is_within else "--"
        lw = 2.5 if is_within else 1.5
        label = f"{PHASE_TEX[pA]}-{PHASE_TEX[pB]}"

        ax.plot(K_VALUES, mean_profile, ls, color=pair_colors[pair_idx],
                lw=lw, label=label)
        ax.fill_between(K_VALUES, mean_profile - sem_profile,
                        mean_profile + sem_profile,
                        alpha=0.15, color=pair_colors[pair_idx])

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=14)
    ax.set_xlabel("Number of clusters (k)")
    ax.set_ylabel("VI (nats)")
    ax.grid(True, alpha=0.3)

axes[0].legend(fontsize=7, loc="upper left", ncol=2)
fig.suptitle("Multiscale Variation of Information Profiles\n"
             "(4-phase patients, solid=within, dashed=cross)",
             fontsize=15, fontweight="bold")
plt.tight_layout()
fig.savefig(OUTPUT_DIR / "fig3_vi_profiles.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTPUT_DIR / 'fig3_vi_profiles.png'}")

# ── Figure 4: Grand average within vs cross (VI as primary) ───────────
print("Generating Figure 4: Within vs cross condition (VI)...")

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Panel A: VI (primary)
ax = axes[0]
for condition_label, pair_list, color, ls in [
    ("Within-condition", WITHIN_PAIRS, "#2196F3", "-"),
    ("Cross-condition", CROSS_PAIRS, "#F44336", "--"),
]:
    all_profiles = []
    for patient in FOUR_PHASE_PATIENTS:
        for band in BANDS:
            for pA, pB in pair_list:
                if (pA, pB) in results[patient][band]:
                    all_profiles.append(results[patient][band][(pA, pB)]["vi"])

    if not all_profiles:
        continue

    mean_profile = np.nanmean(all_profiles, axis=0)
    sem_profile = np.nanstd(all_profiles, axis=0) / np.sqrt(len(all_profiles))

    ax.plot(K_VALUES, mean_profile, ls, color=color, lw=2.5, label=condition_label)
    ax.fill_between(K_VALUES, mean_profile - sem_profile,
                    mean_profile + sem_profile, alpha=0.2, color=color)

ax.set_xlabel("Number of clusters (k)", fontsize=13)
ax.set_ylabel("VI (nats)", fontsize=13)
ax.set_title("Variation of Information", fontsize=14, fontweight="bold")
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)

# Panel B: Mixing entropy
ax = axes[1]
for condition_label, pair_list, color, ls in [
    ("Within-condition", WITHIN_PAIRS, "#2196F3", "-"),
    ("Cross-condition", CROSS_PAIRS, "#F44336", "--"),
]:
    all_profiles = []
    for patient in FOUR_PHASE_PATIENTS:
        for band in BANDS:
            for pA, pB in pair_list:
                if (pA, pB) in results[patient][band]:
                    all_profiles.append(results[patient][band][(pA, pB)]["entropy"])

    if not all_profiles:
        continue

    mean_profile = np.nanmean(all_profiles, axis=0)
    sem_profile = np.nanstd(all_profiles, axis=0) / np.sqrt(len(all_profiles))

    ax.plot(K_VALUES, mean_profile, ls, color=color, lw=2.5, label=condition_label)
    ax.fill_between(K_VALUES, mean_profile - sem_profile,
                    mean_profile + sem_profile, alpha=0.2, color=color)

ax.set_xlabel("Number of clusters (k)", fontsize=13)
ax.set_ylabel("Mixing entropy H(B|A) [bits]", fontsize=13)
ax.set_title("Mixing Entropy", fontsize=14, fontweight="bold")
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)

fig.suptitle("Grand Average: Within vs Cross-Condition\n"
             "(averaged across bands and 4-phase patients)",
             fontsize=15, fontweight="bold")
plt.tight_layout()
fig.savefig(OUTPUT_DIR / "fig4_within_vs_cross_grand_average.png", dpi=200,
            bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTPUT_DIR / 'fig4_within_vs_cross_grand_average.png'}")

# ── Figure 5: VI-based reorganization specificity ratio per band ──────
print("Generating Figure 5: VI reorganization specificity ratio...")

fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True)
axes = axes.flatten()

ratio_data = {}  # band -> ratio profile

for idx, band in enumerate(BANDS):
    ax = axes[idx]

    within_profiles = []
    cross_profiles = []

    for patient in FOUR_PHASE_PATIENTS:
        for pA, pB in WITHIN_PAIRS:
            if (pA, pB) in results[patient][band]:
                within_profiles.append(results[patient][band][(pA, pB)]["vi"])
        for pA, pB in CROSS_PAIRS:
            if (pA, pB) in results[patient][band]:
                cross_profiles.append(results[patient][band][(pA, pB)]["vi"])

    if within_profiles and cross_profiles:
        mean_within = np.nanmean(within_profiles, axis=0)
        mean_cross = np.nanmean(cross_profiles, axis=0)

        # Ratio: cross / within (>1 means cross-condition partitions are
        # more different from each other, as expected)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(
                mean_within > 0.01,
                mean_cross / mean_within,
                np.nan,
            )

        ratio_data[band] = ratio

        ax.plot(K_VALUES, ratio, "k-", lw=2)
        ax.axhline(1.0, color="gray", ls=":", lw=1)

        # Highlight peak
        valid_mask = ~np.isnan(ratio)
        if valid_mask.any():
            peak_idx = np.nanargmax(ratio)
            ax.axvline(K_VALUES[peak_idx], color="red", ls="--", alpha=0.7)
            ax.scatter([K_VALUES[peak_idx]], [ratio[peak_idx]], color="red",
                       s=80, zorder=5)
            ax.annotate(f"k={K_VALUES[peak_idx]}\nratio={ratio[peak_idx]:.2f}",
                        xy=(K_VALUES[peak_idx], ratio[peak_idx]),
                        xytext=(5, 10), textcoords="offset points", fontsize=9,
                        color="red")

    ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=14)
    ax.set_xlabel("Number of clusters (k)")
    ax.set_ylabel("Cross/Within VI ratio")
    ax.grid(True, alpha=0.3)

fig.suptitle("Reorganization Specificity Ratio (Cross / Within VI)\n"
             "Peak = scale of maximum task-sensitivity",
             fontsize=15, fontweight="bold")
plt.tight_layout()
fig.savefig(OUTPUT_DIR / "fig5_reorganization_specificity_ratio.png", dpi=200,
            bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTPUT_DIR / 'fig5_reorganization_specificity_ratio.png'}")

# ── Figure 6: Heatmap of peak reorganization scale ────────────────────
print("Generating Figure 6: Peak reorganization heatmap...")

peak_k_matrix = np.full((len(FOUR_PHASE_PATIENTS), len(BANDS)), np.nan)

for pi, patient in enumerate(FOUR_PHASE_PATIENTS):
    for bi, band in enumerate(BANDS):
        within_profiles = []
        cross_profiles = []

        for pA, pB in WITHIN_PAIRS:
            if (pA, pB) in results[patient][band]:
                within_profiles.append(results[patient][band][(pA, pB)]["vi"])
        for pA, pB in CROSS_PAIRS:
            if (pA, pB) in results[patient][band]:
                cross_profiles.append(results[patient][band][(pA, pB)]["vi"])

        if within_profiles and cross_profiles:
            mean_within = np.nanmean(within_profiles, axis=0)
            mean_cross = np.nanmean(cross_profiles, axis=0)
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(
                    mean_within > 0.01,
                    mean_cross / mean_within,
                    np.nan,
                )
            valid = ~np.isnan(ratio)
            if valid.any():
                peak_k_matrix[pi, bi] = K_VALUES[np.nanargmax(ratio)]

fig, ax = plt.subplots(figsize=(10, 5))
im = ax.imshow(peak_k_matrix, aspect="auto", cmap="viridis",
               vmin=K_VALUES[0], vmax=K_VALUES[-1])
ax.set_xticks(range(len(BANDS)))
ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=12)
ax.set_yticks(range(len(FOUR_PHASE_PATIENTS)))
ax.set_yticklabels(FOUR_PHASE_PATIENTS, fontsize=12)
ax.set_xlabel("Frequency Band", fontsize=13)
ax.set_ylabel("Patient", fontsize=13)
ax.set_title("Scale of Maximum Reorganization Specificity (VI-based)\n"
             "(k where cross/within VI ratio peaks)",
             fontsize=14, fontweight="bold")

# Annotate cells
for pi in range(len(FOUR_PHASE_PATIENTS)):
    for bi in range(len(BANDS)):
        val = peak_k_matrix[pi, bi]
        if not np.isnan(val):
            ax.text(bi, pi, f"{int(val)}", ha="center", va="center",
                    fontsize=12, fontweight="bold", color="white")

plt.colorbar(im, ax=ax, label="k (number of clusters)")
plt.tight_layout()
fig.savefig(OUTPUT_DIR / "fig6_peak_reorganization_heatmap.png", dpi=200,
            bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTPUT_DIR / 'fig6_peak_reorganization_heatmap.png'}")

# ── Figure 7: Per-patient alpha VI profiles ───────────────────────────
print("Generating Figure 7: Per-patient alpha band VI profiles...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

for pi, patient in enumerate(FOUR_PHASE_PATIENTS):
    ax = axes[pi]
    band = "alpha"

    for pair_idx, (pA, pB) in enumerate(ALL_PAIRS):
        if (pA, pB) not in results[patient][band]:
            continue

        vi_prof = results[patient][band][(pA, pB)]["vi"]
        is_within = (pA, pB) in WITHIN_PAIRS
        ls = "-" if is_within else "--"
        lw = 2.5 if is_within else 1.5
        label = f"{PHASE_TEX[pA]}-{PHASE_TEX[pB]}"

        ax.plot(K_VALUES, vi_prof, ls, color=pair_colors[pair_idx], lw=lw,
                label=label)

    ax.set_title(f"{patient} — {BRAIN_BAND_TEX_DICT['alpha']}", fontsize=13)
    ax.set_xlabel("Number of clusters (k)")
    ax.set_ylabel("VI (nats)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc="upper left", ncol=2)

fig.suptitle("Per-Patient VI Profiles (Alpha Band)\n"
             "solid=within-condition, dashed=cross-condition",
             fontsize=15, fontweight="bold")
plt.tight_layout()
fig.savefig(OUTPUT_DIR / "fig7_per_patient_alpha_vi_profiles.png", dpi=200,
            bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {OUTPUT_DIR / 'fig7_per_patient_alpha_vi_profiles.png'}")

# ── CSV Export ─────────────────────────────────────────────────────────
print("\nExporting CSV results...")

rows = []
for patient in ALL_PATIENTS:
    for band in BANDS:
        for pA, pB in ALL_PAIRS:
            if (pA, pB) not in results[patient][band]:
                continue
            d = results[patient][band][(pA, pB)]
            condition = "within" if (pA, pB) in WITHIN_PAIRS else "cross"
            for i, k in enumerate(K_VALUES):
                rows.append({
                    "patient": patient,
                    "band": band,
                    "phase_A": pA,
                    "phase_B": pB,
                    "condition": condition,
                    "k": int(k),
                    "mixing_entropy": d["entropy"][i],
                    "nmi": d["nmi"][i],
                    "vi": d["vi"][i],
                    "n_nodes": d["n_nodes"],
                })

df = pd.DataFrame(rows)
csv_path = OUTPUT_DIR / "multiscale_community_flow.csv"
df.to_csv(csv_path, index=False, float_format="%.6f")
print(f"  Saved: {csv_path}  ({len(df)} rows)")

# ── Summary statistics for FINDINGS ───────────────────────────────────
print("\n" + "=" * 70)
print("KEY FINDINGS")
print("=" * 70)

findings_lines = []

# 1. Grand average within vs cross (VI)
within_all = []
cross_all = []
for patient in FOUR_PHASE_PATIENTS:
    for band in BANDS:
        for pA, pB in WITHIN_PAIRS:
            if (pA, pB) in results[patient][band]:
                within_all.append(results[patient][band][(pA, pB)]["vi"])
        for pA, pB in CROSS_PAIRS:
            if (pA, pB) in results[patient][band]:
                cross_all.append(results[patient][band][(pA, pB)]["vi"])

within_mean = np.nanmean(within_all, axis=0)
cross_mean = np.nanmean(cross_all, axis=0)

# Scale of maximum separation (cross VI > within VI expected)
diff = cross_mean - within_mean
max_sep_k = K_VALUES[np.nanargmax(diff)]
max_sep_val = np.nanmax(diff)

msg = (f"1. SCALE OF MAXIMUM WITHIN-CROSS VI SEPARATION: k={max_sep_k} "
       f"(delta={max_sep_val:.4f} nats)")
print(msg)
findings_lines.append(msg)

# 1b. Average VI at peak scale
ki_peak = np.argmin(np.abs(K_VALUES - max_sep_k))
msg = (f"   At k={max_sep_k}: within VI={within_mean[ki_peak]:.4f}, "
       f"cross VI={cross_mean[ki_peak]:.4f}")
print(msg)
findings_lines.append(msg)

# 2. Per-band peak reorganization ratio (VI-based)
msg = "\n2. PEAK REORGANIZATION SPECIFICITY RATIO BY BAND (VI-based):"
print(msg)
findings_lines.append(msg)
for band in BANDS:
    if band in ratio_data:
        r = ratio_data[band]
        valid = ~np.isnan(r)
        if valid.any():
            peak_idx = np.nanargmax(r)
            msg = (f"   {BRAIN_BAND_TEX_DICT[band]:>20s}: peak at k={K_VALUES[peak_idx]:2d}, "
                   f"ratio={r[peak_idx]:.3f}")
            print(msg)
            findings_lines.append(msg)

# 3. Cross-patient consistency
msg = "\n3. CROSS-PATIENT CONSISTENCY (peak k across patients, alpha band):"
print(msg)
findings_lines.append(msg)
alpha_peaks = peak_k_matrix[:, BANDS.index("alpha")]
msg = (f"   Alpha band peak k values: {[int(x) if not np.isnan(x) else 'N/A' for x in alpha_peaks]}")
print(msg)
findings_lines.append(msg)
valid_peaks = alpha_peaks[~np.isnan(alpha_peaks)]
if len(valid_peaks) > 1:
    msg = f"   Mean={np.mean(valid_peaks):.1f}, Std={np.std(valid_peaks):.1f}"
    print(msg)
    findings_lines.append(msg)

# 4. Mesoscale specificity check
msg = "\n4. MESOSCALE SPECIFICITY (VI-based):"
print(msg)
findings_lines.append(msg)

for band in BANDS:
    if band in ratio_data:
        r = ratio_data[band]
        valid = ~np.isnan(r)
        if valid.any():
            peak_k = K_VALUES[np.nanargmax(r)]
            # Check if peak is at intermediate scales (5-20) vs extremes
            if 5 <= peak_k <= 20:
                specificity = "MESOSCALE"
            elif peak_k < 5:
                specificity = "COARSE"
            else:
                specificity = "FINE"
            msg = f"   {band:>12s}: peak k={peak_k:2d} -> {specificity}"
            print(msg)
            findings_lines.append(msg)

# 5. Average VI within vs cross at different scales
msg = "\n5. VI WITHIN vs CROSS AT KEY SCALES:"
print(msg)
findings_lines.append(msg)

for k_target in [3, 5, 10, 20]:
    ki = np.argmin(np.abs(K_VALUES - k_target))
    within_vi = []
    cross_vi = []
    for patient in FOUR_PHASE_PATIENTS:
        for band in BANDS:
            for pA, pB in WITHIN_PAIRS:
                if (pA, pB) in results[patient][band]:
                    within_vi.append(results[patient][band][(pA, pB)]["vi"][ki])
            for pA, pB in CROSS_PAIRS:
                if (pA, pB) in results[patient][band]:
                    cross_vi.append(results[patient][band][(pA, pB)]["vi"][ki])

    w_mean = np.nanmean(within_vi)
    c_mean = np.nanmean(cross_vi)
    gap = c_mean - w_mean  # positive = cross is more different (expected)
    msg = (f"   k={k_target:2d}: within VI={w_mean:.4f}, cross VI={c_mean:.4f}, "
           f"gap={gap:+.4f} nats")
    print(msg)
    findings_lines.append(msg)

# 6. Most reorganized phase pair (highest avg VI = most different)
msg = "\n6. MOST REORGANIZED PHASE PAIR (highest avg VI = most different partitions):"
print(msg)
findings_lines.append(msg)

pair_avg_vi = {}
for pA, pB in ALL_PAIRS:
    vals = []
    for patient in FOUR_PHASE_PATIENTS:
        for band in BANDS:
            if (pA, pB) in results[patient][band]:
                vals.append(np.nanmean(results[patient][band][(pA, pB)]["vi"]))
    if vals:
        pair_avg_vi[(pA, pB)] = np.mean(vals)

for pair, val in sorted(pair_avg_vi.items(), key=lambda x: -x[1]):
    cond = "WITHIN" if pair in WITHIN_PAIRS else "CROSS"
    msg = f"   {pair[0]:>10s} - {pair[1]:<10s} ({cond}): {val:.4f} nats"
    print(msg)
    findings_lines.append(msg)

# ── Write FINDINGS.md ──────────────────────────────────────────────────
findings_path = OUTPUT_DIR / "FINDINGS.md"
with open(findings_path, "w") as f:
    f.write("# Cross-Phase Multiscale Community Flow - Findings\n\n")
    f.write("## Method\n\n")
    f.write("For each pair of phases, at each hierarchical level k (2-30), "
            "we cut both LRG dendrograms and compute:\n")
    f.write("- **Variation of Information (VI)**: true metric on partitions, "
            "VI(P,Q) = H(P|Q) + H(Q|P). VI=0 means identical; higher = more "
            "different.\n")
    f.write("- **Mixing entropy** H(B|A): conditional entropy of phase B "
            "partition given phase A partition\n")
    f.write("- **NMI**: normalized mutual information between partitions\n\n")
    f.write("VI is the primary partition comparison metric (information-theoretic, "
            "true metric, decomposes into directional conditional entropies).\n\n")
    f.write("The VI profile across scales k is a **multiscale "
            "reorganization fingerprint**.\n\n")
    f.write("Within-condition pairs: rest_pre-rest_post, task_learn-task_test\n")
    f.write("Cross-condition pairs: rest_pre-task_learn, rest_pre-task_test, "
            "task_learn-rest_post, task_test-rest_post\n\n")
    f.write(f"FC method: {FC_METHOD}, 4-phase patients: "
            f"{', '.join(FOUR_PHASE_PATIENTS)}\n\n")
    f.write("## Results\n\n")
    f.write("```\n")
    for line in findings_lines:
        f.write(line + "\n")
    f.write("```\n\n")
    f.write("## Figures\n\n")
    f.write("- `fig1_mixing_entropy_profiles.png` - Mixing entropy vs k, "
            "per band, all phase pairs\n")
    f.write("- `fig2_nmi_profiles.png` - NMI vs k, per band, all phase pairs\n")
    f.write("- `fig3_vi_profiles.png` - VI vs k, per band, all phase pairs "
            "(primary metric)\n")
    f.write("- `fig4_within_vs_cross_grand_average.png` - Grand average "
            "within vs cross-condition (VI + mixing entropy)\n")
    f.write("- `fig5_reorganization_specificity_ratio.png` - Cross/within "
            "VI ratio per band\n")
    f.write("- `fig6_peak_reorganization_heatmap.png` - Peak k across "
            "patients x bands (VI-based)\n")
    f.write("- `fig7_per_patient_alpha_vi_profiles.png` - Individual patient "
            "VI profiles for alpha\n")

print(f"\n  Saved: {findings_path}")
print("\nDone!")
