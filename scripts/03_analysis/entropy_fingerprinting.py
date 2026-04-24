#!/usr/bin/env python3
"""Entropy Curve Fingerprinting Analysis.

Investigates whether LRG entropy curves (1-S and C) capture
task-specific multiscale reorganization patterns in EEG data.
"""

import warnings
from pathlib import Path
from itertools import combinations

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, mannwhitneyu

from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS, BRAIN_BAND_TEX_DICT, PATIENTS_4PHASE
from lrg_eegfc.config.paths import FIGURES_ROOT

# ── Configuration ──────────────────────────────────────────────────────
FC_METHOD = "msc"
ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]
FOUR_PHASE_PATIENTS = PATIENTS_4PHASE  # exclude Pat_06, Pat_07
BANDS = list(BRAIN_BANDS.keys())
PHASES = list(PHASE_LABELS)  # rest_pre, task_learn, task_test, rest_post

OUTPUT_DIR = FIGURES_ROOT / "multiscale_investigation" / "entropy_fingerprints"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Phase categories
REST_PHASES = ["rest_pre", "rest_post"]
TASK_PHASES = ["task_learn", "task_test"]

PHASE_COLORS = {
    "rest_pre": "#1f77b4",      # blue
    "rest_post": "#6baed6",     # light blue
    "task_learn": "#d62728",  # red
    "task_test": "#ff7f7f",   # light red
}
PHASE_DISPLAY = {
    "rest_pre": "Rest Pre",
    "rest_post": "Rest Post",
    "task_learn": "Task Learn",
    "task_test": "Task Test",
}

# ── 1. Load all LRG results ───────────────────────────────────────────
print("=" * 60)
print("ENTROPY CURVE FINGERPRINTING ANALYSIS")
print("=" * 60)

results = {}  # results[(patient, band, phase)] = LRGResult
loaded, missing = 0, 0

for patient in ALL_PATIENTS:
    phases_for_patient = PHASES if patient in FOUR_PHASE_PATIENTS else TASK_PHASES
    for band in BANDS:
        for phase in phases_for_patient:
            r = load_lrg_result(patient, phase, band, FC_METHOD)
            if r is not None:
                results[(patient, band, phase)] = r
                loaded += 1
            else:
                missing += 1
                print(f"  MISSING: {patient} {band} {phase}")

print(f"\nLoaded {loaded} LRG results, {missing} missing.")

# ── 2. Compute curve metrics ──────────────────────────────────────────
print("\nComputing curve metrics...")

rows = []
for (patient, band, phase), r in results.items():
    tau = r.entropy_tau
    one_minus_S = r.entropy_1_minus_S
    C = r.entropy_C

    # Area under 1-S curve (trapezoidal integration over log-tau)
    auc_1mS = np.trapz(one_minus_S, np.log10(tau))

    # C peak
    C_peak_idx = np.argmax(C)
    C_peak_val = C[C_peak_idx]
    # tau for C: midpoints between tau values (C has len-1)
    tau_mid = 0.5 * (tau[:-1] + tau[1:])
    C_peak_tau = tau_mid[C_peak_idx]

    rows.append({
        "patient": patient,
        "band": band,
        "phase": phase,
        "auc_1mS": auc_1mS,
        "C_peak_val": C_peak_val,
        "C_peak_tau": C_peak_tau,
        "n_nodes": r.n_nodes,
    })

df_metrics = pd.DataFrame(rows)
print(f"  Metrics table: {len(df_metrics)} entries")

# ── 3. Pairwise curve distances ──────────────────────────────────────
print("\nComputing pairwise curve distances...")

pair_rows = []
phase_pairs = list(combinations(PHASES, 2))

for patient in FOUR_PHASE_PATIENTS:
    for band in BANDS:
        # Collect curves for this patient/band
        curves = {}
        for phase in PHASES:
            key = (patient, band, phase)
            if key in results:
                curves[phase] = results[key]

        if len(curves) < 4:
            continue

        for p1, p2 in phase_pairs:
            r1, r2 = curves[p1], curves[p2]

            # Ensure same length (should be 400 for 1-S, 399 for C)
            min_len_S = min(len(r1.entropy_1_minus_S), len(r2.entropy_1_minus_S))
            s1 = r1.entropy_1_minus_S[:min_len_S]
            s2 = r2.entropy_1_minus_S[:min_len_S]

            min_len_C = min(len(r1.entropy_C), len(r2.entropy_C))
            c1 = r1.entropy_C[:min_len_C]
            c2 = r2.entropy_C[:min_len_C]

            # L2 distance for 1-S
            l2_S = np.sqrt(np.mean((s1 - s2) ** 2))
            # Pearson correlation for 1-S
            corr_S, pval_S = pearsonr(s1, s2)

            # L2 distance for C
            l2_C = np.sqrt(np.mean((c1 - c2) ** 2))
            # Pearson correlation for C
            corr_C, pval_C = pearsonr(c1, c2)

            # Classify pair
            if (p1 in REST_PHASES and p2 in REST_PHASES):
                condition = "within_rest"
            elif (p1 in TASK_PHASES and p2 in TASK_PHASES):
                condition = "within_task"
            elif ((p1 in REST_PHASES and p2 in TASK_PHASES) or
                  (p1 in TASK_PHASES and p2 in REST_PHASES)):
                condition = "cross"
            else:
                condition = "other"

            pair_rows.append({
                "patient": patient,
                "band": band,
                "phase1": p1,
                "phase2": p2,
                "l2_S": l2_S,
                "corr_S": corr_S,
                "pval_S": pval_S,
                "l2_C": l2_C,
                "corr_C": corr_C,
                "pval_C": pval_C,
                "condition": condition,
            })

df_pairs = pd.DataFrame(pair_rows)
print(f"  Pairwise table: {len(df_pairs)} entries")

# Aggregate condition type
within = df_pairs[df_pairs["condition"].isin(["within_rest", "within_task"])]
cross = df_pairs[df_pairs["condition"] == "cross"]

print(f"\n  Within-condition pairs: {len(within)} ({len(within[within['condition']=='within_rest'])} rest, {len(within[within['condition']=='within_task'])} task)")
print(f"  Cross-condition pairs: {len(cross)}")

# ── 4. Statistical tests ─────────────────────────────────────────────
print("\n--- Statistical Summary ---")

for metric, label in [("l2_S", "L2 dist (1-S)"), ("l2_C", "L2 dist (C)"),
                       ("corr_S", "Corr (1-S)"), ("corr_C", "Corr (C)")]:
    w_vals = within[metric].values
    c_vals = cross[metric].values
    stat, pval = mannwhitneyu(w_vals, c_vals, alternative="two-sided")
    print(f"  {label}: within={np.mean(w_vals):.4f}+/-{np.std(w_vals):.4f}, "
          f"cross={np.mean(c_vals):.4f}+/-{np.std(c_vals):.4f}, "
          f"U={stat:.0f}, p={pval:.4e}")

# Per-band analysis
print("\n--- Per-Band Analysis (L2 distance, 1-S) ---")
band_effects = []
for band in BANDS:
    b_within = within[within["band"] == band]["l2_S"]
    b_cross = cross[cross["band"] == band]["l2_S"]
    if len(b_within) > 0 and len(b_cross) > 0:
        stat, pval = mannwhitneyu(b_within.values, b_cross.values, alternative="two-sided")
        effect = np.mean(b_cross.values) - np.mean(b_within.values)
        print(f"  {band:12s}: within={np.mean(b_within):.4f}, cross={np.mean(b_cross):.4f}, "
              f"diff={effect:+.4f}, p={pval:.4e}")
        band_effects.append({
            "band": band,
            "within_mean_l2S": np.mean(b_within),
            "cross_mean_l2S": np.mean(b_cross),
            "effect_l2S": effect,
            "pval_l2S": pval,
        })

        # Also for C
        b_within_C = within[within["band"] == band]["l2_C"]
        b_cross_C = cross[cross["band"] == band]["l2_C"]
        stat_C, pval_C = mannwhitneyu(b_within_C.values, b_cross_C.values, alternative="two-sided")
        effect_C = np.mean(b_cross_C.values) - np.mean(b_within_C.values)
        band_effects[-1].update({
            "within_mean_l2C": np.mean(b_within_C),
            "cross_mean_l2C": np.mean(b_cross_C),
            "effect_l2C": effect_C,
            "pval_l2C": pval_C,
        })

df_band_effects = pd.DataFrame(band_effects)

# ── 5. FIGURES ────────────────────────────────────────────────────────
print("\n\nGenerating figures...")

# ── Fig 1: Overlay of 1-S curves ─────────────────────────────────────
for patient in FOUR_PHASE_PATIENTS:
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle(f"Entropy Curves (1-S) across Phases — {patient}", fontsize=16, y=0.98)

    for idx, band in enumerate(BANDS):
        ax = axes.flat[idx]
        ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]} ({band})", fontsize=12)

        for phase in PHASES:
            key = (patient, band, phase)
            if key in results:
                r = results[key]
                ax.plot(
                    np.log10(r.entropy_tau),
                    r.entropy_1_minus_S,
                    color=PHASE_COLORS[phase],
                    label=PHASE_DISPLAY[phase],
                    linewidth=1.8,
                    alpha=0.85,
                )

        ax.set_xlabel(r"$\log_{10}(\tau)$")
        ax.set_ylabel(r"$1 - S(\tau)$")
        ax.set_ylim(-0.05, 1.05)
        if idx == 0:
            ax.legend(fontsize=9, loc="best")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fpath = OUTPUT_DIR / f"fig1_entropy_overlay_{patient}.png"
    fig.savefig(fpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fpath}")

# Also do C(tau) overlay
for patient in FOUR_PHASE_PATIENTS:
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle(f"Spectral Complexity C(τ) across Phases — {patient}", fontsize=16, y=0.98)

    for idx, band in enumerate(BANDS):
        ax = axes.flat[idx]
        ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]} ({band})", fontsize=12)

        for phase in PHASES:
            key = (patient, band, phase)
            if key in results:
                r = results[key]
                tau_mid = 0.5 * (r.entropy_tau[:-1] + r.entropy_tau[1:])
                ax.plot(
                    np.log10(tau_mid),
                    r.entropy_C,
                    color=PHASE_COLORS[phase],
                    label=PHASE_DISPLAY[phase],
                    linewidth=1.8,
                    alpha=0.85,
                )

        ax.set_xlabel(r"$\log_{10}(\tau)$")
        ax.set_ylabel(r"$C(\tau)$")
        if idx == 0:
            ax.legend(fontsize=9, loc="best")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fpath = OUTPUT_DIR / f"fig1_complexity_overlay_{patient}.png"
    fig.savefig(fpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fpath}")

# ── Fig 2: Heatmap of pairwise L2 distances (averaged) ──────────────
# Average across all patients and bands
phase_order = PHASES
n_phases = len(phase_order)

for curve_type, metric in [("1-S", "l2_S"), ("C", "l2_C")]:
    avg_matrix = np.zeros((n_phases, n_phases))
    count_matrix = np.zeros((n_phases, n_phases))

    for _, row in df_pairs.iterrows():
        i = phase_order.index(row["phase1"])
        j = phase_order.index(row["phase2"])
        avg_matrix[i, j] += row[metric]
        avg_matrix[j, i] += row[metric]
        count_matrix[i, j] += 1
        count_matrix[j, i] += 1

    # Avoid division by zero on diagonal
    mask = count_matrix > 0
    avg_matrix[mask] /= count_matrix[mask]

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(avg_matrix, cmap="YlOrRd", aspect="equal")
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label(f"Mean L2 distance ({curve_type})", fontsize=11)

    phase_labels_display = [PHASE_DISPLAY[p] for p in phase_order]
    ax.set_xticks(range(n_phases))
    ax.set_xticklabels(phase_labels_display, rotation=45, ha="right")
    ax.set_yticks(range(n_phases))
    ax.set_yticklabels(phase_labels_display)

    # Annotate values
    for i in range(n_phases):
        for j in range(n_phases):
            val = avg_matrix[i, j]
            if i != j:
                ax.text(j, i, f"{val:.4f}", ha="center", va="center",
                        fontsize=10, color="black" if val < np.max(avg_matrix)*0.7 else "white")

    ax.set_title(f"Pairwise Entropy Curve Distance ({curve_type})\nAveraged across patients & bands", fontsize=13)
    plt.tight_layout()
    fpath = OUTPUT_DIR / f"fig2_distance_heatmap_{curve_type.replace('-','')}.png"
    fig.savefig(fpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fpath}")

# ── Fig 3: Boxplot within vs cross condition ─────────────────────────
fig, axes = plt.subplots(1, 4, figsize=(18, 5))

for idx, (metric, label) in enumerate([
    ("l2_S", "L2 dist (1-S)"),
    ("l2_C", "L2 dist (C)"),
    ("corr_S", "Pearson r (1-S)"),
    ("corr_C", "Pearson r (C)"),
]):
    ax = axes[idx]

    within_rest = df_pairs[df_pairs["condition"] == "within_rest"][metric].values
    within_task = df_pairs[df_pairs["condition"] == "within_task"][metric].values
    cross_vals = df_pairs[df_pairs["condition"] == "cross"][metric].values

    data = [within_rest, within_task, cross_vals]
    labels_box = ["Rest\n(rest_pre-rest_post)", "Task\n(Learn-Test)", "Cross\n(rest-task)"]
    colors = ["#6baed6", "#ff7f7f", "#9467bd"]

    bp = ax.boxplot(data, patch_artist=True, labels=labels_box, widths=0.6)
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Overlay individual points
    for i, d in enumerate(data):
        jitter = np.random.RandomState(42).uniform(-0.12, 0.12, len(d))
        ax.scatter(np.full(len(d), i + 1) + jitter, d, alpha=0.35, s=15,
                   color=colors[i], edgecolors="none")

    ax.set_title(label, fontsize=12)
    ax.set_ylabel(label)

    # Add significance annotation
    stat, pval = mannwhitneyu(
        np.concatenate([within_rest, within_task]),
        cross_vals, alternative="two-sided"
    )
    sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "n.s."
    ymax = max(np.max(d) for d in data if len(d) > 0) * 1.05
    ax.text(2, ymax, f"within vs cross: {sig}\n(p={pval:.3e})", ha="center",
            fontsize=9, style="italic")

fig.suptitle("Entropy Curve Distances: Within-Condition vs Cross-Condition",
             fontsize=14, y=1.02)
plt.tight_layout()
fpath = OUTPUT_DIR / "fig3_within_vs_cross_boxplot.png"
fig.savefig(fpath, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {fpath}")

# ── Fig 4: Band-specific fingerprint effect ──────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for idx, (metric_w, metric_c, metric_eff, metric_p, curve_label) in enumerate([
    ("within_mean_l2S", "cross_mean_l2S", "effect_l2S", "pval_l2S", "1-S"),
    ("within_mean_l2C", "cross_mean_l2C", "effect_l2C", "pval_l2C", "C"),
]):
    ax = axes[idx]
    x = np.arange(len(BANDS))
    width = 0.35

    bars1 = ax.bar(x - width/2, df_band_effects[metric_w], width,
                   label="Within-condition", color="#6baed6", alpha=0.8)
    bars2 = ax.bar(x + width/2, df_band_effects[metric_c], width,
                   label="Cross-condition", color="#ff7f7f", alpha=0.8)

    # Significance stars
    for i, (_, row) in enumerate(df_band_effects.iterrows()):
        pval = row[metric_p]
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
        if sig:
            ymax = max(row[metric_w], row[metric_c])
            ax.text(i, ymax * 1.05, sig, ha="center", fontsize=12, fontweight="bold")

    band_tex = [BRAIN_BAND_TEX_DICT[b] for b in BANDS]
    ax.set_xticks(x)
    ax.set_xticklabels(band_tex, fontsize=11)
    ax.set_ylabel(f"Mean L2 Distance ({curve_label})")
    ax.set_title(f"Band-Specific Fingerprint Effect ({curve_label})", fontsize=13)
    ax.legend()

plt.tight_layout()
fpath = OUTPUT_DIR / "fig4_band_specific_effect.png"
fig.savefig(fpath, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {fpath}")

# ── Additional: Per-patient band effect heatmap ──────────────────────
fig, axes = plt.subplots(1, 2, figsize=(16, 5))

for idx, (metric, label) in enumerate([("l2_S", "L2(1-S)"), ("l2_C", "L2(C)")]):
    ax = axes[idx]
    effect_mat = np.zeros((len(FOUR_PHASE_PATIENTS), len(BANDS)))

    for pi, patient in enumerate(FOUR_PHASE_PATIENTS):
        for bi, band in enumerate(BANDS):
            w = df_pairs[(df_pairs["patient"] == patient) &
                         (df_pairs["band"] == band) &
                         (df_pairs["condition"].isin(["within_rest", "within_task"]))][metric]
            c = df_pairs[(df_pairs["patient"] == patient) &
                         (df_pairs["band"] == band) &
                         (df_pairs["condition"] == "cross")][metric]
            if len(w) > 0 and len(c) > 0:
                effect_mat[pi, bi] = np.mean(c.values) - np.mean(w.values)

    im = ax.imshow(effect_mat, cmap="RdBu_r", aspect="auto",
                   vmin=-np.max(np.abs(effect_mat)), vmax=np.max(np.abs(effect_mat)))
    fig.colorbar(im, ax=ax, shrink=0.8, label=f"Effect (cross - within)")

    band_tex = [BRAIN_BAND_TEX_DICT[b] for b in BANDS]
    ax.set_xticks(range(len(BANDS)))
    ax.set_xticklabels(band_tex, fontsize=11)
    ax.set_yticks(range(len(FOUR_PHASE_PATIENTS)))
    ax.set_yticklabels(FOUR_PHASE_PATIENTS)
    ax.set_title(f"Fingerprint Effect by Patient & Band ({label})", fontsize=12)

    # Annotate
    for pi in range(len(FOUR_PHASE_PATIENTS)):
        for bi in range(len(BANDS)):
            val = effect_mat[pi, bi]
            ax.text(bi, pi, f"{val:+.3f}", ha="center", va="center", fontsize=8,
                    color="white" if abs(val) > np.max(np.abs(effect_mat)) * 0.5 else "black")

plt.tight_layout()
fpath = OUTPUT_DIR / "fig5_patient_band_effect_heatmap.png"
fig.savefig(fpath, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {fpath}")

# ── 6. Save CSVs ─────────────────────────────────────────────────────
print("\nSaving CSV files...")

df_metrics.to_csv(OUTPUT_DIR / "entropy_curve_metrics.csv", index=False)
print(f"  Saved: {OUTPUT_DIR / 'entropy_curve_metrics.csv'}")

df_pairs.to_csv(OUTPUT_DIR / "pairwise_distances.csv", index=False)
print(f"  Saved: {OUTPUT_DIR / 'pairwise_distances.csv'}")

df_band_effects.to_csv(OUTPUT_DIR / "band_specific_effects.csv", index=False)
print(f"  Saved: {OUTPUT_DIR / 'band_specific_effects.csv'}")

# ── 7. Write FINDINGS.md ─────────────────────────────────────────────
print("\nWriting findings summary...")

# Compute summary stats for the report
overall_within_l2S = within["l2_S"].mean()
overall_cross_l2S = cross["l2_S"].mean()
stat_overall, pval_overall = mannwhitneyu(within["l2_S"].values, cross["l2_S"].values, alternative="two-sided")

overall_within_l2C = within["l2_C"].mean()
overall_cross_l2C = cross["l2_C"].mean()
stat_overall_C, pval_overall_C = mannwhitneyu(within["l2_C"].values, cross["l2_C"].values, alternative="two-sided")

overall_within_corrS = within["corr_S"].mean()
overall_cross_corrS = cross["corr_S"].mean()

# Best discriminating band
best_band_row = df_band_effects.loc[df_band_effects["effect_l2S"].idxmax()]
best_band = best_band_row["band"]

# Correlation summary
overall_within_corrC = within["corr_C"].mean()
overall_cross_corrC = cross["corr_C"].mean()

findings = f"""# Entropy Curve Fingerprinting — Findings

**Date:** 2026-03-10
**FC method:** MSC (Magnitude-Squared Coherence)
**Patients (4-phase):** {', '.join(FOUR_PHASE_PATIENTS)}
**Frequency bands:** {', '.join(BANDS)}

---

## Summary

This analysis tested whether LRG entropy curves — specifically the normalized
entropy 1-S(tau) and spectral complexity C(tau) — carry condition-specific
"fingerprints" that distinguish rest from task states. If the multiscale network
organization is modulated by cognitive state, we expect entropy curves from the
same condition (rest-rest or task-task) to be more similar than curves from
different conditions (rest-task).

## Key Results

### Overall Distances (1-S curves)

| Comparison | Mean L2 Distance | Std |
|-----------|-----------------|-----|
| Within-condition (rest+task) | {overall_within_l2S:.4f} | {within['l2_S'].std():.4f} |
| Cross-condition (rest-task) | {overall_cross_l2S:.4f} | {cross['l2_S'].std():.4f} |
| **Mann-Whitney U** | **{stat_overall:.0f}** | **p = {pval_overall:.4e}** |

### Overall Distances (C curves)

| Comparison | Mean L2 Distance | Std |
|-----------|-----------------|-----|
| Within-condition | {overall_within_l2C:.4f} | {within['l2_C'].std():.4f} |
| Cross-condition | {overall_cross_l2C:.4f} | {cross['l2_C'].std():.4f} |
| **Mann-Whitney U** | **{stat_overall_C:.0f}** | **p = {pval_overall_C:.4e}** |

### Correlation Summary

| Metric | Within-condition | Cross-condition |
|--------|-----------------|-----------------|
| Pearson r (1-S) | {overall_within_corrS:.4f} | {overall_cross_corrS:.4f} |
| Pearson r (C) | {overall_within_corrC:.4f} | {overall_cross_corrC:.4f} |

### Band-Specific Effects (L2 distance, 1-S)

| Band | Within | Cross | Effect | p-value |
|------|--------|-------|--------|---------|
"""

for _, row in df_band_effects.iterrows():
    sig = "***" if row["pval_l2S"] < 0.001 else "**" if row["pval_l2S"] < 0.01 else "*" if row["pval_l2S"] < 0.05 else ""
    findings += f"| {row['band']} | {row['within_mean_l2S']:.4f} | {row['cross_mean_l2S']:.4f} | {row['effect_l2S']:+.4f} | {row['pval_l2S']:.4e} {sig} |\n"

findings += f"""
### Interpretation

"""

# Determine direction of effect
if overall_cross_l2S > overall_within_l2S:
    direction = "larger"
    supports = "supports"
else:
    direction = "smaller"
    supports = "does not support"

sig_label = "significant" if pval_overall < 0.05 else "not significant"

findings += f"""- Cross-condition L2 distances are **{direction}** than within-condition distances
  for 1-S curves (p = {pval_overall:.4e}, {sig_label}).
- This {supports} the hypothesis that entropy curves capture task-specific
  multiscale reorganization.
- The band with the strongest fingerprint effect is **{best_band}**
  (cross-within difference = {best_band_row['effect_l2S']:+.4f}, p = {best_band_row['pval_l2S']:.4e}).

## Figures

- `fig1_entropy_overlay_PatXX.png` — 1-S(tau) and C(tau) curves per phase/band
- `fig2_distance_heatmap_*.png` — Average pairwise L2 distances
- `fig3_within_vs_cross_boxplot.png` — Within vs cross-condition comparison
- `fig4_band_specific_effect.png` — Band-resolved fingerprint effect
- `fig5_patient_band_effect_heatmap.png` — Patient x band effect matrix

## Data Files

- `entropy_curve_metrics.csv` — Per-condition curve summary metrics
- `pairwise_distances.csv` — All pairwise distance computations
- `band_specific_effects.csv` — Band-specific statistical results
"""

with open(OUTPUT_DIR / "FINDINGS.md", "w") as f:
    f.write(findings)
print(f"  Saved: {OUTPUT_DIR / 'FINDINGS.md'}")

print("\n" + "=" * 60)
print("ANALYSIS COMPLETE")
print("=" * 60)
