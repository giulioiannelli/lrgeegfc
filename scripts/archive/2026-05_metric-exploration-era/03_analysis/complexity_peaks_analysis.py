"""Spectral Complexity Peaks as Natural Scale Markers.

Analyzes C(tau) curves from LRG results to identify peaks that indicate
natural scales of brain network hierarchical organization.
"""

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.signal import find_peaks
from scipy.integrate import trapezoid

from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.config.const import (
    BRAIN_BANDS, PHASE_LABELS, BRAIN_BAND_TEX_DICT, BRAIN_BANDS_NAMES,
    PATIENTS_4PHASE,
)
from lrg_eegfc.config.paths import FIGURES_ROOT

# ---------- Configuration ----------
FC_METHOD = "msc"
ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]
# Pat_06/Pat_07 only have task_learn/task_test
CROSS_PHASE_PATIENTS = PATIENTS_4PHASE
BANDS = list(BRAIN_BANDS.keys())
PHASES = list(PHASE_LABELS)

OUTPUT_DIR = FIGURES_ROOT / "multiscale_investigation" / "complexity_peaks"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Peak detection params
PROMINENCE_THRESHOLD = 0.005
HEIGHT_THRESHOLD = 0.01
DISTANCE = 5  # minimum distance between peaks in samples

# Phase colors
PHASE_COLORS = {
    "rest_pre": "#2196F3",
    "task_learn": "#FF9800",
    "task_test": "#F44336",
    "rest_post": "#4CAF50",
}
PHASE_NICE = {
    "rest_pre": "Rest Pre",
    "task_learn": "Task Learn",
    "task_test": "Task Test",
    "rest_post": "Rest Post",
}

# ---------- 1. Load all LRG results and extract peak info ----------
print("=" * 60)
print("Loading LRG results and analyzing C(tau) peaks...")
print("=" * 60)

records = []

for patient in ALL_PATIENTS:
    phases_for_patient = PHASES if patient in CROSS_PHASE_PATIENTS else ["task_learn", "task_test"]
    for phase in phases_for_patient:
        for band in BANDS:
            result = load_lrg_result(patient, phase, band, FC_METHOD)
            if result is None:
                print(f"  MISSING: {patient} {phase} {band}")
                continue

            tau = result.entropy_tau
            C = result.entropy_C

            # C has 399 values, tau has 400 -- C is derivative-like, use midpoints
            tau_mid = 0.5 * (tau[:-1] + tau[1:])

            # Find peaks
            peaks, props = find_peaks(
                C,
                prominence=PROMINENCE_THRESHOLD,
                height=HEIGHT_THRESHOLD,
                distance=DISTANCE,
            )

            n_peaks = len(peaks)
            area = float(trapezoid(C, tau_mid))

            # Main peak info
            if n_peaks > 0:
                prominences = props["prominences"]
                heights = props["peak_heights"]
                main_idx = np.argmax(prominences)
                main_peak_tau = float(tau_mid[peaks[main_idx]])
                main_peak_height = float(heights[main_idx])
                main_peak_prominence = float(prominences[main_idx])

                # FWHM of main peak (approximate)
                half_max = main_peak_height / 2.0
                peak_pos = peaks[main_idx]
                # walk left
                left = peak_pos
                while left > 0 and C[left] > half_max:
                    left -= 1
                # walk right
                right = peak_pos
                while right < len(C) - 1 and C[right] > half_max:
                    right += 1
                fwhm_tau = float(tau_mid[min(right, len(tau_mid)-1)] - tau_mid[max(left, 0)])

                all_peak_taus = [float(tau_mid[p]) for p in peaks]
                all_peak_heights = [float(h) for h in heights]
            else:
                main_peak_tau = np.nan
                main_peak_height = np.nan
                main_peak_prominence = np.nan
                fwhm_tau = np.nan
                all_peak_taus = []
                all_peak_heights = []

            is_rest = phase in ("rest_pre", "rest_post")
            condition = "rest" if is_rest else "task"

            records.append({
                "patient": patient,
                "phase": phase,
                "band": band,
                "condition": condition,
                "n_peaks": n_peaks,
                "main_peak_tau": main_peak_tau,
                "main_peak_height": main_peak_height,
                "main_peak_prominence": main_peak_prominence,
                "main_peak_fwhm": fwhm_tau,
                "area_under_C": area,
                "all_peak_taus": all_peak_taus,
                "all_peak_heights": all_peak_heights,
                "n_nodes": result.n_nodes,
            })

df = pd.DataFrame(records)
print(f"\nLoaded {len(df)} LRG results across {df['patient'].nunique()} patients")
print(f"Peaks detected: mean={df['n_peaks'].mean():.1f}, range=[{df['n_peaks'].min()}, {df['n_peaks'].max()}]")

# Save CSV (without list columns)
csv_df = df.drop(columns=["all_peak_taus", "all_peak_heights"])
csv_df.to_csv(OUTPUT_DIR / "complexity_peaks_results.csv", index=False)
print(f"CSV saved to {OUTPUT_DIR / 'complexity_peaks_results.csv'}")

# ---------- Figure 1: C(tau) curves for Pat_02, one subplot per band ----------
print("\n--- Figure 1: C(tau) curves for Pat_02 ---")

fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharex=True, sharey=False)
axes = axes.flatten()

for i, band in enumerate(BANDS):
    ax = axes[i]
    for phase in PHASES:
        result = load_lrg_result("Pat_02", phase, band, FC_METHOD)
        if result is None:
            continue
        tau = result.entropy_tau
        C = result.entropy_C
        tau_mid = 0.5 * (tau[:-1] + tau[1:])
        log_tau_mid = np.log10(tau_mid)
        ax.plot(log_tau_mid, C, label=PHASE_NICE[phase], color=PHASE_COLORS[phase],
                linewidth=1.5, alpha=0.85)

        # Mark peaks
        peaks, props = find_peaks(C, prominence=PROMINENCE_THRESHOLD,
                                  height=HEIGHT_THRESHOLD, distance=DISTANCE)
        if len(peaks) > 0:
            ax.scatter(log_tau_mid[peaks], C[peaks], color=PHASE_COLORS[phase],
                      s=30, zorder=5, edgecolors="k", linewidths=0.5)

    ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]} ({band})", fontsize=13)
    ax.set_xlabel(r"$\log_{10}(\tau)$")
    ax.set_ylabel(r"$C(\tau)$")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3)

fig.suptitle("Pat_02 — Spectral Complexity C(τ) Across Phases\n(peaks marked with dots)",
             fontsize=15, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(OUTPUT_DIR / "fig1_ctau_curves_pat02.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("  Saved fig1_ctau_curves_pat02.png")


# ---------- Figure 2: Number of peaks, rest vs task (boxplot) ----------
print("\n--- Figure 2: Number of peaks, rest vs task ---")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# 2a: Overall rest vs task
ax = axes[0]
rest_peaks = df[df["condition"] == "rest"]["n_peaks"]
task_peaks = df[df["condition"] == "task"]["n_peaks"]
bp = ax.boxplot([rest_peaks, task_peaks], labels=["Rest", "Task"],
                patch_artist=True, widths=0.5)
bp["boxes"][0].set_facecolor("#90CAF9")
bp["boxes"][1].set_facecolor("#FFCC80")
ax.set_ylabel("Number of Peaks")
ax.set_title("Number of Complexity Peaks:\nRest vs Task (all patients/bands)")
ax.grid(True, alpha=0.3, axis="y")

# Add individual points
for j, (data, x_pos) in enumerate([(rest_peaks, 1), (task_peaks, 2)]):
    jitter = np.random.default_rng(42).uniform(-0.12, 0.12, size=len(data))
    ax.scatter(np.full(len(data), x_pos) + jitter, data,
              alpha=0.3, s=15, color="black", zorder=3)

# stats annotation
from scipy.stats import mannwhitneyu
stat, pval = mannwhitneyu(rest_peaks, task_peaks, alternative="two-sided")
ax.text(0.5, 0.95, f"Mann-Whitney U p={pval:.4f}",
        transform=ax.transAxes, ha="center", va="top", fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.5))

# 2b: Per band
ax = axes[1]
band_data = []
band_labels_list = []
for band in BANDS:
    for cond in ["rest", "task"]:
        subset = df[(df["band"] == band) & (df["condition"] == cond)]["n_peaks"]
        band_data.append(subset.values)
        band_labels_list.append(f"{BRAIN_BAND_TEX_DICT[band]}\n{cond}")

positions = []
for i in range(len(BANDS)):
    positions.extend([i * 3 + 1, i * 3 + 2])

bp = ax.boxplot(band_data, positions=positions, patch_artist=True, widths=0.6)
for j, box in enumerate(bp["boxes"]):
    box.set_facecolor("#90CAF9" if j % 2 == 0 else "#FFCC80")

ax.set_xticks([(i * 3 + 1.5) for i in range(len(BANDS))])
ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=11)
ax.set_ylabel("Number of Peaks")
ax.set_title("Number of Peaks per Band:\nRest (blue) vs Task (orange)")
ax.grid(True, alpha=0.3, axis="y")

fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig2_peaks_rest_vs_task.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("  Saved fig2_peaks_rest_vs_task.png")


# ---------- Figure 3: Main peak tau position across phases ----------
print("\n--- Figure 3: Main peak tau position across phases ---")

fig, axes = plt.subplots(2, 3, figsize=(16, 10))
axes = axes.flatten()

# Use only cross-phase patients
df_cross = df[df["patient"].isin(CROSS_PHASE_PATIENTS)]

for i, band in enumerate(BANDS):
    ax = axes[i]
    for j, phase in enumerate(PHASES):
        subset = df_cross[(df_cross["band"] == band) & (df_cross["phase"] == phase)]
        vals = subset["main_peak_tau"].dropna()
        if len(vals) == 0:
            continue
        mean_val = vals.mean()
        std_val = vals.std()
        ax.bar(j, mean_val, yerr=std_val, color=PHASE_COLORS[phase],
               capsize=5, edgecolor="black", linewidth=0.5, alpha=0.8,
               label=PHASE_NICE[phase])
        # individual patient dots
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(np.full(len(vals), j) + jitter, vals,
                  color="black", s=20, alpha=0.5, zorder=5)

    ax.set_xticks(range(len(PHASES)))
    ax.set_xticklabels([PHASE_NICE[p] for p in PHASES], fontsize=9, rotation=15)
    ax.set_ylabel(r"$\tau$ of main peak")
    ax.set_title(f"{BRAIN_BAND_TEX_DICT[band]}", fontsize=13)
    ax.grid(True, alpha=0.3, axis="y")

fig.suptitle("Main Peak Position (τ) Across Phases\n(cross-phase patients, mean ± SD)",
             fontsize=15, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(OUTPUT_DIR / "fig3_main_peak_tau_phases.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("  Saved fig3_main_peak_tau_phases.png")


# ---------- Figure 4: Peak height heatmap (phases x bands) ----------
print("\n--- Figure 4: Peak height heatmap ---")

fig, ax = plt.subplots(figsize=(10, 5))

# Average main peak height across cross-phase patients
heatmap_data = np.full((len(PHASES), len(BANDS)), np.nan)
for i, phase in enumerate(PHASES):
    for j, band in enumerate(BANDS):
        subset = df_cross[(df_cross["phase"] == phase) & (df_cross["band"] == band)]
        vals = subset["main_peak_height"].dropna()
        if len(vals) > 0:
            heatmap_data[i, j] = vals.mean()

im = ax.imshow(heatmap_data, aspect="auto", cmap="YlOrRd", interpolation="nearest")
ax.set_xticks(range(len(BANDS)))
ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=12)
ax.set_yticks(range(len(PHASES)))
ax.set_yticklabels([PHASE_NICE[p] for p in PHASES], fontsize=11)

# Annotate cells
for i in range(len(PHASES)):
    for j in range(len(BANDS)):
        val = heatmap_data[i, j]
        if not np.isnan(val):
            ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                    fontsize=9, color="black" if val < 0.5 * np.nanmax(heatmap_data) else "white")

cbar = fig.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label("Mean Dominant Peak Height", fontsize=11)
ax.set_title("Dominant Peak Height — Phase × Band\n(averaged across patients)",
             fontsize=14, fontweight="bold")
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig4_peak_height_heatmap.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("  Saved fig4_peak_height_heatmap.png")


# ---------- Figure 5: Total spectral complexity (AUC) ----------
print("\n--- Figure 5: Total spectral complexity (area under C) ---")

fig, ax = plt.subplots(figsize=(10, 5))

auc_data = np.full((len(PHASES), len(BANDS)), np.nan)
for i, phase in enumerate(PHASES):
    for j, band in enumerate(BANDS):
        subset = df_cross[(df_cross["phase"] == phase) & (df_cross["band"] == band)]
        vals = subset["area_under_C"].dropna()
        if len(vals) > 0:
            auc_data[i, j] = vals.mean()

im = ax.imshow(auc_data, aspect="auto", cmap="viridis", interpolation="nearest")
ax.set_xticks(range(len(BANDS)))
ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in BANDS], fontsize=12)
ax.set_yticks(range(len(PHASES)))
ax.set_yticklabels([PHASE_NICE[p] for p in PHASES], fontsize=11)

for i in range(len(PHASES)):
    for j in range(len(BANDS)):
        val = auc_data[i, j]
        if not np.isnan(val):
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=9, color="white" if val > 0.5 * np.nanmax(auc_data) else "black")

cbar = fig.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label("Mean Area Under C(τ)", fontsize=11)
ax.set_title("Total Spectral Complexity — Phase × Band\n(averaged across patients)",
             fontsize=14, fontweight="bold")
fig.tight_layout()
fig.savefig(OUTPUT_DIR / "fig5_total_complexity_heatmap.png", dpi=200, bbox_inches="tight")
plt.close(fig)
print("  Saved fig5_total_complexity_heatmap.png")


# ---------- Statistical summary and FINDINGS.md ----------
print("\n" + "=" * 60)
print("ANALYSIS SUMMARY")
print("=" * 60)

# 1. Rest vs task peak counts
rest_df = df[df["condition"] == "rest"]
task_df = df[df["condition"] == "task"]
print(f"\n1. Number of natural scales (peaks):")
print(f"   Rest:  mean={rest_df['n_peaks'].mean():.2f} ± {rest_df['n_peaks'].std():.2f}")
print(f"   Task:  mean={task_df['n_peaks'].mean():.2f} ± {task_df['n_peaks'].std():.2f}")
stat, pval = mannwhitneyu(rest_df["n_peaks"], task_df["n_peaks"], alternative="two-sided")
print(f"   Mann-Whitney U test: U={stat:.1f}, p={pval:.4f}")

# 2. Main peak position consistency
print(f"\n2. Main peak τ position (cross-phase patients):")
for band in BANDS:
    band_df = df_cross[df_cross["band"] == band]
    for phase in PHASES:
        vals = band_df[band_df["phase"] == phase]["main_peak_tau"].dropna()
        if len(vals) > 0:
            print(f"   {band:12s} {phase:12s}: τ={vals.mean():.3f} ± {vals.std():.3f}")

# 3. Band-specific effects
print(f"\n3. Band-specific peak count changes (rest→task):")
for band in BANDS:
    r = df_cross[(df_cross["band"] == band) & (df_cross["condition"] == "rest")]["n_peaks"]
    t = df_cross[(df_cross["band"] == band) & (df_cross["condition"] == "task")]["n_peaks"]
    if len(r) > 0 and len(t) > 0:
        diff = t.mean() - r.mean()
        print(f"   {band:12s}: rest={r.mean():.1f}, task={t.mean():.1f}, delta={diff:+.1f}")

# 4. Total complexity
print(f"\n4. Total spectral complexity (AUC):")
for phase in PHASES:
    vals = df_cross[df_cross["phase"] == phase]["area_under_C"]
    print(f"   {phase:12s}: {vals.mean():.3f} ± {vals.std():.3f}")

# 5. Peak height
print(f"\n5. Dominant peak height by phase:")
for phase in PHASES:
    vals = df_cross[df_cross["phase"] == phase]["main_peak_height"].dropna()
    print(f"   {phase:12s}: {vals.mean():.4f} ± {vals.std():.4f}")

# 6. FWHM
print(f"\n6. Main peak FWHM by phase:")
for phase in PHASES:
    vals = df_cross[df_cross["phase"] == phase]["main_peak_fwhm"].dropna()
    print(f"   {phase:12s}: {vals.mean():.3f} ± {vals.std():.3f}")

# 7. Per-patient consistency
print(f"\n7. Per-patient main peak τ (averaged over bands):")
for patient in CROSS_PHASE_PATIENTS:
    for phase in PHASES:
        vals = df[(df["patient"] == patient) & (df["phase"] == phase)]["main_peak_tau"].dropna()
        print(f"   {patient} {phase:12s}: τ={vals.mean():.3f} ± {vals.std():.3f}  (n={len(vals)})")

# ---------- Write FINDINGS.md ----------
findings_lines = []
findings_lines.append("# Spectral Complexity Peaks as Natural Scale Markers\n")
findings_lines.append(f"**Date:** 2026-03-10  ")
findings_lines.append(f"**FC method:** {FC_METHOD}  ")
findings_lines.append(f"**Patients (cross-phase):** {', '.join(CROSS_PHASE_PATIENTS)}  ")
findings_lines.append(f"**All patients:** {', '.join(ALL_PATIENTS)}  ")
findings_lines.append(f"**Peak detection:** prominence>{PROMINENCE_THRESHOLD}, height>{HEIGHT_THRESHOLD}, distance>{DISTANCE}\n")

findings_lines.append("## Key Findings\n")

# Finding 1
rest_mean = rest_df['n_peaks'].mean()
task_mean = task_df['n_peaks'].mean()
findings_lines.append(f"### 1. Number of Natural Scales (Peaks)\n")
findings_lines.append(f"- **Rest:** {rest_mean:.2f} +/- {rest_df['n_peaks'].std():.2f} peaks")
findings_lines.append(f"- **Task:** {task_mean:.2f} +/- {task_df['n_peaks'].std():.2f} peaks")
findings_lines.append(f"- Mann-Whitney U p={pval:.4f}")
if pval < 0.05:
    direction = "more" if task_mean > rest_mean else "fewer"
    findings_lines.append(f"- **Significant difference**: task conditions show {direction} hierarchical scales than rest.")
else:
    findings_lines.append(f"- No significant difference in number of scales between rest and task.\n")

# Finding 2
findings_lines.append(f"\n### 2. Main Peak Position Consistency\n")
for band in BANDS:
    taus_by_phase = {}
    for phase in PHASES:
        vals = df_cross[(df_cross["band"] == band) & (df_cross["phase"] == phase)]["main_peak_tau"].dropna()
        if len(vals) > 0:
            taus_by_phase[phase] = (vals.mean(), vals.std())
    if taus_by_phase:
        mean_vals = [v[0] for v in taus_by_phase.values()]
        spread = max(mean_vals) - min(mean_vals)
        findings_lines.append(f"- **{band}**: τ range across phases = {spread:.3f}")

# Finding 3
findings_lines.append(f"\n### 3. Band-Specific Effects\n")
for band in BANDS:
    r = df_cross[(df_cross["band"] == band) & (df_cross["condition"] == "rest")]["n_peaks"]
    t = df_cross[(df_cross["band"] == band) & (df_cross["condition"] == "task")]["n_peaks"]
    if len(r) > 0 and len(t) > 0:
        diff = t.mean() - r.mean()
        findings_lines.append(f"- **{band}**: delta peaks = {diff:+.2f} (rest={r.mean():.1f}, task={t.mean():.1f})")

# Finding 4
findings_lines.append(f"\n### 4. Total Spectral Complexity (AUC)\n")
for phase in PHASES:
    vals = df_cross[df_cross["phase"] == phase]["area_under_C"]
    findings_lines.append(f"- **{PHASE_NICE[phase]}**: {vals.mean():.3f} +/- {vals.std():.3f}")

# Finding 5
findings_lines.append(f"\n### 5. Dominant Peak Height\n")
for phase in PHASES:
    vals = df_cross[df_cross["phase"] == phase]["main_peak_height"].dropna()
    findings_lines.append(f"- **{PHASE_NICE[phase]}**: {vals.mean():.4f} +/- {vals.std():.4f}")

# Finding 6
findings_lines.append(f"\n### 6. Main Peak Width (FWHM)\n")
for phase in PHASES:
    vals = df_cross[df_cross["phase"] == phase]["main_peak_fwhm"].dropna()
    findings_lines.append(f"- **{PHASE_NICE[phase]}**: {vals.mean():.3f} +/- {vals.std():.3f}")

findings_lines.append(f"\n## Figures\n")
findings_lines.append("1. `fig1_ctau_curves_pat02.png` — C(τ) curves for Pat_02, all phases, one subplot per band")
findings_lines.append("2. `fig2_peaks_rest_vs_task.png` — Peak count boxplot: rest vs task")
findings_lines.append("3. `fig3_main_peak_tau_phases.png` — Main peak τ position across phases per band")
findings_lines.append("4. `fig4_peak_height_heatmap.png` — Dominant peak height heatmap (phase × band)")
findings_lines.append("5. `fig5_total_complexity_heatmap.png` — Total spectral complexity heatmap (phase × band)")

findings_lines.append(f"\n## Data\n")
findings_lines.append(f"- `complexity_peaks_results.csv` — Full results table ({len(df)} rows)")

findings_text = "\n".join(findings_lines) + "\n"
(OUTPUT_DIR / "FINDINGS.md").write_text(findings_text)
print(f"\nFINDINGS.md saved to {OUTPUT_DIR / 'FINDINGS.md'}")
print("\nDone!")
