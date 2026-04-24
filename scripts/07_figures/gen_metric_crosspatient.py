#!/usr/bin/env python3
"""Section 5.2b — Cross-patient metric consistency analysis (beta band).

For each of the 9 reorganization metrics, computes the phase-pair distance
vector for every patient, then measures cross-patient agreement via pairwise
Spearman rank correlation. Identifies which metric(s) produce the most
reproducible reorganization signal.

Outputs:
  - fig_metric_crosspatient_beta.pdf   (consistency bar + patient profiles)
  - metric_crosspatient_beta.csv       (raw values)

Run: python scripts/gen_metric_crosspatient.py
"""
from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT

# ── Config ────────────────────────────────────────────────────────────
BAND = "beta"
BAND_TEX = BRAIN_BAND_TEX_DICT[BAND]

from lrg_eegfc.config.paths import METRIC_CONCORDANCE_CACHE, FIGURES_ROOT
DATA_DIR = METRIC_CONCORDANCE_CACHE
OUTPUT_DIR = FIGURES_ROOT / "metric_concordance"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

METRIC_DISPLAY = {
    "matrix_distance":      "Frobenius",
    "scaled_distance":      "Scaled (log)",
    "rank_distance":        "Rank (1−Spearman)",
    "quantile_rmse":        "Quantile RMSE",
    "permutation_robust":   "Permutation-robust",
    "tree_robinson_foulds": "Robinson−Foulds",
    "tree_cophenetic_corr": "Cophenetic corr.",
    "tree_baker_gamma":     "Baker gamma",
    "tree_fowlkes_mallows": "Fowlkes−Mallows",
}

PAIR_SHORT = {
    "rsPre_vs_taskLearn": "rs₁–tL",
    "rsPre_vs_taskTest":  "rs₁–tT",
    "rsPre_vs_rsPost":    "rs₁–rs₂",
    "taskLearn_vs_taskTest":  "tL–tT",
    "taskLearn_vs_rsPost":    "tL–rs₂",
    "taskTest_vs_rsPost":     "tT–rs₂",
}

PATIENT_COLORS = {
    "Pat_02": "#e6194b",
    "Pat_03": "#3cb44b",
    "Pat_05": "#4363d8",
    "Pat_07": "#f58231",
    "Pat_08": "#911eb4",
}

# ── Load data ─────────────────────────────────────────────────────────
df_all = pd.read_csv(DATA_DIR / "metric_concordance_all_patients.csv")
df = df_all[df_all["band"] == BAND].copy()

# Only patients with all 6 phase pairs (4 phases)
pair_counts = df.groupby("patient")["phase_pair"].nunique()
patients_full = sorted(pair_counts[pair_counts == 6].index)
df = df[df["patient"].isin(patients_full)]

metrics = list(METRIC_DISPLAY.keys())
phase_pairs = sorted(df["phase_pair"].unique())
n_metrics = len(metrics)
n_pairs = len(phase_pairs)
n_patients = len(patients_full)

print(f"Band: {BAND}, {n_patients} patients with 4 phases: {patients_full}")
print(f"{n_metrics} metrics, {n_pairs} phase pairs\n")

# ── Compute cross-patient Spearman for each metric ────────────────────
consistency = {}
for m in metrics:
    sub = df[df["metric"] == m]
    pivot = sub.pivot_table(index="patient", columns="phase_pair", values="value")
    pivot = pivot[phase_pairs]

    rhos = []
    for pi, pj in combinations(pivot.index, 2):
        vi = pivot.loc[pi].values
        vj = pivot.loc[pj].values
        valid = np.isfinite(vi) & np.isfinite(vj)
        if valid.sum() >= 3:
            rho, _ = spearmanr(vi[valid], vj[valid])
            rhos.append(rho)

    consistency[m] = {
        "mean_rS": np.mean(rhos) if rhos else np.nan,
        "min_rS": np.min(rhos) if rhos else np.nan,
        "max_rS": np.max(rhos) if rhos else np.nan,
        "std_rS": np.std(rhos) if rhos else np.nan,
        "n_comparisons": len(rhos),
    }
    print(f"  {METRIC_DISPLAY[m]:22s}  mean rS={consistency[m]['mean_rS']:+.3f}  "
          f"[{consistency[m]['min_rS']:+.3f}, {consistency[m]['max_rS']:+.3f}]")

# ── Figure: top panel = consistency bar, bottom = patient profiles ────
print("\nGenerating figure...")

fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(2, 1, height_ratios=[0.35, 1], hspace=0.35)

# --- Top panel: consistency bar chart ---
ax_bar = fig.add_subplot(gs[0])

mean_rs = [consistency[m]["mean_rS"] for m in metrics]
min_rs = [consistency[m]["min_rS"] for m in metrics]
max_rs = [consistency[m]["max_rS"] for m in metrics]

# Sort by mean consistency
sort_idx = np.argsort(mean_rs)[::-1]
metrics_sorted = [metrics[i] for i in sort_idx]
mean_sorted = [mean_rs[i] for i in sort_idx]
min_sorted = [min_rs[i] for i in sort_idx]
max_sorted = [max_rs[i] for i in sort_idx]

x_bar = np.arange(n_metrics)
bar_colors = ["#4363d8" if v > 0.2 else "#aaa" for v in mean_sorted]
bars = ax_bar.bar(x_bar, mean_sorted, color=bar_colors, edgecolor="white",
                  linewidth=0.5, zorder=3)

# Error bars showing min-max range
for i in range(n_metrics):
    ax_bar.plot([i, i], [min_sorted[i], max_sorted[i]],
               color="0.3", linewidth=1.5, zorder=4)

ax_bar.set_xticks(x_bar)
ax_bar.set_xticklabels([METRIC_DISPLAY[m] for m in metrics_sorted],
                        rotation=35, ha="right", fontsize=8)
ax_bar.set_ylabel("Mean pairwise\nSpearman $r_S$", fontsize=10)
ax_bar.set_title(f"Cross-patient consistency ({BAND_TEX} band, {n_patients} patients)",
                 fontsize=11, fontweight="bold")
ax_bar.axhline(0, color="0.5", linewidth=0.8, zorder=1)
ax_bar.axhline(0.2, color="0.7", linewidth=0.8, linestyle="--", zorder=1)
ax_bar.set_ylim(-0.6, 1.0)
ax_bar.grid(axis="y", alpha=0.3)

# Annotate values
for i, v in enumerate(mean_sorted):
    ax_bar.text(i, v + 0.04, f"{v:+.2f}", ha="center", va="bottom",
               fontsize=7, fontweight="bold")

# --- Bottom panel: 3x3 grid of patient profiles ---
gs_inner = gs[1].subgridspec(3, 3, hspace=0.45, wspace=0.3)

for idx, m in enumerate(metrics_sorted):
    ax = fig.add_subplot(gs_inner[idx // 3, idx % 3])
    sub = df[df["metric"] == m]

    for patient in patients_full:
        psub = sub[sub["patient"] == patient].set_index("phase_pair")
        vals = [psub.loc[pp, "value"] if pp in psub.index else np.nan
                for pp in phase_pairs]
        # Normalize to [0, 1] range for comparison across patients
        vmin, vmax = np.nanmin(vals), np.nanmax(vals)
        if vmax > vmin:
            vals_norm = [(v - vmin) / (vmax - vmin) for v in vals]
        else:
            vals_norm = [0.5] * len(vals)

        ax.plot(range(n_pairs), vals_norm,
                "-o", color=PATIENT_COLORS[patient], markersize=3.5,
                linewidth=1.2, alpha=0.85, label=patient)

    ax.set_xticks(range(n_pairs))
    ax.set_xticklabels([PAIR_SHORT.get(pp, pp) for pp in phase_pairs],
                        fontsize=6, rotation=30, ha="right")
    ax.set_title(METRIC_DISPLAY[m], fontsize=8, fontweight="bold")
    ax.set_ylabel("norm. dist.", fontsize=6)
    ax.tick_params(labelsize=6)
    ax.set_ylim(-0.1, 1.1)
    ax.grid(axis="y", alpha=0.2)

    # Show mean rS in corner
    rs = consistency[m]["mean_rS"]
    ax.text(0.97, 0.95, f"$\\bar{{r}}_S$={rs:+.2f}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=7, bbox=dict(boxstyle="round,pad=0.2",
                                   facecolor="white", alpha=0.8))

# Legend
handles = [plt.Line2D([0], [0], color=PATIENT_COLORS[p], marker="o",
                       markersize=4, linewidth=1.2, label=p)
           for p in patients_full]
fig.legend(handles=handles, loc="lower center", ncol=n_patients,
           fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.01))

fig_path = OUTPUT_DIR / f"fig_metric_crosspatient_{BAND}.pdf"
fig.savefig(fig_path, bbox_inches="tight", dpi=300)
plt.close(fig)
print(f"  Saved: {fig_path}")

# ── Save CSV ──────────────────────────────────────────────────────────
rows = []
for m in metrics:
    c = consistency[m]
    rows.append({
        "metric": m,
        "display_name": METRIC_DISPLAY[m],
        "band": BAND,
        "n_patients": n_patients,
        "mean_spearman": c["mean_rS"],
        "min_spearman": c["min_rS"],
        "max_spearman": c["max_rS"],
        "std_spearman": c["std_rS"],
    })
csv_path = DATA_DIR / f"metric_crosspatient_{BAND}.csv"
pd.DataFrame(rows).to_csv(csv_path, index=False)
print(f"  Saved: {csv_path}")

print(f"\nBest metric: {METRIC_DISPLAY[metrics_sorted[0]]} "
      f"(mean rS = {mean_sorted[0]:+.3f})")
