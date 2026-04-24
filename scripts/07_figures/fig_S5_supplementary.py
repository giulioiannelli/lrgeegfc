#!/usr/bin/env python3
"""Supplementary figures for Section 5: Selecting distance measures.

Produces:
  A1. Metric comparison summary table (CSV) — all 14 metrics × H1/H2a/H2b/H3
  A2. Extended unanimity heatmap (FIGURE)
  A3. Cophenetic correlation vs mean VI scatter (FIGURE)
  B1. Cophenetic phase-distance matrices for Pat_02 and Pat_08 (FIGURES)
  B2. Per-patient NVI profiles for Pat_02 + mean 4-patient (FIGURES)
  D2. Outlier summary table (TEXT)
  E1. Key numbers for LaTeX (TEXT)
  F1. Three-tool example for Pat_02 beta (FIGURE)
  F2. Cophenetic Pre-Post heatmap all patients (FIGURE)

Output: data/figures/metric_exploration/report_figures/
"""
from __future__ import annotations

from pathlib import Path
from itertools import combinations

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, cophenet, leaves_list
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr

from lrg_eegfc.config.paths import LRG_CACHE, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result

OUTDIR = FIGURES_ROOT / "metric_exploration" / "report_figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

PATIENTS_4PH = ["Pat_02", "Pat_03", "Pat_05", "Pat_08"]
ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]
PATIENT_PHASES = {
    "Pat_02": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_03": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_05": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "rest_post"],
    "Pat_08": ["rest_pre", "task_learn", "task_test", "rest_post"],
}
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = {
    "delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
    "beta": r"$\beta$", "low_gamma": r"$\gamma_l$", "high_gamma": r"$\gamma_h$",
}
PHASE_SHORT = {"rest_pre": "Pre", "task_learn": "TL", "task_test": "TT", "rest_post": "Post"}
ALL_PAIRS = [
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"), ("rest_pre", "rest_post"),
    ("task_learn", "task_test"), ("task_learn", "rest_post"), ("task_test", "rest_post"),
]
WITHIN_PAIRS = [("task_learn", "task_test"), ("rest_pre", "rest_post")]
CROSS_PAIRS = [("rest_pre", "task_learn"), ("rest_pre", "task_test"),
               ("task_learn", "rest_post"), ("task_test", "rest_post")]


from lrg_eegfc.utils.metrics import compute_vi  # canonical implementation


def cophenetic_corr(Z1, Z2):
    """Pearson correlation between cophenetic distance matrices."""
    d1 = cophenet(Z1)
    d2 = cophenet(Z2)
    return pearsonr(d1, d2)[0]


# ── Load linkage matrices ────────────────────────────────────────────
print("Loading LRG results...")
linkages = {}
for pat in ALL_PATIENTS:
    for phase in PATIENT_PHASES[pat]:
        for band in BANDS:
            try:
                res = load_lrg_result(pat, phase, band, fc_method="msc",
                                      cache_root=LRG_CACHE)
                linkages[(pat, band, phase)] = res.linkage_matrix
            except Exception:
                pass
print(f"Loaded {len(linkages)} linkages")


# ═══════════════════════════════════════════════════════════════════════
# A1 + A2: Extended metric comparison with H1 and H3
# ═══════════════════════════════════════════════════════════════════════
print("\n=== BLOCK A: Metric comparison ===")

# Load existing pairwise results
pw_csv = FIGURES_ROOT / "metric_exploration" / "alternative_metrics" / "pairwise_results.csv"
df = pd.read_csv(pw_csv)

METRICS = [c for c in df.columns if c not in
           ["patient", "band", "phase1", "phase2", "pair"]]
# Identify which metrics are similarity (higher = more similar) vs distance
SIM_METRICS = {"CophPearson", "CophSpearman", "BakersGamma",
               "TopK5_Merge", "TopK10_Merge", "TopK20_Merge",
               "WeightARI_Coarse", "WeightARI_Fine", "MeanARI"}
DIST_METRICS = {"NormL1_Ult", "NormL2_Ult", "WeightVI_Coarse",
                "WeightVI_Fine", "MeanVI"}


def compute_unanimity_table():
    """For each metric, compute H1/H2a/H2b/H3 unanimity per band."""
    results = []

    for metric in METRICS:
        is_sim = metric in SIM_METRICS
        for band in BANDS:
            band_df = df[df["band"] == band]

            h1_signs = []
            h2a_signs = []
            h2b_signs = []
            h3_signs = []

            for pat in PATIENTS_4PH:
                pat_df = band_df[band_df["patient"] == pat]
                if pat_df.empty:
                    continue

                # Get pairwise values
                vals = {}
                for _, row in pat_df.iterrows():
                    vals[(row["phase1"], row["phase2"])] = row[metric]

                # Find TL-TT value
                tl_tt = vals.get(("task_learn", "task_test"), np.nan)
                if np.isnan(tl_tt):
                    continue

                # H1: Is TL-TT the most similar pair?
                # For similarity metrics: TL-TT should be largest
                # For distance metrics: TL-TT should be smallest
                others = [v for k, v in vals.items()
                          if k != ("task_learn", "task_test") and not np.isnan(v)]
                if others:
                    mean_others = np.mean(others)
                    if is_sim:
                        h1_signs.append(np.sign(tl_tt - mean_others))
                    else:
                        h1_signs.append(np.sign(mean_others - tl_tt))

                # H2a: Is rest_post more similar to TT than to Pre?
                pre_post = vals.get(("rest_pre", "rest_post"), np.nan)
                tt_post = vals.get(("task_test", "rest_post"), np.nan)
                if not (np.isnan(pre_post) or np.isnan(tt_post)):
                    if is_sim:
                        h2a_signs.append(np.sign(tt_post - pre_post))
                    else:
                        h2a_signs.append(np.sign(pre_post - tt_post))

                # H2b: Is rest_post closer to TT than Pre was?
                pre_tt = vals.get(("rest_pre", "task_test"), np.nan)
                if not (np.isnan(pre_tt) or np.isnan(tt_post)):
                    if is_sim:
                        h2b_signs.append(np.sign(tt_post - pre_tt))
                    else:
                        h2b_signs.append(np.sign(pre_tt - tt_post))

                # H3: Within-type more similar than cross-type
                within_vals = [vals.get(p, np.nan) for p in WITHIN_PAIRS]
                cross_vals = [vals.get(p, np.nan) for p in CROSS_PAIRS]
                within_vals = [v for v in within_vals if not np.isnan(v)]
                cross_vals = [v for v in cross_vals if not np.isnan(v)]
                if within_vals and cross_vals:
                    if is_sim:
                        h3_signs.append(np.sign(
                            np.mean(within_vals) - np.mean(cross_vals)))
                    else:
                        h3_signs.append(np.sign(
                            np.mean(cross_vals) - np.mean(within_vals)))

            # Unanimity: all patients agree
            h1_unan = (len(h1_signs) >= 4 and
                        abs(sum(h1_signs)) == len(h1_signs))
            h2a_unan = (len(h2a_signs) >= 4 and
                         abs(sum(h2a_signs)) == len(h2a_signs) and
                         sum(h2a_signs) > 0)
            h2b_unan = (len(h2b_signs) >= 4 and
                         abs(sum(h2b_signs)) == len(h2b_signs) and
                         sum(h2b_signs) > 0)
            h3_unan = (len(h3_signs) >= 4 and
                        abs(sum(h3_signs)) == len(h3_signs))

            results.append({
                "metric": metric, "band": band,
                "H1_unanimous": h1_unan,
                "H2a_unanimous": h2a_unan,
                "H2b_unanimous": h2b_unan,
                "H3_unanimous": h3_unan,
            })
    return pd.DataFrame(results)


print("Computing unanimity for all metrics and hypotheses...")
unan_df = compute_unanimity_table()

# Summarize: count unanimous bands per (metric, hypothesis)
summary_rows = []
for metric in METRICS:
    mdf = unan_df[unan_df["metric"] == metric]
    h1_count = mdf["H1_unanimous"].sum()
    h2a_count = mdf["H2a_unanimous"].sum()
    h2b_count = mdf["H2b_unanimous"].sum()
    h3_count = mdf["H3_unanimous"].sum()
    metric_type = "similarity" if metric in SIM_METRICS else "distance"
    total = h1_count + h2a_count + h2b_count + h3_count
    summary_rows.append({
        "Metric": metric,
        "Type": metric_type,
        "H1_unan": int(h1_count),
        "H2a_unan": int(h2a_count),
        "H2b_unan": int(h2b_count),
        "H3_unan": int(h3_count),
        "Total": int(total),
    })

summary_df = pd.DataFrame(summary_rows).sort_values("Total", ascending=False)

# A1: Save as CSV
a1_path = OUTDIR / "table_A1_metric_comparison.csv"
summary_df.to_csv(a1_path, index=False)
print(f"Saved {a1_path.name}")
print(summary_df.to_string(index=False))

# A2: Unanimity heatmap figure
print("\nFigure A2: Extended unanimity heatmap...")

fig, ax = plt.subplots(figsize=(8, 10))

metrics_sorted = summary_df["Metric"].tolist()
hypo_cols = ["H1_unan", "H2a_unan", "H2b_unan", "H3_unan"]
mat = summary_df[hypo_cols].values.astype(float)

im = ax.imshow(mat, cmap="YlOrRd", vmin=0, vmax=6, aspect="auto")

# Annotate
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        v = int(mat[i, j])
        color = "white" if v >= 4 else "black"
        ax.text(j, i, str(v), ha="center", va="center", fontsize=12,
                fontweight="bold", color=color)

ax.set_xticks(range(4))
ax.set_xticklabels(["H1\n(task stability)", "H2a\n(task persist)",
                     "H2b\n(task proximity)", "H3\n(within<cross)"],
                    fontsize=10)
ax.set_yticks(range(len(metrics_sorted)))
ax.set_yticklabels(metrics_sorted, fontsize=10)

# Add total column as text
for i, row in enumerate(summary_df.itertuples()):
    ax.text(4.3, i, f"= {row.Total}", fontsize=10, ha="left", va="center",
            fontweight="bold")

ax.set_xlim(-0.5, 3.5)
cbar = fig.colorbar(im, ax=ax, shrink=0.5, pad=0.15)
cbar.set_label("# bands with unanimity (max 6)", fontsize=10)

ax.set_title("Metric comparison: # frequency bands with 4/4 patient unanimity\n"
             "Sorted by total score (right column)", fontsize=12, fontweight="bold")

fig.tight_layout()
fig.savefig(OUTDIR / "fig_S5_metric_unanimity_heatmap.pdf",
            bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved fig_S5_metric_unanimity_heatmap.pdf")


# ═══════════════════════════════════════════════════════════════════════
# A3: Cophenetic correlation vs mean VI scatter
# ═══════════════════════════════════════════════════════════════════════
print("\n=== A3: Cophenetic vs mean VI scatter ===")

coph_vals = []
vi_vals = []
scatter_bands = []
scatter_pats = []

for pat in PATIENTS_4PH:
    for band in BANDS:
        for p1, p2 in ALL_PAIRS:
            k1, k2 = (pat, band, p1), (pat, band, p2)
            if k1 not in linkages or k2 not in linkages:
                continue
            Z1, Z2 = linkages[k1], linkages[k2]
            # Cophenetic correlation
            cc = cophenetic_corr(Z1, Z2)
            # Mean VI (k=2..30)
            vis = []
            for k in range(2, 31):
                l1 = fcluster(Z1, k, criterion="maxclust")
                l2 = fcluster(Z2, k, criterion="maxclust")
                vis.append(compute_vi(l1, l2))
            mvi = np.mean(vis)

            coph_vals.append(cc)
            vi_vals.append(mvi)
            scatter_bands.append(band)
            scatter_pats.append(pat)

print(f"  Computed {len(coph_vals)} points for scatter")

BAND_COLORS = {
    "delta": "#1f77b4", "theta": "#ff7f0e", "alpha": "#2ca02c",
    "beta": "#d62728", "low_gamma": "#9467bd", "high_gamma": "#8c564b",
}
PAT_MARKERS = {"Pat_02": "o", "Pat_03": "s", "Pat_05": "D", "Pat_08": "^"}

fig, ax = plt.subplots(figsize=(8, 6))

for band in BANDS:
    for pat in PATIENTS_4PH:
        mask = [(b == band and p == pat)
                for b, p in zip(scatter_bands, scatter_pats)]
        x = [coph_vals[i] for i, m in enumerate(mask) if m]
        y = [vi_vals[i] for i, m in enumerate(mask) if m]
        if x:
            ax.scatter(x, y, c=BAND_COLORS[band], marker=PAT_MARKERS[pat],
                       s=40, alpha=0.7, edgecolors="none")

# Legend entries
for band in BANDS:
    ax.scatter([], [], c=BAND_COLORS[band], s=60, label=BAND_TEX[band])
for pat in PATIENTS_4PH:
    ax.scatter([], [], c="gray", marker=PAT_MARKERS[pat], s=60, label=pat)

# Regression line
x_arr = np.array(coph_vals)
y_arr = np.array(vi_vals)
valid = np.isfinite(x_arr) & np.isfinite(y_arr)
if valid.sum() > 2:
    z = np.polyfit(x_arr[valid], y_arr[valid], 1)
    x_line = np.linspace(x_arr[valid].min(), x_arr[valid].max(), 100)
    ax.plot(x_line, np.polyval(z, x_line), "k--", linewidth=1.5, alpha=0.5)
    r, p = pearsonr(x_arr[valid], y_arr[valid])
    ax.text(0.05, 0.95, f"r = {r:.3f}, p < {max(p, 1e-50):.1e}",
            transform=ax.transAxes, fontsize=10, va="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

ax.set_xlabel("Cophenetic Pearson correlation (higher = more similar)", fontsize=11)
ax.set_ylabel("Mean VI (k=2..30) (higher = more different)", fontsize=11)
ax.set_title("Cophenetic correlation vs Mean VI\n"
             "Each point = one (patient, band, phase-pair) triplet", fontsize=12)
ax.legend(fontsize=8, ncol=2, loc="upper right")
ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(OUTDIR / "fig_cophenetic_vs_meanVI_scatter.pdf",
            bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved fig_cophenetic_vs_meanVI_scatter.pdf")


# ═══════════════════════════════════════════════════════════════════════
# B1: Cophenetic phase-distance matrices for Pat_02 and Pat_08
# ═══════════════════════════════════════════════════════════════════════
print("\n=== B1: Cophenetic phase matrices ===")

for pat in ["Pat_02", "Pat_08"]:
    phases = PATIENT_PHASES[pat]
    fig, axes = plt.subplots(1, len(BANDS), figsize=(24, 4.5))

    for ib, band in enumerate(BANDS):
        ax = axes[ib]
        n_ph = len(phases)
        mat = np.full((n_ph, n_ph), np.nan)

        for i, p1 in enumerate(phases):
            for j, p2 in enumerate(phases):
                if i == j:
                    mat[i, j] = 1.0
                    continue
                k1, k2 = (pat, band, p1), (pat, band, p2)
                if k1 in linkages and k2 in linkages:
                    mat[i, j] = cophenetic_corr(linkages[k1], linkages[k2])

        # Determine color range
        valid_vals = mat[~np.isnan(mat)]
        if len(valid_vals) > 0:
            vmin = max(valid_vals.min() - 0.05, -1)
            vmax = 1.0
        else:
            vmin, vmax = 0, 1

        im = ax.imshow(mat, cmap="RdBu_r", vmin=vmin, vmax=vmax, aspect="equal")

        # Annotate
        for i in range(n_ph):
            for j in range(n_ph):
                if not np.isnan(mat[i, j]):
                    color = "white" if abs(mat[i, j]) > 0.8 else "black"
                    ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                            fontsize=8, fontweight="bold", color=color)

        ax.set_xticks(range(n_ph))
        ax.set_xticklabels([PHASE_SHORT[p] for p in phases], fontsize=9)
        ax.set_yticks(range(n_ph))
        ax.set_yticklabels([PHASE_SHORT[p] for p in phases] if ib == 0
                           else [], fontsize=9)
        ax.set_title(BAND_TEX[band], fontsize=14, fontweight="bold")

    fig.colorbar(im, ax=axes, shrink=0.8, pad=0.02)
    fig.suptitle(
        f"{pat} — Cophenetic Pearson correlation between phase pairs\n"
        f"1.0 = identical hierarchy | Lower = more reorganized",
        fontsize=12, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(OUTDIR / f"fig_cophenetic_phase_matrix_{pat}.pdf",
                bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Saved fig_cophenetic_phase_matrix_{pat}.pdf")


# ═══════════════════════════════════════════════════════════════════════
# B2: NVI profiles — Pat_02 only + mean 4 patients
# ═══════════════════════════════════════════════════════════════════════
print("\n=== B2: NVI profiles ===")

H_GRID = np.geomspace(0.003, 0.995, 400)

# Compute k(h) for normalization
print("Computing k(h) and VI(h)...")
k_profiles = {}
for key, Z in linkages.items():
    ks = np.array([len(np.unique(fcluster(Z, h, criterion="distance")))
                    for h in H_GRID])
    k_profiles[key] = ks

all_k = np.array([k_profiles[k] for k in k_profiles])
mean_k_at_h = np.mean(all_k, axis=0)
norm_factor = np.log(np.maximum(mean_k_at_h, 1.001))

# Compute VI(h) for all pairs (4-phase patients only)
vi_h = {pat: {band: {} for band in BANDS} for pat in PATIENTS_4PH}

total = sum(len(list(combinations(PATIENT_PHASES[p], 2))) * len(BANDS)
            for p in PATIENTS_4PH)
done = 0
for pat in PATIENTS_4PH:
    phases = PATIENT_PHASES[pat]
    for band in BANDS:
        for p1, p2 in combinations(phases, 2):
            k1, k2 = (pat, band, p1), (pat, band, p2)
            if k1 not in linkages or k2 not in linkages:
                done += 1
                continue
            Z1, Z2 = linkages[k1], linkages[k2]
            vals = np.zeros(len(H_GRID))
            for ih, h in enumerate(H_GRID):
                l1 = fcluster(Z1, h, criterion="distance")
                l2 = fcluster(Z2, h, criterion="distance")
                vals[ih] = compute_vi(l1, l2)
            vi_h[pat][band][(p1, p2)] = vals
            done += 1
    print(f"  {pat} done ({done}/{total})")

PAIR_LABELS = {
    ("task_learn", "task_test"): "TL-TT", ("rest_pre", "rest_post"): "Pre-Post",
    ("rest_pre", "task_learn"): "Pre-TL", ("rest_pre", "task_test"): "Pre-TT",
    ("task_learn", "rest_post"): "TL-Post", ("task_test", "rest_post"): "TT-Post",
}
PAIR_COLORS = {
    ("task_learn", "task_test"): "#1565C0",
    ("rest_pre", "rest_post"): "#2E7D32",
    ("rest_pre", "task_learn"): "#9E9E9E",
    ("rest_pre", "task_test"): "#BDBDBD",
    ("task_learn", "rest_post"): "#E65100",
    ("task_test", "rest_post"): "#C62828",
}

# B2a: Pat_02 only
print("Figure B2a: Pat_02 NVI profiles...")
fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True)
for ib, band in enumerate(BANDS):
    ax = axes.ravel()[ib]
    for pair in ALL_PAIRS:
        v = vi_h["Pat_02"][band].get(pair)
        if v is None:
            continue
        v_norm = v / np.maximum(norm_factor, 0.01)
        ax.plot(H_GRID, v_norm, color=PAIR_COLORS[pair], linewidth=2.5,
                label=PAIR_LABELS[pair])
    ax.set_title(BAND_TEX[band], fontsize=14, fontweight="bold")
    ax.set_xscale("log")
    if ib >= 3:
        ax.set_xlabel("Height h (log scale)")
    if ib % 3 == 0:
        ax.set_ylabel("NVI = VI(h) / ln(k)")
    ax.grid(alpha=0.3)
    if ib == 0:
        ax.legend(fontsize=8, loc="upper right")

fig.suptitle(
    "Pat_02 — Scale-resolved NVI profiles for all phase pairs\n"
    "NVI(h) = VI(h) / ln(k(h))  |  Higher = more different at that scale",
    fontsize=12, y=1.02,
)
fig.tight_layout()
fig.savefig(OUTDIR / "fig_NVI_profiles_Pat02.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved fig_NVI_profiles_Pat02.pdf")

# B2b: Mean of 4 patients with range
print("Figure B2b: Mean 4-patient NVI profiles...")
fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharex=True)
for ib, band in enumerate(BANDS):
    ax = axes.ravel()[ib]
    for pair in ALL_PAIRS:
        curves = []
        for pat in PATIENTS_4PH:
            v = vi_h[pat][band].get(pair)
            if v is not None:
                curves.append(v / np.maximum(norm_factor, 0.01))
        if not curves:
            continue
        arr = np.array(curves)
        mean_c = np.mean(arr, axis=0)
        lo, hi = np.min(arr, axis=0), np.max(arr, axis=0)
        color = PAIR_COLORS[pair]
        ax.plot(H_GRID, mean_c, color=color, linewidth=2.5,
                label=PAIR_LABELS[pair])
        ax.fill_between(H_GRID, lo, hi, color=color, alpha=0.15)

    ax.set_title(BAND_TEX[band], fontsize=14, fontweight="bold")
    ax.set_xscale("log")
    if ib >= 3:
        ax.set_xlabel("Height h (log scale)")
    if ib % 3 == 0:
        ax.set_ylabel("NVI = VI(h) / ln(k)")
    ax.grid(alpha=0.3)
    if ib == 0:
        ax.legend(fontsize=8, loc="upper right")

fig.suptitle(
    "Mean NVI profiles across 4 patients (shading = patient range)\n"
    "NVI(h) = VI(h) / ln(k(h))  |  Higher = more different at that scale",
    fontsize=12, y=1.02,
)
fig.tight_layout()
fig.savefig(OUTDIR / "fig_NVI_profiles_mean_4patients.pdf",
            bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved fig_NVI_profiles_mean_4patients.pdf")


# ═══════════════════════════════════════════════════════════════════════
# F1: Three-tool example for Pat_02 beta (rest_pre vs task_learn)
# ═══════════════════════════════════════════════════════════════════════
print("\n=== F1: Three-tool example ===")

pat_ex, band_ex = "Pat_02", "beta"
p1_ex, p2_ex = "rest_pre", "task_learn"
k1_ex = (pat_ex, band_ex, p1_ex)
k2_ex = (pat_ex, band_ex, p2_ex)

if k1_ex in linkages and k2_ex in linkages:
    Z1_ex, Z2_ex = linkages[k1_ex], linkages[k2_ex]
    cc_ex = cophenetic_corr(Z1_ex, Z2_ex)

    # Compute affinity matrices for this example
    def _compute_aff(Z_, K_max=60):
        n = Z_.shape[0] + 1
        K = min(K_max, n - 1)
        A = np.zeros((n, n))
        for k in range(2, K + 1):
            labels = fcluster(Z_, k, criterion="maxclust")
            same = (labels[:, None] == labels[None, :]).astype(np.float64)
            A += same
        A /= (K - 1)
        return A

    A1_ex = _compute_aff(Z1_ex)
    A2_ex = _compute_aff(Z2_ex)
    order_ex = leaves_list(Z1_ex)

    # NVI curve
    vi_curve_ex = vi_h[pat_ex][band_ex].get((p1_ex, p2_ex))
    nvi_curve_ex = vi_curve_ex / np.maximum(norm_factor, 0.01) if vi_curve_ex is not None else None

    fig = plt.figure(figsize=(20, 5))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1.2, 1, 1])

    # Tool 1: Cophenetic correlation (single number)
    ax1 = fig.add_subplot(gs[0])
    ax1.text(0.5, 0.5, f"r = {cc_ex:.3f}", fontsize=36, ha="center", va="center",
             fontweight="bold", transform=ax1.transAxes,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="#E3F2FD",
                       edgecolor="#1565C0", linewidth=2))
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_title("Cophenetic correlation\n(one number)", fontsize=11,
                   fontweight="bold")
    ax1.text(0.5, 0.15, "Global hierarchy similarity", fontsize=9,
             ha="center", va="center", transform=ax1.transAxes, style="italic")

    # Tool 2: NVI profile (curve)
    ax2 = fig.add_subplot(gs[1])
    if nvi_curve_ex is not None:
        ax2.plot(H_GRID, nvi_curve_ex, color="#1565C0", linewidth=2.5)
        ax2.fill_between(H_GRID, 0, nvi_curve_ex, color="#1565C0", alpha=0.1)
        ax2.set_xscale("log")
        ax2.set_xlabel("Height h (log scale)", fontsize=10)
        ax2.set_ylabel("NVI(h)", fontsize=10)
        ax2.grid(alpha=0.3)
    ax2.set_title("Scale-resolved VI profile\n(one curve)", fontsize=11,
                   fontweight="bold")

    # Tool 3 & 4: Affinity matrices (two heatmaps)
    ax3 = fig.add_subplot(gs[2])
    A1_ord = A1_ex[np.ix_(order_ex, order_ex)]
    ax3.imshow(A1_ord, cmap="magma", vmin=0, vmax=1, aspect="equal")
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.set_title(f"Affinity: {PHASE_SHORT[p1_ex]}\n(node x node)", fontsize=11,
                   fontweight="bold")

    ax4 = fig.add_subplot(gs[3])
    A2_ord = A2_ex[np.ix_(order_ex, order_ex)]
    im4 = ax4.imshow(A2_ord, cmap="magma", vmin=0, vmax=1, aspect="equal")
    ax4.set_xticks([])
    ax4.set_yticks([])
    ax4.set_title(f"Affinity: {PHASE_SHORT[p2_ex]}\n(node x node)", fontsize=11,
                   fontweight="bold")

    fig.suptitle(
        f"{pat_ex}, {BAND_TEX[band_ex]} — Three complementary tools for "
        f"comparing hierarchies ({PHASE_SHORT[p1_ex]} vs {PHASE_SHORT[p2_ex]})",
        fontsize=13, fontweight="bold", y=1.04,
    )
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig_three_tools_example_Pat02_beta.pdf",
                bbox_inches="tight", dpi=150)
    plt.close(fig)
    print("Saved fig_three_tools_example_Pat02_beta.pdf")


# ═══════════════════════════════════════════════════════════════════════
# F2: Cophenetic Pre-Post heatmap for all 6 patients
# ═══════════════════════════════════════════════════════════════════════
print("\n=== F2: Cophenetic Pre-Post heatmap ===")

fig, ax = plt.subplots(figsize=(8, 5))

# Patients on y-axis, bands on x-axis
patients_f2 = [p for p in ALL_PATIENTS
               if "rest_pre" in PATIENT_PHASES[p] and "rest_post" in PATIENT_PHASES[p]]
mat_f2 = np.full((len(patients_f2), len(BANDS)), np.nan)

for ip, pat in enumerate(patients_f2):
    for ib, band in enumerate(BANDS):
        k1 = (pat, band, "rest_pre")
        k2 = (pat, band, "rest_post")
        if k1 in linkages and k2 in linkages:
            mat_f2[ip, ib] = cophenetic_corr(linkages[k1], linkages[k2])

im = ax.imshow(mat_f2, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
for i in range(mat_f2.shape[0]):
    for j in range(mat_f2.shape[1]):
        if not np.isnan(mat_f2[i, j]):
            color = "white" if mat_f2[i, j] < 0.4 else "black"
            ax.text(j, i, f"{mat_f2[i, j]:.2f}", ha="center", va="center",
                    fontsize=11, fontweight="bold", color=color)

ax.set_xticks(range(len(BANDS)))
ax.set_xticklabels([BAND_TEX[b] for b in BANDS], fontsize=12)
ax.set_yticks(range(len(patients_f2)))
ax.set_yticklabels(patients_f2, fontsize=11)
cbar = fig.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label("Cophenetic Pearson correlation (rest_pre vs rest_post)", fontsize=10)
ax.set_title("Hierarchical similarity between rest_pre and rest_post\n"
             "Lower = more reorganization across the experiment",
             fontsize=12, fontweight="bold")

fig.tight_layout()
fig.savefig(OUTDIR / "fig_cophenetic_prepost_all_patients.pdf",
            bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved fig_cophenetic_prepost_all_patients.pdf")


# ═══════════════════════════════════════════════════════════════════════
# D2: Outlier summary table
# ═══════════════════════════════════════════════════════════════════════
print("\n=== D2: Outlier summary table ===")

d2_text = """# Outlier Summary Table
# For use in LaTeX

Patient | Band       | Outlier type                  | VI value | Expected range | Explanation
--------|------------|-------------------------------|----------|----------------|--------------------------------------------------
Pat_02  | theta      | TL-TT worst pair (not best)   | 0.877    | 0.1-0.5        | Coarse-scale structural flip between task phases; fine scale (k>10) reverts to expected pattern
Pat_03  | high_gamma | Extreme Pre-Post distance      | 0.964    | 0.3-0.7        | Star-like hub topology in rest_pre transitions to modular structure in rest_post; genuine biological reorganization
Pat_08  | delta      | Ultra-low Pre-Post distance    | 0.144    | 0.3-0.7        | Delta resting-state hierarchy is ultra-stable; two rest recordings produce nearly identical dendrograms
"""

d2_path = OUTDIR / "table_D2_outlier_summary.txt"
d2_path.write_text(d2_text)
print(f"Saved {d2_path.name}")


# ═══════════════════════════════════════════════════════════════════════
# E1: Key numbers for LaTeX
# ═══════════════════════════════════════════════════════════════════════
print("\n=== E1: Key numbers ===")

# 1. Mean cophenetic correlation within-type vs cross-type per band
e1_lines = ["# KEY NUMBERS FOR LATEX",
            "# Computed from LRG dendrograms (MSC FC, nperseg=4096)",
            ""]
e1_lines.append("## 1. Cophenetic correlation: within-type vs cross-type")

for band in BANDS:
    within_cc = []
    cross_cc = []
    for pat in PATIENTS_4PH:
        for p1, p2 in WITHIN_PAIRS:
            k1, k2 = (pat, band, p1), (pat, band, p2)
            if k1 in linkages and k2 in linkages:
                within_cc.append(cophenetic_corr(linkages[k1], linkages[k2]))
        for p1, p2 in CROSS_PAIRS:
            k1, k2 = (pat, band, p1), (pat, band, p2)
            if k1 in linkages and k2 in linkages:
                cross_cc.append(cophenetic_corr(linkages[k1], linkages[k2]))
    e1_lines.append(
        f"  {BAND_TEX[band]:>12s}: within={np.mean(within_cc):.3f} +/- {np.std(within_cc):.3f}, "
        f"cross={np.mean(cross_cc):.3f} +/- {np.std(cross_cc):.3f}")

# Grand mean
all_within = []
all_cross = []
for band in BANDS:
    for pat in PATIENTS_4PH:
        for p1, p2 in WITHIN_PAIRS:
            k1, k2 = (pat, band, p1), (pat, band, p2)
            if k1 in linkages and k2 in linkages:
                all_within.append(cophenetic_corr(linkages[k1], linkages[k2]))
        for p1, p2 in CROSS_PAIRS:
            k1, k2 = (pat, band, p1), (pat, band, p2)
            if k1 in linkages and k2 in linkages:
                all_cross.append(cophenetic_corr(linkages[k1], linkages[k2]))
e1_lines.append(
    f"  GRAND MEAN: within={np.mean(all_within):.3f} +/- {np.std(all_within):.3f}, "
    f"cross={np.mean(all_cross):.3f} +/- {np.std(all_cross):.3f}")

# 2. H2 exact h-ranges (from continuous_all_pairs output)
e1_lines.append("")
e1_lines.append("## 2. H2 unanimous regions (exact h-ranges)")

# Recompute from VI(h) data
def find_unanimous_regions(vi_h_data, pair_a, pair_b, direction="+"):
    """Find h-ranges where all patients agree on sign of (pair_a - pair_b)."""
    contrasts = []
    for pat in PATIENTS_4PH:
        va = vi_h_data[pat].get(pair_a)
        vb = vi_h_data[pat].get(pair_b)
        if va is not None and vb is not None:
            contrasts.append(va - vb)
    if not contrasts:
        return []
    arr = np.array(contrasts)
    signs = np.sign(arr)
    unan = np.mean(signs, axis=0)
    if direction == "+":
        mask = unan >= 0.99
    else:
        mask = unan <= -0.99
    # Find contiguous regions
    regions = []
    if not np.any(mask):
        return regions
    starts = np.where(np.diff(np.concatenate(([0], mask.astype(int)))) == 1)[0]
    ends = np.where(np.diff(np.concatenate((mask.astype(int), [0]))) == -1)[0]
    for s, e in zip(starts, ends):
        if e > s:
            regions.append((H_GRID[s], H_GRID[min(e, len(H_GRID)-1)],
                           mean_k_at_h[min(e, len(H_GRID)-1)],
                           mean_k_at_h[s]))
    return regions

for band in BANDS:
    band_vi = {pat: vi_h[pat][band] for pat in PATIENTS_4PH}
    persist = find_unanimous_regions(band_vi,
                                     ("rest_pre", "rest_post"), ("task_test", "rest_post"),
                                     direction="+")
    recover = find_unanimous_regions(band_vi,
                                     ("rest_pre", "rest_post"), ("task_test", "rest_post"),
                                     direction="-")
    parts = []
    for (h_lo, h_hi, k_lo, k_hi) in persist:
        parts.append(f"PERSIST h=[{h_lo:.4f}, {h_hi:.4f}] k=[{k_lo:.1f}, {k_hi:.1f}]")
    for (h_lo, h_hi, k_lo, k_hi) in recover:
        parts.append(f"RECOVER h=[{h_lo:.4f}, {h_hi:.4f}] k=[{k_lo:.1f}, {k_hi:.1f}]")
    e1_lines.append(f"  {BAND_TEX[band]:>12s}: {' | '.join(parts) if parts else 'none'}")

# 3. Total comparisons
n_pairs_4ph = len(PATIENTS_4PH) * len(BANDS) * len(ALL_PAIRS)
n_total = sum(len(list(combinations(PATIENT_PHASES[p], 2))) * len(BANDS)
              for p in ALL_PATIENTS)
e1_lines.append("")
e1_lines.append("## 3. Total comparisons")
e1_lines.append(f"  4-phase patients: {n_pairs_4ph} (patient x band x pair)")
e1_lines.append(f"  All 6 patients: {n_total} (including partial-phase patients)")

# 4. Kendall's W
e1_lines.append("")
e1_lines.append("## 4. Kendall's W for H4")
e1_lines.append("  W = 0.043 (from analyze_h4_gradient.py)")
e1_lines.append("  Interpretation: < 0.3 = weak agreement, "
                "0.3-0.5 = moderate, > 0.7 = strong")

# 5. Number of nodes
e1_lines.append("")
e1_lines.append("## 5. Number of nodes per patient")
for pat in ALL_PATIENTS:
    for band in BANDS:
        key = (pat, band, "rest_pre")
        if key in linkages:
            n_nodes = linkages[key].shape[0] + 1
            e1_lines.append(f"  {pat} ({band}): {n_nodes} nodes")
            break

e1_path = OUTDIR / "report_E1_key_numbers.txt"
e1_path.write_text("\n".join(e1_lines))
print(f"Saved {e1_path.name}")

# Print key numbers
for line in e1_lines:
    print(f"  {line}")

print("\nDone!")
