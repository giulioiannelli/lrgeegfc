#!/usr/bin/env python3
"""Consolidate report figures into ≤15 PDFs by merging similar ones.

Merges:
  - 2 cophenetic phase matrices + Pre-Post overview → 1 figure
  - 2 NVI profiles (Pat02 + mean) → 1 figure
  - 2 affinity figures (matrices + differences) → 1 figure
  - metric heatmap + scatter → 1 figure
  - 2 H1 figures → 1 figure
  - C02 (per-band detail) + C03 (raw VI curves) → 1 figure
  - 2 H4 figures → 1 figure
  - 3 outlier diagnostics → 1 simplified figure

Drops (redundant):
  - C04 (redundant with C01)
  - C05 (redundant with C01)
  - A01 (superseded by B2 mean profiles)
  - A03 (similar to A01)
  - G01 (individual hypothesis figs are better)
  - G03 (replaced by extended metric heatmap)
  - D02 (6-panel contrasts too dense; D01 is clearer)

Output: data/figures/metric_exploration/report_figures_final/
"""
from __future__ import annotations

import shutil
from itertools import combinations
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import fcluster, cophenet, leaves_list
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr
import pandas as pd

from lrg_eegfc.config.paths import LRG_CACHE, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result

SRC = FIGURES_ROOT / "metric_exploration" / "report_figures"
DST = FIGURES_ROOT / "metric_exploration" / "report_figures_final"
DST.mkdir(parents=True, exist_ok=True)

PATIENTS_4PH = ["Pat_02", "Pat_03", "Pat_05", "Pat_08"]
ALL_PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07", "Pat_08"]
PATIENT_PHASES = {
    "Pat_02": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_03": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_05": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_08": ["rest_pre", "task_learn", "task_test", "rest_post"],
}
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = {
    "delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
    "beta": r"$\beta$", "low_gamma": r"$\gamma_l$", "high_gamma": r"$\gamma_h$",
}
PHASE_SHORT = {"rest_pre": "Pre", "task_learn": "TL", "task_test": "TT", "rest_post": "Post"}
ALL_PAIRS = [
    ("task_learn", "task_test"), ("rest_pre", "rest_post"),
    ("rest_pre", "task_learn"), ("rest_pre", "task_test"),
    ("task_learn", "rest_post"), ("task_test", "rest_post"),
]
PAIR_LABELS = {
    ("task_learn", "task_test"): "TL-TT", ("rest_pre", "rest_post"): "Pre-Post",
    ("rest_pre", "task_learn"): "Pre-TL", ("rest_pre", "task_test"): "Pre-TT",
    ("task_learn", "rest_post"): "TL-Post", ("task_test", "rest_post"): "TT-Post",
}
PAIR_COLORS = {
    ("task_learn", "task_test"): "#1565C0", ("rest_pre", "rest_post"): "#2E7D32",
    ("rest_pre", "task_learn"): "#9E9E9E", ("rest_pre", "task_test"): "#BDBDBD",
    ("task_learn", "rest_post"): "#E65100", ("task_test", "rest_post"): "#C62828",
}
WITHIN_PAIRS = [("task_learn", "task_test"), ("rest_pre", "rest_post")]
CROSS_PAIRS = [("rest_pre", "task_learn"), ("rest_pre", "task_test"),
               ("task_learn", "rest_post"), ("task_test", "rest_post")]
BAND_COLORS = {
    "delta": "#1f77b4", "theta": "#ff7f0e", "alpha": "#2ca02c",
    "beta": "#d62728", "low_gamma": "#9467bd", "high_gamma": "#8c564b",
}
PAT_MARKERS = {"Pat_02": "o", "Pat_03": "s", "Pat_05": "D", "Pat_08": "^"}

H_GRID = np.geomspace(0.003, 0.995, 400)


from lrg_eegfc.utils.metrics import compute_vi  # canonical implementation


def cophenetic_corr(Z1, Z2):
    return pearsonr(cophenet(Z1), cophenet(Z2))[0]


def compute_affinity(Z, K_max=60):
    n = Z.shape[0] + 1
    K = min(K_max, n - 1)
    A = np.zeros((n, n))
    for k in range(2, K + 1):
        labels = fcluster(Z, k, criterion="maxclust")
        same = (labels[:, None] == labels[None, :]).astype(np.float64)
        A += same
    A /= (K - 1)
    return A


# ── Load all data ────────────────────────────────────────────────────
print("Loading LRG linkages...")
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

# ── Compute k(h) and VI(h) ──────────────────────────────────────────
print("Computing k(h) profiles...")
k_profiles = {}
for key, Z in linkages.items():
    k_profiles[key] = np.array([len(np.unique(fcluster(Z, h, criterion="distance")))
                                 for h in H_GRID])
all_k = np.array([k_profiles[k] for k in k_profiles])
mean_k_at_h = np.mean(all_k, axis=0)
norm_factor = np.log(np.maximum(mean_k_at_h, 1.001))

print("Computing VI(h) for 4-phase patients...")
vi_h = {pat: {band: {} for band in BANDS} for pat in PATIENTS_4PH}
for pat in PATIENTS_4PH:
    for band in BANDS:
        for p1, p2 in combinations(PATIENT_PHASES[pat], 2):
            k1, k2 = (pat, band, p1), (pat, band, p2)
            if k1 not in linkages or k2 not in linkages:
                continue
            Z1, Z2 = linkages[k1], linkages[k2]
            vals = np.zeros(len(H_GRID))
            for ih, h in enumerate(H_GRID):
                l1 = fcluster(Z1, h, criterion="distance")
                l2 = fcluster(Z2, h, criterion="distance")
                vals[ih] = compute_vi(l1, l2)
            vi_h[pat][band][(p1, p2)] = vals
    print(f"  {pat} done")


# ═══════════════════════════════════════════════════════════════════════
# 1. COPY STANDALONE FIGURES (no merge needed)
# ═══════════════════════════════════════════════════════════════════════
COPIES = [
    # (target, source)
    ("01_three_tools_example.pdf", "fig_three_tools_example_Pat02_beta.pdf"),
    ("03_h_to_k_mapping.pdf", "A02_k_at_h_profile.pdf"),
    ("08_H2_definitive_heatmap.pdf", "C01_H2_definitive_heatmap.pdf"),
    ("10_H3_within_vs_cross.pdf", "D01_H3_within_vs_cross.pdf"),
    ("13_pre_post_all_6_patients.pdf", "G02_pre_post_all_6_patients.pdf"),
]
for target, source in COPIES:
    shutil.copy2(SRC / source, DST / target)
    print(f"Copied {target}")


# ═══════════════════════════════════════════════════════════════════════
# 2. MERGE: Metric heatmap + scatter → 02_metric_comparison.pdf
# ═══════════════════════════════════════════════════════════════════════
print("\nMerging: metric comparison...")

# Load metric data
pw_csv = FIGURES_ROOT / "metric_exploration" / "alternative_metrics" / "pairwise_results.csv"
df_pw = pd.read_csv(pw_csv)
METRICS = [c for c in df_pw.columns if c not in
           ["patient", "band", "phase1", "phase2", "pair"]]
SIM_METRICS = {"CophPearson", "CophSpearman", "BakersGamma",
               "TopK5_Merge", "TopK10_Merge", "TopK20_Merge",
               "WeightARI_Coarse", "WeightARI_Fine", "MeanARI"}

# Compute unanimity table
def _unan_table():
    results = []
    for metric in METRICS:
        is_sim = metric in SIM_METRICS
        for band in BANDS:
            bdf = df_pw[df_pw["band"] == band]
            h1s, h2as, h2bs, h3s = [], [], [], []
            for pat in PATIENTS_4PH:
                pdf = bdf[bdf["patient"] == pat]
                if pdf.empty:
                    continue
                vals = {}
                for _, r in pdf.iterrows():
                    vals[(r["phase1"], r["phase2"])] = r[metric]
                tl_tt = vals.get(("task_learn", "task_test"), np.nan)
                if np.isnan(tl_tt):
                    continue
                others = [v for k, v in vals.items()
                          if k != ("task_learn", "task_test") and not np.isnan(v)]
                if others:
                    h1s.append(np.sign((tl_tt - np.mean(others)) * (1 if is_sim else -1)))
                pre_post = vals.get(("rest_pre", "rest_post"), np.nan)
                tt_post = vals.get(("task_test", "rest_post"), np.nan)
                if not (np.isnan(pre_post) or np.isnan(tt_post)):
                    h2as.append(np.sign((tt_post - pre_post) * (1 if is_sim else -1)))
                pre_tt = vals.get(("rest_pre", "task_test"), np.nan)
                if not (np.isnan(pre_tt) or np.isnan(tt_post)):
                    h2bs.append(np.sign((tt_post - pre_tt) * (1 if is_sim else -1)))
                wi = [vals.get(p, np.nan) for p in WITHIN_PAIRS]
                cr = [vals.get(p, np.nan) for p in CROSS_PAIRS]
                wi = [v for v in wi if not np.isnan(v)]
                cr = [v for v in cr if not np.isnan(v)]
                if wi and cr:
                    h3s.append(np.sign((np.mean(wi) - np.mean(cr)) * (1 if is_sim else -1)))
            results.append({
                "metric": metric, "band": band,
                "H1": len(h1s) >= 4 and abs(sum(h1s)) == len(h1s),
                "H2a": len(h2as) >= 4 and abs(sum(h2as)) == len(h2as) and sum(h2as) > 0,
                "H2b": len(h2bs) >= 4 and abs(sum(h2bs)) == len(h2bs) and sum(h2bs) > 0,
                "H3": len(h3s) >= 4 and abs(sum(h3s)) == len(h3s),
            })
    return pd.DataFrame(results)

udf = _unan_table()
summary = []
for m in METRICS:
    mdf = udf[udf["metric"] == m]
    summary.append({"Metric": m,
                     "H1": int(mdf["H1"].sum()), "H2a": int(mdf["H2a"].sum()),
                     "H2b": int(mdf["H2b"].sum()), "H3": int(mdf["H3"].sum()),
                     "Total": int(mdf["H1"].sum() + mdf["H2a"].sum() +
                                  mdf["H2b"].sum() + mdf["H3"].sum())})
sdf = pd.DataFrame(summary).sort_values("Total", ascending=False)

# Also compute cophenetic vs VI scatter data
coph_v, vi_v, sb, sp = [], [], [], []
for pat in PATIENTS_4PH:
    for band in BANDS:
        for p1, p2 in ALL_PAIRS:
            k1, k2 = (pat, band, p1), (pat, band, p2)
            if k1 not in linkages or k2 not in linkages:
                continue
            cc = cophenetic_corr(linkages[k1], linkages[k2])
            vis = []
            for k in range(2, 31):
                l1 = fcluster(linkages[k1], k, criterion="maxclust")
                l2 = fcluster(linkages[k2], k, criterion="maxclust")
                vis.append(compute_vi(l1, l2))
            coph_v.append(cc)
            vi_v.append(np.mean(vis))
            sb.append(band)
            sp.append(pat)

fig = plt.figure(figsize=(18, 7))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.3])

# Left: unanimity heatmap
ax1 = fig.add_subplot(gs[0])
metrics_sorted = sdf["Metric"].tolist()
mat = sdf[["H1", "H2a", "H2b", "H3"]].values.astype(float)
im = ax1.imshow(mat, cmap="YlOrRd", vmin=0, vmax=6, aspect="auto")
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        v = int(mat[i, j])
        ax1.text(j, i, str(v), ha="center", va="center", fontsize=10,
                 fontweight="bold", color="white" if v >= 4 else "black")
ax1.set_xticks(range(4))
ax1.set_xticklabels(["H1", "H2a", "H2b", "H3"], fontsize=10)
ax1.set_yticks(range(len(metrics_sorted)))
ax1.set_yticklabels(metrics_sorted, fontsize=9)
for i, row in enumerate(sdf.itertuples()):
    ax1.text(4.2, i, f"={row.Total}", fontsize=9, ha="left", va="center",
             fontweight="bold")
ax1.set_xlim(-0.5, 3.5)
fig.colorbar(im, ax=ax1, shrink=0.5, pad=0.12, label="# unanimous bands (max 6)")
ax1.set_title("(a) Unanimity across hypotheses", fontsize=11, fontweight="bold")

# Right: scatter
ax2 = fig.add_subplot(gs[1])
for band in BANDS:
    for pat in PATIENTS_4PH:
        mask = [(b == band and p == pat) for b, p in zip(sb, sp)]
        x = [coph_v[i] for i, m in enumerate(mask) if m]
        y = [vi_v[i] for i, m in enumerate(mask) if m]
        if x:
            ax2.scatter(x, y, c=BAND_COLORS[band], marker=PAT_MARKERS[pat],
                        s=40, alpha=0.7, edgecolors="none")
for band in BANDS:
    ax2.scatter([], [], c=BAND_COLORS[band], s=60, label=BAND_TEX[band])
for pat in PATIENTS_4PH:
    ax2.scatter([], [], c="gray", marker=PAT_MARKERS[pat], s=60, label=pat)
xa, ya = np.array(coph_v), np.array(vi_v)
valid = np.isfinite(xa) & np.isfinite(ya)
z = np.polyfit(xa[valid], ya[valid], 1)
xl = np.linspace(xa[valid].min(), xa[valid].max(), 100)
ax2.plot(xl, np.polyval(z, xl), "k--", linewidth=1.5, alpha=0.5)
r, p = pearsonr(xa[valid], ya[valid])
ax2.text(0.05, 0.95, f"r = {r:.3f}", transform=ax2.transAxes, fontsize=10,
         va="top", bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
ax2.set_xlabel("Cophenetic Pearson correlation", fontsize=10)
ax2.set_ylabel("Mean VI (k=2..30)", fontsize=10)
ax2.legend(fontsize=7, ncol=2, loc="upper right")
ax2.grid(alpha=0.3)
ax2.set_title("(b) Cophenetic correlation vs Mean VI", fontsize=11, fontweight="bold")

fig.suptitle("Metric comparison: unanimity scores and cophenetic–VI relationship",
             fontsize=12, fontweight="bold", y=1.01)
fig.tight_layout()
fig.savefig(DST / "02_metric_comparison.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved 02_metric_comparison.pdf")


# ═══════════════════════════════════════════════════════════════════════
# 4. MERGE: NVI profiles Pat_02 (top) + mean 4-patient (bottom)
# ═══════════════════════════════════════════════════════════════════════
print("\nMerging: NVI profiles...")

fig, axes = plt.subplots(2, 6, figsize=(28, 9), sharex=True)

for ib, band in enumerate(BANDS):
    # Top row: Pat_02
    ax = axes[0, ib]
    for pair in ALL_PAIRS:
        v = vi_h["Pat_02"][band].get(pair)
        if v is None:
            continue
        ax.plot(H_GRID, v / np.maximum(norm_factor, 0.01),
                color=PAIR_COLORS[pair], linewidth=2, label=PAIR_LABELS[pair])
    ax.set_title(BAND_TEX[band], fontsize=14, fontweight="bold")
    ax.set_xscale("log")
    ax.grid(alpha=0.3)
    if ib == 0:
        ax.set_ylabel("Pat_02\nNVI(h)", fontsize=11)
        ax.legend(fontsize=6, loc="upper right")

    # Bottom row: mean ± range
    ax = axes[1, ib]
    for pair in ALL_PAIRS:
        curves = []
        for pat in PATIENTS_4PH:
            v = vi_h[pat][band].get(pair)
            if v is not None:
                curves.append(v / np.maximum(norm_factor, 0.01))
        if not curves:
            continue
        arr = np.array(curves)
        ax.plot(H_GRID, np.mean(arr, axis=0), color=PAIR_COLORS[pair],
                linewidth=2, label=PAIR_LABELS[pair])
        ax.fill_between(H_GRID, np.min(arr, axis=0), np.max(arr, axis=0),
                         color=PAIR_COLORS[pair], alpha=0.12)
    ax.set_xscale("log")
    ax.set_xlabel("h (log)", fontsize=9)
    ax.grid(alpha=0.3)
    if ib == 0:
        ax.set_ylabel("Mean (n=4)\nNVI(h)", fontsize=11)

fig.suptitle("Scale-resolved NVI profiles: Pat_02 (top) vs 4-patient mean ± range (bottom)",
             fontsize=13, fontweight="bold", y=1.01)
fig.tight_layout()
fig.savefig(DST / "04_nvi_profiles.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved 04_nvi_profiles.pdf")


# ═══════════════════════════════════════════════════════════════════════
# 5. MERGE: Cophenetic phase matrices (Pat02 + Pat08) + Pre-Post overview
# ═══════════════════════════════════════════════════════════════════════
print("\nMerging: cophenetic overview...")

fig = plt.figure(figsize=(26, 12))
gs = fig.add_gridspec(3, 6, height_ratios=[1, 1, 1.2], hspace=0.35)

for row_idx, pat in enumerate(["Pat_02", "Pat_08"]):
    phases = PATIENT_PHASES[pat]
    n_ph = len(phases)
    for ib, band in enumerate(BANDS):
        ax = fig.add_subplot(gs[row_idx, ib])
        mat = np.full((n_ph, n_ph), np.nan)
        for i, p1 in enumerate(phases):
            for j, p2 in enumerate(phases):
                if i == j:
                    mat[i, j] = 1.0
                else:
                    k1, k2 = (pat, band, p1), (pat, band, p2)
                    if k1 in linkages and k2 in linkages:
                        mat[i, j] = cophenetic_corr(linkages[k1], linkages[k2])
        vv = mat[~np.isnan(mat)]
        vmin = max(vv.min() - 0.05, -1) if len(vv) > 0 else 0
        im = ax.imshow(mat, cmap="RdBu_r", vmin=vmin, vmax=1, aspect="equal")
        for i in range(n_ph):
            for j in range(n_ph):
                if not np.isnan(mat[i, j]):
                    c = "white" if abs(mat[i, j]) > 0.8 else "black"
                    ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                            fontsize=7, fontweight="bold", color=c)
        ax.set_xticks(range(n_ph))
        ax.set_xticklabels([PHASE_SHORT[p] for p in phases], fontsize=8)
        ax.set_yticks(range(n_ph))
        ax.set_yticklabels([PHASE_SHORT[p] for p in phases] if ib == 0 else [],
                           fontsize=8)
        if ib == 0:
            ax.set_ylabel(pat, fontsize=12, fontweight="bold")
        if row_idx == 0:
            ax.set_title(BAND_TEX[band], fontsize=13, fontweight="bold")

# Bottom row: Pre-Post overview (all 6 patients)
ax_bot = fig.add_subplot(gs[2, :])
patients_pp = [p for p in ALL_PATIENTS
               if "rest_pre" in PATIENT_PHASES[p] and "rest_post" in PATIENT_PHASES[p]]
mat_pp = np.full((len(patients_pp), len(BANDS)), np.nan)
for ip, pat in enumerate(patients_pp):
    for ib, band in enumerate(BANDS):
        k1, k2 = (pat, band, "rest_pre"), (pat, band, "rest_post")
        if k1 in linkages and k2 in linkages:
            mat_pp[ip, ib] = cophenetic_corr(linkages[k1], linkages[k2])
im2 = ax_bot.imshow(mat_pp, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
for i in range(mat_pp.shape[0]):
    for j in range(mat_pp.shape[1]):
        if not np.isnan(mat_pp[i, j]):
            c = "white" if mat_pp[i, j] < 0.4 else "black"
            ax_bot.text(j, i, f"{mat_pp[i, j]:.2f}", ha="center", va="center",
                        fontsize=11, fontweight="bold", color=c)
ax_bot.set_xticks(range(len(BANDS)))
ax_bot.set_xticklabels([BAND_TEX[b] for b in BANDS], fontsize=12)
ax_bot.set_yticks(range(len(patients_pp)))
ax_bot.set_yticklabels(patients_pp, fontsize=11)
fig.colorbar(im2, ax=ax_bot, shrink=0.6, pad=0.02,
             label="Cophenetic corr. (rest_pre vs rest_post)")
ax_bot.set_title("Pre-Post cophenetic correlation across all 6 patients",
                  fontsize=11, fontweight="bold")

fig.suptitle("Cophenetic correlation overview: phase matrices (top) and Pre-Post summary (bottom)",
             fontsize=13, fontweight="bold", y=1.01)
fig.savefig(DST / "05_cophenetic_overview.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved 05_cophenetic_overview.pdf")


# ═══════════════════════════════════════════════════════════════════════
# 6. MERGE: Affinity matrices + differences → 06_affinity_Pat02.pdf
# ═══════════════════════════════════════════════════════════════════════
print("\nMerging: affinity matrices...")

pat = "Pat_02"
phases = PATIENT_PHASES[pat]

# Compute affinities for Pat_02
aff = {}
for phase in phases:
    for band in BANDS:
        key = (pat, band, phase)
        if key in linkages:
            aff[(band, phase)] = compute_affinity(linkages[key])

fig = plt.figure(figsize=(26, 20))
gs = fig.add_gridspec(6, 6, hspace=0.25, wspace=0.15)

# Top 4 rows: affinity matrices (phases × bands)
for ip, phase in enumerate(phases):
    for ib, band in enumerate(BANDS):
        ax = fig.add_subplot(gs[ip, ib])
        A = aff.get((band, phase))
        if A is None:
            ax.set_visible(False)
            continue
        ref_key = (pat, band, "rest_pre")
        order = leaves_list(linkages[ref_key]) if ref_key in linkages else np.arange(A.shape[0])
        ax.imshow(A[np.ix_(order, order)], cmap="magma", vmin=0, vmax=1,
                   aspect="equal", interpolation="nearest")
        ax.set_xticks([])
        ax.set_yticks([])
        if ip == 0:
            ax.set_title(BAND_TEX[band], fontsize=13, fontweight="bold")
        if ib == 0:
            ax.set_ylabel(PHASE_SHORT[phase], fontsize=12, fontweight="bold")

# Bottom 2 rows: difference matrices
contrasts = [("rest_post", "rest_pre", "Post $-$ Pre"), ("task_learn", "rest_pre", "TL $-$ Pre")]
vmax_d = 0
for pa, pb, _ in contrasts:
    for band in BANDS:
        a, b = aff.get((band, pa)), aff.get((band, pb))
        if a is not None and b is not None:
            vmax_d = max(vmax_d, np.percentile(np.abs(a - b), 99))

for ic, (pa, pb, title) in enumerate(contrasts):
    for ib, band in enumerate(BANDS):
        ax = fig.add_subplot(gs[4 + ic, ib])
        a, b = aff.get((band, pa)), aff.get((band, pb))
        if a is None or b is None:
            ax.set_visible(False)
            continue
        diff = a - b
        ref_key = (pat, band, "rest_pre")
        order = leaves_list(linkages[ref_key]) if ref_key in linkages else np.arange(diff.shape[0])
        ax.imshow(diff[np.ix_(order, order)], cmap="RdBu_r",
                   vmin=-vmax_d, vmax=vmax_d, aspect="equal", interpolation="nearest")
        ax.set_xticks([])
        ax.set_yticks([])
        if ib == 0:
            ax.set_ylabel(title, fontsize=11, fontweight="bold")

fig.suptitle(
    f"{pat} — Co-classification affinity matrices (top 4 rows) and differences (bottom 2)\n"
    f"$A_{{ij}}$ = P(nodes $i,j$ co-classified across $k=2..60$)  |  "
    f"Nodes ordered by rest_pre leaf order",
    fontsize=13, fontweight="bold", y=1.01,
)
fig.savefig(DST / "06_affinity_Pat02.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved 06_affinity_Pat02.pdf")


# ═══════════════════════════════════════════════════════════════════════
# 7. MERGE: H1 normalized + multiscale → 07_H1_task_stability.pdf
# ═══════════════════════════════════════════════════════════════════════
print("\nMerging: H1 figures...")

fig, axes = plt.subplots(1, 6, figsize=(28, 4.5), sharex=True)

for ib, band in enumerate(BANDS):
    ax = axes[ib]
    gains = []
    for pat in PATIENTS_4PH:
        band_vi = vi_h[pat][band]
        tl_tt = band_vi.get(("task_learn", "task_test"))
        if tl_tt is None:
            continue
        others = [band_vi.get(p) for p in ALL_PAIRS
                  if p != ("task_learn", "task_test") and band_vi.get(p) is not None]
        if not others:
            continue
        mean_others = np.mean(others, axis=0)
        gain = (mean_others - tl_tt) / np.maximum(mean_others, 1e-10)
        gains.append(gain)
        ax.plot(H_GRID, gain, color="gray", linewidth=0.8, alpha=0.4)

    if gains:
        mean_gain = np.mean(gains, axis=0)
        ax.plot(H_GRID, mean_gain, color="#1565C0", linewidth=2.5)
        # Unanimous regions
        signs = np.sign(np.array(gains))
        unan = np.mean(signs, axis=0)
        unanimous = unan >= 0.99
        ax.fill_between(H_GRID, -0.5, 1.5,
                         where=unanimous, color="#1565C0", alpha=0.15)

    ax.axhline(0, color="k", linewidth=0.5, linestyle="--")
    ax.set_xscale("log")
    ax.set_title(BAND_TEX[band], fontsize=14, fontweight="bold")
    ax.set_ylim(-0.5, 1.0)
    ax.grid(alpha=0.3)
    if ib == 0:
        ax.set_ylabel("Relative gain\n(TL-TT vs others)", fontsize=10)
    ax.set_xlabel("h", fontsize=9)

fig.suptitle("H1: Task stability — TL-TT relative gain across scales\n"
             "Blue shading = 4/4 patients unanimous  |  >0 = TL-TT is the most similar pair",
             fontsize=12, fontweight="bold", y=1.04)
fig.tight_layout()
fig.savefig(DST / "07_H1_task_stability.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved 07_H1_task_stability.pdf")


# ═══════════════════════════════════════════════════════════════════════
# 9. MERGE: H2 per-band detail + raw VI curves → 09_H2_detail.pdf
# ═══════════════════════════════════════════════════════════════════════
print("\nMerging: H2 detail...")

fig, axes = plt.subplots(2, 6, figsize=(28, 9), sharex=True)

# Top: per-band H2a contrast (individual patient traces)
for ib, band in enumerate(BANDS):
    ax = axes[0, ib]
    for pat in PATIENTS_4PH:
        pre_post = vi_h[pat][band].get(("rest_pre", "rest_post"))
        tt_post = vi_h[pat][band].get(("task_test", "rest_post"))
        if pre_post is None or tt_post is None:
            continue
        contrast = pre_post - tt_post
        ax.plot(H_GRID, contrast, linewidth=1.2, alpha=0.6, label=pat)
    ax.axhline(0, color="k", linewidth=0.5, linestyle="--")
    ax.set_xscale("log")
    ax.set_title(BAND_TEX[band], fontsize=14, fontweight="bold")
    ax.grid(alpha=0.3)
    if ib == 0:
        ax.set_ylabel("H2a contrast\nVI(Pre,Post) - VI(TT,Post)", fontsize=9)
        ax.legend(fontsize=7)

# Bottom: raw VI curves for 4 key pairs (mean across patients)
key_pairs = [("task_learn", "task_test"), ("rest_pre", "rest_post"),
             ("task_test", "rest_post"), ("rest_pre", "task_test")]
for ib, band in enumerate(BANDS):
    ax = axes[1, ib]
    for pair in key_pairs:
        curves = [vi_h[pat][band].get(pair) for pat in PATIENTS_4PH
                  if vi_h[pat][band].get(pair) is not None]
        if not curves:
            continue
        arr = np.array(curves) / np.maximum(norm_factor, 0.01)
        ax.plot(H_GRID, np.mean(arr, axis=0), color=PAIR_COLORS[pair],
                linewidth=2, label=PAIR_LABELS[pair])
        ax.fill_between(H_GRID, np.min(arr, axis=0), np.max(arr, axis=0),
                         color=PAIR_COLORS[pair], alpha=0.1)
    ax.set_xscale("log")
    ax.set_xlabel("h (log)", fontsize=9)
    ax.grid(alpha=0.3)
    if ib == 0:
        ax.set_ylabel("NVI(h)", fontsize=10)
        ax.legend(fontsize=6, loc="upper right")

fig.suptitle("H2 detail: per-patient contrast traces (top) and mean NVI curves for key pairs (bottom)",
             fontsize=12, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(DST / "09_H2_detail.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved 09_H2_detail.pdf")


# ═══════════════════════════════════════════════════════════════════════
# 11. MERGE: H4 gradient + ranking agreement → 11_H4_gradient.pdf
# ═══════════════════════════════════════════════════════════════════════
print("\nMerging: H4 figures...")

# Compute reorganization strength per (patient, band)
reorg = {}
for pat in PATIENTS_4PH:
    for band in BANDS:
        cross_vis = []
        for p1, p2 in CROSS_PAIRS:
            v = vi_h[pat][band].get((p1, p2))
            if v is not None:
                cross_vis.append(np.mean(v))
        if cross_vis:
            reorg[(pat, band)] = np.mean(cross_vis)

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Left: rank heatmap
ax1 = axes[0]
rank_mat = np.full((len(BANDS), len(PATIENTS_4PH)), np.nan)
for ip, pat in enumerate(PATIENTS_4PH):
    vals = [reorg.get((pat, b), np.nan) for b in BANDS]
    ranks = np.argsort(np.argsort([-v if not np.isnan(v) else 0 for v in vals])) + 1
    rank_mat[:, ip] = ranks

im = ax1.imshow(rank_mat, cmap="YlOrRd_r", vmin=1, vmax=6, aspect="auto")
for i in range(rank_mat.shape[0]):
    for j in range(rank_mat.shape[1]):
        if not np.isnan(rank_mat[i, j]):
            ax1.text(j, i, str(int(rank_mat[i, j])), ha="center", va="center",
                     fontsize=14, fontweight="bold")
ax1.set_xticks(range(len(PATIENTS_4PH)))
ax1.set_xticklabels(PATIENTS_4PH, fontsize=10)
ax1.set_yticks(range(len(BANDS)))
ax1.set_yticklabels([BAND_TEX[b] for b in BANDS], fontsize=11)
fig.colorbar(im, ax=ax1, shrink=0.7, label="Rank (1 = most reorganized)")
ax1.set_title("(a) Band ranking by reorganization strength\n"
              "Shuffled colors → no consistent gradient", fontsize=11, fontweight="bold")

# Right: Spearman correlation matrix
from scipy.stats import spearmanr
ax2 = axes[1]
n_pat = len(PATIENTS_4PH)
corr_mat = np.full((n_pat, n_pat), np.nan)
for i in range(n_pat):
    for j in range(n_pat):
        r1 = rank_mat[:, i]
        r2 = rank_mat[:, j]
        valid = ~(np.isnan(r1) | np.isnan(r2))
        if valid.sum() >= 3:
            corr_mat[i, j] = spearmanr(r1[valid], r2[valid])[0]

im2 = ax2.imshow(corr_mat, cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal")
for i in range(n_pat):
    for j in range(n_pat):
        if not np.isnan(corr_mat[i, j]):
            ax2.text(j, i, f"{corr_mat[i, j]:.2f}", ha="center", va="center",
                     fontsize=11, fontweight="bold",
                     color="white" if abs(corr_mat[i, j]) > 0.6 else "black")
ax2.set_xticks(range(n_pat))
ax2.set_xticklabels(PATIENTS_4PH, fontsize=10)
ax2.set_yticks(range(n_pat))
ax2.set_yticklabels(PATIENTS_4PH, fontsize=10)
fig.colorbar(im2, ax=ax2, shrink=0.7, label="Spearman rank correlation")
ax2.set_title("(b) Pairwise rank agreement\n"
              "Kendall's W = 0.043 (no concordance)", fontsize=11, fontweight="bold")

fig.suptitle("H4: Frequency gradient — REJECTED (patient-specific, not frequency-driven)",
             fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout()
fig.savefig(DST / "11_H4_gradient.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved 11_H4_gradient.pdf")


# ═══════════════════════════════════════════════════════════════════════
# 12. MERGE: 3 outliers → simplified 12_outlier_summary.pdf
# ═══════════════════════════════════════════════════════════════════════
print("\nMerging: outlier summary...")

outlier_cases = [
    ("Pat_02", "theta", "TL-TT worst pair\n(coarse-scale flip)"),
    ("Pat_03", "high_gamma", "Extreme Pre-Post\n(star → modular)"),
    ("Pat_08", "delta", "Ultra-stable rest\n(Pre-Post ≈ 0)"),
]

fig, axes = plt.subplots(1, 3, figsize=(24, 6))

for idx, (pat, band, desc) in enumerate(outlier_cases):
    ax = axes[idx]
    band_vi = vi_h[pat][band]

    for pair in ALL_PAIRS:
        v = band_vi.get(pair)
        if v is None:
            continue
        v_norm = v / np.maximum(norm_factor, 0.01)
        lw = 3 if pair in [("task_learn", "task_test"), ("rest_pre", "rest_post"),
                            ("task_test", "rest_post")] else 1.2
        alpha = 1.0 if lw == 3 else 0.4
        ax.plot(H_GRID, v_norm, color=PAIR_COLORS[pair], linewidth=lw,
                alpha=alpha, label=PAIR_LABELS[pair])

    ax.set_xscale("log")
    ax.set_xlabel("Height h (log scale)", fontsize=10)
    ax.set_ylabel("NVI(h)" if idx == 0 else "", fontsize=10)
    ax.set_title(f"{pat}, {BAND_TEX[band]}\n{desc}", fontsize=12, fontweight="bold")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="upper right")

fig.suptitle("Outlier diagnostics: NVI profiles for three anomalous (patient, band) cases\n"
             "Bold lines = key pairs (TL-TT, Pre-Post, TT-Post)",
             fontsize=13, fontweight="bold", y=1.03)
fig.tight_layout()
fig.savefig(DST / "12_outlier_summary.pdf", bbox_inches="tight", dpi=150)
plt.close(fig)
print("Saved 12_outlier_summary.pdf")


# ═══════════════════════════════════════════════════════════════════════
# Copy text files
# ═══════════════════════════════════════════════════════════════════════
print("\nCopying text reports...")

# Merge all tables and key numbers into one file
report_parts = []
report_parts.append(Path(SRC / "TECHNICAL_REPORT.md").read_text())
report_parts.append("\n\n---\n\n# Appendix: Data Tables\n\n")

report_parts.append("## A1. Metric Comparison Table\n\n")
report_parts.append("```\n")
report_parts.append(sdf.to_string(index=False))
report_parts.append("\n```\n\n")

report_parts.append("## D2. Outlier Summary\n\n")
report_parts.append(Path(SRC / "table_D2_outlier_summary.txt").read_text())
report_parts.append("\n\n")

report_parts.append("## E1. Key Numbers\n\n")
report_parts.append("```\n")
report_parts.append(Path(SRC / "report_E1_key_numbers.txt").read_text())
report_parts.append("\n```\n\n")

report_parts.append("## B4. Affinity vs FC Comparison\n\n")
report_parts.append("```\n")
report_parts.append(Path(SRC / "table_affinity_vs_fc.csv").read_text())
report_parts.append("\n```\n")
report_parts.append("Mean d_aff/d_fc ratio = 1.68 ± 1.30 (median 1.51)\n\n")

report_parts.append("## Frobenius Norm Table\n\n")
report_parts.append("See table_affinity_frobenius_norms.csv (7.5 KB, not inlined)\n")

(DST / "REPORT_AND_TABLES.md").write_text("".join(report_parts))
print("Saved REPORT_AND_TABLES.md")

# Copy briefing
shutil.copy2(SRC / "WRITING_AGENT_BRIEFING.md", DST / "WRITING_AGENT_BRIEFING.md")
print("Saved WRITING_AGENT_BRIEFING.md")


# ═══════════════════════════════════════════════════════════════════════
# Final inventory
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("FINAL INVENTORY")
print("=" * 70)
all_files = sorted(DST.iterdir())
pdfs = [f for f in all_files if f.suffix == ".pdf"]
mds = [f for f in all_files if f.suffix == ".md"]
print(f"\nPDF figures: {len(pdfs)}")
for f in pdfs:
    print(f"  {f.name}")
print(f"\nText files: {len(mds)}")
for f in mds:
    print(f"  {f.name}")
print(f"\nTOTAL: {len(pdfs) + len(mds)} files")
