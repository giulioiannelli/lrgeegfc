#!/usr/bin/env python3
"""BLOCK B3 + B4: Multiscale co-classification affinity matrices.

Computes and visualizes:
  B3a. Affinity matrices for Pat_02 (4 phases × 6 bands)
  B3b. Affinity difference matrices for Pat_02 (2 contrasts × 6 bands)
  B3c. Frobenius norm table for all patients
  B4.  Affinity vs raw FC comparison (ratio table)

Output: data/figures/metric_exploration/report_figures/
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import fcluster, leaves_list

from lrg_eegfc.config.paths import LRG_CACHE, MSC_CACHE, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.workflow.msc import load_msc_matrix

OUTDIR = FIGURES_ROOT / "metric_exploration" / "report_figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

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


# ═══════════════════════════════════════════════════════════════════════
# Co-classification affinity matrix
# ═══════════════════════════════════════════════════════════════════════
def compute_affinity(Z, K_max=60):
    """Co-classification affinity: A_ij = P(nodes i,j co-classified across k=2..K)."""
    n = Z.shape[0] + 1
    K = min(K_max, n - 1)
    A = np.zeros((n, n))
    for k in range(2, K + 1):
        labels = fcluster(Z, k, criterion="maxclust")
        same = (labels[:, None] == labels[None, :]).astype(np.float64)
        A += same
    A /= (K - 1)  # number of levels: k=2..K
    return A


# ── Load all linkage matrices ────────────────────────────────────────
print("Loading LRG linkage matrices...")
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


# ── Compute affinity matrices ────────────────────────────────────────
print("Computing affinity matrices...")
affinities = {}
for key, Z in linkages.items():
    affinities[key] = compute_affinity(Z, K_max=60)
    # Print progress every 20
    if len(affinities) % 20 == 0:
        print(f"  {len(affinities)}/{len(linkages)} done")
print(f"Computed {len(affinities)} affinity matrices")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE B3a: Affinity matrices for Pat_02 (4 phases × 6 bands)
# ═══════════════════════════════════════════════════════════════════════
print("\nFigure B3a: Pat_02 affinity matrices...")

pat = "Pat_02"
phases = PATIENT_PHASES[pat]

# Get node ordering from rest_pre (first band available) for consistent display
ref_key = (pat, BANDS[0], "rest_pre")
if ref_key in linkages:
    node_order = leaves_list(linkages[ref_key])
else:
    node_order = None

fig, axes = plt.subplots(len(phases), len(BANDS), figsize=(24, 16))

for ip, phase in enumerate(phases):
    for ib, band in enumerate(BANDS):
        ax = axes[ip, ib]
        key = (pat, band, phase)
        if key not in affinities:
            ax.set_visible(False)
            continue
        A = affinities[key]
        # Reorder by rest_pre leaf order for this specific band
        band_ref = (pat, band, "rest_pre")
        if band_ref in linkages:
            order = leaves_list(linkages[band_ref])
        elif node_order is not None:
            order = node_order
        else:
            order = np.arange(A.shape[0])
        A_ord = A[np.ix_(order, order)]

        im = ax.imshow(A_ord, cmap="magma", vmin=0, vmax=1, aspect="equal",
                        interpolation="nearest")
        ax.set_xticks([])
        ax.set_yticks([])

        if ip == 0:
            ax.set_title(BAND_TEX[band], fontsize=14, fontweight="bold")
        if ib == 0:
            ax.set_ylabel(PHASE_SHORT[phase], fontsize=14, fontweight="bold")

# Shared colorbar
cbar = fig.colorbar(im, ax=axes, shrink=0.4, pad=0.02, location="right")
cbar.set_label("Co-classification affinity $A_{ij}$", fontsize=12)

fig.suptitle(
    f"{pat} — Multiscale co-classification affinity matrices\n"
    f"$A_{{ij}}$ = P(nodes $i$,$j$ co-classified across $k=2..60$)  |  "
    f"Nodes ordered by rest_pre dendrogram leaf order",
    fontsize=13, y=1.01,
)
fig.savefig(OUTDIR / "fig_affinity_matrices_Pat02.pdf", bbox_inches="tight",
            dpi=150)
plt.close(fig)
print("Saved fig_affinity_matrices_Pat02.pdf")


# ═══════════════════════════════════════════════════════════════════════
# FIGURE B3b: Affinity difference matrices for Pat_02
# ═══════════════════════════════════════════════════════════════════════
print("\nFigure B3b: Pat_02 affinity differences...")

contrasts = [
    ("rest_post", "rest_pre", "Post $-$ Pre (task trace)"),
    ("task_learn", "rest_pre", "TL $-$ Pre (task onset)"),
]

fig, axes = plt.subplots(len(contrasts), len(BANDS), figsize=(24, 9))

vmax_diff = 0
# Pre-compute vmax
for phase_a, phase_b, _ in contrasts:
    for band in BANDS:
        ka = (pat, band, phase_a)
        kb = (pat, band, phase_b)
        if ka in affinities and kb in affinities:
            diff = affinities[ka] - affinities[kb]
            vmax_diff = max(vmax_diff, np.percentile(np.abs(diff), 99))

for ic, (phase_a, phase_b, title) in enumerate(contrasts):
    for ib, band in enumerate(BANDS):
        ax = axes[ic, ib]
        ka = (pat, band, phase_a)
        kb = (pat, band, phase_b)
        if ka not in affinities or kb not in affinities:
            ax.set_visible(False)
            continue

        diff = affinities[ka] - affinities[kb]
        # Reorder
        band_ref = (pat, band, "rest_pre")
        if band_ref in linkages:
            order = leaves_list(linkages[band_ref])
        else:
            order = np.arange(diff.shape[0])
        diff_ord = diff[np.ix_(order, order)]

        im = ax.imshow(diff_ord, cmap="RdBu_r", vmin=-vmax_diff, vmax=vmax_diff,
                        aspect="equal", interpolation="nearest")
        ax.set_xticks([])
        ax.set_yticks([])

        if ic == 0:
            ax.set_title(BAND_TEX[band], fontsize=14, fontweight="bold")
        if ib == 0:
            ax.set_ylabel(title, fontsize=11, fontweight="bold")

cbar = fig.colorbar(im, ax=axes, shrink=0.5, pad=0.02, location="right")
cbar.set_label("$\\Delta A_{ij}$ (affinity change)", fontsize=12)

fig.suptitle(
    f"{pat} — Affinity difference matrices\n"
    f"Red = nodes became more co-classified | Blue = nodes became less co-classified",
    fontsize=13, y=1.01,
)
fig.savefig(OUTDIR / "fig_affinity_difference_Pat02.pdf", bbox_inches="tight",
            dpi=150)
plt.close(fig)
print("Saved fig_affinity_difference_Pat02.pdf")


# ═══════════════════════════════════════════════════════════════════════
# B3c: Frobenius norm table for all patients
# ═══════════════════════════════════════════════════════════════════════
print("\nComputing Frobenius norm table...")

from itertools import combinations

frob_lines = ["# Frobenius norm of affinity difference ||Delta_A||_F",
              "# for all (patient, band, phase_pair) combinations",
              "",
              "patient,band,phase1,phase2,frob_norm,frob_norm_relative"]

for pat_ in ALL_PATIENTS:
    phases_ = PATIENT_PHASES[pat_]
    for band in BANDS:
        for p1, p2 in combinations(phases_, 2):
            k1 = (pat_, band, p1)
            k2 = (pat_, band, p2)
            if k1 not in affinities or k2 not in affinities:
                continue
            diff = affinities[k1] - affinities[k2]
            frob = np.linalg.norm(diff, "fro")
            frob_ref = np.linalg.norm(affinities[k1], "fro")
            frob_rel = frob / frob_ref if frob_ref > 0 else 0
            frob_lines.append(
                f"{pat_},{band},{p1},{p2},{frob:.4f},{frob_rel:.4f}")

frob_path = OUTDIR / "table_affinity_frobenius_norms.csv"
frob_path.write_text("\n".join(frob_lines))
print(f"Saved {frob_path.name}")


# ═══════════════════════════════════════════════════════════════════════
# B4: Affinity vs raw FC comparison
# ═══════════════════════════════════════════════════════════════════════
print("\nComputing affinity vs raw FC comparison (B4)...")

b4_lines = ["# Affinity vs raw FC: relative Frobenius distance comparison",
            "# d_aff = ||A_post - A_pre||_F / ||A_pre||_F",
            "# d_fc  = ||W_post - W_pre||_F / ||W_pre||_F",
            "",
            "patient,band,d_fc,d_aff,ratio"]

ratios = []
for pat_ in ALL_PATIENTS:
    if "rest_pre" not in PATIENT_PHASES[pat_] or "rest_post" not in PATIENT_PHASES[pat_]:
        continue
    for band in BANDS:
        # Affinity distance
        ka_pre = (pat_, band, "rest_pre")
        ka_post = (pat_, band, "rest_post")
        if ka_pre not in affinities or ka_post not in affinities:
            continue
        d_aff = np.linalg.norm(affinities[ka_post] - affinities[ka_pre], "fro")
        ref_aff = np.linalg.norm(affinities[ka_pre], "fro")
        d_aff_rel = d_aff / ref_aff if ref_aff > 0 else 0

        # FC distance — load MSC matrices
        try:
            W_pre = load_msc_matrix(pat_, "rest_pre", band, cache_root=MSC_CACHE,
                                     sparsify="soft")
            W_post = load_msc_matrix(pat_, "rest_post", band, cache_root=MSC_CACHE,
                                      sparsify="soft")
            d_fc = np.linalg.norm(W_post - W_pre, "fro")
            ref_fc = np.linalg.norm(W_pre, "fro")
            d_fc_rel = d_fc / ref_fc if ref_fc > 0 else 0
        except Exception:
            continue

        ratio = d_aff_rel / d_fc_rel if d_fc_rel > 0 else float("nan")
        ratios.append(ratio)
        b4_lines.append(f"{pat_},{band},{d_fc_rel:.4f},{d_aff_rel:.4f},{ratio:.2f}")

b4_path = OUTDIR / "table_affinity_vs_fc.csv"
b4_path.write_text("\n".join(b4_lines))
print(f"Saved {b4_path.name}")

if ratios:
    mean_ratio = np.nanmean(ratios)
    std_ratio = np.nanstd(ratios)
    print(f"  Mean ratio d_aff/d_fc = {mean_ratio:.2f} +/- {std_ratio:.2f}")

print("\nDone!")
