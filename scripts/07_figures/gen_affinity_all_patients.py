#!/usr/bin/env python3
"""Generate 06_affinity_PatXX.pdf for Pat_03, Pat_05, Pat_07, Pat_08.

Same format as 06_affinity_Pat02.pdf: affinity matrices (phases × bands)
on top, difference matrices (Post−Pre, TL−Pre) on bottom.
"""
from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from scipy.cluster.hierarchy import fcluster, leaves_list

from lrg_eegfc.config.paths import LRG_CACHE, FIGURES_ROOT
from lrg_eegfc.workflow.lrg import load_lrg_result

DST = FIGURES_ROOT / "metric_exploration" / "report_figures_final"
DST.mkdir(parents=True, exist_ok=True)

PATIENT_PHASES = {
    "Pat_02": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_03": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_05": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_07": ["rest_pre", "task_learn", "task_test", "rest_post"],
    "Pat_08": ["rest_pre", "task_learn", "task_test", "rest_post"],
}
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_TEX = {
    "delta": r"$\delta$", "theta": r"$\theta$", "alpha": r"$\alpha$",
    "beta": r"$\beta$", "low_gamma": r"$\gamma_l$", "high_gamma": r"$\gamma_h$",
}
PHASE_SHORT = {"rest_pre": "Pre", "task_learn": "TL", "task_test": "TT", "rest_post": "Post"}


def compute_affinity(Z, K_max=60):
    """Co-classification affinity: A_ij = P(nodes i,j co-classified across k=2..K)."""
    n = Z.shape[0] + 1
    K = min(K_max, n - 1)
    A = np.zeros((n, n))
    for k in range(2, K + 1):
        labels = fcluster(Z, k, criterion="maxclust")
        same = (labels[:, None] == labels[None, :]).astype(np.float64)
        A += same
    A /= (K - 1)
    return A


def generate_affinity_figure(pat, phases, linkages):
    """Generate and save one 06_affinity_PatXX.pdf figure."""
    # Compute affinities
    aff = {}
    for phase in phases:
        for band in BANDS:
            key = (pat, band, phase)
            if key in linkages:
                aff[(band, phase)] = compute_affinity(linkages[key])

    contrasts = [("rest_post", "rest_pre", "Post $-$ Pre"), ("task_learn", "rest_pre", "TL $-$ Pre")]
    n_rows = len(phases) + len(contrasts)
    n_cols = len(BANDS)

    # Compact figure: ~3.2 in per cell + small cbar column
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3.2 + 0.8, n_rows * 3.0))
    fig.subplots_adjust(hspace=0.12, wspace=0.08)

    # ── Top rows: affinity matrices (phases × bands) ──
    row_ims_aff = {}  # row_idx -> last imshow for colorbar
    for ip, phase in enumerate(phases):
        for ib, band in enumerate(BANDS):
            ax = axes[ip, ib]
            A = aff.get((band, phase))
            if A is None:
                ax.set_visible(False)
                continue
            ref_key = (pat, band, "rest_pre")
            order = leaves_list(linkages[ref_key]) if ref_key in linkages else np.arange(A.shape[0])
            im = ax.imshow(A[np.ix_(order, order)], cmap="magma", vmin=0, vmax=1,
                           aspect="equal", interpolation="nearest")
            ax.set_xticks([])
            ax.set_yticks([])
            if ip == 0:
                ax.set_title(BAND_TEX[band], fontsize=13, fontweight="bold")
            if ib == 0:
                ax.set_ylabel(PHASE_SHORT[phase], fontsize=12, fontweight="bold")
            row_ims_aff[ip] = im

    # Per-row colorbar on the last column for affinity rows
    for ip in sorted(row_ims_aff):
        ax_last = axes[ip, -1]
        div = make_axes_locatable(ax_last)
        cax = div.append_axes("right", size="5%", pad=0.05)
        cb = fig.colorbar(row_ims_aff[ip], cax=cax)
        if ip == 0:
            cb.set_label("$A_{ij}$", fontsize=10)
        else:
            cb.set_ticks([0, 0.5, 1])

    # ── Compute vmax for difference colorscale ──
    vmax_d = 0
    for pa, pb, _ in contrasts:
        for band in BANDS:
            a, b = aff.get((band, pa)), aff.get((band, pb))
            if a is not None and b is not None:
                vmax_d = max(vmax_d, np.percentile(np.abs(a - b), 99))

    # ── Bottom rows: difference matrices ──
    row_ims_diff = {}
    for ic, (pa, pb, title) in enumerate(contrasts):
        row_idx = len(phases) + ic
        for ib, band in enumerate(BANDS):
            ax = axes[row_idx, ib]
            a, b = aff.get((band, pa)), aff.get((band, pb))
            if a is None or b is None:
                ax.set_visible(False)
                continue
            diff = a - b
            ref_key = (pat, band, "rest_pre")
            order = leaves_list(linkages[ref_key]) if ref_key in linkages else np.arange(diff.shape[0])
            im = ax.imshow(diff[np.ix_(order, order)], cmap="RdBu_r",
                           vmin=-vmax_d, vmax=vmax_d, aspect="equal",
                           interpolation="nearest")
            ax.set_xticks([])
            ax.set_yticks([])
            if ib == 0:
                ax.set_ylabel(title, fontsize=11, fontweight="bold")
            row_ims_diff[row_idx] = im

    # Per-row colorbar for difference rows
    for row_idx in sorted(row_ims_diff):
        ax_last = axes[row_idx, -1]
        div = make_axes_locatable(ax_last)
        cax = div.append_axes("right", size="5%", pad=0.05)
        cb = fig.colorbar(row_ims_diff[row_idx], cax=cax)
        if row_idx == len(phases):
            cb.set_label("$\\Delta A_{ij}$", fontsize=10)

    outpath = DST / f"06_affinity_{pat.replace('_', '')}.pdf"
    fig.savefig(outpath, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Saved {outpath.name}")


# ── Load linkage matrices for target patients ────────────────────────
print("Loading LRG linkage matrices...")
linkages = {}
for pat in PATIENT_PHASES:
    for phase in PATIENT_PHASES[pat]:
        for band in BANDS:
            try:
                res = load_lrg_result(pat, phase, band, fc_method="msc",
                                      cache_root=LRG_CACHE)
                linkages[(pat, band, phase)] = res.linkage_matrix
            except Exception:
                pass
print(f"Loaded {len(linkages)} linkages")

# ── Generate figures ─────────────────────────────────────────────────
for pat, phases in PATIENT_PHASES.items():
    print(f"\nGenerating affinity figure for {pat}...")
    generate_affinity_figure(pat, phases, linkages)

print("\nDone!")
