#!/usr/bin/env python3
r"""fig_reinstatement_scales_measure -- the SAME 4-phase reinstatement row, read at THREE scales.

Measure-slide variant of ``fig_reinstatement_learn_mst020``. The bottom cohort forest and
the shared per-leaf colourbar are STRIPPED; what remains is the row of four circular
dendrograms (rest_pre -> task_learn -> task_test -> rest_post), each coloured by per-leaf
cophenetic preservation to task_test and topped by a ``rho^coph``-to-test bar -- and that
row is stacked at three diffusion scales s (fine / report / coarse).

The bar anchor is the SAME structure at every scale: full bar = identical to the task_test
hierarchy (rho^coph = 1). So a fill level reads the same way across the three rows -- the
point of the measure slide: ONE cophenetic concordance, read at several scales, and the
trace (pre pale -> learn warms -> post stays warm) survives all of them.

Scales are chosen inside the scale-robust band (Pat_08 beta): pre ~ 0, learn ~ 0.3,
post ~ 0.4 at each -- so the reading is directly interpretable, not a coarse-scale artifact.

Reads : imcoh_abs FC (Pat_08, four phases, via load_fc_matrix + mst@0.20 backbone).
Writes: data/outputs/figures/talk/slide_measure_reinstatement_scales.{pdf,png}
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.stats import spearmanr

import _common as C
from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.workflow.fc import load_fc_matrix

C.use_lrg_style()

PATIENT, BAND, FC = "Pat_08", "beta", "imcoh_abs"
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
TITLES = {"rest_pre": r"$\mathrm{rest_{pre}}$", "task_learn": r"$\mathrm{task_{learn}}$",
          "task_test": r"$\mathrm{task_{test}}$", "rest_post": r"$\mathrm{rest_{post}}$"}
# three scales inside the robust band (see the 16-scale sweep): fine / report / coarse
SCALE_IDX = [0, 5, 10]                             # SGRID -> 1.00, 5.65, 31.88
SCALES = [float(C.SGRID[i]) for i in SCALE_IDX]
SCALE_TAG = ["fine", "report", "coarse"]

PRES_CM = LinearSegmentedColormap.from_list("pres", ["#d7d7d7", "#7cc47f", "#136b32"])
NORM = Normalize(0.0, 0.8)
C_PRE, C_LEARN, C_POST = "0.55", "#6f8fb0", "#3a9a4f"


def phase_fc(phase):
    """imcoh_abs FC with the canonical load_phase preprocessing (diag 0, clip [0,1], symm)."""
    W = np.asarray(load_fc_matrix(PATIENT, phase, BAND, fc_method=FC), float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


def per_leaf_preservation(A, B):
    """For each leaf i, Spearman of its cophenetic-distance row (A vs B), excluding self."""
    n = A.shape[0]
    out = np.zeros(n)
    for i in range(n):
        m = np.ones(n, bool)
        m[i] = False
        out[i] = spearmanr(A[i, m], B[i, m]).correlation
    return out


def _sim_bar(ax, val, color, caption):
    """Bar filling: background = the reference structure (full = rho^coph 1)."""
    ax.barh([0], [1.0], color="0.90", height=0.55, zorder=1)
    if val is not None:
        ax.barh([0], [max(val, 0.0)], color=color, height=0.55, zorder=2)
    ax.set_xlim(0, 1)
    ax.set_ylim(-1.5, 0.7)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(0.5, -0.9, caption, ha="center", va="top", fontsize=9.5, color="0.25")


def main():
    from lrg_eegfc.visuals.lrg import plot_circular_dendrogram

    FCm = {ph: phase_fc(ph) for ph in PHASES}

    fig = plt.figure(figsize=(14.6, 12.2))
    gs = fig.add_gridspec(6, 4, height_ratios=[3.4, 0.5, 3.4, 0.5, 3.4, 0.5],
                          hspace=0.18, wspace=0.06, left=0.055, right=0.99,
                          top=0.965, bottom=0.02)
    node_s = None
    print(f"fig_reinstatement_scales_measure ({PATIENT} {BAND}, mst@0.20):")
    for r, (si, s, tag) in enumerate(zip(SCALE_IDX, SCALES, SCALE_TAG)):
        D = {ph: C.coph_square_at_scale(FCm[ph], s) for ph in PHASES}
        Z = {ph: C.tree_at_scale(FCm[ph], s) for ph in PHASES}
        pl = {ph: per_leaf_preservation(D[ph], D["task_test"]) for ph in PHASES}
        sim = {ph: float(spearmanr(D[ph][np.triu_indices(D[ph].shape[0], 1)],
                                   D["task_test"][np.triu_indices(D[ph].shape[0], 1)]).correlation)
               for ph in PHASES}
        if node_s is None:
            node_s = max(5.0, 760.0 / (Z["rest_pre"].shape[0] + 1))
        print(f"  s={s:6.2f} ({tag:6s}): "
              + "  ".join(f"{ph.split('_')[-1]} {sim[ph]:+.2f}" for ph in PHASES))
        for j, ph in enumerate(PHASES):
            ax = fig.add_subplot(gs[2 * r, j])
            leaf_colors = [PRES_CM(NORM(np.clip(np.nan_to_num(v, nan=0.0), 0, 0.8)))
                           for v in pl[ph]]
            plot_circular_dendrogram(ax, C.log_height_linkage(Z[ph]),
                                     leaf_colors=leaf_colors, node_size=node_s,
                                     line_width=1.05, r_inner=0.10)
            if r == 0:
                ax.set_title(TITLES[ph], fontsize=17, pad=8)
            if j == 0:
                ax.text(-0.06, 0.5, rf"$s = {s:.1f}$", transform=ax.transAxes,
                        rotation=90, ha="right", va="center", fontsize=15,
                        fontweight="bold", color="0.2")
            axb = fig.add_subplot(gs[2 * r + 1, j])
            if ph == "task_test":
                _sim_bar(axb, 1.0, "#3a9a4f", "(reference)")
            else:
                col = C_POST if ph == "rest_post" else (C_PRE if ph == "rest_pre" else C_LEARN)
                _sim_bar(axb, sim[ph], col, rf"$\rho^{{\mathrm{{coph}}}}$ to test = {sim[ph]:.2f}")

    outdir = FIGURES_ROOT / "talk"
    outdir.mkdir(parents=True, exist_ok=True)
    stem = outdir / "slide_measure_reinstatement_scales"
    fig.savefig(f"{stem}.pdf", transparent=True, bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {stem}.pdf and {stem}.png")


if __name__ == "__main__":
    main()
