#!/usr/bin/env python3
r"""fig_reinstatement_contrast_coloring -- CONCEPT TEST: colour rest_pre too, via the pre<->post contrast.

Alternative to the "rest_pre = grey baseline" colouring. Instead of colouring each phase by its
own resemblance to task (raw, confounded by the shared anatomical skeleton) or by rho_sym
(which makes rest_pre 0 by construction), colour every contact by the per-contact REINSTATEMENT
CONTRAST -- how much MORE its cophenetic row resembles task in rest_post than in rest_pre:

    sim_X(i)  = Spearman( D_X[i, :], D_test[i, :] )          (leaf i's row resemblance to task)
    Delta(i)  = sim_post(i) - sim_pre(i)                     (the shared skeleton CANCELS)

Painted:
    rest_post  ->  +Delta   (green  = contact moved TOWARD the task hierarchy)
    rest_pre   ->  -Delta   (purple = contact STARTED below task, by that much)   <-- pre now carries signal
    task_learn ->  sim_learn - sim_pre   (its own progress from baseline)
    task_test  ->  sim_test  - sim_pre = 1 - sim_pre   (the reference: full gain)

So the pre->post colour swing per contact IS the trace. NOTE pre's purple = "below task relative
to post", NOT an anti-correlation. The BARS stay the locked pre-referenced gate rho_sym.

Single scale s=5.6 (the reporting scale), Pat_08 beta.
Writes: data/outputs/figures/talk/slide_measure_reinstatement_contrast.{pdf,png}
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.stats import spearmanr

import _common as C
from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.workflow.fc import load_fc_matrix
from audit_150_rho_sym_gate import load_phase

C.use_lrg_style()

PATIENT, BAND, FC = "Pat_08", "beta", "imcoh_abs"
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
TITLES = {"rest_pre": r"$\mathrm{rest_{pre}}$", "task_learn": r"$\mathrm{task_{learn}}$",
          "task_test": r"$\mathrm{task_{test}}$", "rest_post": r"$\mathrm{rest_{post}}$"}
# contrast map kept OFF the bar-green: red = below/different from task, blue = reinstated/same
DIV_CM = LinearSegmentedColormap.from_list("div", ["#88314b", "#c9c9c9", "#343171"])
NORM = Normalize(-0.5, 0.5)
C_LEARN, C_TEST, C_POST_BAR, C_PRE = "#6f8fb0", "#8fca90", "#136b32", "0.72"


def _fc(W):
    W = np.asarray(W, float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


def per_leaf_resemblance(DX, DT):
    """Per-leaf Spearman of leaf i's cophenetic row: phase X vs task_test (excluding self)."""
    n = DX.shape[0]
    out = np.zeros(n)
    for i in range(n):
        m = np.ones(n, bool)
        m[i] = False
        out[i] = np.nan_to_num(spearmanr(DX[i, m], DT[i, m]).correlation)
    return out


def _bar(ax, val, rel, null, color, caption):
    ax.barh([0], [rel], color="0.9", height=0.55, zorder=1)
    if val is not None and val > 0:
        ax.barh([0], [min(val, rel)], color=color, height=0.55, zorder=2)
    ax.plot([null, null], [-0.30, 0.30], color="white", lw=1.5, ls=(0, (3.5, 2.8)),
            zorder=3, solid_capstyle="butt")
    ax.set_xlim(0, rel)
    ax.set_ylim(-1.5, 0.7)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.text(rel * 0.5, -0.82, caption, ha="center", va="top", fontsize=9.5, color="0.25")


def main():
    from lrg_eegfc.visuals.lrg import plot_circular_dendrogram

    s0 = float(C.SGRID[C.I_REPORT])
    FCm = {ph: _fc(load_fc_matrix(PATIENT, ph, BAND, fc_method=FC)) for ph in PHASES}
    WA, WB = _fc(load_phase(PATIENT, "A", BAND)), _fc(load_phase(PATIENT, "B", BAND))
    Z = {ph: C.tree_at_scale(FCm[ph]) for ph in PHASES}
    D = {ph: C.coph_square_at_scale(FCm[ph]) for ph in PHASES}
    DA, DB, DT = C.coph_square_at_scale(WA), C.coph_square_at_scale(WB), D["task_test"]
    iu = np.triu_indices(DT.shape[0], 1)

    # --- per-contact reinstatement contrast (colour) ---
    sim = {ph: per_leaf_resemblance(D[ph], DT) for ph in PHASES}
    delta = sim["rest_post"] - sim["rest_pre"]
    cval = {"rest_pre": -delta,                              # mirror: started below task
            "task_learn": sim["task_learn"] - sim["rest_pre"],
            "task_test": sim["task_test"] - sim["rest_pre"],  # = 1 - sim_pre (reference)
            "rest_post": delta}
    leaf_colors = {ph: [DIV_CM(NORM(np.clip(v, -0.5, 0.5))) for v in cval[ph]] for ph in PHASES}

    # --- gate rho_sym for the bars (locked, pre-referenced) ---
    rel = float(spearmanr(DA[iu], DB[iu]).correlation)
    pp = pd.read_csv(C.MS / "per_patient_scale.csv")
    xrow = pp[(pp.patient == PATIENT) & (pp.band == BAND) & np.isclose(pp.s, s0)].iloc[0]
    null = float(xrow.surr_p95)

    def rho_sym(DX):
        a = spearmanr(DT[iu] - DA[iu], DX[iu] - DB[iu]).correlation
        b = spearmanr(DT[iu] - DB[iu], DX[iu] - DA[iu]).correlation
        return 0.5 * (a + b)

    bar = {"rest_pre": 0.0, "task_learn": rho_sym(D["task_learn"]),
           "task_test": rel, "rest_post": float(xrow.obs_rho)}

    fig = plt.figure(figsize=(14.6, 5.3))
    gs = fig.add_gridspec(2, 4, height_ratios=[3.4, 0.5], hspace=0.06, wspace=0.06,
                          left=0.02, right=0.99, top=0.9, bottom=0.17)
    node_s = max(5.0, 820.0 / (Z["rest_pre"].shape[0] + 1))
    for j, ph in enumerate(PHASES):
        ax = fig.add_subplot(gs[0, j])
        plot_circular_dendrogram(ax, C.log_height_linkage(Z[ph]), leaf_colors=leaf_colors[ph],
                                 node_size=node_s, line_width=1.05, r_inner=0.10)
        ax.set_title(TITLES[ph], fontsize=17, pad=8)
        axb = fig.add_subplot(gs[1, j])
        if ph == "task_test":
            _bar(axb, rel, rel, null, C_TEST, "(reference)")
        elif ph == "rest_pre":
            _bar(axb, 0.0, rel, null, C_PRE, "baseline")
        else:
            col = C_POST_BAR if ph == "rest_post" else C_LEARN
            _bar(axb, bar[ph], rel, null, col, rf"$\rho_{{\mathrm{{sym}}}} = {bar[ph]:.2f}$")

    cax = fig.add_axes([0.32, 0.055, 0.36, 0.02])
    cb = fig.colorbar(ScalarMappable(norm=NORM, cmap=DIV_CM), cax=cax, orientation="horizontal")
    cb.set_label(r"per-contact reinstatement  $\Delta = \mathrm{sim}(\mathrm{post}) - "
                 r"\mathrm{sim}(\mathrm{pre})$ to task   "
                 r"($\mathrm{pre}=-\Delta$: below $\leftrightarrow$ post$=+\Delta$: reinstated)",
                 fontsize=9)
    cb.set_ticks([-0.5, 0.0, 0.5])
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_visible(False)

    print(f"fig_reinstatement_contrast_coloring ({PATIENT} {BAND}, s={s0:.2f}):")
    print(f"  per-leaf sim->task:  pre {sim['rest_pre'].mean():+.2f}  learn "
          f"{sim['task_learn'].mean():+.2f}  post {sim['rest_post'].mean():+.2f}")
    print(f"  contrast Delta mean {delta.mean():+.2f}  (frac>0 {np.mean(delta > 0):.2f})")

    outdir = FIGURES_ROOT / "talk"
    outdir.mkdir(parents=True, exist_ok=True)
    stem = outdir / "slide_measure_reinstatement_contrast"
    fig.savefig(f"{stem}.pdf", transparent=True, bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {stem}.pdf and {stem}.png")


if __name__ == "__main__":
    main()
