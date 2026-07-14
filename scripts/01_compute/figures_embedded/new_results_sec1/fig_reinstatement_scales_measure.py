#!/usr/bin/env python3
r"""fig_reinstatement_scales_measure -- the rho_sym reinstatement, read at THREE scales.

Measure-slide figure. Four phase dendrograms per scale row (rest_pre -> task_learn ->
task_test -> rest_post), stacked at three diffusion scales s (fine / report / coarse).
Everything is the SYMMETRIC, split-baseline-subtracted cophenetic concordance rho_sym
(== rho^coph; the raw non-symmetric phase->test form is retired, 2026-07-13):

    rho_sym(X <-> test) = 1/2[ Spearman(D_test-D_preA, D_X-D_preB)
                             + Spearman(D_test-D_preB, D_X-D_preA) ]

with A,B = the two rest_pre split-halves. Per column:
  - rest_pre   = the subtracted BASELINE (no shift -> rho_sym is 0 by construction);
  - task_learn = rho_sym(learn-shift <-> test-shift)   (encoding concordance);
  - task_test  = the REFERENCE (its bar fills to the reachable ceiling);
  - rest_post  = rho_sym(post-shift <-> test-shift)    (the reinstatement GATE value,
                 identical to data/sparsified_arc/ms_mst020/per_patient_scale.csv obs_rho).

Each bar is cut at the REACHABLE ceiling rel = Spearman(D_A, D_B) (the model-free
split-half reliability, ~0.5 for beta) -- NEVER 1 -- and marks the matched-strength
null (dashed white line at the surrogate p95). The unreachable region is not drawn.

Per-contact COLOUR = the reinstatement CONTRAST (so rest_pre is not blank). With per-leaf
resemblance sim_X(i) = Spearman(D_X[i,:], D_test[i,:]) and Delta(i) = sim_post(i) - sim_pre(i)
(the shared skeleton cancels): rest_post -> +Delta (green, moved toward task); rest_pre ->
-Delta (purple, started below task); task_learn -> sim_learn - sim_pre; task_test -> 1 - sim_pre
(reference). The pre->post colour flip per contact IS the trace. pre-purple = "below task
relative to post", not an anti-correlation.

Reads : imcoh_abs FC (Pat_08, four phases + rest_pre split-halves A,B) via load_fc_matrix
        / load_phase + mst@0.20 backbone; per_patient_scale.csv for obs_rho + surr_p95.
Writes: data/outputs/figures/talk/slide_measure_reinstatement_scales.{pdf,png}
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
COLS = ["rest_pre", "task_learn", "task_test", "rest_post"]
TITLES = {"rest_pre": r"$\mathrm{rest_{pre}}$", "task_learn": r"$\mathrm{task_{learn}}$",
          "task_test": r"$\mathrm{task_{test}}$", "rest_post": r"$\mathrm{rest_{post}}$"}
SCALE_IDX = [0, 5, 10]                              # SGRID -> 1.00, 5.65, 31.88
SCALES = [float(C.SGRID[i]) for i in SCALE_IDX]

# per-leaf rho_sym: anti (purple) <-> 0 (grey) <-> trace (green)
# contrast map kept OFF the bar-green: red = below/different from task, blue = reinstated/same
DIV_CM = LinearSegmentedColormap.from_list("div", ["#88314b", "#c9c9c9", "#343171"])
NORM = Normalize(-0.5, 0.5)
C_PRE, C_LEARN, C_TEST, C_POST = "0.72", "#6f8fb0", "#8fca90", "#136b32"


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
    """Horizontal bar cut at the reachable ceiling rel; dashed white null; no unreachable end."""
    ax.barh([0], [rel], color="0.9", height=0.55, zorder=1)          # reachable region -> ends at ceiling
    if val is not None and val > 0:
        ax.barh([0], [min(val, rel)], color=color, height=0.55, zorder=2)
    ax.plot([null, null], [-0.30, 0.30], color="white", lw=1.5, ls=(0, (3.5, 2.8)),
            zorder=3, solid_capstyle="butt")                        # dashed white null (matched-strength p95)
    ax.set_xlim(0, rel)
    ax.set_ylim(-1.5, 0.7)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.text(rel * 0.5, -0.82, caption, ha="center", va="top", fontsize=9.5, color="0.25")


def main():
    from lrg_eegfc.visuals.lrg import plot_circular_dendrogram

    pp = pd.read_csv(C.MS / "per_patient_scale.csv")
    FCm = {ph: _fc(load_fc_matrix(PATIENT, ph, BAND, fc_method=FC)) for ph in COLS}
    WA, WB = _fc(load_phase(PATIENT, "A", BAND)), _fc(load_phase(PATIENT, "B", BAND))

    fig = plt.figure(figsize=(14.6, 12.4))
    gs = fig.add_gridspec(6, 4, height_ratios=[3.4, 0.5, 3.4, 0.5, 3.4, 0.5],
                          hspace=0.22, wspace=0.06, left=0.055, right=0.94,
                          top=0.965, bottom=0.025)
    node_s, iu = None, None
    print(f"fig_reinstatement_scales_measure ({PATIENT} {BAND}, mst@0.20, rho_sym):")
    for r, s in enumerate(SCALES):
        D = {ph: C.coph_square_at_scale(FCm[ph], s) for ph in COLS}
        DA, DB, DT = (C.coph_square_at_scale(W, s) for W in (WA, WB, FCm["task_test"]))
        Z = {ph: C.tree_at_scale(FCm[ph], s) for ph in COLS}
        if iu is None:
            iu = np.triu_indices(DT.shape[0], 1)
            node_s = max(5.0, 760.0 / (Z["rest_pre"].shape[0] + 1))
        rel = float(spearmanr(DA[iu], DB[iu]).correlation)
        row = pp[(pp.patient == PATIENT) & (pp.band == BAND) & np.isclose(pp.s, s)].iloc[0]
        null = float(row.surr_p95)

        def rho_sym(DX):
            a = spearmanr(DT[iu] - DA[iu], DX[iu] - DB[iu]).correlation
            b = spearmanr(DT[iu] - DB[iu], DX[iu] - DA[iu]).correlation
            return 0.5 * (a + b)

        vals = {"rest_pre": 0.0,                        # baseline: no shift by construction
                "task_learn": rho_sym(D["task_learn"]),
                "task_test": rel,                       # reference -> fills to ceiling
                "rest_post": float(row.obs_rho)}        # gate value (== rho_sym(D["rest_post"]))
        # per-contact reinstatement contrast (colour): Delta = sim(post) - sim(pre) to task,
        # skeleton cancels; pre painted -Delta (below), post +Delta (reinstated).
        sim = {ph: per_leaf_resemblance(D[ph], DT) for ph in COLS}
        delta = sim["rest_post"] - sim["rest_pre"]
        cval = {"rest_pre": -delta,
                "task_learn": sim["task_learn"] - sim["rest_pre"],
                "task_test": sim["task_test"] - sim["rest_pre"],   # = 1 - sim_pre (reference)
                "rest_post": delta}
        leafc = {ph: [DIV_CM(NORM(np.clip(v, -0.5, 0.5))) for v in cval[ph]] for ph in COLS}
        print(f"  s={s:6.2f}: rel={rel:.2f} null_p95={null:.2f}  "
              f"learn {vals['task_learn']:+.2f}  post {vals['rest_post']:+.2f}")

        for j, ph in enumerate(COLS):
            ax = fig.add_subplot(gs[2 * r, j])
            plot_circular_dendrogram(ax, C.log_height_linkage(Z[ph]), leaf_colors=leafc[ph],
                                     node_size=node_s, line_width=1.05, r_inner=0.10)
            if r == 0:
                ax.set_title(TITLES[ph], fontsize=17, pad=8)
            if j == 0:
                ax.text(-0.06, 0.5, rf"$s = {s:.1f}$", transform=ax.transAxes, rotation=90,
                        ha="right", va="center", fontsize=15, fontweight="bold", color="0.2")
            axb = fig.add_subplot(gs[2 * r + 1, j])
            if ph == "task_test":
                _bar(axb, rel, rel, null, C_TEST, "(reference)")
            elif ph == "rest_pre":
                _bar(axb, 0.0, rel, null, C_PRE, "baseline")
            else:
                col = C_POST if ph == "rest_post" else C_LEARN
                _bar(axb, vals[ph], rel, null, col, rf"$\rho_{{\mathrm{{sym}}}} = {vals[ph]:.2f}$")

    # vertical diverging colourbar (no label)
    cax = fig.add_axes([0.955, 0.30, 0.014, 0.40])
    cb = fig.colorbar(ScalarMappable(norm=NORM, cmap=DIV_CM), cax=cax, orientation="vertical")
    cb.set_ticks([-0.5, 0.0, 0.5])
    cb.ax.tick_params(labelsize=9)
    cb.outline.set_visible(False)

    outdir = FIGURES_ROOT / "talk"
    outdir.mkdir(parents=True, exist_ok=True)
    stem = outdir / "slide_measure_reinstatement_scales"
    fig.savefig(f"{stem}.pdf", transparent=True, bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {stem}.pdf and {stem}.png")


if __name__ == "__main__":
    main()
