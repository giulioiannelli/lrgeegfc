#!/usr/bin/env python3
r"""fig_reinstatement_scales_baseline -- the rho_sym reinstatement at THREE scales, FAITHFUL colouring.

Fall-back sibling of fig_reinstatement_scales_measure.py (the contrast-coloured slide). SAME
layout -- four phase dendrograms per scale row (rest_pre -> task_learn -> task_test -> rest_post),
stacked at three diffusion scales s (fine / report / coarse) -- but here the per-contact COLOUR
IS the measure itself: the per-leaf projection of the SYMMETRIC, split-baseline-subtracted
cophenetic concordance rho_sym (== rho^coph; raw non-symmetric phase->test form retired 2026-07-13):

    rho_sym(X <-> test) = 1/2[ Spearman(D_test-D_preA, D_X-D_preB)
                             + Spearman(D_test-D_preB, D_X-D_preA) ]

with A,B = the two rest_pre split-halves. Per-leaf, the same expression is evaluated on leaf i's
cophenetic ROW (i,:) instead of the full upper triangle -> a per-contact rho_sym in [-1,1] whose
whole-graph average is the bar value. rest_pre is the SUBTRACTED baseline (no shift -> 0 by
construction), rendered a clean neutral grey; task_learn / task_test / rest_post light up as their
contacts' hierarchy rows come to resemble the task-shift. The picture and the on-slide formula are
the SAME object -- nothing to reconcile.

Per column bars (identical to the contrast sibling; the locked pre-referenced gate rho_sym):
  - rest_pre   = baseline (0, grey);
  - task_learn = rho_sym(learn-shift <-> test-shift)   (encoding concordance);
  - task_test  = the REFERENCE (fills to the reachable ceiling);
  - rest_post  = rho_sym(post-shift <-> test-shift)    (the reinstatement GATE value,
                 == data/sparsified_arc/ms_mst020/per_patient_scale.csv obs_rho).
Each bar is cut at the REACHABLE ceiling rel = Spearman(D_A, D_B) (model-free split-half
reliability, ~0.5 for beta) -- NEVER 1 -- and marks the matched-strength null (dashed white line
at surrogate p95). The unreachable region is not drawn. Per-contact colour = a SEQUENTIAL grey
-> task-green ramp on the per-leaf rho_sym: grey = baseline (0; anti folds into baseline), green
= trace (reinstated toward the task hierarchy). Scale limits [0, 0.75] fit the per-leaf spread
(medians learn ~0.4 / post ~0.5 / task_test ~0.65) so the greens resolve rather than saturate.

Reads : imcoh_abs FC (Pat_08, four phases + rest_pre split-halves A,B) via load_fc_matrix
        / load_phase + mst@0.20 backbone; per_patient_scale.csv for obs_rho + surr_p95.
Writes: data/outputs/figures/talk/slide_measure_reinstatement_scales_baseline.{pdf,png}
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

# per-leaf rho_sym: sequential grey (baseline / anti) -> task-green (trace). Endpoints reuse
# the bar greens: #8fca90 = task_test reference, #136b32 = rest_post (strongest reinstatement).
SEQ_CM = LinearSegmentedColormap.from_list("seq", ["#c9c9c9", "#8fca90", "#136b32"])
NORM = Normalize(0.0, 0.75)
C_PRE, C_LEARN, C_TEST, C_POST = "0.72", "#6f8fb0", "#8fca90", "#136b32"


def _fc(W):
    W = np.asarray(W, float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


def per_leaf_rho_sym(DX, DT, DA, DB):
    """Per-leaf symmetric rho_sym projection: leaf i's task-shift row vs test-shift row."""
    n = DX.shape[0]
    out = np.zeros(n)
    for i in range(n):
        m = np.ones(n, bool)
        m[i] = False
        a = spearmanr(DT[i, m] - DA[i, m], DX[i, m] - DB[i, m]).correlation
        b = spearmanr(DT[i, m] - DB[i, m], DX[i, m] - DA[i, m]).correlation
        out[i] = 0.5 * (np.nan_to_num(a) + np.nan_to_num(b))
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
    print(f"fig_reinstatement_scales_baseline ({PATIENT} {BAND}, mst@0.20, faithful rho_sym):")
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
        # per-contact COLOUR = the measure itself, resolved per leaf. rest_pre = subtracted
        # baseline -> uniform grey (0); other phases light up as their rows resemble task-shift.
        pl = {ph: (np.zeros(DT.shape[0]) if ph == "rest_pre"
                   else per_leaf_rho_sym(D[ph], DT, DA, DB)) for ph in COLS}
        leafc = {ph: [SEQ_CM(NORM(np.clip(np.nan_to_num(v), 0.0, 0.75))) for v in pl[ph]]
                 for ph in COLS}
        print(f"  s={s:6.2f}: rel={rel:.2f} null_p95={null:.2f}  "
              f"learn {vals['task_learn']:+.2f} (leaf {np.nanmean(pl['task_learn']):+.2f})  "
              f"post {vals['rest_post']:+.2f} (leaf {np.nanmean(pl['rest_post']):+.2f})")

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

    # vertical sequential colourbar (grey baseline -> task-green trace; no label), full height
    cax = fig.add_axes([0.955, 0.025, 0.014, 0.94])
    cb = fig.colorbar(ScalarMappable(norm=NORM, cmap=SEQ_CM), cax=cax, orientation="vertical")
    cb.set_ticks([0.0, 0.25, 0.5, 0.75])
    cb.ax.tick_params(labelsize=9)
    cb.outline.set_visible(False)

    outdir = FIGURES_ROOT / "talk"
    outdir.mkdir(parents=True, exist_ok=True)
    stem = outdir / "slide_measure_reinstatement_scales_baseline"
    fig.savefig(f"{stem}.pdf", transparent=True, bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {stem}.pdf and {stem}.png")


if __name__ == "__main__":
    main()
