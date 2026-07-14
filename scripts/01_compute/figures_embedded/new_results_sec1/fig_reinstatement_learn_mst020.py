#!/usr/bin/env python3
r"""fig_reinstatement_learn_mst020 -- the resting hierarchy becomes task-like AFTER the task,
across the full rest -> learn -> test -> rest arc (mst@0.20).

HEAD: the post-task resting cophenetic hierarchy re-expresses the task hierarchy; the pre-task
one does not -- rebuilt on the **mst@0.20** backbone at the reporting scale **s=5.6**, resolved
across the FULL paradigm rest_pre -> task_learn -> task_test -> rest_post. Everything is the
SYMMETRIC, split-baseline-subtracted concordance rho_sym (== rho^coph; the raw non-symmetric
phase->test form is retired, 2026-07-13):

    rho_sym(X <-> test) = 1/2[ Spearman(D_test-D_preA, D_X-D_preB)
                             + Spearman(D_test-D_preB, D_X-D_preA) ]

with A,B = the two rest_pre split-halves. Per column: rest_pre = the subtracted BASELINE (no
shift -> 0 by construction, rendered neutral grey); task_learn = rho_sym(learn-shift <-> test-
shift), the hierarchy already forming; task_test = the REFERENCE (its bar fills to the reachable
ceiling); rest_post = rho_sym(post-shift <-> test-shift), the reinstatement GATE value. Each bar
is cut at the REACHABLE ceiling rel = Spearman(D_A, D_B) (the model-free split-half reliability
~0.5 for beta) -- NEVER 1 -- and marks the matched-strength null (dashed white line at surrogate
p95). Per-contact colour = the per-leaf rho_sym projection (diverging: anti <-> trace).

Exemplar Pat_08 (beta, a strong tracer: whole-phase beta rho_sym = +0.48 at s=5.6). The cohort
claim (bottom) is the load-bearing one: per-patient beta rho_sym (task_test <-> rest_post) at
s=5.6 against the matched-strength null -- the same gate that clears 16/16 scales at p=0.001.
Read on the mst@0.20 multiscale backbone, not the degenerate dense graph.

Reads : imcoh_abs FC (Pat_08, four phases + rest_pre split-halves A,B) via load_fc_matrix /
        load_phase + mst@0.20 backbone; per_patient_scale.csv (cohort beta rho_sym @ s=5.6).
Writes: data/preprint/figures/new_results_sec1/fig_reinstatement_learn_mst020.pdf
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
from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.visuals.lrg import plot_circular_dendrogram
from lrg_eegfc.workflow.fc import load_fc_matrix
from audit_150_rho_sym_gate import load_phase

C.use_lrg_style()

PATIENT, BAND, FC = "Pat_08", "beta", "imcoh_abs"
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
TITLES = {"rest_pre": r"$\mathrm{rest_{pre}}$", "task_learn": r"$\mathrm{task_{learn}}$",
          "task_test": r"$\mathrm{task_{test}}$", "rest_post": r"$\mathrm{rest_{post}}$"}
# per-leaf rho_sym: anti (purple) <-> 0 (grey) <-> trace (green)
DIV_CM = LinearSegmentedColormap.from_list("div", ["#6a3d9a", "#c9c9c9", "#136b32"])
NORM = Normalize(-0.5, 0.5)
C_LEARN, C_TEST, C_POST_BAR, C_PRE = "#6f8fb0", "#8fca90", "#136b32", "0.72"
C_DOT, C_ANTI = "#3a9a4f", "#d1352b"          # cohort-panel dot colours
OUT = C.FIGDIR / "fig_reinstatement_learn_mst020.pdf"


def _fc(W):
    """imcoh_abs FC with canonical preprocessing (diag 0, clip [0,1], symmetric)."""
    W = np.asarray(W, float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


def phase_fc(phase):
    return _fc(load_fc_matrix(PATIENT, phase, BAND, fc_method=FC))


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
    ax.text(rel * 0.5, -0.82, caption, ha="center", va="top", fontsize=10.0, color="0.25")


def draw(target):
    s0 = float(C.SGRID[C.I_REPORT])
    FCm = {ph: phase_fc(ph) for ph in PHASES}
    WA, WB = _fc(load_phase(PATIENT, "A", BAND)), _fc(load_phase(PATIENT, "B", BAND))
    Z = {ph: C.tree_at_scale(FCm[ph]) for ph in PHASES}
    D = {ph: C.coph_square_at_scale(FCm[ph]) for ph in PHASES}
    DA, DB, DT = C.coph_square_at_scale(WA), C.coph_square_at_scale(WB), D["task_test"]
    iu = np.triu_indices(DT.shape[0], 1)
    rel = float(spearmanr(DA[iu], DB[iu]).correlation)

    pp = pd.read_csv(C.MS / "per_patient_scale.csv")
    xrow = pp[(pp.patient == PATIENT) & (pp.band == BAND) & np.isclose(pp.s, s0)].iloc[0]
    null = float(xrow.surr_p95)

    def rho_sym(DX):
        a = spearmanr(DT[iu] - DA[iu], DX[iu] - DB[iu]).correlation
        b = spearmanr(DT[iu] - DB[iu], DX[iu] - DA[iu]).correlation
        return 0.5 * (a + b)

    vals = {"rest_pre": 0.0,                          # baseline: no shift by construction
            "task_learn": rho_sym(D["task_learn"]),
            "task_test": rel,                         # reference -> fills to ceiling
            "rest_post": float(xrow.obs_rho)}         # gate value (== rho_sym(D["rest_post"]))
    pl = {ph: (np.zeros(DT.shape[0]) if ph == "rest_pre"
               else per_leaf_rho_sym(D[ph], DT, DA, DB)) for ph in PHASES}
    leaf_colors = {ph: [DIV_CM(NORM(np.clip(np.nan_to_num(v), -0.5, 0.5))) for v in pl[ph]]
                   for ph in PHASES}

    rdf = load_channel_regions(PATIENT)
    rdf = rdf.loc[:, ~rdf.columns.duplicated()]
    leaf_labels = rdf["label"].astype(str).tolist()

    gs = target.add_gridspec(4, 4, height_ratios=[3.4, 0.42, 0.5, 1.15], hspace=0.28,
                             wspace=0.05, left=0.03, right=0.965, top=0.975, bottom=0.085)
    node_s = max(6.0, 900.0 / (Z["rest_pre"].shape[0] + 1))
    for j, ph in enumerate(PHASES):
        ax = target.add_subplot(gs[0, j])
        # DISPLAY at log-height: mst@0.20 cophenetic heights span ~4-8 decades, so a
        # linear-radius fan dendrogram crushes ~95-99% of merges onto the rim (only 2-3
        # splits show). log_height_linkage spreads the hierarchy over the radius; topology
        # + rho_sym colours are unchanged (see _common.log_height_linkage).
        plot_circular_dendrogram(ax, C.log_height_linkage(Z[ph]), leaf_colors=leaf_colors[ph],
                                 leaf_labels=leaf_labels, node_size=node_s,
                                 line_width=1.05, r_inner=0.10,
                                 label_fontsize=3.0, label_offset=1.03)
        ax.set_title(TITLES[ph], fontsize=16, pad=6)
        axb = target.add_subplot(gs[1, j])
        if ph == "task_test":
            _bar(axb, rel, rel, null, C_TEST, "(reference)")
        elif ph == "rest_pre":
            _bar(axb, 0.0, rel, null, C_PRE, "baseline")
        else:
            col = C_POST_BAR if ph == "rest_post" else C_LEARN
            _bar(axb, vals[ph], rel, null, col, rf"$\rho_{{\mathrm{{sym}}}} = {vals[ph]:.2f}$")

    target.text(0.03, 0.955, r"$\mathbf{a}$", fontsize=17, va="bottom", ha="left",
                fontweight="bold")

    # horizontal colourbar (leaf colour = per-contact rho_sym to task_test)
    band = gs[2, :].get_position(target)
    cb_w, cb_h = 0.34, 0.022
    cax = target.add_axes([0.5 - cb_w / 2, band.y0 + 0.5 * band.height, cb_w, cb_h])
    cb = target.colorbar(ScalarMappable(norm=NORM, cmap=DIV_CM), cax=cax,
                         orientation="horizontal")
    cb.set_label(r"per-contact $\rho_{\mathrm{sym}}$ to $\mathrm{task_{test}}$"
                 r"  (anti $\leftrightarrow$ trace)", fontsize=10)
    cb.set_ticks([-0.5, 0.0, 0.5])
    cb.ax.tick_params(labelsize=9)
    cb.solids.set_rasterized(False)

    # -- bottom: cohort beta rho_sym (task<->post) at s=5.6 vs matched-strength null --
    axs = target.add_subplot(gs[3, :])
    x = pp[(pp.band == BAND) & np.isclose(pp.s, s0)].set_index("patient")
    obs = np.array([float(x.loc[p, "obs_rho"]) for p in C.COHORT])
    p95 = np.array([float(x.loc[p, "surr_p95"]) for p in C.COHORT])
    clears = np.array([float(x.loc[p, "p"]) < 0.05 for p in C.COHORT])
    y = np.arange(len(C.COHORT))
    axs.axvline(0.0, color="0.55", lw=1.0, ls="--", zorder=1)
    for yi, o, t, cl in zip(y, obs, p95, clears):
        axs.plot([0, o], [yi, yi], color=C_DOT if o >= 0 else C_ANTI, lw=1.3, zorder=2)
        axs.plot([t], [yi], "|", color="0.5", ms=12, mew=1.8, zorder=3)   # own null p95
        axs.scatter(o, yi, s=70, color=C_DOT if o >= 0 else C_ANTI,
                    edgecolor="black" if cl else "white", lw=1.2 if cl else 0.6, zorder=4)
    npos = int((obs > 0).sum())
    axs.text(0.985, 0.90, rf"$\beta$: {npos}/10 patients reinstate  ·  cohort gate "
             r"$p=0.001$ (16/16 scales)", transform=axs.transAxes, ha="right", va="center",
             fontsize=12, color=C_DOT, fontweight="bold")
    axs.text(0.015, 0.90, r"cohort, $s=5.6$", transform=axs.transAxes, ha="left",
             va="center", fontsize=10, color="0.4")
    axs.set_yticks(y)
    axs.set_yticklabels([p.replace("Pat_", "") for p in C.COHORT], fontsize=8)
    axs.set_ylim(-0.6, len(C.COHORT) - 0.4)
    axs.set_xlabel(r"$\mathrm{task_{test}}\!\leftrightarrow\!\mathrm{rest_{post}}$ cophenetic "
                   r"concordance  $\rho_{\mathrm{sym}}$  (vs matched-strength "
                   r"$p_{95}$, $|$)", fontsize=12)
    axs.tick_params(axis="y", length=0)
    for s in ("top", "right"):
        axs.spines[s].set_visible(False)
    target.text(0.03, band.y0 - 0.10, r"$\mathbf{b}$", fontsize=17, va="top", ha="left",
                fontweight="bold")

    print(f"fig_reinstatement_learn ({PATIENT} {BAND}, mst@0.20 s={s0:.2f}, rho_sym):")
    print(f"  reachable ceiling rel={rel:.2f}  matched-strength null p95={null:.2f}")
    for ph in PHASES:
        print(f"  {ph:10s}: bar rho_sym={vals[ph]:+.2f}   per-leaf mean {np.nanmean(pl[ph]):+.2f}")
    print(f"  cohort beta reinstatement: {npos}/10 obs>0 @ s=5.6")


def main():
    fig = plt.figure(figsize=(15.0, 8.9))
    draw(fig)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
