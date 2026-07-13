#!/usr/bin/env python3
r"""fig_reinstatement_learn_mst020 -- the resting hierarchy becomes task-like AFTER the task,
across the full rest -> learn -> test -> rest arc (mst@0.20).

HEAD: the post-task resting cophenetic hierarchy re-expresses the task hierarchy; the pre-task
one does not -- rebuilt on the **mst@0.20** backbone at the reporting scale **s=5.6**, and now
resolved across the FULL paradigm: rest_pre -> task_learn -> task_test -> rest_post. The colour
IS the trace: every contact is shaded by its per-leaf cophenetic preservation to task_test (how
well it keeps its cophenetic distances to all other contacts relative to the test phase -- the
per-contact term of rho^coph). rest_pre renders pale (it has not seen the task); task_learn
already warms (the hierarchy is forming); task_test is the reference; and rest_post stays warm
-- the offline resting hierarchy holds the task geometry it did not have before. The whole-tree
rho^coph-to-test bars climb pre -> learn -> (test) -> post.

Exemplar Pat_08 (beta, a strong tracer: whole-phase beta rho_sym = +0.48 at s=5.6). The cohort
claim (bottom) is the load-bearing one: per-patient beta rho_sym^coph (task_test <-> rest_post
cophenetic concordance) at s=5.6 against the matched-strength null -- the same gate that clears
16/16 scales at p=0.001. Read on the mst@0.20 multiscale backbone, not the degenerate dense graph.

Why not colour by a fixed partition: rho^coph is a pairwise cophenetic-distance correlation, not
clade identity, so the faithful leaf-level colour is per-leaf preservation, not a module palette.

Reads : imcoh_abs FC (Pat_08, four phases, via load_fc_matrix + mst@0.20 backbone);
        data/sparsified_arc/ms_mst020/per_patient_scale.csv (cohort beta rho_sym @ s=5.6).
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

C.use_lrg_style()

PATIENT, BAND, FC = "Pat_08", "beta", "imcoh_abs"
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
TITLES = {"rest_pre": r"$\mathrm{rest_{pre}}$", "task_learn": r"$\mathrm{task_{learn}}$",
          "task_test": r"$\mathrm{task_{test}}$", "rest_post": r"$\mathrm{rest_{post}}$"}
PRES_CM = LinearSegmentedColormap.from_list("pres", ["#d7d7d7", "#7cc47f", "#136b32"])
NORM = Normalize(0.0, 0.8)
C_POST, C_PRE, C_ANTI = "#3a9a4f", "#9a9a9a", "#d1352b"
OUT = C.FIGDIR / "fig_reinstatement_learn_mst020.pdf"


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
    ax.barh([0], [1.0], color="0.90", height=0.55, zorder=1)
    if val is not None:
        ax.barh([0], [max(val, 0.0)], color=color, height=0.55, zorder=2)
    ax.set_xlim(0, 1)
    ax.set_ylim(-1.5, 0.7)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(0.5, -0.9, caption, ha="center", va="top", fontsize=10.0, color="0.25")


def draw(target):
    Z = {ph: C.tree_at_scale(phase_fc(ph)) for ph in PHASES}
    D = {ph: C.coph_square_at_scale(phase_fc(ph)) for ph in PHASES}
    pl = {ph: per_leaf_preservation(D[ph], D["task_test"]) for ph in PHASES}
    leaf_colors = {ph: [PRES_CM(NORM(np.clip(np.nan_to_num(v, nan=0.0), 0, 0.8)))
                        for v in pl[ph]] for ph in PHASES}
    sim = {ph: float(spearmanr(D[ph][np.triu_indices(D[ph].shape[0], 1)],
                               D["task_test"][np.triu_indices(D[ph].shape[0], 1)]).correlation)
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
        # + preservation colours are unchanged (see _common.log_height_linkage).
        plot_circular_dendrogram(ax, C.log_height_linkage(Z[ph]), leaf_colors=leaf_colors[ph],
                                 leaf_labels=leaf_labels, node_size=node_s,
                                 line_width=1.05, r_inner=0.10,
                                 label_fontsize=3.0, label_offset=1.03)
        ax.set_title(TITLES[ph], fontsize=16, pad=6)
        axb = target.add_subplot(gs[1, j])
        if ph == "task_test":
            _sim_bar(axb, 1.0, "0.85", "(reference)")
        else:
            col = C_POST if ph == "rest_post" else (C_PRE if ph == "rest_pre" else "#6f8fb0")
            _sim_bar(axb, sim[ph], col,
                     rf"$\rho^{{\mathrm{{coph}}}}$ to test = {sim[ph]:.2f}")

    target.text(0.03, 0.955, r"$\mathbf{a}$", fontsize=17, va="bottom", ha="left",
                fontweight="bold")

    # horizontal colourbar (leaf colour = per-contact preservation to task_test)
    band = gs[2, :].get_position(target)
    cb_w, cb_h = 0.40, 0.024
    cax = target.add_axes([0.5 - cb_w / 2, band.y0 + 0.5 * band.height, cb_w, cb_h])
    cb = target.colorbar(ScalarMappable(norm=NORM, cmap=PRES_CM), cax=cax,
                         orientation="horizontal")
    cb.set_label(r"per-contact $\rho^{\mathrm{coph}}$ to $\mathrm{task_{test}}$", fontsize=10)
    cb.set_ticks([0.0, 0.4, 0.8])
    cb.ax.tick_params(labelsize=9)
    cb.solids.set_rasterized(False)

    # -- bottom: cohort beta rho_sym (task<->post) at s=5.6 vs matched-strength null --
    axs = target.add_subplot(gs[3, :])
    pp = pd.read_csv(C.MS / "per_patient_scale.csv")
    s0 = C.SGRID[C.I_REPORT]
    x = pp[(pp.band == BAND) & np.isclose(pp.s, s0)].set_index("patient")
    obs = np.array([float(x.loc[p, "obs_rho"]) for p in C.COHORT])
    p95 = np.array([float(x.loc[p, "surr_p95"]) for p in C.COHORT])
    clears = np.array([float(x.loc[p, "p"]) < 0.05 for p in C.COHORT])
    y = np.arange(len(C.COHORT))
    axs.axvline(0.0, color="0.55", lw=1.0, ls="--", zorder=1)
    for yi, o, t, cl in zip(y, obs, p95, clears):
        axs.plot([0, o], [yi, yi], color=C_POST if o >= 0 else C_ANTI, lw=1.3, zorder=2)
        axs.plot([t], [yi], "|", color="0.5", ms=12, mew=1.8, zorder=3)   # own null p95
        axs.scatter(o, yi, s=70, color=C_POST if o >= 0 else C_ANTI,
                    edgecolor="black" if cl else "white", lw=1.2 if cl else 0.6, zorder=4)
    npos = int((obs > 0).sum())
    axs.text(0.985, 0.90, rf"$\beta$: {npos}/10 patients reinstate  ·  cohort gate "
             r"$p=0.001$ (16/16 scales)", transform=axs.transAxes, ha="right", va="center",
             fontsize=12, color=C_POST, fontweight="bold")
    axs.text(0.015, 0.90, r"cohort, $s=5.6$", transform=axs.transAxes, ha="left",
             va="center", fontsize=10, color="0.4")
    axs.set_yticks(y)
    axs.set_yticklabels([p.replace("Pat_", "") for p in C.COHORT], fontsize=8)
    axs.set_ylim(-0.6, len(C.COHORT) - 0.4)
    axs.set_xlabel(r"$\mathrm{task_{test}}\!\leftrightarrow\!\mathrm{rest_{post}}$ cophenetic "
                   r"concordance  $\rho_{\mathrm{sym}}^{\mathrm{coph}}$  (vs matched-strength "
                   r"$p_{95}$, $|$)", fontsize=12)
    axs.tick_params(axis="y", length=0)
    for s in ("top", "right"):
        axs.spines[s].set_visible(False)
    target.text(0.03, band.y0 - 0.10, r"$\mathbf{b}$", fontsize=17, va="top", ha="left",
                fontweight="bold")

    print(f"fig_reinstatement_learn ({PATIENT} {BAND}, mst@0.20 s={s0:.2f}):")
    for ph in PHASES:
        print(f"  per-leaf pres->test {ph:10s}: mean {np.nanmean(pl[ph]):.2f}  "
              f"frac>0.5 {np.mean(pl[ph] > 0.5):.2f}   whole-tree rho^coph {sim[ph]:+.2f}")
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
