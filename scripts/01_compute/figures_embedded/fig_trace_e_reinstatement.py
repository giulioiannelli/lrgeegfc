#!/usr/bin/env python3
r"""fig:trace_e — the resting hierarchy becomes task-like after the task (Results §1, para e).

CORE MESSAGE: the post-task resting cophenetic hierarchy re-expresses the task hierarchy;
the pre-task one does not. The colour IS the trace: every contact is shaded by its per-leaf
cophenetic preservation to task_test — how well it keeps its cophenetic distances to all
other contacts, relative to the task phase (the per-contact term of rho^coph). So rest_pre
renders grey (per-leaf preservation mean 0.20, 17% of contacts > 0.5) and rest_post renders
green (mean 0.50, 77% > 0.5); 95/120 contacts are greener in post than pre. The
size-weighted branches make a preserved sub-tree glow green. Exemplar Pat_08 (beta,
strongest tracer): whole-tree rho^coph to task rises 0.05 (pre) -> 0.58 (post), shown by
the bars. Cohort: a window-level measure agrees, post-task rest closer to task than pre for
10/10 patients (bottom). Read at the single working scale tau_min = 1/lambda_max.

Why not colour by a fixed partition: rho^coph is a pairwise cophenetic-distance correlation,
not clade identity (task and post share only 3 exact clades despite rho=0.58) — so the
faithful leaf-level colour is per-leaf preservation, not an arbitrary module palette.

Reads : LRG cache (Pat_08 beta imcoh_abs, three phases)
        data/audit/replay_states/sustained_reinstatement_per_patient.csv
Writes: data/preprint/figures/results_section1/fig_trace_e_reinstatement.pdf
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.lrg import plot_circular_dendrogram
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.lrg import load_lrg_result

ROOT = setup_script_env()
use_lrg_style()

PP = ROOT / "data/audit/replay_states/sustained_reinstatement_per_patient.csv"
OUT = ROOT / "data/preprint/figures/results_section1/fig_trace_e_reinstatement.pdf"

PATIENT, BAND, FC = "Pat_08", "beta", "imcoh_abs"
PHASES = ["rest_pre", "task_test", "rest_post"]
TITLES = {"rest_pre": r"$\mathrm{rest_{pre}}$", "task_test": r"$\mathrm{task_{test}}$",
          "rest_post": r"$\mathrm{rest_{post}}$"}

# per-contact preservation colormap: grey (not preserved) -> green (preserved)
PRES_CM = LinearSegmentedColormap.from_list("pres", ["#d7d7d7", "#7cc47f", "#136b32"])
NORM = Normalize(0.0, 0.8)
C_POST, C_PRE, C_ANTI = "#3a9a4f", "#9a9a9a", "#d1352b"


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
    ax.text(0.5, -0.85, caption, ha="center", va="top", fontsize=10.5, color="0.25")


def _touch_node_size(target, cell_spec, n_leaves, r_outer=1.0, frac=0.86):
    """Scatter s (pts^2) so equiangular leaf markers just touch (frac<1 leaves a hair gap).
    Sized from the cell's PHYSICAL inches, so it holds standalone AND in the compound tile
    (a fixed pts^2 would overlap once the tile shrinks the circle)."""
    dpi = getattr(target, "dpi", None) or target.get_figure().dpi
    pos = cell_spec.get_position(target)                    # cell as fraction of the tile
    used_in = min(pos.width * target.bbox.width,
                  pos.height * target.bbox.height) / dpi     # equal-aspect circle fits smaller side
    ring_r_in = (used_in / 2.26) * r_outer                  # data half-range 1.13 -> span 2.26
    arc_in = 2.0 * np.pi * ring_r_in / max(n_leaves, 1)     # arc between adjacent leaves
    return np.pi / 4.0 * (frac * arc_in * 72.0) ** 2


def draw(target):
    res = {ph: load_lrg_result(PATIENT, ph, BAND, FC) for ph in PHASES}
    Z = {ph: res[ph].linkage_matrix for ph in PHASES}
    D = {ph: squareform(res[ph].ultrametric_matrix) for ph in PHASES}

    # colour = per-contact cophenetic preservation to task (the leaf-level term of rho^coph)
    pl = {ph: per_leaf_preservation(D[ph], D["task_test"]) for ph in PHASES}
    leaf_colors = {ph: [PRES_CM(NORM(np.clip(np.nan_to_num(v, nan=0.0), 0, 0.8)))
                        for v in pl[ph]] for ph in PHASES}

    rdf = load_channel_regions(PATIENT)
    rdf = rdf.loc[:, ~rdf.columns.duplicated()]
    leaf_labels = rdf["label"].astype(str).tolist()

    sim = {ph: float(spearmanr(res[ph].ultrametric_matrix,
                               res["task_test"].ultrametric_matrix).correlation)
           for ph in ("rest_pre", "rest_post")}

    # more vertical stack: dendrograms / whole-tree bars / horizontal colourbar / window plot.
    # right runs to 0.965 (the old right-side vertical colourbar is gone) so the trees grow.
    gs = target.add_gridspec(4, 3, height_ratios=[3.5, 0.42, 0.55, 1.15], hspace=0.26,
                             wspace=0.05, left=0.03, right=0.965, top=0.975, bottom=0.075)
    node_s = _touch_node_size(target, gs[0, 0], int(Z["rest_pre"].shape[0]) + 1)

    for j, ph in enumerate(PHASES):
        ax = target.add_subplot(gs[0, j])
        plot_circular_dendrogram(ax, Z[ph], leaf_colors=leaf_colors[ph],
                                 leaf_labels=leaf_labels, node_size=node_s,
                                 line_width=1.05, r_inner=0.10,
                                 label_fontsize=3.0, label_offset=1.03)
        ax.set_title(TITLES[ph], fontsize=17, pad=6)
        axb = target.add_subplot(gs[1, j])
        if ph == "task_test":
            _sim_bar(axb, 1.0, "0.78", "(reference)")
        else:
            col = C_POST if ph == "rest_post" else C_PRE
            _sim_bar(axb, sim[ph], col,
                     r"whole-tree $\rho^{\mathrm{coph}}$ to $\mathrm{task_{test}}$"
                     f" = {sim[ph]:.2f}")

    target.text(0.03, 0.955, r"$\mathbf{d}$", fontsize=17, va="bottom", ha="left",
                fontweight="bold")

    # horizontal colourbar (leaf colour = per-contact preservation to task), centred in its
    # own band below the whole-tree bars and above the window-level plot.
    band = gs[2, :].get_position(target)
    cb_w, cb_h = 0.40, 0.024
    cax = target.add_axes([0.5 - cb_w / 2, band.y0 + 0.55 * band.height, cb_w, cb_h])
    cb = target.colorbar(ScalarMappable(norm=NORM, cmap=PRES_CM), cax=cax,
                         orientation="horizontal")
    cb.set_label(r"per-contact $\rho^{\mathrm{coph}}$ to $\mathrm{task_{test}}$", fontsize=10)
    cb.set_ticks([0.0, 0.4, 0.8])
    cb.ax.tick_params(labelsize=9)
    cax.xaxis.set_ticks_position("bottom")
    cax.xaxis.set_label_position("bottom")

    # -- bottom: cohort window-level reinstatement (post vs pre proximity to task) --
    axs = target.add_subplot(gs[3, :])
    pp = pd.read_csv(PP)
    v = pp["A1_level_post_minus_pre"].to_numpy(float)
    z, _ = wilcoxon_z(v)  # normal-approx z, diagnostic only
    # Exact one-sided signed-rank p (n=10): the normal approximation overstates
    # p at this sample size (0.0025->"0.003"); the manuscript reports the exact
    # test (1/2^10 = 0.001, LOO 0.002). Keep figure and prose on the same test.
    p = float(wilcoxon(v, alternative="greater", method="exact").pvalue)
    rng = np.random.default_rng(0)
    yj = 0.11 * rng.standard_normal(len(v))
    axs.axvline(0.0, color="0.55", lw=1.1, ls="--", zorder=1)
    axs.scatter(v, yj, s=135, color=[C_POST if x > 0 else C_ANTI for x in v],
                edgecolor="white", linewidth=0.9, zorder=4)
    med = float(np.median(v))
    axs.scatter([med], [0], marker="D", s=150, color="#1f7a34", edgecolor="white",
                linewidth=1.2, zorder=6)
    axs.annotate("cohort median", xy=(med, 0), xytext=(med, 0.42), fontsize=9.5,
                 color="#1f7a34", ha="center", va="bottom")
    npos = int((v > 0).sum())
    axs.text(0.985, 0.82, f"{npos}/10 patients closer to task   $p={p:.3f}$",
             transform=axs.transAxes, ha="right", va="center", fontsize=12.5,
             color=C_POST, fontweight="bold")
    axs.text(0.015, 0.82, "window-level",
             transform=axs.transAxes, ha="left", va="center", fontsize=10, color="0.4")
    axs.set_xlim(-0.06, 0.20)
    axs.set_ylim(-0.6, 0.6)
    axs.set_yticks([])
    axs.set_xlabel(r"$\mathrm{rest_{post}}$ closer to $\mathrm{task_{test}}$ than "
                   r"$\mathrm{rest_{pre}}$   "
                   r"(window-level, $A_1^{\mathrm{post}}\!-\!A_1^{\mathrm{pre}}$)",
                   fontsize=12.5)
    axs.tick_params(axis="x", labelsize=11)
    for s in ("top", "left", "right"):
        axs.spines[s].set_visible(False)

    print(f"fig:trace_e — reinstatement ({PATIENT} {BAND}), per-leaf preservation colour\n")
    for ph in PHASES:
        print(f"  per-leaf pres to task, {ph:10s}: mean {np.nanmean(pl[ph]):.2f} "
              f"frac>0.5 {np.mean(pl[ph] > 0.5):.2f}")
    print(f"  whole-tree rho^coph: pre {sim['rest_pre']:.3f} -> post {sim['rest_post']:.3f}")
    print(f"  cohort window: {npos}/10 post>pre, z={z:.2f} p={p:.4f}")


def main():
    fig = plt.figure(figsize=(12.2, 8.7))
    draw(fig)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
