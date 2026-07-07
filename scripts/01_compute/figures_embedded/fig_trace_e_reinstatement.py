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
10/10 patients (bottom), scale-robust across ten diffusion tau.

Why not colour by a fixed partition: rho^coph is a pairwise cophenetic-distance correlation,
not clade identity (task and post share only 3 exact clades despite rho=0.58) — so the
faithful leaf-level colour is per-leaf preservation, not an arbitrary module palette.

Reads : LRG cache (Pat_08 beta imcoh_abs, three phases)
        data/audit/replay_states/sustained_reinstatement_per_patient.csv
Writes: data/reports/results_section1/fig_trace_e_reinstatement.pdf
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
from scipy.stats import spearmanr

from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.lrg import plot_circular_dendrogram
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.lrg import load_lrg_result

ROOT = setup_script_env()
use_lrg_style()

PP = ROOT / "data/audit/replay_states/sustained_reinstatement_per_patient.csv"
OUT = ROOT / "data/reports/results_section1/fig_trace_e_reinstatement.pdf"

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


def main():
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

    fig = plt.figure(figsize=(12.8, 7.9))
    gs = fig.add_gridspec(3, 3, height_ratios=[3.05, 0.5, 1.1], hspace=0.30,
                          wspace=0.06, left=0.03, right=0.93, top=0.95, bottom=0.10)

    for j, ph in enumerate(PHASES):
        ax = fig.add_subplot(gs[0, j])
        plot_circular_dendrogram(ax, Z[ph], leaf_colors=leaf_colors[ph],
                                 leaf_labels=leaf_labels, node_size=12.0,
                                 line_width=1.05, r_inner=0.10,
                                 label_fontsize=3.0, label_offset=1.03)
        ax.set_title(TITLES[ph], fontsize=17, pad=6)
        axb = fig.add_subplot(gs[1, j])
        if ph == "task_test":
            _sim_bar(axb, 1.0, "0.78", "(reference)")
        else:
            col = C_POST if ph == "rest_post" else C_PRE
            _sim_bar(axb, sim[ph], col,
                     r"whole-tree $\rho^{\mathrm{coph}}$ to $\mathrm{task_{test}}$"
                     f" = {sim[ph]:.2f}")

    fig.text(0.03, 0.955, r"$\mathbf{e}$", fontsize=17, va="bottom", ha="left",
             fontweight="bold")

    # colourbar defining the leaf colour = per-contact preservation to task
    cax = fig.add_axes([0.945, 0.52, 0.013, 0.34])
    cb = fig.colorbar(ScalarMappable(norm=NORM, cmap=PRES_CM), cax=cax)
    cb.set_label(r"per-contact $\rho^{\mathrm{coph}}$ to $\mathrm{task_{test}}$", fontsize=10)
    cb.set_ticks([0.0, 0.4, 0.8])
    cb.ax.tick_params(labelsize=9)

    # -- bottom: cohort window-level reinstatement (post vs pre proximity to task) --
    axs = fig.add_subplot(gs[2, :])
    pp = pd.read_csv(PP)
    v = pp["A1_level_post_minus_pre"].to_numpy(float)
    z, p = wilcoxon_z(v)
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
    axs.text(0.015, 0.82, "window-level; scale-robust across ten $\\tau$",
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

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight")
    plt.close(fig)

    print(f"fig:trace_e — reinstatement ({PATIENT} {BAND}), per-leaf preservation colour\n")
    for ph in PHASES:
        print(f"  per-leaf pres to task, {ph:10s}: mean {np.nanmean(pl[ph]):.2f} "
              f"frac>0.5 {np.mean(pl[ph] > 0.5):.2f}")
    print(f"  whole-tree rho^coph: pre {sim['rest_pre']:.3f} -> post {sim['rest_post']:.3f}")
    print(f"  cohort window: {npos}/10 post>pre, z={z:.2f} p={p:.4f}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
