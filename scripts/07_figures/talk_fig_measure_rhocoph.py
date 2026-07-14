#!/usr/bin/env python3
r"""talk_fig_measure_rhocoph -- the MEASURE slide (S2): what does rho^coph read?

Shows the measure by putting the colour ON the dendrogram: a reference hierarchy, then
two comparison hierarchies shaded by per-contact cophenetic preservation to it (green =
preserved, grey = reshuffled), each captioned with the whole-tree rho^coph. The dynamic
range is anchored by two neutral, pre-result examples -- NOT the trace:
  (1) the OTHER split-half of the SAME resting phase  -> high rho^coph (green): the
      hierarchy is reproducible within one recording;
  (2) a STRENGTH-MATCHED shuffle of the same phase     -> low  rho^coph (grey): destroy
      the wiring and the hierarchy scrambles (this also foreshadows the sole null).
No task, no cohort gate, no 4-fold taxonomy -- those are S3/S7. Style matches
fig_reinstatement_learn_mst020 (the figure the measure visual should look like).

Writes: data/outputs/figures/talk/fig_measure_rhocoph.png  (transparent, Canva-ready)
"""
from __future__ import annotations
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.stats import spearmanr

# reuse the mst@0.20 backbone/scale layer + the circular-dendrogram engine
sys.path.insert(0, str(Path(__file__).resolve().parents[1] /
                       "01_compute" / "figures_embedded" / "new_results_sec1"))
import _common as C  # noqa: E402
from lrg_eegfc.visuals.lrg import plot_circular_dendrogram  # noqa: E402
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle  # noqa: E402

C.use_lrg_style()

PAT, BAND = "Pat_08", "beta"
PRES_CM = LinearSegmentedColormap.from_list("pres", ["#d7d7d7", "#7cc47f", "#136b32"])
NORM = Normalize(0.0, 0.8)
CAP = "0.55"                    # caption grey (reads on light + mid backgrounds)
OUT = C.ROOT / "data" / "outputs" / "figures" / "talk" / "fig_measure_rhocoph.png"


def per_leaf_preservation(A, B):
    """For each leaf i: Spearman of its cophenetic-distance row (A vs B), excluding self."""
    n = A.shape[0]
    out = np.zeros(n)
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        out[i] = spearmanr(A[i, m], B[i, m]).correlation
    return out


def _triu(D):
    return D[np.triu_indices(D.shape[0], 1)]


def _sim_bar(ax, val, color, caption):
    ax.barh([0], [1.0], color="0.90", height=0.55, zorder=1)
    if val is not None:
        ax.barh([0], [max(val, 0.0)], color=color, height=0.55, zorder=2)
    ax.set_xlim(0, 1); ax.set_ylim(-1.6, 0.7)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(0.5, -0.95, caption, ha="center", va="top", fontsize=11.0, color=CAP)


def main():
    WA = C.load_phase(PAT, "A", BAND)          # rest_pre split-half A  (reference)
    WB = C.load_phase(PAT, "B", BAND)          # rest_pre split-half B  (same phase)
    N = WA.shape[0]
    rng = np.random.default_rng(2026)
    n_swaps = 20 * N * (N - 1) // 2
    WS = matched_strength_shuffle(WA, n_swaps, rng, 1.0)   # strength-matched shuffle

    Zref = C.tree_at_scale(WA); Dref = C.coph_square_at_scale(WA)
    ZB = C.tree_at_scale(WB);   DB = C.coph_square_at_scale(WB)
    ZS = C.tree_at_scale(WS);   DS = C.coph_square_at_scale(WS)

    rhoB = float(spearmanr(_triu(DB), _triu(Dref)).correlation)
    rhoS = float(spearmanr(_triu(DS), _triu(Dref)).correlation)
    colB = [PRES_CM(NORM(np.clip(np.nan_to_num(v), 0, 0.8))) for v in per_leaf_preservation(DB, Dref)]
    colS = [PRES_CM(NORM(np.clip(np.nan_to_num(v), 0, 0.8))) for v in per_leaf_preservation(DS, Dref)]
    colRef = ["#c9c9c9"] * N                    # reference: neutral

    panels = [
        (Zref, colRef, None, "0.85", "a hierarchy  (reference)"),
        (ZB,   colB,   rhoB, "#3a9a4f", rf"same rest, other half   $\rho^{{\mathrm{{coph}}}}={rhoB:.2f}$"),
        (ZS,   colS,   rhoS, "0.6",     rf"strength-matched shuffle   $\rho^{{\mathrm{{coph}}}}={rhoS:.2f}$"),
    ]

    fig = plt.figure(figsize=(12.6, 5.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[3.4, 0.5], hspace=0.10, wspace=0.06,
                          left=0.02, right=0.98, top=0.99, bottom=0.16)
    node_s = max(6.0, 900.0 / (N + 1))
    for j, (Z, colors, val, barcol, cap) in enumerate(panels):
        ax = fig.add_subplot(gs[0, j])
        plot_circular_dendrogram(ax, C.log_height_linkage(Z), leaf_colors=colors,
                                 leaf_labels=None, node_size=node_s, line_width=1.05,
                                 r_inner=0.10)
        axb = fig.add_subplot(gs[1, j])
        _sim_bar(axb, val, barcol, cap)

    # shared grey->green colourbar
    cax = fig.add_axes([0.34, 0.075, 0.32, 0.022])
    cb = fig.colorbar(ScalarMappable(norm=NORM, cmap=PRES_CM), cax=cax, orientation="horizontal")
    cb.set_label(r"per-contact preservation  $\rho^{\mathrm{coph}}$", fontsize=10, color=CAP)
    cb.set_ticks([0.0, 0.4, 0.8]); cb.ax.tick_params(labelsize=9, colors=CAP)
    cb.solids.set_rasterized(False)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"{PAT} {BAND}  N={N}   rho^coph(same-half)={rhoB:+.2f}   rho^coph(shuffle)={rhoS:+.2f}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
