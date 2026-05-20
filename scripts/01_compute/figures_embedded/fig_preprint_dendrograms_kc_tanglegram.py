"""KC λ=0 tanglegram TEST: 4 dendros + per-leaf polylines below.

For each leaf L, draw one polyline in figure coordinates threading from
its display x in rPre → tLearn → tTest → rPost. Polylines occupy a band
below the dendrograms.

Reading rule:
- Shallow (near-horizontal) segment between two phases = leaf stayed at
  the same display position = trees grouped the same way locally.
- Steep segment = leaf moved a lot = trees disagree on that leaf.
- Cluster of near-horizontal, non-crossing lines between two levels =
  two trees similar in that region.

Line color: leaf-resolved KC λ=0 trace strength
    s_leaf[L] = mean_j s(L, j),   s(i,j) defined as in the per-subtree
script. Red = trace-strong leaf, blue = anti-trace, white ≈ neutral.

Single-patient test only — filename prefixed TEST_ so it doesn't
collide with the batch PDFs.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.metrics.tree_distance import kc_vectors
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.lrg import load_lrg_result

use_lrg_style()


PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
PHASE_TITLES = {
    "rest_pre": r"$\mathrm{rPre}$",
    "task_learn": r"$\mathrm{tLearn}$",
    "task_test": r"$\mathrm{tTest}$",
    "rest_post": r"$\mathrm{rPost}$",
}
FC = "imcoh_abs"

ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "data/reports/preprint/figure1_beta_trace/dendrograms_kc_trace"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def leaves_under_each_node(Z: np.ndarray) -> dict[int, frozenset[int]]:
    n = Z.shape[0] + 1
    out: dict[int, frozenset[int]] = {i: frozenset({i}) for i in range(n)}
    for i in range(n - 1):
        a = int(Z[i, 0])
        b = int(Z[i, 1])
        out[n + i] = out[a] | out[b]
    return out


def m_vec_to_pair_matrix(m_vec: np.ndarray, n: int) -> np.ndarray:
    M = np.zeros((n, n), dtype=float)
    M[np.triu_indices(n, k=1)] = m_vec
    return M + M.T


def per_subtree_score(
    leaves_dict: dict[int, frozenset[int]],
    s_pair: np.ndarray,
) -> dict[int, float]:
    out: dict[int, float] = {}
    for nid, leaves in leaves_dict.items():
        if len(leaves) < 2:
            continue
        arr = np.fromiter(leaves, dtype=int)
        sub = s_pair[np.ix_(arr, arr)]
        out[nid] = float(sub[np.triu_indices(len(arr), k=1)].mean())
    return out


def main(patient: str, band: str) -> None:
    out_pdf = OUT_DIR / f"TEST_kc_tanglegram_{patient}_{band}.pdf"
    Zs = {p: load_lrg_result(patient, p, band, FC).linkage_matrix for p in PHASES}
    n = Zs[PHASES[0]].shape[0] + 1

    m_vec = {p: kc_vectors(Z)[0] for p, Z in Zs.items()}
    M = {p: m_vec_to_pair_matrix(m, n) for p, m in m_vec.items()}

    s = 0.5 * (np.abs(M["rest_pre"] - M["task_test"]) +
               np.abs(M["rest_pre"] - M["rest_post"])) \
        - np.abs(M["task_test"] - M["rest_post"])

    s_off = s.copy()
    np.fill_diagonal(s_off, np.nan)
    s_leaf = np.nanmean(s_off, axis=1)

    leaves_per_phase = {p: leaves_under_each_node(Z) for p, Z in Zs.items()}
    scores_per_phase = {
        p: per_subtree_score(leaves_per_phase[p], s) for p in PHASES
    }

    all_scores = np.fromiter(
        (v for d in scores_per_phase.values() for v in d.values()),
        dtype=float,
    )
    vmax_sub = float(np.max(np.abs(all_scores)))
    cmap = plt.get_cmap("RdBu_r")
    norm_sub = Normalize(vmin=-vmax_sub, vmax=vmax_sub)

    vmax_leaf = float(np.max(np.abs(s_leaf)))
    norm_leaf = Normalize(vmin=-vmax_leaf, vmax=vmax_leaf)
    leaf_color = [mcolors.to_hex(cmap(norm_leaf(s_leaf[L]))) for L in range(n)]

    def color_for_phase(phase: str):
        scores = scores_per_phase[phase]

        def cf(node_id: int) -> str:
            sc = scores.get(node_id)
            if sc is None or not np.isfinite(sc):
                return "#888888"
            return mcolors.to_hex(cmap(norm_sub(sc)))

        return cf

    tmin = min(Z[0, 2] for Z in Zs.values()) * 0.8
    tmax = max(Z[-1, 2] for Z in Zs.values()) * 1.05

    fig = plt.figure(figsize=(13.5, 5.5))
    gs = fig.add_gridspec(
        2, 4, height_ratios=[8, 2.5],
        hspace=0.05, wspace=0.10,
        left=0.05, right=0.97, top=0.84, bottom=0.05,
    )
    axes_top = [fig.add_subplot(gs[0, 0])]
    for j in range(1, 4):
        axes_top.append(fig.add_subplot(gs[0, j], sharey=axes_top[0]))
    ax_tangle = fig.add_subplot(gs[1, :])

    d_per_phase = {}
    for col, phase in enumerate(PHASES):
        ax = axes_top[col]
        ax.set_title(PHASE_TITLES[phase])
        Z = Zs[phase]
        d = dendrogram(
            Z, ax=ax,
            no_labels=True,
            link_color_func=color_for_phase(phase),
            above_threshold_color="#888888",
        )
        d_per_phase[phase] = d
        ax.set_ylim(tmin, tmax)
        ax.set_xticks([])
        ax.set_xlabel("")
        if col == 0:
            ax.set_ylabel(r"merge height  $\hat{D}(\tau)$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    ax_tangle.set_axis_off()
    ax_tangle.set_xlim(0, 1)
    ax_tangle.set_ylim(0, 1)

    fig.canvas.draw()

    leaf_fig_x: dict[str, dict[int, float]] = {}
    for phase, ax in zip(PHASES, axes_top):
        d = d_per_phase[phase]
        pos_of_leaf = {L: pos for pos, L in enumerate(d["leaves"])}
        leaf_fig_x[phase] = {}
        for L in range(n):
            data_x = 5 + 10 * pos_of_leaf[L]
            px, _ = ax.transData.transform((data_x, 0))
            fx, _ = fig.transFigure.inverted().transform((px, 0))
            leaf_fig_x[phase][L] = float(fx)

    bbox_tangle = ax_tangle.get_position()
    y_top = bbox_tangle.y1
    y_bot = bbox_tangle.y0
    y_levels = np.linspace(y_top, y_bot, len(PHASES))

    for L in range(n):
        xs = [leaf_fig_x[p][L] for p in PHASES]
        ys = list(y_levels)
        line = Line2D(
            xs, ys,
            color=leaf_color[L], lw=0.5, alpha=0.55,
            transform=fig.transFigure, zorder=2,
        )
        fig.add_artist(line)

    for y, phase in zip(y_levels, PHASES):
        fig.text(0.975, y, PHASE_TITLES[phase],
                 fontsize=6, va="center", ha="left", color="#444",
                 transform=fig.transFigure)

    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    fig.text(
        0.05, 0.95,
        rf"{patient}  |  {band_tex}  |  KC$_{{\lambda=0}}$ tanglegram TEST  |  "
        rf"$N = {n}$  |  $s_{{\mathrm{{leaf}}}}$ $\pm{vmax_leaf:.2f}$  |  "
        rf"$\bar{{s}}_{{\mathrm{{sub}}}}$ $\pm{vmax_sub:.2f}$",
        fontweight="bold",
    )
    fig.text(
        0.05, 0.92,
        "U-shapes: subtree-mean KC trace (RdBu_r).  Below: one polyline per "
        "leaf threading from rPre to rPost; line color = leaf-resolved KC "
        "trace strength. Crossings between two adjacent levels = trees "
        "reorder leaves locally; near-horizontal bundle = trees agree.",
        fontsize=7, color="#444",
    )

    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"DONE: {out_pdf}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--patient", default="Pat_06")
    parser.add_argument("--band", default="beta")
    args = parser.parse_args()
    main(args.patient, args.band)
