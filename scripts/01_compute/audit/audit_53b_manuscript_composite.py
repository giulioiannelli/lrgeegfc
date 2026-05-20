#!/usr/bin/env python3
"""Audit 53b — manuscript composite figure for KC leaf topology trace.

Single 4-row x 3-column PDF for §5 of the manuscript:
  Row 1: Pat_07 β  focal leaf 62
  Row 2: Pat_05 β  focal leaf 66
  Row 3: Pat_10 β  focal leaf 37
  Row 4: Pat_14 γ_h (sentinel; predicate empty)

No f_i histograms; column headers on top row only; per-row label on the
left margin. Reuses helpers from
``audit_53_kc_leaf_topological_memory``.

Output:
    data/audit/kc_leaf_topological_memory/figures/manuscript_kc_leaf_topology.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgba
import numpy as np
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.metrics.tree_distance import (
    _build_parent_height_depth,
    kc_vectors,
)
from lrg_eegfc.workflow.lrg import load_lrg_result

sys.path.insert(0, str(Path(__file__).parent))
from audit_53_kc_leaf_topological_memory import (  # noqa: E402
    GRAY_LINK,
    PHASES,
    PHASE_LABELS,
    compute_leaf_x_positions,
    expand_per_leaf_score,
    find_mrca,
    mrca_horizontal,
    pair_index,
    pair_mrca_subtree_sizes,
    path_segments_to_mrca,
)

OUT_DIR = ROOT / "data/audit/kc_leaf_topological_memory"

CELLS = [
    ("Pat_07", "beta",       62),
    ("Pat_05", "beta",       66),
    ("Pat_10", "beta",       37),
    ("Pat_14", "high_gamma", None),  # sentinel
]

DELTA         = 1
TOP_M         = 8
MIN_M         = 2
MAX_MRCA_FRAC = 0.50

ACCENT       = "#d62728"
ACCENT_PART  = "#fa8a87"


def compute_cell(
    patient: str,
    band: str,
    focal_override: int | None,
) -> dict:
    Z_by_phase = {}
    for ph in PHASES:
        res = load_lrg_result(patient, ph, band, fc_method="imcoh_abs")
        if res is None or res.linkage_matrix is None:
            raise FileNotFoundError(f"missing LRG result {patient} {ph} {band}")
        Z_by_phase[ph] = np.asarray(res.linkage_matrix)
    n = int(Z_by_phase["rest_pre"].shape[0]) + 1

    m_pre  = kc_vectors(Z_by_phase["rest_pre"])[0]
    m_tt   = kc_vectors(Z_by_phase["task_test"])[0]
    m_post = kc_vectors(Z_by_phase["rest_post"])[0]
    size_tt   = pair_mrca_subtree_sizes(Z_by_phase["task_test"])
    size_post = pair_mrca_subtree_sizes(Z_by_phase["rest_post"])

    s_max = int(np.ceil(MAX_MRCA_FRAC * n))
    matched_tt_post = np.abs(m_tt - m_post) <= DELTA
    differs_pre_tt  = np.abs(m_pre - m_tt) > DELTA
    nontrivial_m    = (m_tt >= MIN_M) & (m_post >= MIN_M)
    fine_grained    = (size_tt <= s_max) & (size_post <= s_max)
    T = matched_tt_post & differs_pre_tt & nontrivial_m & fine_grained

    f = expand_per_leaf_score(T, n)

    if focal_override is not None:
        focal = int(focal_override)
    else:
        focal = int(np.argmax(f))
    partners = []
    for j in range(n):
        if j == focal:
            continue
        k = pair_index(n, focal, j)
        if T[k]:
            partners.append(j)
    partners.sort(key=lambda j: -int(m_tt[pair_index(n, focal, j)]))
    top_M = partners[:TOP_M]

    coord_by_phase = {}
    for ph in PHASES:
        Z = Z_by_phase[ph]
        parent, height, _ = _build_parent_height_depth(Z)
        x_pos, leaf_order = compute_leaf_x_positions(Z)
        coord_by_phase[ph] = {
            "Z": Z, "parent": parent, "height": height,
            "x_pos": x_pos, "leaf_order": leaf_order,
        }

    return {
        "patient": patient,
        "band": band,
        "n": n,
        "focal": focal,
        "top_M": top_M,
        "f_focal": float(f[focal]),
        "n_partners": len(partners),
        "coord_by_phase": coord_by_phase,
    }


def _draw_dendrogram_with_bundle(ax, info: dict, ph: str) -> tuple[float, float]:
    coord  = info["coord_by_phase"][ph]
    Z      = coord["Z"]
    parent = coord["parent"]
    height = coord["height"]
    x_pos  = coord["x_pos"]
    n      = info["n"]
    focal  = info["focal"]
    top_M  = info["top_M"]

    dendrogram(
        Z, ax=ax, no_labels=True,
        color_threshold=0,
        above_threshold_color=GRAY_LINK,
        link_color_func=lambda _id: GRAY_LINK,
    )
    merge_heights = sorted(set(Z[:, 2]))
    tmin = 0.0
    tmax = merge_heights[-1] * 1.05
    ax.set_ylim(tmin, tmax)
    ax.set_xticks([])
    ax.tick_params(axis="y", labelsize=7)
    for sn in ("top", "right", "bottom"):
        ax.spines[sn].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)

    if top_M:
        rgba_path = to_rgba(ACCENT, 0.55)
        rgba_mrca = to_rgba(ACCENT, 0.80)
        segs: list = []
        cols: list = []
        for j in top_M:
            mrca = find_mrca(parent, focal, j)
            if mrca < n:
                continue
            for s in path_segments_to_mrca(focal, mrca, parent, height, x_pos):
                segs.append(s); cols.append(rgba_path)
            for s in path_segments_to_mrca(j, mrca, parent, height, x_pos):
                segs.append(s); cols.append(rgba_path)
            segs.append(mrca_horizontal(Z, mrca, n, x_pos, height))
            cols.append(rgba_mrca)
        if segs:
            ax.add_collection(
                LineCollection(segs, colors=cols, linewidths=1.6, zorder=4)
            )

        for lf in range(n):
            if lf == focal:
                ax.vlines(x_pos[lf], tmin, height[int(parent[lf])],
                          color=ACCENT, lw=1.9, zorder=5)
            elif lf in top_M:
                ax.vlines(x_pos[lf], tmin, height[int(parent[lf])],
                          color=ACCENT_PART, lw=1.4, zorder=4)

        ax.scatter([x_pos[focal]], [tmax * 0.985], marker="v",
                   color=ACCENT, s=46, zorder=6,
                   edgecolors="white", linewidths=0.6)

    return tmin, tmax


def main() -> int:
    cells_info = []
    for patient, band, focal in CELLS:
        info = compute_cell(patient, band, focal)
        cells_info.append(info)
        print(
            f"  {patient} {band:11s}  n={info['n']:3d}  "
            f"focal={info['focal']:3d}  f*={info['f_focal']:.3f}  "
            f"|P|={info['n_partners']:3d}  drawn={len(info['top_M'])}"
        )

    n_rows = len(cells_info)
    fig = plt.figure(figsize=(11.5, 2.5 * n_rows + 0.6))
    left, right = 0.135, 0.99
    top, bottom = 0.945, 0.045
    gs = fig.add_gridspec(
        n_rows, 3,
        left=left, right=right, top=top, bottom=bottom,
        hspace=0.28, wspace=0.06,
    )

    first_axes = []
    for r, info in enumerate(cells_info):
        first_ax = None
        for c, ph in enumerate(PHASES):
            sharey = first_ax
            ax = fig.add_subplot(gs[r, c], sharey=sharey)
            if first_ax is None:
                first_ax = ax
            _draw_dendrogram_with_bundle(ax, info, ph)
            if r == 0:
                ax.set_title(PHASE_LABELS[ph], fontsize=11, pad=4)
            if c == 0:
                ax.set_ylabel("merge height", fontsize=8)
            else:
                ax.tick_params(axis="y", labelleft=False)
                ax.spines["left"].set_visible(False)
        first_axes.append(first_ax)

    fig.canvas.draw()
    for r, (info, first_ax) in enumerate(zip(cells_info, first_axes)):
        bbox = first_ax.get_position()
        band_tex = BRAIN_BAND_TEX_DICT.get(info["band"], info["band"])
        line1 = rf"{info['patient']}  {band_tex}"
        line2 = rf"$f^\star={info['f_focal']:.2f}$"
        line3 = rf"$|P|={info['n_partners']}$"
        y_center = bbox.y0 + bbox.height / 2.0
        fig.text(0.012, y_center + 0.022, line1,
                 ha="left", va="center", fontsize=11, fontweight="bold")
        fig.text(0.012, y_center,         line2,
                 ha="left", va="center", fontsize=9.5)
        fig.text(0.012, y_center - 0.022, line3,
                 ha="left", va="center", fontsize=9.5)


    out_path = OUT_DIR / "figures" / "manuscript_kc_leaf_topology.pdf"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    print(f"\nwrote {out_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
