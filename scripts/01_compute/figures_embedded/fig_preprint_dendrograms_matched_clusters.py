"""tTest-anchored cluster trace visualization on the 4-phase dendrogram.

Cut the **tTest** tree at K clusters via ``fcluster(maxclust=K)`` and
colour every leaf by its tTest-cluster id. The *same* per-leaf colours
are then drawn under all four phases. Every leaf is coloured (no
"unmatched" → no leaf is grey).

Trace reads as **spatial coherence** of each colour across phases:

- On **tTest**: the K colours form K contiguous blocks by construction
  (the cut defines them).
- On **rPost**: each tTest-cluster's leaves are either contiguous
  (the clade survives as a clade in rPost) or scattered (the clade
  broke up). Likewise the U-shapes that envelope a single tTest-cluster
  inherit its colour; mixed U-shapes are grey.
- On **rPre**: usually heavy fragmentation — the same colours appear
  scattered across the bottom strip because rPre's tree groups leaves
  differently.
- On **tLearn**: intermediate.

This is the bias-free version: no per-subtree score, no asymmetric
formula — just one cluster cut on tTest, propagated to every phase.

Usage:
    python … --patient Pat_06 --band beta --k 20
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, fcluster

from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.lrg import load_lrg_result

use_lrg_style()


COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
PHASE_TITLES = {
    "rest_pre": r"$\mathrm{rPre}$",
    "task_learn": r"$\mathrm{tLearn}$",
    "task_test": r"$\mathrm{tTest}$",
    "rest_post": r"$\mathrm{rPost}$",
}
FC = "imcoh_abs"
UNMATCHED_COLOR = "#dddddd"
HETEROGENEOUS_COLOR = "black"
ANCHOR_PHASE = "rest_post"

ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "data/reports/preprint/figure1_beta_trace/dendrograms_matched_clusters"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def leaves_under_each_node(Z: np.ndarray) -> dict[int, frozenset[int]]:
    n = Z.shape[0] + 1
    out: dict[int, frozenset[int]] = {i: frozenset({i}) for i in range(n)}
    for i in range(n - 1):
        a = int(Z[i, 0])
        b = int(Z[i, 1])
        out[n + i] = out[a] | out[b]
    return out


def render_patient_page(patient: str, band: str, k: int, pdf: PdfPages) -> None:
    try:
        Zs = {p: load_lrg_result(patient, p, band, FC).linkage_matrix for p in PHASES}
    except Exception as err:
        print(f"  [{patient} {band}] skip ({err})")
        return
    n = Zs[PHASES[0]].shape[0] + 1

    # Cut the ANCHOR_PHASE tree at K — these clusters define the per-leaf
    # colors propagated to every panel.
    anchor_labels = fcluster(Zs[ANCHOR_PHASE], t=k, criterion="maxclust")
    unique_clusters = np.unique(anchor_labels)
    n_clusters = unique_clusters.size
    sizes = [int((anchor_labels == c).sum()) for c in unique_clusters]
    print(f"  [{patient} {band}] {ANCHOR_PHASE} cut K={k} → "
          f"{n_clusters} clusters, sizes={sorted(sizes, reverse=True)}")

    # Re-rank cluster ids by leaf count (largest first).
    order = sorted(
        unique_clusters,
        key=lambda c: (-int((anchor_labels == c).sum()),
                       int(np.where(anchor_labels == c)[0].min())),
    )
    cluster_to_color: dict[int, int] = {int(c): r for r, c in enumerate(order)}

    leaf_color_id = np.array(
        [cluster_to_color[int(c)] for c in anchor_labels], dtype=int
    )

    # Palette: tab20 has 20 distinguishable colors. For K > 20, cycle.
    palette = plt.get_cmap("tab20")
    color_for_id: dict[int, str] = {
        cid: mcolors.to_hex(palette(cid % palette.N))
        for cid in range(max(n_clusters, 1))
    }

    leaf_colors = [color_for_id[int(c)] for c in leaf_color_id]

    # Per-node leaves for each phase (for link homogeneity test).
    leaves_per_phase = {p: leaves_under_each_node(Z) for p, Z in Zs.items()}

    def make_link_color_func(phase: str):
        leaves_dict = leaves_per_phase[phase]

        def cf(node_id: int) -> str:
            leaves = leaves_dict[node_id]
            color_set = {int(leaf_color_id[L]) for L in leaves}
            if len(color_set) == 1:
                return color_for_id[color_set.pop()]
            return HETEROGENEOUS_COLOR

        return cf

    # Layout: GridSpec 2 rows (dendro, leaf strip), 4 cols (phases).
    fig = plt.figure(figsize=(13.5, 4.2))
    gs = fig.add_gridspec(
        2, 4,
        height_ratios=[14, 0.6],
        hspace=0.02, wspace=0.10,
        left=0.05, right=0.99, top=0.84, bottom=0.06,
    )
    axes_top = [fig.add_subplot(gs[0, 0])]
    for j in range(1, 4):
        axes_top.append(fig.add_subplot(gs[0, j], sharey=axes_top[0]))
    axes_bot = [fig.add_subplot(gs[1, j], sharex=axes_top[j]) for j in range(4)]

    tmin = min(Z[0, 2] for Z in Zs.values()) * 0.8
    tmax = max(Z[-1, 2] for Z in Zs.values()) * 1.05

    for col, phase in enumerate(PHASES):
        ax_top = axes_top[col]
        ax_bot = axes_bot[col]
        Z = Zs[phase]

        ax_top.set_title(PHASE_TITLES[phase])
        d = dendrogram(
            Z, ax=ax_top,
            no_labels=True,
            link_color_func=make_link_color_func(phase),
            above_threshold_color=HETEROGENEOUS_COLOR,
        )
        # Per-leaf colored stub at the bottom of each leaf's leg, so
        # singletons (whose first U-shape is heterogeneous → black) are
        # still visible at the leaf level in the dendrogram itself.
        # scipy places leaf positions at x = 5, 15, …, 5 + 10·(n-1).
        leaf_order = d["leaves"]
        leaf_xs = [5 + 10 * i for i in range(len(leaf_order))]
        stub_top = tmin + 0.025 * (tmax - tmin)
        for x, L in zip(leaf_xs, leaf_order):
            ax_top.plot(
                [x, x], [tmin, stub_top],
                color=leaf_colors[L], linewidth=2.0,
                solid_capstyle="butt", zorder=20,
            )
        ax_top.set_ylim(tmin, tmax)
        ax_top.set_xticks([])
        ax_top.set_xlabel("")
        if col == 0:
            ax_top.set_ylabel(r"merge height  $\hat{D}(\tau)$")
        ax_top.spines["top"].set_visible(False)
        ax_top.spines["right"].set_visible(False)

        # Leaf strip: each leaf in display order at x = 5, 15, 25, ...
        leaves_in_order = d["leaves"]
        for i, L in enumerate(leaves_in_order):
            ax_bot.axvspan(10 * i, 10 * (i + 1), color=leaf_colors[L], lw=0)
        ax_bot.set_xlim(0, 10 * n)
        ax_bot.set_xticks([])
        ax_bot.set_yticks([])
        for spine in ax_bot.spines.values():
            spine.set_visible(False)

    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    fig.text(
        0.05, 0.95,
        rf"{patient}  |  {band_tex}  |  {ANCHOR_PHASE}-anchored clusters "
        rf"($K = {k}$, {n_clusters} actual)  |  $N = {n}$",
        fontweight="bold",
    )
    fig.text(
        0.05, 0.91,
        f"every leaf colored by its {ANCHOR_PHASE}-cluster id; same colors "
        "used in all 4 panels; U-shapes painted black when descendants span "
        ">1 anchor-cluster",
        fontsize=7, color="#444",
    )

    pdf.savefig(fig)
    plt.close(fig)


def build_band_pdf(band: str, k: int) -> Path:
    anchor_tag = ANCHOR_PHASE.replace("_", "")
    out_pdf = OUT_DIR / f"matched_{anchor_tag}_K{k}_{band}_all_patients.pdf"
    print(f"\n=== {band} → {out_pdf.name} ===")
    with PdfPages(out_pdf) as pdf:
        for pat in COHORT:
            render_patient_page(pat, band, k, pdf)
    print(f"DONE: {out_pdf}")
    return out_pdf


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--band", default=None,
                        help="single band (default: all 6 bands)")
    parser.add_argument("--k", type=int, default=20)
    args = parser.parse_args()
    bands = [args.band] if args.band else list(BRAIN_BANDS.keys())
    for band in bands:
        build_band_pdf(band, args.k)


if __name__ == "__main__":
    main()
