"""Pairwise tanglegrams against rPost.

Three textbook tanglegrams per page, one per non-anchor phase:

    rPre     ↔  rPost
    tLearn   ↔  rPost
    tTest    ↔  rPost

Each tanglegram = two dendrograms facing each other (left tree's leaves
on its right edge, right tree's leaves on its left edge) with a colored
connection line per leaf bridging the gap. **Fewer line crossings =
more topologically similar trees** (= fewer leaf-order disagreements).
Crossings count printed as the headline number per row.

Reading: the row with the fewest crossings to rPost is the phase whose
tree is most leaf-permutation-similar to rPost. The persistence/trace
question reduces to: does tTest↔rPost have visibly fewer crossings
than rPre↔rPost?

Leaf orderings: scipy's natural orderings on each tree, no optimization.
A tanglegram-untangling pass would reduce crossings by sibling-swapping
internal nodes; not done here. The relative crossings count is still
meaningful because the same heuristic is applied to all trees.

Line color: probe-shaft identity (same as the existing matched-trace
barcode). Each color = one sEEG shaft.

Output:
    data/reports/preprint/figure1_beta_trace/tanglegrams_vs_rpost/
        tanglegrams_vs_rpost_{band}_all_patients.pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.io.patient import load_channel_labels
from lrg_eegfc.utils.probe import extract_probe_labels
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.lrg import load_lrg_result

use_lrg_style()


COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
PHASES_LEFT = ["rest_pre", "task_learn", "task_test"]
PHASE_TITLES = {
    "rest_pre": r"$\mathrm{rPre}$",
    "task_learn": r"$\mathrm{tLearn}$",
    "task_test": r"$\mathrm{tTest}$",
    "rest_post": r"$\mathrm{rPost}$",
}
FC = "imcoh_abs"
ANCHOR_PHASE = "rest_post"

ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "data/reports/preprint/figure1_beta_trace/tanglegrams_vs_rpost"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def count_crossings(pairs: list[tuple[float, float]]) -> int:
    """Number of crossings between lines connecting (y_left_i, y_right_i)."""
    n_cross = 0
    for i in range(len(pairs)):
        yli, yri = pairs[i]
        for j in range(i + 1, len(pairs)):
            ylj, yrj = pairs[j]
            if (yli - ylj) * (yri - yrj) < 0:
                n_cross += 1
    return n_cross


def reorder_linkage_to_target(
    Z: np.ndarray,
    target_leaf_pos: dict[int, float],
) -> np.ndarray:
    """Sibling-swap each internal node of Z so the visual leaf order is
    pulled toward `target_leaf_pos` (= scipy display position of each
    leaf in the reference tree). Greedy bottom-up by mean target
    position of descendants. Topology is preserved; only sibling
    orientation changes."""
    n = Z.shape[0] + 1
    leaves_under: dict[int, list[int]] = {i: [i] for i in range(n)}
    for i in range(n - 1):
        a, b = int(Z[i, 0]), int(Z[i, 1])
        leaves_under[n + i] = leaves_under[a] + leaves_under[b]

    mean_target: dict[int, float] = {
        nid: float(np.mean([target_leaf_pos[L] for L in leaves]))
        for nid, leaves in leaves_under.items()
    }

    Z_new = Z.copy()
    for i in range(n - 1):
        a, b = int(Z[i, 0]), int(Z[i, 1])
        if mean_target[a] > mean_target[b]:
            Z_new[i, 0] = b
            Z_new[i, 1] = a
    return Z_new


def render_tanglegram_row(
    fig,
    gs_row,
    Z_left: np.ndarray,
    Z_right: np.ndarray,
    n: int,
    probe_color_per_leaf: list[str],
    label_left: str,
    label_right: str,
) -> int:
    ax_left = fig.add_subplot(gs_row[0])
    ax_mid = fig.add_subplot(gs_row[1])
    ax_right = fig.add_subplot(gs_row[2])

    d_left = dendrogram(
        Z_left, ax=ax_left,
        no_labels=True,
        orientation="left",
        link_color_func=lambda _k: "#222222",
        above_threshold_color="#222222",
    )
    d_right = dendrogram(
        Z_right, ax=ax_right,
        no_labels=True,
        orientation="right",
        link_color_func=lambda _k: "#222222",
        above_threshold_color="#222222",
    )

    # scipy places leaves at y = 5, 15, ... regardless of left/right.
    y_left_of = {int(L): i * 10 + 5 for i, L in enumerate(d_left["leaves"])}
    y_right_of = {int(L): i * 10 + 5 for i, L in enumerate(d_right["leaves"])}

    y_pairs = [(y_left_of[L], y_right_of[L]) for L in range(n)]
    n_cross = count_crossings(y_pairs)

    ax_mid.set_xlim(0, 1)
    ax_mid.set_ylim(0, n * 10)
    for L in range(n):
        yl, yr = y_pairs[L]
        ax_mid.plot(
            [0.02, 0.98], [yl, yr],
            color=probe_color_per_leaf[L],
            lw=0.6, alpha=0.7,
            solid_capstyle="butt",
        )
    ax_mid.axis("off")

    ax_left.set_title(label_left, fontsize=8, loc="left")
    ax_right.set_title(label_right, fontsize=8, loc="right")
    for ax in (ax_left, ax_right):
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

    # Crossings headline at the top of the middle gap.
    ax_mid.text(
        0.5, 1.02,
        f"crossings: {n_cross}",
        transform=ax_mid.transAxes,
        ha="center", va="bottom",
        fontsize=8, fontweight="bold",
        color="#222222",
    )
    return n_cross


def render_patient_page(patient: str, band: str, pdf: PdfPages) -> None:
    try:
        Zs = {p: load_lrg_result(patient, p, band, FC).linkage_matrix
              for p in PHASES_LEFT + [ANCHOR_PHASE]}
    except Exception as err:
        print(f"  [{patient} {band}] skip ({err})")
        return
    n = Zs[ANCHOR_PHASE].shape[0] + 1

    # Anchor (rPost) keeps its scipy-default leaf order. Every other
    # phase's tree gets its linkage matrix reordered (sibling-swaps
    # only — topology preserved) so its visual leaf order is pulled
    # toward rPost's. This minimizes connection-line crossings within
    # the constraints imposed by each tree's topology; the residual
    # crossings count is then a fair similarity comparison across the
    # 3 non-anchor phases.
    d_anchor = dendrogram(Zs[ANCHOR_PHASE], no_plot=True)
    anchor_pos = {int(L): i * 10 + 5 for i, L in enumerate(d_anchor["leaves"])}
    for phase in PHASES_LEFT:
        Zs[phase] = reorder_linkage_to_target(Zs[phase], anchor_pos)

    channel_labels = load_channel_labels(patient)
    probe_per_leaf = extract_probe_labels(channel_labels)
    unique_probes = sorted(set(probe_per_leaf))
    palette = plt.get_cmap("tab20")
    probe_to_color = {
        pr: mcolors.to_hex(palette(i % palette.N))
        for i, pr in enumerate(unique_probes)
    }
    probe_color_per_leaf = [probe_to_color[pr] for pr in probe_per_leaf]

    fig = plt.figure(figsize=(13.5, 10.5))
    gs = fig.add_gridspec(
        3, 3,
        width_ratios=[4, 2, 4],
        hspace=0.30, wspace=0.04,
        left=0.04, right=0.99, top=0.92, bottom=0.04,
    )

    crossings: dict[str, int] = {}
    for row, phase in enumerate(PHASES_LEFT):
        gs_row = [gs[row, c] for c in range(3)]
        crossings[phase] = render_tanglegram_row(
            fig, gs_row,
            Zs[phase], Zs[ANCHOR_PHASE], n,
            probe_color_per_leaf,
            PHASE_TITLES[phase], PHASE_TITLES[ANCHOR_PHASE],
        )

    most_similar = min(crossings, key=crossings.get)
    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    headline = (
        rf"{patient}  |  {band_tex}  |  tanglegrams vs rPost  |  $N={n}$"
        rf"  |  crossings — "
        rf"rPre:{crossings['rest_pre']}  "
        rf"tLearn:{crossings['task_learn']}  "
        rf"tTest:{crossings['task_test']}  "
        rf"|  most similar to rPost: {PHASE_TITLES[most_similar]}"
    )
    fig.text(0.04, 0.97, headline, fontweight="bold")
    fig.text(
        0.04, 0.945,
        "Fewer line crossings between two trees ⇒ leaves end up in similar "
        "positions in both ⇒ more topologically similar.  Line color: probe "
        "(shaft) identity.",
        fontsize=7, color="#444",
    )

    print(f"  [{patient} {band}]  crossings  rPre={crossings['rest_pre']}  "
          f"tLearn={crossings['task_learn']}  "
          f"tTest={crossings['task_test']}  →  most-similar={most_similar}")

    pdf.savefig(fig)
    plt.close(fig)


def build_band_pdf(band: str) -> Path:
    out_pdf = OUT_DIR / f"tanglegrams_vs_rpost_{band}_all_patients.pdf"
    print(f"\n=== {band} → {out_pdf.name} ===")
    with PdfPages(out_pdf) as pdf:
        for pat in COHORT:
            render_patient_page(pat, band, pdf)
    print(f"DONE: {out_pdf}")
    return out_pdf


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--band", default=None,
                        help="single band (default: all 6 bands)")
    args = parser.parse_args()
    bands = [args.band] if args.band else list(BRAIN_BANDS.keys())
    for band in bands:
        build_band_pdf(band)


if __name__ == "__main__":
    main()
