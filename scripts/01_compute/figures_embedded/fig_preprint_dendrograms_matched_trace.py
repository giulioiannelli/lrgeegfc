"""Matched-trace-clades coloring on the 4-phase dendrogram.

Identifies *preserved clades* between tTest and rPost — pairs
(S_test, S_post) of internal subtrees with Jaccard ≥ TAU — and gives
each pair a **unique color**. The same color appears on the matching
links in both the tTest and rPost panels → visual pairing by eye.

For each preserved clade, asks the further question: does the same
leaf-set ALSO form a clade in rPre (Jaccard ≥ TAU_PRE)?
- If yes → **anchor**: rPre panel inherits the same color (muted).
- If no  → **trace**: rPre panel stays gray for those leaves; the
  *absence* of the color in rPre is the visual signal of the CHANGE.

tLearn is colored by the same scheme (homogeneity rule against the
test↔post clade assignment), so you can see when the consolidation
clade already exists during learning.

No continuous Jaccard scalar is plotted. The matching uses set overlap
internally to decide kept/not-kept (one binary decision per candidate
pair), but the figure shows only the resulting *colored clade
identities*.

Output:
    data/reports/preprint/figure1_beta_trace/dendrograms_matched_trace/
        matched_trace_{band}_all_patients.pdf
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
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
PHASE_TITLES = {
    "rest_pre": r"$\mathrm{rPre}$",
    "task_learn": r"$\mathrm{tLearn}$",
    "task_test": r"$\mathrm{tTest}$",
    "rest_post": r"$\mathrm{rPost}$",
}
FC = "imcoh_abs"
MIN_CLADE_SIZE = 3   # discard pair clades (too small to be meaningful)
UNMATCHED = "#000000"

ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "data/reports/preprint/figure1_beta_trace/dendrograms_matched_trace"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def leaves_under_each_node(Z: np.ndarray) -> dict[int, frozenset[int]]:
    n = Z.shape[0] + 1
    out: dict[int, frozenset[int]] = {i: frozenset({i}) for i in range(n)}
    for i in range(n - 1):
        a = int(Z[i, 0])
        b = int(Z[i, 1])
        out[n + i] = out[a] | out[b]
    return out


def find_preserved_clades(
    leaves_test: dict[int, frozenset[int]],
    leaves_post: dict[int, frozenset[int]],
) -> list[dict]:
    """Exact-equality clade matching between test and post.

    A clade is preserved iff a subtree of rPost has the *identical*
    leaf-set as some subtree of tTest (set equality, not Jaccard ≥ τ).
    Greedy non-overlapping selection, largest first, so nested matches
    keep the largest representative.

    Returns a list of dicts with keys: nid_test, nid_post, leaves
    (frozenset, identical in both trees by construction).
    """
    n_all = max(max(S) for S in leaves_test.values()) + 1
    # Cap at half the tree so each bipartition is represented only by
    # its smaller side. Without this cap, "everything minus a few
    # outlier leaves" (e.g., 110 of 113) registers as a preserved
    # clade — that's the near-trivial root split, not a real clade.
    max_size = n_all // 2
    post_index: dict[frozenset[int], int] = {
        T: nid_t
        for nid_t, T in leaves_post.items()
        if MIN_CLADE_SIZE <= len(T) <= max_size
    }

    cands: list[tuple] = []
    for nid_s, S in leaves_test.items():
        if not (MIN_CLADE_SIZE <= len(S) <= max_size):
            continue
        if S in post_index:
            cands.append((len(S), nid_s, post_index[S], S))
    # Largest first → nested exact matches keep the outer clade.
    cands.sort(key=lambda c: -c[0])

    claimed: set[int] = set()
    kept: list[dict] = []
    for sz, nid_s, nid_t, S in cands:
        if S & claimed:
            continue
        claimed |= S
        kept.append({"nid_test": nid_s, "nid_post": nid_t, "leaves": S})
    return kept


def clade_in_phase(canonical: frozenset[int],
                    leaves_other: dict[int, frozenset[int]]) -> bool:
    """True iff some subtree of `leaves_other` has the exact leaf-set."""
    return any(canonical == V for V in leaves_other.values())


def lighten(hex_color: str, frac: float = 0.55) -> str:
    """Mix `hex_color` with white at `frac` weight (anchor=lighter version)."""
    r, g, b = mcolors.to_rgb(hex_color)
    return mcolors.to_hex((
        r + frac * (1.0 - r),
        g + frac * (1.0 - g),
        b + frac * (1.0 - b),
    ))


def render_patient_page(patient: str, band: str, pdf: PdfPages) -> None:
    try:
        Zs = {p: load_lrg_result(patient, p, band, FC).linkage_matrix for p in PHASES}
    except Exception as err:
        print(f"  [{patient} {band}] skip ({err})")
        return
    n = Zs[PHASES[0]].shape[0] + 1
    leaves_per_phase = {p: leaves_under_each_node(Z) for p, Z in Zs.items()}

    kept = find_preserved_clades(
        leaves_per_phase["task_test"],
        leaves_per_phase["rest_post"],
    )

    # Anchor iff the same leaf-set is also an exact clade in rPre.
    for c in kept:
        c["anchor"] = clade_in_phase(c["leaves"], leaves_per_phase["rest_pre"])

    palette = (
        [mcolors.to_hex(c) for c in plt.cm.tab20.colors]
        + [mcolors.to_hex(c) for c in plt.cm.tab20b.colors]
        + [mcolors.to_hex(c) for c in plt.cm.tab20c.colors]
    )
    for i, c in enumerate(kept):
        c["color"] = palette[i % len(palette)]

    n_anchor = sum(1 for c in kept if c["anchor"])
    n_trace = sum(1 for c in kept if not c["anchor"])
    sizes = sorted([len(c["leaves"]) for c in kept], reverse=True)
    print(f"  [{patient} {band}] exact-match clades: {len(kept)}  "
          f"(trace: {n_trace},  anchor: {n_anchor})  sizes: {sizes}")

    # Leaves of a preserved clade are identical in test and post by
    # construction → one leaf→clade map suffices, applied uniformly to
    # every panel. Same color for the same leaf across all 4 panels.
    leaf_to_clade = np.full(n, -1, dtype=int)
    for i, c in enumerate(kept):
        for L in c["leaves"]:
            leaf_to_clade[L] = i

    def make_color_func(phase: str):
        # Uniform homogeneity rule across ALL 4 panels: a link is
        # colored iff all its descendant leaves share one clade id
        # (= the leaves of that one preserved test↔post clade also
        # form a clade in THIS phase's tree). Otherwise black.
        # In rPre, trace clades' leaves are scattered → no rPre
        # subtree is homogeneous for them → those leaves' clade color
        # never appears on rPre links. Anchor clades still do.
        leaves_dict = leaves_per_phase[phase]

        def cf(node_id: int) -> str:
            leaves_under = leaves_dict[node_id]
            clade_ids = {int(leaf_to_clade[L]) for L in leaves_under}
            if len(clade_ids) != 1:
                return UNMATCHED
            cid = clade_ids.pop()
            if cid == -1:
                return UNMATCHED
            return kept[cid]["color"]

        return cf

    tmin = min(Z[0, 2] for Z in Zs.values()) * 0.8
    tmax = max(Z[-1, 2] for Z in Zs.values()) * 1.05

    fig = plt.figure(figsize=(13.5, 4.0))
    gs = fig.add_gridspec(
        2, 4,
        height_ratios=[14, 0.6],
        hspace=0.02, wspace=0.10,
        left=0.06, right=0.99, top=0.84, bottom=0.06,
    )
    axes_top = [fig.add_subplot(gs[0, 0])]
    for j in range(1, 4):
        axes_top.append(fig.add_subplot(gs[0, j], sharey=axes_top[0]))
    axes_bot = [fig.add_subplot(gs[1, j], sharex=axes_top[j]) for j in range(4)]

    BARCODE_UNASSIGNED = "#dddddd"

    for col, phase in enumerate(PHASES):
        ax_top = axes_top[col]
        ax_bot = axes_bot[col]
        ax_top.set_title(PHASE_TITLES[phase])
        Z = Zs[phase]
        d = dendrogram(
            Z, ax=ax_top,
            no_labels=True,
            link_color_func=make_color_func(phase),
            above_threshold_color=UNMATCHED,
        )
        leaf_order = d["leaves"]
        leaf_xs = [5 + 10 * i for i in range(len(leaf_order))]
        # Overlay colored leg for every assigned leaf up to its first
        # merge height — so even isolated assignments show on the tree.
        leaf_parent_height = np.zeros(n, dtype=float)
        for i in range(Z.shape[0]):
            for child in (int(Z[i, 0]), int(Z[i, 1])):
                if child < n:
                    leaf_parent_height[child] = Z[i, 2]
        for x, L in zip(leaf_xs, leaf_order):
            cid = int(leaf_to_clade[L])
            if cid < 0:
                continue
            color = kept[cid]["color"]
            ax_top.plot(
                [x, x], [tmin, leaf_parent_height[L]],
                color=color, linewidth=1.2,
                solid_capstyle="butt", zorder=20,
            )
        ax_top.set_ylim(tmin, tmax)
        ax_top.set_xticks([])
        ax_top.set_xlabel("")
        if col == 0:
            ax_top.set_ylabel(r"merge height  $\hat{D}(\tau)$")
        ax_top.spines["top"].set_visible(False)
        ax_top.spines["right"].set_visible(False)

        # Colorbar-style barcode aligned with leaves below the
        # dendrogram. Each leaf in its display order is one cell.
        for i, L in enumerate(leaf_order):
            cid = int(leaf_to_clade[L])
            color = kept[cid]["color"] if cid >= 0 else BARCODE_UNASSIGNED
            ax_bot.axvspan(10 * i, 10 * (i + 1), color=color, lw=0)
        ax_bot.set_xlim(0, 10 * n)
        ax_bot.set_xticks([])
        ax_bot.set_yticks([])
        for spine in ax_bot.spines.values():
            spine.set_visible(False)

    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    fig.text(
        0.06, 0.95,
        rf"{patient}  |  {band_tex}  |  exact test↔post clade match  |  "
        rf"$N={n}$  |  preserved={len(kept)} "
        rf"(trace={n_trace}, anchor={n_anchor})",
        fontweight="bold",
    )
    fig.text(
        0.06, 0.92,
        "Colors = preserved (test↔post) clades.  Same color = same leaf "
        "across all panels.  Link colored iff all its descendant leaves "
        "share one clade in THIS phase's tree → in rPre, trace-clade leaves "
        "scatter → those colors only appear on the stubs, not on the links "
        "(visible CHANGE).  Default link color: black.",
        fontsize=7, color="#444",
    )

    pdf.savefig(fig)
    plt.close(fig)


def build_band_pdf(band: str) -> Path:
    out_pdf = OUT_DIR / f"matched_trace_{band}_all_patients.pdf"
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
