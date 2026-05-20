"""Per-patient hierarchical trees (cut at k=20) across the 4 phases, all bands.

One PDF per band, 10 pages per PDF, each page = 1×4 subplots showing the
LRG dendrogram at rPre / tLearn / tTest / rPost.  At each phase, the tree
is cut into K_CUT=20 flat clusters via ``fcluster(maxclust=20)`` and the
cluster ids are re-ranked by leaf count (largest → palette[0], smallest →
palette[19]) so the same colour means "the k-th largest cluster" in every
phase.  Internal links above the cut (whose descendants span >1 cluster)
are drawn black.

Output:
    data/reports/preprint/figure1_beta_trace/dendrograms_k20/
        dendrograms_k20_<band>_all_patients.pdf  (×6, one per band)

Preliminary view — later figures in a separate folder will try to make
the "memory" blocks / leaf groups emerge more cleanly.
"""
from pathlib import Path

import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, fcluster

from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.visuals.styles import use_lrg_style

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
BANDS = list(BRAIN_BANDS.keys())
FC = "imcoh_abs"
K_CUT = 20
ABOVE_COLOR = "black"

TAB20_HEX = [mcolors.to_hex(c) for c in plt.cm.tab20.colors]

ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "data/reports/preprint/figure1_beta_trace/dendrograms_k20"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _leaves_under_each_node(Z: np.ndarray) -> dict[int, frozenset[int]]:
    """Map every node id (leaf or merge) → frozenset of descendant leaf ids."""
    n = Z.shape[0] + 1
    leaves: dict[int, frozenset[int]] = {i: frozenset({i}) for i in range(n)}
    for i in range(n - 1):
        a = int(Z[i, 0])
        b = int(Z[i, 1])
        leaves[n + i] = leaves[a] | leaves[b]
    return leaves


def _make_link_color_func(Z: np.ndarray, k: int):
    """Build a scipy ``link_color_func`` that colours each link by its
    descendants' fcluster id, **re-ranked by cluster size** so that
    palette[0] is always the largest cluster, palette[1] the 2nd largest,
    etc.  Ties broken by smallest leaf id (deterministic).  This makes the
    colour sequence comparable across phases of the same patient.
    """
    flat = fcluster(Z, t=k, criterion="maxclust")  # 1-indexed
    leaves = _leaves_under_each_node(Z)

    # Rank original cluster ids (1..k_actual) by descending leaf count.
    unique_ids, counts = np.unique(flat, return_counts=True)
    min_leaf = {int(c): int(np.where(flat == c)[0].min()) for c in unique_ids}
    order = sorted(
        unique_ids.tolist(),
        key=lambda c: (-int(counts[unique_ids.tolist().index(c)]), min_leaf[int(c)]),
    )
    rank: dict[int, int] = {int(orig): r for r, orig in enumerate(order)}

    def color_for(node_id: int) -> str:
        cluster_ids = {flat[leaf] for leaf in leaves[node_id]}
        if len(cluster_ids) == 1:
            orig = int(cluster_ids.pop())
            return TAB20_HEX[rank[orig] % len(TAB20_HEX)]
        return ABOVE_COLOR

    return color_for


def _band_glyph(band: str) -> str:
    tex = BRAIN_BAND_TEX_DICT.get(band, band)
    return tex if isinstance(tex, str) else str(tex)


def page_for_patient(patient: str, band: str, pdf: PdfPages) -> None:
    Z_per_phase: dict[str, np.ndarray | None] = {}
    for phase in PHASES:
        try:
            lrg = load_lrg_result(patient, phase, band, FC)
            Z_per_phase[phase] = lrg.linkage_matrix
        except Exception as err:
            print(f"  [{patient} {band} {phase}] missing: {err}")
            Z_per_phase[phase] = None

    available = {p: Z for p, Z in Z_per_phase.items() if Z is not None}
    if not available:
        print(f"[{patient} {band}] no data, skipping")
        return

    # Per-rule: tmin = merge_heights[0]*0.8, tmax = merge_heights[-1]*1.05.
    # Shared y-axis across the 4 phases ⇒ take cohort min/max across phases.
    tmin = min(Z[0, 2] for Z in available.values()) * 0.8
    tmax = max(Z[-1, 2] for Z in available.values()) * 1.05

    fig, axes = plt.subplots(
        1, 4, figsize=(13.5, 3.6),
        sharey=True,
    )
    fig.subplots_adjust(left=0.06, right=0.99, top=0.86, bottom=0.06, wspace=0.10)

    for col, (ax, phase) in enumerate(zip(axes, PHASES)):
        ax.set_title(PHASE_TITLES[phase])
        Z = Z_per_phase[phase]
        if Z is None:
            ax.text(0.5, 0.5, "missing", ha="center", va="center",
                    transform=ax.transAxes, color="#888")
            ax.set_xticks([])
            ax.set_yticks([])
            continue
        color_func = _make_link_color_func(Z, K_CUT)
        dendrogram(
            Z,
            ax=ax,
            no_labels=True,
            link_color_func=color_func,
            above_threshold_color=ABOVE_COLOR,
        )
        ax.set_ylim(tmin, tmax)
        ax.set_xticks([])
        ax.set_xlabel("")
        if col == 0:
            ax.set_ylabel(r"merge height  $\hat{D}(\tau)$")
        # Tighten spines: keep left + bottom only on leftmost; bottom on all.
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    n = next(iter(available.values())).shape[0] + 1
    fig.text(
        0.06, 0.94,
        rf"{patient}  |  {_band_glyph(band)}  |  cut $k = {K_CUT}$  |  "
        rf"$\tau = 1/\lambda_{{\max}}$  |  $N = {n}$",
        fontweight="bold",
    )

    pdf.savefig(fig)
    plt.close(fig)


def build_band_pdf(band: str) -> Path:
    out_pdf = OUT_DIR / f"dendrograms_k20_{band}_all_patients.pdf"
    print(f"\n=== {band} → {out_pdf.name} ===")
    with PdfPages(out_pdf) as pdf:
        for pat in COHORT:
            print(f"[{pat}]")
            page_for_patient(pat, band, pdf)
    print(f"DONE: {out_pdf}")
    return out_pdf


def main() -> None:
    for band in BANDS:
        build_band_pdf(band)


if __name__ == "__main__":
    main()
