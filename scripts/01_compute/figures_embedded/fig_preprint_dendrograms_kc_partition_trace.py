"""KC λ-mixed partition trace coloring (per-merge-node aggregation).

Difference from ``fig_preprint_dendrograms_kc_trace`` (per-subtree mean,
overlapping): here each leaf pair ``(i,j)`` is credited to *exactly* its
MRCA node in the reference phase's tree, and never to ancestor subtrees.
The cross-product ``L_left(N) × L_right(N)`` = pairs first co-merging at
N → strict partition of all leaf pairs across all merge nodes. Therefore
the per-merge-node mean is a faithful decomposition of the per-pair KC
trace sum, with no upward smearing into ancestor U-shapes.

Per-pair KC λ-mixed score (tree-invariant, normalised pooled-max):
    v_X(i,j)  =  (1-λ)·m_X(i,j)/m_max + λ·M_X(i,j)/M_max
    s(i,j)    =  ½(|v_pre − v_test| + |v_pre − v_post|)
                 −  |v_test − v_post|

Per-merge-node score (phase-specific, partition):
    trace_X(N) = mean of s(i,j) over (i,j) ∈ L_left(N) × L_right(N)

Visualization: 4 phase panels (rPre / tLearn / tTest / rPost), U-shapes
coloured by trace_X(N) (RdBu_r divergent, symmetric vlim pooled across
phases for the page). Probe-shaft barcode below each panel for leaf
identity. Each panel's reference tree is its own phase's linkage; the
per-pair s is shared (tree-invariant).

Scope report: `.agents/guides/task-persistence-investigation/
2026-05-11_kc-partition-merge-node-trace.md`.

Output:
    data/reports/preprint/figure1_beta_trace/dendrograms_kc_partition_trace/
        kc_partition_trace_lambda{lam}_{band}_all_patients.pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.io.patient import load_channel_labels
from lrg_eegfc.utils.metrics.tree_distance import kc_vectors
from lrg_eegfc.utils.probe import extract_probe_labels
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrgsglib.plotlib import imshow_colorbar_caxdivider

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
UNCOLORED = "#888888"

ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "data/reports/preprint/figure1_beta_trace/dendrograms_kc_partition_trace"
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


def per_merge_node_score(Z: np.ndarray, s_mat: np.ndarray) -> dict[int, float]:
    """Partition decomposition: each pair credited to its MRCA only."""
    n = Z.shape[0] + 1
    leaves = leaves_under_each_node(Z)
    out: dict[int, float] = {}
    for i in range(n - 1):
        L = list(leaves[int(Z[i, 0])])
        R = list(leaves[int(Z[i, 1])])
        sub = s_mat[np.ix_(L, R)]
        out[n + i] = float(sub.mean())
    return out


def _lam_tag(lam: float) -> str:
    return f"{lam:g}".replace(".", "p")


def render_patient_page(patient: str, band: str, lam: float, pdf: PdfPages) -> None:
    try:
        Zs = {p: load_lrg_result(patient, p, band, FC).linkage_matrix for p in PHASES}
    except Exception as err:
        print(f"  [{patient} {band}] skip ({err})")
        return
    n = Zs[PHASES[0]].shape[0] + 1

    # KC λ-mixed per-pair vectors with pooled normalisation.
    m_vec = {p: kc_vectors(Z)[0].astype(float) for p, Z in Zs.items()}
    M_vec = {p: kc_vectors(Z)[1].astype(float) for p, Z in Zs.items()}
    mmax = max(v.max() for v in m_vec.values()) or 1.0
    Mmax = max(v.max() for v in M_vec.values()) or 1e-12
    v_pair = {
        p: (1.0 - lam) * m_vec[p] / mmax + lam * M_vec[p] / Mmax
        for p in PHASES
    }
    Vmat = {p: m_vec_to_pair_matrix(v_pair[p], n) for p in PHASES}

    V_pre = Vmat["rest_pre"]
    V_test = Vmat["task_test"]
    V_post = Vmat["rest_post"]
    s = 0.5 * (np.abs(V_pre - V_test) + np.abs(V_pre - V_post)) \
        - np.abs(V_test - V_post)

    # Per-merge-node trace, one map per phase.
    scores_per_phase = {p: per_merge_node_score(Zs[p], s) for p in PHASES}

    all_scores = np.fromiter(
        (v for d in scores_per_phase.values() for v in d.values()),
        dtype=float,
    )
    vmax = max(float(np.max(np.abs(all_scores))), 0.05)
    vmin = -vmax
    cmap = plt.get_cmap("RdBu_r")
    norm = Normalize(vmin=vmin, vmax=vmax)

    print(f"  [{patient} {band}] λ={lam:.2f}  s pair-range "
          f"[{s.min():+.3f}, {s.max():+.3f}]  "
          f"per-merge-node vlim ±{vmax:.3f}")

    def color_for_phase(phase: str):
        scores = scores_per_phase[phase]

        def cf(node_id: int) -> str:
            sc = scores.get(node_id)
            if sc is None or not np.isfinite(sc):
                return UNCOLORED
            return mcolors.to_hex(cmap(norm(sc)))

        return cf

    # Probe-shaft leaf identity for the barcode.
    channel_labels = load_channel_labels(patient)
    probe_per_leaf = extract_probe_labels(channel_labels)
    unique_probes = sorted(set(probe_per_leaf))
    stub_palette = plt.get_cmap("tab20")
    probe_to_color = {
        pr: mcolors.to_hex(stub_palette(i % stub_palette.N))
        for i, pr in enumerate(unique_probes)
    }
    leaf_stub_color = [probe_to_color[pr] for pr in probe_per_leaf]

    tmin = min(Z[0, 2] for Z in Zs.values()) * 0.8
    tmax = max(Z[-1, 2] for Z in Zs.values()) * 1.05

    fig = plt.figure(figsize=(13.5, 4.0))
    gs = fig.add_gridspec(
        2, 4,
        height_ratios=[14, 0.6],
        hspace=0.02, wspace=0.10,
        left=0.06, right=0.94, top=0.84, bottom=0.06,
    )
    axes_top = [fig.add_subplot(gs[0, 0])]
    for j in range(1, 4):
        axes_top.append(fig.add_subplot(gs[0, j], sharey=axes_top[0]))
    axes_bot = [fig.add_subplot(gs[1, j], sharex=axes_top[j]) for j in range(4)]

    for col, phase in enumerate(PHASES):
        ax_top = axes_top[col]
        ax_bot = axes_bot[col]
        ax_top.set_title(PHASE_TITLES[phase])
        Z = Zs[phase]
        d = dendrogram(
            Z, ax=ax_top,
            no_labels=True,
            link_color_func=color_for_phase(phase),
            above_threshold_color=UNCOLORED,
        )
        ax_top.set_ylim(tmin, tmax)
        ax_top.set_xticks([])
        ax_top.set_xlabel("")
        if col == 0:
            ax_top.set_ylabel(r"merge height  $\hat{D}(\tau)$")
        ax_top.spines["top"].set_visible(False)
        ax_top.spines["right"].set_visible(False)

        # Probe-shaft barcode.
        leaf_order = d["leaves"]
        for i, L in enumerate(leaf_order):
            ax_bot.axvspan(10 * i, 10 * (i + 1),
                            color=leaf_stub_color[L], lw=0)
        ax_bot.set_xlim(0, 10 * n)
        ax_bot.set_xticks([])
        ax_bot.set_yticks([])
        for spine in ax_bot.spines.values():
            spine.set_visible(False)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    _, _, cb = imshow_colorbar_caxdivider(sm, axes_top[-1], size="4%", pad=0.10)
    cb.set_label(
        r"per-merge-node KC trace $\bar{s}_N$"
        "\n"
        r"(partition: each pair credited to its MRCA only)",
        fontsize=6,
    )

    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    fig.text(
        0.06, 0.94,
        rf"{patient}  |  {band_tex}  |  KC$_{{\lambda={lam:g}}}$ partition "
        rf"trace  |  $\tau = 1/\lambda_{{\max}}$  |  $N = {n}$  |  "
        rf"sym vlim $\pm{vmax:.2f}$",
        fontweight="bold",
    )
    fig.text(
        0.06, 0.91,
        "U-shapes: per-merge-node KC trace (RdBu_r, partition of per-pair "
        "sum).  Bottom barcode: probe-shaft identity (same color = same "
        "leaf across panels).",
        fontsize=7, color="#444",
    )

    pdf.savefig(fig)
    plt.close(fig)


def build_band_pdf(band: str, lam: float) -> Path:
    out_pdf = OUT_DIR / f"kc_partition_trace_lambda{_lam_tag(lam)}_{band}_all_patients.pdf"
    print(f"\n=== {band} λ={lam:g} → {out_pdf.name} ===")
    with PdfPages(out_pdf) as pdf:
        for pat in COHORT:
            render_patient_page(pat, band, lam, pdf)
    print(f"DONE: {out_pdf}")
    return out_pdf


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--band", default=None,
                        help="single band (default: all 6 bands)")
    parser.add_argument("--lam", type=float, default=0.0,
                        help="KC λ (0 = topology, 0.5 = mixed, 1 = heights)")
    args = parser.parse_args()
    bands = [args.band] if args.band else list(BRAIN_BANDS.keys())
    for band in bands:
        build_band_pdf(band, args.lam)


if __name__ == "__main__":
    main()
