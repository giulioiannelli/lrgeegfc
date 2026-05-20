"""KC λ=0 per-subtree trace coloring on the 4-phase dendrogram (batch).

Per-pair quantity (KC λ=0):
    m_phase[i,j] = #edges from root to MRCA(i,j)  in phase's dendrogram
                 = ``kc_vectors(Z)[0]`` reshaped into a symmetric matrix

Per-pair 3-phase trace strength (tree-invariant scalar — same value used
on all four panels):
    s(i,j) = ½(|Δm_{pre,test}| + |Δm_{pre,post}|) − |Δm_{test,post}|

Per-subtree score (for each internal node S of *each* phase's dendrogram):
    s̄(S) = mean of s(i,j) over leaf pairs (i,j) inside S

Each U-shape link is coloured by s̄ via a divergent RdBu_r centred at 0,
symmetric vlim per page. Red = subtree concentrates trace-strong pairs;
blue = anti-trace pairs.

**Caveat (important for reading the figure):** s(i,j) is identical
across all four panels — what varies is which subtrees the pair falls
into. A red rPre subtree does *not* mean rPre is similar to rPost; it
means rPre happens to group trace-strong pairs into that clade. To
compare two trees pair-by-pair, use the matched-cluster figure
(``fig_preprint_dendrograms_matched_clusters.py``).

One PDF per band, one page per patient × 4 phases.

Output:
    data/reports/preprint/figure1_beta_trace/dendrograms_kc_trace/
        kc_lambda0_trace_{band}_all_patients.pdf
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
    """Upper-triangle scipy-ordered vector → symmetric (n, n) matrix."""
    M = np.zeros((n, n), dtype=float)
    M[np.triu_indices(n, k=1)] = m_vec
    return M + M.T


def per_subtree_score(
    leaves_dict: dict[int, frozenset[int]],
    s_pair: np.ndarray,
) -> dict[int, float]:
    out: dict[int, float] = {}
    for node_id, leaves in leaves_dict.items():
        if len(leaves) < 2:
            continue
        arr = np.fromiter(leaves, dtype=int)
        sub = s_pair[np.ix_(arr, arr)]
        out[node_id] = float(sub[np.triu_indices(len(arr), k=1)].mean())
    return out


def render_patient_page(patient: str, band: str, lam: float, pdf: PdfPages) -> None:
    try:
        Zs = {p: load_lrg_result(patient, p, band, FC).linkage_matrix for p in PHASES}
    except Exception as err:
        print(f"  [{patient} {band}] skip ({err})")
        return
    n = Zs[PHASES[0]].shape[0] + 1

    # KC λ-mixed per-pair vectors. m = MRCA-edge-count (topology),
    # M = MRCA merge height. Normalize each by per-patient pooled max
    # across all 4 phases so the two components blend on comparable
    # scales (Kendall & Colijn 2016, kc_distance with normalize=True).
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

    leaves_per_phase = {p: leaves_under_each_node(Z) for p, Z in Zs.items()}
    scores_per_phase = {
        p: per_subtree_score(leaves_per_phase[p], s) for p in PHASES
    }

    all_scores = np.fromiter(
        (v for d in scores_per_phase.values() for v in d.values()),
        dtype=float,
    )
    vmax = float(np.max(np.abs(all_scores))) if all_scores.size else 1.0
    vmin = -vmax

    print(f"  [{patient} {band}] λ={lam:.2f}  per-pair s: min={s.min():+.3f} "
          f"max={s.max():+.3f}  per-subtree vmax (sym)={vmax:.3f}")

    cmap = plt.get_cmap("RdBu_r")
    norm = Normalize(vmin=vmin, vmax=vmax)

    def color_for_phase(phase: str):
        scores = scores_per_phase[phase]

        def cf(node_id: int) -> str:
            score = scores.get(node_id)
            if score is None or not np.isfinite(score):
                return "#888888"
            return mcolors.to_hex(cmap(norm(score)))

        return cf

    # Cross-panel leaf identity anchor: probe (shaft) of each electrode.
    # All leaves on the same sEEG shaft get the same color; ~8-15 probes
    # per patient means tab20 covers them cleanly with no color reuse.
    channel_labels = load_channel_labels(patient)
    probe_per_leaf = extract_probe_labels(channel_labels)
    unique_probes = sorted(set(probe_per_leaf))
    palette = plt.get_cmap("tab20")
    probe_to_color = {
        pr: mcolors.to_hex(palette(i % palette.N))
        for i, pr in enumerate(unique_probes)
    }
    leaf_anchor_color = [probe_to_color[pr] for pr in probe_per_leaf]
    print(f"  [{patient} {band}] probes ({len(unique_probes)}): {unique_probes}")

    tmin = min(Z[0, 2] for Z in Zs.values()) * 0.8
    tmax = max(Z[-1, 2] for Z in Zs.values()) * 1.05

    fig, axes = plt.subplots(
        1, 4, figsize=(13.5, 3.8),
        sharey=True,
        gridspec_kw={"width_ratios": [1, 1, 1, 1.10]},
    )
    fig.subplots_adjust(left=0.06, right=0.94, top=0.84, bottom=0.06, wspace=0.10)

    for col, (ax, phase) in enumerate(zip(axes, PHASES)):
        ax.set_title(PHASE_TITLES[phase])
        Z = Zs[phase]
        d = dendrogram(
            Z, ax=ax,
            no_labels=True,
            link_color_func=color_for_phase(phase),
            above_threshold_color="#888888",
        )
        # Per-leaf identity stub at the bottom of each leaf's leg.
        # scipy places leaves at x = 5, 15, …, 5 + 10·(n-1).
        leaf_order = d["leaves"]
        leaf_xs = [5 + 10 * i for i in range(len(leaf_order))]
        stub_top = tmin + 0.025 * (tmax - tmin)
        for x, L in zip(leaf_xs, leaf_order):
            ax.plot(
                [x, x], [tmin, stub_top],
                color=leaf_anchor_color[L], linewidth=2.0,
                solid_capstyle="butt", zorder=20,
            )
        ax.set_ylim(tmin, tmax)
        ax.set_xticks([])
        ax.set_xlabel("")
        if col == 0:
            ax.set_ylabel(r"merge height  $\hat{D}(\tau)$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    _, _, cb = imshow_colorbar_caxdivider(sm, axes[-1], size="4%", pad=0.10)
    cb.set_label(
        r"subtree KC$_{\lambda=0}$ trace score $\bar{s}$"
        "\n"
        r"$\frac{1}{2}(|\Delta m_{\mathrm{pre,test}}|+|\Delta m_{\mathrm{pre,post}}|)"
        r"\,-\,|\Delta m_{\mathrm{test,post}}|$",
        fontsize=7,
    )

    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    fig.text(
        0.06, 0.94,
        rf"{patient}  |  {band_tex}  |  KC$_{{\lambda={lam:g}}}$ trace coloring  "
        rf"|  $\tau = 1/\lambda_{{\max}}$  |  $N = {n}$  |  "
        rf"sym vlim $\pm{vmax:.2f}$",
        fontweight="bold",
    )
    fig.text(
        0.06, 0.91,
        r"U-shapes: subtree-mean KC trace $\bar{s}$ (RdBu_r).  "
        "Bottom-of-leg stubs: probe (shaft) identity — "
        "leaves on the same sEEG shaft share a stub color across all panels.",
        fontsize=7, color="#444",
    )

    pdf.savefig(fig)
    plt.close(fig)


def _lam_tag(lam: float) -> str:
    return f"{lam:g}".replace(".", "p")


def build_band_pdf(band: str, lam: float) -> Path:
    out_pdf = OUT_DIR / f"kc_lambda{_lam_tag(lam)}_trace_{band}_all_patients.pdf"
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
                        help="KC λ (default 0 = topology; 0.5 = mixed; "
                             "1 = pure heights)")
    args = parser.parse_args()
    bands = [args.band] if args.band else list(BRAIN_BANDS.keys())
    for band in bands:
        build_band_pdf(band, args.lam)


if __name__ == "__main__":
    main()
