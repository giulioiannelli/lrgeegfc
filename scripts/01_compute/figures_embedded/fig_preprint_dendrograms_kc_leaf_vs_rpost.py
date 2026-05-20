"""Per-leaf KC contribution to rPost, as a barcode below each phase's
dendrogram. Visualises KC tree-similarity directly on the dendrograms.

For each leaf L and each phase X, define:

    kc_leaf(L, X) = Σ_{j ≠ L}  |v_X^λ(L, j) − v_rPost^λ(L, j)|

where `v_X^λ(i,j) = (1-λ)·m_X(i,j)/m_max + λ·M_X(i,j)/M_max` is the
Kendall-Colijn λ-mixed pair vector (KC `(m, M)` pair, pooled
normalization à la `kc_distance(normalize=True)`).

Properties:
- **Partition** — `Σ_L kc_leaf(L, X) = 2 · KC^{L1}(X, rPost)` (every
  pair `|v_X − v_rPost|` is counted once for each endpoint). So the
  barcode total *is* the KC distance.
- **rPost** — `kc_leaf(L, rPost) = 0` for every leaf, so its barcode
  is uniformly zero. This is the visual anchor for "perfectly similar
  to rPost".

Visualisation:
- Four dendrograms side-by-side (rPre, tLearn, tTest, rPost), black
  links (topology only).
- Colorbar-style barcode strip below each, one cell per leaf in that
  panel's display order, colored by `kc_leaf(L, X)` using a sequential
  cmap with shared `vmin = 0`, `vmax = max across the 3 non-anchor
  phases for this page` (so the 4 barcodes are visually comparable).
- Per-panel headline: total `KC(X, rPost)`.

Reading: the phase whose barcode is **overall coolest** (closest to
rPost's uniform zero) is the most KC-similar to rPost. Hot stripes
mark the leaves driving the disagreement.

Output:
    data/reports/preprint/figure1_beta_trace/dendrograms_kc_leaf_vs_rpost/
        kc_leaf_vs_rpost_lambda{lam}_{band}_all_patients.pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.metrics.tree_distance import kc_vectors
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrgsglib.plotlib import imshow_colorbar_caxdivider

use_lrg_style()


COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
NON_ANCHOR = ["rest_pre", "task_learn", "task_test"]
ANCHOR_PHASE = "rest_post"
PHASE_TITLES = {
    "rest_pre": r"$\mathrm{rPre}$",
    "task_learn": r"$\mathrm{tLearn}$",
    "task_test": r"$\mathrm{tTest}$",
    "rest_post": r"$\mathrm{rPost}$",
}
FC = "imcoh_abs"

ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "data/reports/preprint/figure1_beta_trace/dendrograms_kc_leaf_vs_rpost"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _lam_tag(lam: float) -> str:
    return f"{lam:g}".replace(".", "p")


def m_vec_to_pair_matrix(m_vec: np.ndarray, n: int) -> np.ndarray:
    M = np.zeros((n, n), dtype=float)
    M[np.triu_indices(n, k=1)] = m_vec
    return M + M.T


def render_patient_page(patient: str, band: str, lam: float, pdf: PdfPages) -> None:
    try:
        Zs = {p: load_lrg_result(patient, p, band, FC).linkage_matrix for p in PHASES}
    except Exception as err:
        print(f"  [{patient} {band}] skip ({err})")
        return
    n = Zs[PHASES[0]].shape[0] + 1

    m = {p: kc_vectors(Z)[0].astype(float) for p, Z in Zs.items()}
    M = {p: kc_vectors(Z)[1].astype(float) for p, Z in Zs.items()}
    mmax = max(v.max() for v in m.values()) or 1.0
    Mmax = max(v.max() for v in M.values()) or 1e-12
    v_pair = {
        p: (1.0 - lam) * m[p] / mmax + lam * M[p] / Mmax
        for p in PHASES
    }
    Vmat = {p: m_vec_to_pair_matrix(v_pair[p], n) for p in PHASES}
    V_post = Vmat[ANCHOR_PHASE]

    kc_per_leaf: dict[str, np.ndarray] = {}
    total_kc: dict[str, float] = {}
    for p in PHASES:
        diff = np.abs(Vmat[p] - V_post)
        np.fill_diagonal(diff, 0.0)
        kc_per_leaf[p] = diff.sum(axis=1)
        total_kc[p] = float(kc_per_leaf[p].sum() / 2.0)

    vmax = max(kc_per_leaf[p].max() for p in NON_ANCHOR)
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0
    norm = Normalize(vmin=0.0, vmax=float(vmax))
    cmap = plt.get_cmap("inferno")

    tmin = min(Z[0, 2] for Z in Zs.values()) * 0.8
    tmax = max(Z[-1, 2] for Z in Zs.values()) * 1.05

    fig = plt.figure(figsize=(13.5, 4.2))
    gs = fig.add_gridspec(
        2, 4,
        height_ratios=[12, 0.8],
        hspace=0.04, wspace=0.10,
        left=0.06, right=0.93, top=0.80, bottom=0.06,
    )
    axes_top = [fig.add_subplot(gs[0, 0])]
    for j in range(1, 4):
        axes_top.append(fig.add_subplot(gs[0, j], sharey=axes_top[0]))
    axes_bot = [fig.add_subplot(gs[1, j], sharex=axes_top[j]) for j in range(4)]

    for col, phase in enumerate(PHASES):
        ax_top = axes_top[col]
        ax_bot = axes_bot[col]
        ax_top.set_title(
            f"{PHASE_TITLES[phase]}"
            f"\n$KC(\\cdot,\\,\\mathrm{{rPost}}) = {total_kc[phase]:.1f}$",
            fontsize=8,
        )
        Z = Zs[phase]
        d = dendrogram(
            Z, ax=ax_top,
            no_labels=True,
            link_color_func=lambda _k: "#222222",
            above_threshold_color="#222222",
        )
        ax_top.set_ylim(tmin, tmax)
        ax_top.set_xticks([])
        ax_top.set_xlabel("")
        if col == 0:
            ax_top.set_ylabel(r"merge height  $\hat{D}(\tau)$")
        ax_top.spines["top"].set_visible(False)
        ax_top.spines["right"].set_visible(False)

        leaf_order = d["leaves"]
        for i, L in enumerate(leaf_order):
            color = cmap(norm(float(kc_per_leaf[phase][L])))
            ax_bot.axvspan(10 * i, 10 * (i + 1), color=color, lw=0)
        ax_bot.set_xlim(0, 10 * n)
        ax_bot.set_xticks([])
        ax_bot.set_yticks([])
        for spine in ax_bot.spines.values():
            spine.set_visible(False)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    _, _, cb = imshow_colorbar_caxdivider(sm, axes_top[-1], size="4%", pad=0.10)
    cb.set_label(
        r"$kc_{\mathrm{leaf}}(L, X) = \sum_{j} |v^\lambda_X(L,j) - v^\lambda_{\mathrm{rPost}}(L,j)|$",
        fontsize=6,
    )

    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    best = min(NON_ANCHOR, key=lambda p: total_kc[p])
    fig.text(
        0.06, 0.94,
        rf"{patient}  |  {band_tex}  |  per-leaf KC vs rPost  |  "
        rf"$\lambda = {lam:g}$  |  $N = {n}$  |  closest to rPost: "
        rf"{PHASE_TITLES[best]} ($KC = {total_kc[best]:.1f}$)",
        fontweight="bold",
    )
    fig.text(
        0.06, 0.905,
        "Each leaf cell = its share of $KC(X, \\mathrm{rPost})$. "
        "rPost's barcode is uniformly 0 by construction (visual anchor for "
        "'matches rPost'). Phase with overall coolest barcode = most "
        "KC-similar to rPost.",
        fontsize=7, color="#444",
    )

    print(f"  [{patient} {band}] λ={lam:g}  "
          f"KC(rPre,post)={total_kc['rest_pre']:.1f}  "
          f"KC(tLearn,post)={total_kc['task_learn']:.1f}  "
          f"KC(tTest,post)={total_kc['task_test']:.1f}  →  closest={best}")

    pdf.savefig(fig)
    plt.close(fig)


def build_band_pdf(band: str, lam: float) -> Path:
    out_pdf = OUT_DIR / f"kc_leaf_vs_rpost_lambda{_lam_tag(lam)}_{band}_all_patients.pdf"
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
