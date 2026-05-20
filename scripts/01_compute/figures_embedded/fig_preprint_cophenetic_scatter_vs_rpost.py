"""Cophenetic scatter: per-pair MRCA-depth (KC m-vector) of each phase
plotted against rPost. Three scatter panels per patient page, one for
each non-anchor phase. The tighter the points hug the diagonal, the
more topologically similar that tree is to rPost in the KC sense.

Per-pair KC λ-mixed value (Kendall & Colijn 2016, pooled normalisation):

    v_X^λ(i,j) = (1-λ) * m_X(i,j) / m_max + λ * M_X(i,j) / M_max

For each phase X in {rPre, tLearn, tTest} the panel shows the
hexbin density of {(v_X(i,j), v_post(i,j)) : i<j}, the identity line
y=x, Pearson r, Spearman ρ, and L1 distance (= KC trace contribution
to rPost at this λ). Phase with highest r (= tightest scatter) is the
most KC-similar to rPost; headline per page names it.

Visual reads at a glance because the *shape* of the cloud encodes the
agreement: a thin line along the diagonal = perfect agreement; a fat
oval = random agreement. No clade matching, no thresholds, no
untangling heuristics — direct cophenetic comparison.

Output:
    data/reports/preprint/figure1_beta_trace/cophenetic_scatters_vs_rpost/
        cophenetic_scatter_lambda{lam}_{band}_all_patients.pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import spearmanr

from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.metrics.tree_distance import kc_vectors
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
OUT_DIR = ROOT / "data/reports/preprint/figure1_beta_trace/cophenetic_scatters_vs_rpost"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _lam_tag(lam: float) -> str:
    return f"{lam:g}".replace(".", "p")


def compute_v_lambda(Zs: dict[str, np.ndarray], lam: float) -> dict[str, np.ndarray]:
    """Return per-pair KC λ-mixed vector per phase (pooled normalization)."""
    m = {p: kc_vectors(Z)[0].astype(float) for p, Z in Zs.items()}
    M = {p: kc_vectors(Z)[1].astype(float) for p, Z in Zs.items()}
    mmax = max(v.max() for v in m.values()) or 1.0
    Mmax = max(v.max() for v in M.values()) or 1e-12
    return {p: (1.0 - lam) * m[p] / mmax + lam * M[p] / Mmax for p in m}


def render_patient_page(
    patient: str, band: str, lam: float, pdf: PdfPages
) -> dict[str, tuple[float, float, float]]:
    try:
        Zs = {p: load_lrg_result(patient, p, band, FC).linkage_matrix
              for p in PHASES_LEFT + [ANCHOR_PHASE]}
    except Exception as err:
        print(f"  [{patient} {band}] skip ({err})")
        return {}

    v = compute_v_lambda(Zs, lam)
    v_post = v[ANCHOR_PHASE]

    # Global axis range so all 3 panels are visually comparable.
    all_vals = np.concatenate([v[p] for p in PHASES_LEFT + [ANCHOR_PHASE]])
    lo = float(all_vals.min())
    hi = float(all_vals.max())
    pad = 0.04 * (hi - lo) if hi > lo else 1.0

    fig, axes = plt.subplots(
        1, 3, figsize=(13.5, 4.5),
        sharex=True, sharey=True,
    )
    fig.subplots_adjust(left=0.06, right=0.99, top=0.82, bottom=0.10, wspace=0.08)

    # Jitter scale: a fraction of the smallest non-zero spacing in v_post.
    # For λ=0 the v values are integer/m_max so spacing is 1/m_max → jitter
    # ≈ 0.4·spacing spreads integer pile-ups visibly. For λ=0.5 the
    # values are mostly distinct floats so jitter is essentially a tiny
    # nudge.
    unique_v = np.unique(v_post)
    if unique_v.size >= 2:
        diffs = np.diff(unique_v)
        spacing = float(diffs[diffs > 0].min())
    else:
        spacing = 0.01
    jitter_sigma = 0.4 * spacing
    rng = np.random.default_rng(seed=42)

    metrics: dict[str, tuple[float, float, float]] = {}
    for ax, phase in zip(axes, PHASES_LEFT):
        v_phase = v[phase]
        x = v_phase + rng.normal(0.0, jitter_sigma, size=v_phase.size)
        y = v_post + rng.normal(0.0, jitter_sigma, size=v_post.size)
        ax.scatter(
            x, y,
            s=4, c="#1f77b4",
            alpha=0.18, linewidths=0,
            rasterized=False,
        )
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
                color="#888888", linestyle="--", linewidth=0.8, zorder=10)
        r = float(np.corrcoef(v_phase, v_post)[0, 1])
        rho = float(spearmanr(v_phase, v_post).statistic)
        l1 = float(np.abs(v_phase - v_post).sum())
        metrics[phase] = (r, rho, l1)
        ax.set_title(
            f"{PHASE_TITLES[phase]}"
            f"\n$r = {r:.3f}$,  $\\rho = {rho:.3f}$,  $L_1 = {l1:.1f}$",
            fontsize=8,
        )
        ax.set_xlabel(rf"$v^{{\lambda={lam:g}}}_{{\mathrm{{{phase.replace('_','')}}}}}(i,j)$")
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_aspect("equal")

    axes[0].set_ylabel(rf"$v^{{\lambda={lam:g}}}_{{\mathrm{{rPost}}}}(i,j)$")

    best_r = max(metrics, key=lambda p: metrics[p][0])
    best_l1 = min(metrics, key=lambda p: metrics[p][2])
    n = Zs[ANCHOR_PHASE].shape[0] + 1
    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    fig.text(
        0.06, 0.94,
        rf"{patient}  |  {band_tex}  |  cophenetic vs rPost  |  "
        rf"$\lambda = {lam:g}$  |  $N = {n}$  |  pairs $= {len(v_post)}$  |  "
        rf"closest to rPost — by $r$: {PHASE_TITLES[best_r]},  "
        rf"by $L_1$: {PHASE_TITLES[best_l1]}",
        fontweight="bold",
    )
    fig.text(
        0.06, 0.91,
        "Each point = one leaf pair $(i,j)$. Hexbin density on log scale; "
        "dashed line = identity. Tighter cloud around the diagonal $\\Rightarrow$ "
        "more KC-similar to rPost (smaller $L_1$, higher $r$).",
        fontsize=7, color="#444",
    )

    print(f"  [{patient} {band}] λ={lam:g}  "
          f"r: rPre={metrics['rest_pre'][0]:.3f}  "
          f"tLearn={metrics['task_learn'][0]:.3f}  "
          f"tTest={metrics['task_test'][0]:.3f}  →  best_r={best_r}  "
          f"|  L1: rPre={metrics['rest_pre'][2]:.1f}  "
          f"tLearn={metrics['task_learn'][2]:.1f}  "
          f"tTest={metrics['task_test'][2]:.1f}  →  best_L1={best_l1}")

    pdf.savefig(fig)
    plt.close(fig)
    return metrics


def build_band_pdf(band: str, lam: float) -> Path:
    out_pdf = OUT_DIR / f"cophenetic_scatter_lambda{_lam_tag(lam)}_{band}_all_patients.pdf"
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
