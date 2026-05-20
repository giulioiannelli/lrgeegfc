"""Per-patient D̂(τ) matrices across the 4 phases, probe-view ticks.

One PDF per band (6 bands × 1 PDF each = 6 PDFs), 10 pages per PDF,
each page = 1×4 subplots (D̂_rPre, D̂_tLearn, D̂_tTest, D̂_rPost) on a
shared log color scale.

Output: data/reports/preprint/figure1_beta_trace/laplacian_distances/
        D_tau_<band>_all_patients.pdf for band ∈ {delta, theta, alpha,
        beta, gamma_l, gamma_h}.

D̂(τ) = 1/ρ̂(τ) at τ = 1/λ_max (LRG communication distance, ultrametric form
stored in the imcoh_abs LRG cache). Probe-view ticks (chnames mode) put one
major tick per probe family at the family midpoint, with light minor ticks
at probe boundaries — matches `fc_templates / row_per_phase` style.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.spatial.distance import squareform

from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.visuals.fc_templates import (
    _apply_factored_sci_format,
    draw_probe_outlines,
    plot_fc_adjacency,
)
from lrg_eegfc.visuals.styles import use_lrg_style
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
BANDS = list(BRAIN_BANDS.keys())
FC = "imcoh_abs"

ROOT = Path(__file__).resolve().parents[3]
RAW_ROOT = ROOT / "data/raw/stereoeeg_patients"
OUT_DIR = ROOT / "data/reports/preprint/figure1_beta_trace/laplacian_distances"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def load_channel_labels(patient: str) -> list[str]:
    csv_path = RAW_ROOT / patient / "channel_labels.csv"
    return pd.read_csv(csv_path)["label"].astype(str).tolist()


def load_D_per_phase(patient: str, band: str) -> dict[str, np.ndarray | None]:
    out: dict[str, np.ndarray | None] = {}
    for phase in PHASES:
        try:
            lrg = load_lrg_result(patient, phase, band, FC)
        except Exception as err:
            print(f"  [{patient} {band} {phase}] missing: {err}")
            out[phase] = None
            continue
        out[phase] = squareform(lrg.ultrametric_matrix)
    return out


def _band_glyph(band: str) -> str:
    tex = BRAIN_BAND_TEX_DICT.get(band, band)
    return tex if isinstance(tex, str) else str(tex)


def page_for_patient(patient: str, band: str, pdf: PdfPages) -> None:
    D_per_phase = load_D_per_phase(patient, band)
    available = [(p, D) for p, D in D_per_phase.items() if D is not None]
    if not available:
        print(f"[{patient} {band}] no data, skipping")
        return

    n = available[0][1].shape[0]
    try:
        labels = load_channel_labels(patient)
        if len(labels) != n:
            print(f"  [{patient} {band}] label/n mismatch ({len(labels)} vs {n}); generic ticks")
            labels = None
            tick_mode = "generic"
        else:
            tick_mode = "chnames"
    except FileNotFoundError:
        print(f"  [{patient} {band}] no channel_labels.csv; generic ticks")
        labels = None
        tick_mode = "generic"

    flat = np.concatenate([
        np.asarray(D, dtype=float)[~np.eye(D.shape[0], dtype=bool)].ravel()
        for _, D in available
    ])
    flat = flat[np.isfinite(flat) & (flat > 0)]
    vmin = float(np.percentile(flat, 1))
    vmax = float(np.percentile(flat, 99))

    width_ratios = [1.0] * 3 + [1.06]
    fig, axes = plt.subplots(
        1, 4, figsize=(13.5, 3.6),
        sharey=True, gridspec_kw={"width_ratios": width_ratios},
    )
    fig.subplots_adjust(left=0.06, right=0.94, top=0.86, bottom=0.16, wspace=0.10)

    last_im = None
    for col, (ax, phase) in enumerate(zip(axes, PHASES)):
        D = D_per_phase[phase]
        ax.set_title(PHASE_TITLES[phase])
        if D is None:
            ax.text(0.5, 0.5, "missing", ha="center", va="center",
                    transform=ax.transAxes, color="#888")
            ax.set_xticks([]); ax.set_yticks([])
            continue
        _, _, im = plot_fc_adjacency(
            D, ax=ax,
            vmin=vmin, vmax=vmax,
            zero_diagonal=True,
            log_scale=True,
            tick_labels=tick_mode,
            channel_labels=labels,
            colorbar=False,
            is_first_col=(col == 0),
            is_last_row=True,
        )
        if labels is not None:
            draw_probe_outlines(ax, labels, lw=0.6, alpha=0.8)
        last_im = im

    if last_im is not None:
        _, _, clb = imshow_colorbar_caxdivider(last_im, axes[-1], size="4%", pad=0.06)
        clb.set_label(r"$\hat{D}(\tau) = 1/\hat{\rho}(\tau)$")
        _apply_factored_sci_format(clb, axis_orientation="vertical")

    fig.text(
        0.06, 0.94,
        rf"{patient}  |  {_band_glyph(band)}  |  "
        rf"$\tau = 1/\lambda_{{\max}}$  |  $N = {n}$",
        fontweight="bold",
    )

    pdf.savefig(fig)
    plt.close(fig)


def build_band_pdf(band: str) -> Path:
    out_pdf = OUT_DIR / f"D_tau_{band}_all_patients.pdf"
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
