"""mosaic_band_phase — phases × bands grid for one patient (horizontal).

A 4 × 6 mosaic of FC adjacency matrices: rows iterate over phases
(rest$_{\\mathrm{pre}}$, task$_{\\mathrm{learn}}$, task$_{\\mathrm{test}}$,
rest$_{\\mathrm{post}}$), columns iterate over bands (δ, θ, α, β,
γ$_l$, γ$_h$).  Default colorbar mode is **per-col** so each band
column has its own scale (γ would be crushed otherwise); horizontal
cbs at the bottom carry the band glyph.

Usage
-----
    conda activate lapbrain
    python .agents/guides/05_plotting/fc_templates/mosaic_band_phase.py \\
        --patient Pat_05 --fc-method imcoh_abs
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT, BRAIN_BANDS_NAMES, PHASE_LABELS
from lrg_eegfc.config.paths import FIGURES_ROOT, SEEG_DATAPATH
from lrg_eegfc.visuals.fc_templates import plot_fc_adjacency_grid
from lrg_eegfc.visuals.layout import add_provenance_footer
from lrg_eegfc.workflow.fc import load_fc_matrix


PHASE_TEX = {
    "rest_pre": r"rest$_{\mathrm{pre}}$",
    "task_learn": r"task$_{\mathrm{learn}}$",
    "task_test": r"task$_{\mathrm{test}}$",
    "rest_post": r"rest$_{\mathrm{post}}$",
}


def _load_channel_labels(patient: str) -> list[str]:
    csv = SEEG_DATAPATH / patient / "channel_labels.csv"
    return pd.read_csv(csv)["label"].astype(str).tolist()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--patient", default="Pat_05")
    parser.add_argument(
        "--fc-method", default="imcoh_abs",
        choices=["corr", "msc", "imcoh_abs", "imcoh_sq"],
    )
    parser.add_argument(
        "--tick-labels", default="chnames",
        choices=["generic", "index", "chnames"],
        help="Default: chnames (probe-family midpoint labels).",
    )
    parser.add_argument("--watermark", action="store_true")
    parser.add_argument(
        "--log-scale", action=argparse.BooleanOptionalAction, default=True,
        help="Logarithmic colour norm; default on (use --no-log-scale "
             "for a linear scale).",
    )
    parser.add_argument(
        "--colorbar-mode", default="per_col",
        choices=["shared", "per_row", "per_col"],
        help="Colour-scale policy.  Default: per_col (each band column "
             "its own scale, with horizontal cbs at the bottom).  Use "
             "--colorbar-mode shared to compare absolute amplitudes "
             "across bands on a single global scale.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(FIGURES_ROOT / "fc_templates" / "mosaic_band_phase"),
    )
    args = parser.parse_args()

    # Horizontal layout: rows = phases (4), cols = bands (6).
    matrices = []
    for phase in PHASE_LABELS:
        row = []
        for band in BRAIN_BANDS_NAMES:
            M = load_fc_matrix(
                patient=args.patient, phase=phase,
                band=band, fc_method=args.fc_method,
            )
            if M is None:
                raise SystemExit(
                    f"No cached FC for {args.patient}/{band}/{phase}/"
                    f"{args.fc_method}.  Compute it first."
                )
            row.append(M)
        matrices.append(row)

    chnames = (
        _load_channel_labels(args.patient)
        if args.tick_labels == "chnames" else None
    )

    fig, _axes = plot_fc_adjacency_grid(
        matrices,
        row_titles=[PHASE_TEX[p] for p in PHASE_LABELS],
        col_titles=[BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
        fc_method=args.fc_method,
        col_bands=BRAIN_BANDS_NAMES,
        tick_labels=args.tick_labels,
        channel_labels=chnames,
        log_scale=args.log_scale,
        colorbar_mode=args.colorbar_mode,
    )

    if args.watermark:
        add_provenance_footer(
            fig, f"mosaic_band_phase · {args.patient} · {args.fc_method}",
        )

    scale_tag = "log" if args.log_scale else "lin"
    out = (
        Path(args.out_dir)
        / f"{args.patient}_{args.fc_method}_"
          f"{args.tick_labels}_{scale_tag}_{args.colorbar_mode}.pdf"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
