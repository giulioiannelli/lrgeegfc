"""mosaic_patient_phase — patients × phases grid at a fixed band.

Rows iterate over a patient subset, columns iterate over the four
phases.  Default colorbar mode is **per-row**: each patient has its
own scale, since cohort-wide amplitude varies.  Use --colorbar-mode
shared if you want to compare absolute amplitudes across patients.

Usage
-----
    conda activate lapbrain
    python .agents/guides/05_plotting/fc_templates/mosaic_patient_phase.py \\
        --patients Pat_05,Pat_06,Pat_08,Pat_13 --band beta
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from lrg_eegfc.config.const import PHASE_LABELS
from lrg_eegfc.config.paths import FIGURES_ROOT, SEEG_DATAPATH
from lrg_eegfc.visuals.fc_templates import plot_fc_adjacency_grid
from lrg_eegfc.visuals.layout import add_provenance_footer
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.fc import load_fc_matrix


PHASE_TEX = {
    "rest_pre": r"rest$_{\mathrm{pre}}$",
    "task_learn": r"task$_{\mathrm{learn}}$",
    "task_test": r"task$_{\mathrm{test}}$",
    "rest_post": r"rest$_{\mathrm{post}}$",
}

DEFAULT_PATIENTS = "Pat_05,Pat_06,Pat_08,Pat_13"


def _load_channel_labels(patient: str) -> list[str]:
    csv = SEEG_DATAPATH / patient / "channel_labels.csv"
    return pd.read_csv(csv)["label"].astype(str).tolist()


def main() -> None:
    use_lrg_style()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--patients", default=DEFAULT_PATIENTS,
        help="Comma-separated patient list. Default: %(default)s",
    )
    parser.add_argument("--band", default="beta")
    parser.add_argument(
        "--fc-method", default="imcoh_abs",
        choices=["corr", "msc", "imcoh_abs", "imcoh_sq"],
    )
    parser.add_argument(
        "--tick-labels", default="chnames",
        choices=["generic", "index", "chnames"],
        help="Default: chnames (per-row probe-family labels).",
    )
    parser.add_argument(
        "--colorbar-label", default=None,
        help="Override the auto imcoh/band colorbar label with a "
             "method-neutral string (e.g. '$\\mathrm{FC}_{ij}$'); pass '' "
             "to drop it. Default: canonical <|imcoh|>_band.",
    )
    parser.add_argument("--watermark", action="store_true")
    parser.add_argument(
        "--log-scale", action=argparse.BooleanOptionalAction, default=True,
        help="Logarithmic colour norm; default on (use --no-log-scale "
             "for a linear scale).",
    )
    parser.add_argument(
        "--colorbar-mode", default="per_row",
        choices=["shared", "per_row", "per_col"],
        help="Default: per_row (each patient its own scale).",
    )
    parser.add_argument(
        "--out-dir",
        default=str(FIGURES_ROOT / "fc_templates" / "mosaic_patient_phase"),
    )
    args = parser.parse_args()

    patients = [p.strip() for p in args.patients.split(",") if p.strip()]
    if not patients:
        raise SystemExit("--patients must list at least one patient")

    matrices = []
    for patient in patients:
        row = []
        for phase in PHASE_LABELS:
            M = load_fc_matrix(
                patient=patient, phase=phase,
                band=args.band, fc_method=args.fc_method,
            )
            if M is None:
                raise SystemExit(
                    f"No cached FC for {patient}/{args.band}/{phase}/"
                    f"{args.fc_method}.  Compute it first."
                )
            row.append(M)
        matrices.append(row)

    chnames_per_row = (
        [_load_channel_labels(p) for p in patients]
        if args.tick_labels == "chnames" else None
    )

    fig, _axes = plot_fc_adjacency_grid(
        matrices,
        row_titles=patients,
        col_titles=[PHASE_TEX[p] for p in PHASE_LABELS],
        fc_method=args.fc_method,
        band=args.band,  # all rows share this band → every cb labelled with it
        colorbar_label=args.colorbar_label,
        tick_labels=args.tick_labels,
        channel_labels_per_row=chnames_per_row,
        log_scale=args.log_scale,
        colorbar_mode=args.colorbar_mode,
    )

    if args.watermark:
        add_provenance_footer(
            fig, f"mosaic_patient_phase · n={len(patients)} · "
                 f"{args.band} · {args.fc_method}",
        )

    scale_tag = "log" if args.log_scale else "lin"
    pat_tag = f"n{len(patients)}" if len(patients) > 1 else patients[0]
    out = (
        Path(args.out_dir)
        / f"{pat_tag}_{args.band}_{args.fc_method}_"
          f"{args.tick_labels}_{scale_tag}_{args.colorbar_mode}.pdf"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
