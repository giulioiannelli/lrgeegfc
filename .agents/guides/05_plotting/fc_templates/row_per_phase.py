"""row_per_phase — one patient × one band × the four phases.

Thin wrapper around ``lrg_eegfc.visuals.plot_fc_adjacency_row``.
Defaults to a shared colour scale across the row + one colorbar at
the right (phases share a band, so amplitudes are directly
comparable).

Usage
-----
    conda activate lapbrain
    python .agents/guides/05_plotting/fc_templates/row_per_phase.py \\
        --patient Pat_06 --band beta --fc-method imcoh_abs \\
        --tick-labels generic
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from lrg_eegfc.config.const import PHASE_LABELS
from lrg_eegfc.config.paths import FIGURES_ROOT, SEEG_DATAPATH
from lrg_eegfc.visuals.fc_templates import plot_fc_adjacency_row
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
    parser.add_argument("--patient", default="Pat_06")
    parser.add_argument("--band", default="beta")
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
        "--out-dir",
        default=str(FIGURES_ROOT / "fc_templates" / "row_per_phase"),
    )
    args = parser.parse_args()

    matrices = []
    for phase in PHASE_LABELS:
        M = load_fc_matrix(
            patient=args.patient, phase=phase,
            band=args.band, fc_method=args.fc_method,
        )
        if M is None:
            raise SystemExit(
                f"No cached FC for {args.patient}/{args.band}/"
                f"{phase}/{args.fc_method}.  Compute it first."
            )
        matrices.append(M)

    chnames = (
        _load_channel_labels(args.patient)
        if args.tick_labels == "chnames" else None
    )

    fig, _axes = plot_fc_adjacency_row(
        matrices,
        titles=[PHASE_TEX[p] for p in PHASE_LABELS],
        fc_method=args.fc_method,
        band=args.band,
        tick_labels=args.tick_labels,
        channel_labels=chnames,
        log_scale=args.log_scale,
        shared_scale=True,
    )

    if args.watermark:
        add_provenance_footer(
            fig,
            f"row_per_phase · {args.patient} · {args.band} · "
            f"{args.fc_method}",
        )

    scale_tag = "log" if args.log_scale else "lin"
    out = (
        Path(args.out_dir)
        / f"{args.patient}_{args.band}_{args.fc_method}_"
          f"{args.tick_labels}_{scale_tag}.pdf"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
