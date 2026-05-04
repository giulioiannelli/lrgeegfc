"""single_adjacency — canonical FC adjacency-matrix figure (one panel).

Thin wrapper around ``lrg_eegfc.visuals.plot_fc_adjacency``.  Loads
one (patient, band, phase, fc_method) FC matrix from cache, renders
it as a vector PDF, no watermark by default.

Usage
-----
    conda activate lapbrain
    python .agents/guides/05_plotting/fc_templates/single_adjacency.py \\
        --patient Pat_02 --band beta --phase rest_pre \\
        --fc-method imcoh_abs --tick-labels generic

``--tick-labels`` ∈ {``generic`` (default, math ``$i$``/``$j$``),
``index`` (numeric ticks), ``chnames`` (probe-family midpoints)}.
``--watermark`` opts in to the grey provenance footer (off by
default).

See sibling ``single_adjacency.md`` for the style sheet.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import FIGURES_ROOT, SEEG_DATAPATH
from lrg_eegfc.visuals.fc_templates import plot_fc_adjacency
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
    parser.add_argument("--patient", default="Pat_02")
    parser.add_argument("--band", default="beta")
    parser.add_argument("--phase", default="rest_pre")
    parser.add_argument(
        "--fc-method", default="imcoh_abs",
        choices=["corr", "msc", "imcoh_abs", "imcoh_sq"],
    )
    parser.add_argument(
        "--tick-labels", default="chnames",
        choices=["generic", "index", "chnames"],
        help="Axis-tick mode: generic (i/j), index (numeric), chnames "
             "(probe-family midpoints).  Default: chnames.",
    )
    parser.add_argument(
        "--watermark", action="store_true",
        help="Add a small grey provenance string bottom-right (off by default).",
    )
    parser.add_argument(
        "--log-scale", action=argparse.BooleanOptionalAction, default=True,
        help="Logarithmic colour norm (LogNorm); default on, "
             "use --no-log-scale for a linear scale.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(FIGURES_ROOT / "fc_templates" / "single_adjacency"),
    )
    args = parser.parse_args()

    M = load_fc_matrix(
        patient=args.patient, phase=args.phase,
        band=args.band, fc_method=args.fc_method,
    )
    if M is None:
        raise SystemExit(
            f"No cached FC matrix for {args.patient}/{args.band}/"
            f"{args.phase}/{args.fc_method}.  Compute it first."
        )

    chnames = (
        _load_channel_labels(args.patient)
        if args.tick_labels == "chnames" else None
    )

    fig, _ax, _im = plot_fc_adjacency(
        M,
        fc_method=args.fc_method,
        band=args.band,
        tick_labels=args.tick_labels,
        channel_labels=chnames,
        log_scale=args.log_scale,
    )

    if args.watermark:
        add_provenance_footer(
            fig,
            f"single_adjacency · {args.patient} · "
            f"{BRAIN_BAND_TEX_DICT.get(args.band, args.band)} · "
            f"{PHASE_TEX.get(args.phase, args.phase)} · {args.fc_method}",
        )

    # Folder is already named single_adjacency/ — don't repeat in filename.
    scale_tag = "log" if args.log_scale else "lin"
    out = (
        Path(args.out_dir)
        / f"{args.patient}_{args.band}_{args.phase}_"
          f"{args.fc_method}_{args.tick_labels}_{scale_tag}.pdf"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
