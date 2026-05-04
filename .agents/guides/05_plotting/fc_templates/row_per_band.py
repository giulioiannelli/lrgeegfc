"""row_per_band — one patient × one phase × the six bands.

Thin wrapper around ``lrg_eegfc.visuals.plot_fc_adjacency_row``.
Defaults to **per-panel** colour scale (γ matrices have much smaller
amplitude than α/β; a shared linear scale would crush them).  Pass
``--shared-scale`` for a single colorbar with a global ``vmax``, or
``--log-scale`` to compress the dynamic range.

Usage
-----
    conda activate lapbrain
    python .agents/guides/05_plotting/fc_templates/row_per_band.py \\
        --patient Pat_06 --phase rest_pre --fc-method imcoh_abs \\
        --tick-labels generic
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT, BRAIN_BANDS_NAMES
from lrg_eegfc.config.paths import FIGURES_ROOT, SEEG_DATAPATH
from lrg_eegfc.visuals.fc_templates import plot_fc_adjacency_row
from lrg_eegfc.visuals.layout import add_provenance_footer
from lrg_eegfc.workflow.fc import load_fc_matrix


def _load_channel_labels(patient: str) -> list[str]:
    csv = SEEG_DATAPATH / patient / "channel_labels.csv"
    return pd.read_csv(csv)["label"].astype(str).tolist()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--patient", default="Pat_06")
    parser.add_argument("--phase", default="rest_pre")
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
        "--shared-scale", action="store_true",
        help="One global vmax + a single right-side colorbar "
             "(default: per-panel scale + per-panel colorbar).",
    )
    parser.add_argument(
        "--out-dir",
        default=str(FIGURES_ROOT / "fc_templates" / "row_per_band"),
    )
    args = parser.parse_args()

    matrices = []
    for band in BRAIN_BANDS_NAMES:
        M = load_fc_matrix(
            patient=args.patient, phase=args.phase,
            band=band, fc_method=args.fc_method,
        )
        if M is None:
            raise SystemExit(
                f"No cached FC for {args.patient}/{band}/"
                f"{args.phase}/{args.fc_method}.  Compute it first."
            )
        matrices.append(M)

    chnames = (
        _load_channel_labels(args.patient)
        if args.tick_labels == "chnames" else None
    )

    fig, _axes = plot_fc_adjacency_row(
        matrices,
        titles=[BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES],
        fc_method=args.fc_method,
        bands=BRAIN_BANDS_NAMES,
        tick_labels=args.tick_labels,
        channel_labels=chnames,
        log_scale=args.log_scale,
        shared_scale=args.shared_scale,
    )

    if args.watermark:
        add_provenance_footer(
            fig,
            f"row_per_band · {args.patient} · {args.phase} · "
            f"{args.fc_method}",
        )

    scale_tag = "log" if args.log_scale else "lin"
    share_tag = "shared" if args.shared_scale else "perpanel"
    out = (
        Path(args.out_dir)
        / f"{args.patient}_{args.phase}_{args.fc_method}_"
          f"{args.tick_labels}_{scale_tag}_{share_tag}.pdf"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
