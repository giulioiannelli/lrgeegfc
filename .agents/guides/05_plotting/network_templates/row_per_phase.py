"""row_per_phase — 1 row × 4 phases of FC networks at fixed band.

Thin wrapper around ``lrg_eegfc.visuals.plot_fc_network_row``.  Useful
for "how does the network reorganise across rest_pre → task_learn →
task_test → rest_post" at a single band.

Usage
-----
    conda activate lapbrain
    python .agents/guides/05_plotting/network_templates/row_per_phase.py \
        --patient Pat_05 --band beta --layout kk

See sibling ``row_per_phase.md`` for the style sheet.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from lrg_eegfc.config.const import PHASE_LABELS
from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.visuals.network_templates import LAYOUT_REGISTRY, plot_fc_network_row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--patient", default="Pat_05")
    parser.add_argument("--band", default="beta")
    parser.add_argument(
        "--fc-method", default="imcoh_abs",
        choices=["corr", "msc", "imcoh", "imcoh_abs", "imcoh_sq"],
    )
    parser.add_argument(
        "--layout", default="spring",
        choices=[k for k in LAYOUT_REGISTRY.keys() if not k.startswith("lrg_")],
        help="Layout algorithm (LRG-seeded layouts excluded — they need "
             "different community labels per phase, see lrg_seeded template).",
    )
    parser.add_argument(
        "--node-color", default="shaft",
        choices=["shaft", "uniform"],
    )
    parser.add_argument(
        "--coloring", default=None,
        choices=["probe", "signed", "weight"],
    )
    parser.add_argument(
        "--out-dir",
        default=str(FIGURES_ROOT / "network_templates" / "row_per_phase"),
    )
    args = parser.parse_args()

    fig = plot_fc_network_row(
        args.patient, args.band,
        fc_method=args.fc_method,
        phases=PHASE_LABELS,
        layout=args.layout,
        coloring=args.coloring,
        node_color=args.node_color,
    )

    coloring_tag = args.coloring or "auto"
    out = (
        Path(args.out_dir)
        / f"{args.patient}_{args.band}_{args.fc_method}_"
          f"{args.layout}_{args.node_color}_{coloring_tag}.pdf"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
