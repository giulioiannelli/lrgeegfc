"""mosaic_band_phase — 4 phases × 6 bands of FC networks for one patient.

Thin wrapper around ``lrg_eegfc.visuals.plot_fc_network_grid``.  The
network analogue of ``fc_templates/mosaic_band_phase``: a full
spectral × phase decomposition of one patient's FC graph at a glance.

Usage
-----
    conda activate lapbrain
    python .agents/guides/05_plotting/network_templates/mosaic_band_phase.py \
        --patient Pat_05 --layout kk

See sibling ``mosaic_band_phase.md`` for the style sheet.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, PHASE_LABELS
from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.visuals.network_templates import LAYOUT_REGISTRY, plot_fc_network_grid


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--patient", default="Pat_05")
    parser.add_argument(
        "--fc-method", default="imcoh_abs",
        choices=["corr", "msc", "imcoh", "imcoh_abs", "imcoh_sq"],
    )
    parser.add_argument(
        "--layout", default="spring",
        choices=[k for k in LAYOUT_REGISTRY.keys() if not k.startswith("lrg_")],
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
        default=str(FIGURES_ROOT / "network_templates" / "mosaic_band_phase"),
    )
    args = parser.parse_args()

    fig = plot_fc_network_grid(
        args.patient,
        fc_method=args.fc_method,
        phases=PHASE_LABELS,
        bands=BRAIN_BANDS_NAMES,
        layout=args.layout,
        coloring=args.coloring,
        node_color=args.node_color,
    )

    coloring_tag = args.coloring or "auto"
    out = (
        Path(args.out_dir)
        / f"{args.patient}_{args.fc_method}_"
          f"{args.layout}_{args.node_color}_{coloring_tag}.pdf"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
