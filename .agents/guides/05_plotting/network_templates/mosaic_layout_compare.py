"""mosaic_layout_compare — every layout side-by-side for one case.

Thin wrapper around ``lrg_eegfc.visuals.plot_layout_gallery``.  This is
the "test all layouts" template: same FC matrix rendered with every
layout in :data:`DEFAULT_GALLERY` (spring, kk, spectral, laplacian_pca,
community_grouped, backbone_guided, circular_by_shaft, sfdp, arf,
lrg_kk, lrg_sfdp).  LRG-seeded layouts are dropped silently if no
LRG cache is present.

Usage
-----
    conda activate lapbrain
    python .agents/guides/05_plotting/network_templates/mosaic_layout_compare.py \
        --patient Pat_05 --band beta --phase rest_pre

See sibling ``mosaic_layout_compare.md`` for the style sheet.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.visuals.network_templates import DEFAULT_GALLERY, plot_layout_gallery


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--patient", default="Pat_05")
    parser.add_argument("--band", default="beta")
    parser.add_argument("--phase", default="rest_pre")
    parser.add_argument(
        "--fc-method", default="imcoh_abs",
        choices=["corr", "msc", "imcoh", "imcoh_abs", "imcoh_sq"],
    )
    parser.add_argument(
        "--layouts", default=None,
        help="Comma-separated subset of layouts (default: full gallery).",
    )
    parser.add_argument(
        "--node-color", default="shaft",
        choices=["shaft", "community", "uniform"],
    )
    parser.add_argument(
        "--coloring", default=None,
        choices=["probe", "signed", "weight"],
    )
    parser.add_argument(
        "--n-communities", type=int, default=10,
        help="LRG dendrogram cut for LRG-seeded layouts.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(FIGURES_ROOT / "network_templates" / "mosaic_layout_compare"),
    )
    args = parser.parse_args()

    layouts = (
        DEFAULT_GALLERY if args.layouts is None
        else tuple(args.layouts.split(","))
    )

    fig = plot_layout_gallery(
        args.patient, args.band, args.phase,
        fc_method=args.fc_method,
        layouts=layouts,
        n_communities=args.n_communities,
        coloring=args.coloring,
        node_color=args.node_color,
    )

    coloring_tag = args.coloring or "auto"
    out = (
        Path(args.out_dir)
        / f"{args.patient}_{args.band}_{args.phase}_"
          f"{args.fc_method}_{args.node_color}_{coloring_tag}.pdf"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
