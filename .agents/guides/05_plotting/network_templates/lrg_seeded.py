"""lrg_seeded — LRG-driven network drawing (the cluster figure).

Thin wrapper around ``lrg_eegfc.visuals.plot_fc_network_lrg``.  Uses
the LRG dendrogram cut at ``n_communities`` as a "groups" prior to
force modules apart on the canvas.  Two backends:

- ``backend=gt`` (default): graph-tool SFDP with ``groups=`` from the
  LRG cut.  Multilevel, scales well; ``sfdp_gamma`` controls
  module-separation strength.
- ``backend=nx``: NetworkX Kamada-Kawai with inflated inter-community
  distances; ``separation`` controls module separation.

Default ``node_color=community`` because the figure's whole point is
showing LRG modules.  Switch to ``shaft`` to confirm modules aren't
probe-trivial.

Usage
-----
    conda activate lapbrain
    python .agents/guides/05_plotting/network_templates/lrg_seeded.py \
        --patient Pat_05 --band beta --phase rest_post --n-communities 10

See sibling ``lrg_seeded.md`` for the style sheet.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from lrg_eegfc.config.paths import FIGURES_ROOT
from lrg_eegfc.visuals.network_templates import plot_fc_network_lrg


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
        "--n-communities", type=int, default=10,
        help="LRG dendrogram cut (default 10).",
    )
    parser.add_argument(
        "--backend", default="gt", choices=["gt", "nx"],
        help="Layout backend: gt = graph-tool SFDP (default), "
             "nx = NetworkX Kamada-Kawai.",
    )
    parser.add_argument(
        "--separation", type=float, default=5.0,
        help="NetworkX KK inter-community distance multiplier.",
    )
    parser.add_argument(
        "--sfdp-gamma", type=float, default=0.1,
        help="graph-tool SFDP module-separation gamma "
             "(0.03=weak, 0.6=strong).",
    )
    parser.add_argument(
        "--node-color", default="community",
        choices=["community", "shaft", "uniform"],
    )
    parser.add_argument(
        "--coloring", default=None,
        choices=["probe", "signed", "weight"],
    )
    parser.add_argument(
        "--out-dir",
        default=str(FIGURES_ROOT / "network_templates" / "lrg_seeded"),
    )
    args = parser.parse_args()

    fig = plot_fc_network_lrg(
        args.patient, args.band, args.phase,
        fc_method=args.fc_method,
        n_communities=args.n_communities,
        backend=args.backend,
        separation=args.separation,
        sfdp_gamma=args.sfdp_gamma,
        coloring=args.coloring,
        node_color=args.node_color,
    )

    coloring_tag = args.coloring or "auto"
    out = (
        Path(args.out_dir)
        / f"{args.patient}_{args.band}_{args.phase}_"
          f"{args.fc_method}_n{args.n_communities}_"
          f"{args.backend}_{args.node_color}_{coloring_tag}.pdf"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
