"""metastable_sankey_mpl — static matplotlib alluvial of cluster evolution.

Vector-PDF twin of ``metastable_sankey.py``.  Shares the data-prep
pipeline (``compute_clustering_across_tau`` + ``compute_sankey_flows``)
with the Plotly template; only the renderer differs.

Layout:

  * One column per τ at x = i · ``x_step``.
  * Stacked cluster bars within each column; bar height ∝ cluster size,
    bar fill = qualitative-palette colour by cluster ID.
  * Cubic-bezier ribbon polygons between consecutive columns, one per
    ``(c_t → c_{t+1})`` transition.  Ribbon width = transition count,
    ribbon colour = target-bar colour at ``ribbon_alpha`` (default 0.55).
  * Per-column header (``τ = X (k clusters)``) above the top bar.

PDF only, fully vector — never rasterise the ribbon patches.

Usage
-----
    conda activate lapbrain
    python .agents/guides/05_plotting/network_templates/metastable_sankey_mpl.py \\
        --patient Pat_05 --band beta --phase rest_pre

    # Custom τ-schedule (must be monotone-decreasing in cluster count)
    python .agents/guides/05_plotting/network_templates/metastable_sankey_mpl.py \\
        --patient Pat_05 --band beta --phase rest_pre \\
        --tau-pairs 0.1:15 0.3:12 0.5:8 0.8:6 1.0:5 2.0:3 5.0:2

See companion ``metastable_sankey_mpl.md`` for the style sheet.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    PHASE_LABELS,
)
from lrg_eegfc.config.paths import FIGURES_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_patient_metadata
from lrg_eegfc.visuals.metastable import (
    DEFAULT_TAU_NCLUST_PAIRS,
    compute_clustering_across_tau,
    compute_sankey_flows,
    make_metastable_sankey_mpl,
)
from lrg_eegfc.visuals.network_templates import _load_inputs
from lrg_eegfc.visuals.styles import use_lrg_style

try:
    from lrgsglib.core import get_giant_component
except ImportError:  # pragma: no cover
    def get_giant_component(G: nx.Graph) -> nx.Graph:
        nodes = max(nx.connected_components(G), key=len)
        return G.subgraph(nodes).copy()


PHASE_TEX = {
    "rest_pre": r"rest$_{\mathrm{pre}}$",
    "task_learn": r"task$_{\mathrm{learn}}$",
    "task_test": r"task$_{\mathrm{test}}$",
    "rest_post": r"rest$_{\mathrm{post}}$",
}


def _parse_pair(s: str) -> tuple[float, int]:
    t, n = s.split(":")
    return float(t), int(n)


def _parse_contact_label(raw: str) -> str:
    head = str(raw).split(",")[0]
    return "".join(c for c in head if c.isascii() and c.isalnum())


def _load_contact_labels(patient: str) -> list[str]:
    meta = load_patient_metadata(patient, SEEG_DATAPATH)
    if meta is None or "label" not in meta.columns:
        return []
    return [_parse_contact_label(s) for s in meta["label"].tolist()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--patient", default="Pat_05")
    parser.add_argument("--band", default="beta", choices=list(BRAIN_BANDS_NAMES))
    parser.add_argument("--phase", default="rest_pre", choices=PHASE_LABELS)
    parser.add_argument(
        "--fc-method", default="imcoh_abs",
        choices=["imcoh_abs", "imcoh_sq", "msc", "corr"],
    )
    parser.add_argument(
        "--tau-pairs",
        nargs="+",
        type=_parse_pair,
        default=list(DEFAULT_TAU_NCLUST_PAIRS),
        help="Space-separated τ:n_clusters pairs.",
    )
    parser.add_argument(
        "--figsize",
        nargs=2,
        type=float,
        default=(12.0, 6.0),
        help="Figure (width height) in inches.",
    )
    parser.add_argument(
        "--node-width", type=float, default=0.06,
        help="Cluster-bar width in x-axis units.  Default 0.06.",
    )
    parser.add_argument(
        "--ribbon-alpha", type=float, default=0.55,
        help="Ribbon opacity.  Default 0.55 (slightly airier than the "
             "Plotly default of 0.6).",
    )
    parser.add_argument(
        "--no-bar-labels", action="store_true",
        help="Suppress the ``C{c} ({n})`` text on each bar.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(
            FIGURES_ROOT / "network_templates" / "metastable_sankey_mpl"
        ),
    )
    args = parser.parse_args()

    use_lrg_style()

    A, _ = _load_inputs(args.patient, args.band, args.phase, args.fc_method)
    np.fill_diagonal(A, 0.0)
    G = nx.from_numpy_array(A)
    Gcc = get_giant_component(G)

    tau_values = np.array([t for t, _ in args.tau_pairs])
    n_clusters_list = [n for _, n in args.tau_pairs]

    node_labels = _load_contact_labels(args.patient)
    if not node_labels:
        node_labels = [f"Node_{i}" for i in range(Gcc.number_of_nodes())]
    node_labels = node_labels[: Gcc.number_of_nodes()]

    partitions, _ = compute_clustering_across_tau(
        None, tau_values, Gcc, n_clusters_list=n_clusters_list,
    )

    flow = compute_sankey_flows(partitions, tau_values, node_labels)

    band_tex = BRAIN_BAND_TEX_DICT.get(args.band, args.band)
    phase_tex = PHASE_TEX.get(args.phase, args.phase)
    title = (
        f"{args.patient} · {band_tex} · {phase_tex} · {args.fc_method}"
        r"    cluster evolution across $\tau$"
    )

    fig, ax = make_metastable_sankey_mpl(
        flow,
        figsize=tuple(args.figsize),
        node_width=args.node_width,
        ribbon_alpha=args.ribbon_alpha,
        show_bar_labels=not args.no_bar_labels,
        title=title,
    )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tau_tag = "_".join(f"t{t:g}n{n}" for t, n in args.tau_pairs)
    out = (
        out_dir
        / f"{args.patient}_{args.band}_{args.phase}_{args.fc_method}_"
          f"{tau_tag}.pdf"
    )
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")
    print(f"size: {out.stat().st_size / 1024:.0f} KB")


if __name__ == "__main__":
    main()
