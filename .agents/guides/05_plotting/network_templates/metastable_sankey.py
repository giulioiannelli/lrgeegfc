"""metastable_sankey — interactive Plotly Sankey of cluster evolution.

Crystallises the recipe from
``ipynb/06_presentation_figures/05_metastable_sankey.ipynb``:

  1. Build the giant component of the FC graph at one ``(patient, band,
     phase, fc_method)`` cell.
  2. Re-compute the LRG distance ``T_rho(τ)`` at each τ in
     ``DEFAULT_TAU_NCLUST_PAIRS`` (small τ → fine structure, large τ →
     coarse structure).
  3. Cut the τ-specific UPGMA dendrogram at the target cluster count to
     get a partition; collect partitions across τ.
  4. Emit an interactive Plotly Sankey HTML colouring each ribbon by
     its target cluster.

The sibling ``metastable_sankey_mpl.py`` reuses the same data prep and
emits a static PDF.

Usage
-----
    conda activate lapbrain
    python .agents/guides/05_plotting/network_templates/metastable_sankey.py \\
        --patient Pat_05 --band beta --phase rest_pre

    # Different τ-schedule (must be monotone-decreasing in cluster count
    # for a clean alluvial)
    python .agents/guides/05_plotting/network_templates/metastable_sankey.py \\
        --patient Pat_05 --band beta --phase rest_pre \\
        --tau-pairs 0.1:15 0.3:12 0.5:8 1.0:5 2.0:3 5.0:2

See companion ``metastable_sankey.md`` for the style sheet.
"""
from __future__ import annotations

import argparse
from pathlib import Path

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
    make_metastable_sankey_plotly,
)
from lrg_eegfc.visuals.network_templates import _load_inputs

try:  # graph-tool optional at import time; we only need nx here
    from lrgsglib.core import get_giant_component
except ImportError:  # pragma: no cover
    def get_giant_component(G: nx.Graph) -> nx.Graph:
        nodes = max(nx.connected_components(G), key=len)
        return G.subgraph(nodes).copy()


PHASE_TEX = {
    "rest_pre": r"rest_pre",
    "task_learn": r"task_learn",
    "task_test": r"task_test",
    "rest_post": r"rest_post",
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
        help="Magnitude FC only (Laplacian requires non-negative).",
    )
    parser.add_argument(
        "--tau-pairs",
        nargs="+",
        type=_parse_pair,
        default=list(DEFAULT_TAU_NCLUST_PAIRS),
        help="Space-separated τ:n_clusters pairs.  Default: "
             "0.1:15 0.5:8 1.0:5 2.0:3 5.0:2.  Cluster counts MUST be "
             "monotone-decreasing for a clean Sankey layout.",
    )
    parser.add_argument(
        "--width", type=int, default=1400, help="Plotly figure width (px)."
    )
    parser.add_argument(
        "--height", type=int, default=700, help="Plotly figure height (px)."
    )
    parser.add_argument(
        "--out-dir",
        default=str(
            FIGURES_ROOT / "network_templates" / "metastable_sankey"
        ),
    )
    args = parser.parse_args()

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

    partitions, n_clusters = compute_clustering_across_tau(
        None, tau_values, Gcc, n_clusters_list=n_clusters_list,
    )

    flow = compute_sankey_flows(partitions, tau_values, node_labels)

    band_tex = BRAIN_BAND_TEX_DICT.get(args.band, args.band)
    phase_tex = PHASE_TEX.get(args.phase, args.phase)
    title = (
        f"{args.patient} · {band_tex} · {phase_tex} · {args.fc_method}"
        f" — cluster evolution across τ"
    )

    fig = make_metastable_sankey_plotly(
        flow, title=title, width=args.width, height=args.height,
    )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tau_tag = "_".join(f"t{t:g}n{n}" for t, n in args.tau_pairs)
    out = (
        out_dir
        / f"{args.patient}_{args.band}_{args.phase}_{args.fc_method}_"
          f"{tau_tag}.html"
    )
    fig.write_html(str(out))
    print(f"wrote {out}")
    print(f"size: {out.stat().st_size / 1024:.0f} KB")


if __name__ == "__main__":
    main()
