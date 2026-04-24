#!/usr/bin/env python3
"""Section 2 Figure F: network metrics directly on the graph layout.

3 subpanels per method:
  1. Node strength       — sum of |FC| weights on each node.
  2. Distance-weighted betweenness centrality — bridge nodes on
     shortest paths through the inverse-weight distance graph.
  3. Current-flow betweenness centrality — random-walk betweenness;
     surfaces peripheral connector nodes (often anti-correlated with
     strength, so genuinely orthogonal to panel 1).

These three were chosen to be **mutually independent** on dense
weighted FC graphs. The earlier choices (clustering, participation)
turned out to track strength to r ≈ 0.96 / 0.51 and produced
visually-similar panels — replaced.

Each metric encoded by node SIZE (per-panel winsorised at 5/95th
percentile); nodes coloured by electrode shaft (no colorbar — same
shaft mapping as fig_E2).

Run:
  python scripts/10_notes_imcoh/fig_F_metrics_on_graph.py [-v]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import sys
sys.path.insert(0, str(Path(__file__).parent))

import matplotlib.pyplot as plt
import networkx as nx

from lrg_eegfc.utils.metrics.hypothesis import (
    BRAIN_BAND_TEX_DICT,
    REPR_PATIENTS, REPR_BANDS, REPR_PHASES,
    load_channel_labels, extract_probe_labels,
    apply_pub_style, save_fig, SECTION2_ROOT,
    draw_network_edges, spring_k_for, compute_network_layout,
    SECTION2_METHODS as METHODS,
    SECTION2_METHOD_LABELS as METHOD_LABELS
)
from lrg_eegfc.workflow.fc import load_fc_matrix



def fig_f_metrics_on_graph(
    patient: str,
    band: str,
    phase: str,
    edge_pct: float = 0.15,
    output_dir: Path = SECTION2_ROOT / "fig_F",
    verbose: bool = False,
):
    """Generate 4-panel metric-on-graph figures for both MSC and ImCoh."""
    ch = load_channel_labels(patient)
    probes = extract_probe_labels(ch)
    N = len(ch)

    # Shaft colors used uniformly across all panels (single legend mode:
    # node color = shaft, node size = metric).
    unique_shafts = sorted(set(probes))
    shaft_cmap = plt.get_cmap("tab20", len(unique_shafts))
    shaft_colors = [shaft_cmap(unique_shafts.index(p)) for p in probes]

    for method in METHODS:
        mat = load_fc_matrix(patient, phase, band, method)
        if mat is None:
            print(f"  SKIP {method}: no data for {patient} {phase} {band}")
            continue

        A = np.abs(mat.copy())
        np.fill_diagonal(A, 0)

        # fig_F-specific override: force `mds_lrg_continuous` for ALL
        # methods (including MSC). fig_E2 keeps its config-driven choice
        # (MSC = spring) so the contrast there stays — this override is
        # local to fig_F where same-layout-across-methods is what's
        # needed for the 4-panel metric overlay to be comparable.
        pos_arr = compute_network_layout(
            A, fc_method=method,
            patient=patient, phase=phase, band=band,
            layout_method="mds_lrg_continuous",
        )

        # Compute three mutually-independent metrics on the full graph.
        # NOTE: clustering & participation removed — both correlated
        # ≥ 0.95 / 0.51 with strength on dense FC graphs and produced
        # visually-redundant panels.
        strengths = A.sum(axis=1)
        # Distance-weighted graph: high |FC| → short distance.
        D = np.where(A > 0, 1.0 / np.clip(A, 1e-6, None), 0.0)
        np.fill_diagonal(D, 0.0)
        G_dist = nx.from_numpy_array(D)
        bc_w = np.array(
            list(nx.betweenness_centrality(
                G_dist, weight="weight", normalized=True).values())
        )
        cfb = np.array(
            list(nx.current_flow_betweenness_centrality(
                G_dist, weight="weight", normalized=True).values())
        )

        # 3 panels — every panel uses the same encoding:
        #   node color = electrode shaft (no colormap, no colorbar)
        #   node size  = the metric value of that panel
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        titles = ["Node strength",
                  "Betweenness centrality (dist = 1/|FC|)",
                  "Current-flow betweenness centrality"]
        metrics = [strengths, bc_w, cfb]

        for ax, title, metric in zip(axes, titles, metrics):
            # Draw edges with method-aware gamma scaling (weak edges fade)
            draw_network_edges(ax, pos_arr, A, fc_method=method,
                               probe_labels=probes,
                               highlight_same_probe=True)

            # Node size encodes THIS panel's metric (per-panel scale).
            # Winsorise at 5th/95th percentile before rescaling so a
            # single outlier doesn't compress the bulk distribution
            # into a narrow size band. Each panel gets its OWN
            # rescaling — different metrics → genuinely different
            # node-size patterns instead of a flat [30, 230] uniform.
            lo = float(np.percentile(metric, 5))
            hi = float(np.percentile(metric, 95))
            if hi > lo:
                s_norm = np.clip((metric - lo) / (hi - lo), 0.0, 1.0)
                node_sizes = 25 + 275 * s_norm
            else:
                node_sizes = 60

            ax.scatter(pos_arr[:, 0], pos_arr[:, 1],
                       c=shaft_colors, s=node_sizes,
                       edgecolors="white", linewidths=0.5, zorder=5)

            ax.set_title(title, fontsize=12, fontweight="bold")
            ax.axis("off")

            margin = 0.1
            xmin, xmax = pos_arr[:, 0].min(), pos_arr[:, 0].max()
            ymin, ymax = pos_arr[:, 1].min(), pos_arr[:, 1].max()
            span = max(xmax - xmin, ymax - ymin)
            ax.set_xlim(xmin - margin * span, xmax + margin * span)
            ax.set_ylim(ymin - margin * span, ymax + margin * span)

        label = METHOD_LABELS[method]
        band_tex = BRAIN_BAND_TEX_DICT[band]
        fig.tight_layout()
        save_fig(fig, output_dir / f"fig_F_{method}_{patient}_{band}_{phase}")

        if verbose:
            print(f"    {method}: strength [{strengths.min():.3f}, {strengths.max():.3f}], "
                  f"participation [{participation.min():.3f}, {participation.max():.3f}]")


def main():
    parser = argparse.ArgumentParser(description="Section 2 network metrics on graph.")
    parser.add_argument("--output-dir", type=Path, default=SECTION2_ROOT)
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    apply_pub_style()

    for patient in REPR_PATIENTS:
        for band in REPR_BANDS:
            for phase in REPR_PHASES:
                print(f"--- Figure F: {patient} {band} {phase} ---")
                fig_f_metrics_on_graph(
                    patient, band, phase,
                    output_dir=args.output_dir / "fig_F",
                    verbose=args.verbose,
                )

    print("\nDone.")


if __name__ == "__main__":
    main()
