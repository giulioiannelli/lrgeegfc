#!/usr/bin/env python
"""Generate dense-vs-CReMa comparison figures for all patients.

For each patient/band/phase produces a 4-panel figure:
  Top-left:  Dense MSC matrix       Top-right:  CReMa-validated matrix
  Bot-left:  Dense network graph    Bot-right:  CReMa network graph

Follows codebase plotting conventions (viridis, dpi=300, etc.).
"""

from __future__ import annotations

from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT

# -----------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------

PATIENTS_PHASES = {
    "Pat_02": ["rest_pre", "rest_post", "task_learn", "task_test"],
    "Pat_03": ["rest_pre", "rest_post", "task_learn", "task_test"],
    "Pat_05": ["rest_pre", "rest_post", "task_learn", "task_test"],
    "Pat_06": ["rest_pre", "rest_post"],
    "Pat_07": ["rest_pre", "rest_post", "task_learn"],
    "Pat_08": ["rest_pre", "rest_post", "task_learn", "task_test"],
}

BANDS = list(BRAIN_BANDS.keys())
NPERSEG = 4096
from lrg_eegfc.config.paths import FIGURES_ROOT
OUTPUT_ROOT = FIGURES_ROOT / "crema_comparison"

# Edge drawing constants (from visuals/lrg.py)
EDGE_POWER = 2.0
MAX_WIDTH = 4.0


# -----------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------

def load_pair(patient, phase, band):
    """Load dense and CReMa-validated MSC matrices."""
    dense = load_msc_matrix(patient, phase, band, sparsify="none", nperseg=NPERSEG)
    crema = load_msc_matrix(
        patient, phase, band, sparsify="ecm_adaptive", nperseg=NPERSEG,
        ecm_alpha_min=0.01, ecm_alpha_max=0.50,
        ecm_n_ensemble=100, ecm_weight_scale=1000,
    )
    return dense, crema


def graph_stats(A):
    """Quick stats string for a matrix."""
    N = A.shape[0]
    triu = np.triu_indices(N, k=1)
    w = A[triu]
    nz = (w > 0).sum()
    density = nz / len(w)
    mean_w = w[w > 0].mean() if nz > 0 else 0
    return nz, density, mean_w


def draw_network(ax, A, pos, title, color_accent="#4C72B0"):
    """Draw a network on ax using codebase conventions."""
    N = A.shape[0]
    G = nx.Graph()
    G.add_nodes_from(range(N))

    triu_i, triu_j = np.triu_indices(N, k=1)
    edges = []
    weights = []
    for idx in range(len(triu_i)):
        w = A[triu_i[idx], triu_j[idx]]
        if w > 0:
            i, j = int(triu_i[idx]), int(triu_j[idx])
            G.add_edge(i, j, weight=float(w))
            edges.append((i, j))
            weights.append(float(w))

    if weights:
        weights_arr = np.array(weights)
        w_norm = weights_arr / weights_arr.max()
        # Power-law edge scaling (from visuals/lrg.py)
        scaled = np.power(w_norm, EDGE_POWER)
        widths = MAX_WIDTH * scaled
        # Edge colors from viridis with alpha
        cmap = plt.cm.get_cmap("viridis")
        edge_colors = [(*cmap(wn)[:3], max(0.08, a)) for wn, a in zip(w_norm, scaled)]

        nx.draw_networkx_edges(
            G, pos, ax=ax, edgelist=edges,
            width=widths.tolist(), edge_color=edge_colors,
        )

    # Node degree for coloring
    degrees = np.array([G.degree(n) for n in range(N)])
    if degrees.max() > 0:
        deg_norm = degrees / degrees.max()
    else:
        deg_norm = np.zeros(N)

    nx.draw_networkx_nodes(
        G, pos, ax=ax,
        node_size=80,
        node_color=deg_norm,
        cmap=plt.cm.YlOrRd,
        vmin=0, vmax=1,
        edgecolors="black",
        linewidths=0.4,
    )

    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.axis("off")


# -----------------------------------------------------------------------
# Main figure: 4-panel (matrix + network) × (dense + CReMa)
# -----------------------------------------------------------------------

def plot_comparison(patient, phase, band, dense, crema, output_dir):
    """Generate 4-panel dense-vs-CReMa comparison."""
    N = dense.shape[0]
    n_dense, dens_dense, mw_dense = graph_stats(dense)
    n_crema, dens_crema, mw_crema = graph_stats(crema)

    band_name = BRAIN_BAND_TEX_DICT.get(band, band)

    # Shared layout: compute from CReMa backbone (sparser → better layout)
    G_layout = nx.from_numpy_array(crema)
    try:
        pos = nx.spring_layout(G_layout, seed=43, k=0.3, iterations=100)
    except Exception:
        pos = nx.circular_layout(G_layout)

    # Figure
    fig = plt.figure(figsize=(18, 16))
    gs = fig.add_gridspec(2, 2, hspace=0.25, wspace=0.20)

    # --- Top row: matrices ---
    ax_mat_d = fig.add_subplot(gs[0, 0])
    ax_mat_c = fig.add_subplot(gs[0, 1])

    # Dense matrix
    vmax = max(dense.max(), crema.max(), 0.01)
    im0 = ax_mat_d.imshow(dense, cmap="viridis", vmin=0, vmax=vmax,
                           aspect="equal", interpolation="none")
    divider0 = make_axes_locatable(ax_mat_d)
    cax0 = divider0.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im0, cax=cax0)
    ax_mat_d.set_title(
        f"Dense MSC\n{n_dense} edges, density={dens_dense:.3f}, <w>={mw_dense:.4f}",
        fontsize=11, fontweight="bold",
    )
    ax_mat_d.set_xlabel("Channel", fontsize=10)
    ax_mat_d.set_ylabel("Channel", fontsize=10)

    # CReMa matrix
    im1 = ax_mat_c.imshow(crema, cmap="viridis", vmin=0, vmax=vmax,
                           aspect="equal", interpolation="none")
    divider1 = make_axes_locatable(ax_mat_c)
    cax1 = divider1.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im1, cax=cax1)
    pct = 100 * n_crema / n_dense if n_dense > 0 else 0
    ax_mat_c.set_title(
        f"CReMa-validated MSC\n{n_crema} edges ({pct:.1f}%), density={dens_crema:.3f}, <w>={mw_crema:.4f}",
        fontsize=11, fontweight="bold",
    )
    ax_mat_c.set_xlabel("Channel", fontsize=10)
    ax_mat_c.set_ylabel("Channel", fontsize=10)

    # --- Bottom row: networks ---
    ax_net_d = fig.add_subplot(gs[1, 0])
    ax_net_c = fig.add_subplot(gs[1, 1])

    draw_network(ax_net_d, dense, pos, "Dense Network")
    draw_network(ax_net_c, crema, pos, "CReMa-validated Network")

    # Suptitle
    fig.suptitle(
        f"{patient} — {band_name} — {phase}",
        fontsize=14, fontweight="bold", y=0.98,
    )

    # Save
    out_dir = output_dir / patient
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{band}_{phase}_crema_comparison.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


# -----------------------------------------------------------------------
# Band overview: 6 bands in one figure per patient/phase
# -----------------------------------------------------------------------

def plot_band_overview(patient, phase, matrices, output_dir):
    """6-band overview: top row = dense matrices, bottom = CReMa matrices."""
    n_bands = len(BANDS)
    fig, axes = plt.subplots(2, n_bands, figsize=(4 * n_bands, 8))
    fig.suptitle(
        f"{patient} — {phase} — Dense vs CReMa-validated MSC",
        fontsize=14, fontweight="bold", y=1.01,
    )

    for col, band in enumerate(BANDS):
        dense, crema = matrices[band]
        band_name = BRAIN_BAND_TEX_DICT.get(band, band)
        n_d, d_d, _ = graph_stats(dense)
        n_c, d_c, _ = graph_stats(crema)
        vmax = max(dense.max(), crema.max(), 0.01)

        # Dense
        ax = axes[0, col]
        im = ax.imshow(dense, cmap="viridis", vmin=0, vmax=vmax,
                       aspect="equal", interpolation="none")
        ax.set_title(f"{band_name}\n{n_d} edges", fontsize=9, fontweight="bold")
        ax.axis("off")
        if col == 0:
            ax.set_ylabel("Dense", fontsize=11, fontweight="bold")

        # CReMa
        ax = axes[1, col]
        im = ax.imshow(crema, cmap="viridis", vmin=0, vmax=vmax,
                       aspect="equal", interpolation="none")
        pct = 100 * n_c / n_d if n_d > 0 else 0
        ax.set_title(f"{n_c} edges ({pct:.0f}%)", fontsize=9)
        ax.axis("off")
        if col == 0:
            ax.set_ylabel("CReMa", fontsize=11, fontweight="bold")

    plt.tight_layout()
    out_dir = output_dir / patient
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"all_bands_{phase}_crema_overview.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------

def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    total = 0
    failed = 0

    for patient, phases in PATIENTS_PHASES.items():
        for phase in phases:
            band_matrices = {}
            for band in BANDS:
                dense, crema = load_pair(patient, phase, band)
                if dense is None or crema is None:
                    print(f"  SKIP {patient}/{band}/{phase} — missing cache")
                    failed += 1
                    continue

                band_matrices[band] = (dense, crema)

                # Per-band 4-panel figure
                out = plot_comparison(patient, phase, band, dense, crema, OUTPUT_ROOT)
                total += 1
                n_d, _, _ = graph_stats(dense)
                n_c, _, _ = graph_stats(crema)
                print(f"  {patient}/{band}/{phase}: {n_d} -> {n_c} edges  [{out.name}]")

            # Band overview for this patient/phase
            if len(band_matrices) == len(BANDS):
                out = plot_band_overview(patient, phase, band_matrices, OUTPUT_ROOT)
                print(f"  -> overview: {out.name}")

    print(f"\nDone: {total} figures, {failed} skipped")
    print(f"Output: {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()
