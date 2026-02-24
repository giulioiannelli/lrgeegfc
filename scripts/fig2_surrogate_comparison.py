#!/usr/bin/env python3
"""
Figure 2: MSC Validation Comparison - Multiple Surrogate Levels

Creates a 4-row figure comparing:
- Row 1: Dense (no surrogates)
- Row 2: Validated with 50 surrogates
- Row 3: Validated with 100 surrogates
- Row 4: Validated with 200 surrogates

Run in background, saves to data/figures/presentation_figures/
"""

import sys
import os
from pathlib import Path

# Setup paths
ROOT = Path(__file__).resolve().parent.parent
os.chdir(ROOT)
sys.path.insert(0, str(ROOT / "src"))

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import networkx as nx
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

from lrg_eegfc.workflow.msc import load_msc_matrix, compute_msc_matrix

# ============================================================================
# CONFIGURATION
# ============================================================================
PATIENT = "Pat_02"
PHASE = "rsPre"
BAND = "beta"
SURROGATE_LEVELS = [0, 50, 100, 200]  # 0 = dense (no validation)

LAYOUT_K = 0.1
LAYOUT_SEED = 42
EDGE_POWER = 2.0  # Same as Figure 1
MAX_WIDTH = 4.0

MSC_CACHE = Path("data/msc_cache")
OUTPUT_DIR = Path("data/figures/presentation_figures") / PATIENT
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# Load or compute MSC matrices
# ============================================================================
print("=" * 60)
print("FIGURE 2: MSC Validation Comparison")
print("=" * 60)
print(f"Patient: {PATIENT}")
print(f"Phase: {PHASE}")
print(f"Band: {BAND}")
print(f"Surrogate levels: {SURROGATE_LEVELS}")
print()

msc_matrices = {}

for n_surr in SURROGATE_LEVELS:
    print(f"\n[{n_surr:3d} surrogates] ", end="", flush=True)

    if n_surr == 0:
        # Dense matrix
        matrix = load_msc_matrix(
            PATIENT, PHASE, BAND,
            cache_root=MSC_CACHE,
            sparsify="none",
            n_surrogates=0,
            nperseg=1024,
        )
        label = "Dense"
    else:
        # Validated matrix
        matrix = load_msc_matrix(
            PATIENT, PHASE, BAND,
            cache_root=MSC_CACHE,
            sparsify="soft",
            n_surrogates=n_surr,
            nperseg=1024,
        )
        label = f"{n_surr} surrogates"

    if matrix is not None:
        print(f"Loaded from cache - shape {matrix.shape}, mean={matrix.mean():.4f}")
        msc_matrices[n_surr] = matrix
    else:
        print(f"Computing (n_workers=1)...")
        result = compute_msc_matrix(
            PATIENT, PHASE, BAND,
            sparsify="soft" if n_surr > 0 else "none",
            n_surrogates=n_surr,
            verbose=True,
            n_workers=1,  # Sequential to avoid memory issues
        )
        msc_matrices[n_surr] = result.adjacency_matrix
        print(f"  Done - shape {result.adjacency_matrix.shape}, mean={result.adjacency_matrix.mean():.4f}")

print("\n" + "=" * 60)
print("All matrices loaded/computed")
print("=" * 60)

# ============================================================================
# Plotting function
# ============================================================================
def plot_surrogate_comparison(
    matrices: dict,
    surrogate_levels: list,
    patient: str,
    phase: str,
    band: str,
    k: float = 0.1,
    seed: int = 42,
    figsize: tuple = (16, 24),
    edge_power: float = 2.0,
    max_width: float = 4.0,
) -> plt.Figure:
    """Create 4-row figure comparing different surrogate levels.

    Each row: MSC matrix | Network visualization
    """
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(len(surrogate_levels), 2, figure=fig, hspace=0.25, wspace=0.2)

    # Get colormap (same as Figure 1 - magma)
    cmap = "magma"
    cmap_obj = plt.cm.get_cmap(cmap)

    # Compute shared layout based on dense matrix
    G_dense = nx.from_numpy_array(matrices[0])
    shared_pos = nx.spring_layout(G_dense, seed=seed, k=k, iterations=100)

    for row, n_surr in enumerate(surrogate_levels):
        matrix = matrices[n_surr]
        G = nx.from_numpy_array(matrix)

        # Get edge stats
        triu_idx = np.triu_indices_from(matrix, k=1)
        edges_vals = matrix[triu_idx]
        mean_val = edges_vals.mean()
        std_val = edges_vals.std()

        # Label
        if n_surr == 0:
            label = "Dense (No Validation)"
        else:
            label = f"Validated ({n_surr} surrogates)"

        # =================================================================
        # Left column: MSC Matrix
        # =================================================================
        ax_matrix = fig.add_subplot(gs[row, 0])

        im = ax_matrix.imshow(
            matrix,
            cmap=cmap,
            vmin=0,
            vmax=1,
            aspect='equal',
            interpolation='none',
        )

        # Colorbar
        divider = make_axes_locatable(ax_matrix)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)
        cbar.set_label("MSC", fontsize=11)

        ax_matrix.set_title(
            f"{label}\nmean={mean_val:.4f}, std={std_val:.4f}",
            fontsize=12,
            fontweight='bold'
        )
        ax_matrix.set_xlabel("Channel", fontsize=10)
        ax_matrix.set_ylabel("Channel", fontsize=10)

        # =================================================================
        # Right column: Network visualization
        # =================================================================
        ax_network = fig.add_subplot(gs[row, 1])

        # Get edge weights
        edges = list(G.edges(data=True))
        weights = np.array([e[2]['weight'] for e in edges])

        # Power law scaling for width (same as Figure 1)
        scaled = np.power(weights, edge_power)
        widths = max_width * scaled

        # Color edges using same colormap as matrix (magma)
        alphas = np.power(weights, edge_power)
        edge_colors = [(*cmap_obj(w)[:3], a) for w, a in zip(weights, alphas)]

        # Draw network
        nx.draw_networkx_nodes(
            G, shared_pos, ax=ax_network,
            node_size=80,
            node_color='steelblue',
            alpha=0.9,
            edgecolors='white',
            linewidths=0.5,
        )

        nx.draw_networkx_edges(
            G, shared_pos, ax=ax_network,
            width=widths,
            edge_color=edge_colors,
        )

        # Stats
        n_nodes = G.number_of_nodes()
        n_edges = G.number_of_edges()
        n_visible = np.sum(weights > 0.1)

        ax_network.set_title(
            f"{label}\n{n_nodes} nodes, {n_visible} edges with MSC > 0.1",
            fontsize=12,
            fontweight='bold'
        )
        ax_network.axis('off')

    # Overall title
    fig.suptitle(
        f"{patient} | {phase} | {band.upper()} Band\n"
        "MSC Validation: Comparing Surrogate Levels",
        fontsize=16,
        fontweight='bold',
        y=0.995
    )

    return fig


# ============================================================================
# Generate figure
# ============================================================================
print("\nGenerating figure...")

fig = plot_surrogate_comparison(
    msc_matrices,
    SURROGATE_LEVELS,
    PATIENT,
    PHASE,
    BAND,
    k=LAYOUT_K,
    seed=LAYOUT_SEED,
    edge_power=EDGE_POWER,
    max_width=MAX_WIDTH,
)

output_path = OUTPUT_DIR / f"fig2_msc_surrogate_comparison_{BAND}.pdf"
fig.savefig(output_path, dpi=150, bbox_inches='tight')
print(f"\nSaved: {output_path}")

# Also save PNG for quick preview
output_png = OUTPUT_DIR / f"fig2_msc_surrogate_comparison_{BAND}.png"
fig.savefig(output_png, dpi=100, bbox_inches='tight')
print(f"Saved: {output_png}")

plt.close(fig)

print("\n" + "=" * 60)
print("FIGURE 2 GENERATION COMPLETE")
print("=" * 60)
