# Figure Patterns Guide

Templates and patterns for generating publication-quality figures.

---

## Standard Setup

### Required Imports
```python
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Project imports
from lrg_eegfc.notebook import *
move_to_root(pathname="lrgeegfc")

from lrg_eegfc.config import BRAIN_BANDS, PHASE_LABELS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.workflow import load_corr_matrix, load_msc_matrix, load_lrg_result
from lrg_eegfc.visuals import (
    plot_correlation_and_network,
    plot_lrg_full_panel,
    plot_msc_and_network,
)
```

### Standard Constants
```python
# DPI settings
DPI_STANDARD = 300      # Publication quality
DPI_3D = 150            # 3D and interactive

# Figure sizes by type
FIGSIZE_SINGLE = (8, 8)
FIGSIZE_DOUBLE = (16, 8)
FIGSIZE_TRIPLE = (18, 6)
FIGSIZE_LRG_PANEL = (20, 12)
FIGSIZE_PHASE = (34, 20)
FIGSIZE_3D = (10, 8)

# Font sizes
FONTSIZE_TITLE = 14
FONTSIZE_LABEL = 11
FONTSIZE_TICK = 9
FONTSIZE_LEGEND = 10

# Colormaps
CMAP_CORRELATION = "viridis"      # For absolute values [0,1]
CMAP_SIGNED = "coolwarm"          # For signed values [-1,1]
CMAP_MSC = "viridis"              # MSC [0,1]
CMAP_ULTRAMETRIC = "viridis"      # Ultrametric distances
CMAP_DIFFERENCE = "Reds"          # Difference/removed edges
CMAP_DIVERGING = "RdBu_r"         # Diverging comparisons
```

---

## Pattern 1: Single Matrix Heatmap

```python
def plot_matrix(matrix: np.ndarray, output_path: Path, title: str = "Matrix",
                vmin: float = 0, vmax: float = 1, cmap: str = "viridis"):
    """Single matrix heatmap with colorbar."""
    fig, ax = plt.subplots(figsize=(8, 8))

    im = ax.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")
    ax.set_title(title, fontsize=14)
    ax.set_xlabel("Channel", fontsize=11)
    ax.set_ylabel("Channel", fontsize=11)

    # Colorbar with divider
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path
```

---

## Pattern 2: Matrix + Network Side-by-Side

```python
def plot_matrix_network(matrix: np.ndarray, output_path: Path,
                        patient: str, phase: str, band: str):
    """Side-by-side matrix and network graph."""
    import networkx as nx

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    # Left: Matrix
    ax = axes[0]
    im = ax.imshow(np.abs(matrix), cmap="viridis", vmin=0, vmax=1)
    ax.set_title(f"{patient} {phase} {band} - Matrix", fontsize=14)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)

    # Right: Network
    ax = axes[1]
    threshold = np.percentile(np.abs(matrix[matrix != 0]), 75)
    adj = np.where(np.abs(matrix) > threshold, np.abs(matrix), 0)
    G = nx.from_numpy_array(adj)
    G.remove_edges_from(nx.selfloop_edges(G))

    pos = nx.spring_layout(G, k=0.3, seed=42)
    edges = G.edges()
    weights = [G[u][v]['weight'] for u, v in edges]

    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=80, node_color="lightblue")
    nx.draw_networkx_edges(G, pos, ax=ax, width=[w * 2 for w in weights], alpha=0.6)
    ax.set_title(f"Network (threshold={threshold:.2f})", fontsize=14)
    ax.axis("off")

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path
```

---

## Pattern 3: Full LRG Panel (5-panel)

```python
def plot_lrg_panel(patient: str, phase: str, band: str, fc_method: str = "msc"):
    """5-panel LRG analysis using built-in function."""
    from lrg_eegfc.visuals.lrg import plot_lrg_full_panel

    output_path = Path(f"data/figures/lrg/{patient}/{band}_{phase}_{fc_method}_full.png")

    return plot_lrg_full_panel(
        patient=patient,
        phase=phase,
        band=band,
        fc_method=fc_method,
        cache_root=Path("data/lrg_cache"),
        dataset_root=Path("data/stereoeeg_patients"),
        output_path=output_path,
        figsize=(20, 12),
    )
```

---

## Pattern 4: Multi-Phase Comparison

```python
def plot_phase_comparison(patient: str, band: str, fc_method: str = "msc"):
    """Compare LRG across all phases."""
    from lrg_eegfc.visuals.reorganization import plot_phase_reorganization

    output_path = Path(f"data/figures/reorganization/{patient}/{band}_{fc_method}_phases.png")

    return plot_phase_reorganization(
        patient=patient,
        band=band,
        fc_method=fc_method,
        phases=["rest_pre", "task_learn", "task_test", "rest_post"],
        cache_root=Path("data/lrg_cache"),
        output_path=output_path,
        figsize=(34, 20),
    )
```

---

## Pattern 5: GridSpec Complex Layout

```python
def plot_complex_panel(patient: str, phase: str, band: str):
    """Custom GridSpec layout for complex figures."""
    fig = plt.figure(figsize=(20, 12))
    gs = gridspec.GridSpec(4, 6, figure=fig, hspace=0.3, wspace=0.4)

    # Top-left: 2x2 matrix
    ax_matrix = fig.add_subplot(gs[0:2, 0:2])

    # Top-middle: 2x2 entropy
    ax_entropy = fig.add_subplot(gs[0:2, 2:4])

    # Right: Full-height dendrogram
    ax_dendro = fig.add_subplot(gs[0:4, 4:6])

    # Bottom-left: 2x3 network
    ax_network = fig.add_subplot(gs[2:4, 0:3])

    # Bottom-middle: 2x1 PSI
    ax_psi = fig.add_subplot(gs[2:4, 3:4])

    # ... populate axes ...

    output_path = Path(f"data/figures/custom/{patient}/{band}_{phase}_panel.png")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path
```

---

## Pattern 6: 3D Spatial Network

```python
def plot_3d_brain(patient: str, phase: str, band: str, fc_method: str = "msc"):
    """Interactive 3D brain network."""
    from lrg_eegfc.visuals.spatial import plot_spatial_network_3d, view_brain_connectome

    # Plotly interactive
    fig_plotly = plot_spatial_network_3d(
        patient=patient,
        phase=phase,
        band=band,
        fc_method=fc_method,
        edge_threshold=0.3,
        max_edges=200,
        node_size=8,
        colorby="cluster",
    )

    # Nilearn glass brain
    output_path = Path(f"data/figures/spatial/{patient}/{band}_{phase}_brain.png")
    view_brain_connectome(
        patient=patient,
        phase=phase,
        band=band,
        fc_method=fc_method,
        edge_threshold=0.3,
    )

    return fig_plotly
```

---

## Pattern 7: Batch Figure Generation

```python
def generate_patient_figures(patient: str):
    """Generate all standard figures for a patient."""
    from lrg_eegfc.config import BRAIN_BANDS_NAMES, PHASE_LABELS
    from lrg_eegfc.visuals import (
        plot_correlation_and_network,
        plot_msc_and_network,
        plot_lrg_full_panel,
    )

    figures = []

    for phase in PHASE_LABELS:
        for band in BRAIN_BANDS_NAMES:
            # Correlation figure
            try:
                fig = plot_correlation_and_network(
                    patient=patient, phase=phase, band=band,
                    output_path=Path(f"data/figures/correlation/{patient}/{band}_{phase}.png")
                )
                figures.append(fig)
            except FileNotFoundError:
                print(f"Skip {patient}/{phase}/{band} - no cache")

            # LRG figure
            try:
                fig = plot_lrg_full_panel(
                    patient=patient, phase=phase, band=band, fc_method="msc",
                    output_path=Path(f"data/figures/lrg/{patient}/{band}_{phase}_msc.png")
                )
                figures.append(fig)
            except FileNotFoundError:
                print(f"Skip LRG {patient}/{phase}/{band}")

    return figures
```

---

## Pattern 8: Entropy Curves

```python
def plot_entropy_curves(tau: np.ndarray, entropy_1_minus_S: np.ndarray,
                        entropy_C: np.ndarray, output_path: Path):
    """Dual entropy curves with log scale."""
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(tau, entropy_1_minus_S, 'b-', linewidth=2, label=r'$1 - S$')
    ax.plot(tau, entropy_C, 'r--', linewidth=2, label=r'$C$')

    ax.set_xscale('log')
    ax.set_xlabel(r'$\tau$', fontsize=12)
    ax.set_ylabel('Entropy', fontsize=12)
    ax.set_title('LRG Entropy Curves', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path
```

---

## Pattern 9: Dendrogram with Threshold

```python
from scipy.cluster.hierarchy import dendrogram

def plot_dendrogram_custom(linkage_matrix: np.ndarray, threshold: float,
                           labels: list, output_path: Path):
    """Dendrogram with optimal threshold line."""
    fig, ax = plt.subplots(figsize=(10, 8))

    dendro = dendrogram(
        linkage_matrix,
        labels=labels,
        orientation='right',
        color_threshold=threshold,
        above_threshold_color='gray',
        leaf_font_size=8,
        ax=ax,
    )

    ax.axvline(x=threshold, color='red', linestyle='--', linewidth=2,
               label=f'Threshold = {threshold:.3f}')
    ax.set_xscale('log')
    ax.set_xlabel('Ultrametric Distance', fontsize=12)
    ax.set_title('Hierarchical Clustering', fontsize=14)
    ax.legend(fontsize=10)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path
```

---

## Output Path Conventions

```
data/figures/
├── correlation/
│   └── Pat_XX/
│       └── {band}_{phase}_corr.png
├── msc/
│   └── Pat_XX/
│       └── {band}_{phase}_msc.png
├── lrg/
│   └── Pat_XX/
│       └── {band}_{phase}_{fc_method}_full.png
├── reorganization/
│   └── Pat_XX/
│       └── {band}_{fc_method}_phases.png
├── spatial/
│   └── Pat_XX/
│       └── {band}_{phase}_brain.png
└── comparison/
    └── Pat_XX/
        └── {band}_{metric}.png
```

---

## Quick Checklist

Before saving any figure:
- [ ] `output_path.parent.mkdir(parents=True, exist_ok=True)`
- [ ] `fig.savefig(..., dpi=300, bbox_inches="tight")`
- [ ] `plt.close(fig)`
- [ ] Use cache loaders, never recompute
- [ ] Include patient/phase/band in title or filename
