# Instructions: Implement LRG Analysis Visualizations

**Date:** 2025-12-11
**Task:** Implement LRG (Laplacian Renormalization Group) analysis visualization suite
**Goal:** Create comprehensive 4-panel visualization for both Correlation-based and MSC-based FC networks

---

## Overview

You will implement LRG visualization functions that analyze functional connectivity networks using:
1. **Laplacian-based ultrametric distances** - Hierarchical structure from graph Laplacian
2. **Optimal partitioning** - Data-driven threshold selection for community detection
3. **Partition stability analysis** - PSI (Partition Stability Index) to assess quality
4. **Thermodynamic observables** - Entropy and specific heat from diffusion dynamics

The visualization will be applied to **both** correlation and MSC adjacency matrices.

---

## Reference Implementation

### Existing Notebooks

**Primary reference:** `ipynb/FIGMNTGN01.ipynb`
- 4-panel layout with custom gridspec
- Implements PSI calculation
- Shows entropy/specific heat plots
- Network visualization with partition colors

**Secondary reference:** `ipynb/FIGMNTGN02.ipynb`
- Multi-phase comparison layout
- Consistent dendrogram coloring approach

### Key Visualization from FIGMNTGN01

```
┌─────────────────┬─────────────────┬─────────────────────────┐
│  (a) Matrix     │  (b) Entropy    │  (c) Dendrogram         │
│  Correlation/   │  C (speC heat)  │  Hierarchical tree      │
│  MSC Heatmap    │  S (entropy)    │  Color threshold        │
│                 │  vs log(tau)    │  Optimal partition      │
├─────────────────┴─────────────────┤                         │
│  (d) Network                      │                         │
│  Graph with partition colors      │                         │
│  from dendrogram clustering       │                         │
│                                   │                         │
└───────────────────────────────────┴─────────────────────────┘
                  (e) PSI Plot
              Partition Stability Index
```

---

## Implementation Requirements

### File Structure

Create new file: `src/lrg_eegfc/visuals/lrg.py`

```python
"""LRG analysis visualizations for functional connectivity networks."""

from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import networkx as nx
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster
from scipy.spatial.distance import squareform

__all__ = [
    "plot_lrg_analysis",
    "compute_partition_stability_index",
]
```

### Dependencies from lrgsglib

The implementation requires these functions from `lrgsglib`:

```python
from lrgsglib.utils.lrg.clustering import compute_normalized_linkage, compute_optimal_threshold
from lrgsglib.utils.lrg.spectral import compute_laplacian_properties
from lrgsglib.utils.lrg.observables import entropy
```

**Important:** These are external library functions that MUST be available.

---

## Core Helper Functions

### Function 1: Partition Stability Index (PSI)

```python
def compute_partition_stability_index(linkage_matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the Partition Stability Index from hierarchical clustering.

    The PSI quantifies the stability of a partition with n communities by comparing
    the ultrametric gap when merging from n+1 to n communities versus merging from
    n to n-1 communities. Large PSI values indicate stable partitions.

    Parameters
    ----------
    linkage_matrix : np.ndarray
        Linkage matrix from hierarchical clustering (N-1, 4) where N is number of nodes.
        Column 2 contains the ultrametric distances (gaps) at each merge step.

    Returns
    -------
    psi_values : np.ndarray
        Partition Stability Index values for each partition level.
        Normalized to [-1, 1] range where peaks indicate stable partitions.
    n_communities : np.ndarray
        Number of communities corresponding to each PSI value.

    Notes
    -----
    The PSI is computed as:
        PSI(n) = N * (log10(Δ_n) - log10(Δ_{n+1}))

    where:
        - N is the total number of nodes
        - Δ_n is the ultrametric gap when merging from n+1 to n communities
        - Δ_{n+1} is the gap when merging from n to n-1 communities

    Peaks in PSI indicate partition levels where the hierarchical structure
    is particularly stable (large gap before merge, small gap after).

    Examples
    --------
    >>> psi, n_comm = compute_partition_stability_index(linkage_matrix)
    >>> optimal_n = n_comm[np.argmax(-psi)]  # Find most stable partition
    >>> print(f"Optimal partition has {optimal_n} communities")
    """
    N = linkage_matrix.shape[0] + 1  # Number of original nodes

    # Extract ultrametric gaps (distances at each merge step)
    # linkage_matrix[:, 2] contains distances, sorted from smallest to largest
    deltas = linkage_matrix[:, 2]

    # Initialize arrays
    psi_values = []
    n_communities = []

    # Compute PSI for each partition level n
    # At step i, we have N-i communities (before merging to N-i-1)
    for i in range(len(deltas) - 1):
        n = N - i - 1  # Number of communities after this merge

        delta_n = deltas[i]          # Gap when going from n+1 to n communities
        delta_n_plus_1 = deltas[i + 1]  # Gap when going from n to n-1 communities

        # Avoid log(0) by adding small epsilon if needed
        if delta_n <= 0:
            delta_n = 1e-10
        if delta_n_plus_1 <= 0:
            delta_n_plus_1 = 1e-10

        # Compute Partition Stability Index
        psi = N * (np.log10(delta_n) - np.log10(delta_n_plus_1))

        psi_values.append(psi)
        n_communities.append(n)

    # Normalize PSI to [-1, 1] for visualization
    psi_array = np.array(psi_values)
    psi_normalized = -psi_array / np.max(np.abs(psi_array))

    return psi_normalized, np.array(n_communities)
```

---

## Main Visualization Function

### Function Signature

```python
def plot_lrg_analysis(
    patient: str,
    phase: str,
    band: str,
    fc_method: str = "correlation",
    cache_root_corr: Path = Path("data/corr_cache"),
    cache_root_msc: Path = Path("data/msc_cache"),
    dataset_root: Path = Path("data/stereoeeg_patients"),
    output_path: Optional[Path] = None,
    nperseg: int = 1024,
    scaling_factor: float = 0.98,
    figsize: tuple = (16, 10),
    tau_max: float = 3.0,
    verbose: bool = False,
) -> Path:
    """Create comprehensive LRG analysis visualization.

    Generates a 4-panel layout analyzing functional connectivity using
    Laplacian-based ultrametric distances and hierarchical clustering:

    Layout:
    ┌─────────┬─────────┬───────────┐
    │ Matrix  │ Entropy │ Dendrogram│
    │         │ + Heat  │           │
    ├─────────┴─────────┤           │
    │ Network           │           │
    │ (partitioned)     │           │
    └───────────────────┴───────────┘
              PSI Plot

    Parameters
    ----------
    patient : str
        Patient identifier (e.g., "Pat_02")
    phase : str
        Recording phase (e.g., "rsPre", "taskLearn")
    band : str
        Frequency band (e.g., "beta", "alpha")
    fc_method : str
        FC method: "correlation" or "msc" (default: "correlation")
    cache_root_corr : Path
        Root directory for correlation cache
    cache_root_msc : Path
        Root directory for MSC cache
    dataset_root : Path
        Root directory for patient data (to get channel labels)
    output_path : Path, optional
        Output file path. If None, uses default location
    nperseg : int
        Window length for MSC (must match cached files)
    scaling_factor : float
        Scaling factor for optimal threshold selection (default: 0.98)
        Higher values → more conservative clustering (more communities)
    figsize : tuple
        Figure size (width, height)
    tau_max : float
        Maximum diffusion time for entropy calculation (log10 scale)
    verbose : bool
        Print progress information

    Returns
    -------
    Path
        Path to saved figure

    Raises
    ------
    ImportError
        If lrgsglib is not available
    FileNotFoundError
        If cache files or patient data not found

    Notes
    -----
    This function requires the lrgsglib package for:
    - compute_laplacian_properties: Extract graph Laplacian and ultrametric distances
    - compute_normalized_linkage: Hierarchical clustering from distance matrix
    - compute_optimal_threshold: Data-driven threshold selection
    - entropy: Compute thermodynamic observables (entropy, specific heat)

    The optimal partition is determined by the compute_optimal_threshold
    function, which finds the threshold maximizing partition stability.

    Examples
    --------
    >>> # Correlation-based LRG analysis
    >>> output = plot_lrg_analysis("Pat_02", "rsPre", "beta", fc_method="correlation")

    >>> # MSC-based LRG analysis
    >>> output = plot_lrg_analysis("Pat_02", "rsPre", "beta", fc_method="msc")
    """
```

---

## Implementation Steps

### Step 1: Import and Validation

```python
# Check lrgsglib availability
try:
    from lrgsglib.utils.lrg.clustering import (
        compute_normalized_linkage,
        compute_optimal_threshold
    )
    from lrgsglib.utils.lrg.spectral import compute_laplacian_properties
    from lrgsglib.utils.lrg.observables import entropy
except ImportError as e:
    raise ImportError(
        "lrgsglib package required for LRG analysis. "
        "Install with: pip install ./lrgsglib"
    ) from e

# Load appropriate FC matrix
if fc_method == "correlation":
    from lrg_eegfc.workflow_corr import load_corr_matrix

    fc_matrix = load_corr_matrix(
        patient, phase, band,
        cache_root=cache_root_corr,
        use_cleaned=True  # Use MP-cleaned correlation
    )

    if fc_matrix is None:
        raise FileNotFoundError(
            f"Correlation matrix not found: {patient} {phase} {band}\n"
            f"Run: python src/compute_corr_matrices.py --patient {patient}"
        )

    method_label = "Correlation"

elif fc_method == "msc":
    from lrg_eegfc.workflow_msc import load_msc_matrix

    fc_matrix = load_msc_matrix(
        patient, phase, band,
        cache_root=cache_root_msc,
        sparsify="none",  # Use dense MSC
        n_surrogates=0,
        nperseg=nperseg
    )

    if fc_matrix is None:
        raise FileNotFoundError(
            f"MSC matrix not found: {patient} {phase} {band}\n"
            f"Run: python src/compute_msc_matrices.py --patient {patient}"
        )

    method_label = "MSC"

else:
    raise ValueError(f"Unknown fc_method: {fc_method}. Use 'correlation' or 'msc'")

if verbose:
    print(f"Loaded {method_label} matrix: {fc_matrix.shape}")
```

### Step 2: Build Graph and Compute LRG Properties

```python
# Build NetworkX graph from FC matrix
G = nx.from_numpy_array(fc_matrix)

# Extract giant component
if not nx.is_connected(G):
    largest_cc = max(nx.connected_components(G), key=len)
    G = G.subgraph(largest_cc).copy()
    G = nx.convert_node_labels_to_integers(G)
    if verbose:
        print(f"Extracted giant component: {G.number_of_nodes()} nodes")

# Compute Laplacian-based ultrametric distances
# This computes: spectrum, Laplacian, ρ(tau), T_ρ(tau), tau
spectrum, L, rho, Trho, tau = compute_laplacian_properties(G, tau=None)

if verbose:
    print(f"Computed Laplacian properties")
    print(f"  Spectrum range: [{spectrum.min():.4f}, {spectrum.max():.4f}]")

# Convert to condensed distance matrix for clustering
dists = squareform(Trho)

# Compute hierarchical clustering (Ward linkage on ultrametric distances)
linkage_matrix, label_list, _ = compute_normalized_linkage(
    dists, G, method='ward'
)

if verbose:
    print(f"Computed normalized linkage matrix")

# Find optimal clustering threshold
FlatClusteringTh, optTh, stabs, opt_idx = compute_optimal_threshold(
    linkage_matrix,
    scaling_factor=scaling_factor
)

if verbose:
    print(f"Optimal threshold: {FlatClusteringTh:.6f}")

# Get optimal partition
optimal_clusters = fcluster(
    linkage_matrix,
    t=FlatClusteringTh,
    criterion='distance'
)

n_clusters = len(np.unique(optimal_clusters))
if verbose:
    print(f"Optimal partition: {n_clusters} communities")
```

### Step 3: Compute Thermodynamic Observables

```python
# Compute entropy and specific heat
# Returns: (entropy_vs_tau, specific_heat_vs_tau, ???, tau_values)
net_entropy = entropy(G, t1=-3, t2=tau_max)

# Extract components
tau_scale = net_entropy[-1]  # Diffusion time values (tau)
entropy_values = net_entropy[0]  # S(tau)
specific_heat = net_entropy[1]   # C(tau) = dS/d(log tau)

if verbose:
    print(f"Computed thermodynamic observables")
    print(f"  Tau range: [{tau_scale[0]:.4f}, {tau_scale[-1]:.4f}]")
```

### Step 4: Compute PSI

```python
# Compute Partition Stability Index
psi_values, n_communities = compute_partition_stability_index(linkage_matrix)

# Find most stable partition
max_psi_idx = np.argmax(-psi_values)  # Negative because we normalized
optimal_n_communities = n_communities[max_psi_idx]
max_psi_value = psi_values[max_psi_idx]

if verbose:
    print(f"PSI analysis:")
    print(f"  Most stable partition: {optimal_n_communities} communities")
    print(f"  PSI value: {max_psi_value:.4f}")
```

### Step 5: Load Channel Labels

```python
# Load channel labels for visualization
from lrg_eegfc.utils.datamanag.patient_robust import load_patient_dataset_robust

try:
    dataset = load_patient_dataset_robust(patient, dataset_root, phases=[phase])
    recording = dataset[phase]

    # Try to get channel labels
    if hasattr(recording, 'channel_labels') and recording.channel_labels is not None:
        channel_labels = {i: label for i, label in enumerate(recording.channel_labels)}
    else:
        channel_labels = {i: f"Ch{i}" for i in range(G.number_of_nodes())}

except Exception as e:
    if verbose:
        print(f"Could not load channel labels: {e}")
    channel_labels = {i: f"Ch{i}" for i in range(G.number_of_nodes())}
```

### Step 6: Create Figure Layout

```python
# Create custom grid layout matching FIGMNTGN01
fig = plt.figure(figsize=figsize)
gs = gridspec.GridSpec(4, 6, figure=fig, hspace=0.3, wspace=0.4)

# Define subplot areas
ax_matrix = fig.add_subplot(gs[0:2, 0:2])      # Top-left: Matrix
ax_entropy = fig.add_subplot(gs[0:2, 2:4])     # Top-mid: Entropy/Heat
ax_dendro = fig.add_subplot(gs[0:4, 4:6])      # Right: Dendrogram (full height)
ax_network = fig.add_subplot(gs[2:4, 0:3])     # Bottom-left: Network
ax_psi = fig.add_subplot(gs[2:4, 3:4])         # Bottom-mid-right: PSI

# Overall title
fig.suptitle(
    f"{patient} {phase} {band.upper()} - LRG Analysis ({method_label})",
    fontsize=16,
    fontweight='bold',
    y=0.98
)
```

### Step 7: Panel (a) - FC Matrix Heatmap

```python
# Plot FC matrix
if fc_method == "correlation":
    # Show absolute correlation
    plot_matrix = np.abs(fc_matrix)
    vmin, vmax = 0, 1
    cmap = 'viridis'
else:
    # MSC already in [0, 1]
    plot_matrix = fc_matrix
    vmin, vmax = 0, 1
    cmap = 'viridis'

im_matrix = ax_matrix.imshow(
    plot_matrix,
    cmap=cmap,
    vmin=vmin,
    vmax=vmax,
    origin='upper',
    aspect='auto'
)

ax_matrix.set_title(f'(a) {method_label} Matrix', fontsize=12, fontweight='bold')
ax_matrix.set_xlabel('Channel', fontsize=10)
ax_matrix.set_ylabel('Channel', fontsize=10)

# Colorbar
cbar_matrix = plt.colorbar(im_matrix, ax=ax_matrix, fraction=0.046, pad=0.04)
cbar_matrix.set_label(method_label, fontsize=10)
```

### Step 8: Panel (b) - Entropy and Specific Heat

```python
# Create twin axis for entropy and specific heat
ax_entropy_twin = ax_entropy.twinx()

# Plot specific heat on left y-axis (blue)
line_heat = ax_entropy.plot(
    tau_scale[1:], specific_heat, '-',
    label='C (Specific Heat)', color='blue', linewidth=2
)
ax_entropy.set_ylabel('C (Specific Heat)', color='blue', fontsize=11)
ax_entropy.tick_params(axis='y', labelcolor='blue')

# Plot entropy on right y-axis (red)
line_entropy = ax_entropy_twin.plot(
    tau_scale, entropy_values, '-',
    label='S (Entropy)', color='red', linewidth=2
)
ax_entropy_twin.set_ylabel('S (Entropy)', color='red', fontsize=11)
ax_entropy_twin.tick_params(axis='y', labelcolor='red')

# Set x-axis (log scale for tau)
ax_entropy.set_xscale('log')
ax_entropy.set_xlabel(r'$\tau$ (Diffusion Time)', fontsize=11)
ax_entropy.set_title('(b) Thermodynamic Observables', fontsize=12, fontweight='bold')

# Combined legend
lines_heat, labels_heat = ax_entropy.get_legend_handles_labels()
lines_entropy, labels_entropy = ax_entropy_twin.get_legend_handles_labels()
ax_entropy_twin.legend(
    lines_heat + lines_entropy,
    labels_heat + labels_entropy,
    loc='best',
    fontsize=9
).set_zorder(200)

ax_entropy.grid(alpha=0.3)
```

### Step 9: Panel (c) - Dendrogram with Optimal Threshold

```python
# Plot dendrogram (right orientation to save space)
node_list = list(G.nodes())
labels_for_dendro = [channel_labels.get(n, f"Ch{n}") for n in node_list]

# Generate color palette for communities
n_colors = min(n_clusters + 3, 20)  # At least n_clusters colors
palette = plt.cm.tab20(np.linspace(0, 1, n_colors))
palette = [plt.matplotlib.colors.to_hex(c) for c in palette]

from scipy.cluster import hierarchy
hierarchy.set_link_color_palette(palette)

# Create dendrogram
dendro = dendrogram(
    linkage_matrix,
    ax=ax_dendro,
    color_threshold=FlatClusteringTh,
    labels=labels_for_dendro,
    above_threshold_color='k',
    leaf_font_size=5,
    orientation='right'
)

# Set log scale and threshold line
tmin = linkage_matrix[:, 2][0] * 0.8
tmax = linkage_matrix[:, 2][-1] * 1.01
ax_dendro.set_xscale('log')
ax_dendro.axvline(FlatClusteringTh, color='b', linestyle='--', linewidth=2,
                  label=r'$\mathcal{D}_{\rm th}$')
ax_dendro.set_xlim(tmin, tmax)
ax_dendro.set_xlabel(r'$\mathcal{D}/\mathcal{D}_{\max}$', fontsize=11)
ax_dendro.set_title('(c) Hierarchical Tree', fontsize=12, fontweight='bold')
ax_dendro.legend(fontsize=9)
```

### Step 10: Panel (d) - Network with Partition Colors

```python
# Extract node colors from dendrogram
leaf_label_colors = {
    lbl: col for lbl, col in zip(dendro['ivl'], dendro['leaves_color_list'])
}

# Map colors to nodes
node_colors = [
    leaf_label_colors.get(channel_labels.get(n, f"Ch{n}"), 'gray')
    for n in node_list
]

# Compute network layout
pos = nx.spring_layout(G, seed=43, scale=1, k=0.27, iterations=50)

# Scale edge widths by weight
edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
wmin, wmax = min(edge_weights), max(edge_weights)
widths_scaled = [0.05 + 0.3 * (w - wmin) / (wmax - wmin) for w in edge_weights]

# Draw network
nx.draw(
    G,
    pos=pos,
    ax=ax_network,
    width=widths_scaled,
    node_color=node_colors,
    node_size=140,
    with_labels=True,
    labels=channel_labels,
    font_size=7,
    font_color='k'
)

ax_network.set_title(
    f'(d) Network (Optimal Partition: {n_clusters} communities)',
    fontsize=12,
    fontweight='bold'
)

# Set axis limits with small margin
x, y = np.array(list(pos.values())).T
ax_network.set_xlim(x.min() - 0.05, x.max() + 0.05)
ax_network.set_ylim(y.min() - 0.05, y.max() + 0.05)
```

### Step 11: Panel (e) - PSI Plot

```python
# Plot Partition Stability Index
ax_psi.plot(n_communities, psi_values, '-o', color='green', linewidth=2, markersize=4)

# Mark optimal partition
ax_psi.axvline(
    optimal_n_communities,
    ls='--',
    c='red',
    linewidth=2,
    label=f'Optimal: n={optimal_n_communities}'
)

# Also mark the partition we actually used (from compute_optimal_threshold)
ax_psi.axvline(
    n_clusters,
    ls=':',
    c='blue',
    linewidth=2,
    label=f'Used: n={n_clusters}'
)

ax_psi.set_xlabel(r'$n$ (Number of Communities)', fontsize=11)
ax_psi.set_ylabel(r'$\Psi(n, \tau)$ (PSI)', fontsize=11)
ax_psi.set_title('(e) Partition Stability', fontsize=12, fontweight='bold')
ax_psi.legend(fontsize=9)
ax_psi.grid(alpha=0.3)
ax_psi.set_xlim(1, min(20, max(n_communities)))
```

### Step 12: Save Figure

```python
if output_path is None:
    output_dir = Path("data/figures/lrg") / fc_method / patient
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{band}_{phase}_lrg_{fc_method}.png"

plt.tight_layout()
plt.savefig(output_path, dpi=150, bbox_inches='tight')
plt.close(fig)

if verbose:
    print(f"Saved LRG analysis: {output_path}")

return output_path
```

---

## CLI Integration

### Create Script: `src/visualize_lrg.py`

```python
#!/usr/bin/env python3
"""CLI script for LRG analysis visualizations.

Examples
--------
# Correlation-based LRG analysis
python src/visualize_lrg.py --patient Pat_02 --phase rsPre --band beta --fc-method correlation --verbose

# MSC-based LRG analysis
python src/visualize_lrg.py --patient Pat_02 --phase rsPre --band beta --fc-method msc --verbose

# Batch mode: all bands for correlation
python src/visualize_lrg.py --patient Pat_02 --phase rsPre --fc-method correlation --batch --verbose

# Both FC methods for one band/phase
python src/visualize_lrg.py --patient Pat_02 --phase rsPre --band beta --fc-method both --verbose
"""

import argparse
import sys
from pathlib import Path

from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS
from lrg_eegfc.visuals.lrg import plot_lrg_analysis


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate LRG analysis visualizations for functional connectivity"
    )

    # Required arguments
    parser.add_argument("--patient", required=True, help="Patient ID")
    parser.add_argument("--phase", help="Recording phase")
    parser.add_argument("--band", help="Frequency band")

    # FC method selection
    parser.add_argument(
        "--fc-method",
        choices=["correlation", "msc", "both"],
        default="correlation",
        help="FC method to analyze (default: correlation)"
    )

    # LRG parameters
    parser.add_argument(
        "--scaling-factor",
        type=float,
        default=0.98,
        help="Threshold scaling factor (default: 0.98)"
    )

    parser.add_argument(
        "--tau-max",
        type=float,
        default=3.0,
        help="Maximum diffusion time for entropy (default: 3.0)"
    )

    parser.add_argument(
        "--nperseg",
        type=int,
        default=1024,
        help="MSC window length (must match cache). Default: 1024"
    )

    # Paths
    parser.add_argument(
        "--cache-root-corr",
        type=Path,
        default=Path("data/corr_cache"),
        help="Correlation cache directory"
    )

    parser.add_argument(
        "--cache-root-msc",
        type=Path,
        default=Path("data/msc_cache"),
        help="MSC cache directory"
    )

    parser.add_argument(
        "--dataset-root",
        type=Path,
        default=Path("data/stereoeeg_patients"),
        help="Patient data directory"
    )

    # Batch mode
    parser.add_argument(
        "--batch",
        action="store_true",
        help="Generate for all bands (requires --phase and --fc-method)"
    )

    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    # Validation
    if args.batch and not args.phase:
        parser.error("--batch requires --phase")

    if not args.batch and (not args.phase or not args.band):
        parser.error("--phase and --band required unless using --batch")

    # Determine combinations
    if args.batch:
        bands = list(BRAIN_BANDS.keys())
        phases = [args.phase]
    else:
        bands = [args.band]
        phases = [args.phase]

    # Determine FC methods
    if args.fc_method == "both":
        fc_methods = ["correlation", "msc"]
    else:
        fc_methods = [args.fc_method]

    # Process
    for band in bands:
        for phase in phases:
            for fc_method in fc_methods:
                if args.verbose:
                    print(f"\nProcessing LRG analysis: {args.patient} {phase} {band} ({fc_method})")

                try:
                    output_path = plot_lrg_analysis(
                        patient=args.patient,
                        phase=phase,
                        band=band,
                        fc_method=fc_method,
                        cache_root_corr=args.cache_root_corr,
                        cache_root_msc=args.cache_root_msc,
                        dataset_root=args.dataset_root,
                        nperseg=args.nperseg,
                        scaling_factor=args.scaling_factor,
                        tau_max=args.tau_max,
                        verbose=args.verbose,
                    )

                    if args.verbose:
                        print(f"  ✓ Saved: {output_path}")

                except ImportError as e:
                    print(f"  ERROR: {e}")
                    print(f"  Install lrgsglib: pip install ./lrgsglib")
                    return 1
                except FileNotFoundError as e:
                    print(f"  ERROR: {e}")
                except Exception as e:
                    print(f"  ERROR: Failed for {args.patient} {phase} {band} ({fc_method}): {e}")
                    if args.verbose:
                        import traceback
                        traceback.print_exc()

    return 0


if __name__ == "__main__":
    sys.exit(main())
```

---

## Testing Instructions

### Test 1: Single LRG Analysis (Correlation)

```bash
# Ensure correlation matrix exists
python src/compute_corr_matrices.py --patient Pat_02 --verbose
python src/clean_correlation_matrices.py --patient Pat_02 --verbose

# Generate LRG analysis
python src/visualize_lrg.py --patient Pat_02 --phase rsPre --band beta --fc-method correlation --verbose
```

Expected output: `data/figures/lrg/correlation/Pat_02/beta_rsPre_lrg_correlation.png`

### Test 2: Single LRG Analysis (MSC)

```bash
# Ensure MSC matrix exists
python src/compute_msc_matrices.py --patient Pat_02 --verbose

# Generate LRG analysis
python src/visualize_lrg.py --patient Pat_02 --phase rsPre --band beta --fc-method msc --verbose
```

Expected output: `data/figures/lrg/msc/Pat_02/beta_rsPre_lrg_msc.png`

### Test 3: Both FC Methods

```bash
python src/visualize_lrg.py --patient Pat_02 --phase rsPre --band beta --fc-method both --verbose
```

Expected: 2 figures (one for correlation, one for MSC)

### Test 4: Batch Mode

```bash
# All bands for correlation
python src/visualize_lrg.py --patient Pat_02 --phase rsPre --fc-method correlation --batch --verbose
```

Expected: 6 correlation LRG figures (one per band)

---

## Integration with Pipeline

### Update `scripts/run_full_analysis.sh`

Add after comparison visualization (after Step 5c):

```bash
# -------------------------------------------------------------------------
# Step 6: Generate LRG analysis visualizations
# -------------------------------------------------------------------------
echo -e "${BLUE}[6/8]${NC} Generating LRG analysis visualizations..."

# Correlation-based LRG
if [ -f "src/visualize_lrg.py" ]; then
  if python src/visualize_lrg.py --patient "$PATIENT" --batch --fc-method correlation --verbose 2>/dev/null; then
    echo -e "${GREEN}✓${NC} Correlation LRG analysis generated"
  else
    echo -e "${YELLOW}⚠${NC} Correlation LRG analysis skipped or failed for $PATIENT"
  fi
else
  echo -e "${YELLOW}⚠${NC} LRG visualization script not available yet"
fi

# MSC-based LRG
if [ -f "src/visualize_lrg.py" ]; then
  if python src/visualize_lrg.py --patient "$PATIENT" --batch --fc-method msc --verbose 2>/dev/null; then
    echo -e "${GREEN}✓${NC} MSC LRG analysis generated"
  else
    echo -e "${YELLOW}⚠${NC} MSC LRG analysis skipped or failed for $PATIENT"
  fi
fi
```

---

## Expected Figure Output

### Panel Descriptions

**(a) FC Matrix:**
- Heatmap of correlation (absolute) or MSC matrix
- Colorbar: viridis, [0, 1]

**(b) Thermodynamic Observables:**
- Blue line (left y-axis): C (specific heat) vs log(τ)
- Red line (right y-axis): S (entropy) vs log(τ)
- Shows diffusion dynamics and scale-dependent structure

**(c) Dendrogram:**
- Hierarchical tree (right orientation)
- Log-scale x-axis: ultrametric distances
- Blue dashed line: optimal threshold
- Colors indicate communities

**(d) Network:**
- Spring layout graph
- Node colors from dendrogram communities
- Edge widths scaled by FC strength
- Channel labels displayed

**(e) PSI Plot:**
- Green line: PSI vs number of communities
- Red dashed line: optimal partition (from PSI)
- Blue dotted line: partition used (from compute_optimal_threshold)
- Shows partition stability landscape

---

## Additional Requirements

### Export Functions

Add to `src/lrg_eegfc/__init__.py`:

```python
from .visuals.lrg import plot_lrg_analysis, compute_partition_stability_index
```

### Documentation String for Module

```python
"""
LRG Analysis Visualization Module

This module provides visualization tools for Laplacian Renormalization Group (LRG)
analysis of functional connectivity networks. The LRG approach reveals hierarchical
community structure through ultrametric distances computed from the graph Laplacian.

Key features:
- Laplacian-based ultrametric distances (multi-scale structure)
- Optimal partition detection via data-driven threshold selection
- Partition Stability Index (PSI) for quality assessment
- Thermodynamic observables (entropy, specific heat)

The visualization can be applied to both correlation-based and MSC-based FC networks.

References
----------
.. [1] Laplacian Renormalization Group for heterogeneous networks
       De Bacco et al., Physical Review E (2016)
"""
```

---

## Validation Checklist

Before marking as complete:

- [ ] Function created in `src/lrg_eegfc/visuals/lrg.py`
- [ ] PSI function implemented and tested
- [ ] CLI script created: `src/visualize_lrg.py`
- [ ] Exports added to `__init__.py`
- [ ] Test with Pat_02 rsPre beta correlation
- [ ] Test with Pat_02 rsPre beta MSC
- [ ] Test batch mode (all bands)
- [ ] Verify 4-panel layout matches specification
- [ ] Verify dendrogram coloring is consistent with network
- [ ] Verify PSI peaks indicate stable partitions
- [ ] Verify entropy/specific heat plots display correctly
- [ ] Integration with pipeline script
- [ ] Works for both fs=1024 and fs=2048 data

---

## Key Design Decisions

1. **Use cleaned correlation matrices**: MP spectral cleaning improves signal
2. **Use dense MSC matrices**: Compare raw methods before sparsification
3. **Ward linkage**: Best for compact, spherical clusters
4. **Scaling factor = 0.98**: Conservative default (from FIGMNTGN01)
5. **Right-oriented dendrogram**: Saves horizontal space
6. **Spring layout**: Good default for brain networks
7. **Log scale for tau and distances**: Natural for multi-scale analysis
8. **PSI normalization**: Easier interpretation in [-1, 1] range

---

## Notes for Implementation

- Follow existing code style in `src/lrg_eegfc/visuals/corr.py` and `msc.py`
- Use consistent color schemes (viridis for matrices, tab20 for communities)
- Function should be ~400-500 lines including docstring
- CLI should follow same pattern as other visualization scripts
- Gracefully handle lrgsglib import errors with clear instructions
- Provide verbose output for debugging LRG computation steps

---

## Example Expected Output Statistics

For Pat_02 rsPre beta correlation:

```
Loaded Correlation matrix: (117, 117)
Extracted giant component: 117 nodes
Computed Laplacian properties
  Spectrum range: [0.0034, 3.2156]
Computed normalized linkage matrix
Optimal threshold: 0.049269
Optimal partition: 9 communities
Computed thermodynamic observables
  Tau range: [0.0010, 1000.0000]
PSI analysis:
  Most stable partition: 9 communities
  PSI value: -0.1234
Saved LRG analysis: data/figures/lrg/correlation/Pat_02/beta_rsPre_lrg_correlation.png
```

Values will vary by patient, phase, and band.

---

**End of Instructions**

When implementing this, you should be able to run:
```bash
python src/visualize_lrg.py --patient Pat_02 --phase rsPre --band beta --fc-method both --verbose
```

And get two publication-quality 4-panel LRG analysis figures (one for correlation, one for MSC).
