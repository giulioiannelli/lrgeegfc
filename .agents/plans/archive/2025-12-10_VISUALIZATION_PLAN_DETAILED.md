# Detailed Implementation Plan: Phases 3-7

This document provides comprehensive implementation details for the remaining phases of the visualization pipeline.

---

## PHASE 3: MSC Visualization Suite (DETAILED)

### Goal
Create visualization scripts for MSC analysis (matrices and networks only, no cleaning)

### Why MSC Doesn't Need Cleaning
- **User requirement**: "for the msc there is no need of thresholding or removing spectral noise"
- MSC is computed in frequency domain via Welch's method
- Already robust to temporal noise through spectral averaging
- Circular shift surrogates provide validation if needed (separate from cleaning)

### Implementation Details

#### File 1: `src/lrg_eegfc/visuals/msc.py`

**Complete implementation specifications:**

##### Function 1: plot_msc_heatmap()

```python
def plot_msc_heatmap(
    msc_matrix: np.ndarray,
    output_path: Path,
    title: str = "MSC Matrix",
    vmin: float = 0.0,  # MSC range is [0, 1] unlike correlation [-1, 1]
    vmax: float = 1.0,
    channel_labels: Optional[List[str]] = None,
    figsize: tuple = (8, 8),
    cmap: str = "viridis",
) -> Path:
    """Plot MSC matrix heatmap.

    Implementation pattern (adapted from correlation version):

    1. Create figure and axis
    2. Plot using imshow with viridis colormap
    3. Add colorbar using mpl_toolkits.axes_grid1.make_axes_locatable
    4. Set channel labels if provided (and if n_channels <= 50)
    5. Set title
    6. Save to output_path
    7. Close figure
    8. Return output_path

    Key differences from correlation:
    - vmin=0.0 (MSC is always non-negative)
    - vmax=1.0 (MSC is bounded by 1)
    - No need for symmetric colormap (RdBu not needed)

    Code pattern:
        fig, ax = plt.subplots(figsize=figsize)
        im = ax.imshow(msc_matrix, cmap=cmap, vmin=vmin, vmax=vmax,
                       interpolation='none', aspect='auto')

        # Add colorbar
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)

        # Labels
        if channel_labels is not None and len(channel_labels) <= 50:
            ax.set_xticks(range(len(channel_labels)))
            ax.set_yticks(range(len(channel_labels)))
            ax.set_xticklabels(channel_labels, rotation=90, fontsize=8)
            ax.set_yticklabels(channel_labels, fontsize=8)
        else:
            ax.set_xticks([])
            ax.set_yticks([])

        ax.set_title(title, fontsize=14)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(fig)

        return output_path
    """
```

##### Function 2: plot_msc_and_network()

```python
def plot_msc_and_network(
    patient: str,
    phase: str,
    band: str,
    cache_root: Path = Path("data/msc_cache"),
    output_path: Optional[Path] = None,
    sparsify: str = "none",
    n_surrogates: int = 0,
    nperseg: int = 256,
    figsize: tuple = (16, 8),
    dataset_root: Path = Path("data/stereoeeg_patients"),
) -> Path:
    """Side-by-side: MSC matrix + network graph with edge weights.

    Implementation steps:

    1. Load MSC matrix from cache
       - Use workflow_msc.load_msc_matrix() or get_msc_cache_path()
       - Handle both dense and validated versions based on sparsify param

    2. Load channel labels if available
       - Try loading from dataset_root/{patient}/channel_labels.mat
       - Fallback to integer labels if not found

    3. Create 2-panel figure
       - Left: MSC heatmap
       - Right: Network graph

    4. Left panel: MSC matrix heatmap
       - Call plot_msc_heatmap() internally or reuse code
       - vmin=0, vmax=1, viridis colormap

    5. Right panel: Network visualization
       - Build graph: G = nx.from_numpy_array(msc_matrix)
       - Extract edge weights: widths = [G[u][v]['weight'] for u, v in G.edges()]
       - Scale widths to visual range [0.1, 3.0]:
           wmin, wmax = min(widths), max(widths)
           widths_scaled = [0.1 + 2.9 * (w - wmin) / (wmax - wmin + 1e-10) for w in widths]
       - Choose layout: pos = nx.spring_layout(G, seed=42) or nx.kamada_kawai_layout(G)
       - Draw with networkx:
           nx.draw(G, pos=pos, ax=ax[1],
                   width=widths_scaled,
                   node_size=100,
                   node_color='lightblue',
                   edge_color='gray',
                   with_labels=True if n_nodes <= 30 else False,
                   labels=label_dict,
                   font_size=8)

    6. Add titles and annotations
       - Left: "MSC Matrix - {band} {phase}"
       - Right: "MSC Network - {band} {phase}"
       - Add info text: sparsify method, n_surrogates if applicable

    7. Save figure

    Cache filename patterns:
    - Dense: {band}_{phase}_msc_sparsify-none_nperseg-{nperseg}.npy
    - Validated: {band}_{phase}_msc_sparsify-soft_nsurr-{n_surrogates}_nperseg-{nperseg}.npy

    Code pattern:
        from lrg_eegfc import load_msc_matrix, get_msc_cache_path

        # Load MSC matrix
        cache_path = get_msc_cache_path(patient, phase, band, cache_root,
                                        sparsify=sparsify, nperseg=nperseg)
        if sparsify == "soft":
            # Adjust cache path for n_surrogates
            cache_path = cache_root / patient / f"{band}_{phase}_msc_sparsify-soft_nsurr-{n_surrogates}_nperseg-{nperseg}.npy"

        if not cache_path.exists():
            raise FileNotFoundError(f"MSC matrix not found: {cache_path}")

        msc_matrix = np.load(cache_path)

        # Load channel labels
        try:
            label_data = scipy.io.loadmat(dataset_root / patient / "channel_labels.mat")
            channel_labels = [str(label[0]) for label in label_data['channel_labels'].flatten()]
        except:
            channel_labels = [str(i) for i in range(msc_matrix.shape[0])]

        # Create figure
        fig, ax = plt.subplots(1, 2, figsize=figsize)

        # Left: Heatmap
        im = ax[0].imshow(msc_matrix, cmap='viridis', vmin=0, vmax=1,
                          interpolation='none', aspect='auto')
        divider = make_axes_locatable(ax[0])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
        ax[0].set_title(f"MSC Matrix - {band} {phase}")

        # Right: Network
        G = nx.from_numpy_array(msc_matrix)
        widths = [G[u][v]['weight'] for u, v in G.edges()]
        if len(widths) > 0:
            wmin, wmax = min(widths), max(widths)
            widths_scaled = [0.1 + 2.9 * (w - wmin) / (wmax - wmin + 1e-10) for w in widths]
        else:
            widths_scaled = []

        pos = nx.spring_layout(G, seed=42)
        label_dict = {i: channel_labels[i] for i in range(len(channel_labels))}

        nx.draw(G, pos=pos, ax=ax[1],
                width=widths_scaled,
                node_size=100,
                node_color='lightblue',
                edge_color='gray',
                with_labels=(len(channel_labels) <= 30),
                labels=label_dict if len(channel_labels) <= 30 else {},
                font_size=8)

        ax[1].set_title(f"MSC Network - {band} {phase}")

        # Save
        if output_path is None:
            output_dir = Path("data/figures/msc") / patient
            output_dir.mkdir(parents=True, exist_ok=True)
            suffix = "dense" if sparsify == "none" else f"validated_nsurr{n_surrogates}"
            output_path = output_dir / f"{band}_{phase}_msc_{suffix}_network.png"

        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(fig)

        return output_path
    """
```

##### Function 3: plot_msc_comparison_dense_vs_validated()

```python
def plot_msc_comparison_dense_vs_validated(
    patient: str,
    phase: str,
    band: str,
    n_surrogates: int = 200,
    cache_root: Path = Path("data/msc_cache"),
    output_path: Optional[Path] = None,
    figsize: tuple = (18, 6),
) -> Path:
    """Compare dense vs validated MSC side-by-side.

    3-panel layout:
    Panel 0: Dense MSC matrix
    Panel 1: Validated MSC matrix
    Panel 2: Difference (dense - validated) showing removed connections

    Implementation:
    1. Load both versions of MSC
    2. Create 3-panel figure
    3. Plot each with appropriate colormaps
    4. Add statistics annotations:
       - Mean MSC (dense)
       - Mean MSC (validated)
       - % edges retained
       - Mean threshold from surrogate analysis

    This is useful for understanding what surrogate validation removes.
    """
```

#### File 2: `src/visualize_msc.py`

**Complete CLI script specification:**

```python
#!/usr/bin/env python3
"""Visualize MSC-based functional connectivity networks.

This script generates visualizations for MSC (Magnitude-Squared Coherence)
matrices without requiring recomputation. Reads from cached MSC matrices.

Usage:
    # Single patient, single band, single phase
    python src/visualize_msc.py --patient Pat_02 --phase rsPre --band beta

    # All plot types for a patient/phase/band
    python src/visualize_msc.py --patient Pat_02 --phase rsPre --band beta --plot-type all

    # Compare dense vs validated MSC
    python src/visualize_msc.py --patient Pat_02 --phase rsPre --band beta \\
        --plot-type comparison --sparsify soft --n-surrogates 200

    # Batch mode: all bands for a patient/phase
    python src/visualize_msc.py --patient Pat_02 --phase rsPre --batch
"""

import argparse
from pathlib import Path
from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS
from lrg_eegfc.visuals.msc import (
    plot_msc_heatmap,
    plot_msc_and_network,
    plot_msc_comparison_dense_vs_validated,
)

def main():
    parser = argparse.ArgumentParser(
        description="Visualize MSC functional connectivity networks"
    )

    # Required arguments
    parser.add_argument("--patient", required=True, help="Patient ID")
    parser.add_argument("--phase", help="Recording phase")
    parser.add_argument("--band", help="Frequency band")

    # Plot type selection
    parser.add_argument(
        "--plot-type",
        choices=["matrix", "network", "comparison", "all"],
        default="network",
        help="Type of plot to generate"
    )

    # MSC version selection
    parser.add_argument(
        "--sparsify",
        choices=["none", "soft"],
        default="none",
        help="MSC version: none (dense) or soft (validated)"
    )
    parser.add_argument(
        "--n-surrogates",
        type=int,
        default=200,
        help="Number of surrogates (if sparsify='soft')"
    )
    parser.add_argument(
        "--nperseg",
        type=int,
        default=256,
        help="Window length for MSC computation"
    )

    # Cache and output
    parser.add_argument(
        "--cache-root",
        type=Path,
        default=Path("data/msc_cache"),
        help="MSC cache directory"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/figures/msc"),
        help="Output directory for figures"
    )

    # Batch mode
    parser.add_argument(
        "--batch",
        action="store_true",
        help="Generate plots for all bands (requires --phase)"
    )

    # Verbosity
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    # Validation
    if args.batch and not args.phase:
        parser.error("--batch requires --phase")

    if not args.batch:
        if not args.phase or not args.band:
            parser.error("--phase and --band required unless using --batch")

    # Determine combinations to process
    if args.batch:
        bands = list(BRAIN_BANDS.keys())
        phases = [args.phase]
    else:
        bands = [args.band]
        phases = [args.phase]

    # Process each combination
    for band in bands:
        for phase in phases:
            if args.verbose:
                print(f"Processing {args.patient} {phase} {band}...")

            output_dir = args.output_dir / args.patient
            output_dir.mkdir(parents=True, exist_ok=True)

            try:
                if args.plot_type in ["matrix", "all"]:
                    # Just heatmap
                    # (Would need to load matrix first, then call plot_msc_heatmap)
                    pass

                if args.plot_type in ["network", "all"]:
                    # Side-by-side matrix + network
                    output_path = plot_msc_and_network(
                        args.patient, phase, band,
                        cache_root=args.cache_root,
                        sparsify=args.sparsify,
                        n_surrogates=args.n_surrogates,
                        nperseg=args.nperseg,
                    )
                    if args.verbose:
                        print(f"  ✓ Saved: {output_path}")

                if args.plot_type == "comparison":
                    # Dense vs validated comparison
                    output_path = plot_msc_comparison_dense_vs_validated(
                        args.patient, phase, band,
                        n_surrogates=args.n_surrogates,
                        cache_root=args.cache_root,
                    )
                    if args.verbose:
                        print(f"  ✓ Saved: {output_path}")

            except Exception as e:
                print(f"  ✗ Failed {args.patient} {phase} {band}: {e}")

    print("\n✓ MSC visualization complete!")

if __name__ == "__main__":
    main()
```

### Output Organization

```
data/figures/msc/
├── Pat_02/
│   ├── beta_rsPre_msc_dense_network.png
│   ├── beta_rsPre_msc_validated_nsurr200_network.png
│   ├── beta_rsPre_msc_comparison.png
│   ├── beta_taskLearn_msc_dense_network.png
│   └── ...
├── Pat_03/
│   └── ...
```

### Testing Checklist for Phase 3

- [ ] Load dense MSC matrix successfully
- [ ] Load validated MSC matrix successfully
- [ ] Handle missing cache files gracefully
- [ ] Channel labels display correctly (or fallback to numbers)
- [ ] Edge widths scale appropriately
- [ ] Network layout is readable
- [ ] Colorbar shows correct range [0, 1]
- [ ] Output files created in correct directory structure
- [ ] Batch mode processes all bands
- [ ] Comparison plot shows meaningful differences

---

## PHASE 4: LRG Visualization Suite (DETAILED)

### Goal
Create comprehensive LRG visualization (entropy, dendrograms, ultrametrics)

### Background: What is LRG Analysis?

**LRG (Laplacian Renormalization Group)** extracts hierarchical structure from networks:
1. Compute graph Laplacian L = D - A
2. Compute resistance distances via matrix exponential: ρ = exp(-τL) / tr(exp(-τL))
3. Invert to get ultrametric distances: d[i,j] = 1/ρ[i,j]
4. Hierarchical clustering on ultrametric distances
5. Optimal threshold selection via partition stability index
6. Entropy analysis across scales

**Cached LRG fields** (from workflow_lrg.py):
- `ultrametric_matrix`: Condensed distance matrix (n*(n-1)/2,)
- `linkage_matrix`: Hierarchical clustering linkage (n-1, 4)
- `entropy_tau`: τ values for entropy computation (n_steps,)
- `entropy_1_minus_S`: Normalized entropy 1-S (n_steps,)
- `entropy_C`: Spectral complexity C (n_steps,)
- `optimal_threshold`: Dendrogram cutting threshold (float)

### Implementation Details

#### File 1: `src/lrg_eegfc/visuals/lrg.py`

##### Function 1: plot_lrg_entropy_curves()

```python
def plot_lrg_entropy_curves(
    patient: str,
    phase: str,
    band: str,
    fc_method: str,  # "corr" or "msc"
    cache_root: Path = Path("data/lrg_cache"),
    output_path: Optional[Path] = None,
    figsize: tuple = (10, 6),
) -> Path:
    """Plot entropy curves (1-S) and spectral complexity (C) vs tau.

    Pattern from: src/lrg_eegfc/plotting.py plot_entropy()

    Implementation:

    1. Load LRG result from cache
       from lrg_eegfc import load_lrg_result
       result = load_lrg_result(patient, phase, band, fc_method, cache_root)

    2. Extract entropy data:
       tau = result.entropy_tau
       S_norm = result.entropy_1_minus_S  # 1-S (normalized entropy)
       C = result.entropy_C  # Spectral complexity

    3. Create figure with log-scale x-axis:
       fig, ax = plt.subplots(figsize=figsize)
       ax.set_xscale('log')

    4. Plot both curves:
       ax.plot(tau, S_norm, label='1-S (Normalized Entropy)', color='blue', lw=2)
       ax.plot(tau, C, label='C (Spectral Complexity)', color='red', lw=2)

    5. Formatting:
       ax.set_xlabel(r'$\\tau$ (Scale Parameter)', fontsize=12)
       ax.set_ylabel('Entropy Metrics', fontsize=12)
       ax.set_title(f'LRG Entropy - {patient} {phase} {band} ({fc_method})', fontsize=14)
       ax.legend(fontsize=10)
       ax.grid(True, alpha=0.3)
       ax.set_ylim(0, 1)  # Both metrics normalized to [0, 1]

    6. Save and return path

    Interpretation:
    - 1-S increases as network becomes more hierarchical
    - C measures complexity of hierarchical organization
    - Crossings and plateaus indicate scale transitions
    """
```

##### Function 2: plot_lrg_dendrogram()

```python
def plot_lrg_dendrogram(
    patient: str,
    phase: str,
    band: str,
    fc_method: str,
    cache_root: Path = Path("data/lrg_cache"),
    output_path: Optional[Path] = None,
    orientation: str = "top",
    figsize: tuple = (12, 8),
    optimal_leaf_order: bool = True,
    show_labels: bool = True,
    max_labels: int = 50,
) -> Path:
    """Plot hierarchical dendrogram with optimal threshold line.

    Pattern from: ipynb/NEW_distance_of_distances.ipynb

    Implementation:

    1. Load LRG result:
       result = load_lrg_result(patient, phase, band, fc_method, cache_root)
       linkage = result.linkage_matrix
       ultrametric_condensed = result.ultrametric_matrix
       optimal_th = result.optimal_threshold

    2. Optional: Optimal leaf ordering for better visualization
       if optimal_leaf_order:
           from scipy.cluster.hierarchy import optimal_leaf_ordering
           linkage = optimal_leaf_ordering(linkage, ultrametric_condensed)

    3. Load channel labels:
       try:
           # Load from patient data
           labels = load_channel_labels(patient)
       except:
           labels = [str(i) for i in range(n_nodes)]

    4. Create dendrogram:
       from scipy.cluster.hierarchy import dendrogram

       fig, ax = plt.subplots(figsize=figsize)

       dendro = dendrogram(
           linkage,
           ax=ax,
           orientation=orientation,
           labels=labels if (show_labels and len(labels) <= max_labels) else None,
           no_labels=(len(labels) > max_labels),
           color_threshold=optimal_th,
           above_threshold_color='black',  # Clusters above threshold are black
           leaf_font_size=8 if len(labels) <= max_labels else 5,
       )

    5. Add optimal threshold line:
       if orientation == "top":
           ax.axhline(optimal_th, color='blue', linestyle='--', lw=2,
                      label=f'Optimal Threshold = {optimal_th:.3f}')
       elif orientation == "right":
           ax.axvline(optimal_th, color='blue', linestyle='--', lw=2,
                      label=f'Optimal Threshold = {optimal_th:.3f}')

    6. Formatting:
       ax.set_yscale('log')  # Log scale for distances
       if orientation == "top":
           ax.set_ylabel('Ultrametric Distance (log scale)', fontsize=12)
           ax.set_xlabel('Channel Index', fontsize=12)
       else:
           ax.set_xlabel('Ultrametric Distance (log scale)', fontsize=12)
           ax.set_ylabel('Channel Index', fontsize=12)

       ax.set_title(f'LRG Dendrogram - {patient} {phase} {band} ({fc_method})',
                    fontsize=14)
       ax.legend(fontsize=10)

    7. Save and return

    Returned dendrogram dict contains:
    - 'icoord': x-coordinates of dendrogram segments
    - 'dcoord': y-coordinates (distances)
    - 'ivl': leaf labels in plotted order
    - 'leaves': original leaf indices in plotted order
    - 'color_list': colors assigned to each link
    - 'leaves_color_list': colors for leaves

    This can be used for network visualization coloring
    """
```

##### Function 3: plot_ultrametric_heatmap()

```python
def plot_ultrametric_heatmap(
    patient: str,
    phase: str,
    band: str,
    fc_method: str,
    cache_root: Path = Path("data/lrg_cache"),
    output_path: Optional[Path] = None,
    figsize: tuple = (10, 10),
    cmap: str = "viridis",
) -> Path:
    """Plot ultrametric distance matrix as heatmap.

    Implementation:

    1. Load LRG result:
       result = load_lrg_result(patient, phase, band, fc_method, cache_root)
       ultrametric_condensed = result.ultrametric_matrix

    2. Convert condensed to square form:
       from scipy.spatial.distance import squareform
       ultrametric_square = squareform(ultrametric_condensed)

    3. Load channel labels:
       labels = load_channel_labels(patient) or range(ultrametric_square.shape[0])

    4. Plot symmetric heatmap:
       fig, ax = plt.subplots(figsize=figsize)

       im = ax.imshow(ultrametric_square, cmap=cmap, interpolation='none',
                      aspect='auto')

       # Colorbar
       from mpl_toolkits.axes_grid1 import make_axes_locatable
       divider = make_axes_locatable(ax)
       cax = divider.append_axes("right", size="5%", pad=0.05)
       cbar = plt.colorbar(im, cax=cax)
       cbar.set_label('Ultrametric Distance', fontsize=12)

    5. Add labels (if not too many):
       if len(labels) <= 50:
           ax.set_xticks(range(len(labels)))
           ax.set_yticks(range(len(labels)))
           ax.set_xticklabels(labels, rotation=90, fontsize=8)
           ax.set_yticklabels(labels, fontsize=8)
       else:
           ax.set_xticks([])
           ax.set_yticks([])

    6. Title and save:
       ax.set_title(f'Ultrametric Distance - {patient} {phase} {band} ({fc_method})',
                    fontsize=14)

       fig.savefig(output_path, dpi=300, bbox_inches='tight')
       plt.close(fig)
       return output_path

    Interpretation:
    - Diagonal is zero (distance to self)
    - Symmetric matrix
    - Ultrametric property: d(i,k) <= max(d(i,j), d(j,k)) for all i,j,k
    - Darker colors = closer in hierarchy
    - Block diagonal structure indicates hierarchical clusters
    """
```

##### Function 4: plot_lrg_full_panel()

```python
def plot_lrg_full_panel(
    patient: str,
    phase: str,
    band: str,
    fc_method: str,
    cache_root: Path = Path("data/lrg_cache"),
    output_path: Optional[Path] = None,
    figsize: tuple = (20, 14),
) -> Path:
    """4-panel comprehensive LRG visualization.

    Layout:
    ┌─────────────┬─────────────┐
    │  Entropy    │  Dendrogram │
    │  (1-S, C)   │  (log-scale)│
    ├─────────────┼─────────────┤
    │  Ultrametric│  Network    │
    │  Heatmap    │  (colored)  │
    └─────────────┴─────────────┘

    Implementation:

    1. Load all necessary data:
       result = load_lrg_result(...)
       fc_matrix = load FC matrix (corr or msc)
       labels = load_channel_labels(...)

    2. Create 2x2 subplot layout:
       fig, ax = plt.subplots(2, 2, figsize=figsize)

    3. Panel [0, 0]: Entropy curves
       - Plot 1-S and C vs tau (log x-axis)
       - Legend, grid, labels

    4. Panel [0, 1]: Dendrogram
       - Hierarchical dendrogram with optimal threshold
       - Color-coded clusters
       - Log y-axis

    5. Panel [1, 0]: Ultrametric heatmap
       - Square-form ultrametric distance matrix
       - Viridis colormap with colorbar

    6. Panel [1, 1]: Network colored by dendrogram
       - Build graph G from FC matrix
       - Get dendrogram leaf colors
       - Map node indices to dendrogram leaves
       - Color nodes by cluster
       - Edge widths proportional to weights
       - Spring layout or Kamada-Kawai

       Code pattern for coloring:
           dendro = dendrogram(linkage, ax=temp_ax, no_plot=True,
                              color_threshold=optimal_th)

           # Map leaves to colors
           leaf_to_color = {}
           for leaf_idx, color in zip(dendro['leaves'], dendro['leaves_color_list']):
               leaf_to_color[leaf_idx] = color

           # Color nodes
           node_colors = [leaf_to_color.get(i, 'gray') for i in G.nodes()]

           # Draw network
           pos = nx.spring_layout(G, seed=42)
           widths = [G[u][v]['weight'] for u, v in G.edges()]
           nx.draw(G, pos=pos, ax=ax[1,1],
                   node_color=node_colors,
                   width=[w*5 for w in widths],
                   node_size=100,
                   with_labels=False)

    7. Overall title and tight layout:
       fig.suptitle(f'LRG Analysis - {patient} {phase} {band} ({fc_method})',
                    fontsize=16, y=0.98)
       fig.tight_layout(rect=[0, 0, 1, 0.97])

    8. Save high-resolution figure:
       fig.savefig(output_path, dpi=300, bbox_inches='tight')
       plt.close(fig)

    This provides a complete overview of hierarchical organization in one figure.
    """
```

#### File 2: `src/visualize_lrg.py`

**Complete CLI specification:**

```python
#!/usr/bin/env python3
"""Visualize LRG (Laplacian Renormalization Group) analysis results.

This script generates visualizations for LRG hierarchical network analysis
without requiring recomputation. Reads from cached LRG results.

Usage:
    # Entropy curves for one analysis
    python src/visualize_lrg.py --patient Pat_02 --phase rsPre --band beta \\
        --fc-method corr --plot-type entropy

    # Full 4-panel visualization
    python src/visualize_lrg.py --patient Pat_02 --phase rsPre --band beta \\
        --fc-method corr --plot-type full

    # All plot types
    python src/visualize_lrg.py --patient Pat_02 --phase rsPre --band beta \\
        --fc-method corr --plot-type all

    # Batch mode: all bands for a patient/phase
    python src/visualize_lrg.py --patient Pat_02 --phase rsPre --fc-method corr --batch
"""

import argparse
from pathlib import Path
from lrg_eegfc.config.const import BRAIN_BANDS
from lrg_eegfc.visuals.lrg import (
    plot_lrg_entropy_curves,
    plot_lrg_dendrogram,
    plot_ultrametric_heatmap,
    plot_lrg_full_panel,
)

def main():
    parser = argparse.ArgumentParser(
        description="Visualize LRG hierarchical network analysis"
    )

    # Required
    parser.add_argument("--patient", required=True)
    parser.add_argument("--phase", help="Recording phase")
    parser.add_argument("--band", help="Frequency band")
    parser.add_argument("--fc-method", required=True, choices=["corr", "msc"])

    # Plot type
    parser.add_argument(
        "--plot-type",
        choices=["entropy", "dendrogram", "ultrametric", "full", "all"],
        default="full",
        help="Type of plot to generate"
    )

    # Dendrogram options
    parser.add_argument(
        "--orientation",
        choices=["top", "right", "bottom", "left"],
        default="top",
        help="Dendrogram orientation"
    )

    # Paths
    parser.add_argument(
        "--cache-root",
        type=Path,
        default=Path("data/lrg_cache"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/figures/lrg"),
    )

    # Batch mode
    parser.add_argument("--batch", action="store_true")
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    # Validation
    if args.batch and not args.phase:
        parser.error("--batch requires --phase")
    if not args.batch and (not args.phase or not args.band):
        parser.error("--phase and --band required unless using --batch")

    # Determine combinations
    bands = list(BRAIN_BANDS.keys()) if args.batch else [args.band]

    # Process
    for band in bands:
        if args.verbose:
            print(f"Processing {args.patient} {args.phase} {band} ({args.fc_method})...")

        output_dir = args.output_dir / args.patient
        output_dir.mkdir(parents=True, exist_ok=True)

        try:
            if args.plot_type in ["entropy", "all"]:
                output_path = plot_lrg_entropy_curves(
                    args.patient, args.phase, band, args.fc_method,
                    cache_root=args.cache_root
                )
                if args.verbose:
                    print(f"  ✓ Entropy: {output_path}")

            if args.plot_type in ["dendrogram", "all"]:
                output_path = plot_lrg_dendrogram(
                    args.patient, args.phase, band, args.fc_method,
                    cache_root=args.cache_root,
                    orientation=args.orientation
                )
                if args.verbose:
                    print(f"  ✓ Dendrogram: {output_path}")

            if args.plot_type in ["ultrametric", "all"]:
                output_path = plot_ultrametric_heatmap(
                    args.patient, args.phase, band, args.fc_method,
                    cache_root=args.cache_root
                )
                if args.verbose:
                    print(f"  ✓ Ultrametric: {output_path}")

            if args.plot_type in ["full", "all"]:
                output_path = plot_lrg_full_panel(
                    args.patient, args.phase, band, args.fc_method,
                    cache_root=args.cache_root
                )
                if args.verbose:
                    print(f"  ✓ Full panel: {output_path}")

        except Exception as e:
            print(f"  ✗ Failed: {e}")

    print("\n✓ LRG visualization complete!")

if __name__ == "__main__":
    main()
```

### Testing Checklist for Phase 4

- [ ] Load LRG result from cache successfully
- [ ] Entropy curves plotted with correct axes (log x-axis)
- [ ] Dendrogram shows correct hierarchical structure
- [ ] Optimal threshold line visible on dendrogram
- [ ] Ultrametric heatmap symmetric
- [ ] Network colored by dendrogram clusters
- [ ] Full panel layout correct (2x2)
- [ ] All visualizations handle missing labels gracefully
- [ ] Batch mode processes all bands

---

## PHASE 5: Comparison Visualization Suite (DETAILED)

### Goal
Create visualizations for comparing FC methods and phases (reorganization analysis)

### Background: Network Reorganization Analysis

**User requirement:** "correlation between cluster belonging of each node or measure of distance of ultrametric matrices or even topological measure of similarity hierarchical tree to hierarchical tree. all of these have been explored in the notebooks"

**Key concepts:**
1. **MSC vs Correlation comparison**: How do two FC methods produce different hierarchies?
2. **Phase reorganization**: How does network structure change across experimental phases (rsPre → taskLearn → taskTest → rsPost)?
3. **Memory effects**: Persistent changes in network organization after task performance

**Available metrics** (from compare.py):
- matrix_distance (Frobenius norm)
- multiscale_distance
- quantile_rmse
- rank_correlation
- tree_similarity
- scale_profile_distance

### Implementation Details

#### File 1: `src/lrg_eegfc/visuals/comparisons.py`

##### Function 1: compute_cluster_membership_correlation()

```python
def compute_cluster_membership_correlation(
    linkage_1: np.ndarray,
    linkage_2: np.ndarray,
    threshold_1: Optional[float] = None,
    threshold_2: Optional[float] = None,
    method: str = "spearman",
) -> float:
    """Compute correlation of cluster memberships between two hierarchies.

    From user requirement: "correlation between cluster belonging of each node"

    Algorithm:
    1. Cut dendrograms at thresholds to get flat clusters
       - Use scipy.cluster.hierarchy.fcluster()
       - criterion='distance' with given threshold
       - If threshold not provided, use t=0.5*max_distance

    2. Get cluster labels for each node
       - labels_1 = fcluster(linkage_1, threshold_1, criterion='distance')
       - labels_2 = fcluster(linkage_2, threshold_2, criterion='distance')

    3. Compute correlation of label assignments
       - Spearman: scipy.stats.spearmanr(labels_1, labels_2)
       - Pearson: scipy.stats.pearsonr(labels_1, labels_2)
       - Kendall: scipy.stats.kendalltau(labels_1, labels_2)

    Returns: correlation coefficient [-1, 1]
    - 1.0: Identical clustering (perfect agreement)
    - 0.0: No relationship between clusterings
    - -1.0: Inverse clustering

    Implementation:
        from scipy.cluster.hierarchy import fcluster
        from scipy.stats import spearmanr, pearsonr, kendalltau

        # Determine thresholds
        if threshold_1 is None:
            max_dist_1 = linkage_1[:, 2].max()
            threshold_1 = 0.5 * max_dist_1

        if threshold_2 is None:
            max_dist_2 = linkage_2[:, 2].max()
            threshold_2 = 0.5 * max_dist_2

        # Cut dendrograms
        labels_1 = fcluster(linkage_1, threshold_1, criterion='distance')
        labels_2 = fcluster(linkage_2, threshold_2, criterion='distance')

        # Compute correlation
        if method == "spearman":
            corr, pval = spearmanr(labels_1, labels_2)
        elif method == "pearson":
            corr, pval = pearsonr(labels_1, labels_2)
        elif method == "kendall":
            corr, pval = kendalltau(labels_1, labels_2)
        else:
            raise ValueError(f"Unknown method: {method}")

        return corr

    Note: This assumes both dendrograms have the same number of leaves
    """
```

##### Function 2: plot_msc_vs_corr_comparison()

```python
def plot_msc_vs_corr_comparison(
    patient: str,
    phase: str,
    band: str,
    cache_root: Path = Path("data/lrg_cache"),
    output_path: Optional[Path] = None,
    figsize: tuple = (20, 12),
) -> Path:
    """Side-by-side comparison of MSC and Correlation LRG results.

    2x3 grid layout:
    ┌─────────────┬─────────────┬─────────────┐
    │ MSC Entropy │ MSC Dendro  │ MSC Ultra   │
    ├─────────────┼─────────────┼─────────────┤
    │ Corr Entropy│ Corr Dendro │ Corr Ultra  │
    └─────────────┴─────────────┴─────────────┘

    Implementation:

    1. Load both LRG results:
       result_msc = load_lrg_result(patient, phase, band, "msc", cache_root)
       result_corr = load_lrg_result(patient, phase, band, "corr", cache_root)

    2. Create 2x3 subplot layout:
       fig, ax = plt.subplots(2, 3, figsize=figsize)

    3. Row 0: MSC visualizations
       ax[0, 0]: Entropy curves (1-S and C vs tau)
       ax[0, 1]: Dendrogram with optimal threshold
       ax[0, 2]: Ultrametric heatmap

    4. Row 1: Correlation visualizations
       ax[1, 0]: Entropy curves
       ax[1, 1]: Dendrogram
       ax[1, 2]: Ultrametric heatmap

    5. Compute comparison metrics:
       - Load ultrametric matrices from both
       - Compute matrix distance, tree similarity, cluster correlation
       - Add text annotations showing metrics

    6. Add row labels:
       ax[0, 0].text(-0.15, 0.5, 'MSC', transform=ax[0, 0].transAxes,
                     fontsize=16, fontweight='bold', va='center', rotation=90)
       ax[1, 0].text(-0.15, 0.5, 'Correlation', transform=ax[1, 0].transAxes,
                     fontsize=16, fontweight='bold', va='center', rotation=90)

    7. Overall title and metrics:
       fig.suptitle(f'MSC vs Correlation - {patient} {phase} {band}',
                    fontsize=18)

       # Add comparison metrics textbox
       from lrg_eegfc.compare import compare_ultrametric_matrices
       comparison = compare_ultrametric_matrices(
           result_msc.ultrametric_matrix,
           result_corr.ultrametric_matrix,
           result_msc.linkage_matrix,
           result_corr.linkage_matrix
       )

       metrics_text = f'''Comparison Metrics:
       Matrix Distance: {comparison.matrix_distance:.4f}
       Tree Similarity: {comparison.tree_similarity:.4f}
       Rank Correlation: {comparison.rank_correlation:.4f}
       Multiscale Distance: {comparison.multiscale_distance:.4f}
       '''

       fig.text(0.02, 0.02, metrics_text, fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    8. Save high-resolution figure

    This allows direct visual comparison of network structures
    derived from two different FC methods.
    """
```

##### Function 3: plot_phase_reorganization()

```python
def plot_phase_reorganization(
    patient: str,
    band: str,
    fc_method: str,
    phases: List[str] = ["rsPre", "taskLearn", "taskTest", "rsPost"],
    cache_root: Path = Path("data/lrg_cache"),
    output_path: Optional[Path] = None,
    figsize: tuple = (24, 16),
) -> Path:
    """Visualize network reorganization across experimental phases.

    Pattern from: ipynb/NEW_distance_of_distances.ipynb

    3-row layout:
    ┌────────┬────────┬────────┬────────┐
    │ rsPre  │taskLrn │taskTst │rsPost  │  Row 0: Dendrograms
    ├────────┼────────┼────────┼────────┤
    │ rsPre  │taskLrn │taskTst │rsPost  │  Row 1: Ultrametric heatmaps
    ├────────┴────────┴────────┴────────┤
    │  Cluster Membership Correlation   │  Row 2: Phase-to-phase correlation matrix
    └───────────────────────────────────┘

    Implementation:

    1. Load all LRG results for specified phases:
       results = {}
       for phase in phases:
           results[phase] = load_lrg_result(patient, phase, band, fc_method, cache_root)

    2. Create 3-row subplot layout:
       n_phases = len(phases)
       fig = plt.figure(figsize=figsize)
       gs = fig.add_gridspec(3, n_phases, height_ratios=[2, 2, 1.5])

       # Row 0: Dendrograms
       ax_dendro = [fig.add_subplot(gs[0, i]) for i in range(n_phases)]

       # Row 1: Ultrametric heatmaps
       ax_ultra = [fig.add_subplot(gs[1, i]) for i in range(n_phases)]

       # Row 2: Correlation matrix (spans all columns)
       ax_corr = fig.add_subplot(gs[2, :])

    3. Row 0: Side-by-side dendrograms
       - Plot all dendrograms with SAME y-axis scale for comparison
       - Use optimal leaf ordering for each
       - Color-code by phase (different colors for different phases)
       - Mark optimal thresholds

       Code pattern:
           max_height = max(results[p].linkage_matrix[:, 2].max() for p in phases)

           for i, phase in enumerate(phases):
               linkage = results[phase].linkage_matrix
               optimal_th = results[phase].optimal_threshold

               dendro = dendrogram(linkage, ax=ax_dendro[i],
                                  no_labels=True,
                                  color_threshold=optimal_th)

               ax_dendro[i].set_title(phase, fontsize=14)
               ax_dendro[i].set_ylim(0, max_height * 1.1)  # Same scale
               ax_dendro[i].set_yscale('log')
               ax_dendro[i].axhline(optimal_th, color='b', linestyle='--', lw=2)

    4. Row 1: Ultrametric heatmaps for each phase
       - Convert condensed to square form
       - Plot with same colormap and scale
       - Shared colorbar

       Code pattern:
           from scipy.spatial.distance import squareform

           # Find global min/max for consistent colorscale
           all_ultras = [squareform(results[p].ultrametric_matrix) for p in phases]
           vmin = min(u.min() for u in all_ultras)
           vmax = max(u.max() for u in all_ultras)

           for i, phase in enumerate(phases):
               ultra = all_ultras[i]
               im = ax_ultra[i].imshow(ultra, cmap='viridis', vmin=vmin, vmax=vmax,
                                       interpolation='none', aspect='auto')
               ax_ultra[i].set_title(phase, fontsize=12)
               ax_ultra[i].set_xticks([])
               ax_ultra[i].set_yticks([])

           # Shared colorbar
           fig.colorbar(im, ax=ax_ultra, orientation='horizontal',
                       fraction=0.05, pad=0.05, label='Ultrametric Distance')

    5. Row 2: Cluster membership correlation matrix across phases
       - Compute pairwise correlations between all phase pairs
       - Create symmetric correlation matrix
       - Plot as heatmap with annotations

       Code pattern:
           n = len(phases)
           corr_matrix = np.zeros((n, n))

           for i, phase1 in enumerate(phases):
               for j, phase2 in enumerate(phases):
                   if i == j:
                       corr_matrix[i, j] = 1.0
                   else:
                       corr = compute_cluster_membership_correlation(
                           results[phase1].linkage_matrix,
                           results[phase2].linkage_matrix,
                           results[phase1].optimal_threshold,
                           results[phase2].optimal_threshold
                       )
                       corr_matrix[i, j] = corr

           # Plot correlation matrix
           im = ax_corr.imshow(corr_matrix, cmap='RdBu_r', vmin=-1, vmax=1,
                               interpolation='none', aspect='auto')

           ax_corr.set_xticks(range(n))
           ax_corr.set_yticks(range(n))
           ax_corr.set_xticklabels(phases, fontsize=12)
           ax_corr.set_yticklabels(phases, fontsize=12)

           # Annotate with correlation values
           for i in range(n):
               for j in range(n):
                   text = ax_corr.text(j, i, f'{corr_matrix[i, j]:.2f}',
                                      ha="center", va="center", color="black",
                                      fontsize=10)

           ax_corr.set_title('Cluster Membership Correlation Across Phases', fontsize=14)
           fig.colorbar(im, ax=ax_corr, orientation='horizontal', pad=0.1)

    6. Overall title and layout:
       fig.suptitle(f'Phase Reorganization - {patient} {band} ({fc_method})',
                    fontsize=18, y=0.98)
       fig.tight_layout(rect=[0, 0, 1, 0.97])

    7. Save figure

    Interpretation:
    - Dendrogram changes show structural reorganization
    - Ultrametric changes quantify distance modifications
    - Correlation matrix shows stability/change in cluster assignments
    - High correlation = stable clustering across phases
    - Low correlation = significant reorganization
    """
```

##### Function 4: plot_reorganization_metrics()

```python
def plot_reorganization_metrics(
    patient: str,
    band: str,
    fc_method: str,
    comparison_csv: Path,  # From compare_fc_methods.py --mode memory output
    output_path: Optional[Path] = None,
    figsize: tuple = (18, 12),
) -> Path:
    """Plot reorganization metrics from comparison analysis.

    From user requirement: "measure of distance of ultrametric matrices or
    even topological measure of similarity hierarchical tree to hierarchical tree"

    4-panel layout showing metric evolution across phase pairs:
    ┌─────────────┬─────────────┐
    │ Matrix Dist │ Tree Sim    │
    ├─────────────┼─────────────┤
    │ Multiscale  │ Rank Corr   │
    └─────────────┴─────────────┘

    Implementation:

    1. Load comparison results from CSV:
       import pandas as pd
       df = pd.read_csv(comparison_csv)

       # Filter for this patient/band/fc_method
       df_filtered = df[
           (df['patient'] == patient) &
           (df['band'] == band) &
           (df['fc_method'] == fc_method)
       ]

    2. Extract phase pairs and metrics:
       phase_pairs = df_filtered['phase_pair'].values
       # Typical pairs: "rsPre-rsPost", "rsPre-taskLearn", etc.

       metrics = {
           'matrix_distance': df_filtered['matrix_distance'].values,
           'tree_similarity': df_filtered['tree_similarity'].values,
           'multiscale_distance': df_filtered['multiscale_distance'].values,
           'rank_correlation': df_filtered['rank_correlation'].values,
       }

    3. Create 2x2 subplot layout:
       fig, ax = plt.subplots(2, 2, figsize=figsize)

    4. Panel [0, 0]: Matrix Distance evolution
       - Bar plot or line plot showing matrix distance for each phase pair
       - Higher values = more reorganization

       ax[0, 0].bar(range(len(phase_pairs)), metrics['matrix_distance'])
       ax[0, 0].set_xticks(range(len(phase_pairs)))
       ax[0, 0].set_xticklabels(phase_pairs, rotation=45, ha='right')
       ax[0, 0].set_ylabel('Matrix Distance (Frobenius)')
       ax[0, 0].set_title('Ultrametric Matrix Distance')
       ax[0, 0].grid(True, alpha=0.3)

    5. Panel [0, 1]: Tree Similarity
       - Higher = more similar tree structure
       - Inverted from distance metrics

       ax[0, 1].bar(range(len(phase_pairs)), metrics['tree_similarity'],
                    color='green')
       ax[0, 1].set_xticks(range(len(phase_pairs)))
       ax[0, 1].set_xticklabels(phase_pairs, rotation=45, ha='right')
       ax[0, 1].set_ylabel('Tree Similarity')
       ax[0, 1].set_title('Hierarchical Tree Topology Similarity')
       ax[0, 1].set_ylim(0, 1)
       ax[0, 1].grid(True, alpha=0.3)

    6. Panel [1, 0]: Multiscale Distance
       - Quantifies changes across hierarchical scales

       ax[1, 0].bar(range(len(phase_pairs)), metrics['multiscale_distance'],
                    color='orange')
       ax[1, 0].set_xticks(range(len(phase_pairs)))
       ax[1, 0].set_xticklabels(phase_pairs, rotation=45, ha='right')
       ax[1, 0].set_ylabel('Multiscale Distance')
       ax[1, 0].set_title('Multiscale Reorganization')
       ax[1, 0].grid(True, alpha=0.3)

    7. Panel [1, 1]: Rank Correlation
       - Preservation of distance ordering
       - High correlation = preserved ranking

       ax[1, 1].bar(range(len(phase_pairs)), metrics['rank_correlation'],
                    color='purple')
       ax[1, 1].set_xticks(range(len(phase_pairs)))
       ax[1, 1].set_xticklabels(phase_pairs, rotation=45, ha='right')
       ax[1, 1].set_ylabel('Rank Correlation')
       ax[1, 1].set_title('Distance Ranking Preservation')
       ax[1, 1].set_ylim(-1, 1)
       ax[1, 1].axhline(0, color='black', linestyle='--', lw=1)
       ax[1, 1].grid(True, alpha=0.3)

    8. Add classification of reorganization strength:
       - Compute mean of all distance metrics
       - Classify as: weak, moderate, strong reorganization
       - Add text annotation

       mean_dist = (metrics['matrix_distance'].mean() +
                    metrics['multiscale_distance'].mean() +
                    (1 - metrics['tree_similarity'].mean())) / 3

       if mean_dist < 0.3:
           strength = "Weak"
           color = "green"
       elif mean_dist < 0.6:
           strength = "Moderate"
           color = "orange"
       else:
           strength = "Strong"
           color = "red"

       fig.text(0.5, 0.02, f'Overall Reorganization Strength: {strength}',
                ha='center', fontsize=14, fontweight='bold',
                bbox=dict(boxstyle='round', facecolor=color, alpha=0.3))

    9. Overall title and save:
       fig.suptitle(f'Reorganization Metrics - {patient} {band} ({fc_method})',
                    fontsize=16)
       fig.tight_layout(rect=[0, 0.05, 1, 0.95])

    This provides quantitative assessment of network reorganization
    across experimental phases.
    """
```

#### File 2: `src/visualize_comparisons.py`

**Complete CLI specification:**

```python
#!/usr/bin/env python3
"""Visualize comparisons between FC methods and across phases.

This script generates comparison visualizations:
1. MSC vs Correlation (different FC methods)
2. Phase reorganization (network changes across experimental phases)
3. Reorganization metrics (quantitative analysis)

Usage:
    # Compare MSC vs Correlation for one case
    python src/visualize_comparisons.py --patient Pat_02 --phase rsPre --band beta \\
        --comparison-type msc-vs-corr

    # Phase reorganization analysis
    python src/visualize_comparisons.py --patient Pat_02 --band beta \\
        --fc-method corr --comparison-type phase-reorganization

    # Reorganization metrics from comparison CSV
    python src/visualize_comparisons.py --patient Pat_02 --band beta \\
        --fc-method corr --comparison-type metrics \\
        --comparison-csv results/memory_corr.csv
"""

import argparse
from pathlib import Path
from lrg_eegfc.config.const import BRAIN_BANDS
from lrg_eegfc.visuals.comparisons import (
    plot_msc_vs_corr_comparison,
    plot_phase_reorganization,
    plot_reorganization_metrics,
)

def main():
    parser = argparse.ArgumentParser(
        description="Visualize FC method and phase comparisons"
    )

    # Required
    parser.add_argument("--patient", required=True)
    parser.add_argument("--band", help="Frequency band (required for most types)")

    # Comparison type
    parser.add_argument(
        "--comparison-type",
        required=True,
        choices=["msc-vs-corr", "phase-reorganization", "metrics"],
        help="Type of comparison visualization"
    )

    # Type-specific args
    parser.add_argument("--phase", help="Phase (required for msc-vs-corr)")
    parser.add_argument(
        "--fc-method",
        choices=["corr", "msc"],
        help="FC method (required for phase-reorganization and metrics)"
    )
    parser.add_argument(
        "--comparison-csv",
        type=Path,
        help="Comparison CSV file (required for metrics type)"
    )

    # Paths
    parser.add_argument(
        "--cache-root",
        type=Path,
        default=Path("data/lrg_cache"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/figures/comparisons"),
    )

    # Options
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    # Validation
    if args.comparison_type == "msc-vs-corr":
        if not args.phase or not args.band:
            parser.error("--phase and --band required for msc-vs-corr")

    elif args.comparison_type == "phase-reorganization":
        if not args.band or not args.fc_method:
            parser.error("--band and --fc-method required for phase-reorganization")

    elif args.comparison_type == "metrics":
        if not args.band or not args.fc_method or not args.comparison_csv:
            parser.error("--band, --fc-method, and --comparison-csv required for metrics")

    # Create output directory
    output_dir = args.output_dir / args.patient
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate visualization
    try:
        if args.comparison_type == "msc-vs-corr":
            if args.verbose:
                print(f"Comparing MSC vs Correlation for {args.patient} {args.phase} {args.band}...")

            output_path = plot_msc_vs_corr_comparison(
                args.patient, args.phase, args.band,
                cache_root=args.cache_root
            )
            print(f"✓ Saved: {output_path}")

        elif args.comparison_type == "phase-reorganization":
            if args.verbose:
                print(f"Analyzing phase reorganization for {args.patient} {args.band} ({args.fc_method})...")

            output_path = plot_phase_reorganization(
                args.patient, args.band, args.fc_method,
                cache_root=args.cache_root
            )
            print(f"✓ Saved: {output_path}")

        elif args.comparison_type == "metrics":
            if args.verbose:
                print(f"Plotting reorganization metrics for {args.patient} {args.band} ({args.fc_method})...")

            output_path = plot_reorganization_metrics(
                args.patient, args.band, args.fc_method,
                comparison_csv=args.comparison_csv
            )
            print(f"✓ Saved: {output_path}")

    except Exception as e:
        print(f"✗ Failed: {e}")
        raise

    print("\n✓ Comparison visualization complete!")

if __name__ == "__main__":
    main()
```

### Testing Checklist for Phase 5

- [ ] Cluster membership correlation computed correctly
- [ ] MSC vs Corr comparison shows both methods side-by-side
- [ ] Phase reorganization dendrograms aligned vertically
- [ ] Ultrametric heatmaps use consistent colorscale
- [ ] Cluster correlation matrix symmetric and properly annotated
- [ ] Reorganization metrics load from CSV correctly
- [ ] All metrics plotted with appropriate scales
- [ ] Reorganization strength classification sensible
- [ ] All comparison types handle missing data gracefully

---

## PHASE 6: Documentation and Integration (DETAILED)

### Goal
Create comprehensive VISUALS.md guide and integrate with existing pipeline

### VISUALS.md Structure (Complete)