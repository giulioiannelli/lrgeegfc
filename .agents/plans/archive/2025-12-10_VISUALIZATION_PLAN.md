---
name: visualization-plan
type: plan
era: MSC
status: superseded
created: 2025-12-10
updated: 2026-04-24
pointers: []
---

# Visualization and Cleaning Pipeline Implementation Plan

## Overview

Implement a comprehensive visualization and cleaning pipeline for functional connectivity analysis. This includes:
1. Correlation matrix cleaning (Marchenko-Pastur spectral filtering + percolation thresholding)
2. Visualization suite (matrices, networks, LRG analysis, comparisons)
3. Structured workflow with separate visualization scripts
4. Comprehensive documentation

## User Requirements (Original Request)

> "ok for the analysis but now i need visuals, like for the correlation matrices and so on. so in order:
> 0. do the program for the pipeline for cleaning matrices (you can find code in the notebooks) correlation by first removing the spectral noise (marchenko pastur method) and then thresholding up to just before the first node gets detached. create the program then maybe also store the output in the casched as cleaned matrix.
> 1. for the msc there is no need of thresholding or removing spectral noise then here we dont need anything.
> 3. create visuals: the logic of the visuals is "not having to compute data" so they have to leverage dthe fils coming out from the analysis pipeline. for the Msc and correlation plot side by side the correlation matrix and the result netowork (edges proportional to the weight) while only for correlation matrix analysis also show the spectral cleaning of the marchenko pastur (i.e. the comparison of the two distributions) and the thresholding analysis i.e. P_inf and E_inf vs \theta (where \theta is threshold)
> 4. create a "visuals.md" similar to what you did for the analysis pipeline to keep everithing under track.
> 5. now the same goes on with the LRG pipeline and then with the comparison of topologies
> 6. i think these are a lot of things so maybe create and structure a plan for this.
> 7. if you think to other fancy plots that we can do for MSC or kind of analysis feel free to suggests, but add them at the end of the whole implementation"

### Cleaning Pipeline Requirements
1. **Marchenko-Pastur spectral cleaning** - Remove noise from correlation eigenvalues
2. **Percolation thresholding** - Threshold just before first node detaches from giant component
3. **Caching strategy** - Store cleaned matrices with `_cleaned.npy` suffix in same cache directory
4. **No MSC cleaning** - MSC matrices don't need spectral cleaning or thresholding

### Visualization Requirements
1. **Separate visualization scripts** - Keep computation and plotting separate ("not having to compute data")
2. **Leverage cached data** - All plots read from existing analysis outputs
3. **For Correlation**:
   - Side-by-side: correlation matrix + network graph (edge weights proportional)
   - Marchenko-Pastur comparison (eigenvalue histogram with MP curve overlay)
   - Percolation analysis (P_inf and E_inf vs threshold θ)
4. **For MSC**:
   - Side-by-side: MSC matrix + network graph (edge weights proportional)
   - No spectral cleaning or percolation plots
5. **For LRG**:
   - Entropy curves (1-S and C vs tau)
   - Dendrograms with optimal threshold
   - Ultrametric distance heatmaps
   - Phase comparison (reorganization analysis: cluster correlation, ultrametric distance, tree topology)
6. **Documentation** - Create VISUALS.md guide

## Critical Files Identified

### Existing Code to Reuse
- `src/lrg_eegfc/plotting.py` - plot_correlation_matrix, plot_entropy, plot_dendrogram, plot_graph
- `src/lrg_eegfc/utils/corrmat/base.py` - clean_correlation_matrix()
- `src/lrg_eegfc/utils/corrmat/thresholds.py` - find_exact_detachment_threshold(), find_threshold_jumps()
- `lrgsglib/.../thresholding.py` - compute_threshold_stats(), compute_threshold_stats_fast()
- `lrgsglib/.../probability.py` - marchenko_pastur()
- `lrgsglib/.../clustering.py` - compute_normalized_linkage(), compute_optimal_threshold()

### Notebook Patterns to Follow

**Percolation Curves (P_inf, E_inf vs θ):**
- `ipynb/poster01.ipynb` (cell ad38bec7) - Complete percolation visualization with jump annotations

**Marchenko-Pastur Visualization:**
- `ipynb/pat02_corrnet_bands.ipynb` (cell 49e51f2b) - 3-panel: original matrix | cleaned matrix | eigenvalue comparison

**Side-by-Side Matrix + Network:**
- `ipynb/pat02_corrnet_bands.ipynb` (cell f9bd5d3c) - 3-row layout with correlation matrices, networks, giant components

**Complete LRG + Comparison Workflow:**
- `ipynb/NEW_distance_of_distances.ipynb` - Full dendrogram, ultrametric, and phase comparison workflow

---

## PHASE 1: Correlation Cleaning Workflow ✅ COMPLETE

### Goal
Create workflow for cleaning correlation matrices with Marchenko-Pastur + percolation thresholding

### Implementation

**Create:** `src/lrg_eegfc/workflow_cleaning.py` ✅

**Key functions:**
```python
@dataclass
class CleanedCorrResult:
    """Result from correlation cleaning workflow."""
    cleaned_matrix: np.ndarray  # Cleaned correlation matrix
    original_matrix: np.ndarray  # Original for comparison
    threshold: float  # Percolation threshold (just before detachment)
    eigenvalues: np.ndarray
    lambda_min: float  # MP lower bound
    lambda_max: float  # MP upper bound
    signal_eigenvalues: np.ndarray
    patient: str
    phase: str
    band: str
    n_signal_components: int
    n_channels: int
    detachment_info: Dict  # Info about first detachment

def clean_correlation_matrix_full(
    patient: str,
    phase: str,
    band: str,
    dataset_root: Path = Path("data/stereoeeg_patients"),
    cache_root: Path = Path("data/corr_cache"),
    *,
    filter_order: int = 4,
    sample_rate: float = 2048.0,
    use_cache: bool = True,
    overwrite_cache: bool = False,
    verbose: bool = False,
) -> CleanedCorrResult:
    """Complete cleaning workflow: bandpass → correlate → MP clean → threshold.

    Steps:
    1. Load timeseries and bandpass filter
    2. Compute correlation matrix
    3. Apply Marchenko-Pastur spectral cleaning
    4. Find percolation threshold (just before first detachment)
    5. Apply threshold to cleaned matrix
    6. Cache with suffix: {band}_{phase}_corr_cleaned.npy
    """

def compute_cleaned_corr_for_patient(patient, **kwargs):
    """Batch process all bands/phases for a patient."""
```

**Caching:**
- Cleaned matrices: `data/corr_cache/{patient}/{band}_{phase}_corr_cleaned.npy`
- Metadata: `data/corr_cache/{patient}/{band}_{phase}_corr_cleaned_meta.npz` containing:
  - `threshold`, `eigenvalues`, `lambda_min`, `lambda_max`
  - `signal_eigenvalues`, `n_signal_components`
  - `detachment_info` (percolation analysis details)

**Create:** `src/clean_correlation_matrices.py` ✅
- CLI script to run cleaning pipeline on all patients
- Arguments: `--patients`, `--overwrite`, `--verbose`
- Outputs both cleaned matrices and metadata

**Modify:** `src/lrg_eegfc/__init__.py` ✅
- Added: `from .workflow_cleaning import CleanedCorrResult, clean_correlation_matrix_full, load_cleaned_corr_matrix, get_cleaned_corr_cache_path, compute_cleaned_corr_for_patient`

---

## PHASE 2: Correlation Visualization Suite

### Goal
Create visualization scripts for correlation analysis (matrices, networks, spectral cleaning, percolation)

### Implementation

**Create:** `src/lrg_eegfc/visuals/` module directory

**Create:** `src/lrg_eegfc/visuals/correlation.py`

**Key functions with notebook patterns:**

#### 1. plot_correlation_heatmap()
```python
def plot_correlation_heatmap(
    corr_matrix: np.ndarray,
    output_path: Path,
    title: str = "Correlation Matrix",
    vmin: float = -1.0,
    vmax: float = 1.0,
    channel_labels: Optional[List[str]] = None,
    figsize: tuple = (8, 8),
) -> Path:
    """Plot correlation matrix heatmap with colorbar.

    Pattern from: pat02_corrnet_bands.ipynb
    - Uses viridis colormap
    - Adds colorbar with divider
    - No axis labels if >50 channels
    """
```

#### 2. plot_correlation_and_network()
```python
def plot_correlation_and_network(
    patient: str,
    phase: str,
    band: str,
    cache_root: Path = Path("data/corr_cache"),
    output_path: Optional[Path] = None,
    cleaned: bool = False,  # Use cleaned or original matrix
    figsize: tuple = (16, 8),
) -> Path:
    """Side-by-side: correlation matrix + network graph with edge weights.

    Pattern from: pat02_corrnet_bands.ipynb (cell f9bd5d3c)

    Left panel: Heatmap of correlation matrix
    Right panel: Network graph with:
        - Edge widths proportional to correlation strength
        - widths = [G[u][v]['weight'] for u, v in G.edges()]
        - Scaled to visual range [0.05, 0.35]
        - Node labels from channel mapping
    """
```

#### 3. plot_marchenko_pastur_comparison()
```python
def plot_marchenko_pastur_comparison(
    patient: str,
    phase: str,
    band: str,
    cache_root: Path = Path("data/corr_cache"),
    output_path: Optional[Path] = None,
    figsize: tuple = (18, 6),
) -> Path:
    """Plot eigenvalue histogram vs Marchenko-Pastur distribution.

    Pattern from: pat02_corrnet_bands.ipynb (cell 49e51f2b)

    3-panel layout:
    Panel 0: Original correlation matrix heatmap
    Panel 1: Cleaned correlation matrix heatmap
    Panel 2: Eigenvalue comparison:
        - Empirical eigenvalue histogram (original)
        - Theoretical MP curve: marchenko_pastur(bins, gamma)
        - Cleaned eigenvalue histogram
        - Vertical lines at lambda_min, lambda_max (red/blue dashed)
        - Shaded gray region between lambda_min and lambda_max
        - Log scale on y-axis, symlog on x-axis

    Key calculations:
    - gamma = time_steps / n_channels
    - lambda_min = (1 - sqrt(1/gamma))^2
    - lambda_max = (1 + sqrt(1/gamma))^2
    - MP_dist = marchenko_pastur(bins, gamma)
    """
```

#### 4. plot_percolation_curves()
```python
def plot_percolation_curves(
    patient: str,
    phase: str,
    band: str,
    cache_root: Path = Path("data/corr_cache"),
    output_path: Optional[Path] = None,
    cleaned: bool = False,
    figsize: tuple = (10, 7),
) -> Path:
    """Plot P_inf and E_inf vs threshold θ.

    Pattern from: poster01.ipynb (cell ad38bec7)

    Shows:
    - P_inf (fraction of nodes in giant component) vs θ
        - Plotted with 'kH' markers (black hexagons), ms=15, hollow
    - E_inf (fraction of edges in giant component) vs θ
        - Plotted with blue line, lw=3
    - Vertical red dashed lines at percolation jumps (node detachments)
    - Annotations showing which nodes detached at each jump
        - Rotated 90°, bbox with rounded corners
        - Arrow pointing to jump location

    Algorithm:
    1. Load correlation matrix (cleaned or original)
    2. Build network G = nx.from_numpy_array(C)
    3. Compute Th, Einf, Pinf = compute_threshold_stats(G)
    4. Find jumps: np.where(np.diff(Pinf) != 0)[0]
    5. For each jump, identify detached nodes:
        - Filter graph at Th[jump]
        - Find giant component
        - Nodes not in giant = detached nodes
    6. Annotate with channel labels
    """
```

**Create:** `src/visualize_correlation.py`
- CLI script for all correlation visualizations
- Arguments:
  - `--patient`, `--phase`, `--band`
  - `--plot-type` (matrix, network, mp, percolation, all)
  - `--cleaned` (use cleaned matrices)
  - `--output-dir` (default: data/figures/correlation/)
- Batch mode: `--batch` generates all plots for patient

---

## PHASE 3: MSC Visualization Suite

### Goal
Create visualization scripts for MSC analysis (matrices and networks only, no cleaning)

### Implementation

**Create:** `src/lrg_eegfc/visuals/msc.py`

**Key functions:**

#### 1. plot_msc_heatmap()
```python
def plot_msc_heatmap(
    msc_matrix: np.ndarray,
    output_path: Path,
    title: str = "MSC Matrix",
    vmin: float = 0.0,  # MSC range is [0, 1] unlike correlation [-1, 1]
    vmax: float = 1.0,
    channel_labels: Optional[List[str]] = None,
    figsize: tuple = (8, 8),
) -> Path:
    """Plot MSC matrix heatmap (note: vmin=0, vmax=1 for MSC unlike correlation)."""
```

#### 2. plot_msc_and_network()
```python
def plot_msc_and_network(
    patient: str,
    phase: str,
    band: str,
    cache_root: Path = Path("data/msc_cache"),
    output_path: Optional[Path] = None,
    sparsify: str = "none",  # Which MSC version to visualize
    n_surrogates: int = 0,
    nperseg: int = 256,
    figsize: tuple = (16, 8),
) -> Path:
    """Side-by-side: MSC matrix + network graph with edge weights.

    Similar to correlation version but for MSC.
    Supports both:
    - Dense (sparsify='none'): All MSC connections
    - Validated (sparsify='soft', n_surrogates=200): Significant connections only

    Pattern follows correlation visualization but NO spectral cleaning
    """
```

**Create:** `src/visualize_msc.py`
- CLI script for MSC visualizations
- Arguments:
  - `--patient`, `--phase`, `--band`
  - `--sparsify` (none, soft)
  - `--n-surrogates` (if sparsify='soft')
  - `--plot-type` (matrix, network, all)
  - `--output-dir` (default: data/figures/msc/)

---

## PHASE 4: LRG Visualization Suite

### Goal
Create comprehensive LRG visualization (entropy, dendrograms, ultrametrics)

### Implementation

**Create:** `src/lrg_eegfc/visuals/lrg.py`

**Key functions:**

#### 1. plot_lrg_entropy_curves()
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

    Pattern from: plotting.py plot_entropy()

    Uses cached LRG result fields:
    - entropy_tau: τ values (log-spaced)
    - entropy_1_minus_S: Normalized entropy (1-S)
    - entropy_C: Spectral complexity (C)

    Visualization:
    - Log-scale x-axis (tau)
    - Normalized y-axis [0, 1]
    - Two curves: 1-S (blue) and C (red)
    - Legend with labels
    """
```

#### 2. plot_lrg_dendrogram()
```python
def plot_lrg_dendrogram(
    patient: str,
    phase: str,
    band: str,
    fc_method: str,
    cache_root: Path = Path("data/lrg_cache"),
    output_path: Optional[Path] = None,
    orientation: str = "top",  # "top", "right", "bottom", "left"
    figsize: tuple = (12, 8),
) -> Path:
    """Plot hierarchical dendrogram with optimal threshold line.

    Pattern from: NEW_distance_of_distances.ipynb

    Uses cached fields:
    - linkage_matrix: Hierarchical clustering linkage
    - optimal_threshold: Cutting threshold

    Features:
    - Optimal leaf ordering: scipy.cluster.hierarchy.optimal_leaf_ordering(linkage, condensed_dists)
    - Color threshold: above_threshold_color='k', color below threshold varies
    - Log-scale y-axis
    - Threshold line: ax.axhline(optimal_threshold, color='b', linestyle='--')
    - No labels if n_nodes > 50 (too crowded)
    """
```

#### 3. plot_ultrametric_heatmap()
```python
def plot_ultrametric_heatmap(
    patient: str,
    phase: str,
    band: str,
    fc_method: str,
    cache_root: Path = Path("data/lrg_cache"),
    output_path: Optional[Path] = None,
    figsize: tuple = (10, 10),
) -> Path:
    """Plot ultrametric distance matrix as heatmap.

    Uses cached ultrametric_matrix (condensed form)

    Algorithm:
    1. Convert condensed to square: scipy.spatial.distance.squareform()
    2. Plot symmetric heatmap
    3. Colorbar with distance interpretation
    4. Channel labels if available
    """
```

#### 4. plot_lrg_full_panel()
```python
def plot_lrg_full_panel(
    patient: str,
    phase: str,
    band: str,
    fc_method: str,
    cache_root: Path = Path("data/lrg_cache"),
    output_path: Optional[Path] = None,
    figsize: tuple = (20, 12),
) -> Path:
    """4-panel summary: Entropy | Dendrogram | Ultrametric heatmap | Network

    Layout:
    ┌─────────────┬─────────────┐
    │  Entropy    │  Dendrogram │
    ├─────────────┼─────────────┤
    │  Ultrametric│  Network    │
    └─────────────┴─────────────┘

    Comprehensive single-figure visualization
    Network colored by dendrogram clusters
    """
```

**Create:** `src/visualize_lrg.py`
- CLI script for LRG visualizations
- Arguments:
  - `--patient`, `--phase`, `--band`, `--fc-method`
  - `--plot-type` (entropy, dendrogram, ultrametric, full, all)
  - `--output-dir` (default: data/figures/lrg/)

---

## PHASE 5: Comparison Visualization Suite

### Goal
Create visualizations for comparing FC methods and phases (reorganization analysis)

### Implementation

**Create:** `src/lrg_eegfc/visuals/comparisons.py`

**Key functions:**

#### 1. plot_msc_vs_corr_comparison()
```python
def plot_msc_vs_corr_comparison(
    patient: str,
    phase: str,
    band: str,
    cache_root: Path = Path("data/lrg_cache"),
    output_path: Optional[Path] = None,
    figsize: tuple = (18, 12),
) -> Path:
    """Side-by-side comparison of MSC and Correlation LRG results.

    2x3 grid layout:
    ┌─────────────┬─────────────┬─────────────┐
    │ MSC Entropy │ MSC Dendro  │ MSC Ultra   │
    ├─────────────┼─────────────┼─────────────┤
    │ Corr Entropy│ Corr Dendro │ Corr Ultra  │
    └─────────────┴─────────────┴─────────────┘

    Allows visual comparison of network structure between methods
    """
```

#### 2. plot_phase_reorganization()
```python
def plot_phase_reorganization(
    patient: str,
    band: str,
    fc_method: str,
    phases: List[str] = ["rsPre", "taskLearn", "taskTest", "rsPost"],
    cache_root: Path = Path("data/lrg_cache"),
    output_path: Optional[Path] = None,
    figsize: tuple = (20, 12),
) -> Path:
    """Visualize network reorganization across phases.

    Pattern from: NEW_distance_of_distances.ipynb

    Shows:
    - Top row: Dendrograms for all phases side-by-side
        - Same y-axis scale for comparison
        - Color-coded by phase
    - Middle row: Ultrametric distance heatmaps
    - Bottom row: Cluster membership correlation matrix across phases
        - Shows how cluster assignments change between phases

    Metrics displayed:
    - Tree topology similarity scores (from compare.py)
    - Matrix distance evolution
    - Cluster membership correlation heatmap
    """
```

#### 3. plot_reorganization_metrics()
```python
def plot_reorganization_metrics(
    patient: str,
    band: str,
    fc_method: str,
    comparison_csv: Path,  # From compare_fc_methods.py output
    output_path: Optional[Path] = None,
    figsize: tuple = (16, 10),
) -> Path:
    """Plot reorganization metrics from comparison analysis.

    From user requirement: "correlation between cluster belonging of each node
    or measure of distance of ultrametric matrices or even topological measure
    of similarity hierarchical tree to hierarchical tree"

    4-panel layout:
    ┌─────────────┬─────────────┐
    │ Matrix Dist │ Tree Sim    │
    ├─────────────┼─────────────┤
    │ Cluster Corr│ Multiscale  │
    └─────────────┴─────────────┘

    Reads from comparison CSV columns:
    - matrix_distance, multiscale_distance
    - tree_similarity, rank_correlation
    - quantile_rmse, scale_profile_distance

    Shows evolution across phase pairs
    """
```

#### 4. compute_cluster_membership_correlation()
```python
def compute_cluster_membership_correlation(
    linkage_1: np.ndarray,
    linkage_2: np.ndarray,
    threshold: Optional[float] = None,
) -> float:
    """Compute correlation of cluster memberships between two hierarchies.

    Algorithm:
    1. Cut dendrograms at threshold (or optimal threshold from linkage)
        - scipy.cluster.hierarchy.fcluster(linkage, threshold, criterion='distance')
    2. Get cluster labels for each node
    3. Compute Spearman correlation of label vectors
        - scipy.stats.spearmanr(labels_1, labels_2)

    Returns: correlation coefficient [-1, 1]
    - 1.0: Identical clustering
    - 0.0: No relationship
    - -1.0: Inverse clustering
    """
```

**Create:** `src/visualize_comparisons.py`
- CLI script for comparison visualizations
- Arguments:
  - `--patient`, `--band`
  - `--comparison-type` (msc-vs-corr, phase-reorganization, metrics)
  - `--fc-method` (for phase-reorganization and metrics)
  - `--comparison-csv` (for metrics mode)
  - `--output-dir` (default: data/figures/comparisons/)

---

## PHASE 6: Documentation and Integration

### Goal
Create comprehensive VISUALS.md guide and integrate with existing pipeline

### Implementation

**Create:** `VISUALS.md`

**Structure:**
```markdown
# Visualization Pipeline Guide

## Overview
- Visualization philosophy: separate from computation
- Reading from cached data
- Output organization
- Philosophy: "not having to compute data"

## Cleaning Pipeline
- Marchenko-Pastur spectral filtering
- Percolation thresholding (just before first node detachment)
- When to use cleaned vs uncleaned matrices
- Cache naming conventions

## Correlation Visualizations
- Matrix heatmaps
- Network graphs with weighted edges
- Marchenko-Pastur comparison (3-panel)
- Percolation curves (P_inf, E_inf vs θ)
- Command examples for each plot type

## MSC Visualizations
- Matrix heatmaps (MSC range 0-1)
- Network graphs (dense and validated)
- Why MSC doesn't need cleaning
- Command examples

## LRG Visualizations
- Entropy curves (1-S and C vs tau)
- Dendrograms with optimal threshold
- Ultrametric heatmaps
- Full panel summaries
- Command examples

## Comparison Visualizations
- MSC vs Correlation side-by-side
- Phase reorganization analysis
- Memory effects visualization
- Reorganization metrics
- Command examples

## Complete Workflow Example
1. Clean correlation matrices
2. Generate all visualizations
3. Compare MSC vs Correlation
4. Analyze phase reorganization
5. Interpret results

## Output Organization
data/figures/
├── correlation/
│   ├── Pat_02/
│   │   ├── beta_rsPre_matrix.png
│   │   ├── beta_rsPre_network.png
│   │   ├── beta_rsPre_mp_comparison.png
│   │   └── beta_rsPre_percolation.png
├── msc/
│   ├── Pat_02/
│   │   ├── beta_rsPre_matrix.png
│   │   └── beta_rsPre_network.png
├── lrg/
│   ├── Pat_02/
│   │   ├── beta_rsPre_corr_entropy.png
│   │   ├── beta_rsPre_corr_dendrogram.png
│   │   ├── beta_rsPre_corr_ultrametric.png
│   │   └── beta_rsPre_corr_full_panel.png
└── comparisons/
    ├── Pat_02/
    │   ├── beta_rsPre_msc_vs_corr.png
    │   ├── beta_corr_phase_reorganization.png
    │   └── beta_corr_reorganization_metrics.png

## Key Visualization Patterns
- Edge width scaling: widths = [G[u][v]['weight'] for u, v in G.edges()]
- Percolation jump detection: jumps = np.where(np.diff(Pinf) != 0)[0]
- Marchenko-Pastur bounds: lambda_min/max = (1 ± sqrt(1/Q))^2
- Optimal leaf ordering for dendrograms
- Cluster coloring from dendrogram
```

**Modify:** `ANALYSIS_COMMANDS.md` ✅
- Added Step 2C: Clean Correlation Matrices
- Integration with cleaning pipeline

**Create:** `src/lrg_eegfc/visuals/__init__.py`
- Export all visualization functions
- Organized imports by category (correlation, msc, lrg, comparisons)

---

## PHASE 7: Testing and Validation

### Implementation

**Create:** `tests/test_visuals.py`
- Test cleaning workflow with mock data
- Test visualization functions produce valid outputs
- Test file creation and naming conventions
- Validate plot dimensions and content

**Create:** `tests/test_cleaning.py`
- Test Marchenko-Pastur bounds calculation
- Test percolation threshold detection
- Test cleaned matrix properties (symmetry, range)
- Test caching and loading

**Create:** Example notebook: `ipynb/EXAMPLE_visualization_workflow.ipynb`
- Section 1: Cleaning correlation matrices
- Section 2: Correlation visualizations (all types)
- Section 3: MSC visualizations
- Section 4: LRG visualizations
- Section 5: Comparison visualizations
- Section 6: Interpretation guide
- Uses Pat_02 beta rsPre as running example

---

## Implementation Order

### Week 1: Cleaning and Core Visualizations
1. ✅ **Days 1-2**: Phase 1 - Cleaning workflow (workflow_cleaning.py, clean_correlation_matrices.py)
2. **Days 3-4**: Phase 2 - Correlation visualizations (visuals/correlation.py, visualize_correlation.py)
3. **Day 5**: Testing cleaning and correlation visuals

### Week 2: MSC and LRG Visualizations
4. **Day 1**: Phase 3 - MSC visualizations (visuals/msc.py, visualize_msc.py)
5. **Days 2-3**: Phase 4 - LRG visualizations (visuals/lrg.py, visualize_lrg.py)
6. **Days 4-5**: Testing MSC and LRG visuals

### Week 3: Comparisons and Documentation
7. **Days 1-2**: Phase 5 - Comparison visualizations (visuals/comparisons.py, visualize_comparisons.py)
8. **Days 3-4**: Phase 6 - Documentation (VISUALS.md, integration)
9. **Day 5**: Phase 7 - Final testing and example notebook

---

## Files to Create (15 new files)

### Core Workflow
1. ✅ `src/lrg_eegfc/workflow_cleaning.py` - Correlation cleaning workflow
2. ✅ `src/clean_correlation_matrices.py` - CLI for cleaning pipeline

### Visualization Modules
3. `src/lrg_eegfc/visuals/__init__.py` - Module exports
4. `src/lrg_eegfc/visuals/correlation.py` - Correlation visualizations
5. `src/lrg_eegfc/visuals/msc.py` - MSC visualizations
6. `src/lrg_eegfc/visuals/lrg.py` - LRG visualizations
7. `src/lrg_eegfc/visuals/comparisons.py` - Comparison visualizations

### CLI Scripts
8. `src/visualize_correlation.py` - Correlation visualization CLI
9. `src/visualize_msc.py` - MSC visualization CLI
10. `src/visualize_lrg.py` - LRG visualization CLI
11. `src/visualize_comparisons.py` - Comparison visualization CLI

### Documentation
12. `VISUALS.md` - Comprehensive visualization guide

### Testing
13. `tests/test_visuals.py` - Visualization test suite
14. `tests/test_cleaning.py` - Cleaning workflow tests

### Examples
15. `ipynb/EXAMPLE_visualization_workflow.ipynb` - Example notebook

---

## Files to Modify (3 existing files)

1. ✅ `src/lrg_eegfc/__init__.py` - Add cleaning and visualization exports
2. ✅ `ANALYSIS_COMMANDS.md` - Add visualization commands
3. `README.md` - Mention visualization pipeline

---

## Key Design Decisions

### Why separate visualization scripts?
- **User requirement**: "not having to compute data"
- Faster iteration on plot styling
- Can regenerate plots without recomputing expensive analysis
- Cleaner separation of concerns

### Why clean correlation but not MSC?
- **Correlation**: Affected by temporal correlations and finite-size effects → needs spectral cleaning
- **MSC**: Already frequency-specific, spectral estimation is robust to noise
- **User requirement**: "for the msc there is no need of thresholding or removing spectral noise"

### Why suffix instead of separate directory?
- **User preference**: "Same cache with different filename suffix"
- Easier to find corresponding cleaned/uncleaned versions
- Less directory nesting
- Clear naming: `beta_rsPre_corr.npy` vs `beta_rsPre_corr_cleaned.npy`

### Percolation threshold strategy
- **Goal**: "just before the first node gets detached"
- **Method**: Use `find_threshold_jumps()` to detect jumps in P_inf
- **Choose**: Threshold at `jumps[0] - 1` (just before first detachment)
  - `jumps[0]` is WHERE first detachment occurs
  - `jumps[0] - 1` is just before it
- **Validate**: Check that `nx.number_connected_components == 1` at chosen threshold

### Visualization output organization
```
data/figures/
├── correlation/
│   └── {patient}/
│       └── {band}_{phase}_{plot_type}.png
├── msc/
│   └── {patient}/
│       └── {band}_{phase}_{plot_type}.png
├── lrg/
│   └── {patient}/
│       └── {band}_{phase}_{fc_method}_{plot_type}.png
└── comparisons/
    └── {patient}/
        └── {comparison_type}_{band}.png
```

---

## User-Requested Additional Visualizations

From user: "if you think to other fancy plots that we can do for MSC or kind of analysis feel free to suggests, but add them at the end of the whole implementation"

### Potential Additional Plots (Post-Implementation)
1. **Coherence spectra evolution**
   - MSC(f) across frequencies for key channel pairs
   - How coherence changes across phases

2. **Network modularity analysis**
   - Community detection on correlation/MSC networks
   - Modularity evolution across phases

3. **Topological features**
   - Degree distribution
   - Clustering coefficient
   - Small-world metrics

4. **Multiscale analysis**
   - MSC at different nperseg values
   - Sensitivity to window length

5. **Cross-frequency coupling**
   - Phase-amplitude coupling visualization
   - Cross-frequency coherence matrices

**Implementation strategy for additional plots:**
- Add as separate module: `src/lrg_eegfc/visuals/advanced.py`
- Only after core visualizations are complete and tested
- Document in separate section of VISUALS.md

---

## Expected Outcomes

After complete implementation:

✅ **Cleaning Pipeline**: Marchenko-Pastur + percolation thresholding for correlation
✅ **Correlation Visuals**: Matrix, network, MP comparison, percolation curves
✅ **MSC Visuals**: Matrix, network (dense and validated)
✅ **LRG Visuals**: Entropy, dendrograms, ultrametric heatmaps
✅ **Comparison Visuals**: MSC vs Corr, phase reorganization, memory effects
✅ **CLI Tools**: 4 visualization scripts for different analysis types
✅ **Documentation**: VISUALS.md comprehensive guide
✅ **Testing**: Unit tests for cleaning and visualization
✅ **Example**: Jupyter notebook demonstrating complete workflow

---

## Critical Success Factors

1. **Reuse existing code** - Leverage plotting.py, thresholding.py, clustering.py patterns
2. **Follow notebook patterns** - Use exact plotting code from working notebooks
3. **Cache strategy** - Cleaned matrices with `_cleaned.npy` suffix, metadata in separate .npz
4. **Separation of concerns** - Computation scripts vs visualization scripts
5. **Comprehensive documentation** - VISUALS.md should be self-contained guide
6. **User validation** - Ensure reorganization analysis matches user's research questions

---

## Phase Comparison Metrics (From User Requirement)

User wants: "correlation between cluster belonging of each node or measure of distance of ultrametric matrices or even topological measure of similarity hierarchical tree to hierarchical tree. all of these have been explored in the notebooks"

### Metrics to Implement

1. **Cluster Membership Correlation**
   - Cut dendrograms at optimal threshold
   - Get cluster labels for each node
   - Compute Spearman correlation of label vectors across phases

2. **Ultrametric Distance Measures** (Already in compare.py)
   - Matrix distance (Frobenius)
   - Multiscale distance
   - Quantile RMSE
   - Rank correlation
   - Scale profile distance

3. **Tree Topology Similarity**
   - Tree similarity from `compare_ultrametric_trees()` (already implemented)
   - Cophenetic correlation
   - Robinson-Foulds distance (if available in lrgsglib)

4. **Visualization of Reorganization Strength**
   - Classify as: no reorganization, weak, strong based on metric thresholds
   - Color-code phases in plots
   - Summary statistics across patients

All these are already computed in the comparison pipeline - visualization just needs to read from CSV and plot effectively.

---

## Implementation Progress

### ✅ Completed
- **PHASE 1: Correlation Cleaning Workflow**
  - ✅ `src/lrg_eegfc/workflow_cleaning.py` - CleanedCorrResult dataclass, cleaning functions
  - ✅ `src/clean_correlation_matrices.py` - CLI script for batch cleaning
  - ✅ `src/lrg_eegfc/__init__.py` - Exports added
  - ✅ `ANALYSIS_COMMANDS.md` - Step 2C added

### 🚧 In Progress
- None

### ⏳ Pending
- PHASE 2: Correlation Visualization Suite
- PHASE 3: MSC Visualization Suite
- PHASE 4: LRG Visualization Suite
- PHASE 5: Comparison Visualization Suite
- PHASE 6: Documentation and Integration
- PHASE 7: Testing and Validation
