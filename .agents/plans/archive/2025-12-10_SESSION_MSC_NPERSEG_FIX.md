# Session: MSC nperseg Parameter Fix and Visualization Plan

**Date:** 2025-12-10
**Context:** Fixing MSC network visualization and preparing MSC summary plots

---

## Issue Discovered

### Problem: MSC Networks Were Fully Connected

When checking MSC visualizations, discovered that MSC networks were showing as complete hairballs (fully connected graphs):

```python
# Checking Pat_02 beta rsPre with nperseg=256
m = np.load('data/msc_cache/Pat_02/beta_rsPre_msc_sparsify-none_nperseg-256.npy')
G = nx.from_numpy_array(m)

# Results:
Number of nodes: 117
Number of edges: 6786  # <-- FULLY CONNECTED (117*116/2 = 6786)
Density: 1.0
Edge weight min: 2.191637018401794e-05
Edge weight max: 0.9725248022062299
```

**Root cause:**
- `nx.from_numpy_array()` creates an edge for every non-zero value
- MSC values are never exactly zero (smallest ~2e-05)
- With `nperseg=256`, frequency resolution is too coarse → all channel pairs have non-zero MSC
- Result: Complete graph with 6786 edges (useless for visualization)

### Verification Against Notebook

Checked `ipynb/UTILS-FC_COMPARISON_PIPELINE.ipynb` and found they use:

```python
msc_matrices = coherence_fc_pipeline(
    X,
    fs,
    bands=BRAIN_BANDS,
    sparsify="none",
    nperseg=512,  # <-- Different from our cache (256)!
)
```

And their results show naturally sparse matrices:
```
delta       : 117 nodes, ~   0 edges (threshold=0.01)
theta       : 117 nodes, ~4163 edges (threshold=0.01)
alpha       : 117 nodes, ~4163 edges (threshold=0.01)
beta        : 117 nodes, ~4660 edges (threshold=0.01)
low_gamma   : 117 nodes, ~3894 edges (threshold=0.01)
high_gamma  : 117 nodes, ~5106 edges (threshold=0.01)
```

**Conclusion:** Higher `nperseg` → Better frequency resolution → Naturally sparser matrices

---

## Analysis: nperseg Parameter Constraints

### Frequency Resolution Requirements

Key constraint: **Delta band needs good frequency resolution**

```
BRAIN_BANDS frequency ranges:
delta       :   0.53-  4.00 Hz (BW=  3.47 Hz) → need df < 1.16 Hz  ← BOTTLENECK
theta       :   4.00-  8.00 Hz (BW=  4.00 Hz) → need df < 1.33 Hz
alpha       :   8.00- 13.00 Hz (BW=  5.00 Hz) → need df < 1.67 Hz
beta        :  13.00- 30.00 Hz (BW= 17.00 Hz) → need df < 5.67 Hz
low_gamma   :  30.00- 80.00 Hz (BW= 50.00 Hz) → need df < 16.67 Hz
high_gamma  :  80.00-300.00 Hz (BW=220.00 Hz) → need df < 73.33 Hz
```

### Sampling Rates in Data

- **fs=2048 Hz** (most common)
- **fs=1024 Hz** (some recordings)

### nperseg Performance Analysis

For **fs=2048 Hz**:

| nperseg | Window(s) | Freq Res (Hz) | n_segments | Status |
|---------|-----------|---------------|------------|---------|
| 256     | 0.125     | 8.00          | ~10336     | ❌ Too coarse for delta (df=8 Hz) |
| 512     | 0.250     | 4.00          | ~5167      | ⚠️ Marginal for delta (df=4 Hz) |
| **1024** | **0.500** | **2.00**      | **~2583**  | **✓ Good balance** |
| 2048    | 1.000     | 1.00          | ~1291      | ✓ Best for delta, but fewer segments |
| 4096    | 2.000     | 0.50          | ~645       | ✓ Excellent resolution, but few segments |

For **fs=1024 Hz**:

| nperseg | Window(s) | Freq Res (Hz) | n_segments | Status |
|---------|-----------|---------------|------------|---------|
| 512     | 0.500     | 2.00          | ~2583      | ✓ Good |
| **1024** | **1.000** | **1.00**      | **~1291**  | **✓ Excellent** |

### Recommendation: nperseg=1024

**Rationale:**
1. ✓ Adequate frequency resolution for delta band (df=2 Hz at fs=2048)
2. ✓ Reasonable window length (0.5s - good for EEG stationarity)
3. ✓ Plenty of segments for averaging (~2580 at fs=2048)
4. ✓ Naturally produces sparser matrices (better for visualization)
5. ✓ Works well for both fs=1024 and fs=2048

**Limits:**
- **Lower bound:** nperseg ≥ 512 (avoid too coarse resolution)
- **Upper bound:** nperseg ≤ 4096 (avoid too few segments)
- **Sweet spot:** 1024-2048

---

## Changes Made

### File 1: `src/compute_msc_matrices.py`

**Line 56:**
```python
# BEFORE
default=256,
help="Window length for Welch's method (default: 256)"

# AFTER
default=1024,
help="Window length for Welch's method (default: 1024)"
```

### File 2: `src/visualize_msc.py`

**Line 152-155:**
```python
# BEFORE
parser.add_argument(
    "--nperseg",
    type=int,
    default=256,
    help="Window length for Welch's method. Default: 256",
)

# AFTER
parser.add_argument(
    "--nperseg",
    type=int,
    default=1024,
    help="Window length for Welch's method. Default: 1024",
)
```

### File 3: `src/lrg_eegfc/workflow_msc.py`

Updated **3 function signatures**:

**1. `get_msc_cache_path()` - Line 72:**
```python
# BEFORE
nperseg: int = 256,

# AFTER
nperseg: int = 1024,
```

**2. `load_msc_matrix()` - Line 117:**
```python
# BEFORE
nperseg: int = 256,

# AFTER
nperseg: int = 1024,
```

**3. `compute_msc_matrix()` - Line 162:**
```python
# BEFORE
nperseg: int = 256,

# AFTER
nperseg: int = 1024,
```

**Updated docstring - Line 197-198:**
```python
# BEFORE
nperseg : int, optional
    Window length for Welch's method (default: 256)

# AFTER
nperseg : int, optional
    Window length for Welch's method (default: 1024)
```

### File 4: `src/lrg_eegfc/visuals/msc.py`

Updated **2 function signatures**:

**1. `plot_msc_and_network()` - Line 109:**
```python
# BEFORE
nperseg: int = 256,

# AFTER
nperseg: int = 1024,
```

**2. `plot_msc_comparison_dense_vs_validated()` - Line 254:**
```python
# BEFORE
nperseg: int = 256,

# AFTER
nperseg: int = 1024,
```

---

## Impact

### Before (nperseg=256):
- Frequency resolution: df = 4.0 Hz (at fs=2048)
- MSC matrices: Fully connected (6786 edges)
- Network visualization: Unusable hairball
- Delta band: Inadequate resolution (needs < 1.16 Hz)

### After (nperseg=1024):
- Frequency resolution: df = 2.0 Hz (at fs=2048)
- MSC matrices: Naturally sparse (~4000-5000 edges)
- Network visualization: Clear structure visible
- Delta band: Adequate resolution

---

## Next Steps

### 1. Recompute MSC Matrices

**REQUIRED:** Recompute with new nperseg=1024 default:

```bash
# Single patient
python src/compute_msc_matrices.py --patient Pat_02 --overwrite --verbose

# All patients
python src/compute_msc_matrices.py --overwrite --verbose
```

**Cache files will be:**
```
data/msc_cache/Pat_XX/{band}_{phase}_msc_sparsify-none_nperseg-1024.npy
```

### 2. Create MSC Summary Visualization

**User request:** Create MSC summary plot similar to correlation summary, with 3-panel horizontal layout:

```
┌─────────────┬─────────────┬─────────────┐
│ MSC Matrix  │ Full Network│ Percolation │
│ (already    │ (unthresh.) │   Curves    │
│  positive)  │             │ (info only) │
└─────────────┴─────────────┴─────────────┘
```

**Key differences from correlation:**
- **Panel 1:** MSC matrix (already [0,1], no abs() needed)
- **Panel 2:** Network from FULL unthresholded MSC (but with nperseg=1024, will be naturally sparse)
- **Panel 3:** Percolation curves showing how modules detach (informational - threshold NOT applied to network)

**Implementation needed:**
- Add `plot_msc_summary()` function to `src/lrg_eegfc/visuals/msc.py`
- Add `--plot-type summary` option to `src/visualize_msc.py`
- Update `scripts/run_full_analysis.sh` to generate MSC summaries

### 3. Update Pipeline Script

Modify `scripts/run_full_analysis.sh` to include MSC summary visualization step after MSC computation.

---

## Technical Notes

### Why nperseg Affects Sparsity

**Welch's Method:**
1. Divide signal into `nperseg`-length windows
2. Compute FFT for each window
3. Average power/cross-power across windows
4. Longer windows → Better frequency resolution

**Frequency Resolution:**
```
df = fs / nperseg
```

**Why Better Resolution → Sparser Matrices:**
- Finer frequency bins → More precise band averaging
- Precise averaging → Weak coherence approaches zero more accurately
- More near-zero values → Sparser adjacency matrix when converted to graph

### Cache Filename Convention

MSC cache files encode all parameters:
```
{band}_{phase}_msc_sparsify-{method}_nperseg-{nperseg}.npy

# Dense (no sparsification)
beta_rsPre_msc_sparsify-none_nperseg-1024.npy

# Validated (with surrogates)
beta_rsPre_msc_sparsify-soft_nsurr-200_nperseg-1024.npy
```

**Important:** Old cache files with `nperseg-256` will remain but won't be used by default. Can be deleted after verification.

---

## Comparison: Correlation vs MSC

### Workflow Differences

| Aspect | Correlation | MSC |
|--------|-------------|-----|
| **Input** | Bandpass filtered time series | Raw time series |
| **Method** | Pearson correlation | Welch's coherence |
| **Range** | [-1, 1] | [0, 1] |
| **Cleaning** | Marchenko-Pastur spectral cleaning | Not needed (already freq-specific) |
| **Thresholding** | Percolation-based (applied) | Percolation (shown but not applied) |
| **Sparsity** | From thresholding | Natural (from nperseg) |

### Visualization Differences

| Plot | Correlation | MSC |
|------|-------------|-----|
| **Matrix** | Abs value, viridis | Direct, viridis |
| **Network** | Thresholded | Full (but naturally sparse with nperseg=1024) |
| **Cleaning** | MP comparison | Not applicable |
| **Percolation** | Shows threshold selection | Shows for info only |

---

## File Status Summary

### ✅ Updated Files (nperseg: 256 → 1024)
- `src/compute_msc_matrices.py`
- `src/visualize_msc.py`
- `src/lrg_eegfc/workflow_msc.py`
- `src/lrg_eegfc/visuals/msc.py`

### 🚧 Pending Implementation
- `plot_msc_summary()` function (3-panel layout)
- MSC summary CLI integration
- Pipeline script update

### ⚠️ Requires Recomputation
- All MSC cache files (currently have nperseg=256)

---

## User Context

**User has:**
- Fixed correlation visualization code with Codex
- Generated correlation summary visualizations successfully
- Modified `scripts/run_full_analysis.sh` to use `--plot-type summary`
- Now wants equivalent MSC summary visualization

**User wants:**
- MSC summary with 3 panels (Matrix | Network | Percolation)
- No thresholding applied to MSC network (unlike correlation)
- Percolation curves shown for informational/analytical purposes only
- Clean, interpretable network visualization (enabled by nperseg=1024)

**User has confirmed:**
- nperseg=1024 is appropriate default
- Understands frequency resolution constraints
- Ready to proceed with MSC summary implementation

---

## Commands to Run After Recomputation

```bash
# Recompute MSC with nperseg=1024
python src/compute_msc_matrices.py --patient Pat_02 --overwrite --verbose

# Test MSC network visualization (should now be sparse)
python src/visualize_msc.py --patient Pat_02 --phase rsPre --band beta --plot-type network --verbose

# After implementing MSC summary:
python src/visualize_msc.py --patient Pat_02 --phase rsPre --band beta --plot-type summary --verbose
```

---

**Session saved:** 2025-12-10
**Ready to compact conversation and continue with MSC summary visualization implementation**
