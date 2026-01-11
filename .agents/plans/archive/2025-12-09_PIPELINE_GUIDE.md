# Complete Analysis Pipeline Guide

This guide walks through the complete functional connectivity analysis pipeline from raw EEG data to network comparisons.

## Overview

The pipeline consists of 4 main stages:

1. **Data Inspection** - Validate and inspect patient data
2. **FC Matrix Computation** - Compute MSC and correlation matrices
3. **LRG Analysis** - Extract ultrametric distances via renormalization group
4. **Comparison & Aggregation** - Compare networks and analyze memory effects

## Prerequisites

- Patient data in standardized format: `data/stereoeeg_patients/Pat_XX/{phase}.mat`
- Python environment with lrg_eegfc and lrgsglib installed

## Stage 1: Data Inspection

**Check data completeness and divergences:**

```bash
python src/inspect_patient_data.py
```

This generates `PATIENT_DATA_REPORT.txt` with:
- Missing phase files
- Missing sampling rates (fs)
- Missing channel labels
- Data format inconsistencies

**Review reports:**
- `PATIENT_DATA_REPORT.txt` - Detailed technical report
- `DATA_STATUS.md` - Summary for collaborators

**Key findings to check:**
- Sampling rate divergences (Pat_02: 2048 Hz, Pat_03: 1024 Hz)
- Missing fs parameters (Pat_05-08)
- Missing phase files (Pat_06, Pat_07)

## Stage 2: FC Matrix Computation

### 2A: Correlation-Based FC

**Compute correlation matrices for all patients:**

```bash
python src/compute_corr_matrices.py --verbose
```

**Parameters:**
- `--filter-type` - Filter type: "abs" (default), "pos", "neg", "none"
- `--zero-diagonal` - Set diagonal to zero (default: True)
- `--filter-order` - Bandpass filter order (default: 4)
- `--patients` - Specific patients (default: all)

**Output:** `data/corr_cache/{patient}/{band}_{phase}_corr_ftype-{type}_zdiag-{bool}.npy`

**Example - positive correlations only:**
```bash
python src/compute_corr_matrices.py --filter-type pos --patients Pat_02 Pat_03
```

### 2B: MSC-Based FC

**Compute MSC matrices (dense, no validation):**

```bash
python src/compute_msc_matrices.py --verbose
```

**Parameters:**
- `--sparsify` - Sparsification: "none" (default, dense MSC) or "soft"
- `--n-surrogates` - Number of surrogates for validation (if sparsify="soft")
- `--nperseg` - Window length for Welch's method (default: 256)
- `--patients` - Specific patients (default: all)

**Output:** `data/msc_cache/{patient}/{band}_{phase}_msc_sparsify-{method}_nsurr-{n}_nperseg-{nperseg}.npy`

**Example - validated MSC with surrogates:**
```bash
python src/compute_msc_matrices.py --sparsify soft --n-surrogates 200 --patients Pat_02
```

### Sampling Rate Handling

**Important:** The pipeline automatically extracts sampling rates from patient data:
- Pat_02: 2048 Hz (from data)
- Pat_03: 1024 Hz (from data)
- Pat_05-08: 2048 Hz (default, with warning if fs missing)

No manual fs specification needed!

## Stage 3: LRG Analysis

**Compute ultrametric distances via Laplacian Renormalization Group:**

### 3A: LRG from Correlation Matrices

```bash
python src/compute_lrg_analysis.py --fc-method corr --verbose
```

### 3B: LRG from MSC Matrices

```bash
python src/compute_lrg_analysis.py --fc-method msc --verbose
```

**Parameters:**
- `--fc-method` - Required: "msc" or "corr"
- `--entropy-steps` - Steps for entropy computation (default: 400)
- `--entropy-t1` - Tau range start (log scale, default: -3.0)
- `--entropy-t2` - Tau range end (log scale, default: 5.0)
- `--patients` - Specific patients (default: all)

**Output:** `data/lrg_cache/{patient}/{band}_{phase}_lrg_{fc_method}.npz`

**NPZ file contents:**
- `ultrametric_matrix` - Condensed ultrametric distance matrix
- `linkage_matrix` - Hierarchical clustering linkage
- `entropy_tau` - Tau values for entropy curve
- `entropy_1_minus_S` - 1-S (normalized entropy)
- `entropy_C` - C (spectral complexity)
- `optimal_threshold` - Dendrogram cutting threshold
- Metadata: patient, phase, band, fc_method, n_nodes

## Stage 4: Comparison & Aggregation

### 4A: Compare MSC vs Correlation

**Compare FC methods across all bands/phases:**

```bash
python src/compare_fc_methods.py --mode fc-methods --output results/msc_vs_corr.csv
```

**Output metrics:**
- `matrix_distance` - Frobenius distance between ultrametric matrices
- `multiscale_distance` - Multi-scale hierarchical comparison
- `quantile_rmse` - RMSE of quantile distributions
- `rank_correlation` - Spearman rank correlation
- `scale_profile_distance` - Distance between scale profiles
- `scaled_distance` - Normalized distance metric
- `tree_similarity` - Tree topology similarity (0-1)

### 4B: Memory Effects Analysis

**Compare pre vs post phases to detect memory effects:**

```bash
python src/compare_fc_methods.py \\
    --mode memory \\
    --fc-method corr \\
    --output results/memory_effects.csv
```

This compares:
- `rsPre` vs `rsPost` (resting state before/after)
- `taskLearn` vs `taskTest` (learning vs testing)

**Per patient/band analysis:**

```bash
python src/compare_fc_methods.py \\
    --mode memory \\
    --fc-method corr \\
    --patients Pat_02 \\
    --bands beta alpha
```

## Complete Pipeline Example

**Full workflow for Pat_02 and Pat_03:**

```bash
# 1. Inspect data
python src/inspect_patient_data.py

# 2. Compute FC matrices
python src/compute_corr_matrices.py --patients Pat_02 Pat_03 --verbose
python src/compute_msc_matrices.py --patients Pat_02 Pat_03 --verbose

# 3. Compute LRG analysis
python src/compute_lrg_analysis.py --fc-method corr --patients Pat_02 Pat_03 --verbose
python src/compute_lrg_analysis.py --fc-method msc --patients Pat_02 Pat_03 --verbose

# 4. Compare networks
python src/compare_fc_methods.py --patients Pat_02 Pat_03 --output results/comparison.csv

# 5. Memory effects
python src/compare_fc_methods.py --mode memory --fc-method corr --output results/memory.csv
```

## Python API Usage

### Correlation FC

```python
from lrg_eegfc import compute_corr_matrix

# Compute correlation matrix
result = compute_corr_matrix(
    "Pat_02", "rsPre", "beta",
    filter_type="abs",
    zero_diagonal=True
)

print(f"Shape: {result.adjacency_matrix.shape}")
print(f"Mean correlation: {result.mean_corr:.4f}")
```

### MSC FC

```python
from lrg_eegfc import compute_msc_matrix

# Dense MSC (no validation)
result = compute_msc_matrix(
    "Pat_02", "rsPre", "beta",
    sparsify="none"
)

# Validated MSC with surrogates
result = compute_msc_matrix(
    "Pat_02", "rsPre", "beta",
    sparsify="soft",
    n_surrogates=200
)

print(f"Mean MSC: {result.mean_msc:.4f}")
print(f"Sparsify: {result.sparsify}, Surrogates: {result.n_surrogates}")
```

### LRG Analysis

```python
from lrg_eegfc import compute_corr_matrix, compute_lrg_analysis

# Get FC matrix
fc_result = compute_corr_matrix("Pat_02", "rsPre", "beta")

# Compute LRG analysis
lrg_result = compute_lrg_analysis(
    fc_result.adjacency_matrix,
    "Pat_02", "rsPre", "beta", "corr"
)

print(f"Ultrametric matrix shape: {lrg_result.ultrametric_matrix.shape}")
print(f"Optimal threshold: {lrg_result.optimal_threshold:.4f}")
print(f"Giant component nodes: {lrg_result.n_nodes}")
```

### Comparisons

```python
from lrg_eegfc import compare_fc_methods, compare_phases

# Compare MSC vs correlation
comp = compare_fc_methods("Pat_02", "rsPre", "beta")
if comp:
    print(f"Matrix distance: {comp.matrix_distance:.4f}")
    print(f"Tree similarity: {comp.tree_similarity:.4f}")

# Compare pre vs post (memory effects)
comp = compare_phases("Pat_02", "rsPre", "rsPost", "beta", "corr")
if comp:
    print(f"Pre-Post distance: {comp.matrix_distance:.4f}")
```

### Batch Processing

```python
from lrg_eegfc import batch_compare_fc_methods, batch_compare_phases
from lrg_eegfc.config.const import BRAIN_BANDS, PHASE_LABELS

# Compare MSC vs corr across all combinations
df = batch_compare_fc_methods(
    ["Pat_02", "Pat_03"],
    list(PHASE_LABELS),
    list(BRAIN_BANDS.keys())
)

# Analyze by band
print(df.groupby("band")["matrix_distance"].mean())

# Memory effects analysis
phase_pairs = [("rsPre", "rsPost"), ("taskLearn", "taskTest")]
df_memory = batch_compare_phases(
    ["Pat_02", "Pat_03"],
    phase_pairs,
    ["beta", "alpha"],
    "corr"
)

print(df_memory.groupby("phase_pair")["matrix_distance"].mean())
```

## Data Organization

```
project/
├── data/
│   ├── stereoeeg_patients/          # Raw data
│   │   ├── Pat_02/
│   │   │   ├── rsPre.mat
│   │   │   ├── taskLearn.mat
│   │   │   ├── taskTest.mat
│   │   │   ├── rsPost.mat
│   │   │   └── channel_labels.csv
│   │   ├── Pat_03/
│   │   └── ...
│   ├── corr_cache/                  # Correlation matrices
│   │   ├── Pat_02/
│   │   │   └── beta_rsPre_corr_ftype-abs_zdiag-True.npy
│   │   └── ...
│   ├── msc_cache/                   # MSC matrices
│   │   ├── Pat_02/
│   │   │   └── beta_rsPre_msc_sparsify-none_nperseg-256.npy
│   │   └── ...
│   └── lrg_cache/                   # LRG analysis results
│       ├── Pat_02/
│       │   ├── beta_rsPre_lrg_corr.npz
│       │   └── beta_rsPre_lrg_msc.npz
│       └── ...
└── results/                          # Comparison outputs
    ├── msc_vs_corr.csv
    └── memory_effects.csv
```

## Frequency Bands

Default brain frequency bands:

- **delta**: 0.5 - 4 Hz
- **theta**: 4 - 8 Hz
- **alpha**: 8 - 13 Hz
- **beta**: 13 - 30 Hz
- **low_gamma**: 30 - 70 Hz
- **high_gamma**: 70 - 150 Hz

## Phases

Default recording phases:

- **rsPre**: Resting state before task
- **taskLearn**: Learning task
- **taskTest**: Testing task
- **rsPost**: Resting state after task

## Troubleshooting

### No data found errors

**Problem:** `FileNotFoundError` for patient data

**Solution:**
1. Check data is in correct format: `data/stereoeeg_patients/Pat_XX/{phase}.mat`
2. Run `python src/inspect_patient_data.py` to validate data
3. Review `PATIENT_DATA_REPORT.txt` for missing files

### Sampling rate warnings

**Problem:** `WARNING: No sampling rate found in data, using provided value: 2048 Hz`

**Solution:**
- For Pat_05-08: Add 'fs' parameter to .mat files (request from collaborators)
- Pipeline will use 2048 Hz default - verify this is correct for your data

### No cached FC matrices

**Problem:** LRG analysis fails with "No cached FC matrix"

**Solution:**
1. Compute FC matrices first:
   ```bash
   python src/compute_corr_matrices.py --patients Pat_02
   python src/compute_msc_matrices.py --patients Pat_02
   ```
2. Then run LRG analysis:
   ```bash
   python src/compute_lrg_analysis.py --fc-method corr --patients Pat_02
   ```

### Empty comparison results

**Problem:** `compare_fc_methods.py` returns no results

**Solution:**
- Ensure both MSC and corr LRG analyses are cached
- Run both:
  ```bash
  python src/compute_lrg_analysis.py --fc-method corr
  python src/compute_lrg_analysis.py --fc-method msc
  ```

## Performance Tips

### Parallel Processing

Currently scripts process sequentially. For large batches, process patients in parallel:

```bash
# Terminal 1
python src/compute_corr_matrices.py --patients Pat_02 &

# Terminal 2
python src/compute_corr_matrices.py --patients Pat_03 &
```

### Caching Strategy

- FC matrices (~1MB each): Cached per patient/phase/band/parameters
- LRG results (~100KB each): Cached per patient/phase/band/fc_method
- Comparisons: Not cached (fast recomputation from LRG cache)

To recompute with different parameters:
```bash
# Different filter for correlations
python src/compute_corr_matrices.py --filter-type pos --overwrite

# Different nperseg for MSC
python src/compute_msc_matrices.py --nperseg 512 --overwrite
```

## Next Steps

After completing the pipeline:

1. **Visualization** - Plot entropy curves, dendrograms, heatmaps
2. **Statistical Testing** - Test significance of MSC vs corr differences
3. **Memory Effects** - Statistical analysis of pre-post changes
4. **Cross-Patient Analysis** - Aggregate patterns across patients
5. **Validation** - Implement surrogate-based sparsification for MSC

## References

- Laplacian Renormalization Group (LRG) methodology
- Magnitude-Squared Coherence (MSC) for frequency-specific connectivity
- Hierarchical clustering and ultrametric distances
