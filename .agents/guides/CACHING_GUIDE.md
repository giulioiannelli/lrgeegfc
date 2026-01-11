# Caching Guide (LRG EEG FC)

This guide explains how cached outputs are produced and reused to avoid
reloading raw data or recomputing heavy steps. Use this workflow in all
notebooks and scripts.

## Why caching matters
- Raw SEEG data are large and slow to process.
- FC matrices (corr/MSC) and LRG results are expensive to recompute.
- Visualizations should read cached outputs and never recompute data.

## Cache directories (default)

- Correlation FC: `data/corr_cache/Pat_XX/`
- Cleaned correlation: `data/corr_cache/Pat_XX/` (+ `_meta.npz`)
- MSC FC: `data/msc_cache/Pat_XX/`
- LRG analysis: `data/lrg_cache/Pat_XX/`
- Figures: `data/figures/<category>/Pat_XX/`
- CLI (legacy): `data/correlations/Pat_XX/` from `lrg-eegfc-corr`
- Time-window MSC: `data/msc_cache_windows/Pat_XX/<phase>/<band>/`
- Time-window LRG: `data/lrg_cache_windows/Pat_XX/<phase>/<band>/`

## Cache keys and filenames

Correlation:
- `data/corr_cache/Pat_XX/{band}_{phase}_corr_ftype-{filter_type}_zdiag-{bool}.npy`

MSC:
- `data/msc_cache/Pat_XX/{band}_{phase}_msc_sparsify-{mode}_nsurr-{N}_nperseg-{n}.npy`

Time-window MSC (proposed):
- `data/msc_cache_windows/Pat_XX/<phase>/<band>/win_<idx>_msc_sparsify-{mode}_nsurr-{N}_nperseg-{n}.npy`

Cleaned correlation:
- `data/corr_cache/Pat_XX/{band}_{phase}_corr_cleaned.npy`
- `data/corr_cache/Pat_XX/{band}_{phase}_corr_cleaned_meta.npz`

LRG:
- `data/lrg_cache/Pat_XX/{band}_{phase}_lrg_{fc_method}.npz`

Time-window LRG (proposed):
- `data/lrg_cache_windows/Pat_XX/<phase>/<band>/win_<idx>_lrg_{fc_method}.npz`

## Load-first pattern (preferred)

Always try a cached load before computing:

```python
from lrg_eegfc import load_corr_matrix, compute_corr_matrix

corr = load_corr_matrix(patient, phase, band)
if corr is None:
    corr = compute_corr_matrix(patient, phase, band)
```

Repeat the same pattern for MSC and LRG:

```python
from lrg_eegfc import load_msc_matrix, compute_msc_matrix

msc = load_msc_matrix(patient, phase, band)
if msc is None:
    msc = compute_msc_matrix(patient, phase, band)
```

```python
from lrg_eegfc import load_lrg_result, compute_lrg_analysis

lrg = load_lrg_result(patient, phase, band, fc_method)
if lrg is None:
    lrg = compute_lrg_analysis(fc_matrix, patient, phase, band, fc_method)
```

## Compute vs visualize separation

- Compute scripts should generate caches and stop.
- Visualization scripts should fail if caches are missing.
- Notebooks should load caches first and compute only if needed.

## Cache invalidation

- Use `overwrite_cache=True` in compute functions when you must recompute.
- Keep parameters consistent (filter type, nperseg, surrogates) or you will
  write new cache files with different suffixes.

## Recommended usage in this repo

- Batch compute: `src/compute_corr_matrices.py`, `src/compute_msc_matrices.py`
- LRG compute: `src/compute_lrg_analysis.py`
- Visualizations: `src/visualize_*.py` (should read caches only)

## Notes

- `config/const.py` reads the filesystem at import time; prefer `constants.py`
  for lightweight imports in utilities.
- The CLI `lrg-eegfc-corr` uses a different output root
  (`data/correlations/Pat_XX/`). Use it for quick checks only.
