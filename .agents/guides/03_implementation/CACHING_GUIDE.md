# Caching Guide (LRG EEG FC)

This guide explains how cached outputs are produced, located, and reused.

## NOTICE: data/ folder has been reorganized

The `data/` folder now has **only 4 top-level categories**: `raw/`, `cache/`,
`reports/`, `outputs/`. Old flat paths (`data/corr_cache/`, `data/msc_cache/`,
`data/figures/`, etc.) **no longer exist** — they were moved into the new
structure. Any script still using `Path("data/corr_cache")` will break.

**Always use `config.paths` constants:**
```python
from lrg_eegfc.config.paths import (
    CORR_CACHE, MSC_CACHE, LRG_CACHE,
    IMCOH_CACHE, IMCOH_LRG_CACHE,
    SEEG_DATAPATH, FIGURES_ROOT, TABLES_ROOT,
)
```

## Why caching matters
- Raw SEEG data are large and slow to process.
- FC matrices (corr/MSC/ImCoh) and LRG results are expensive to recompute.
- Visualizations should read cached outputs and never recompute data.

## Cache directories

All cache directories live under `data/cache/` (symlinked from old flat layout):

| FC method | Cache | Path constant |
|-----------|-------|---------------|
| Correlation | `data/cache/corr/Pat_XX/` | `CORR_CACHE` |
| MSC | `data/cache/msc/Pat_XX/` | `MSC_CACHE` |
| MSC (dev) | `data/cache/msc_dev/Pat_XX/` | `MSC_DEV_CACHE` |
| ImCoh | `data/cache/imcoh/Pat_XX/` | `IMCOH_CACHE` |
| Cleaned corr | `data/cache/cleaned_corr/Pat_XX/` | `CLEANED_CORR_CACHE` |
| LRG (from MSC/corr) | `data/cache/lrg/Pat_XX/` | `LRG_CACHE` |
| LRG (from ImCoh) | `data/cache/imcoh_lrg/Pat_XX/` | `IMCOH_LRG_CACHE` |
| LRG (CREMA) | `data/cache/lrg_crema/Pat_XX/` | `LRG_CREMA_CACHE` |
| FC figure cache | `data/cache/fc_fig/` | `FC_FIG_CACHE` |
| Surrogate validation | `data/cache/surrogate_validation/` | `SURROGATE_VALIDATION_CACHE` |
| Bipolar (experimental) | `data/cache/bipolar/` | `BIPOLAR_CACHE` |
| Rescaled (experimental) | `data/cache/rescaled/` | `RESCALED_CACHE` |

Reports and investigation outputs live under `data/reports/`:

| Investigation | Path |
|--------------|------|
| Metric exploration | `data/reports/metric_exploration/` |
| Clinical application | `data/reports/clinical_application/` |
| ImCoh figures | `data/reports/imcoh/` |
| ImCoh unanimity | `data/reports/imcoh_unanimity/` |
| ImCoh VI profiles | `data/reports/imcoh_vi/` |

## Cache keys and filenames

Correlation:
- `{band}_{phase}_corr_ftype-{filter_type}_zdiag-{bool}.npy`

MSC:
- `{band}_{phase}_msc_sparsify-{mode}_nsurr-{N}_nperseg-{n}.npy`

ImCoh:
- `{band}_{phase}_imcoh_sparsify-none_nperseg-{n}.npy`
- Note: ImCoh should always use `sparsify="none"` (surrogates are wrong for ImCoh)

LRG:
- `{band}_{phase}_lrg_{fc_method}.npz`
- For ImCoh: `{band}_{phase}_lrg_imcoh.npz` (stored under `imcoh_lrg/`)

Cleaned correlation:
- `{band}_{phase}_corr_cleaned.npy`
- `{band}_{phase}_corr_cleaned_meta.npz`

## Load-first pattern (preferred)

Always try a cached load before computing:

```python
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.config.paths import MSC_CACHE, LRG_CACHE, IMCOH_CACHE, IMCOH_LRG_CACHE

# MSC
msc = load_msc_matrix(patient, phase, band, cache_root=MSC_CACHE)

# ImCoh (volume-conduction immune — preferred for community analysis)
imcoh = load_msc_matrix(patient, phase, band, cache_root=IMCOH_CACHE,
                        sparsify="none", n_surrogates=0)

# LRG from MSC
lrg = load_lrg_result(patient, phase, band, "msc", cache_root=LRG_CACHE)

# LRG from ImCoh
lrg_imcoh = load_lrg_result(patient, phase, band, "imcoh", cache_root=IMCOH_LRG_CACHE)
```

## Compute vs visualize separation

- Compute scripts generate caches and stop.
- Visualization scripts fail if caches are missing.
- Notebooks load caches first and compute only if needed.

## Cache invalidation

- Use `overwrite_cache=True` in compute functions when you must recompute.
- Keep parameters consistent (filter type, nperseg, surrogates) or you will
  write new cache files with different suffixes.
- Dev runs with `filter_time` automatically redirect to `*_dev` caches.
