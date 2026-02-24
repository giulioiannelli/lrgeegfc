# Agent Notes - lrgeegfc

Canonical agent guidance lives in `.agents/README.md`.
Root `AGENTS.md` and `CLAUDE.md` should remain identical.

---

## Quick Reference Map

Use this section to quickly locate functions and modules.

### Skills Available

Invoke with `/skillname`:
- `/figure` - Generate figure code with templates
- `/patient <id>` - Run full pipeline for patient
- `/lrg` - LRG analysis with customization
- `/cache` - Manage cache files
- `/validate` - Run tests and checks
- `/style` - Apply plot styling
- `/notebook` - Create/execute notebooks
- `/data` - Inspect patient data

### Figure Generation

| Type | Module | Key Function |
|------|--------|--------------|
| Correlation matrix | `visuals.correlation` | `plot_correlation_and_network()` |
| MSC matrix | `visuals.msc` | `plot_msc_and_network()` |
| LRG full panel | `visuals.lrg` | `plot_lrg_full_panel()` |
| Dendrogram | `visuals.lrg` | `plot_lrg_dendrogram()` |
| Entropy curves | `visuals.lrg` | `plot_lrg_entropy_curves()` |
| 3D brain (Plotly) | `visuals.spatial` | `plot_spatial_network_3d()` |
| 3D brain (nilearn) | `visuals.spatial` | `view_brain_connectome()` |
| Phase comparison | `visuals.reorganization` | `plot_phase_reorganization()` |

### Data Loading (always use cache)

| Task | Module | Function |
|------|--------|----------|
| Load correlation | `workflow.corr` | `load_corr_matrix(patient, phase, band, cache_root)` |
| Load MSC | `workflow.msc` | `load_msc_matrix(patient, phase, band, cache_root, sparsify, n_surrogates, nperseg)` |
| Load LRG result | `workflow.lrg` | `load_lrg_result(patient, phase, band, fc_method, cache_root)` |
| Load timeseries | `utils.io.patient` | `load_timeseries(patient, phase, root_path)` |
| Inspect data | `utils.io.inspect` | `inspect_patient(patient, root_path)` |

### Computing (creates cache)

| Task | Module | Function |
|------|--------|----------|
| Compute correlation | `workflow.corr` | `compute_corr_matrix(patient, phase, band, ...)` |
| Compute MSC | `workflow.msc` | `compute_msc_matrix(patient, phase, band, ...)` |
| Compute LRG | `workflow.lrg` | `compute_lrg_analysis(fc_matrix, ...)` |
| Full patient | `workflow.corr` | `compute_corr_for_patient(patient, ...)` |

### Configuration

| Item | Location | Access |
|------|----------|--------|
| Frequency bands | `config.const` | `BRAIN_BANDS`, `BRAIN_BANDS_NAMES` |
| Phase labels | `config.const` | `PHASE_LABELS` |
| LaTeX band names | `config.const` | `BRAIN_BAND_TEX_DICT` |
| Sample rate | `config.const` | `DEFAULT_SAMPLE_RATE` (2048 Hz) |
| Surrogates default | `config.const` | `DEFAULT_N_SURROGATES` (200) |

---

## Cache Locations

```
data/corr_cache/Pat_XX/{band}_{phase}_*.npy
data/msc_cache/Pat_XX/{band}_{phase}_*.npy
data/lrg_cache/Pat_XX/{band}_{phase}_{fc_method}_*.npz
data/figures/{category}/Pat_XX/*.png
```

**Naming sensitivity:**
- MSC: includes `sparsify-{none|soft}`, `nperseg-{value}`, `n_surrogates-{value}`
- Correlation: includes `filter-{none|abs}`, `zero_diag-{true|false}`

---

## Critical Invariants

1. **Never recompute in visualization** - always load from cache
2. **Standard notebook header:**
   ```python
   from lrg_eegfc.notebook import *
   move_to_root(pathname="lrgeegfc")
   ```
3. **Exclude patients for cross-phase:** Pat_06, Pat_07 (missing phases)
4. **Always close figures:** `plt.close(fig)` after saving
5. **Create output dirs:** `output_path.parent.mkdir(parents=True, exist_ok=True)`

---

## Standard Imports

```python
# Notebook header
from lrg_eegfc.notebook import *
move_to_root(pathname="lrgeegfc")

# Config
from lrg_eegfc.config import BRAIN_BANDS, PHASE_LABELS, BRAIN_BAND_TEX_DICT

# Workflows
from lrg_eegfc.workflow import (
    load_corr_matrix, compute_corr_matrix,
    load_msc_matrix, compute_msc_matrix,
    load_lrg_result, compute_lrg_analysis,
)

# Visualization
from lrg_eegfc.visuals import (
    plot_correlation_and_network,
    plot_msc_and_network,
    plot_lrg_full_panel,
    plot_phase_reorganization,
    plot_spatial_network_3d,
)

# Data loading
from lrg_eegfc.utils.io import load_timeseries, inspect_patient
```

---

## Directory Structure

```
src/lrg_eegfc/
├── __init__.py          # Package exports
├── cli.py               # Command-line interface
├── notebook.py          # Notebook utilities
├── config/              # Constants, plotting config
├── workflow/            # High-level compute + cache
│   ├── corr.py          # Correlation workflow
│   ├── msc.py           # MSC workflow
│   ├── lrg.py           # LRG analysis workflow
│   └── cleaning.py      # Marchenko-Pastur cleaning
├── visuals/             # All visualization
│   ├── correlation.py   # Correlation plots
│   ├── msc.py           # MSC plots
│   ├── lrg.py           # LRG/dendrogram plots
│   ├── spatial.py       # 3D brain visualization
│   └── reorganization.py # Phase comparison
└── utils/               # Low-level utilities
    ├── io/              # Data loading
    ├── fc/              # FC computation primitives
    ├── lrg/             # LRG utilities
    └── metrics/         # Comparison metrics
```

---

## Detailed Guides

- `.agents/guides/FUNCTION_MAP.md` - Complete function reference
- `.agents/guides/FIGURE_PATTERNS.md` - Figure templates
- `.agents/guides/CACHING_GUIDE.md` - Cache management
- `.agents/guides/MSC_METHOD_GUIDE.md` - MSC methodology
- `.agents/guides/AGENT_PLAYBOOK.md` - Session workflow
