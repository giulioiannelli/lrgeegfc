---
name: refactor-plan-updated
type: plan
era: MSC
status: superseded
created: 2025-12-11
updated: 2026-04-24
pointers: []
---

# Notebook Refactoring Plan (UPDATED)
**Date:** 2025-12-11
**Status:** Ready to implement
**User Decisions Incorporated:** ✅

---

## User Decisions

1. **Phase 1 Scope:** Extract all possible, but check existing code first (no duplication)
2. **Old Notebooks:** Keep only interactive notebooks, discard rest after code extraction
3. **Top 3 Priorities:**
   - UTILS-COHERENCE_NETWORKS.ipynb
   - DSTCMP_all_distance_measures.ipynb
   - FIGMNTGN03.ipynb
4. **Organization:** Simplified structure - `tutorials/`, `notebooks/`, `dev/`

---

## Audit Results: What Already Exists ✅

### Existing src/lrg_eegfc/ Infrastructure

**Core Modules:**
- `plotting.py` - Core plotting (correlation matrix, dendrogram, entropy, graph)
- `notebook.py` - Notebook convenience imports (needs enhancement)
- `workflow.py`, `workflow_corr.py`, `workflow_msc.py`, `workflow_lrg.py` - High-level workflows
- `compare.py` - FC method comparison

**Visualization Package (`src/lrg_eegfc/visuals/`):**
- ✅ `correlation.py` - Correlation heatmaps, network plots, Marchenko-Pastur, percolation
- ✅ `msc.py` - MSC/coherence heatmaps and networks
- ✅ `lrg.py` - LRG entropy, dendrograms, ultrametric heatmaps, full panels
- ✅ `metastable.py` - **Sankey diagrams, cluster evolution tracking**
- ✅ `reorganization.py` - Phase reorganization visualization
- ✅ `compare.py` - FC comparison plots

**Data Management (`src/lrg_eegfc/utils/datamanag/`):**
- ✅ `loaders.py` - Data loading (load_data_dict)
- ✅ `patient.py` - Patient data structures

**Correlation FC (`src/lrg_eegfc/utils/corrmat/`):**
- ✅ `base.py` - Core correlation computations
- ✅ `thresholds.py` - Percolation thresholds
- ✅ `network.py` - Network processing
- ✅ `bands.py` - Per-band analysis
- ✅ `structures.py` - Batch processing (compute_structures_for_patient)

**Coherence FC (`src/lrg_eegfc/utils/coherence/`):**
- ✅ `msc.py` - Magnitude-squared coherence
- ✅ `surrogates.py` - Circular shift surrogates
- ✅ `sparsify.py` - Soft sparsification

---

## What's Missing (Need to Add) ❌

### 1. Clustering Utilities (NEW: `src/lrg_eegfc/utils/clustering/`)

**File: `hierarchical.py`**
```python
def fcluster_with_outliers(linkage_matrix, threshold, outlier_sensitivity=1.5, min_cluster_size=2):
    """Enhanced fcluster that identifies outliers (from FIGMNTGN03)"""

def get_dendrogram_consistent_clusters(linkage_matrix, dendro_result, threshold):
    """Cluster assignments matching dendrogram colors (from FIGMNTGN03)"""

def compute_optimal_clusters_auto(linkage_matrix, method='auto'):
    """Wrapper with multiple methods for optimal clustering"""
```

### 2. Distance Comparison Utilities (NEW: `src/lrg_eegfc/utils/distances/`)

**File: `comparison.py`**
```python
def compute_phase_distance_matrix(pat, band, phase_labels, compute_fn):
    """Compute distance matrix across phases (from DSTCMP_all)"""

def compute_cross_patient_consistency(M1, M2):
    """Pearson, Spearman, normalized diff (from DSTCMP_all)"""

def rank_distance_measures(measure_results, patients):
    """Rank measures by cross-patient consistency (from DSTCMP_all)"""
```

**File: `plotting.py`** (extend existing)
```python
def plot_phase_distance_heatmap(distance_matrix, phase_labels, **kwargs):
    """Heatmap for phase-phase distances"""

def plot_multi_patient_heatmaps(measure_results, patients, bands, **kwargs):
    """Multi-patient comparison grid (from DSTCMP_all)"""
```

### 3. Enhance `notebook.py`

Currently very minimal - enhance with convenient imports:
```python
# Keep existing
from .config import *
from .utils import *

# ADD:
# Visualization
from .visuals import (
    plot_correlation_heatmap,
    plot_msc_heatmap,
    plot_lrg_dendrogram,
    plot_lrg_full_panel,
    create_sankey_diagram,
    plot_fc_comparison,
)

# Workflows
from .workflow_corr import compute_band_connectivity
from .compare import compare_fc_methods

# Data loading
from .utils.datamanag.loaders import load_data_dict

# Constants
from .config.const import BRAIN_BANDS, PHASE_LABELS

# Utility: Notebook setup
def setup_notebook(figure_dir='figures', **kwargs):
    """One-line notebook setup"""
    from pathlib import Path
    path_figs = Path('data') / figure_dir
    path_figs.mkdir(parents=True, exist_ok=True)
    print(f'✓ Notebook setup complete')
    print(f'  Figure output: {path_figs}')
    return path_figs
```

---

## Revised Implementation Plan

### Phase 1: Add Missing Functions (1-2 days)

#### Step 1.1: Create `src/lrg_eegfc/utils/clustering/`
- `hierarchical.py` - Add outlier-aware clustering functions from FIGMNTGN03
- `__init__.py` - Export public API

#### Step 1.2: Create `src/lrg_eegfc/utils/distances/`
- `comparison.py` - Add distance comparison utilities from DSTCMP_all
- `plotting.py` - Add distance heatmap functions
- `__init__.py` - Export public API

#### Step 1.3: Enhance `src/lrg_eegfc/notebook.py`
- Add comprehensive imports from existing modules
- Add `setup_notebook()` convenience function
- Update docstring with usage examples

#### Step 1.4: Update `src/lrg_eegfc/visuals/__init__.py`
- Ensure all useful functions are exported
- Add any missing re-exports

---

### Phase 2: Reorganize Notebooks (1-2 days)

#### Step 2.1: Archive Current Notebooks
```bash
mkdir -p ipynb/archive_2025-12-11
# Move ALL current notebooks (preserve everything!)
mv ipynb/*.ipynb ipynb/archive_2025-12-11/
# Keep only the archive folder
```

#### Step 2.2: Create New Structure
```
ipynb/
├── archive_2025-12-11/     # All current notebooks preserved
├── tutorials/               # Learning materials
│   ├── 01_data_loading.ipynb
│   ├── 02_correlation_fc.ipynb
│   ├── 03_coherence_fc.ipynb
│   ├── 04_fc_comparison.ipynb
│   └── 05_lrg_analysis.ipynb
├── notebooks/              # All scientific analyses
│   ├── distance_measures_comparison.ipynb
│   ├── time_windows_analysis.ipynb
│   ├── phase_reorganization.ipynb
│   ├── metastable_states_sankey.ipynb
│   └── interactive_single_patient.ipynb  # User requested
└── dev/                    # Scratch work
    └── scratch.ipynb
```

---

### Phase 3: Migrate Priority Notebooks (2-3 days)

#### Priority 1: UTILS-COHERENCE_NETWORKS.ipynb → `tutorials/03_coherence_fc.ipynb`
**Status:** Already excellent, minimal changes needed
- ✅ Uses `coherence_fc_pipeline()` from src
- ✅ Well-documented with markdown
- ✅ Clear examples
- **Changes:** Add CONFIG section, use `setup_notebook()`, minor cleanup

#### Priority 2: DSTCMP_all_distance_measures.ipynb → `notebooks/distance_measures_comparison.ipynb`
**Key extractions:**
- Move `compute_phase_distance_matrix()` → `src/lrg_eegfc/utils/distances/comparison.py`
- Move `compute_cross_patient_consistency()` → same
- Move `plot_measure_across_patients()` → `src/lrg_eegfc/utils/distances/plotting.py`
- **Result:** Notebook becomes ~100 lines of high-level calls

#### Priority 3: FIGMNTGN03.ipynb → `notebooks/time_windows_analysis.ipynb`
**Key extractions:**
- Move `fcluster_with_outliers()` → `src/lrg_eegfc/utils/clustering/hierarchical.py`
- Move `get_dendrogram_consistent_clusters()` → same
- Move Sankey diagram code → already in `visuals/metastable.py`! ✅
- **Result:** Notebook uses existing `create_sankey_diagram()` function

---

### Phase 4: Create Additional Notebooks (1-2 days)

#### Tutorial Notebooks (from best practices across archive)

**`tutorials/01_data_loading.ipynb`** (NEW)
- Source: Common patterns + UTILS-FILEREADER.ipynb
- Show: `load_data_dict()`, inspect structure, access metadata

**`tutorials/02_correlation_fc.ipynb`** (NEW)
- Source: Best correlation practices from archive
- Show: Full correlation pipeline using `workflow_corr.py` functions

**`tutorials/04_fc_comparison.ipynb`** (MIGRATE)
- Source: UTILS-FC_COMPARISON_PIPELINE.ipynb (already good)
- Minor cleanup: CONFIG section, `setup_notebook()`

**`tutorials/05_lrg_analysis.ipynb`** (NEW)
- Source: Common LRG patterns
- Show: Complete LRG workflow using `plot_lrg_full_panel()`

#### Analysis Notebooks

**`notebooks/metastable_states_sankey.ipynb`** (NEW)
- Source: FIGMNTGN03 Sankey code
- Uses: `create_sankey_diagram()` from `visuals/metastable.py`

**`notebooks/phase_reorganization.ipynb`** (NEW)
- Source: Reorganization patterns from archive
- Uses: `visuals/reorganization.py` functions

**`notebooks/interactive_single_patient.ipynb`** (MIGRATE)
- Source: TEST_interactive_single_patient.ipynb
- Keep interactive widgets, use src functions

---

## Standard Notebook Template

Every new/migrated notebook follows this structure:

```python
# %% [markdown]
# # Notebook Title
#
# **Purpose:** One-sentence description
# **Inputs:** Data requirements
# **Outputs:** Generated files/figures
# **Date:** 2025-12-11

# %% Configuration
CONFIG = {
    'patients': ['Pat_02', 'Pat_03'],
    'phases': ['rsPre', 'taskLearn', 'taskTest', 'rsPost'],
    'bands': ['delta', 'theta', 'alpha', 'beta', 'low_gamma', 'high_gamma'],
    'output_dir': 'data/figures/notebook_name',
    # ... all parameters here
}

# %% Setup
from lrgsglib import move_to_rootf
move_to_rootf(pathname='lrg_eegfc')

from lrg_eegfc.notebook import *
path_figs = setup_notebook(CONFIG['output_dir'])

# %% Load Data
data_dict, int_label_map = load_data_dict(pat_list=CONFIG['patients'])
print(f"✓ Loaded {len(data_dict)} patients")

# %% Analysis
# High-level function calls only, < 200 lines total

# %% Save
# Save figures with descriptive names
```

---

## Timeline

### Week 1 (Days 1-2): Add Missing Functions ⚡ HIGH PRIORITY
- Day 1: Create `utils/clustering/` and `utils/distances/`
- Day 2: Enhance `notebook.py`, test new functions

### Week 1 (Days 3-5): Migrate Priority Notebooks ⚡ HIGH PRIORITY
- Day 3: Archive all current notebooks, create new structure
- Day 4: Migrate UTILS-COHERENCE → `tutorials/03_coherence_fc.ipynb`
- Day 5: Migrate DSTCMP_all → `notebooks/distance_measures_comparison.ipynb`

### Week 2 (Days 6-8): Complete Migration 📊 MEDIUM PRIORITY
- Day 6: Migrate FIGMNTGN03 → `notebooks/time_windows_analysis.ipynb`
- Day 7: Create `tutorials/01_data_loading.ipynb` and `02_correlation_fc.ipynb`
- Day 8: Create remaining tutorial notebooks

### Week 2 (Days 9-10): Polish & Documentation 📝 LOWER PRIORITY
- Day 9: Create additional analysis notebooks, test all
- Day 10: Update documentation, final review

---

## Success Criteria

✅ **Code Organization:**
- All reusable functions in `src/lrg_eegfc/`
- No code duplication across notebooks
- Notebooks < 200 lines (high-level only)

✅ **Notebook Quality:**
- All have CONFIG section at top
- All use `setup_notebook()` for initialization
- Clear markdown documentation
- Run without errors

✅ **Structure:**
- `tutorials/` - 5 step-by-step learning notebooks
- `notebooks/` - Clean scientific analysis notebooks
- `dev/` - Scratch work
- `archive_2025-12-11/` - All original work preserved

✅ **Documentation:**
- Updated README with notebook organization
- Each notebook has purpose/inputs/outputs header
- Code is self-documenting with clear variable names

---

## Risk Mitigation

1. **No Data Loss:** All current notebooks archived (never deleted)
2. **Incremental:** Migrate one notebook at a time, test each
3. **Reversible:** Git tag before changes: `git tag notebook-archive-2025-12-11`
4. **Validation:** Each migrated notebook must run end-to-end

---

## Next Immediate Actions

**Ready to implement Phase 1, Step 1.1:**

1. Create `src/lrg_eegfc/utils/clustering/hierarchical.py`
2. Extract `fcluster_with_outliers()` from FIGMNTGN03
3. Add tests
4. Update imports

**Shall I proceed?** 🚀
