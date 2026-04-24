---
name: notebook-refactor-plan
type: plan
era: MSC
status: superseded
created: 2025-12-11
updated: 2026-04-24
pointers: []
---

# Notebook Refactoring Plan
**Date:** 2025-12-11
**Project:** LRG EEG Functional Connectivity
**Objective:** Comprehensive refactoring of ipynb/ folder to create clean, modular, maintainable notebooks

---

## Current State Analysis

### Notebook Inventory (32 notebooks total)

#### ✅ Recent & Working (Dec 2025 - 5 notebooks)
1. **UTILS-COHERENCE_NETWORKS.ipynb** - Coherence-based FC tutorial (Dec 11)
2. **UTILS-FC_COMPARISON_PIPELINE.ipynb** - FC methods comparison (Dec 11)
3. **FIGMNTGN03.ipynb** - Time windows analysis figures (Dec 11)
4. **FIGMNTGN04.ipynb** - Additional figures (Dec 11)
5. **tmp.ipynb** - Temporary work (Dec 10)

#### ⚠️ Medium Age - Likely Working (Oct-Nov 2025 - 10 notebooks)
6. **DSTCMP_all_distance_measures.ipynb** - Comprehensive distance comparison (Nov 25)
7. **UTILS-FILEREADER.ipynb** - Simple data loader demo (Nov 25)
8. **DSTCMP_permutation_robust.ipynb** - Permutation-robust distance (Oct 22)
9. **NEW_standard_disXpat.ipynb** - Standard distance patient comparison (Oct 22)
10. **NEW_ultrametric_scaled_distance_pat_comparison.ipynb** (Oct 22)
11. **NEW_ultrametric_quantile_rmse_pat_comparison.ipynb** (Oct 22)
12. **NEW_ultrametric_rank_correlation_pat_comparison.ipynb** (Oct 22)
13. **NEW_tree_measures_pat_comparison.ipynb** (Oct 22)
14. **NEW_distance_of_distances.ipynb** - Distance of distances analysis (Oct 22)
15. **FIGMNTGN02.ipynb** - Figure generation (Oct 8)

#### ❌ Old - Likely Outdated (May-Sep 2025 - 17 notebooks)
16. **FIGMNTGN01.ipynb** - Early figure generation (Oct 8)
17. **TEST_distance_of_distances.ipynb** - Exploratory (Oct 8-9)
18. **TEST_distance_of_distances_2.ipynb** - Exploratory v2 (Oct 9)
19. **TEST_interactive_single_patient.ipynb** - Interactive exploration (Oct 8)
20. **single_patient_experiments.ipynb** - Single patient analysis (Oct 8)
21. **distance_of_distances.ipynb** - Original distance analysis (Oct 3)
22. **poster02.ipynb** - Poster figures (Sep 25)
23. **TEST_improve_speed.ipynb** - Performance testing (Sep 25)
24. **TEST_per_patient_analysis.ipynb** - Per-patient exploration (Sep 25)
25. **TEST_per_patient_time_windows.ipynb** - Time window analysis (Sep 10)
26. **MLB-F01.ipynb** - Unknown analysis (Sep 10)
27. **poster01.ipynb** - Poster figures (Jul 22)
28. **tests.ipynb** - General testing (Jun 24)
29. **pat02_corrnet_bands.ipynb** - Very exploratory band analysis (Jun 19)
30. **band_pca.ipynb** - PCA analysis (Jun 11)
31. **specific_heat_animation.ipynb** - Animation experiments (Jun 10)
32. **pat02_corrnet_emd.ipynb** - EMD experiments (May 14)

---

## Key Issues Identified

### 1. **Code Duplication**
- Data loading code repeated in 20+ notebooks
- Graph plotting functions (dendrograms, networks) duplicated everywhere
- Network processing pipelines copy-pasted across notebooks
- LRG analysis workflows reimplemented multiple times

### 2. **Helper Functions in Notebooks**
Critical reusable functions found in notebooks that should be in `src/`:

**From TEST_per_patient_analysis.ipynb:**
- `plot_dendrogram()` - dendrogram plotting with consistent styling
- `plot_network()` - network visualization with dendrogram colors
- `handle_empty_phase()` - graceful handling of empty networks
- `process_network_for_phase()` - complete phase processing pipeline

**From FIGMNTGN03.ipynb:**
- `plot_metastable_sankey()` - Sankey diagrams for metastable evolution
- `get_dendrogram_consistent_clusters()` - cluster assignments matching dendrogram
- `fcluster_with_outliers()` - outlier-aware hierarchical clustering
- Clustering consistency analysis functions

**From pat02_corrnet_bands.ipynb:**
- `plot_graph_analysis()` - comprehensive graph analysis figure (dendrogram + network + entropy + adjacency)
- Marchenko-Pastur distribution utilities
- Correlation matrix cleaning functions

**From DSTCMP_all_distance_measures.ipynb:**
- `compute_phase_distance_matrix()` - distance matrix computation wrapper
- `plot_measure_across_patients()` - multi-patient heatmap plotting
- `compute_cross_patient_consistency()` - cross-patient comparison metrics

### 3. **Poor Organization**
- No clear naming convention (UTILS-, FIGMNTGN, TEST_, NEW_ all mixed)
- Notebooks serve multiple purposes (exploration + analysis + figure generation)
- No separation between tutorials, analyses, and paper figures

### 4. **Configuration Issues**
- Parameters hardcoded throughout cells
- No centralized configuration section
- Jump indices, thresholds, band selections scattered

### 5. **Import Inconsistencies**
- Some use `from lrg_eegfc.notebook import *`
- Others manually import individual modules
- Inconsistent use of `move_to_rootf()`

---

## Refactoring Strategy

### Phase 1: Code Extraction to src/ (Priority: HIGH)

Create new modules in `src/lrg_eegfc/` to house reusable notebook code:

#### A. `src/lrg_eegfc/visuals/` (NEW PACKAGE)
**Purpose:** Consolidate all visualization functions

**Files to create:**
1. **`src/lrg_eegfc/visuals/__init__.py`**
   - Re-export main plotting functions

2. **`src/lrg_eegfc/visuals/dendrograms.py`**
   ```python
   def plot_dendrogram(ax, linkage_matrix, threshold, labels, **kwargs)
   def plot_dendrogram_with_inset_network(linkage_matrix, graph, labels, **kwargs)
   def get_dendrogram_node_colors(dendrogram_result, label_map)
   ```

3. **`src/lrg_eegfc/visuals/networks.py`**
   ```python
   def plot_network(ax, graph, labels, colors=None, layout='spring', **kwargs)
   def plot_network_colored_by_dendrogram(ax, graph, dendrogram_result, labels, **kwargs)
   def plot_3d_network_with_positions(graph, positions, labels, **kwargs)
   ```

4. **`src/lrg_eegfc/visuals/combined.py`**
   ```python
   def plot_graph_analysis_panel(graph, labels, tau_range=(-3, 5), **kwargs)
   # Combines: dendrogram + network + entropy curves + adjacency matrix
   ```

5. **`src/lrg_eegfc/visuals/heatmaps.py`**
   ```python
   def plot_phase_distance_heatmap(distance_matrix, phase_labels, **kwargs)
   def plot_multi_patient_heatmaps(measure_results, patients, bands, **kwargs)
   def plot_correlation_matrix(corr_matrix, labels, **kwargs)
   ```

6. **`src/lrg_eegfc/visuals/sankey.py`**
   ```python
   def plot_metastable_sankey(graph, tau_list, Trho_list, **kwargs)
   # For visualizing metastable state evolution
   ```

7. **`src/lrg_eegfc/visuals/diagnostics.py`**
   ```python
   def plot_marchenko_pastur(eigenvalues, gamma, **kwargs)
   def plot_threshold_statistics(graph, thresholds, **kwargs)
   def plot_entropy_curves(graph, tau_range, **kwargs)
   ```

#### B. `src/lrg_eegfc/workflow_*.py` (EXTEND EXISTING)
**Purpose:** High-level analysis pipelines

**Enhancements:**
1. **`src/lrg_eegfc/workflow_corr.py`** (already exists, enhance)
   - Add `process_phase_correlation_network()` - full pipeline for one phase
   - Add `process_patient_all_phases()` - iterate over phases
   - Add `process_multi_patient_analysis()` - cross-patient workflows

2. **`src/lrg_eegfc/workflow_comparison.py`** (NEW)
   ```python
   def compare_fc_methods(data, fs, bands, methods=['correlation', 'coherence'])
   def compute_all_distance_measures(ultra_dict, linkage_dict, condensed_dict)
   def cross_patient_consistency_analysis(measure_results, patients)
   ```

3. **`src/lrg_eegfc/workflow_time_windows.py`** (NEW)
   ```python
   def split_into_time_windows(data, n_windows)
   def compute_fc_per_window(data, fs, band, n_windows, method='correlation')
   def analyze_temporal_evolution(data, fs, bands, n_windows)
   ```

#### C. `src/lrg_eegfc/utils/clustering/` (NEW PACKAGE)
**Purpose:** Clustering and community detection utilities

**Files:**
1. **`src/lrg_eegfc/utils/clustering/hierarchical.py`**
   ```python
   def fcluster_with_outliers(linkage_matrix, threshold, **kwargs)
   def get_dendrogram_consistent_clusters(linkage_matrix, dendrogram_result, threshold)
   def compute_optimal_clusters(linkage_matrix, method='auto')
   ```

2. **`src/lrg_eegfc/utils/clustering/consistency.py`**
   ```python
   def align_partitions_across_tau(partitions)
   def compute_cluster_stability(clustering_results)
   def track_node_trajectories(clusterings)
   ```

#### D. `src/lrg_eegfc/utils/distances/` (NEW PACKAGE)
**Purpose:** Distance measure computation and comparison

**Files:**
1. **`src/lrg_eegfc/utils/distances/comparison.py`**
   ```python
   def compute_all_distance_measures(U1, U2, Z1, Z2, D1, D2, labels)
   def compute_phase_distance_matrix(data_dict, patient, band, measure_fn)
   def cross_patient_consistency(M1, M2)
   ```

2. **`src/lrg_eegfc/utils/distances/metrics.py`**
   ```python
   # Re-export from lrgsglib or implement wrappers
   ```

#### E. `src/lrg_eegfc/notebook.py` (ENHANCE EXISTING)
**Purpose:** Notebook convenience imports and setup

**Enhancements:**
```python
# Current exports (keep)
from lrg_eegfc.config.const import *
from lrg_eegfc.utils.datamanag.loaders import load_data_dict
# ... existing imports ...

# ADD NEW:
# Visualization
from lrg_eegfc.visuals.dendrograms import plot_dendrogram
from lrg_eegfc.visuals.networks import plot_network
from lrg_eegfc.visuals.combined import plot_graph_analysis_panel
from lrg_eegfc.visuals.heatmaps import plot_phase_distance_heatmap

# Workflows
from lrg_eegfc.workflow_corr import process_phase_correlation_network
from lrg_eegfc.workflow_comparison import compare_fc_methods

# Clustering
from lrg_eegfc.utils.clustering.hierarchical import fcluster_with_outliers

# Helper for notebook setup
def setup_notebook_environment(figure_dir='figures', **kwargs):
    """One-line notebook setup: imports, paths, matplotlib config"""
    import matplotlib.pyplot as plt
    from pathlib import Path
    %matplotlib inline

    # Create figure directory
    path_figs = Path('data') / figure_dir
    path_figs.mkdir(parents=True, exist_ok=True)

    # Configure matplotlib
    plt.rcParams['figure.dpi'] = kwargs.get('dpi', 150)
    plt.rcParams['figure.figsize'] = kwargs.get('figsize', (10, 6))

    return path_figs
```

---

### Phase 2: Notebook Reorganization (Priority: HIGH)

#### Step 1: Archive Current Notebooks
```bash
mkdir -p ipynb/archive_2025-12-11
mv ipynb/*.ipynb ipynb/archive_2025-12-11/
```

#### Step 2: Create New Organized Structure
```
ipynb/
├── tutorials/               # How to use the package
│   ├── 01_data_loading.ipynb
│   ├── 02_correlation_fc.ipynb
│   ├── 03_coherence_fc.ipynb
│   ├── 04_fc_comparison.ipynb
│   └── 05_lrg_analysis.ipynb
│
├── analyses/               # Scientific analyses
│   ├── distance_measures/
│   │   ├── all_measures_comparison.ipynb
│   │   ├── permutation_robust.ipynb
│   │   └── cross_patient_consistency.ipynb
│   │
│   ├── temporal/
│   │   ├── time_windows_analysis.ipynb
│   │   └── phase_transitions.ipynb
│   │
│   └── single_patient/
│       ├── patient_overview.ipynb
│       └── band_specific_analysis.ipynb
│
├── figures/                # Paper/presentation figures
│   ├── figure_01_fc_overview.ipynb
│   ├── figure_02_time_windows.ipynb
│   ├── figure_03_distance_comparison.ipynb
│   └── figure_04_metastable_states.ipynb
│
└── dev/                    # Development/testing
    └── scratch.ipynb       # Temporary work
```

#### Step 3: Create New Polished Notebooks

**All notebooks should follow this template:**

```python
# %% [markdown]
# # Notebook Title
# **Purpose:** Clear one-sentence purpose
# **Inputs:** List required data/parameters
# **Outputs:** List generated outputs/figures
# **Date:** YYYY-MM-DD

# %% Configuration
# All parameters in one place at the top
CONFIG = {
    'patients': ['Pat_02', 'Pat_03'],
    'phases': ['rsPre', 'taskLearn', 'taskTest', 'rsPost'],
    'bands': ['delta', 'theta', 'alpha', 'beta', 'low_gamma', 'high_gamma'],
    'jump_index': 0,
    'n_surrogates': 50,
    'output_dir': 'data/figures/notebook_name',
}

# %% Setup
from lrgsglib import move_to_rootf
move_to_rootf(pathname='lrg_eegfc')

from lrg_eegfc.notebook import *
from pathlib import Path

# Create output directory
path_figs = Path(CONFIG['output_dir'])
path_figs.mkdir(parents=True, exist_ok=True)

# %% Load Data
data_dict, int_label_map = load_data_dict(pat_list=CONFIG['patients'])
print(f"✓ Loaded {len(data_dict)} patients")

# %% Analysis
# High-level function calls, minimal code
# All complex logic in src/

# %% Visualization
# Clean plotting code using src/lrg_eegfc/visuals functions

# %% Save Results
# Save figures, export data if needed
```

---

### Phase 3: Specific Notebook Creation (Priority: MEDIUM)

#### Tutorials (ipynb/tutorials/)

**1. `01_data_loading.ipynb`**
- Source: UTILS-FILEREADER.ipynb + common loading patterns
- Show: Load patient data, inspect structure, access metadata

**2. `02_correlation_fc.ipynb`**
- Source: Best practices from recent notebooks
- Show: Full correlation-based FC pipeline with percolation thresholds

**3. `03_coherence_fc.ipynb`**
- Source: UTILS-COHERENCE_NETWORKS.ipynb (already excellent)
- Minor cleanup: use new visualization functions

**4. `04_fc_comparison.ipynb`**
- Source: UTILS-FC_COMPARISON_PIPELINE.ipynb (already excellent)
- Show: Side-by-side comparison of correlation vs coherence

**5. `05_lrg_analysis.ipynb`**
- Source: Common LRG patterns from multiple notebooks
- Show: Laplacian → distances → clustering → visualization

#### Analyses (ipynb/analyses/)

**distance_measures/all_measures_comparison.ipynb**
- Source: DSTCMP_all_distance_measures.ipynb
- Cleanup: Extract helper functions to src/, use CONFIG section
- Focus: Compare all 9 distance measures across patients/bands

**distance_measures/permutation_robust.ipynb**
- Source: DSTCMP_permutation_robust.ipynb
- Focus: Deep dive into permutation-robust distance measure

**distance_measures/cross_patient_consistency.ipynb**
- Source: Consistency analysis from DSTCMP_all_distance_measures.ipynb
- Focus: Which measures are most reliable across patients

**temporal/time_windows_analysis.ipynb**
- Source: FIGMNTGN03.ipynb + TEST_per_patient_time_windows.ipynb
- Focus: Sliding window FC analysis, temporal evolution

**temporal/phase_transitions.ipynb**
- Source: TEST_per_patient_analysis.ipynb patterns
- Focus: Pre/post task comparisons, phase transitions

**single_patient/patient_overview.ipynb**
- Source: single_patient_experiments.ipynb + patterns from TEST notebooks
- Focus: Complete single-patient analysis workflow

**single_patient/band_specific_analysis.ipynb**
- Source: pat02_corrnet_bands.ipynb (cleaned up)
- Focus: Deep dive into band-specific network properties

#### Figures (ipynb/figures/)

**figure_01_fc_overview.ipynb**
- Source: FIGMNTGN01.ipynb (if relevant, otherwise new)
- Purpose: Publication-ready overview figure

**figure_02_time_windows.ipynb**
- Source: FIGMNTGN03.ipynb
- Purpose: Time window analysis figure for paper

**figure_03_distance_comparison.ipynb**
- Source: Best from DSTCMP notebooks
- Purpose: Distance measure comparison figure

**figure_04_metastable_states.ipynb**
- Source: Sankey diagram code from FIGMNTGN03.ipynb
- Purpose: Metastable state evolution figure

---

### Phase 4: Testing & Validation (Priority: MEDIUM)

#### Test Each New Notebook
1. Run from top to bottom without errors
2. Verify all outputs are generated
3. Check figure quality and consistency
4. Ensure no hardcoded paths

#### Create Notebook Test Suite
```python
# tests/test_notebooks.py
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

def test_tutorial_notebooks():
    """Ensure all tutorial notebooks execute successfully"""
    notebooks = [
        'ipynb/tutorials/01_data_loading.ipynb',
        'ipynb/tutorials/02_correlation_fc.ipynb',
        # ...
    ]
    for nb_path in notebooks:
        # Execute notebook and check for errors
        pass
```

---

### Phase 5: Documentation (Priority: LOW)

#### Update README.md
Add section:
```markdown
## Notebooks

### Tutorials (`ipynb/tutorials/`)
Learn how to use the package with step-by-step examples:
- `01_data_loading.ipynb` - Load and inspect SEEG data
- `02_correlation_fc.ipynb` - Correlation-based FC networks
- ...

### Analyses (`ipynb/analyses/`)
Scientific analysis notebooks organized by topic:
- `distance_measures/` - Distance measure comparisons
- `temporal/` - Time-dependent analyses
- ...

### Figures (`ipynb/figures/`)
Notebooks that generate publication figures.

### Archive (`ipynb/archive_2025-12-11/`)
Old notebooks kept for reference.
```

#### Create `ipynb/README.md`
Document:
- Notebook organization
- How to create new notebooks (follow template)
- Coding standards for notebooks

---

## Implementation Priority

### Week 1: Foundation (HIGH PRIORITY)
1. ✅ Create `src/lrg_eegfc/visuals/` package with all visualization modules
2. ✅ Extract helper functions from notebooks → src/
3. ✅ Enhance `src/lrg_eegfc/notebook.py` with new imports
4. ✅ Test new modules in a scratch notebook

### Week 2: Notebook Migration (HIGH PRIORITY)
5. ✅ Archive current notebooks to `ipynb/archive_2025-12-11/`
6. ✅ Create new directory structure
7. ✅ Create tutorial notebooks (priority: 03, 04 are already good)
8. ✅ Create 2-3 key analysis notebooks

### Week 3: Figures & Polish (MEDIUM PRIORITY)
9. ✅ Create figure generation notebooks
10. ✅ Test all notebooks end-to-end
11. ✅ Update documentation

### Week 4: Cleanup (LOW PRIORITY)
12. ⚠️ Review archived notebooks for any missed code
13. ⚠️ Delete truly obsolete notebooks
14. ⚠️ Final documentation pass

---

## Code Patterns to Standardize

### Configuration Section
```python
CONFIG = {
    'patients': ['Pat_02', 'Pat_03'],
    'phases': PHASE_LABELS,  # Use constants
    'bands': BRAIN_BANDS.keys(),
    'correlation_protocol': dict(filter_type='abs', spectral_cleaning=False),
    'jump_index': 0,
    'output_dir': 'data/figures/notebook_name',
}
```

### Data Loading
```python
# Standard pattern
data_dict, int_label_map = load_data_dict(pat_list=CONFIG['patients'])
```

### Figure Saving
```python
# Use pathlib, descriptive names
fig.savefig(path_figs / f'{patient}_{band}_{phase}_analysis.pdf',
            dpi=300, bbox_inches='tight')
print(f'✓ Saved: {path_figs / f"{patient}_{band}_{phase}_analysis.pdf"}')
```

### Progress Indicators
```python
# Use checkmarks for clarity
print('✓ Data loaded')
print('✓ Networks computed')
print('✓ Figures saved')
```

---

## Success Criteria

### By End of Refactoring:
1. ✅ All reusable code moved to `src/lrg_eegfc/`
2. ✅ Notebooks are < 200 lines on average (high-level only)
3. ✅ Clear separation: tutorials / analyses / figures
4. ✅ All notebooks have CONFIG section at top
5. ✅ No code duplication across notebooks
6. ✅ All notebooks run without errors
7. ✅ Figure output locations consistent
8. ✅ Documentation updated

### Maintenance Going Forward:
- New analyses → new notebook in appropriate folder
- New reusable code → add to src/ first, then use in notebook
- Keep tutorials updated as API evolves
- Archive old notebooks rather than delete

---

## Risk Mitigation

### Preserve Existing Work
- ✅ Archive ALL current notebooks (never delete)
- ✅ Keep archive accessible for reference
- ✅ Document what was migrated where

### Incremental Migration
- ✅ Don't try to migrate everything at once
- ✅ Start with most important notebooks
- ✅ Test each new notebook before moving to next

### Version Control
- ✅ Commit after each major step
- ✅ Tag before archive: `git tag notebook-archive-2025-12-11`
- ✅ Can always revert if needed

---

## Notes

### Notebooks to Definitely Keep (Migrate)
1. UTILS-COHERENCE_NETWORKS.ipynb → tutorials/03_coherence_fc.ipynb
2. UTILS-FC_COMPARISON_PIPELINE.ipynb → tutorials/04_fc_comparison.ipynb
3. DSTCMP_all_distance_measures.ipynb → analyses/distance_measures/all_measures_comparison.ipynb
4. FIGMNTGN03.ipynb → figures/figure_02_time_windows.ipynb
5. FIGMNTGN04.ipynb → figures/figure_04_*.ipynb (depending on content)

### Notebooks to Maybe Merge/Consolidate
- NEW_* distance notebooks → Single comprehensive distance analysis
- TEST_distance_of_distances* → Keep insights, discard exploration
- poster*.ipynb → Extract good figures to figures/ folder

### Notebooks to Likely Discard (After Code Extraction)
- tests.ipynb (general testing, use pytest instead)
- tmp.ipynb (temporary)
- pat02_corrnet_emd.ipynb (very old, likely superseded)
- specific_heat_animation.ipynb (exploratory, extract if useful)
- band_pca.ipynb (extract if PCA still relevant)

---

## Open Questions

1. **PCA Analysis**: Is band_pca.ipynb still relevant? If yes, migrate to analyses/
2. **EMD Analysis**: Should EMD be kept? (pat02_corrnet_emd.ipynb)
3. **Poster Figures**: Keep poster notebooks or just extract useful code?
4. **MLB-F01.ipynb**: What is this? Should it be kept?
5. **Interactive Notebooks**: TEST_interactive_single_patient.ipynb - is interactivity needed?

---

## End of Plan

**Next Steps:**
1. Review and approve this plan
2. Begin Phase 1: Code extraction to src/
3. Test extracted modules
4. Proceed with Phase 2: Notebook reorganization
