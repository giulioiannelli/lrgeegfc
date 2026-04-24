---
name: agent-tasks
type: guide
era: CROSS_ERA
status: current
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

# Agent Task Reference

Common tasks agents should know how to perform autonomously.

---

## Task: Generate Figures for Patient

**Goal:** Create all standard figures for a patient across bands/phases.

**Steps:**
1. Check data exists: `data/stereoeeg_patients/Pat_XX/`
2. Verify caches: `data/corr_cache/Pat_XX/`, `data/msc_cache/Pat_XX/`, `data/lrg_cache/Pat_XX/`
3. Generate figures using `lrg_eegfc.visuals` functions
4. Save to `data/figures/{category}/Pat_XX/`

**Key functions:**
```python
from lrg_eegfc.visuals import (
    plot_correlation_and_network,
    plot_msc_and_network,
    plot_lrg_full_panel,
    plot_phase_reorganization,
)
```

**Skip if:** Cache files missing (report which ones)

---

## Task: Run Full Pipeline for New Patient

**Goal:** Compute all FC matrices and LRG analysis from raw data.

**Steps:**
1. Inspect data: `inspect_patient(patient, dataset_root)`
2. Compute correlation: `compute_corr_for_patient(patient, ...)`
3. Compute MSC: `compute_msc_for_patient(patient, ...)`
4. Compute LRG for both methods:
   - `compute_lrg_for_patient(patient, fc_method="corr", ...)`
   - `compute_lrg_for_patient(patient, fc_method="msc", ...)`
5. Generate standard figures

**Caveats:**
- Pat_06, Pat_07: Skip cross-phase analysis (missing phases)
- Default sample rate: 2048 Hz
- Default surrogates: 200

---

## Task: Compare FC Methods

**Goal:** Compare correlation vs MSC functional connectivity.

**Steps:**
1. Load LRG results for both methods
2. Use `compare_fc_methods(patient, phase, band, lrg_cache_root)`
3. Generate comparison plots

**Output:** UltrametricComparison dataclass with multiple distance metrics

---

## Task: Phase Reorganization Analysis

**Goal:** Analyze how brain networks change across experimental phases.

**Steps:**
1. Load LRG results for all 4 phases
2. Use `plot_phase_reorganization(patient, band, fc_method, ...)`
3. Compute distance matrices with `compute_phase_distance_matrix(...)`

**Phases:** rest_pre → task_learn → task_test → rest_post

---

## Task: 3D Brain Visualization

**Goal:** Create interactive/static 3D brain network plots.

**Requirements:**
- Electrode coordinates in `Implant_pat_XX.csv`
- LRG results (for cluster coloring)
- FC matrix (for edge weights)

**Functions:**
- Interactive: `plot_spatial_network_3d()` (Plotly)
- Static: `plot_spatial_network_3d_mpl()` (Matplotlib)
- Glass brain: `view_brain_connectome()` (nilearn)

---

## Task: Cache Status Check

**Goal:** Report what's computed vs missing.

**Check locations:**
```
data/corr_cache/Pat_XX/     # Correlation matrices
data/msc_cache/Pat_XX/      # MSC matrices
data/lrg_cache/Pat_XX/      # LRG results
```

**Expected files per patient:**
- 4 phases × 6 bands = 24 files per FC method
- 24 corr + 24 msc + 48 lrg = 96 cache files total

---

## Task: Validate Notebooks

**Goal:** Ensure notebooks run without errors.

**Steps:**
1. List notebooks: `ipynb/**/*.ipynb`
2. Execute with: `jupyter nbconvert --execute --inplace`
3. Report failures with traceback

**Notebook categories:**
- `00_intake/` - Data loading
- `01_preprocessing/` - Signal processing
- `02_fc_msc/` - Functional connectivity
- `03_lrg/` - LRG analysis
- `04_reorganization/` - Phase comparisons
- `05_figures/` - Publication figures

---

## Task: Add New Visualization Function

**Goal:** Create a new plot type following project patterns.

**Requirements:**
1. Place in appropriate `visuals/*.py` module
2. Follow signature pattern:
   ```python
   def plot_X(
       patient: str,
       phase: str,
       band: str,
       cache_root: Path = Path("data/X_cache"),
       output_path: Optional[Path] = None,
       figsize: tuple = (8, 8),
   ) -> Path:
   ```
3. Use cache loaders, never recompute
4. Always close figure: `plt.close(fig)`
5. Create output dirs: `output_path.parent.mkdir(parents=True, exist_ok=True)`
6. Export in `visuals/__init__.py`

---

## Task: Debug Missing Data

**Goal:** Diagnose why analysis fails for a patient/phase/band.

**Check order:**
1. Raw data exists: `data/stereoeeg_patients/Pat_XX/{phase}.mat`
2. MAT file readable: `inspect_mat_file(path)`
3. Cache exists: `get_X_cache_path(...)`
4. Cache loadable: `load_X_matrix(...)`
5. Parameters match cache naming

**Common issues:**
- Wrong `sparsify` parameter for MSC
- Missing phases for Pat_06/Pat_07
- Corrupted .mat files (try h5py fallback)
