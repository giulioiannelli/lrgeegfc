---
name: checkpoint-visualization-implementation
type: plan
era: MSC
status: superseded
created: 2025-12-10
updated: 2026-04-24
pointers: []
---

# Visualization Implementation Checkpoint

**Date:** 2025-12-10
**Status:** Phases 1-4 Complete, Phase 5 Pending

---

## Summary

We have successfully implemented a comprehensive visualization pipeline for LRG EEG Functional Connectivity analysis. The implementation covers correlation-based FC, MSC-based FC, and LRG hierarchical analysis, with complete CLI tools for batch processing.

---

## ✅ COMPLETED WORK

### Phase 1: Correlation Cleaning Workflow
**Status:** ✅ Complete and tested

**Files created:**
- `src/lrg_eegfc/workflow_cleaning.py` - Marchenko-Pastur noise removal
- `src/clean_correlation_matrices.py` - CLI for batch cleaning

**What it does:**
- Applies spectral cleaning to correlation matrices using Marchenko-Pastur law
- Removes random matrix noise from temporal correlations
- Saves cleaned matrices with metadata (gamma, lambda_min, lambda_max)

**Test command:**
```bash
python src/clean_correlation_matrices.py --patient Pat_02 --verbose
```

**Output:**
- `data/corr_cache/Pat_XX/{band}_{phase}_corr_cleaned.npy` - Cleaned matrices
- `data/corr_cache/Pat_XX/{band}_{phase}_corr_cleaned_meta.npz` - Metadata

---

### Phase 2: Correlation Visualization Suite
**Status:** ✅ Complete and tested

**Files created:**
- `src/lrg_eegfc/visuals/correlation.py` - 4 visualization functions
  - `plot_correlation_heatmap()` - Matrix heatmap
  - `plot_correlation_and_network()` - Side-by-side matrix + network
  - `plot_marchenko_pastur_comparison()` - 3-panel: original | cleaned | eigenvalue comparison
  - `plot_percolation_curves()` - Threshold selection curves (P_inf, E_inf)
- `src/visualize_correlation.py` - CLI script with batch mode

**What it does:**
- Visualizes correlation matrices (cleaned and uncleaned)
- Shows network graphs with edge weights
- Compares original vs MP-cleaned eigenvalue distributions
- Displays percolation curves for threshold selection

**Test commands:**
```bash
# Single visualization (all plot types)
python src/visualize_correlation.py --patient Pat_02 --phase rsPre --band beta --plot-type all --verbose

# Batch mode (all bands for a phase)
python src/visualize_correlation.py --patient Pat_02 --phase rsPre --batch --verbose

# Cleaned matrices
python src/visualize_correlation.py --patient Pat_02 --phase rsPre --batch --cleaned --verbose
```

**Output:**
- `data/figures/correlation/Pat_XX/{band}_{phase}_matrix.png`
- `data/figures/correlation/Pat_XX/{band}_{phase}_network.png`
- `data/figures/correlation/Pat_XX/{band}_{phase}_mp_comparison.png`
- `data/figures/correlation/Pat_XX/{band}_{phase}_percolation.png`

**Verified working:** ✅ Tested with Pat_02 beta rsPre

---

### Phase 3: MSC Visualization Suite
**Status:** ✅ Complete and tested

**Files created:**
- `src/lrg_eegfc/visuals/msc.py` - 3 visualization functions
  - `plot_msc_heatmap()` - MSC matrix heatmap (range [0,1])
  - `plot_msc_and_network()` - Side-by-side matrix + network
  - `plot_msc_comparison_dense_vs_validated()` - Compare dense vs surrogate-validated MSC
- `src/visualize_msc.py` - CLI script with batch mode

**What it does:**
- Visualizes MSC matrices (no cleaning needed - already frequency-specific)
- Shows dense MSC (all connections) vs validated MSC (surrogate-based)
- Network graphs with edge weights scaled to [0.1, 3.0]

**Test commands:**
```bash
# Dense MSC network
python src/visualize_msc.py --patient Pat_02 --phase rsPre --band beta --plot-type network --verbose

# Validated MSC (with 200 surrogates)
python src/visualize_msc.py --patient Pat_02 --phase rsPre --band beta --plot-type network --sparsify soft --n-surrogates 200 --verbose

# Comparison plot
python src/visualize_msc.py --patient Pat_02 --phase rsPre --band beta --plot-type comparison --n-surrogates 200 --verbose

# Batch mode
python src/visualize_msc.py --patient Pat_02 --phase rsPre --batch --verbose
```

**Output:**
- `data/figures/msc/Pat_XX/{band}_{phase}_msc_dense_network.png`
- `data/figures/msc/Pat_XX/{band}_{phase}_msc_validated_nsurr200_network.png`
- `data/figures/msc/Pat_XX/{band}_{phase}_msc_comparison_nsurr200.png`

**Verified working:** ✅ Tested with Pat_02 rsPre (all bands)

---

### Phase 4: LRG Visualization Suite
**Status:** ✅ Complete (awaiting LRG cache data for testing)

**Files created:**
- `src/lrg_eegfc/visuals/lrg.py` - 4 visualization functions
  - `plot_lrg_entropy_curves()` - 1-S and C vs tau (scale parameter)
  - `plot_lrg_dendrogram()` - Hierarchical dendrogram with optimal threshold
  - `plot_ultrametric_heatmap()` - Ultrametric distance matrix
  - `plot_lrg_full_panel()` - 4-panel comprehensive view (2x2 grid)
- `src/visualize_lrg.py` - CLI script with batch mode

**What it does:**
- Visualizes LRG hierarchical analysis results
- Shows entropy evolution across scales
- Displays hierarchical clustering with ultrametric distances
- Comprehensive 4-panel view: entropy | dendrogram | ultrametric | cluster-colored network

**Test commands:**
```bash
# Entropy curves (correlation-based FC)
python src/visualize_lrg.py --patient Pat_02 --phase rsPre --band beta --fc-method corr --plot-type entropy --verbose

# Full 4-panel visualization
python src/visualize_lrg.py --patient Pat_02 --phase rsPre --band beta --fc-method corr --plot-type full --verbose

# All plot types
python src/visualize_lrg.py --patient Pat_02 --phase rsPre --band beta --fc-method corr --plot-type all --verbose

# Batch mode (all bands)
python src/visualize_lrg.py --patient Pat_02 --phase rsPre --fc-method corr --batch --verbose

# MSC-based FC
python src/visualize_lrg.py --patient Pat_02 --phase rsPre --band beta --fc-method msc --plot-type full --verbose
```

**Output:**
- `data/figures/lrg/Pat_XX/{band}_{phase}_lrg_{fc_method}_entropy.png`
- `data/figures/lrg/Pat_XX/{band}_{phase}_lrg_{fc_method}_dendrogram.png`
- `data/figures/lrg/Pat_XX/{band}_{phase}_lrg_{fc_method}_ultrametric.png`
- `data/figures/lrg/Pat_XX/{band}_{phase}_lrg_{fc_method}_full.png`

**Testing status:** Ready to test once LRG cache is generated by `src/compute_lrg_analysis.py`

---

### Complete Analysis Pipeline Script
**Status:** ✅ Complete

**File created:**
- `scripts/run_full_analysis.sh` - End-to-end pipeline for all patients

**What it does:**
Runs the complete analysis pipeline for all available patients:
1. Computes correlation matrices
2. Cleans correlation matrices (Marchenko-Pastur)
3. Generates correlation visualizations (uncleaned)
4. Generates correlation visualizations (cleaned)
5. Computes MSC matrices (if script exists)
6. Compares FC methods (if script exists)
7. Computes LRG analysis (if script exists)

**Execute:**
```bash
bash scripts/run_full_analysis.sh
```

**Features:**
- Color-coded output (green=success, yellow=warning, red=error)
- Timing information per patient
- Graceful handling of missing scripts (phases 5-7)
- Continues if one patient fails

**Current configuration:** Processes Pat_02 and Pat_03 (can be edited to add more)

---

## 🚧 PENDING WORK

### Phase 5: Comparison Visualization Suite
**Status:** ❌ Not started

**What needs to be implemented:**

#### File 1: `src/lrg_eegfc/visuals/comparison.py`

**Functions needed:**
1. `compute_cluster_membership_correlation()` - Correlation of cluster assignments
2. `plot_msc_vs_corr_comparison()` - 2x3 grid comparing MSC vs Correlation LRG results
3. `plot_phase_reorganization()` - Network reorganization across experimental phases

**User requirements addressed:**
- "correlation between cluster belonging of each node"
- "measure of distance of ultrametric matrices" (uses existing `compare.py`)
- "topological measure of similarity hierarchical tree to hierarchical tree"
- Phase reorganization analysis (rsPre → task → rsPost)

#### File 2: `src/visualize_comparison.py`

**CLI script for:**
- Comparing MSC vs Correlation methods
- Tracking reorganization across phases
- Generating comparison reports

**Specifications available in:** `VISUALIZATION_PLAN_DETAILED.md` (lines 923-1542)

---

### Phase 6: Documentation
**Status:** ❌ Not started

**Tasks:**
- Create `VISUALS.md` documenting all visualization functions
- Update README.md with visualization examples
- Document interpretation of each plot type

---

### Phase 7: Testing & Validation
**Status:** ⚠️ Partially complete

**Completed:**
- ✅ Phase 2 tested with Pat_02 beta rsPre
- ✅ Phase 3 tested with Pat_02 rsPre (all bands)

**Pending:**
- ❌ Phase 4 testing (need to generate LRG cache first)
- ❌ Phase 5 testing
- ❌ Full pipeline end-to-end test on all patients

---

## 📁 File Organization

```
src/
├── lrg_eegfc/
│   ├── visuals/
│   │   ├── __init__.py          ✅ Updated with all exports
│   │   ├── correlation.py       ✅ Phase 2
│   │   ├── msc.py              ✅ Phase 3
│   │   ├── lrg.py              ✅ Phase 4
│   │   └── comparison.py        ❌ Phase 5 (pending)
│   ├── workflow_cleaning.py     ✅ Phase 1
│   ├── workflow_corr.py         ✅ Pre-existing
│   ├── workflow_msc.py          ✅ Pre-existing
│   ├── workflow_lrg.py          ✅ Pre-existing
│   └── compare.py               ✅ Pre-existing (used by Phase 5)
├── visualize_correlation.py     ✅ Phase 2 CLI
├── visualize_msc.py            ✅ Phase 3 CLI
├── visualize_lrg.py            ✅ Phase 4 CLI
├── visualize_comparison.py      ❌ Phase 5 CLI (pending)
├── clean_correlation_matrices.py ✅ Phase 1 CLI
├── compute_corr_matrices.py    ✅ Pre-existing
├── compute_msc_matrices.py     ✅ Pre-existing
├── compute_lrg_analysis.py     ✅ Pre-existing
└── compare_fc_methods.py       ✅ Pre-existing

scripts/
└── run_full_analysis.sh         ✅ Complete pipeline

data/
├── corr_cache/Pat_XX/           ✅ Correlation matrices (cleaned & uncleaned)
├── msc_cache/Pat_XX/            ✅ MSC matrices
├── lrg_cache/Pat_XX/            ⚠️ Will be created by pipeline
└── figures/
    ├── correlation/Pat_XX/      ✅ Generated
    ├── msc/Pat_XX/              ✅ Generated
    ├── lrg/Pat_XX/              ⚠️ Will be generated
    └── comparison/Pat_XX/       ❌ Pending Phase 5
```

---

## 🔧 Available Patients

**Full dataset (4 phases):**
- Pat_02 ✅
- Pat_03 ✅
- Pat_05 ✅
- Pat_07 ✅
- Pat_08 ✅

**Partial dataset (2 phases: rsPre, rsPost):**
- Pat_06 ⚠️

**Phases:** `rsPre`, `taskLearn`, `taskTest`, `rsPost`

**Bands:** `delta`, `theta`, `alpha`, `beta`, `low_gamma`, `high_gamma`

---

## 🚀 Quick Start Commands

### Test Individual Components

```bash
# Phase 2: Correlation visualizations
python src/visualize_correlation.py --patient Pat_02 --phase rsPre --band beta --plot-type all --verbose

# Phase 3: MSC visualizations
python src/visualize_msc.py --patient Pat_02 --phase rsPre --band beta --plot-type network --verbose

# Phase 4: LRG visualizations (after generating LRG cache)
python src/visualize_lrg.py --patient Pat_02 --phase rsPre --band beta --fc-method corr --plot-type full --verbose
```

### Run Full Pipeline

```bash
# Complete analysis for all patients
bash scripts/run_full_analysis.sh

# Just correlation workflow (fastest)
for PATIENT in Pat_02 Pat_03; do
  python src/compute_corr_matrices.py --patient $PATIENT --verbose && \
  python src/clean_correlation_matrices.py --patient $PATIENT --verbose && \
  python src/visualize_correlation.py --patient $PATIENT --batch --verbose
done
```

---

## 🐛 Known Issues & Fixes

### Issue 1: Missing metadata in MP comparison
**Status:** ✅ Fixed
**Fix:** Added `metadata = None` initialization in exception handler (correlation.py:295)

### Issue 2: Type checks for None values
**Status:** ✅ Fixed
**Fix:** Added `metadata is not None` checks before dictionary access (correlation.py:313, 343)

### Issue 3: LRG cache not yet generated
**Status:** ⚠️ Expected
**Resolution:** Will be created when `compute_lrg_analysis.py` runs in full pipeline

---

## 📊 Expected Outputs After Full Pipeline

After running `bash scripts/run_full_analysis.sh`, you should have:

```
data/figures/
├── correlation/
│   ├── Pat_02/
│   │   ├── {band}_{phase}_matrix.png                  (24 files: 6 bands × 4 phases)
│   │   ├── {band}_{phase}_network.png                 (24 files)
│   │   ├── {band}_{phase}_mp_comparison.png           (24 files)
│   │   └── {band}_{phase}_percolation.png             (24 files)
│   ├── Pat_03/ ... (same structure)
│   └── ...
├── msc/
│   ├── Pat_02/
│   │   ├── {band}_{phase}_msc_dense_network.png       (24 files)
│   │   └── ... (if validated: _validated_nsurr200_network.png)
│   └── ...
└── lrg/
    ├── Pat_02/
    │   ├── {band}_{phase}_lrg_corr_entropy.png        (24 files)
    │   ├── {band}_{phase}_lrg_corr_dendrogram.png     (24 files)
    │   ├── {band}_{phase}_lrg_corr_ultrametric.png    (24 files)
    │   ├── {band}_{phase}_lrg_corr_full.png           (24 files)
    │   └── ... (same for _lrg_msc_* if MSC-based LRG computed)
    └── ...
```

**Total per patient:** ~96-192 figures (depending on MSC LRG inclusion)

---

## 📝 Next Session Tasks

### Priority 1: Complete Phase 5
1. Implement `src/lrg_eegfc/visuals/comparison.py`
   - `compute_cluster_membership_correlation()`
   - `plot_msc_vs_corr_comparison()`
   - `plot_phase_reorganization()`
2. Implement `src/visualize_comparison.py` CLI
3. Test comparison visualizations

### Priority 2: Run Full Pipeline
1. Execute `bash scripts/run_full_analysis.sh`
2. Verify all outputs generated correctly
3. Check for any errors or missing data

### Priority 3: Documentation
1. Create `VISUALS.md` with:
   - Function documentation
   - Interpretation guides
   - Example gallery
2. Update README.md with visualization section

### Priority 4: Code Organization
1. Clean up markdown files in root directory (per ToDo.md)
2. Move documentation files to `docs/`
3. Archive session notes

---

## 🔗 Key Reference Files

- `VISUALIZATION_PLAN.md` - Original user requirements and high-level plan
- `VISUALIZATION_PLAN_DETAILED.md` - Detailed implementation specs for Phases 3-7
- `NEXT_SESSION_START_HERE.md` - Previous checkpoint (now superseded by this file)
- `CLAUDE.md` - Project overview and coding standards
- `ANALYSIS_COMMANDS.md` - Command reference

---

## ⚙️ Technical Notes

### Visualization Design Principles
1. **Separation of concerns:** Visualization reads from cache, never recomputes
2. **Consistent naming:** `{band}_{phase}_{analysis_type}_{variant}_{plot_type}.png`
3. **Batch processing:** All CLIs support `--batch` mode for all bands
4. **Parameter preservation:** Cache filenames encode all computation parameters
5. **Graceful degradation:** Missing channel labels → use indices

### Color Conventions
- **Correlation matrices:** `viridis` (sequential, [min, max])
- **MSC matrices:** `viridis` (sequential, [0, 1])
- **Eigenvalue comparison:** `viridis` + `orange` (MP curve)
- **Percolation curves:** Black hexagons (P_inf) + Blue line (E_inf)
- **Dendrograms:** Automatic coloring by cluster, blue dashed line = optimal threshold
- **Network graphs:** Lightblue nodes (default), cluster-colored (LRG full panel)

### Edge Width Scaling
- Correlation networks: [0.05, 0.35]
- MSC networks: [0.1, 3.0]
- LRG networks: [0.1, 5.0]

---

## 💾 Checkpoint Commands

### Save current state
```bash
# Already saved in this file!
# Location: CHECKPOINT_VISUALIZATION_IMPLEMENTATION.md
```

### Resume from checkpoint
```bash
# Read this file
cat CHECKPOINT_VISUALIZATION_IMPLEMENTATION.md

# Continue with Phase 5 implementation
# OR run full pipeline to generate all visualizations
bash scripts/run_full_analysis.sh
```

---

**Checkpoint created:** 2025-12-10
**Ready to resume:** Phase 5 implementation OR full pipeline execution
**Current status:** Phases 1-4 complete and tested (2-3) or ready for testing (4)
