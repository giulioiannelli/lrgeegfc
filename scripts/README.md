# Pipeline Scripts Guide

## Directory Structure

Scripts are organized into numbered directories by function:

```
scripts/
├── 01_compute/          # FC computation, LRG, ImCoh, bipolar (16 scripts)
├── 02_visualize/        # Per-patient visualization and comparison (19 scripts)
├── 03_analysis/         # Quantitative structural analysis (18 scripts)
├── 04_reorganization/   # Phase reorganization and hypothesis testing (13 scripts)
├── 05_multiscale/       # Multiscale and cross-frequency analysis (11 scripts)
├── 06_metric_sweep/     # Systematic metric exploration (10 scripts)
├── 07_figures/          # Publication figure generation (30 scripts)
├── 08_epileptic/        # Epileptic node analysis (6 scripts)
├── 09_surrogate/        # Surrogate validation and probe bias (9 scripts)
├── wp0/                 # Work package 0: metric exploration framework
├── wp1/                 # Work package 1: clinical application
├── config.sh            # Default pipeline parameters
├── run_step.sh          # Modular step runner
├── run_full_analysis.sh # Full pipeline runner
└── migrate_data.sh      # Data folder reorganization tool
```

## Key Scripts by Category

### 01_compute/ — Data production
Core computation scripts that produce cached results:
- `compute_corr_matrices.py` — Correlation FC matrices
- `compute_msc_matrices.py` — MSC FC matrices
- `compute_lrg_analysis.py` — LRG analysis from FC
- `compute_imcoh_lrg.py` — LRG from ImCoh FC (volume-conduction immune)
- `compute_imcoh_vi.py` — Multiscale VI profiles for ImCoh
- `compute_imcoh_unanimity.py` — Cross-patient unanimity maps

### 02_visualize/ — Visualization
Generate figures from cached results (never recompute):
- `visualize_lrg.py` — LRG dendrograms, entropy, full panels
- `visualize_imcoh_gallery.py` — ImCoh six-panel gallery
- `compare_msc_imcoh.py` — Side-by-side MSC vs ImCoh comparison

### 07_figures/ — Publication figures
Generate final publication-quality figures:
- `gen_*` scripts produce specific figure sets
- `consolidate_report_figures.py` — Merge figures for reports
- `build_overleaf_bundle.py` — Package for LaTeX

## Usage

### Data paths
All scripts should use `config.paths` constants:
```python
from lrg_eegfc.config.paths import MSC_CACHE, LRG_CACHE, IMCOH_CACHE, FIGURES_ROOT
```

### Script helpers
Shared boilerplate is available from `utils.scripting`:
```python
from lrg_eegfc.utils.scripting import setup_script_env, iter_patient_band_phase, save_figure
ROOT = setup_script_env()
```

### Shell wrappers
```bash
# List available steps
bash scripts/run_step.sh --list

# Test a step on one case
bash scripts/run_step.sh --test --step 5b

# Run full pipeline
bash scripts/run_full_analysis.sh
```

## Work Packages

- `wp0/` — Metric exploration framework (7 scripts with shared `_common.py`)
- `wp1/` — Clinical application analysis (7 scripts)

These are self-contained investigation suites with their own shared utilities.
