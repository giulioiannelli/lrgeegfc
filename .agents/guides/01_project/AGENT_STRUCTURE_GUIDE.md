---
name: agent-structure-guide
type: guide
era: CROSS_ERA
status: current
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

# Agent Guide: Structure & Invariants

## Quick Orientation
- Core package lives in `src/lrg_eegfc/`; runnable scripts live in `scripts/py/` (e.g., `compute_corr_matrices.py`, `compute_msc_matrices.py`, `compare_fc_methods.py`).
- Session context and historical decisions are tracked in `.agents/plans/active/` (start with `.agents/START_HERE.md`, then use `.agents/plans/INDEX.md` for the full list). Completed plans live under `.agents/plans/developed/`.
- Data layout is `data/stereoeeg_patients/Pat_XX/{phase}.mat`; caches are kept alongside analysis outputs under `data/`.
- Companion library `lrgsglib` is a git submodule and must be installed for entropy/percolation utilities.

## Repository Layout (what lives where)
- `src/lrg_eegfc/` – Python package:
  - `utils/` – low-level primitives (`io/` loaders, `fc/corr/` correlation building/thresholds, `fc/msc/` MSC + surrogates, `lrg/` clustering helpers, `metrics/` comparisons, `pipelines/` orchestration).
  - `workflow/` – high-level, cached workflows that wrap the utils and return dataclasses.
  - `visuals/` – plotting modules (`correlation.py`, `msc.py`, `lrg.py`, `metastable.py`, `reorganization.py`, `compare.py`, `plotting.py`) used by CLI/nb.
  - `notebook.py` – curated imports + `setup_notebook()` for consistent figure dirs.
  - `config/const.py` – bands/phases and defaults.
- `scripts/py/*.py` – CLIs/batch scripts mirroring workflows; keep them read-from-cache (no recompute inside visualization scripts).
- `docs/` – high-level overview and developer guide (aligned to current module names).
- `ipynb/` – archived notebooks + refactored tutorials/analysis (see `.agents/plan/dev/PHASE2_PROGRESS.md`).

## Core Invariants to Preserve
- **Caching first:** Correlation (`data/corr_cache/{patient}/{band}_{phase}_corr_ftype-…npy`), cleaned correlation (`…_corr_cleaned.npy` + `…_meta.npz`), MSC (`data/msc_cache/{patient}/{band}_{phase}_msc_sparsify-{method}_…npy`), LRG (`data/lrg_cache/{patient}/{band}_{phase}_lrg_{method}.npz`). Visualizations must only read these caches.
- **Deterministic workflows:** All compute functions accept explicit params (fs, filter order, sparsify, jump index). Avoid hidden globals; prefer dataclass returns for metadata.
- **Robust loading:** Use `utils/io/patient_robust.py` when fs/shape may be missing; channel metadata is optional but should be passed through.
- **Band/phase canon:** Use `BRAIN_BANDS`, `PHASE_LABELS` from `config/const.py` (or re-exported by `lrg_eegfc.__init__`). Do not introduce a third source of truth.
- **Figure/output layout:** Keep figures under `data/figures/{type}/{patient}/` with `{band}_{phase}_{plot}.png`; keep comparison CSVs in `results/`.
- **No recomputation in viz:** CLI/visual modules should fail fast if caches are missing, not regenerate data.

## Primary APIs (preferred entry points)
- Correlation: `workflow.corr.compute_corr_matrix` / `compute_corr_for_patient` (caching), `utils.fc.corr.build_corr_network` for low-level.
- Cleaning: `workflow.cleaning.clean_correlation_matrix_full` (+ `load_cleaned_corr_matrix`) for MP + percolation thresholds.
- Coherence/MSC: `workflow.msc.compute_msc_matrix` / `compute_msc_for_patient`; low-level pipeline in `utils.fc.msc.coherence_fc_pipeline`.
- LRG/ultrametric: `workflow.lrg.compute_lrg_analysis` / `compute_lrg_for_patient`; comparisons in `utils.metrics.compare`.
- Visualization: `visuals/*.py` functions; CLI wrappers in `scripts/py/visualize_*.py`.
- Notebooks: `from lrg_eegfc.notebook import *` then `setup_notebook()` to set figure dir and imports.

## Status & Reference Files
- Data readiness: `.agents/plans/active/*DATA_STATUS*.md`, `.agents/plans/active/*PATIENT_DATA_REPORT*.txt`.
- Pipeline commands: `.agents/plans/active/*ANALYSIS_COMMANDS*.md` and `.agents/plans/archive/*PIPELINE_GUIDE*.md`.
- Notebook migration state: `.agents/plans/active/*PHASE2_PROGRESS*.md`.
- Validation/MSCs: `.agents/plans/active/*VALIDATION_PIPELINE*.md`.
