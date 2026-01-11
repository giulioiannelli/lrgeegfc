# Agent Guide: LRG EEG Functional Connectivity

This is the canonical guidance for coding agents working in this repository.
Root `AGENTS.md` and `CLAUDE.md` should remain identical and point here.
Start with `.agents/START_HERE.md`, then use `.agents/plans/INDEX.md`.

## Purpose and scope

This repo is a Python toolkit for building frequency-specific functional
connectivity (FC) networks from stereo-EEG (SEEG) recordings and analyzing
those networks with Laplacian Renormalization Group (LRG) methods. The suite
exists to:
- convert SEEG time series into correlation- or coherence-based FC matrices
- select thresholds and clean matrices for network construction
- compute LRG ultrametric distances, dendrograms, and entropy curves
- compare FC methods and phase changes across patients
- generate publication-grade visualizations

## Source of truth

Use these as the authoritative sources:
- `README.md` for user-facing setup and CLI usage
- `src/lrg_eegfc/` for the current code structure and APIs
- `scripts/` plus `src/*.py` for batch analysis pipelines

Note: `docs/overview.md` reflects the current module layout.

## Repository layout (relevant)

- `src/lrg_eegfc/` Python package
- `src/*.py` pipeline scripts invoked by `scripts/run_step.sh`
- `scripts/` shell wrappers for full analyses
- `data/` local datasets and outputs (not committed)
- `lrgsglib/` submodule dependency (required for LRG operations)
- `.agents/plans/active/` current plans and session notes
- `.agents/plans/developed/` completed plans (ready for reference)
- `.agents/plans/archive/` completed or superseded plans (see `.agents/plans/INDEX.md`)

## Package map (current)

```
src/lrg_eegfc/
|-- __init__.py
|-- cli.py
|-- notebook.py
|-- workflow/             # canonical workflows
|   |-- core.py
|   |-- corr.py
|   |-- msc.py
|   |-- lrg.py
|   `-- cleaning.py
|-- config/
|   `-- const.py
|-- utils/
|   |-- common.py         # shared imports
|   |-- io/               # MAT loading and patient metadata
|   |-- fc/               # FC pipelines (corr + msc)
|   |-- lrg/              # clustering helpers
|   |-- metrics/          # comparison metrics
|   `-- pipelines/        # batch orchestration
|-- visuals/
|   |-- correlation.py
|   |-- msc.py
|   |-- compare.py
|   |-- lrg.py
|   |-- lrg_revised.py    # wrapper
|   |-- lrg_backup.py     # wrapper
|   |-- reorganization.py
|   |-- plotting.py
|   `-- metastable.py
```

## Core workflows and what they produce

1) Correlation FC (quickest path)
- `workflow.core.compute_band_connectivity`
  - bandpass -> correlation -> percolation jumps -> thresholded graph
- `workflow.corr.compute_corr_matrix`
  - caches `data/corr_cache/Pat_XX/<band>_<phase>_corr_*.npy`

2) Coherence / MSC FC
- `workflow.msc.compute_msc_matrix`
  - magnitude-squared coherence, optional surrogate-based sparsification
  - caches `data/msc_cache/Pat_XX/<band>_<phase>_msc_*.npy`

3) Correlation cleaning (Marchenko-Pastur)
- `workflow.cleaning.clean_correlation_matrix_full`
  - bandpass -> correlation -> spectral cleaning -> threshold
  - caches cleaned matrices under `data/corr_cache/Pat_XX/`

4) LRG analysis
- `workflow.lrg.compute_lrg_analysis`
  - LRG entropy curves, ultrametric distances, dendrogram linkage
  - caches `data/lrg_cache/Pat_XX/<band>_<phase>_lrg_<fc>.npz`

5) Comparisons
- `utils.metrics.compare` compares ultrametric matrices across methods or phases
- `visuals/reorganization.py` summarizes phase reorganization and similarity

## Data layout (expected)

```
data/
`-- stereoeeg_patients/
    `-- Pat_XX/
        |-- rsPre.mat
        |-- taskLearn.mat
        |-- taskTest.mat
        |-- rsPost.mat
        |-- Implant_pat_XX.csv
        `-- channel_labels.csv
```

Notes:
- Standard loader expects a `Data` variable in `.mat` files.
- Robust loader (`utils/io/patient_robust.py`) scans common keys
  (`Data`, `data`, `EEG`, `timeseries`, `signal`, `eeg_data`).

## Output locations

- Correlation cache: `data/corr_cache/Pat_XX/`
- MSC cache: `data/msc_cache/Pat_XX/`
- LRG cache: `data/lrg_cache/Pat_XX/`
- Figures: `data/figures/<category>/Pat_XX/`
- CLI outputs: `data/correlations/Pat_XX/` (from `lrg-eegfc-corr`)

## Entry points and scripts

- CLI (correlation only):
  - `lrg-eegfc-corr --patient Pat_01 --phase rsPre --band beta --plot-all`

- Pipeline scripts (used by `scripts/run_step.sh`):
  - `python src/compute_corr_matrices.py`
  - `python src/compute_msc_matrices.py`
  - `python src/compute_lrg_analysis.py`
  - `python src/visualize_*.py`

- Orchestrators:
  - `bash scripts/run_step.sh --list`
  - `bash scripts/run_full_analysis.sh`

## Setup (local dev)

```
# submodule + venv
git submodule update --init --recursive
python -m venv .venv
source .venv/bin/activate

# install package + tooling
pip install -e .[dev]

# ensure lrgsglib is installed (submodule or git)
pip install ./lrgsglib
```

## Agent quickstart checklist

- Confirm dataset root exists (default `data/stereoeeg_patients/`) and that
  `.mat` files expose a usable data variable (`Data` or one of the robust keys).
- Decide the workflow: correlation (`workflow.corr`), MSC (`workflow.msc`), or
  LRG analysis (`workflow.lrg`) after an FC matrix exists.
- Use cache-aware helpers to avoid recomputation during iteration.
- Prefer `constants.py` for lightweight imports; `config/const.py` reads
  the filesystem at import time.
- For visual outputs, use `visuals/` modules or `scripts/run_step.sh`.

## Coding standards

- Python 3.11+ only
- Type hints required; mypy is configured in strict mode
- Format with black and isort (88-char line length)
- Explicit imports; avoid wildcard imports outside re-export modules

## Known mismatches and caveats

- `docs/overview.md` references modules that no longer exist. Use the
  `src/lrg_eegfc/` tree above as the accurate map.
- `config/const.py` reads `data/stereoeeg_patients` at import time to build a
  patient list; it will fail if the dataset is missing. Prefer `constants.py`
  when you need a lightweight import.
- The CLI is correlation-focused; MSC and LRG workflows are accessed via
  Python APIs or the pipeline scripts.

## Quick mental model

SEEG data -> bandpass filter -> FC matrix (corr or MSC) -> threshold/cleaning
-> graph + LRG ultrametrics/entropy -> comparisons + visualizations.
