---
name: readme
type: guide
era: CROSS_ERA
status: current
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

# Agent Guide: LRG EEG Functional Connectivity

This is the canonical guidance for coding agents working in this repository.
Root `AGENTS.md` and `CLAUDE.md` should remain identical and point here.
Start with `.agents/START_HERE.md`, then use `.agents/plans/INDEX.md`.

> **Before citing any number or re-running any script, read
> `.agents/reports/PIPELINE_STATUS.md` — it is the era index that tags
> every artefact as current (`imcoh_abs` post-reset), superseded
> (|ImCoh|² or MSC era), or dead-end. The current H1-H4 VI(k) results
> live in `.agents/reports/H1_H4_VI_RESULTS_POST_RESET.md`.**

---

## NOTICE: Data folder has been reorganized

The `data/` folder now has **only 4 top-level categories**:
- `data/raw/` — raw SEEG patient data
- `data/cache/` — all computation caches (corr, msc, lrg, imcoh, imcoh_lrg, etc.)
- `data/reports/` — investigation outputs grouped by topic
- `data/outputs/` — CLI figures (`outputs/figures/`) and tables (`outputs/tables/`)

The old flat paths (`data/corr_cache/`, `data/msc_cache/`, `data/figures/`, etc.)
**no longer exist** — they were moved into the new structure. Any code still
using `Path("data/corr_cache")` will break. Always use `config.paths` constants.

**For new code, always use `config.paths` constants:**
```python
from lrg_eegfc.config.paths import (
    SEEG_DATAPATH, CORR_CACHE, MSC_CACHE, LRG_CACHE,
    IMCOH_CACHE, IMCOH_LRG_CACHE, FIGURES_ROOT, TABLES_ROOT,
)
```

**For scripts, shared helpers are available:**
```python
from lrg_eegfc.utils.scripting import setup_script_env, iter_patient_band_phase, save_figure
ROOT = setup_script_env()
```

---

## Purpose and scope

This repo is a Python toolkit for building frequency-specific functional
connectivity (FC) networks from stereo-EEG (SEEG) recordings and analyzing
those networks with Laplacian Renormalization Group (LRG) methods:
- Convert SEEG time series into FC matrices (correlation, MSC, **ImCoh**)
- Select thresholds and clean matrices for network construction
- Compute LRG ultrametric distances, dendrograms, and entropy curves
- Compare FC methods and phase changes across patients
- Generate publication-grade visualizations

**ImCoh (Imaginary Coherence)** is the preferred FC method for community-level
analysis — it eliminates volume conduction bias that contaminates MSC.
See `.agents/guides/IMCOH_GUIDE.md`.

## Source of truth

- `README.md` for user-facing setup and CLI usage
- `src/lrg_eegfc/` for the current code structure and APIs
- `src/lrg_eegfc/config/paths.py` for **all data path constants**
- `scripts/` plus `scripts/py/*.py` for batch analysis pipelines

## Repository layout

- `src/lrg_eegfc/` — Python package
- `scripts/py/*.py` — pipeline scripts (being reorganized into numbered dirs)
- `scripts/` — shell wrappers for full analyses
- `data/` — local datasets and outputs (not committed, see layout above)
- `lrgsglib/` — submodule dependency (required for LRG operations)
- `.agents/guides/` — reference guides for agents
- `.agents/plans/` — active, developed, and archived plans

## Package map (current)

```
src/lrg_eegfc/
|-- __init__.py
|-- notebook.py
|-- cli/                  # unified CLI (lrg-eegfc command)
|   |-- __init__.py       #   exports `app`
|   |-- _app.py           #   LazyGroup root + 7 lazy subcommands
|   |-- _common.py        #   shared options, resolvers, CliReporter
|   |-- compute.py        #   6 compute commands
|   |-- plot.py           #   16 plot commands
|   |-- show.py           #   4 show commands (query cached results)
|   |-- data.py           #   3 data inspection commands
|   |-- cache.py          #   4 cache management commands
|   |-- config_cmd.py     #   2 config display commands
|   |-- bundle.py         #   1 publication bundle command
|   `-- _legacy.py        #   old cli.py (lrg-eegfc-corr compat)
|-- workflow/             # canonical workflows
|   |-- corr.py           #   correlation FC
|   |-- msc.py            #   MSC / ImCoh FC
|   |-- lrg.py            #   LRG analysis (accepts fc_method: corr/msc/imcoh)
|   |-- cleaning.py       #   Marchenko-Pastur cleaning
|   |-- diagnostics.py    #   LRG diagnostics
|   `-- time_windows.py   #   sliding-window FC
|-- config/
|   |-- const.py          #   constants (bands, phases, patients, etc.)
|   |-- paths.py          #   **centralized data paths** (all cache/output dirs)
|   `-- plotlib.py        #   matplotlib config
|-- utils/
|   |-- common.py         # shared imports
|   |-- scripting.py      # **shared script helpers** (setup_script_env, save_figure, etc.)
|   |-- io/               # MAT loading and patient metadata
|   |-- fc/               # FC pipelines (corr + msc/imcoh)
|   |-- lrg/              # clustering helpers
|   |-- metrics/          # comparison metrics
|   `-- pipelines/        # batch orchestration
|-- visuals/
|   |-- correlation.py
|   |-- msc.py
|   |-- compare.py
|   |-- lrg.py
|   |-- spatial.py        # 3D brain visualization
|   |-- reorganization.py
|   |-- cross_condition.py
|   `-- metastable.py
```

## Core workflows

1) **Correlation FC** — `workflow.corr.compute_corr_matrix`
2) **MSC FC** — `workflow.msc.compute_msc_matrix`
3) **ImCoh FC** — `workflow.msc.compute_msc_matrix` with `metric="imcoh"` passed through pipeline
4) **Correlation cleaning** — `workflow.cleaning.clean_correlation_matrix_full`
5) **LRG analysis** — `workflow.lrg.compute_lrg_analysis` (fc_method: `"corr"`, `"msc"`, `"imcoh"`)
6) **Comparisons** — `utils.metrics.compare`, `visuals/reorganization.py`

## Data layout

```
data/
├── raw/stereoeeg_patients/Pat_XX/   # .mat files + electrode metadata
├── cache/                           # all computation caches
│   ├── corr/, msc/, lrg/           # standard FC caches
│   ├── imcoh/, imcoh_lrg/          # ImCoh caches
│   └── ...                         # dev, windows, crema, etc.
├── reports/                         # investigation outputs by topic
│   ├── metric_exploration/
│   ├── imcoh/, imcoh_unanimity/, imcoh_vi/
│   └── ...
└── outputs/
    ├── figures/                     # CLI figure output
    └── tables/                      # CSV, tex, npz tables
```

## Agent quickstart checklist

- Use `config.paths` for all data paths — never hardcode `Path("data/...")`
- Use `config.const` for bands, phases, patients — never redeclare locally
- Use `utils.scripting` helpers in standalone scripts
- Confirm dataset root exists: `SEEG_DATAPATH` (resolves via symlink)
- For ImCoh work, use `IMCOH_CACHE` and `IMCOH_LRG_CACHE` from `config.paths`
- For visual outputs, use `visuals/` modules or `lrg-eegfc plot` CLI

## Quick mental model

SEEG data -> bandpass filter -> FC matrix (corr, MSC, or **ImCoh**) ->
threshold/cleaning -> graph + LRG ultrametrics/entropy -> comparisons +
visualizations.
