# Project overview

This document summarises the architecture of the **LRG EEG Functional
Connectivity** toolkit. Post-2026-04-15 the canonical FC metric is
`imcoh_abs` (= `⟨|Im(S_ij)/√(S_ii S_jj)|⟩_f`, Nolte 2004 / Ewald 2012);
MSC is retained as a methodological baseline only. See
`.agents/guides/02_methods/IMCOH_GUIDE.md` for rationale and
`.agents/reports/PIPELINE_STATUS.md` for the era index of artefacts.

## Package layout (current)

```
src/lrg_eegfc/
|-- __init__.py              # Public API exports
|-- cli/                     # Unified CLI (lazy-loaded subcommands)
|-- workflow/                # Canonical workflows
|-- config/                  # Dataset-aware constants
|-- utils/                   # Loaders + FC primitives
`-- visuals/                 # Figure generation
```

### Data flow

1. **Loading** – `utils/io/` handles `.mat` files and channel metadata.
   The robust loader in `patient_robust.py` supports multiple variable names.
2. **Processing (correlation)** – `utils/fc/corr/` builds correlation matrices,
   applies percolation-based thresholds, and supports MP cleaning.
3. **Processing (MSC)** – `utils/fc/msc/` computes magnitude-squared
   coherence with optional surrogate-based sparsification.
4. **Workflows** – `workflow/corr.py`, `workflow/msc.py`, and
   `workflow/cleaning.py` add caching and orchestration.
5. **LRG analysis** – `workflow/lrg.py` computes ultrametrics, dendrograms, and
   entropy curves using `lrgsglib`.
6. **Presentation** – `visuals/plotting.py` supports the CLI; `visuals/` contains
   notebook-ready figures and multi-panel outputs.

## Command line interface

The unified `lrg-eegfc` CLI exposes subcommand groups (`compute`, `plot`,
`show`, `data`, `cache`, `config`, `bundle`) for all workflows. It is lazy
loaded to keep `--help` fast and avoid importing heavy scientific dependencies
until needed.

The legacy `lrg-eegfc-corr` entry point remains for correlation-only
workflows: load a patient/phase recording, band-pass filter it, compute a
correlation matrix, apply a percolation threshold, and optionally produce
plots. All filesystem paths are configurable and the command reuses cached
correlation matrices unless `--overwrite` is passed.

MSC and LRG workflows are accessed through the unified CLI, the Python APIs,
or the pipeline scripts invoked by `scripts/run_step.sh`.

## Dependency notes

* `lrgsglib` supplies the Laplacian-based operators (entropy and percolation
  utilities) and must be installed separately.
* `networkx`, `numpy`, `pandas`, `matplotlib`, `scipy`, `h5py` and `emd` are
  standard scientific Python dependencies already declared in ``pyproject.toml``.

## Extending the toolkit

* New features should be added as separate modules under ``lrg_eegfc`` and
  exported through ``__init__.__all__``.
* Prefer pure functions with explicit parameters; avoid hidden state and module
  level singletons.
* Keep CLI behaviour deterministic – anything random should accept a seed
  parameter with a sensible default.
