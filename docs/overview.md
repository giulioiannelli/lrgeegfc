# Project overview

This document summarises the architecture of the **LRG EEG Functional
Connectivity** toolkit after the 2024 refactor.

## Package layout (current)

```
src/lrg_eegfc/
|-- __init__.py              # Public API exports
|-- cli.py                   # CLI for correlation workflow
|-- compare.py               # Compatibility wrapper
|-- plotting.py              # Compatibility wrapper
|-- io.py                    # Compatibility wrapper
|-- workflow/                # Canonical workflows
|-- workflow.py              # Compatibility wrapper
|-- workflow_corr.py         # Compatibility wrapper
|-- workflow_msc.py          # Compatibility wrapper
|-- workflow_lrg.py          # Compatibility wrapper
|-- workflow_cleaning.py     # Compatibility wrapper
|-- batch_compute.py         # Compatibility wrapper
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

The `lrg-eegfc-corr` entry point in :mod:`lrg_eegfc.cli` targets the
correlation-based workflow: load a patient/phase recording, band-pass filter
it, compute a correlation matrix, apply a percolation threshold, and optionally
produce plots. All filesystem paths are configurable and the command reuses
cached correlation matrices unless `--overwrite` is passed.

MSC and LRG workflows are accessed through the Python APIs or the
`src/*.py` pipeline scripts invoked by `scripts/run_step.sh`.

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
