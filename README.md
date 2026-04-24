# LRG EEG Functional Connectivity

Tools for computing frequency-specific functional connectivity graphs from
stereo-EEG (SEEG) recordings collected in the [Living and Relational Graphs
(LRG)](https://github.com/giulioiannelli) research project.  The repository now
ships as a regular Python package, complete with a command line interface,
documented public APIs and optional plotting helpers.

---

- [Quick Start](#quick-start)
- [Installation](#installation)
  - [Step 1: Clone the repository](#step-1-clone-the-repository)
  - [Step 2: Initialize submodules](#step-2-initialize-submodules)
  - [Step 3: Create the conda environment](#step-3-create-the-conda-environment)
  - [Step 4: Build and configure lrgsglib](#step-4-build-and-configure-lrgsglib)
  - [Step 5: Install lrg-eegfc](#step-5-install-lrg-eegfc)
- [Dataset layout](#dataset-layout)
- [Command line usage](#command-line-usage)
- [Python API overview](#python-api-overview)
- [Plotting utilities](#plotting-utilities)
- [Developer guide](#developer-guide)

---

## Quick Start

```bash
# Clone and setup
git clone https://github.com/giulioiannelli/lrgeegfc.git
cd lrgeegfc
git submodule update --init --recursive

# Create and activate conda environment
conda env create -f lapbrain.yml
conda activate lapbrain

# Build lrgsglib (generates config scripts and conda hooks)
cd lrgsglib
git checkout lrg_eegfc
CONDA_ENV_NAME=lapbrain make all
pip install -e .
cd ..

# Install lrg-eegfc
pip install -e .
```

## Installation

The project targets Python **3.12** (as specified in `lapbrain.yml`). Follow
these steps to set up the complete environment.

### Step 1: Clone the repository

```bash
git clone https://github.com/giulioiannelli/lrgeegfc.git
cd lrgeegfc
```

### Step 2: Initialize submodules

The project depends on [`lrgsglib`](https://github.com/giulioiannelli/lrgsglib),
a companion library that provides Laplacian-based graph operators. It is
included as a Git submodule:

```bash
git submodule update --init --recursive
```

### Step 3: Create the conda environment

The repository includes a `lapbrain.yml` file that defines all conda
dependencies (numpy, scipy, matplotlib, graph-tool, etc.):

```bash
conda env create -f lapbrain.yml
conda activate lapbrain
```

**Custom installation path:** To install the environment in a specific
location (e.g., on a larger disk or shared filesystem), use the `--prefix`
option instead of the default location:

```bash
conda env create -f lapbrain.yml --prefix /path/to/custom/envs/lapbrain
conda activate /path/to/custom/envs/lapbrain
```

If the environment already exists and you want to update it:

```bash
conda env update -f lapbrain.yml --prune
```

### Step 4: Build and configure lrgsglib

The `lrgsglib` submodule must be built to generate environment configuration
scripts and set up conda activation hooks. Switch to the correct branch and
run the build:

```bash
cd lrgsglib
git checkout lrg_eegfc
CONDA_ENV_NAME=lapbrain make all
```

This command:
- Generates `config_env.sh` and `unconfig_env.sh` with paths rooted at the
  current directory
- Creates conda activation/deactivation hooks that automatically export
  environment variables (`LRGSG_ROOT`, `LRGSG_DATA`, etc.) when you activate
  the `lapbrain` environment
- Compiles any required C extensions

**Custom conda prefix:** If you installed the conda environment with
`--prefix`, set the `CONDA_PREFIX` variable before running make:

```bash
CONDA_PREFIX=/path/to/custom/envs/lapbrain CONDA_ENV_NAME=lapbrain make all
```

After building, install the package in editable mode:

```bash
pip install -e .
cd ..
```

### Step 5: Install lrg-eegfc

Install the main package in editable mode:

```bash
pip install -e .
```

This registers the `lrg-eegfc` command line tool (unified CLI) plus the
legacy `lrg-eegfc-corr` entry point, and makes the `lrg_eegfc` module
available for imports.

For development with linting and testing tools:

```bash
pip install -e ".[dev]"
```

### Verify installation

```bash
python -c "import lrg_eegfc; import lrgsglib; print('OK')"
```

### Developer tools

The `dev` extra installs pytest, black, isort, flake8 and mypy. Run the full
quality gate with:

```bash
pytest
black --check src
isort --check src
flake8 src
mypy src
```

## Dataset layout

The command line tools assume the following directory structure by default:

```
└── data/
    └── raw/
        └── stereoeeg_patients/
            ├── Pat_01/
            │   ├── resting/
            │   │   ├── rest_pre.mat
            │   │   └── rest_post.mat
            │   ├── task/
            │   │   ├── task_learn.mat
            │   │   └── task_test.mat
            │   ├── implant/implant_pat_01.xlsx
            │   ├── implant_pat_01.csv
            │   ├── channel_labels.csv
            │   └── provenance.md
            └── …
```

See `.agents/guides/03_implementation/DATA_LAYOUT.md` for the full canonical
spec + per-patient quirks + the `lrg-eegfc data normalize` CLI.

Use ``--dataset-root`` to point to a different directory when running the CLI.
All generated artefacts (correlation matrices and plots) are written to
``data/correlations/<patient>/`` by default.

## Command line usage

The package exposes a unified `lrg-eegfc` CLI with subcommands for compute,
plot, show, data, cache, config, and bundle. The legacy `lrg-eegfc-corr`
entry point remains for backward-compatible correlation workflows.

Unified CLI examples:

```bash
# Compute MSC matrices (cached)
lrg-eegfc compute msc --patients Pat_02 --band alpha --phase rest_pre -v

# Run LRG analysis from cached FC matrices
lrg-eegfc compute lrg --patients Pat_02 --fc-method msc -v

# Plot LRG panels from cache
lrg-eegfc plot lrg --patient Pat_02 --fc-method msc --plot-type full -v
```

Legacy CLI example (correlation-only):

```bash
lrg-eegfc-corr \
    --patient Pat_01 \
    --phase rest_pre \
    --band beta \
    --dataset-root /path/to/data/stereoeeg_patients \
    --plot-all
```

Key options:

* ``--filter-time`` – restrict the analysis to the first *N* samples.
* ``--jump-index`` – choose which percolation jump to use when selecting the
  correlation threshold (0 = single giant component).
* ``--channel-names`` – path to a ``ChannelNames.mat`` file with the
  ``ChannelNames`` variable, used to label the dendrogram/graph plots.
* ``--plot-*`` flags – enable individual plots; ``--plot-all`` toggles every
  available plot.

Run ``lrg-eegfc --help`` for the full list of command groups and options.

Pipeline scripts used by the shell runners live in `scripts/py/` and are
invoked by `scripts/run_step.sh` and `scripts/run_full_analysis.sh`.

## Python API overview

The public API lives in :mod:`lrg_eegfc`.  Highlights include:

| Function | Description |
| --- | --- |
| ``load_timeseries(patient, phase, root_path)`` | Load a SEEG recording into a ``(channels, samples)`` array. |
| ``load_patient_dataset(patient, root_path)`` | Return a dictionary of ``phase -> PatientRecording`` objects. |
| ``build_correlation_network(timeseries, threshold=...)`` | Produce a processed correlation matrix. |
| ``build_band_correlation_matrices(data_ts, fs)`` | Compute per-band correlation matrices with automatic threshold selection. |
| ``compute_band_connectivity(patient, phase, band, dataset_root)`` | Convenience wrapper that ties together loading, filtering and threshold selection. |

Refer to the in-code docstrings for full parameter documentation.

## Plotting utilities

The :mod:`lrg_eegfc.plotting` module contains the functions that back the CLI
plots (`plot_correlation_matrix`, `plot_entropy`, `plot_dendrogram`,
`plot_graph`).  They accept plain ``pathlib.Path`` destinations so they can be
used interactively inside notebooks.

## Developer guide

* New code must include type hints and informative docstrings.
* Keep imports explicit – wildcard imports are intentionally avoided.
* The ``docs/`` directory contains extended documentation for architecture and
  contribution guidelines.
* Use ``pip install -e .[dev]`` to bring in the linting and formatting tools.

Issues and pull requests are welcome!  Please open an issue with a reproducible
example when reporting bugs.
