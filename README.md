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
  - [Step 4: Configure environment activation hook](#step-4-configure-environment-activation-hook)
  - [Step 5: Install lrgsglib](#step-5-install-lrgsglib)
  - [Step 6: Install lrg-eegfc](#step-6-install-lrg-eegfc)
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

# Configure environment activation hook
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
echo "source $(pwd)/lrgsglib/tools/bash/config_env.sh" \
    > "$CONDA_PREFIX/etc/conda/activate.d/custom_env_setup.sh"
chmod +x "$CONDA_PREFIX/etc/conda/activate.d/custom_env_setup.sh"

# Install dependencies in editable mode
cd lrgsglib && git checkout lrg_eegfc && cd ..
pip install -e ./lrgsglib
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

### Step 4: Configure environment activation hook

The `lrgsglib` submodule provides a shell script that sets up useful
environment variables. Configure conda to source it automatically on
activation:

```bash
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
echo "source $(pwd)/lrgsglib/tools/bash/config_env.sh" \
    > "$CONDA_PREFIX/etc/conda/activate.d/custom_env_setup.sh"
chmod +x "$CONDA_PREFIX/etc/conda/activate.d/custom_env_setup.sh"
```

This creates a hook that runs each time you activate the `lapbrain`
environment, exporting paths like `LRGSG_ROOT`, `LRGSG_DATA`, etc.

### Step 5: Install lrgsglib

The `lrgsglib` submodule must be installed in editable mode on the correct
branch:

```bash
cd lrgsglib
git checkout lrg_eegfc
cd ..
pip install -e ./lrgsglib
```

### Step 6: Install lrg-eegfc

Install the main package in editable mode:

```bash
pip install -e .
```

This registers the `lrg-eegfc-corr` command line tool and makes the
`lrg_eegfc` module available for imports.

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
    └── stereoeeg_patients/
        ├── Pat_01/
        │   ├── rsPre.mat
        │   ├── taskLearn.mat
        │   ├── taskTest.mat
        │   ├── rsPost.mat
        │   ├── Implant_pat_01.csv
        │   └── channel_labels.csv
        └── …
```

Use ``--dataset-root`` to point to a different directory when running the CLI.
All generated artefacts (correlation matrices and plots) are written to
``data/correlations/<patient>/`` by default.

## Command line usage

The package exposes the ``lrg-eegfc-corr`` entry point which computes a band-
specific correlation matrix for a patient and optionally creates a set of
plots:

```bash
lrg-eegfc-corr \
    --patient Pat_01 \
    --phase rsPre \
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

Run ``lrg-eegfc-corr --help`` for the full list of options.

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
