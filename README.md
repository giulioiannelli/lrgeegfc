# LRG EEG Functional Connectivity

Tools for computing frequency-specific functional connectivity (FC)
matrices from stereo-EEG (SEEG) recordings and running the Laplacian
Renormalization Group (LRG) hierarchical analysis on the resulting
graphs. Built for the [Living and Relational Graphs (LRG)
project](https://github.com/giulioiannelli), with a Python package, a
unified `lrg-eegfc` CLI, and a publication-grade plotting pipeline.

---

## Current era (locked)

- **FC carrier:** `imcoh_abs` = ⟨|ImCoh(f)|⟩_f — the band-magnitude of
  the Nolte-2004 imaginary coherency (Ewald 2012 convention).
  Volume-conduction immune. MSC is retained for diagnostics only.
- **Cohort:** 10 sEEG patients (`Pat_02, 03, 05, 06, 07, 08, 10, 13,
  14, 15`), locked 2026-04-25.
- **T_d sign convention** (locked 2026-05-26): at every layer (raw FC,
  D_coph, KC, Grassmann) `T_d := d(rest_pre, task) − d(task, rest_post)`,
  **positive = TRACE**. Wilcoxon trace-direction tests use
  `alternative='greater'`.
- **Live preprint hub:** `.agents/preprint/` (see `WRITING_GUIDE.md`).
- **Agent landing page:** `.agents/START_HERE.md`. Era index:
  `.agents/era-map.md`. Project rules: `CLAUDE.md` (mirrored in
  `AGENTS.md`).

---

## Contents

- [Quick start](#quick-start)
- [Installation](#installation)
- [Canonical CLI flow](#canonical-cli-flow)
- [Dataset layout](#dataset-layout)
- [Python API](#python-api)
- [Plotting utilities](#plotting-utilities)
- [Agent + developer guide](#agent--developer-guide)

---

## Quick start

```bash
git clone https://github.com/giulioiannelli/lrgeegfc.git
cd lrgeegfc
git submodule update --init --recursive

conda env create -f lapbrain.yml
conda activate lapbrain

# Build lrgsglib (the submodule provides Laplacian operators)
cd lrgsglib
CONDA_ENV_NAME=lapbrain make all
pip install -e .
cd ..

pip install -e .
```

End-to-end smoke test on one patient (uses cached `imcoh_abs` data):

```bash
lrg-eegfc compute lrg --patients Pat_02 --fc-method imcoh_abs -v
lrg-eegfc show    lrg --patient  Pat_02 --phase rest_pre --fc-method imcoh_abs
lrg-eegfc plot    lrg --patient  Pat_02 --phase rest_pre --band beta --plot-type full
```

The CLI defaults to `imcoh_abs` everywhere as of 2026-05-28; pass
`--fc-method msc` only for explicit diagnostics.

---

## Installation

Python **3.12** (as pinned by `lapbrain.yml`). Steps:

### 1. Clone

```bash
git clone https://github.com/giulioiannelli/lrgeegfc.git
cd lrgeegfc
```

### 2. Initialise submodules

`lrgsglib` (Laplacian + RG operators) is a submodule on `main`:

```bash
git submodule update --init --recursive
```

### 3. Create the conda env

```bash
conda env create -f lapbrain.yml
conda activate lapbrain
```

Custom prefix:

```bash
conda env create -f lapbrain.yml --prefix /path/to/envs/lapbrain
conda activate /path/to/envs/lapbrain
```

Update if it already exists:

```bash
conda env update -f lapbrain.yml --prune
```

### 4. Build `lrgsglib`

```bash
cd lrgsglib
CONDA_ENV_NAME=lapbrain make all
pip install -e .
cd ..
```

### 5. Install `lrg-eegfc`

```bash
pip install -e .
```

`pip install -e .[dev]` adds `pytest`, `black`, `isort`, etc.; the
`jupyter` extra adds JupyterLab; the `docs` extra adds Sphinx.

---

## Canonical CLI flow

The `lrg-eegfc` command exposes seven groups (40 subcommands total
as of 2026-05-28; run `lrg-eegfc --help` for the live tree):

```text
lrg-eegfc compute  ...   # build caches (corr / msc / lrg / cleaning / ...)
lrg-eegfc plot     ...   # generate figures (PDF by default; --png is opt-in)
lrg-eegfc show     ...   # query cached results (no figures)
lrg-eegfc data     ...   # inspect / normalize / compare patient data
lrg-eegfc cache    ...   # list / verify / clean cache files
lrg-eegfc config   ...   # show config + verify paths
lrg-eegfc bundle   ...   # collect figures into Overleaf-ready bundles
```

Common pipeline on `imcoh_abs`:

```bash
# Compute FC + LRG for one patient
lrg-eegfc compute msc --patients Pat_02 --band alpha --phase rest_pre -v
lrg-eegfc compute lrg --patients Pat_02 --fc-method imcoh_abs -v

# Inspect what's cached
lrg-eegfc show lrg --patient Pat_02 --phase rest_pre --fc-method imcoh_abs

# Plot the full LRG panel
lrg-eegfc plot lrg --patient Pat_02 --fc-method imcoh_abs --plot-type full -v
```

Full subcommand reference:
[`.agents/guides/03_implementation/cli-reference.md`](.agents/guides/03_implementation/cli-reference.md).

---

## Dataset layout

Raw + cache + reports + outputs sit under `data/`:

```text
data/
├── raw/stereoeeg_patients/Pat_NN/   # canonical per-patient layout
├── cache/                            # all computation caches
├── reports/                          # scientific outputs (CSV, .md)
└── outputs/{figures,tables}/         # CLI + publication outputs
```

Per-patient quirks (Pat_03 1024 Hz, Pat_10 task row mask, Pat_14
vendor-replaced `task_test.mat`) are documented in
[`.agents/guides/03_implementation/data-layout.md`](.agents/guides/03_implementation/data-layout.md)
and handled at the config layer (`FS_OVERRIDES`,
`PATIENT_CHANNEL_DROP` in `src/lrg_eegfc/config/const.py`). No
analysis-layer special casing.

---

## Python API

Notebook header (canonical as of 2026-05-28):

```python
from lrg_eegfc.notebook import *
move_to_rootf(pathname="lrgeegfc")
```

Most-used entry points:

```python
from lrg_eegfc.config.paths import (
    SEEG_DATAPATH, CORR_CACHE, MSC_CACHE, LRG_CACHE,
    IMCOH_CACHE, IMCOH_LRG_CACHE, FIGURES_ROOT, TABLES_ROOT,
)

# Unified FC loader (handles corr / msc / imcoh / imcoh_abs / imcoh_sq)
from lrg_eegfc.workflow.fc import load_fc_matrix

# Method-specific compute helpers
from lrg_eegfc.workflow import (
    compute_corr_matrix, compute_msc_matrix,
    load_corr_matrix, load_msc_matrix,
    compute_lrg_analysis, load_lrg_result,
)

# Stats (all promoted from script-local copies)
from lrg_eegfc.utils.metrics.hypothesis import (
    wilcoxon_z, rank_biserial, boot_ci_mean, bh_fdr, cluster_stats,
    surrogate_p_value, loo_sensitivity,
)

# Tree / dendrogram helpers
from lrg_eegfc.utils.metrics.tree import (
    cophenet_matrix, induced_linkage, partition_vi_on_subset,
)

# Patient data
from lrg_eegfc.utils.io.patient import (
    PatientRecording, load_timeseries, load_channel_labels,
    load_epileptic_nodes, PatientMasks, build_epi_masks,
)
```

Full function map:
[`.agents/guides/03_implementation/function-map.md`](.agents/guides/03_implementation/function-map.md)
(representative, not exhaustive — source code is canonical).

---

## Plotting utilities

Activate the project mplstyle at the top of every figure script:

```python
from lrg_eegfc.visuals.styles import use_lrg_style
use_lrg_style()
```

Canonical colorbar helper (single-imshow axes only — for row-shared
colorbars keep the explicit `make_axes_locatable` /
`fig.add_axes([...])` pattern):

```python
from lrg_eegfc.visuals import imshow_colorbar_caxdivider
```

Per-class templates live under
[`.agents/guides/05_plotting/`](.agents/guides/05_plotting/README.md)
(FC matrices, dendrograms, network drawings, etc.). Eight enforced
rules at a glance: figure-level legends, canonical colorbar, no
`fig.suptitle` on publication figures, PDF-only output, no
`set_rasterized(True)`, no default watermark, math axis labels (`$i$`,
`$j$`, not "contact" / "channel"), `use_lrg_style()` activated.

---

## Agent + developer guide

The `.agents/` tree is the source of truth for project state, era
tracking, scope reports, plans, preprint writing, and per-day diary
entries. New contributors should read:

- [`.agents/START_HERE.md`](.agents/START_HERE.md) — current state,
  central numerical results, where to begin per role.
- [`.agents/era-map.md`](.agents/era-map.md) — era landmarks +
  what-invalidated-what.
- [`CLAUDE.md`](CLAUDE.md) (mirrored in `AGENTS.md`) — locked project
  rules (renormalization style, terminology, library-first, never /
  always lists).
- [`.agents/guides/04_rules/never-always-list.md`](.agents/guides/04_rules/never-always-list.md)
  — the single source of truth for enforced preferences.

Developer-facing entry points:

- [`docs/overview.md`](docs/overview.md) — architecture.
- [`docs/developer-guide.md`](docs/developer-guide.md) — dev workflow.
- `tests/` — pytest suite. Run `pytest tests/ -q`.

The repository is licensed under **GPL-3.0-or-later**
(see [`LICENSE`](LICENSE)).
