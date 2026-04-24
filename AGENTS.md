# Agent Notes - lrgeegfc

Canonical agent guidance lives in `.agents/README.md`.
Root `AGENTS.md` and `CLAUDE.md` should remain identical.

---

## Quick Reference Map

Use this section to quickly locate functions and modules.

### Skills Available

Invoke with `/skillname`:
- `/figure` - Generate figure code with templates
- `/patient <id>` - Run full pipeline for patient
- `/lrg` - LRG analysis with customization
- `/cache` - Manage cache files
- `/validate` - Run tests and checks
- `/style` - Apply plot styling
- `/notebook` - Create/execute notebooks
- `/data` - Inspect patient data

### Figure Generation

| Type | Module | Key Function |
|------|--------|--------------|
| Correlation matrix | `visuals.correlation` | `plot_correlation_and_network()` |
| MSC matrix | `visuals.msc` | `plot_msc_and_network()` |
| LRG full panel | `visuals.lrg` | `plot_lrg_full_panel()` |
| Dendrogram | `visuals.lrg` | `plot_lrg_dendrogram()` |
| Entropy curves | `visuals.lrg` | `plot_lrg_entropy_curves()` |
| 3D brain (Plotly) | `visuals.spatial` | `plot_spatial_network_3d()` |
| 3D brain (nilearn) | `visuals.spatial` | `view_brain_connectome()` |
| Phase comparison | `visuals.reorganization` | `plot_phase_reorganization()` |

### Data Loading (always use cache)

| Task | Module | Function |
|------|--------|----------|
| **Load any FC matrix** | **`workflow.fc`** | **`load_fc_matrix(patient, phase, band, fc_method)`** |
| Load correlation | `workflow.corr` | `load_corr_matrix(patient, phase, band, cache_root)` |
| Load MSC | `workflow.msc` | `load_msc_matrix(patient, phase, band, cache_root, sparsify, n_surrogates, nperseg)` |
| Load LRG result | `workflow.lrg` | `load_lrg_result(patient, phase, band, fc_method)` |
| Load timeseries | `utils.io.patient` | `load_timeseries(patient, phase, root_path)` |
| Inspect data | `utils.io.inspect` | `inspect_patient(patient, root_path)` |

**FC method routing:** `fc_method` ∈ `{"corr", "msc", "imcoh", "imcoh_abs", "imcoh_sq"}`. Cache directories and LRG cache roots are selected automatically — never hardcode paths.

**ImCoh taxonomy (post-2026-04-15 reset):**
- `imcoh` = signed Nolte-2004 `Im(S)/√(S_ii·S_jj)`, range `[-1, 1]`. Cached on disk as **freq-resolved** `(N, N, F_band)` per (patient, phase, band). **Cannot feed LRG** (Laplacian needs non-negative) — raises `ValueError`.
- `imcoh_abs` = `<|ImCoh|>_f` = `mean(|signed|, axis=F_band)` at load time. Connectivity-strength magnitude (Ewald 2012 / Bastos & Schoffelen 2016). Default for Section 2.
- `imcoh_sq` = `<|ImCoh|²>_f` = `mean(signed², axis=F_band)` at load time. Squared imaginary coherence (Ewald 2012). Matches pre-reset archive.

Order-of-operations matters (Jensen's inequality): `|mean(signed)| ≠ mean(|signed|)` and `mean(signed)² ≠ mean(signed²)`. The freq-resolved cache lets the loader apply the transform **per-frequency-bin first**, then band-average.

### Computing (creates cache)

| Task | Module | Function |
|------|--------|----------|
| Compute correlation | `workflow.corr` | `compute_corr_matrix(patient, phase, band, ...)` |
| Compute MSC | `workflow.msc` | `compute_msc_matrix(patient, phase, band, ...)` |
| Compute LRG | `workflow.lrg` | `compute_lrg_analysis(fc_matrix, ...)` |
| Full patient | `workflow.corr` | `compute_corr_for_patient(patient, ...)` |

### Configuration

| Item | Location | Access |
|------|----------|--------|
| Frequency bands | `config.const` | `BRAIN_BANDS`, `BRAIN_BANDS_NAMES` |
| Phase labels | `config.const` | `PHASE_LABELS` |
| LaTeX band names | `config.const` | `BRAIN_BAND_TEX_DICT` |
| Sample rate | `config.const` | `DEFAULT_SAMPLE_RATE` (2048 Hz) |
| Surrogates default | `config.const` | `DEFAULT_N_SURROGATES` (200) |

---

## Cache Locations

**Data paths are centralized in `config/paths.py`.** Always import from there:
```python
from lrg_eegfc.config.paths import CORR_CACHE, MSC_CACHE, LRG_CACHE, IMCOH_CACHE, IMCOH_LRG_CACHE
```

**`data/` layout (4 top-level categories only):**
```
data/
├── raw/stereoeeg_patients/Pat_NN/    # sEEG raw data — canonical per-patient layout:
│   ├── resting/{rest_pre,rest_post}.mat   + .info.txt (vendor metadata)
│   ├── task/{task_learn,task_test}.mat    + .info.txt (optional)
│   ├── implant/implant_pat_NN.xlsx         # red-font cells = epileptic contacts
│   ├── channel_labels.csv                  # single labels file, no .txt/.mat duplicates
│   ├── implant_pat_NN.csv                  # derived from xlsx by the normalizer
│   └── provenance.md                       # vendor→canonical mapping + sha256
├── cache/                             # all computation caches
│   ├── corr/, msc/, lrg/             # standard FC + LRG
│   ├── imcoh/, imcoh_lrg/            # ImCoh FC + LRG
│   ├── msc_dev/, lrg_crema/          # variants
│   ├── cleaned_corr/, fc_fig/        # cleaned, figure-ready
│   ├── metric_concordance/, surrogate_validation/
│   └── bipolar/, rescaled/           # experimental (failed approaches)
├── reports/                           # investigation outputs by topic
│   ├── metric_exploration/, scalar_vi/, spatial/
│   ├── clinical_application/, coclassification/
│   └── imcoh/, imcoh_unanimity/, imcoh_vi/
└── outputs/
    ├── figures/                       # CLI figure output
    └── tables/                        # CSV, tex, npz tables
```

The per-patient layout is enforced by `lrg-eegfc data normalize` (dry-run by
default; pass `--apply` to execute). Full spec, vendor → canonical mapping
rules, and per-patient quirks are documented in
[`.agents/guides/03_implementation/DATA_LAYOUT.md`](.agents/guides/03_implementation/DATA_LAYOUT.md)
— **read it before touching any file under `data/raw/`**.

**Naming sensitivity:**
- MSC: includes `sparsify-{none|soft}`, `nperseg-{value}`, `n_surrogates-{value}`
- ImCoh: per-band `{band}_{phase}_imcoh_freqresolved_nperseg-{value}.npy` — freq-resolved tensor `(N, N, F_band)`. No surrogates (they scramble phase and break the magnitude of ImCoh).
- ImCoh LRG: `{band}_{phase}_lrg_imcoh-{abs|sq}.npz` — transform encoded because LRG output depends on it.
- Correlation: includes `filter-{none|abs}`, `zero_diag-{true|false}`

---

## Critical Invariants

1. **Never recompute in visualization** - always load from cache
2. **Standard notebook header:**
   ```python
   from lrg_eegfc.notebook import *
   move_to_root(pathname="lrgeegfc")
   ```
3. **Per-patient quirks** are centrally documented in
   [`.agents/guides/03_implementation/DATA_LAYOUT.md`](.agents/guides/03_implementation/DATA_LAYOUT.md) §6.
   As of 2026-04-22 the full n=10 roster has all 4 phases — Pat_06 was
   previously task-less but has been re-completed, so the old "exclude Pat_06
   for cross-phase" rule no longer applies.
3b. **Pat_03 is a documented outlier (negative control):** Pat_03 was recorded at 1024 Hz
   (all others at 2048 Hz). Its MSC coherence is 3× higher than any other patient,
   producing abnormally dense FC and unstable LRG hierarchies. Include Pat_03 in all
   analyses but flag it as outlier. When Pat_03 diverges from the group pattern, this
   is expected and serves as a negative control validating the method's sensitivity to
   data quality. Always show Pat_03 in figures (marked distinctly) and report its values
   separately in results tables.
4. **Same-probe MSC bias — CRITICAL:** Contacts on the same sEEG probe
   (e.g., A1-A12) have trivially high MSC due to volume conduction and
   spatial proximity. Same-probe MSC is 2-8× higher than cross-probe.
   This dominates LRG community structure at coarse scales: at n=30
   communities, 40-80% of same-community pairs are same-probe (vs ~10%
   expected). **Any community-level analysis must be verified after
   zeroing same-probe edges.** Use `load_epileptic_nodes()` from
   `utils.io.patient` to get epileptic labels; always intersect with
   `channel_labels.csv` since many epileptic contacts are not in the
   network. When comparing epileptic vs non-epileptic, use within-probe
   controls (contacts on the same probe) AND verify with debiased FC.
   The VI metric is misleading here — use same-community/same-probe
   enrichment fraction instead.
5. **Always close figures:** `plt.close(fig)` after saving
6. **Create output dirs:** `output_path.parent.mkdir(parents=True, exist_ok=True)`
7. **nperseg must match sampling rate:** Use `nperseg_for_fs(fs)` from `config.const`
   to compute the correct nperseg for any patient. This gives 2-second segments
   (Δf = 0.5 Hz) regardless of sampling rate. Pat_03 was recorded at 1024 Hz
   (all others at 2048 Hz), so Pat_03 uses `nperseg=2048` while others use `nperseg=4096`.
   The cached data under `nperseg-2048` is the correct version for Pat_03.
   ```python
   from lrg_eegfc.config.const import nperseg_for_fs
   nperseg = nperseg_for_fs(fs)  # 4096 at 2048 Hz, 2048 at 1024 Hz
   ```

---

## Standard Imports

```python
# Notebook header
from lrg_eegfc.notebook import *
move_to_root(pathname="lrgeegfc")

# Data paths (always use these, never hardcode "data/...")
from lrg_eegfc.config.paths import (
    SEEG_DATAPATH, CORR_CACHE, MSC_CACHE, LRG_CACHE,
    IMCOH_CACHE, IMCOH_LRG_CACHE, FIGURES_ROOT, TABLES_ROOT,
)

# Config
from lrg_eegfc.config import BRAIN_BANDS, PHASE_LABELS, BRAIN_BAND_TEX_DICT

# Workflows — unified FC loader (preferred)
from lrg_eegfc.workflow.fc import load_fc_matrix

# Workflows — method-specific (when you need extra kwargs)
from lrg_eegfc.workflow import (
    load_corr_matrix, compute_corr_matrix,
    load_msc_matrix, compute_msc_matrix,
    load_lrg_result, compute_lrg_analysis,
)

# Visualization
from lrg_eegfc.visuals import (
    plot_correlation_and_network,
    plot_msc_and_network,
    plot_lrg_full_panel,
    plot_phase_reorganization,
    plot_spatial_network_3d,
)

# Data loading
from lrg_eegfc.utils.io import load_timeseries, inspect_patient
```

---

## Directory Structure

```
src/lrg_eegfc/
├── __init__.py          # Package exports
├── notebook.py          # Notebook utilities
├── cli/                 # Unified CLI (lrg-eegfc command)
│   ├── __init__.py      #   exports `app`
│   ├── _app.py          #   LazyGroup root + 7 lazy subcommands
│   ├── _common.py       #   shared options, resolvers, CliReporter
│   ├── compute.py       #   6 compute commands
│   ├── plot.py          #   16 plot commands
│   ├── show.py          #   4 show commands (query cached results)
│   ├── data.py          #   3 data inspection commands
│   ├── cache.py         #   4 cache management commands
│   ├── config_cmd.py    #   2 config display commands
│   ├── bundle.py        #   1 publication bundle command
│   └── _legacy.py       #   old cli.py (lrg-eegfc-corr compat)
├── config/              # Constants, plotting config
│   ├── const.py         # Band/phase/patient constants
│   ├── paths.py         # **Centralized data paths** (all cache/output dirs)
│   └── plotlib.py       # Matplotlib configuration
├── workflow/            # High-level compute + cache
│   ├── corr.py          # Correlation workflow
│   ├── msc.py           # MSC / ImCoh workflow
│   ├── lrg.py           # LRG analysis (fc_method: corr/msc/imcoh)
│   ├── cleaning.py      # Marchenko-Pastur cleaning
│   ├── diagnostics.py   # LRG diagnostics
│   └── time_windows.py  # Sliding-window FC
├── visuals/             # All visualization
│   ├── correlation.py   # Correlation plots
│   ├── msc.py           # MSC plots
│   ├── lrg.py           # LRG/dendrogram plots
│   ├── spatial.py       # 3D brain visualization
│   ├── reorganization.py # Phase comparison
│   └── metastable.py   # Sankey diagrams
└── utils/               # Low-level utilities
    ├── scripting.py     # **Shared script helpers** (setup_script_env, save_figure, etc.)
    ├── io/              # Data loading
    ├── fc/              # FC computation primitives (corr, msc, imcoh)
    ├── lrg/             # LRG utilities
    ├── metrics/         # Comparison metrics
    └── pipelines/       # Batch orchestration
```

---

## CLI Quick Reference

The `lrg-eegfc` command provides 36 subcommands across 7 groups.
See `.agents/guides/CLI_REFERENCE.md` for complete documentation.

```bash
# Common workflows (--fc-method accepts: corr, msc, imcoh)
lrg-eegfc compute msc --patients Pat_02 --band alpha --phase rest_pre -v
lrg-eegfc compute lrg --patients Pat_02 --fc-method msc -v
lrg-eegfc compute lrg --patients Pat_02 --fc-method imcoh -v   # ImCoh LRG
lrg-eegfc plot lrg --patient Pat_02 --fc-method imcoh --plot-type full -v
lrg-eegfc plot corr --patient Pat_02 --plot-type all -v

# Query cached results (no figures)
lrg-eegfc show lrg --patient Pat_02 --phase rest_pre --fc-method imcoh
lrg-eegfc show msc --patient Pat_02 --band alpha --phase rest_pre

# Inspection and cache management
lrg-eegfc data inspect --patients Pat_02 -v
lrg-eegfc data normalize                        # dry-run: preview per-patient layout actions
lrg-eegfc data normalize --apply                # execute: enforce canonical layout across all patients
lrg-eegfc cache list
lrg-eegfc config show
```

Command groups: `compute` (6), `plot` (16), `show` (4), `data` (4), `cache` (4), `config` (2), `bundle` (1).
Full data-layout spec + per-patient quirks:
[`.agents/guides/03_implementation/DATA_LAYOUT.md`](.agents/guides/03_implementation/DATA_LAYOUT.md).

---

## Detailed Guides

**START HERE:** `.agents/reports/PIPELINE_STATUS.md` — era index flagging
every artefact as current / superseded / dead-end after the
MSC → |ImCoh|² → |ImCoh| transition. Canonical H1-H4 VI(k) results:
`.agents/reports/H1_H4_VI_RESULTS_POST_RESET.md`.

See `.agents/guides/INDEX.md` for the full index. Key references:

**Implementation:**
- `.agents/guides/03_implementation/CLI_REFERENCE.md` - CLI command reference
- `.agents/guides/03_implementation/CACHING_GUIDE.md` - Cache structure and data paths
- `.agents/guides/03_implementation/FUNCTION_MAP.md` - Function lookup table
- `.agents/guides/03_implementation/FIGURE_PATTERNS.md` - Figure templates

**Methods:**
- `.agents/guides/02_methods/IMCOH_GUIDE.md` - ImCoh: volume-conduction-immune FC
- `.agents/guides/02_methods/PROBE_BIAS_GUIDE.md` - **CRITICAL**: same-probe MSC bias
- `.agents/guides/02_methods/MSC_METHOD_GUIDE.md` - MSC methodology

**Reports:**
- `.agents/reports/IMCOH_PROCESS_REPORT.md` - Full discovery report with results
- `.agents/reports/IMCOH_RESULTS_FOR_WRITING.md` - MSC vs ImCoh comparison for writing
- `.agents/reports/IMCOH_VERIFICATION_RESULTS.md` - Verified cell counts and flags

**Project:**
- `.agents/guides/01_project/AGENT_PLAYBOOK.md` - Session workflow
