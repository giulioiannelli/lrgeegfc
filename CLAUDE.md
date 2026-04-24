# Agent Notes — lrgeegfc

**Read this first. Everything you need in one screen.**

- **Current era:** `IMCOH_ABS` × `COHORT_N9` (9 patients, locked 2026-04-22).
  `imcoh_abs = <|ImCoh|>_f`.
- **Core scientific question:** does `task_test` leave a band-specific,
  multiscale structural trace in `rest_post` LRG dendrograms that is
  cohort-wide (≥ 7/9)? Signal is already visible in existing VI(k) /
  `h2_partition_multiscale` / `h2d_coactivation_persistence` artifacts —
  surface it, don't re-test.
- **Start points:**
  - `.agents/START_HERE.md` — current state, 3 entry paths.
  - `.agents/reports/2026-04-24_pipeline-status.md` — era index.
  - `.agents/reports/2026-04-24_multiscale-task-trace.md` — writing handoff.
  - `.agents/era-map.md` — MSC / IMCOH_SQ / IMCOH_ABS / COHORT_N9 landmarks.
  - `.agents/diary/` — what happened, day by day.

`CLAUDE.md` and `AGENTS.md` must stay identical.

---

## Core rules

### Renormalization communication (read `.agents/guides/01_project/renormalization-style.md`)

Every output — report, commit, diary entry, chat reply, memory — leads
with a 1–2 sentence head (the *juice*). Technical detail goes below,
only as deep as asked. High-level plain language first; zoom in only
when the user requests or something diverges from expectation.

Head-first is a **summary contract**: the head must accurately
summarize what follows. A crisp head over a hand-waved body is worse
than an honest long-form draft.

### Library-first (full rules: `.agents/guides/04_rules/coding-rules.md`)

- General code lives in `src/lrg_eegfc/`. Dataset-specific scripts in
  `scripts/`.
- Before writing any helper in a script: `rg src/lrg_eegfc/` first.
- A helper with ≥2 callers is promoted to the library in the same
  commit as its second use.
- No private-copy forks of `load_fc_matrix`, `wilcoxon_z`, plot
  helpers. Import them.
- FC-method-agnostic defaults. Config-driven constants. No hardcoded
  `data/...` paths — use `lrg_eegfc.config.paths`.

### Never / always list (full: `.agents/guides/04_rules/never-always-list.md`)

**Never**
- Never use `lrg.optimal_threshold` as a diffusion time τ.
- Never call `fig.suptitle` on publication figures.
- Never mock FC data in hypothesis-level tests.
- Never run scripts outside the `lapbrain` conda env.
- Never start a new scalar hypothesis test when the signal is visible
  in existing VI(k) / partition-multiscale / H2d-θ artifacts.
- Never delete files that document research history — `git mv` to
  `<parent>/archive/YYYY-MM/`.
- Never pool metrics into a consensus scalar (user forbidden).
- Never skip frontmatter on a new `.agents/` .md file.
- Never invent metric names — cite literature or existing code.

**Always**
- Always show ≥ 3 patients / bands / phases in published figures.
- Always route data loading through `workflow.fc.load_fc_matrix`.
- Always import statistical helpers from
  `lrg_eegfc.utils.metrics.hypothesis` (`wilcoxon_z`, `bh_fdr`,
  `rank_biserial`, `boot_ci_mean`, `cluster_stats`).
- Always set dendrogram y-limits as
  `tmin = merge_heights[0]*0.8, tmax = merge_heights[-1]*1.05`.
- Always zoom nilearn glass-brain panels to electrode bbox.
- Always flag Pat_03 as 1024 Hz outlier (include, mark distinctly).
- Always exclude Pat_14 `task_test` (corrupt file).
- Always drop Pat_10 task rows `[53, 54, 55]` at load.
- Always add frontmatter to new `.agents/` .md files.
- Always write a renormalization-style head before any body.

### Memory meta-rule

When the user says **"never X"** or **"always Y"**, the rule is added
to `.agents/guides/04_rules/never-always-list.md` AND a matching
`feedback_<short>.md` memory is saved on first mention. This list is
the single source of truth for agent behavior.

---

## Quick reference map

### Skills (invoke with `/name`)

Existing:
- `/figure` — Generate figure code with templates
- `/patient <id>` — Run full pipeline for patient
- `/lrg` — LRG analysis with customization
- `/cache` — Manage cache files
- `/validate` — Run tests and checks
- `/style` — Apply plot styling
- `/notebook` — Create / execute notebooks
- `/data` — Inspect patient data

Added 2026-04 (reorg):
- `/era <path>` — Report era + status from frontmatter
- `/diary` — Append a stamped block to today's `.agents/diary/`
- `/surface <claim>` — Read cached results instead of recomputing
- `/audit <script>` — Library-reuse check on a script

### Figure generation

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

### Data loading (always use cache)

| Task | Module | Function |
|------|--------|----------|
| **Load any FC matrix** | **`workflow.fc`** | **`load_fc_matrix(patient, phase, band, fc_method)`** |
| Load correlation | `workflow.corr` | `load_corr_matrix(patient, phase, band, cache_root)` |
| Load MSC | `workflow.msc` | `load_msc_matrix(patient, phase, band, cache_root, sparsify, n_surrogates, nperseg)` |
| Load LRG result | `workflow.lrg` | `load_lrg_result(patient, phase, band, fc_method)` |
| Load timeseries | `utils.io.patient` | `load_timeseries(patient, phase, root_path)` |
| Inspect data | `utils.io.inspect` | `inspect_patient(patient, root_path)` |

**FC method routing:** `fc_method` ∈ `{"corr", "msc", "imcoh", "imcoh_abs", "imcoh_sq"}`. Cache dirs and LRG cache roots are selected automatically — never hardcode paths.

**ImCoh taxonomy (post-2026-04-15 reset):**
- `imcoh` = signed Nolte-2004 `Im(S)/√(S_ii·S_jj)`, range `[-1, 1]`. Cached freq-resolved `(N, N, F_band)`. **Cannot feed LRG** (Laplacian needs non-negative) — raises `ValueError`.
- `imcoh_abs` = `<|ImCoh|>_f = mean(|signed|, axis=F_band)` at load time. Connectivity-strength magnitude (Ewald 2012 / Bastos & Schoffelen 2016). **Default; current era.**
- `imcoh_sq` = `<|ImCoh|²>_f = mean(signed², axis=F_band)` at load time. Superseded by `imcoh_abs` post-reset.

Order-of-operations matters (Jensen's inequality): `|mean(signed)| ≠ mean(|signed|)` and `mean(signed)² ≠ mean(signed²)`. The freq-resolved cache applies the transform **per-frequency-bin first**, then band-average.

### Hypothesis-test helpers (library)

```python
from lrg_eegfc.utils.metrics.hypothesis import (
    wilcoxon_z, rank_biserial, boot_ci_mean, bh_fdr, cluster_stats,
)
```

### Computing (creates cache)

| Task | Module | Function |
|------|--------|----------|
| Compute correlation | `workflow.corr` | `compute_corr_matrix(...)` |
| Compute MSC | `workflow.msc` | `compute_msc_matrix(...)` |
| Compute LRG | `workflow.lrg` | `compute_lrg_analysis(fc_matrix, ...)` |
| Full patient | `workflow.corr` | `compute_corr_for_patient(...)` |

### Configuration

| Item | Location | Access |
|------|----------|--------|
| Frequency bands | `config.const` | `BRAIN_BANDS`, `BRAIN_BANDS_NAMES` |
| Phase labels | `config.const` | `PHASE_LABELS` |
| LaTeX band names | `config.const` | `BRAIN_BAND_TEX_DICT` |
| Sample rate | `config.const` | `DEFAULT_SAMPLE_RATE` (2048 Hz) |
| Surrogates default | `config.const` | `DEFAULT_N_SURROGATES` (200) |

---

## Cache locations

**Data paths are centralized in `config/paths.py`.** Always import from there:
```python
from lrg_eegfc.config.paths import CORR_CACHE, MSC_CACHE, LRG_CACHE, IMCOH_CACHE, IMCOH_LRG_CACHE
```

**`data/` layout (4 top-level categories):**
```
data/
├── raw/stereoeeg_patients/Pat_NN/    # canonical per-patient layout
├── cache/                             # all computation caches
├── reports/                           # scientific output (CSV, MD, figures)
└── outputs/{figures,tables}/          # CLI + publication outputs
```

Full layout + vendor → canonical mapping + per-patient quirks:
[`.agents/guides/03_implementation/data-layout.md`](.agents/guides/03_implementation/data-layout.md)
— **read before touching `data/raw/`**.

**Agent handoffs vs scientific reports:**
- `.agents/reports/` — agent-to-agent briefings + historical handoffs.
- `data/reports/` — scientific outputs (result tables, per-analysis .md).

**Naming sensitivity:**
- MSC: `sparsify-{none|soft}`, `nperseg-{value}`, `n_surrogates-{value}`
- ImCoh: `{band}_{phase}_imcoh_freqresolved_nperseg-{value}.npy`. No surrogates.
- ImCoh LRG: `{band}_{phase}_lrg_imcoh-{abs|sq}.npz`. Transform encoded in name.
- Correlation: `filter-{none|abs}`, `zero_diag-{true|false}`

---

## Technical invariants

1. **Never recompute in visualization** — always load from cache.
2. **Standard notebook header:**
   ```python
   from lrg_eegfc.notebook import *
   move_to_root(pathname="lrgeegfc")
   ```
3. **Per-patient quirks** in [`.agents/guides/03_implementation/data-layout.md`](.agents/guides/03_implementation/data-layout.md) §6.
   Cohort locked at n=9 on 2026-04-22 (Pat_02, 03, 05, 06, 07, 08, 10, 13, 15).
4. **Pat_03 outlier (negative control):** 1024 Hz (others 2048 Hz). MSC 3× higher.
   Include in analyses, mark distinctly in figures, report values separately.
5. **Same-probe MSC bias:** contacts on the same sEEG probe have trivially
   high MSC (2–8× higher). Dominates LRG community structure at coarse
   scales. Any community-level analysis must verify after zeroing same-probe
   edges. Use `load_epileptic_nodes()` from `utils.io.patient`; intersect
   with `channel_labels.csv`.
6. **Always `plt.close(fig)` after saving.**
7. **Create output dirs:** `output_path.parent.mkdir(parents=True, exist_ok=True)`.
8. **nperseg matches sampling rate:**
   ```python
   from lrg_eegfc.config.const import nperseg_for_fs
   nperseg = nperseg_for_fs(fs)  # 4096 at 2048 Hz, 2048 at 1024 Hz (Pat_03)
   ```

---

## Standard imports

```python
# Notebook header
from lrg_eegfc.notebook import *
move_to_root(pathname="lrgeegfc")

# Data paths (always use these)
from lrg_eegfc.config.paths import (
    SEEG_DATAPATH, CORR_CACHE, MSC_CACHE, LRG_CACHE,
    IMCOH_CACHE, IMCOH_LRG_CACHE, FIGURES_ROOT, TABLES_ROOT,
)

# Config
from lrg_eegfc.config import BRAIN_BANDS, PHASE_LABELS, BRAIN_BAND_TEX_DICT

# Workflows — unified FC loader (preferred)
from lrg_eegfc.workflow.fc import load_fc_matrix

# Workflows — method-specific
from lrg_eegfc.workflow import (
    load_corr_matrix, compute_corr_matrix,
    load_msc_matrix, compute_msc_matrix,
    load_lrg_result, compute_lrg_analysis,
)

# Visualization
from lrg_eegfc.visuals import (
    plot_correlation_and_network, plot_msc_and_network,
    plot_lrg_full_panel, plot_phase_reorganization,
    plot_spatial_network_3d,
)

# Statistics (NEW: elevated from scripts/_shared.py on 2026-04-24)
from lrg_eegfc.utils.metrics.hypothesis import (
    wilcoxon_z, rank_biserial, boot_ci_mean, bh_fdr, cluster_stats,
)

# Data loading
from lrg_eegfc.utils.io import load_timeseries, inspect_patient
```

---

## Directory structure

```
src/lrg_eegfc/
├── cli/                 # Unified CLI (lrg-eegfc command, 36 subcommands)
├── config/              # Constants, paths, plotting
├── workflow/            # High-level compute + cache
├── visuals/             # All visualization
├── notebook.py          # Notebook utilities
└── utils/
    ├── scripting.py     # Shared script helpers
    ├── io/              # Data loading
    ├── fc/              # FC computation primitives
    ├── lrg/              # LRG utilities
    ├── metrics/         # Comparison metrics
    │   ├── hypothesis.py   # wilcoxon_z, bh_fdr, rank_biserial, boot_ci_mean, cluster_stats
    │   ├── tree.py         # tree_internal_nodes, h_log_grid, fcluster_at_h_rel
    │   ├── tree_distance.py # KC, Matching-Cluster, weighted RF
    │   └── ...              # vi, reorganization, compare
    └── pipelines/       # Batch orchestration

scripts/
├── 00_prepare/                    # data preparation, tooling (e.g. frontmatter_init.py)
├── 01_compute/
│   ├── hypothesis_tests/          # h2c_*, h2d_*, rigorous_hypothesis_test, etc.
│   ├── diagnostics/               # diag_*, band_k_landscape
│   ├── batch/                     # compute_*, batch_*, report_h1h4_vi
│   └── figures_embedded/          # fig_* scripts producing publication figures
├── 02_visualize/ … 09_surrogate/  # workpackage scripts
└── archive/
    ├── 2026-02_imcoh-dev-notes/
    ├── 2026-04_imcoh-sq-era/
    └── 2026-04_failed-scalar-session/

.agents/
├── README.md / START_HERE.md / era-map.md
├── diary/YYYY-MM-DD.md            # one file per day, append per session
├── guides/
│   ├── 01_project/  02_methods/  03_implementation/
│   └── 04_rules/                  # renormalization, coding, naming, frontmatter, never/always
├── plans/{active,developed,archive}/
└── reports/
    ├── 2026-04-24_*.md            # current
    └── archive/2026-04/           # superseded + dead, with era in frontmatter
```

---

## CLI quick reference

The `lrg-eegfc` command provides 36 subcommands across 7 groups.
See `.agents/guides/03_implementation/cli-reference.md` for full docs.

```bash
# Common workflows (--fc-method accepts: corr, msc, imcoh, imcoh_abs, imcoh_sq)
lrg-eegfc compute msc --patients Pat_02 --band alpha --phase rest_pre -v
lrg-eegfc compute lrg --patients Pat_02 --fc-method imcoh -v
lrg-eegfc plot lrg --patient Pat_02 --fc-method imcoh --plot-type full -v

# Query cached results (no figures)
lrg-eegfc show lrg --patient Pat_02 --phase rest_pre --fc-method imcoh

# Inspection and cache
lrg-eegfc data inspect --patients Pat_02 -v
lrg-eegfc data normalize                # dry-run
lrg-eegfc data normalize --apply
lrg-eegfc cache list
lrg-eegfc config show
```

Groups: `compute` (6), `plot` (16), `show` (4), `data` (4), `cache` (4), `config` (2), `bundle` (1).

---

## Detailed guides

**Start points (renormalization order):**
1. `.agents/START_HERE.md` — current state + entry paths
2. `.agents/era-map.md` — era landmarks with what-invalidated-what
3. `.agents/reports/2026-04-24_pipeline-status.md` — era index of all artifacts
4. `.agents/reports/2026-04-24_multiscale-task-trace.md` — writing handoff
5. `.agents/reports/2026-04-24_h1-h4-vi-results.md` — canonical VI(k)
6. `.agents/reports/2026-04-24_post-mortem-scalar-session.md` — April lessons

**Rules:**
- `.agents/guides/01_project/renormalization-style.md` — communication
- `.agents/guides/04_rules/coding-rules.md` — library-first, no-duplication
- `.agents/guides/04_rules/naming-conventions.md` — kebab + date + frontmatter
- `.agents/guides/04_rules/frontmatter-schema.md` — YAML schema
- `.agents/guides/04_rules/never-always-list.md` — enforced preferences

**Methods:**
- `.agents/guides/02_methods/imcoh-guide.md` — volume-conduction-immune FC
- `.agents/guides/02_methods/probe-bias-guide.md` — **CRITICAL** same-probe bias
- `.agents/guides/02_methods/h2-metrics.md` — H2-family definitions
- `.agents/guides/02_methods/msc-method-guide.md` — MSC (superseded)
- `.agents/guides/02_methods/time-window-guide.md` — sliding-window FC

**Implementation:**
- `.agents/guides/03_implementation/data-layout.md` — per-patient layout + quirks
- `.agents/guides/03_implementation/cli-reference.md` — CLI command reference
- `.agents/guides/03_implementation/function-map.md` — function lookup
- `.agents/guides/03_implementation/figure-patterns.md` — figure templates
- `.agents/guides/03_implementation/caching-guide.md` — cache structure

**Project workflow:**
- `.agents/guides/01_project/agent-playbook.md` — session workflow
- `.agents/guides/01_project/agent-structure-guide.md` — repo layout
- `.agents/guides/01_project/agent-tasks.md` — common task patterns

**Historical (archive — read for context, don't cite):**
- `.agents/reports/archive/2026-04/` — pre-reset |ImCoh|² + scalar-session artifacts
- `.agents/plans/archive/` — MSC-era + |ImCoh|²-era plans
- `scripts/archive/2026-04_failed-scalar-session/` — dead-end stage scripts
