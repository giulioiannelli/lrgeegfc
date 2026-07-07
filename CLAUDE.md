# Agent Notes — lrgeegfc

**Read this first. Everything you need in one screen.**

- **Current era:** `IMCOH_ABS` × `COHORT_N10` (10 patients; n=9 locked 2026-04-22,
  Pat_14 restored 2026-04-25 after vendor task_test replacement).
  `imcoh_abs = <|ImCoh|>_f`.
  Cohort: Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15.
- **Core scientific question:** does `task_test` leave a band-specific,
  multiscale structural trace in `rest_post` LRG dendrograms that is
  cohort-wide (≥ 8/10)? Signal is already visible in existing VI(k) /
  `h2_partition_multiscale` / `h2d_coactivation_persistence` artifacts —
  surface it, don't re-test.
- **Start points:**
  - `.agents/START_HERE.md` — current state, 3 entry paths.
  - `.agents/reports/2026-04-25_task-trace-audit-and-recovery.md` — current writing handoff.
  - `.agents/guides/task-persistence-investigation/2026-04-25_task-trace-canonical.md` — P/T/R/RA reformalization.
  - `.agents/reports/2026-04-24_pipeline-status.md` — era index (n=9 snapshot).
  - `.agents/era-map.md` — MSC / IMCOH_SQ / IMCOH_ABS / COHORT_N9 / COHORT_N10 landmarks.
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

### Terminology — TRACE / ANCHOR / RESET / EMERGENT (full: `.agents/guides/01_project/terminology.md`)

Use the four-way taxonomy when describing how a module behaves
across phases. The bare word "persistence" is ambiguous between two
opposite phenomena (trace vs anchor) and silently flips the reading
every time.

- **trace** — task reorganized AND change persists into RPost. Our
  raw-FC `T_d^(d_S) > 0` finding is a **trace**, never bare "persistence".
- **anchor** — module unchanged across all phases. A "persistent
  module" in the literature is usually an *anchor*, not a trace.
- **reset** — task reorganized AND module reverts in RPost (`T_d < 0`).
- **emergent** — module that did not exist in RPre (LRG / community
  membership only).

Variable names: `n_trace` (not `n_persist`), `n_anchor`, `n_reset`,
`n_emergent`. Figure annotations: "trace: N/10", "trace zone".

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
- **Library names are general, not local-scope.** Module + function
  names under `src/lrg_eegfc/` NEVER reference manuscript-local tokens
  (`section3`, `figureN`, `preprint`, `chapter`, `H2c`, …). They reflect
  general graph / network / statistics / I/O concepts: `network_layouts`,
  `network_drawing`, `tree_metrics`, `surrogate_helpers`, `patient_io`.
  A library helper used by figures-for-Section-3 today must be importable
  by figures-for-Section-7 tomorrow without renaming. Locked 2026-05-28.

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
- Never save figures as both PDF and PNG. **PDF only** is the
  default and only format. PNG is opt-in on explicit user request.
- Never call `im.set_rasterized(True)`. PDFs are **fully vector**
  for every artist (FC matrices, audit heatmaps, scatter plots,
  `pcolormesh`, dendrograms). The earlier "rasterise heavy artists"
  rule is withdrawn. See `.agents/guides/05_plotting/output-and-rasterization.md`.
- Never add a grey provenance watermark / footer by default —
  the file name carries the metadata. Watermark is **opt-in**
  (`watermark=True` kwarg, or
  `lrg_eegfc.visuals.layout.add_provenance_footer(fig, label)`).
- Never label adjacency-matrix axes with "contact" / "channel"
  — use math `$i$`, `$j$`.
- **Never use `imshow_colorbar_caxdivider` for a colorbar that is
  shared across multiple columns in the same row.** That helper
  attaches the colorbar to a *single* axis via `make_axes_locatable`;
  it has no notion of a multi-column shared cbar and will misalign or
  resize the wrong axis. For row-shared colorbars keep the explicit
  `make_axes_locatable` / `fig.add_axes([...])` pattern. The helper is
  the canonical choice for any single-`imshow` axis. Locked 2026-05-28.
- **Never name a library module / function after a manuscript-local
  token** (`section3`, `figureN`, `preprint`, `chapter`, ...). Library
  names reflect general concepts; see the library-first section above.
- Don't confuse the four cross-phase phenomena. **Trace** = task
  changed it AND change stuck (our `T_d > 0` finding). **Anchor** =
  never changed. **Reset** = changed and reverted. **Emergent** =
  never existed before. Use the explicit taxonomy in cross-phase
  taxonomy tables, mixed-band paragraphs, and any context where the
  reader could otherwise read the wrong phenomenon. In unambiguous
  task-trace sections (e.g. β results subsection) bare "persistence"
  is acceptable (softened 2026-05-18). See terminology guide.
- Never frame `d_P = 1 − Pearson(triu A_a, triu A_b)` as "volume +
  topology" or "orthogonal" to `d_S`. It is a magnitude-weighted
  complement; cohort ρ between `T_d^(d_S)` and `T_d^(d_P)` is
  0.85–0.95 per band — strongly correlated, not orthogonal. β is the
  **convergence** cell.
- **Never present any FC-derived cohort claim without first running
  a strength-preserving matched-strength surrogate null.** Within-
  baseline split-half / drift / sampling-jitter nulls are exploratory
  diagnostics, not verification. KC β 10/10 / q=0.006 within-baseline
  headline collapsed on 2026-05-11 because matched-strength reproduced
  ≈50% of the observed tree-distance shift. Matched-strength is the
  *minimum* null required for any cohort-level claim built on a
  connectivity matrix. See `audit_65_kc_matched_strength_verdict.md`,
  `feedback_matched_strength_mandatory.md`.
- **Never produce sycophantic answers or confidence laundering.**
  Default posture is brutal scientific honesty: critical questioning
  of every methodology in play. If a result depends on a control that
  hasn't been run, say so on the first line. If a null is weaker than
  the alternative explanations require, say so. An honest verdict
  of "currently unverified" is preferred over a confident headline
  that gets retracted three iterations later. See
  `feedback_brutal_honesty_no_sycophancy.md`.
- **Never compute a triangle scalar with the OLD T_d sign convention.**
  Locked 2026-05-26: every triangle scalar `T_d` MUST be
  `T_d = d(rest_pre, task) − d(task, rest_post)` so that **T_d > 0 = TRACE**,
  **T_d < 0 = ANTI-TRACE**, **T_d = 0 = NO TRACE**. Applies at every layer
  (raw FC, LRG D_coph, KC, Grassmann). Wilcoxon one-sided trace tests use
  `alternative='greater'`. Per-patient trace counts use `(T > 0).sum()`.
  Surrogate upper-tail p = `mean(s_finite >= obs_T)`. Plot/ylabel/title
  conventions say "positive = trace" (never "negative = trace"). Do not
  reintroduce `d(task, rsPost) − d(rsPre, task)` in any compute site. See
  `feedback_td_sign_convention.md`.

**Always**
- Always show ≥ 3 patients / bands / phases in published figures.
- Always route data loading through `workflow.fc.load_fc_matrix`.
- Always import statistical helpers from
  `lrg_eegfc.utils.metrics.hypothesis` (`wilcoxon_z`, `bh_fdr`,
  `rank_biserial`, `boot_ci_mean`, `cluster_stats`).
- Always set dendrogram y-limits as
  `tmin = merge_heights[0]*0.8, tmax = merge_heights[-1]*1.05`.
- Always zoom nilearn glass-brain panels to electrode bbox.
- Pat_03 is acquired at 1024 Hz (others at 2048 Hz). Sampling-rate
  handling is **config-level only** (`nperseg_for_fs(fs)`, `FS_OVERRIDES`
  in `config/const.py`). Do **not** treat Pat_03 as an outlier, do
  **not** run Pat_03-dropout sensitivity tests, do **not** mark it
  distinctly in figures, do **not** report values separately. Pat_03
  is a full cohort member at n = 10. (Updated 2026-05-18 — the previous
  "1024 Hz outlier / negative control" framing was retired; the
  sampling-rate difference is absorbed at the config layer and does not
  propagate to analysis-level treatment.)
- **Patient-dropout policy: don't drop patients.** Default cohort is
  the full `n = 10`. Two narrow exceptions only: (a) a single
  cohort-level anti-aligned patient at the probe under test, biology-
  driven (canonical: Pat_15 at β LRG, right-hemisphere-only implant);
  (b) genuinely problematic data (vendor corruption, sampling-rate
  handled at config layer — none currently active). Retired dropouts:
  Pat_03 (1024 Hz "outlier"; retired 2026-05-18 am), Pat_07 (substrate
  marginal anti, solidly pro at LRG; retired 2026-05-18 pm), and the
  legacy `n=8` "pro-cohort restriction" replaced by LRG-native `n=9`
  (drop Pat_15 only). See `feedback_no_patient_dropout.md`.
- Pat_14 `task_test` was corrupt at original import; **vendor-replaced 2026-04-25** and is now valid. Cross-phase cohort returns to n=10.
- Always drop Pat_10 task rows `[53, 54, 55]` at load.
- Always add frontmatter to new `.agents/` .md files.
- Always write a renormalization-style head before any body.
- Always file new task-trace / task-persistence investigation tooling
  under `.agents/guides/task-persistence-investigation/` as a
  mathematically rigorous scope report **before writing any code**.
  See that folder's `README.md` for the required structure (notation,
  predicates, formulas, properties, caveats, pseudocode, prior-tool
  connection, open questions).
- Always use the **trace / anchor / reset / emergent** taxonomy
  where cross-phase ambiguity matters (taxonomy tables, mixed-band
  paragraphs, cross-phase summary captions, variable names like
  `n_trace` not `n_persist`). In unambiguous task-trace contexts
  bare "persistence" is acceptable (softened 2026-05-18 — the
  taxonomy is the disambiguation tool, not a vocabulary ban). See
  `.agents/guides/01_project/terminology.md`.
- **Always open any new test / null / methodology with a 5-point
  critical preamble** (script docstring or scope doc, before any code):
  (1) the claim, (2) the null, (3) the strongest plausible alternative
  the null *should* control for, (4) whether the null actually
  controls for it — by mechanism, not vibes, including what it
  *cannot* reject, (5) what would falsify the claim and which
  limitations remain. The point: catch the KC-style mistake at the
  design stage, before two weeks of figures sit on a null that
  doesn't reach the relevant alternative. See
  `feedback_critical_null_preamble.md`.
- **Always state limitations of unverified methods in the first
  paragraph** of any writeup. Any measure not yet matched-strength
  tested is marked **"unverified"** in writeups, memory, and chat,
  and cannot be cited as the load-bearing claim until that control
  is run.
- **Always optimize + time-estimate + progress-surface long compute.**
  Any pipeline that could exceed ~1 min (surrogates, sweeps,
  per-patient×per-band grids, LOO, bootstraps): (1) JIT/vectorize hot
  loops (numba `@njit` with RNG drawn outside the loop → bit-identical;
  cached eigendecompositions) BEFORE launching — never run hours
  unoptimized; (2) time 1–2 units and extrapolate the full runtime,
  stating the estimate BEFORE the full launch; (3) print live
  `[i/N] label elapsed ETA` progress with `flush=True` (block-buffered
  stdout hides un-flushed prints) + a final wall-clock. See
  `feedback_optimize_time_and_surface_progress.md` +
  `feedback_numba_for_surrogates.md`.

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
- `/plotguide` — Open the plotting style guide before producing or
  editing a figure. Read `.agents/guides/05_plotting/` first.

### Plotting (READ BEFORE FIGURES)

Single source of truth: [`.agents/guides/05_plotting/`](.agents/guides/05_plotting/README.md).

Eight rules at a glance:
1. **Shared legends → figure-level**, not axis-level. Use
   `fig.legend(loc="lower center", bbox_to_anchor=(0.5, -0.04),
   ncol=len(handles), frameon=False)` — or
   `lrg_eegfc.visuals.layout.figure_legend(...)`.
2. **Colorbars → `imshow_colorbar_caxdivider`** from
   `lrgsglib.plotlib`. Never `fig.colorbar` on multi-axis grids.
   For LogNorm colorbars end with
   `_apply_factored_sci_format(clb, axis_orientation=...)` (from
   `lrg_eegfc.visuals.fc_templates`) — kills inline `2 × 10ⁿ`
   mantissa labels in both <1.5-decade and ≥1.5-decade regimes.
   **Scope:** use this helper for any *single-imshow* axis. It
   **cannot** serve a colorbar that is shared across multiple columns
   in the same row (it attaches to one axis via `make_axes_locatable`).
   For row-shared cbars keep the explicit `make_axes_locatable` /
   `fig.add_axes([...])` pattern. Locked 2026-05-28.
3. **Multi-axis layout → figure-level decoration.** Titles, legends,
   colorbars, shared axis labels go on the *figure*.
4. **Library-first.** Check `lrgsglib.plotlib` and `lrg_eegfc.visuals`
   before writing a custom plot helper.
5. **PDF only, full vector — never rasterise.** No PNG siblings.
   Do not call `im.set_rasterized(True)` on any artist. Vector is
   sharper and file sizes for typical FC matrices are small.
6. **No `fig.suptitle`** on publication figures. Captions go in a
   sidecar `.md` *only when explicitly asked*; see `captions.md`
   for the plain-language style.
7. **No watermark / provenance footer by default.** The file name
   is the provenance. Watermark is opt-in (`watermark=True` kwarg,
   or `add_provenance_footer(fig, label)`).
8. **Activate the project mplstyle** at the top of every figure
   script:
   `from lrg_eegfc.visuals.styles import use_lrg_style; use_lrg_style()`.
   Single source: `src/lrg_eegfc/visuals/styles/lrg_eegfc.mplstyle`
   (font sizes, tick widths, `pdf.fonttype=42` TrueType embed).
   For one-off overrides use `with rc_context({...}):` — never edit
   the mplstyle for a single figure.

**Per-class templates** (read before writing FC / dendrogram / etc. figures):
- `.agents/guides/05_plotting/fc_templates/` — FC adjacency matrices
  (single, row-per-phase, mosaics). Each template = `.md` style sheet
  + `.py` script. Always check here before plotting a new
  connectivity / `imcoh_abs` / MSC / `corr` matrix.

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
   move_to_rootf(pathname="lrgeegfc")
   ```
3. **Per-patient quirks** in [`.agents/guides/03_implementation/data-layout.md`](.agents/guides/03_implementation/data-layout.md) §6.
   Cohort locked at n=10 on 2026-04-25 (Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15)
   after Pat_14 vendor `task_test` replacement; n=9 snapshot 2026-04-22.
4. **Pat_03 sampling rate:** 1024 Hz (others 2048 Hz). Handled at the
   config layer via `nperseg_for_fs(fs)` and `FS_OVERRIDES` in
   `config/const.py`. Pat_03 is a full cohort member; do not treat as
   outlier and do not run dropout sensitivity tests. (Updated 2026-05-18.)
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
move_to_rootf(pathname="lrgeegfc")

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
│   ├── 04_rules/                  # renormalization, coding, naming, frontmatter, never/always
│   └── task-persistence-investigation/   # canonical home for new task-trace measures (scope reports first, code follows)
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
- `.agents/guides/02_methods/lrg-framework-guide.md` — **CANONICAL** LRG primitive (formulas, codebase mapping, verification snippets, what the Villegas papers do/don't authorize for our outlier case). Read before touching anything LRG.
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

**Task-persistence investigation (the active research question):**
- `.agents/guides/task-persistence-investigation/README.md` — canonical
  home for every multiscale measure built to investigate
  task-induced reorganization in `rest_post`. **All new task-trace
  tooling MUST land here as a mathematically rigorous scope report
  before any code is written.** README defines the required structure
  (notation → predicates → properties → caveats → pseudocode →
  visualization → connection-to-prior-tools → open-questions).
- `.agents/guides/task-persistence-investigation/archive/scalar-vi-era/2026-04-25_module-retention-landscape.md` —
  MRL: `M̄(b, ξ) ∈ [0,1]` cohort field counting `task_test` subtrees
  absent from `rest_pre` and present in `rest_post` (Jaccard match
  threshold τ).

**Historical (archive — read for context, don't cite):**
- `.agents/reports/archive/2026-04/` — pre-reset |ImCoh|² + scalar-session artifacts
- `.agents/plans/archive/` — MSC-era + |ImCoh|²-era plans
- `scripts/archive/2026-04_failed-scalar-session/` — dead-end stage scripts
