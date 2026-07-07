---
name: stage-02-fc-msc-corr
type: plan
era: MSC
status: superseded
created: 2026-01-10
updated: 2026-04-24
pointers: []
---

# Plan (2026-01-10): Stage 02 - FC matrices (MSC + correlation)

## Goals
- Build MSC FC matrices as the main connectivity signal.
- Compute correlation FC matrices only for comparisons.
- Keep everything cache-first and patient-agnostic.

## Inputs
- Outputs from Stage 01 (dataset inventory)
- `lrg_eegfc.workflow.msc`, `lrg_eegfc.workflow.corr`
- Existing notebooks and legacy notebooks for parameter references

## Current status
- Dev notebooks exist for single-patient FC generation.
- `src/compute_corr_matrices.py` and `src/compute_msc_matrices.py` accept
  `--bands`, `--phases`, and `--filter-time`.
- No surrogate scaling test results yet.
- Sanity check showed MSC cache names include `sparsify-*` and `nperseg-*`
  (e.g., `*_msc_sparsify-none_nperseg-1024.npy`); loader defaults must match.
 - Parameter summary draft created (`2026-01-10_fc_params.md`).
 - MSC method guide drafted (`.agents/guides/msc-method-guide.md`).
 - Stage 02 legacy notebooks consolidated into updated single-case notebooks.

## Tasks
1) Parameter consolidation
- Inspect MSC parameter choices in legacy notebooks (e.g. UTILS-COHERENCE).
- Record defaults for development and final runs.
- Extract correlation comparison settings from notebooks (abs/threshold/cleaned).
- Capture default `nperseg`, `noverlap`, and surrogate policy per band.
- Record dev vs final settings in `.agents/plans/archive/setup-era/2026-01-10_fc_params.md`.

2) Surrogate scaling tests
- Create `tests/` (if missing) and add a timing/stability harness for MSC.
- Sweep surrogate counts (e.g. 0, 25, 50, 100, 200, 500) and record runtime.
- Track network stability vs surrogate count (e.g. Frobenius norm diff,
  Pearson/Spearman correlation on upper triangle, degree distribution change).
- Identify a saturation point to justify final surrogate choice.
- Implementation targets:
  - `tests/test_msc_surrogates_scaling.py` (timing + metrics, opt-in)
  - `scripts/msc_surrogate_scaling.py` (interactive CLI sweep)
- Use the saturation test to justify per-window surrogate counts in Stage S.

3) Single-patient dev notebooks
- Maintain single-patient notebooks as templates for quick runs.
- Keep patient name as a variable; do not hard-code.
- Ensure dev defaults target fast runs (filter_time, low surrogates).
- Target notebooks:
  - `ipynb/02_fc_msc/01_build_corr_networks_singlepat.ipynb`
  - `ipynb/02_fc_msc/02_build_msc_networks_singlepat.ipynb`
  - `ipynb/02_fc_msc/04_msc_dense_vs_validated_singlepat.ipynb`
  - `ipynb/02_fc_msc/05_msc_vs_corr_singlepat.ipynb`

4) Batch-ready scripts
- Ensure `src/compute_msc_matrices.py` and `src/compute_corr_matrices.py`
  accept bands/phases/filter_time and remain cache-first.
- Document how to switch between dev (single patient) and full batch.
- Ensure scripts print cache paths and parameter choices to stdout for traceability.

5) Runtime guardrails (dev only)
- Target sub-60s runtime per dev run.
- Use `filter_time` and low surrogate counts in dev.
- Do not add hard timeouts; use small inputs instead.

6) MSC methodology section (documentation)
- Create `.agents/guides/msc-method-guide.md` with:
  - Rationale: MSC variance depends on Welch segment count, not MP.
  - Window-length rule: min cycles at lowest band frequency.
  - Default `nperseg`/overlap and dev vs full settings.
  - Surrogate validation policy and saturation test reference.
  - Clarify that MP noise applies to correlation matrices, not MSC.
  - Note: final time-window analysis uses per-window surrogates (Stage S).
  - Note the 25% overlap default for time-window analyses.

## Compute vs visualize
- Compute: `src/compute_corr_matrices.py`, `src/compute_msc_matrices.py`,
  `compute_corr_matrix`, `compute_msc_matrix` write caches.
- Visualize: notebooks/scripts must call `load_*` first and should not compute
  unless caches are missing (dev-only exceptions).
- Cache roots are fixed under `data/{corr_cache,msc_cache}/` for full runs.
- Dev runs with `filter_time` auto-redirect to `data/{corr_cache_dev,msc_cache_dev}/`.

## Deliverables
- Parameter summary file for MSC and correlation.
- Surrogate scaling report under `tests/` (runtime vs stability).
- Single-patient notebooks for quick validation.
- Batch scripts ready for external machine execution.
- MSC methodology guide under `.agents/guides/`.

## Suggested files
- `.agents/plans/archive/setup-era/2026-01-10_fc_params.md`
- `tests/test_msc_surrogates_scaling.py`
- `.agents/guides/msc-method-guide.md`

## Exit criteria
- MSC matrices cached under `data/msc_cache/Pat_XX/` for dev scope.
- Correlation matrices cached under `data/corr_cache/Pat_XX/` for comparison.
