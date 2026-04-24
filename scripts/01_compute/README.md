# scripts/01_compute

**Compute-stage scripts: hypothesis tests, diagnostics, batch jobs,
embedded-figure producers. Split into four subdirs by role.**

## Layout

- `hypothesis_tests/` — H2-family cross-patient tests
  (`h2c_ultrametric_drift`, `h2d_coactivation_persistence`,
  `h2_partition_multiscale`, etc.) + `rigorous_hypothesis_test`.
- `diagnostics/` — diagnostic scripts (`diag_*`, `band_k_landscape`).
- `batch/` — FC + LRG computation jobs (`compute_*`, `batch_*`,
  `clean_correlation_matrices`, `report_h1h4_vi`).
- `figures_embedded/` — `fig_*` scripts that produce publication
  figures directly from cached compute outputs.

## Imports

All scripts import statistical helpers from
`lrg_eegfc.utils.metrics.hypothesis` (post-3.8 library elevation).
No more `from _shared import ...`.

## Before adding a new script

1. Check `src/lrg_eegfc/` for existing helpers (grep first).
2. Pick the right subdir (if none fits, think twice before adding a
   5th subdir — probably it belongs in the library).
3. If a helper is used by ≥2 scripts, promote to `lrg_eegfc.utils.*`
   in the same commit.

See `.agents/guides/04_rules/coding-rules.md` for the full ruleset.
