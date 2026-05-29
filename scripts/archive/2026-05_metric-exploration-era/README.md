---
name: 2026-05_metric-exploration-era
type: archive-readme
era: MSC × COHORT_N5 / COHORT_N6 (exploratory)
status: archived
created: 2026-05-28
---

# 2026-05 metric-exploration archive

This folder consolidates exploratory script families that were retired
2026-05-28 during the pre-preprint cleanup pass. None of these scripts'
outputs are cited from the current `.agents/preprint/` hub or from
post-2026-04-25 reports under `.agents/reports/` (live tier). The
scientific value of the explored ideas is preserved here for archaeology
+ recovery; the helpers that turned out to be reusable were promoted
into `src/lrg_eegfc/` BEFORE the archive move (see the cleanup TLDR at
`.agents/reports/2026-05-28_pre-preprint-cleanup.md` for the rescue
list).

## What's here

- `03_analysis/` (18 files, ~10K LOC) — MSC-era + n=6 4-phase cohort
  exploratory analyses: VI(k) vs RF distance (`tree_vs_vi_analysis`),
  segregation / modularity / bridge / merge-rank metric families,
  partition entropy, ultrametric distributions, community-size
  scaling, cophenetic fit, entropy fingerprinting. 16/18 hardcode
  `fc_method="msc"` and a 6-patient list (`Pat_02/03/05/06/07/08`);
  the other 2 are exploratory benchmarks with no cohort framing.
- `04_reorganization/` (13 files) — MSC-era reorganization metrics
  (logCosine D(τ), mean ARI, multiscale ARI, normalized reorg, etc.).
  Knowledge condensed into the preprint-locked ρ_split^coph probe;
  scripts themselves not used.
  - `script_slanzarv.py` is a stub launcher for a non-existent
    compute path; truly dead.
  - `plot_reorganization_msari_v2.py` carried a Pat_08 outlier
    diagnostic panel that violates the 2026-05-18 "no per-patient
    outlier framing" rule.
- `05_multiscale_legacy/` (4 files) — `05_multiscale/` scripts that
  read the legacy `data/figures/metric_exploration/partition_multiscale/results.csv`
  (MSC era): `multiscale_community_flow`, `multiscale_affinity_analysis`,
  `multiscale_h2_profile`, `explore_raw_vi`. The other 7 files in
  `scripts/05_multiscale/` remain live (use `imcoh_abs`).
- `06_metric_sweep/` (10 files, ~8K LOC) — systematic dead-end sweeps:
  17 cophenetic-distance metrics, 4 alternative families, 6 metric
  systems, ultrametric variants. None achieved consensus; preserved
  for documentation of the search.
- `wp0/` (7 files + `_common.py`) — metric exploration workpackage
  (4 tasks: inventory, compute, analysis, multiscale). MSC + n=5.
- `wp1/` (7 files) — co-classification workpackage variants
  (`wp1_coclassification`, `wp1_spatial`, `wp1_scalar_vi`, etc.).
  MSC + n=5.
- `audit_round2_section5_v2.py` + `audit_round2_section5_v2_fixes.py`
  — superseded by the live `audit_round3_section5_redo.py` (Phase 4-B
  will split that one into per-figure scripts).
- `gen_cross_phase_comparison_v2.py` — superseded; v1 stays live.
- `gen_brain_connectome_multiscale.py` — MSC variant of
  `gen_brain_connectome_imcoh.py`; the imcoh-named version stays live
  (will be parametrized in Phase 4-B if a unifying flag is desired).

## What was rescued before archiving

The library helpers re-rolled inside these scripts (4 of them) were
promoted in Phase 4-A:

- `surrogate_p_value` → `lrg_eegfc.utils.metrics.hypothesis`.
- `loo_sensitivity` → same module.
- `cophenet_matrix(Z)` → `lrg_eegfc.utils.metrics.tree`.
- `PatientMasks` + `build_epi_masks` → `lrg_eegfc.utils.io.patient`.

`compute_vi` was already in the library at
`lrg_eegfc.utils.metrics.vi`. The `tree_distance_bootstrap` named by
the initial audit doesn't actually exist as a function anywhere —
the bootstrap patterns in `audit_57_kc_topology_anatomy.py`,
`audit_56_anatomy_distribution.py`, etc. are bespoke per script and
deliberately not abstracted.

## Recovery

`git log --follow <path>` works on the archived paths because the
moves used `git mv`. To restore any file to its original location
(for re-analysis or comparison), simply `git mv` it back.
