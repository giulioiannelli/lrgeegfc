---
name: pipeline-status
type: guide
era: CROSS_ERA
status: current
created: 2026-04-24
updated: 2026-04-26
pointers:
  - .agents/reports/2026-04-25_measure-ledger.csv
  - .agents/reports/2026-04-25_task-trace-audit-and-recovery.md
---

# Pipeline Status — post-ImCoh-reset era index

> **2026-04-26 ledger snapshot.** Per-measure status (current /
> superseded / archived / scope-only / stale-numeric) is now tracked
> in `.agents/reports/2026-04-25_measure-ledger.csv` (35 rows). Cohort
> returned to **n=10** with Pat_14 vendor-replacement; canonical
> reformalization (P/T/R/RA at leaf-set level) is in
> `.agents/guides/task-persistence-investigation/2026-04-25_task-trace-canonical.md`.
> Headline `(band × k)` map at n=10:
> `data/outputs/figures/section6/task_trace_band_k_n10_imcoh_abs.pdf`.

**Last updated:** 2026-04-26 (n=10 backfill + ledger + task-trace audit).
2026-04-24 (original): H2 reorganization + n=10 cohort.
**Canonical FC estimator:** `imcoh_abs` = `<|ImCoh|>_f` (signed Nolte-2004 ImCoh,
freq-resolved cache, abs-then-band-average at load time).
**Canonical writing-agent handoff:** `.agents/reports/2026-04-24_multiscale-task-trace.md`.
**Cohort:** 10 sEEG patients; cross-phase tests at n=9 (Pat_14 `task_test`
vendor-corrupt, excluded from cross-phase only); H3/H4 at n=10.

This file is the one-stop index for *what is current, what is superseded,
and what is a dead end* after the MSC → |ImCoh|² → |ImCoh| transitions
AND the 2026-04-24 H2 reorganization (H2a demoted, H2c + H2d promoted).
Nothing here is deleted — this document flags era and prevents
hallucinations from citing old numbers as if they were current.

## Quick rules

- **For the paper / writing agent**, cite
  `.agents/reports/2026-04-24_multiscale-task-trace.md` as the primary
  handoff. It contains the n=9 results for H1, H2c, H2d, H3, H4, the
  band typology, and the regeneration recipes.
- **For the VI(k) unanimity tables**, cite
  `.agents/reports/2026-04-24_h1-h4-vi-results.md`. Note H2a carries a
  stale-at-n=9 callout at the top of its section.
- **For ImCoh methodology**, cite `.agents/guides/02_methods/imcoh-guide.md`
  and `.agents/reports/2026-04-15_imcoh-reset.md`.
- **Scripts under `scripts/05_multiscale/`** accept `--fc-method` and
  default to `imcoh_abs`. New H2c/H2d scripts live under
  `scripts/01_compute/{h2c_ultrametric_drift,h2d_coactivation_persistence,
  h2_band_selectivity,h2_band_typology,rigorous_hypothesis_test,band_k_landscape}.py`.

## Era timeline

| Era | Dates | Summary |
|-----|-------|---------|
| MSC era | ≤ 2026-03 | primary FC = MSC; H1-H4 on n=5 |
| ImCoh pre-reset (\|ImCoh\|² mislabel) | 2026-03 to 2026-04-15 | primary FC = stored-as-`imcoh` but actually `⟨ImCoh²⟩_f` |
| ImCoh reset | 2026-04-15 | canonical = `imcoh_abs = ⟨\|ImCoh\|⟩_f` |
| Cohort expansion | 2026-04-22 | 5 new patients; Pat_06 re-completed |
| Cohort finalized (n=10) + H2 reorganization | 2026-04-24 | Pat_10 channel drop, Pat_13 rest_pre replaced, Pat_14 task_test excluded; H2c + H2d added as primary; H2a demoted (fails FDR at n=9) |
| H2 audit pass (H2a′ + H2e)                   | 2026-04-24+ | Directed conditional-entropy test (H2a′) on existing LRG caches: fails k-averaged FDR m=6; α passes cluster-perm at k=2-4 (p=0.050★); δ/β/γ_h trend correctly at p≈0.07 — supplementary. Split-half drift noise-floor (H2e) run on rest phases only with halved `nperseg`; results under `data/reports/imcoh_vi/h2e_split_half.{md,csv}` and `.agents/guides/02_methods/h2-metrics.md` §2a. `conditional_entropy` added to `src/lrg_eegfc/utils/metrics/vi.py`. |
| Topology & robustness audit                   | 2026-04-24++ | Added H2-RAW (Spearman on raw D), H2-FROB (Frobenius ratio), H2-TOPO (tree bipartition overlap). **All "rpost closer to task than rpre" tests fail at every metric**, confirming the trace is a residual/conditional signal not a global geometry shift. H1-topo passes 6/6 bands and H3-topo passes 5/6 — β strongest, γ_h weakest — giving a multiscale band-heterogeneity anchor independent of ρ. |ImCoh|² LRG comparison: band ranking essentially identical to |ImCoh| — FC metric is not the bottleneck. Noise-sensitivity diagnostic: ρ survives 5% FC perturbation at ρ_noise≈0.9, cross-phase ρ≈0.48 is signal not noise; but ρ is only modestly correlated with bipartition overlap (Spearman 0.31). Scripts: `h2_raw_matrix_correlation.py`, `h2_frobenius_ratio.py`, `h2_topology_directed.py`, `diag_noise_sensitivity.py`, `diag_topology_vs_rho.py`. |
| Stage 1 simple-first rerun (14 scalars)        | 2026-04-24+++ | Rerun of 14 MSC-era tree-comparison scalars (CophPearson/Spearman, BakersGamma, NormL1/L2 on ultrametric, TopK5/10/20 merges, WeightedARI coarse/fine, MeanARI, WeightedVI coarse/fine, MeanVI) under `imcoh_abs` at n=9. **None meet pre-reg criterion** (β AND γ_l q<0.05, ≥7/9 positive, θ/α null, rb(β)-rb(θ)≥0.5). Best metrics reach 2/4 criteria. **Key finding**: β-only directional signal is present (NormL1_Ult, NormL2_Ult: β p=0.033 pre-FDR, 7/9, rb=+0.69; MeanARI, WeightARI_Fine: β p=0.04-0.05, 8/9, rb=+0.64); γ_l does NOT behave as hypothesized (rb often negative). Paper's "β AND γ_l" hypothesis only half-supported at single-scalar resolution. Scripts: `scripts/01_compute/stage0b_dmax_stability.py`, `scripts/01_compute/stage1_14metrics_rerun.py`. Catalog: `data/reports/imcoh_vi/stage0a_metric_catalog.md`. Outputs: `data/reports/imcoh_vi/stage1_14metrics_{rerun.md,raw.csv,contrast.csv,stats.csv,verdict.csv}` + figures. dmax diagnostic: LRG ultrametric saturates at ≈0.99 so `h_rel` cuts are comparable without renormalization. **Next**: Stage 2 literature search for tree-similarity scalars combining heights + topology (the 14 lack this category). |
| Stage 2/3 KC + MC + wRF tree-distance tests   | 2026-04-24++++ | Implemented Kendall-Colijn (5 λ blends topology↔heights), Matching Cluster (unweighted + weighted), weighted Robinson-Foulds. `src/lrg_eegfc/utils/metrics/tree_distance.py` + `scripts/01_compute/stage3_tree_distance.py`. Also installed `dendropy 5.0.8` via pip into `lapbrain` for wRF cross-checking. **None of the 8 variants meet pre-reg**; best at 2/4 same as Stage 1. **Key diagnostic — KC(λ)**: β−θ rank-biserial separation grows from 0.09 (λ=0 topology only) to 1.46 (λ=1 heights only). At pure topology (λ=0), θ actually goes in the wrong direction (rb_θ=+0.38). At pure heights, θ is strongly negative (rb_θ=−0.82) and β is positive (rb_β=+0.64). **Scalar task-trace signal is HEIGHTS-DRIVEN, not topology-driven.** Secondary: wRF picks up an unexpected δ-band signal (rb=+0.82, 7/9 positive, p=0.014 pre-FDR, q=0.085). Stage 2 menu: `.agents/reports/2026-04-24_stage2-literature-menu.md`. Stage 3 outputs: `data/reports/imcoh_vi/stage3_tree_distance.{md,_stats.csv,_verdict.csv,_raw.csv,_contrast.csv}` + figures incl. `kc_lambda.pdf`. **Consequence for Stage 4**: the originally-planned topology-based band×scale map is unlikely to rescue β (topology is null); Stage 4 must be redesigned — either scale-resolved heights (H2c-refinement) or hybrid heights+topology per cell. User decision pending. |



---

## Current (`imcoh_abs`, post-reset)

| artefact | role |
|---|---|
| `data/cache/imcoh/` | raw signed ImCoh, freq-resolved per (patient, phase, band). Canonical. |
| `data/cache/imcoh_lrg/Pat_*/{band}_{phase}_lrg_imcoh-abs.npz` | LRG results on \|ImCoh\|. Complete 120/120 for 5 × 6 × 4. |
| `scripts/01_compute/compute_imcoh_fc.py` | recomputes the freq-resolved signed ImCoh cache |
| `scripts/01_compute/compute_imcoh_lrg.py` | runs LRG on `imcoh_abs` (transform-tagged output) |
| `scripts/01_compute/compute_imcoh_vi.py` | VI(k) profiles + H1/H2a/H2b/H3 contrasts; coverage guard + `--dry-run` |
| `scripts/01_compute/report_h1h4_vi.py` | emits `2026-04-24_h1-h4-vi-results.md` from the CSVs |
| `scripts/05_multiscale/multiscale_all_hypotheses.py` | H1-H4 per k; **now `--fc-method`-aware**, default `imcoh_abs` |
| `scripts/05_multiscale/continuous_multiscale_h2.py` | H2 continuous over threshold height; `--fc-method`-aware |
| `scripts/05_multiscale/definitive_multiscale_h2.py` | H2 unanimity heatmaps; `--fc-method`-aware |
| `scripts/05_multiscale/continuous_all_pairs.py` | all-pairs VI(h); `--fc-method`-aware |
| `scripts/05_multiscale/analyze_h4_gradient.py` | H4 Kendall W / Spearman; reads imcoh CSV when `--fc-method imcoh_*` |
| `src/lrg_eegfc/workflow/lrg.py::compute_lrg_analysis` | guards: signed `imcoh` raises ValueError |
| `src/lrg_eegfc/config/paths.py::lrg_cache_for`, `lrg_filename` | routing + transform-tagged filenames |
| `src/lrg_eegfc/workflow/fc.py::_load_imcoh` | per-frequency abs/sq then band-average (Jensen's inequality) |
| `.agents/reports/2026-04-15_imcoh-reset.md` | post-reset consolidation doc |
| `.agents/reports/2026-04-18_epileptic-imcoh-final.md` | post-reset epileptic-node clustering (paper-ready) |
| `.agents/reports/2026-04-24_h1-h4-vi-results.md` | **this run's authoritative H1-H4 results** |
| `.agents/guides/02_methods/imcoh-guide.md` | method description |
| `.agents/guides/02_methods/probe-bias-guide.md` | probe-geometry caveat (still applies) |
| `data/reports/imcoh_vi/vi_raw_profiles.csv` | per-(patient, band, phase-pair, k) VI |
| `data/reports/imcoh_vi/hypothesis_contrasts.csv` | per-(patient, band, k, hypothesis) signed contrast |

## Deleted root-level scripts (2026-04-22 CLI consolidation)

15 one-off figure scripts at the `scripts/` root level were removed
when commits `ff2916b` (unified `lrg-eegfc` CLI with 36 subcommands) and
`fc2bed1` (plot CLI `--format`/`--dpi` options) shipped. All functionality
is reachable via the CLI. **Do not restore these files** — they duplicate
CLI capabilities and would drift out of sync.

| Deleted | Successor in `lrg-eegfc` CLI |
|---------|------------------------------|
| `c1_psi_overlay.py` | `lrg-eegfc plot lrg --plot-type entropy` (with Ψ overlay) |
| `fig2_surrogate_comparison.py` | removed, diagnostic only |
| `gen_alluvial_grid.py`, `gen_sankey_grid.py` | `lrg-eegfc plot ...` uses `create_sankey_diagram` |
| `gen_corr_vs_msc_figures.py` | `lrg-eegfc plot corr` / `lrg-eegfc plot msc` via `--fc-method` |
| `gen_fc_figures_fast.py` | superseded by `lrg-eegfc plot` commands |
| `gen_metastable_grid.py` | visualisation folded into `lrg-eegfc plot ... --plot-type …` |
| `gen_metric_concordance.py`, `gen_metric_crosspatient.py` | `lrg-eegfc` metric-exploration commands |
| `gen_mslcd_cross_figures.py`, `gen_mslcd_figures.py` | `lrg-eegfc plot` with MSLCD inputs |
| `gen_phase_reorg_figures.py` | `lrg-eegfc plot` uses `plot_phase_reorganization` |
| `gen_ultrametric_figures.py` | `lrg-eegfc plot lrg --plot-type ultrametric` |
| `msc_surrogate_scaling.py`, `precompute_fc_figures.py` | superseded by the cache-first CLI flow |

If you truly need one of these back, recover from git (`git show
fc2bed1~1:scripts/<name>.py`); do not reintroduce the root-level scripts
into the working tree — the CLI is the canonical interface now.

## Superseded (|ImCoh|² era, pre-2026-04-15 reset)

Do **not** cite numeric values from these. They computed `<ImCoh²>_f`
instead of `<|ImCoh|>_f`. Kept for audit/provenance only.

| artefact | replaced by |
|---|---|
| `data/cache/_archive_absimcoh_sq/` | `data/cache/imcoh/` |
| `data/cache/_archive_absimcoh_sq_lrg/` | `data/cache/imcoh_lrg/` |
| `.agents/reports/2026-04-15_imcoh-verification.md` | `2026-04-24_h1-h4-vi-results.md` |
| `.agents/reports/2026-04-24_imcoh-results-for-writing.md` | `2026-04-24_h1-h4-vi-results.md` |
| `.agents/reports/2026-04-15_imcoh-process.md` | `2026-04-15_imcoh-reset.md` |
| `.agents/reports/2026-04-15_imcoh-pat02-controls.md` | superseded by cohort-level analysis in the post-reset reports |
| `.agents/reports/2026-04-15_imcoh-recovery-plan.md` | historical, transition doc only |

## Superseded (MSC era)

MSC suffers from common-reference / same-probe bias (CLAUDE.md invariant 4).
Outputs remain for cross-method comparison but are not the canonical
answer to H1-H4 any more.

| artefact | status |
|---|---|
| `data/cache/msc/`, `data/cache/msc_dev/` | MSC FC caches |
| `data/cache/corr/` | Pearson correlation cache |
| `data/cache/lrg/`, `data/cache/lrg_crema/` | LRG results on MSC/corr FC |
| `data/outputs/figures/metric_exploration/partition_multiscale/results.csv` | wide MSC-era aggregate; `analyze_h4_gradient.py` reads this when `--fc-method msc` |
| `ipynb/02_fc_msc/*` (5 notebooks) | MSC workflows |
| `ipynb/06_presentation_figures/*` (7 notebooks) | MSC-era figure drafts |
| `.agents/reports/2026-02-26_report-figures.md` | pre-reset figure manifest |
| `.agents/reports/2026-04-13_session-reorganize-scripts.md` | pre-reset session snapshot |
| `scripts/01_compute/batch_lrg_crema.py` | MSC variant (lrg_crema cache) |
| `scripts/01_compute/compute_bipolar_lrg.py` | MSC + bipolar (see "dead ends") |
| `scripts/07_figures/hypotheses_final_figures.py` | hardcodes `fc_method="msc"`; **figures chat will rewire** |
| `scripts/07_figures/consolidate_report_figures.py` | same |
| `scripts/07_figures/fig_S5_*.py`, `gen_affinity_all_patients.py` | same |

## Experimental dead ends (documented, not used)

| artefact | reason |
|---|---|
| `data/cache/bipolar/` | bipolar re-referencing experiment; abandoned |
| `data/cache/rescaled/` | rescaling experiment; abandoned |

## Unknown — verify before reuse

| artefact | concern |
|---|---|
| `data/outputs/tables/mslcd_diagnostics_master.csv` | MSLCD diagnostics; FC provenance not tagged |
| `data/outputs/figures/metric_exploration/*/results.csv` (non-partition_multiscale) | mixed provenance |
| `ipynb/03_lrg/*`, `ipynb/04_reorganization/*`, `ipynb/05_figures/*` | may inherit MSC caches silently |

## Infrastructure rewiring done in this session

The following pre-existing MSC-hardcoded scripts were adapted to accept
`--fc-method` (default `imcoh_abs`) and route cache + output directory
accordingly. **No computation logic was changed** — the same multiscale
machinery now serves both FC families through a flag.

- `scripts/05_multiscale/multiscale_all_hypotheses.py`
- `scripts/05_multiscale/continuous_multiscale_h2.py`
- `scripts/05_multiscale/definitive_multiscale_h2.py`
- `scripts/05_multiscale/continuous_all_pairs.py`
- `scripts/05_multiscale/analyze_h4_gradient.py` (reads `vi_raw_profiles.csv`
  and aggregates to the wide schema when `--fc-method imcoh_*`)

Figure-only scripts in `scripts/07_figures/` still hardcode MSC; they will
be addressed in the dedicated figures chat.

## Top hallucination hazards

Prioritized watch-list — items that look current but aren't:

1. `.agents/reports/2026-04-15_imcoh-verification.md` — cell counts look legitimate but are |ImCoh|².
2. `.agents/reports/2026-04-24_imcoh-results-for-writing.md` — "MSC vs ImCoh" table uses |ImCoh|².
3. `scripts/07_figures/hypotheses_final_figures.py` — defaults to MSC; any figure produced is MSC-era.
4. `ipynb/02_fc_msc/*` — will happily run against the still-present `data/cache/msc/` files.
5. `data/cache/lrg/`, `data/cache/lrg_crema/` — LRG from MSC; not the current canonical.
6. `data/outputs/figures/metric_exploration/partition_multiscale/results.csv` — wide H4 input, MSC.
7. Any citation of "the beta gap" numbers without a source tag — pre- and post-reset values differ.
