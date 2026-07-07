---
name: epi-wm-stratification-rhosym-migration
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-06
supersedes_claims:
  - "epi/WM pair-class tissue numbers (R1.4/R1.6) are ρ_split-era (audit_77/83/85)"
pointers:
  - scripts/01_compute/audit/audit_156_epi_wm_stratified_rhosym.py
  - data/audit/epi_wm_stratified_rhosym/pairclass_rhosym_cohort.csv
---

# Epi/WM tissue stratification — ρ_sym migration (R1.4 / R1.6 numbers)

## Head

The epi/WM pair-class tissue results that R1.4 (β rides healthy cortex) and R1.6
(α recruits the seizure core) cite are now re-estimated under **ρ_sym** (audit_156),
so no ρ_split number lands in the paper. **Every load-bearing verdict holds
like-for-like.** The ONE refinement: **wm↔wm β no longer independently clears the
matched-strength null** (ρ_split Wilcoxon p=0.007 → ρ_sym p=0.065), so the
"WM is strength-clean/separated" sub-claim is dropped — WM↔WM β is simply a
below-average contributor (still depleted vs random pairs, p_pair=1.0). The
gray-dominant headline is unchanged.

## Method (faithful re-run, no old audit edited)

`audit_156_epi_wm_stratified_rhosym.py` ports the pair-class statistic of audit_77
(epi matched-strength), audit_83 (WM matched-strength) and audit_85 --mode pairclass
(pair-decimation null) to the symmetric per-class concordance

    rho_sym(class) = 1/2[ spearman((D_task-D_preA)[pk], (D_post-D_preB)[pk])
                        + spearman((D_task-D_preB)[pk], (D_post-D_preA)[pk]) ]

applied identically to the observed trace, every cached matched-strength surrogate
(canonical seed-20260511 full-graph ensemble, **reused — no regeneration**), and
every random-pair draw. arm1 (single-arm) reproduces the locked ρ_split CSVs
**bit-for-bit** (verified: all pair classes, both nulls, d=0). Spearman kernel is
**numba-parallel, bit-identical to scipy** (max|nb−sp|=0.0); full-grid CSVs
byte-identical scipy-vs-numba. Runtime 60 cells ≈ 3.5 s/cell.

## Verdicts under ρ_sym (cohort, cophenetic) — cited cells

| cell | ρ_sym obs | matched-strength Wilcoxon p | pair-decimation p_pair | verdict | ρ_split was |
|---|---|---|---|---|---|
| **gray↔gray β** | +0.235 | **0.042** (clears) | **0.000** | **CONCENTRATED + clears MS** | 0.032 / 0.000 concentrated |
| wm↔wm β | +0.103 | 0.065 (n.s.) | 1.000 | **depleted**, MS n.s. | 0.007 sep / 1.000 depleted |
| **epi↔epi β** | +0.081 | 0.312 (n.s.) | 0.926 | **generic/trending-depleted → core does NOT carry β** | 0.161 / 0.945 generic |
| **epi↔epi α** | +0.410 | **0.005** (clears) | **0.000** | **CONCENTRATED + clears MS → α RECRUITS core** | 0.003 / 0.000 concentrated |
| wm↔wm α | +0.170 | 0.010 | 0.000 | concentrated | 0.032 / 0.000 concentrated |
| epi cross α | +0.130 | 0.003 | 0.016 | concentrated (corroborates α-recruits) | 0.002 / 0.000 concentrated |

Read the matched-strength gate off the **paired one-sided Wilcoxon p** (the test IS
the gate); ms_verdict="intermediate" for all cells only because the strict
`separated` label additionally requires n_above≥8 (descriptive, never the gate).

## What R1.4 / R1.6 should say

- **R1.4 (β):** gray↔gray cortical coupling CARRIES the β trace — concentrated vs
  random pairs (p_pair<0.001) AND clears the strength-matched null (p=0.042). WM↔WM
  and the seizure core are below-average: WM↔WM depleted (p_pair=1.0, MS n.s.); the
  seizure core carries no more than random pairs and does not clear MS (epi_epi
  p_pair=0.93, MS p=0.31) → "core spared" is DIRECTIONAL. DROP "WM strength-clean".
- **R1.6 (α):** epi↔epi α RECRUITS the core — concentrated (p_pair<0.001) AND
  strengthens under the strength-matched null (p=0.005). Opposite to β. Both nulls agree.

Framework purity: cophenetic only; raw = baseline. Node-count confound handled by the
pair-decimation null (no node removal). estimator-invariant except the one wm_wm
refinement noted.
