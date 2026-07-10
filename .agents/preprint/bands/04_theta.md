---
name: preprint-theta-band
era: IMCOH_ABS_COHORT_N10
status: current
kind: preprint-result-report
band: theta
range_hz: [4, 8]
cohort: [Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15]
n_patients: 10
fc_method: imcoh_abs
tau: 1/lambda_max
canonical_lrg_object: "D_coph = cophenet(UPGMA(D(tau_max)))"
canonical_subspace_object: "U_k = span{phi_2, ..., phi_{k+1}}"
verdict_tag: "no trace"
verdict_source: VERDICT_LEDGER.md (locked 2026-05-18, revised 2026-05-19)
verdict_layers:
  substrate_rank: borderline_positive_direction (raw FC ρ_split^raw 7/10 p=0.138, ratio 20.2× — cohort-median direction is positive but not significant)
  rho_split_coph: no_trace_anti_direction (C3 fails decisively: paired Wilcoxon p=0.722, n_above_surrogate 2/10, ratio −10.2× — cohort-median ρ is in the anti-trace direction)
  grassmann: no_trace (cluster-extent permutation p=0.0995, audit_70; 5-cell observed run within null distribution at 95th percentile = 6.0)
  cophenet_sign_inversion: cophenet_step_flips_substrate_sign (raw FC +0.117 / raw D +0.044 → cophenet −0.040; multiscale aggregation reverses cohort-median direction at θ)
  anatomy: not_yet_audited (no θ-specific Desikan-Killiany enrichment run; not relevant under no-trace verdict)
sources:
  - .agents/preprint/locked/CONTROLS.md (locked 5-control battery)
  - .agents/preprint/locked/VERDICT_LEDGER.md (locked verdicts)
  - .agents/preprint/directives/archive/2026-05/methods_revision_2026-05-18_cophenet.md (binding methods directive)
  - data/audit/ctm_triangle/cohort_summary.csv (θ row: rho_split_median, n_above_drift, rho_xprobe, Wilcoxon)
  - data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv (θ row: paired Wilcoxon vs surrogate)
  - data/audit/grassmann_cluster_extent/cohort_summary.csv (θ cluster-extent p=0.0995, audit_70)
  - data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv (θ per-k T_G + surr stats, audit_66)
  - data/audit/grassmann_epi_exclusion/sensitivity.csv (θ per-k under C5 epi-X, audit_67)
  - data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv (θ substrate row)
  - data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv (θ raw D sensitivity)
revision_history:
  - 2026-05-19: initial brief produced from VERDICT_LEDGER.md lockdown
---

# θ band (4–8 Hz) — preprint result report

## Head

θ is the **cleanest negative** in the panel: **no trace on either probe**, and the cohort-median cophenet `ρ_split^coph` sits in the *anti*-trace direction (−0.040, ratio −10.2× to surrogate, paired Wilcoxon p = 0.722, 2/10 patients above their own surrogate). All four primary controls fail on the cophenet probe (C1 p = 0.687, C2 p = 0.278, C3 p = 0.722, C4 negative-sign at 3/10). The Grassmann subspace probe fails the cluster-extent permutation gate at cluster_p = 0.0995 (5-cell observed longest run within the empirical null distribution; null 95th percentile = 6.0). The substrate-to-cophenet trajectory at θ shows a **sign inversion** that the other bands do not exhibit: raw FC cohort median is positive +0.117 (7/10 patients), raw D(τ_max) is positive +0.044 (6/10), but the cophenet step inverts the sign to −0.040 (3/10 trace-direction). Verdict from `locked/VERDICT_LEDGER.md`: **`no trace`** — useful as the band-specificity benchmark against which β, α, γ_l, and δ traces are measured.

## Headline three-layer cohort table (θ-band row, copied from `locked/VERDICT_LEDGER.md`)

```
Layer                                       θ (full n=10)
-------------------------------------      ---------------------
Raw FC          ρ_split^raw  (raw_fc_ms)   7/10 p=.138  ratio 20.2×  ← borderline +
Raw D(τ_max)    ρ_split      (preprint_05) 6/10 p=.097  ratio 11.2×  ← borderline +
Cophenet D_coph ρ_split^coph (audit_63)    2/10 p=.722  ratio −10.2× ← FAILS, sign-inverted
Grassmann       d_G(k)       (audit_70)    5-cell run, cluster_p = 0.0995  ← no trace
```

Sources: `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv`, `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv`, `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv`, `data/audit/grassmann_cluster_extent/cohort_summary.csv`.

**Key reading**: θ exhibits **cohort-median sign inversion across the substrate-to-cophenet layers**. Raw FC and raw D both find weak positive cohort-median ρ (with raw FC reaching 7/10 patients but only at p = 0.138, not significant). The cophenet aggregation step flips the cohort-median sign to negative — only 2/10 patients are in the + direction after the multiscale step. This is the inverse of γ_l (which keeps the + sign but loses cohort consistency) and the inverse of β (which amplifies cohort consistency). At θ, the cophenetic step is actively *de-coheering* the per-pair signal: what little structure existed at the substrate level does not survive multiscale dendrogram aggregation. This is consistent with θ-band coupling carrying no coherent task-related per-pair reorganization at the LRG layer.

## 1. Scientific claim

**Cohort-level question.** Does the θ-band (4–8 Hz) post-task resting state at the LRG-CTM layer carry task-related reorganization on either of the two locked probes?

**Refined biological claim (null direction).** θ coupling **does not** reorganize during the task in a way that is retained in either the per-pair cophenet multiscale geometry or the leading-mode Laplacian subspace. The substrate-level raw FC and raw D(τ_max) borderline-positive cohort-median ρ values (+0.117 and +0.044 respectively) reflect ambient cohort variability rather than task-specific reorganization, as confirmed by the matched-strength C3 gate failing at all three layers (raw FC p = 0.138, raw D p = 0.097, cophenet p = 0.722, Grassmann cluster_p = 0.0995).

**Falsification budget.** This `no trace` verdict would flip to `weak trace` or `strong trace` if any of: (a) cophenet C3 paired Wilcoxon p < 0.05 (currently 0.722, far from); (b) Grassmann cluster-extent permutation p < 0.05 (currently 0.0995, just outside); (c) a re-run with corrected sample preparation produces a positive cohort-median cophenet ρ at 6+/10 patients. None of these are expected given the cohort-wide pattern and the size of the negative-direction cohort signal.

## 2. Cohort, substrate, and library entry points

Same as α, β, γ_l: n=10 patients (Pat_02/03/05/06/07/08/10/13/14/15); `imcoh_abs` band-averaged magnitude of imaginary coherence; θ band `B = [4, 8] Hz`. Pat_03 acquired at 1024 Hz; handled at the config layer (`nperseg_for_fs(fs)`, `FS_OVERRIDES`) — full cohort member at analysis layer (policy updated 2026-05-18, no dropout sensitivity test). Pat_10 task rows [53, 54, 55] dropped at load. Pat_14 vendor-replaced 2026-04-25.

```python
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
W = load_fc_matrix(patient="Pat_02", phase="rest_pre", band="theta", fc_method="imcoh_abs")
lrg = load_lrg_result(patient="Pat_02", phase="rest_pre", band="theta", fc_method="imcoh_abs")
# Probes documented for completeness; both return no trace at θ
```

## 3. The two LRG probes — methodology, results, and provenance

### 3.1 Per-pair multiscale correlation on `D_coph` — **no trace, sign-inverted**

#### Critical preamble (5-point)
1. **Claim.** θ `ρ_split^coph > 0` at cohort level.
2. **Null (C3, mandatory).** Strength-preserving 4-cycle ±δ rewiring, R=200, seed 20260511. Paired Wilcoxon one-sided greater across patients: obs > surrogate.
3. **Strongest plausible alternative.** Within-session drift driving a per-pair distance shift unrelated to task; same-strength edge re-arrangement producing a small positive cohort-median by chance.
4. **Does C3 cover (3).** Yes for the strength-rearrangement alternative; C2 drift-floor addresses within-session drift. Both fail at θ.
5. **Falsification of the verdict (no trace).** Would require C3 paired Wilcoxon p < 0.05. Currently p = 0.722, far from any reasonable threshold.

#### Results (θ, full cohort)

| Statistic | Value | Source |
|---|---|---|
| Cohort-median `ρ_split^coph` | **−0.0492** | `ctm_triangle/cohort_summary.csv` theta row |
| n_trace_split (`ρ > 0`) | **3/10** | same |
| n_above_drift | 6/10 | same |
| C1 Wilcoxon p (`ρ_split > 0`) | **0.6875** | same col `wilcoxon_split_gt_0_p` — ✗ |
| C2 Wilcoxon p (`ρ_split > ρ_drift`) | 0.2783 | same col `wilcoxon_split_gt_drift_p` — ✗ |
| C3 matched-strength: obs_median | **−0.0399** | `matched_strength_surrogate_split_baseline/cohort_summary.csv` theta |
| C3 surr median (median over surrogate medians) | +0.0039 | same |
| C3 ratio obs/surr | **−10.2×** (sign-inverted) | derived |
| C3 paired Wilcoxon z / p | 22.0 / **0.7217** | same — ✗✗ |
| C3 n_above_surrogate (per-patient) | **2/10** | same |
| C4 cross-probe rho_xprobe | −0.0483 | `ctm_triangle/cohort_summary.csv` theta |
| C4 n_trace_xprobe | 3/10 (−sign matches rho_split −0.049) | same — sign-match only; cohort < 6/10 → ✗ |
| C5 epi-X cophenet | **not run** for θ (audit_68 covered α only) | — |

**Every primary control fails at θ on the cophenet probe.** The cohort-median is in the anti-trace direction, only 2 patients sit above their own matched-strength surrogate, and the cross-probe restriction confirms the negative-sign cohort consensus at 3/10. Note that C4 sign-match passes formally (−sign in both `rho_split` and `rho_xprobe`), but the cohort count 3/10 < 6/10 fails the C4 gate.

#### Reading
θ cophenet has **no trace**, and the cohort signal is **actively in the anti-trace direction**. This is the strongest negative reading in the panel. The substrate-to-cophenet layer shows a sign inversion (raw FC and raw D both positive cohort-median; cophenet negative) — the multiscale dendrogram aggregation is not just losing signal, it is actively reorienting the cohort-median direction. Mechanistically this is consistent with the substrate-level + signal being a low-coherence ambient pattern that re-projects into the dendrogram structure unpredictably across patients; the cophenet step thus acts as a *low-pass filter* on θ per-pair structure, with the residual being noise-dominated.

#### Cache + script provenance
- Full cohort: `audit_63_split_baseline_surrogate.py` → `data/audit/matched_strength_surrogate_split_baseline/`
- Within-baseline triangle: `audit_*_ctm_triangle.py` → `data/audit/ctm_triangle/`
- Object built from `lrg.ultrametric_matrix` field in `data/cache/imcoh_lrg/Pat_NN/theta_{phase}_lrg_imcoh-abs.npz`

### 3.2 Global subspace rotation: Grassmann chordal distance `d_G(k)` — **no trace**

#### Critical preamble (5-point)
1. **Claim.** `T_G(k) < 0` at cohort level for some contiguous window of `k`.
2. **Null (C3, audit_70 cluster-extent permutation).** R=200 matched-strength surrogates; cluster_p = (1 + #(null_LR ≥ obs_LR)) / (R + 1).
3. **Strongest plausible alternative.** Same-strength edge rewiring producing a similar 5-cell contiguous cluster by chance.
4. **Does C3 cover (3).** Yes — the cluster-extent null distribution is built from R=200 phantom-surrogate tests on the same k-grid; the 5-cell observation is compared to the empirical distribution of longest-run lengths.
5. **Falsification of the verdict (no trace).** Would require cluster_p < 0.05. Currently 0.0995, just outside.

#### Results (θ)

| Statistic | Value | Source |
|---|---|---|
| Observed longest contiguous-significant run | **5 cells** | `grassmann_cluster_extent/cohort_summary.csv` theta |
| Cluster-extent null mean LR | 2.13 | same |
| Cluster-extent null median LR | 2.0 | same |
| Cluster-extent null 95th LR | **6.0** | same |
| Cluster-extent null max LR | 17 | same |
| Observed cluster mass (Σ −log10p) | 7.11 | same |
| Null mean cluster mass | 3.95 | same |
| Null 95th cluster mass | 11.84 | same |
| **cluster_p_longest_run** | **0.0995** | same — ✗ |
| **cluster_p_cluster_mass** | 0.1443 | same — ✗ |
| C5 epi-X (audit_67): summary | no improvement; θ remains within null distribution under epi-X | `grassmann_epi_exclusion/sensitivity.csv` theta rows |

The 5-cell observed longest run sits **below the null 95th percentile** (6.0) and only marginally above the null mean (2.13). Both cluster-p statistics fail the 0.05 gate. The cluster-mass companion (0.144) is further from significance than the longest-run statistic (0.099). Under C5 epi-X (audit_67), no improvement in θ's Grassmann signature: per-k cohort Wilcoxon p values remain in the 0.3–0.9 range across most of the k-grid (single isolated cells at k=7 with epi-X p=0.032 are not part of a contiguous window).

#### Reading
θ Grassmann has **no trace**. The 5-cell run is consistent with the chance-cluster distribution from matched-strength rewiring. The cluster-mass statistic (0.144) confirms this is not a "longest-run-too-low, cluster-mass-saved-it" boundary case — both statistics fail. C5 epi-X does not rescue the signal.

#### Cache + script provenance
- `audit_66_grassmann_matched_strength_surrogate.py` → `data/audit/grassmann_matched_strength_surrogate/`
- `audit_67_grassmann_epi_exclusion.py` → `data/audit/grassmann_epi_exclusion/`
- `audit_70_grassmann_cluster_extent.py` → `data/audit/grassmann_cluster_extent/`
- Object built from `eigvecs` field in `data/cache/imcoh_lrg/Pat_NN/theta_{phase}_lrg_imcoh-abs.npz` (obs) and `data/cache/matched_strength_surrogate_lrg/Pat_NN/theta_{phase}_R200_swap20_seed20260511_imcoh_abs.npz` (surrogate)

## 4. Multiscale structure of θ (null)

### Per-pair multiscale via `D_coph` (sign-inverted)

The cophenet step at θ produces a **cohort-median sign inversion**. Raw FC θ cohort median = +0.117 (7/10 trace direction by per-patient sign), raw D(τ_max) cohort median = +0.044 (6/10), but cophenet D_coph cohort median = −0.049 (3/10). The dendrogram aggregation is not preserving the substrate-level positive cohort tendency; it is actively reorganizing it into a slightly negative one. The 2/10 patients above their own matched-strength surrogate is the limiting signature: at θ, the per-pair multiscale trace is *anti-correlated* with what an ambient-FC-noise control would produce.

### Subspace multiscale via Grassmann `d_G(k)` (null)

The 5-cell longest contiguous-significant run is consistent with chance under the cluster-extent permutation. No contiguous window above the empirical 95th percentile. The θ Grassmann signature is structurally indistinguishable from a matched-strength surrogate scan.

## 5. Anatomical distribution — not relevant under no-trace verdict

No θ-specific Desikan-Killiany hypergeometric or implant-geometry analysis exists, and none is required under the `no trace` verdict. If a writing-agent needs an anatomical control for θ, the appropriate framing is **"no localizable trace because no global trace"** — the per-pair and subspace probes both reject θ at the cohort level, so anatomical localization is moot.

## 6. Patient-by-patient reading

### θ `ρ_split^coph` per-patient signs (descriptive)
From `ctm_triangle/cohort_summary.csv`: cohort median `ρ_split^coph` = −0.049 with 3/10 in + direction. The 7 patients in the − direction form a clear cohort majority against the trace direction. No single-patient analysis is decisive here because the cohort signal is itself anti-trace; per-patient inspection would only confirm what the cohort statistic shows.

### θ Grassmann subspace trace per-patient (descriptive)
Per-patient `T_G(k)` at θ shows isolated cells of cohort agreement (e.g., k=7: 3/10 below; k=11..14: 4/10 then 3/10 at the kτ window where it might have been highest), none reaching the contiguous-significant criterion. No mode-by-mode pattern emerges.

## 7. Robustness panel (under the locked battery)

| Check | θ verdict | Source |
|---|---|---|
| C1 within-rsPre split-half null | ✗ 3/10 cohort, p = 0.688 | `ctm_triangle/cohort_summary.csv` theta |
| C2 drift-floor `ρ_split > ρ_drift` | ✗ 6/10 above drift, p = 0.278 | same |
| C3 matched-strength surrogate `ρ_split^coph` | ✗✗ ratio −10.2× (anti-direction), 2/10, p = 0.722 | `matched_strength_surrogate_split_baseline` theta |
| C3 matched-strength surrogate Grassmann (cluster-extent) | ✗ cluster_p = 0.0995, 5-cell run | `grassmann_cluster_extent` theta |
| C4 cross-probe restriction | ✗ −0.048 sign-matches rho_split but only 3/10 cohort | `ctm_triangle` theta (rho_xprobe column) |
| C5 epi-X `ρ_split^coph` | **not run** | (audit_68 was α-only) |
| C5 epi-X Grassmann | no improvement; θ remains null under epi-X | `grassmann_epi_exclusion/sensitivity.csv` theta |
| Three-layer cohort table | substrate borderline + → raw D borderline + → cophenet sign-inverts to − — **band-selective sign reversal** | headline section above |

## 8. Sensitivity panel under the locked 5-control battery

| Control | θ status (D_coph) | θ status (Grassmann) | Notes |
|---|---|---|---|
| C1 split | ✗ p = 0.688 | n/a | C1 specific to cophenet |
| C2 drift | ✗ p = 0.278 | n/a | C2 specific to cophenet |
| C3 matched-strength | ✗✗ p = 0.722 (anti-direction, ratio −10.2×) | ✗ cluster_p = 0.0995 | mandatory; both probes fail |
| C4 cross-probe | ✗ −sign at 3/10 (< 6/10) | n/a | C4 specific to cophenet |
| C5 epi-X | **not run** (no audit_68 for θ cophenet) | no improvement under epi-X | sensitivity layer |

**Verdict (locked)**: `no trace` per `locked/VERDICT_LEDGER.md`. No further sensitivity tests owed; both primary probes reject θ at the 0.05 gate.

## 9. Interpretation

**Neuroscientific framing.** θ (4–8 Hz) — the canonical hippocampal-prefrontal memory-encoding rhythm in rodents and the frontal-midline working-memory rhythm in humans (Buzsáki 2002; Cavanagh & Frank 2014; Hsieh & Ranganath 2014) — does **not** carry a detectable post-task reorganization at the LRG layer in this cohort and paradigm. The cohort substrate signal is in the + direction at raw FC (cohort-median +0.117, 7/10 patients in trace direction) but neither raw FC nor raw D(τ_max) clears the matched-strength gate (p = 0.138 and 0.097 respectively, both > 0.05), and the cophenet step inverts the cohort-median sign to negative.

**Why θ is the cleanest negative.** The substrate-to-cophenet sign inversion makes θ qualitatively distinct from γ_l (cophenet step demotes but preserves +sign; raw D passes), and from β (cophenet amplifies +sign and passes). θ is the *only* band where the multiscale aggregation actively flips the cohort-median direction. This is the signature of a band that lacks structured per-pair task-related reorganization at the LRG scale being probed (τ_max).

**Mechanistically.** Two non-exclusive interpretations:
1. **The task does not produce sustained θ-band reorganization in iEEG.** The paradigm may not engage hippocampal-prefrontal θ in a way that leaves a retained post-task structural signature visible at intracranial electrodes — particularly because most patient implants do not densely sample hippocampus and medial prefrontal cortex (the canonical θ generators).
2. **The implant locations bias against θ structural retention.** sEEG implants are clinically motivated (seizure foci) and the resulting electrode coverage is dominated by lateral temporal cortex and lateral frontal cortex — not the deep midline structures where θ generators predominate. Even if θ reorganization existed, the geometry of the recording would suppress it.

**Why θ is the band-specificity benchmark.** Every "trace" band (β, α, γ_l, δ) earns its verdict against the alternative explanation that the LRG pipeline produces spurious cohort-median directions. θ's clean negative — every primary control fails, both probes reject, cohort-median is *anti*-direction at cophenet — establishes that the pipeline is not producing false positives. The β/α/γ_l/δ verdicts are therefore band-selective findings, not pipeline artefacts.

## 10. Source-of-truth references

### Lockdown documents
- `.agents/preprint/locked/CONTROLS.md` — 5-control battery
- `.agents/preprint/locked/VERDICT_LEDGER.md` — locked verdict tag `no trace` for θ
- `.agents/preprint/directives/archive/2026-05/methods_revision_2026-05-18_cophenet.md` — binding methods directive
- `.agents/preprint/bands/01_beta.md`, `.agents/preprint/bands/02_alpha.md`, `.agents/preprint/bands/03_gammalow.md` — companion band briefs

### Audit data (rows the numbers come from)
- `data/audit/ctm_triangle/cohort_summary.csv` — θ row: C1, C2, C4 statistics
- `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` — θ row: C3 cophenet (fails decisively)
- `data/audit/grassmann_cluster_extent/cohort_summary.csv` — θ cluster-extent permutation (audit_70)
- `data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv` — θ per-k T_G + surrogate stats (audit_66)
- `data/audit/grassmann_epi_exclusion/sensitivity.csv` — θ per-k under C5 epi-X (audit_67)
- `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv` — θ substrate row
- `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv` — θ raw D sensitivity

### Audit scripts (provenance)
- `audit_63_split_baseline_surrogate.py` — C3 ρ_split^coph matched-strength (full cohort)
- `audit_66_grassmann_matched_strength_surrogate.py` — Grassmann C3 per-k
- `audit_67_grassmann_epi_exclusion.py` — Grassmann C5 per-k under epi-X
- `audit_70_grassmann_cluster_extent.py` — Grassmann cluster-extent permutation (2026-05-19)

### Companion guides
- `.agents/guides/02_methods/lrg-framework-guide.md` — LRG primitive (formulas + codebase mapping)
- `.agents/guides/02_methods/imcoh-guide.md` — imcoh_abs convention + Nolte 2004 immunity
- `.agents/guides/04_rules/never-always-list.md` — coding rules + no-hardcoded-thresholds feedback
- `.agents/guides/05_plotting/README.md` — figure conventions

## 11. What follows

This brief freezes the θ `no trace` verdict for the preprint manuscript under the 2026-05-18 / 2026-05-19 lockdown. Next steps:

1. **θ figures (optional, deferred)** — θ does not need a dedicated F1–F4 figure suite because its null status is a single-panel statement. If a θ panel is wanted as a companion to the β/α/γ_l positive panels (e.g., in a cross-band verdict matrix), it can show the four-control failure grid and the Grassmann 5-cell run within the empirical null distribution. No anatomy figure is required. **Note**: the previously cited `F_cohort_2_verdict_matrix.pdf` cross-band figure is not yet produced — figure scripts deferred per `bands/00_cohort.md:157`. If a θ panel is wanted, scope it as a standalone descriptive supplement, not as part of an unrealised cross-band suite.
2. **No anatomy audit owed** at θ.
3. **No further sensitivity tests owed.** The C5 epi-X Grassmann audit (audit_67) confirms no improvement. The cophenet C5 epi-X was not run for θ (audit_68 was α-only) and is not required — the cophenet probe already fails decisively under C1/C2/C3/C4.

Verdict ready for writing-agent handoff: **`no trace`**. The headline story is the **cophenet sign inversion at θ** — the cleanest negative in the panel, useful as the band-specificity benchmark establishing that the LRG pipeline does not produce false-positive trace verdicts when no task-related reorganization is present.
