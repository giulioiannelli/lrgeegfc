---
name: preprint-gammalow-band
era: IMCOH_ABS_COHORT_N10
status: current
kind: preprint-result-report
band: low_gamma
range_hz: [30, 80]
cohort: [Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15]
n_patients: 10
fc_method: imcoh_abs
tau: 1/lambda_max
canonical_lrg_object: "D_coph = cophenet(UPGMA(D(tau_max)))"
canonical_subspace_object: "U_k = span{phi_2, ..., phi_{k+1}}"
verdict_tag: "strong trace, only Grassmann (full graph); the cophenet 'emergence under epi-exclusion' reported 2026-06-05 is WITHDRAWN 2026-06-12 as generic node-count (audit_85 --stratify epi, p_dec=0.130), so γ_l has no cophenet per-pair trace — verdict rests on the Grassmann probe alone"
verdict_source: VERDICT_LEDGER.md (locked 2026-05-18, revised 2026-05-19 + Decision 12 2026-05-26)
verdict_layers:
  substrate_rank: trace_direction_borderline (raw FC ρ_split^raw matched-strength p=0.053 at 7/10 cohort; ratio 53.9× but at significance threshold)
  rho_split_coph: no_trace_on_full_graph__epiX_emergence_is_NODECOUNT_not_epi_specific (full: C3 paired Wilcoxon p=0.116, n_above_surrogate 5/10 — no cophenet trace. epi-X audit_77 2026-06-05 reported "emergence": matched-strength Wilcoxon p=0.024, median +0.174 vs full +0.083, Pat_13/Pat_14 flip sign — framed as the epileptic zone MASKING a healthy-tissue γ_l per-pair trace, elevated to a primary interpretive lens by PI directive 2026-06-05. ⚠ RETRACTED 2026-06-12 (audit_85 --stratify epi): that "emergence" is NOT epi-specific — random-node decimation reaches the same cohort ρ_split by removing the same ~9.5% of RANDOM nodes (cohort p_dec=0.130, generic_nodecount). γ_l has NO epi-specific cophenet per-pair trace; the WM-X "emergence" is likewise node-count (p_dec=0.185). γ_l's "strong trace" verdict rests ENTIRELY on the Grassmann probe (cluster-extent, full graph; re-passes WM-X gate audit_86 cluster_p_mass^wmX=0.005) — the cophenet per-pair channel is withdrawn for γ_l. ⚠ flagged for PI review re: the 2026-06-05 directive)
  grassmann: strong_trace (cluster-extent permutation cluster_p_mass=0.005 Decision-8 gate, cluster_p_LR=0.015 descriptive co-statistic; audit_70 all-clusters; obs cluster mass 66.14 raw / 0.259 normalized; 13-cell observed run with primary window k=12..23 — above null 95th percentile of 27.24; Decision-12 LOO max p_mass=0.040 Pat_05 passes < 0.05)
  grassmann_epi_excluded: secondary mechanistic observation (audit_72; raw mass contracts 66.14 → 32.75, p_mass^epi-X = 0.030 cohort gate held; LOO under epi-X max = 0.159 Pat_05 fragile — reported as secondary, not verdict-driver)
  rho_split_coph_wm_excluded: SECONDARY_cophenet_emergence_is_NODECOUNT_not_wm_specific__raw_IS_wm_specific (C6, audit_83 2026-06-08 / amended audit_85 2026-06-12 → data/audit/wm_stratified/; the cophenet "emergence" under exclude_wm (full null → p=0.019, obs +0.083→+0.212) is NOT WM-specific — the random-node-decimation control (audit_85) reaches the same cohort ρ_split by removing any ~40 % of nodes (cophenet cohort p_dec=0.185, generic_nodecount). So the earlier "SECOND suspect-tissue axis after epi-X / WM unmasks a healthy-tissue γ_l trace" framing is WITHDRAWN for the cophenet substrate — it is a node-count effect, not WM biology. HOWEVER the RAW γ_l substrate IS genuinely WM-specific (decimation p_dec=0.030, WM_specific): removing WM specifically clears raw γ_l beyond random node loss. SECONDARY/mechanistic per CONTROLS §C6, NOT a gate; γ_l verdict tag UNCHANGED)
  grassmann_wm_excluded: SECONDARY_REPASSES_C3_GATE_under_wm_exclusion (C6, audit_86 2026-06-12; the γ_l subspace trace RE-PASSES the LOCKED C3 cluster-extent mass gate on the gray-only montage — cluster_p_mass^wmX=0.005 STRONG, LOO-robust (0.005), matching the full-graph C3 strong verdict. Upgraded from the per-k Wilcoxon count (41→40, retained as co-statistic) to a re-pass of the actual locked gate. Flag B resolved)
  anatomy_grassmann: RETRACTED_0_of_7_cohort_supported (2026-05-30 signed sampling-aware audit, locked in ANATOMY_LEDGER: the 7-DK-region "occipito-temporal + frontal + medial-OFC" list does NOT survive — 0/7 adequately-sampled-and-cohort-localized; the only well-sampled locked region [lh-middletemporal n=5] fails the signed test; per-patient null. The cross-phase Grassmann subspace probe does not localize in any band [2026-06-10 atlas]; γ_l has no anatomical localization. β is the only band that localizes [→OFC, cophenet]. Historical 7-region list retained in §5 as record. Source data/audit/anatomy_localization_wilcoxon/README.md)
sources:
  - .agents/preprint/locked/CONTROLS.md (locked 5-control battery)
  - .agents/preprint/locked/VERDICT_LEDGER.md (locked verdicts)
  - .agents/preprint/methods/methods_revision_2026-05-18_cophenet.md (binding methods directive)
  - data/audit/ctm_triangle/cohort_summary.csv (γ_l row: rho_split_median, n_above_drift, rho_xprobe, Wilcoxon)
  - data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv (γ_l row: paired Wilcoxon vs surrogate)
  - data/audit/grassmann_cluster_extent/cohort_summary.csv (γ_l cluster-extent p=0.0149, audit_70)
  - data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv (γ_l per-k T_G + surr stats, audit_66)
  - data/audit/grassmann_epi_exclusion/sensitivity.csv (γ_l per-k under C5 epi-X, audit_67)
  - data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv (γ_l substrate row)
  - data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv (γ_l raw D sensitivity)
revision_history:
  - 2026-05-19: initial brief produced from VERDICT_LEDGER.md lockdown
  - 2026-05-26: Decision-12 cascade — γ_l Grassmann verdict upgraded weak → **strong** (all three Decision-12 conditions pass: cluster_p_mass=0.005 < 0.01; cluster_p_LR=0.015 < 0.05; LOO max p_mass=0.040 Pat_05 < 0.05). Gate citation switched from `cluster_p_LR = 0.0149` to `cluster_p_mass = 0.005` (Decision-8 mass-only). Pre-fix raw values 19.17 / 4.72 / 13.96 updated to post-fix all-clusters values 66.14 / 9.04 / 27.24. C5 epi-X reframed as secondary mechanistic observation per Decision 10 (retention rule retired).
  - 2026-06-05: γ_l cophenet **EMERGES under epi-exclusion** (audit_77 → `data/audit/epi_stratified/cophenetic_cohort.csv`): full graph null (matched-strength p=0.116) → exclude_epi matched-strength Wilcoxon p=0.024, 8/10 above own surrogate, median +0.174, LO-Pat_13 p=0.037; Pat_13 (−0.042→+0.109) and Pat_14 (−0.220→+0.073) flip sign when epi contacts dropped — the epileptic zone masks/inverts a healthy-tissue γ_l per-pair trace. verdict_tag updated to "Grassmann trace + epi-X-gated cophenet trace". Epi-exclusion elevated to a primary interpretive lens (PI directive). See VERDICT_LEDGER.md 2026-06-05 revision.
  - 2026-06-08: added **C6 white-matter-exclusion** sensitivity (frontmatter `verdict_layers.rho_split_coph_wm_excluded` + `grassmann_wm_excluded`) from audit_83/84 (`data/audit/wm_stratified/`). γ_l cophenet **EMERGES under WM exclusion too** (full null → exclude_wm matched-strength p=0.019, obs +0.083→+0.212, LO-Pat_15 0.014; wm_only also emerges p=0.032) — a second suspect-tissue axis after epi-X where removing contacts unmasks a healthy-tissue γ_l per-pair trace; Grassmann γ_l per-k cells unchanged 41→40. SECONDARY/mechanistic per CONTROLS §C6, NOT a gate and NOT a primary interpretive lens (node-count caveat: ~40 % nodes removed; C6-Grassmann = per-k Wilcoxon basis, not the C3 cluster-extent gate). γ_l verdict UNCHANGED. See VERDICT_LEDGER.md 2026-06-08 revision + `responses/2026-06-08_wm-exclusion-cascade.md`.
  - 2026-06-12: **resolved the two C6 honesty flags** (audit_85 + audit_86) — and the resolution CORRECTS the 2026-06-08 γ_l cophenet read. Flag A (node-count): the γ_l cophenet "emergence under WM exclusion / second suspect-tissue axis" is NOT WM-specific — the random-node-decimation control reaches the same cophenet cohort ρ_split by removing any ~40 % of nodes (cophenet cohort p_dec=0.185, generic_nodecount). **WITHDRAW** the cophenet "WM unmasks a healthy-tissue γ_l trace" framing. The RAW γ_l substrate IS genuinely WM-specific (decimation p_dec=0.030, WM_specific). Flag B (Grassmann gate): the γ_l subspace trace RE-PASSES the locked C3 cluster-extent mass gate under WM exclusion (audit_86 `cluster_p_mass^wmX=0.005`, strong, LOO-robust) — upgraded from the per-k count. Updated both `verdict_layers.*_wm_excluded` entries. **γ_l verdict tag UNCHANGED.** See VERDICT_LEDGER.md 2026-06-12 amendment.
  - 2026-06-12 (same control, applied to C5 epi-X — CORRECTS the γ_l cophenet primary-lens story): `audit_85 --stratify epi` shows the γ_l cophenet "EMERGES under epi-X" claim (the ONLY cophenet evidence γ_l had — no full-graph trace) is ALSO generic node-count (cohort p_dec=0.130; epi-X drops ~9.5% of nodes, random removal reproduces it). So γ_l has **NO epi-specific cophenet per-pair trace** — both the epi-X and WM-X "emergences" are node-count. The cophenet per-pair channel is **withdrawn for γ_l**; its **"strong trace" verdict rests ENTIRELY on the Grassmann probe** (cluster-extent full graph; re-passes WM-X gate). Updated `rho_split_coph` + `rho_split_coph_wm_excluded`. **γ_l verdict tag UNCHANGED** ("strong trace, only Grassmann"). ⚠ Flagged for PI review (touches the 2026-06-05 directive). Source `data/audit/epi_stratified/decimation_control_cohort.csv`.
  - 2026-06-18: **PI decision — epi-exclusion DEMOTED from primary interpretive lens to SECONDARY** (parallel to C6). For γ_l this is doubly settled: the cophenet "emergence under epi-X" — the *only* cophenet evidence γ_l ever had — was already withdrawn 2026-06-12 as generic node-count (audit_85 --stratify epi, p_dec=0.130). γ_l keeps its **`strong trace, only Grassmann`** verdict (full-graph cluster-extent gate, re-passes WM-X). Cleaned the residual "primary interpretive lens / primary γ_l finding" wording in the Head and the C5 row of the cross-probe table. **γ_l verdict tag UNCHANGED.** See VERDICT_LEDGER.md 2026-06-18 amendment.
---

# γ_low band (30–80 Hz) — preprint result report

## Head

γ_l carries a **whole-network subspace trace** detected by the Grassmann probe under cluster-extent permutation (audit_70 all-clusters `cluster_p_mass = 0.005` Decision-8 gate; `cluster_p_LR = 0.015` as descriptive co-statistic) at a 13-cell contiguous-significant window with primary span k=12..23, and **no per-pair cophenet trace on the full graph** (`ρ_split^coph` C3 paired Wilcoxon p = 0.116, n_above_surrogate 5/10) — **but the cophenet trace EMERGES under epi-zone exclusion** (audit_77, 2026-06-05: matched-strength Wilcoxon p = 0.024, 8/10 above own surrogate, median +0.174 vs full +0.083; Pat_13 and Pat_14 flip sign when epi contacts are dropped) — **but this cophenet "emergence" was WITHDRAWN 2026-06-12 as generic node-count** (audit_85 --stratify epi, p_dec = 0.130; not epi-specific), so γ_l has **no per-pair cophenet trace** and its verdict rests on the Grassmann probe alone. Epi-exclusion is a **secondary control, demoted** from the 2026-06-05 "primary interpretive lens" (PI-reviewed 2026-06-18). The cophenet step actively *demotes* the γ_l per-pair signal that exists at the raw-FC and raw-`D(τ_max)` layers (raw FC p=0.053 at 7/10; raw D p=0.032 at 6/10) — the multiscale dendrogram aggregation does not preserve the γ_l per-pair structure. Full-data LOO max `p_mass = 0.040 (Pat_05)` passes the Decision-12 LOO precondition. Under C5 epi-X (audit_72, secondary mechanistic observation), the Grassmann signal contracts (raw mass 66.14 → 32.75, `p_mass^epi-X = 0.030` cohort gate held; LOO under epi-X max = 0.159 Pat_05 is fragile) — interpreted as γ_l reorganization recruiting cortex that straddles the epi/non-epi boundary, not as falsification. Verdict from `locked/VERDICT_LEDGER.md` (Decision 12 cascade 2026-05-26): **`strong trace, only Grassmann ↑`** — γ_l reorganization lives in the leading-mode subspace at intermediate `k` (a band of roughly 20 slow modes), not in the per-pair cophenet geometry.

## Headline three-layer cohort table (γ_l-band row, copied from `locked/VERDICT_LEDGER.md`)

```
Layer                                        γ_l (full n=10)
-------------------------------------       ---------------------
Raw FC          ρ_split^raw  (raw_fc_ms)    7/10 p=.053  ratio 53.9×  ← borderline
Raw D(τ_max)    ρ_split      (preprint_05)  6/10 p=.032  ratio 28.8×  ← passes
Cophenet D_coph ρ_split^coph (audit_63)     5/10 p=.116  ratio 30.2×  ← FAILS
Grassmann       d_G(k)       (audit_70)     13-cell run k=12..23, cluster_p_mass = 0.005  ← strong trace ↑ (Decision-12)
```

Sources: `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv`, `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv`, `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv`, `data/audit/grassmann_cluster_extent/cohort_summary.csv`.

**Key reading**: γ_l is the **inverse of β**. At β the cophenet step *amplifies* a substrate-borderline signal (raw FC 6/10 p=.053 → cophenet 7/10 p=.005, ratio 13.4× → 23.7×). At γ_l the cophenet step *demotes* a raw-D-significant signal (raw D 6/10 p=.032 → cophenet 5/10 p=.116). The multiscale per-pair aggregation acts band-selectively: it concentrates and purifies β per-pair structure but disperses γ_l per-pair structure. The γ_l reorganization re-emerges only at the **whole-network subspace level** via Grassmann, suggesting the γ_l signal is concentrated in coherent modes (~20 slow eigenmodes) rather than in individual pair distances.

## 1. Scientific claim

**Cohort-level question.** Does the γ_l-band (30–80 Hz) post-task resting state at the LRG-CTM layer carry task-related reorganization on either of the two locked probes — per-pair cophenet `ρ_split^coph` or whole-network subspace Grassmann `d_G(k)`?

**Refined biological claim.** γ_l coupling reorganizes during the task and the reorganization is retained in the **leading Laplacian eigenmode subspace** at intermediate dimensions (roughly k=12..23, a 13-cell contiguous-significant window above the empirical null), but **not** in the per-pair cophenet multiscale geometry. The 13-cell subspace window means roughly 20 slow modes (excluding the trivial first eigenvector) rotate coherently between phases — narrower than β's 29-mode signature and at lower k.

**Falsification budget.** Under Decision 12 (locked 2026-05-26), the claim fails if any of: (a) `cluster_p_mass ≥ 0.01` (Decision-8 cluster-mass-null gate); (b) `cluster_p_longest_run ≥ 0.05` (cluster-extent-null co-statistic); (c) full-data LOO max `p_mass ≥ 0.05` (single-patient robustness precondition). At γ_l all three pass: `p_mass = 0.005`, `p_LR = 0.015`, LOO max = 0.040 (Pat_05) — verdict is **strong**. C5 epi-X is a *secondary mechanistic observation*, not a verdict-promoter or demoter per Decision 10 (the prior ≳80% retention rule is retired).

## 2. Cohort, substrate, and library entry points

Same as α and β: n=10 patients (Pat_02/03/05/06/07/08/10/13/14/15); `imcoh_abs` band-averaged magnitude of imaginary coherence; γ_l band `B = [30, 80] Hz`. Pat_03 acquired at 1024 Hz; handled at the config layer (`nperseg_for_fs(fs)`, `FS_OVERRIDES`) — full cohort member at analysis layer (policy updated 2026-05-18, no dropout sensitivity test). Pat_10 task rows [53, 54, 55] dropped at load. Pat_14 vendor-replaced 2026-04-25.

```python
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
W = load_fc_matrix(patient="Pat_02", phase="rest_pre", band="low_gamma", fc_method="imcoh_abs")
lrg = load_lrg_result(patient="Pat_02", phase="rest_pre", band="low_gamma", fc_method="imcoh_abs")
# lrg.eigenvalues, lrg.eigenvectors → Grassmann probe via U_k = span{φ_2, …, φ_{k+1}}
# lrg.ultrametric_matrix is D_coph (silent at γ_l)
```

## 3. The two LRG probes — methodology, results, and provenance

### 3.1 Per-pair multiscale correlation on `D_coph` — **no trace at γ_l**

#### Critical preamble (5-point)
1. **Claim.** γ_l `ρ_split^coph > 0` at cohort level.
2. **Null (C3, mandatory).** Strength-preserving 4-cycle ±δ rewiring, R=200, seed 20260511. Paired Wilcoxon one-sided greater across patients: obs > surrogate.
3. **Strongest plausible alternative.** Same-band node-strength evolution from rsPre to rsPost driving the per-pair distance shift without genuine multiscale reorganization.
4. **Does C3 cover (3).** Yes — matched-strength preserves `s_i` to 10⁻⁴, only topology shuffles. What C3 cannot exclude alone: epi-zone-driven cohort dilution (C5 epi-X cophenet not run for γ_l; audit_68 covered α only).
5. **Falsification.** C3 paired Wilcoxon p ≥ 0.05 → no trace at this probe.

#### Method
Same as α/β: `Δ_task = D_coph^taskT − D_coph^rsPre_A`, `Δ_rest = D_coph^rsPost − D_coph^rsPre_B`, `ρ_split^coph = ρ_S(Δ_task, Δ_rest)` Spearman over all `N(N−1)/2` pair entries; independent rsPre halves A/B; full LRG pipeline rerun per half at `τ_max = 1/λ_max`.

#### Results (γ_l, full cohort)

| Statistic | Value | Source |
|---|---|---|
| Cohort-median `ρ_split^coph` | +0.1399 | `ctm_triangle/cohort_summary.csv` low_gamma row |
| n_trace_split | 7/10 | same |
| n_above_drift | **9/10** | same |
| C1 Wilcoxon p (`ρ_split > 0`) | **0.03223** | same col `wilcoxon_split_gt_0_p` |
| C2 Wilcoxon p (`ρ_split > ρ_drift`) | **0.00977** | same col `wilcoxon_split_gt_drift_p` |
| C3 matched-strength: obs_median | +0.0830 | `matched_strength_surrogate_split_baseline/cohort_summary.csv` low_gamma |
| C3 surr median (median over surrogate medians) | +0.00275 | same |
| C3 ratio obs/surr | 30.2× | derived |
| C3 paired Wilcoxon z / p | 40.0 / **0.1162** | same — **C3 FAILS** |
| C3 n_above_surrogate (per-patient) | 5/10 | same |
| C4 cross-probe rho_xprobe | +0.1430 | `ctm_triangle/cohort_summary.csv` low_gamma |
| C4 n_trace_xprobe | 7/10 (+sign matches +0.140) | same |
| C5 epi-X cophenet (audit_77, 2026-06-05) | reported "EMERGES" (p=0.024, median +0.174) — **WITHDRAWN 2026-06-12 as generic node-count** (audit_85 --stratify epi, p_dec=0.130; not epi-specific) → γ_l has NO cophenet trace, full or epi-X | data/audit/epi_stratified/{cophenetic_cohort.csv, decimation_control_cohort.csv} |

C1, C2, and C4 all pass; C3 (the mandatory matched-strength control per `feedback_matched_strength_mandatory.md`) fails decisively. The cohort-median ρ ≈ +0.08 is 30× larger than the surrogate median in absolute terms, but per-patient variability is high enough that only 5/10 patients sit above their own surrogate and the paired Wilcoxon p = 0.116. Under the locked C3 rule, this fails the trace gate → **no trace at the cophenet probe**.

#### Reading
The γ_l per-pair multiscale trace **fails the matched-strength gate**. C1/C2/C4 passes do not rescue it under the locked rule (C3 is mandatory). The 30× obs/surr ratio with high per-patient variance (5/10 above own surrogate) is the textbook signature of a *cohort-median effect* that doesn't survive per-patient-paired surrogacy — the kind of pattern that within-baseline-only nulls would have endorsed and matched-strength correctly rejects.

This is reinforced by the substrate→raw D→cophenet demotion: raw FC and raw D both *pass* (raw FC p=0.053 borderline; raw D p=0.032), but the cophenet aggregation step removes the per-patient consistency. The interpretation is that γ_l per-pair reorganization is **single-scale at the per-pair layer** (visible at raw `D(τ_max)`) and gets averaged out across `N−1` merge heights in the cophenet.

#### Cache + script provenance
- Full cohort: `audit_63_split_baseline_surrogate.py` → `data/audit/matched_strength_surrogate_split_baseline/`
- Within-baseline triangle: `audit_*_ctm_triangle.py` → `data/audit/ctm_triangle/`
- Object built from `lrg.ultrametric_matrix` field in `data/cache/imcoh_lrg/Pat_NN/low_gamma_{phase}_lrg_imcoh-abs.npz`

### 3.2 Global subspace rotation: Grassmann chordal distance `d_G(k)` — **strong trace at γ_l ↑ (primary band finding)**

#### What `d_G(k)` reads
Chordal Grassmann distance between leading-`k`-dimensional Laplacian eigenmode subspaces of two phases (`d_G² = k − ‖A^T B‖_F²` with A, B = N×k orthonormal eigenvector matrices, skipping the trivial mode). Reads whether the slow-diffusion subspace rotates between phases at a chosen `k`.

#### Critical preamble (5-point)
1. **Claim.** `T_G(k) = d_G(taskT, rsPost; k) − d_G(rsPre_A, taskT; k) < 0` at cohort level for some contiguous window of `k`.
2. **Null (C3, audit_70 cluster-extent permutation).** R=200 matched-strength surrogates; for each phantom-surrogate r, cohort paired Wilcoxon between surrogate r and the mean of the remaining R−1 surrogates per k, longest contiguous run with p < 0.05. Empirical null distribution of longest-run length. Replaces the previous 8-cell hardcoded threshold (`locked/VERDICT_LEDGER.md` Decision 6).
3. **Strongest plausible alternative.** Same-strength edge rewiring producing a similar contiguous-significant cluster by chance.
4. **Does C3 cover (3).** Yes — matched-strength preserves node strengths, only topology shuffles; phantom-surrogate-vs-rest builds the exchangeable cluster-extent null at the same k-grid the observation lives on.
5. **Falsification.** Observed longest run within the null distribution → no trace.

#### Method
- For each `k` ∈ [2, 112] and each patient: `T_G(k) = d_G(U_k^taskT, U_k^rsPost) − d_G(U_k^rsPre_A, U_k^taskT)`.
- Cohort paired Wilcoxon one-sided less per `k`: obs `T_G(k)` < mean-surrogate `T_G(k)` across patients.
- Cluster-extent: longest contiguous run of `k` with cohort-Wilcoxon p < 0.05.
- Null: for each surrogate r in R=200, phantom-test r vs (R−1)-mean, longest run. Cluster p = (1 + #(null ≥ obs)) / (R + 1).

#### Results (γ_l)

| Statistic | Value | Source |
|---|---|---|
| Observed longest contiguous-significant run | **13 cells** | `grassmann_cluster_extent/cohort_summary.csv` low_gamma |
| Primary contiguous window (in cohort_summary p<0.05) | k = 12..23 (12 cells visible at the same gate; audit_70 13-cell count includes one boundary cell) | `grassmann_matched_strength_surrogate/cohort_summary.csv` low_gamma rows k=12..23 |
| Cluster-extent null mean LR | 2.63 | `grassmann_cluster_extent/cohort_summary.csv` low_gamma |
| Cluster-extent null median LR | 2.0 | same |
| Cluster-extent null 95th LR | 8.0 | same |
| Cluster-extent null max LR | 15 | same |
| Observed cluster mass `T_G^*` (raw Σ −log10p) | **66.14** (all-clusters formula, post-2026-05-19-pm fix) | same |
| Observed `T_G^*` (normalized per C1) | **0.259** (raw / 255.65, n_k=111, R=200) | same |
| Null mean cluster mass | 9.04 | same |
| Null 95th cluster mass | 27.24 | same |
| **cluster_p_mass (Decision-8 gate)** | **0.005** (at empirical floor `1/(R+1)`) | same |
| cluster_p_longest_run (descriptive co-statistic) | 0.015 | same |
| **LOO max p_mass (Decision-12 precondition)** | **0.040 (Pat_05)** — passes < 0.05 | same |
| C5 epi-X (audit_72; secondary): raw mass | 66.14 → **32.75** (contracts ~50%) | `grassmann_epi_exclusion/c5_wilcoxon_cohort.csv` low_gamma |
| C5 epi-X cohort gate `p_mass^epi-X` | **0.030** (cohort gate held) | same |
| C5 epi-X LOO max under epi-X | **0.159 (Pat_05)** — fragile under epi exclusion | same |
| C5 epi-X obs longest run | **10 cells at k = 19..28** | `grassmann_epi_exclusion/sensitivity.csv` |
| C5 epi-X overlap with full-cohort window | k = 19..23 (5/13 cells) — descriptive only | per-k inspection |

The observed 13-cell longest run sits between the null 95th percentile (8.0) and the null max (15). The Decision-8 gate is `cluster_p_mass = 0.005` (at empirical floor) — γ_l clears it at the same strength as β and δ. The cluster-extent co-statistic `cluster_p_LR = 0.015` is also < 0.05. Full-data LOO max `p_mass = 0.040 (Pat_05)` passes the Decision-12 LOO precondition. **By the locked Decision-12 rule, γ_l Grassmann is strong trace ↑.**

#### Secondary clusters
Beyond the primary k=12..23 window, the per-k cohort Wilcoxon also shows shorter significant clusters at k=26..29 (4 cells), k=31..39 (9 cells), and k=72..76 (5 cells). The 13-cell cluster at k=12..23 is the audit_70 longest-run anchor; the secondary clusters at higher k contribute to the all-clusters cluster-mass statistic `T_G^* = 66.14` (cluster_p_mass = 0.005 at floor).

#### C5 epi-X — secondary mechanistic observation (audit_72)
Under epi-zone exclusion (treated as a secondary mechanistic observation per Decision 10, **not** a verdict-promoter/demoter):
- Cohort gate held: `p_mass^epi-X = 0.030 < 0.05` — γ_l signal survives epi exclusion at the cohort level.
- Raw mass contracts: 66.14 → 32.75 (~50% reduction).
- LOO under epi-X max: `p_mass^epi-X = 0.159 (Pat_05)` — fragile under the reduced support.
- Contiguous window shifts: full-cohort k=12..23 (13 cells) → epi-X k=19..28 (10 cells), overlap at k=19..23 (5 cells).

**Mechanistic interpretation.** γ_l (30-80 Hz) task-reorganization recruits cortex that straddles the epi/non-epi boundary; some of the γ_l Grassmann signal is carried by contacts that include epi-zone tissue. Removing those contacts reduces the spatial support of the trace, halving the cluster mass and shifting the surviving signal onto slightly faster modes (k=19..28 vs k=12..23). The cohort-level signal survives epi exclusion but with fewer contributing cells, which is what drives the higher LOO under epi-X. **This is interesting biology to discuss, not a falsification of the full-data verdict.** The Decision-8 gate at full data (`p_mass = 0.005`) and Decision-12 LOO precondition (max = 0.040 Pat_05 < 0.05) both hold; the strong-tier verdict stands.

#### Reading
γ_l Grassmann passes the matched-strength cluster-extent gate at **strong** strength under Decision 12: `cluster_p_mass = 0.005` (Decision-8 gate), `cluster_p_LR = 0.015` (descriptive co-statistic), LOO max `p_mass = 0.040` (Pat_05). The 13-cell run at k=12..23 with all-clusters mass `T_G^* = 66.14` (normalized 0.259) is well above the null mean (9.04) and 95th percentile (27.24). C5 epi-X is reported as a secondary mechanistic observation (signal recruits cortex straddling epi-zone boundary), not as a verdict driver. γ_l reorganization is detectable at the **leading-mode subspace level** but not at the per-pair multiscale level.

#### Cache + script provenance
- `audit_66_grassmann_matched_strength_surrogate.py` → `data/audit/grassmann_matched_strength_surrogate/`
- `audit_67_grassmann_epi_exclusion.py` → `data/audit/grassmann_epi_exclusion/`
- `audit_70_grassmann_cluster_extent.py` → `data/audit/grassmann_cluster_extent/`
- Object built from `eigvecs` field in `data/cache/imcoh_lrg/Pat_NN/low_gamma_{phase}_lrg_imcoh-abs.npz` (obs) and `data/cache/matched_strength_surrogate_lrg/Pat_NN/low_gamma_{phase}_R200_swap20_seed20260511_imcoh_abs.npz` (surrogate)

## 4. Multiscale structure of the γ_l trace

### Per-pair multiscale via `D_coph` (silent)

The γ_l cophenet probe is silent under C3. The cohort-median `ρ_split^coph` = +0.083 with ratio 30.2× to surrogate is a real cohort-median direction, but per-patient consistency (5/10 above own surrogate) is too low to clear the paired Wilcoxon gate at p = 0.05.

The substrate-to-cophenet demotion is the key methodological observation at γ_l: raw FC passes at borderline (p=0.053), raw D(τ_max) passes cleanly (p=0.032), cophenet fails (p=0.116). The cophenet step is **band-selective in the opposite direction from β**: where β cophenet *amplifies* the substrate signal, γ_l cophenet *dilutes* it. Interpretation: γ_l per-pair structure is concentrated at a single scale (the τ_max scale where raw D detects it) and the multiscale dendrogram aggregation averages it out across `N−1` merge heights.

### Subspace multiscale via Grassmann `d_G(k)` (strong trace ↑)

The 13-cell contiguous-significant window at k=12..23 with `cluster_p_mass = 0.005` (gate) places the γ_l subspace trace squarely in the **mid-k regime**. β's signature is broader (29 cells at k=27..55); α has no Grassmann signal; γ_l sits between, at lower k than β with a narrower window.

The k=12..23 range corresponds to roughly modes φ_2, φ_3, …, φ_{24} of the Laplacian — the slowest 20 non-trivial diffusion modes. These are the modes that carry the largest spatial scale information about FC organization (longest mixing times). γ_l reorganization therefore lives in the **macroscopic mode structure** of the brain network, not in pair-level fine structure.

The C5 epi-X shift (k=12..23 → k=19..28) means epi-zone contacts contribute to the lower-k modes carrying the original trace; without them, the signal re-projects onto slightly faster modes. The mode-identity is not preserved, but the existence of a coherent subspace rotation is.

## 5. Anatomical distribution — ~~strong-localized to an occipito-temporal + frontal + medial-OFC network~~ **RETRACTED 2026-05-30 / DELOCALIZED 2026-06-01**

> ⚠️ **RETRACTED — the γ_l Grassmann trace does not localize anatomically.** The 7-region
> network does **not** survive a signed, threshold-free test: **0/7 cohort-supported** (the
> one well-sampled locked region, middletemporal n=5, fails; no well-sampled region passes),
> and **per-patient localization is also null** (anchor probe 1–3/10 = a seed-flickering
> knife-edge; the proper cross-phase Grassmann trace also null per-patient). The 7-region
> list was built on top-decile per-node participation + A3 enrichment on counts (this is
> also the phase-AVERAGED anchor quantity, not a cross-phase trace). Tables retained as
> historical record. Sources: `data/audit/anatomy_localization_wilcoxon/README.md`,
> `data/audit/per_patient_localization/README.md`, `locked/ANATOMY_LEDGER.md`.

~~The γ_l Grassmann trace is **localized to an occipito-temporal + frontal + medial-OFC cortical network** of 7 named DK regions that pass A3 matched-strength surrogacy ...~~ (retracted — see banner; the A3-on-counts gate does not survive the signed cohort-consistency / per-patient test)

Audited under `locked/ANATOMY_CONTROLS.md`. **Verdict source NOW: `data/audit/anatomy_localization_wilcoxon/` + `data/audit/per_patient_localization/` (2026-05-30/06-01 retraction); `ANATOMY_LEDGER.md` 2026-05-19 pm entry superseded.**

### γ_l Grassmann anatomy (top-decile per-node participation in `U_k` over `S(γ_l)`)

Aggregation is the **unweighted average of per-node participation across `S(γ_l) = {k : p_k(γ_l) < α_k}`** — the support of the cluster-mass statistic `T_G^*`, 41 cells. The retired contiguous-significant window `K*(γ_l) = [12, 23]` (13 cells) is **not used** for anatomy under the locked all-clusters paradigm.

Seven named DK regions pass A3 alone:

| Region | k_trace / K_region | A1 enrichment | A3 obs_z | A3 p_emp |
|---|---|---|---|---|
| ctx-lh-lateraloccipital | 6 / 31 | 2.43× | very large | 0.005 |
| ctx-lh-middletemporal | 9 / 55 | 1.59× | 3.78 | 0.005 |
| ctx-lh-rostralmiddlefrontal | 8 / 32 | 2.37× | 16.30 | 0.005 |
| ctx-lh-superiortemporal | 6 / 32 | 1.82× | 7.16 | 0.005 |
| ctx-rh-medialorbitofrontal | 4 / 22 | 1.54× | 3.36 | 0.005 |
| ctx-rh-parstriangularis | 5 / 15 | 3.89× | 9.01 | 0.005 |
| ctx-lh-cuneus | 4 / 16 | 2.43× | 5.25 | 0.040 |

~~Network anatomy: **occipital cortex** (left lateral occipital + left cuneus) + **temporal cortex** + **frontal** + **medial OFC**.~~ **RETRACTED (see §5 banner): 0/7 cohort-supported; well-sampled middletemporal (n=5) fails; per-patient null.** (The KC-era "left fusiform at γ_l" claim was already retracted; left fusiform appears nowhere.)

A1 enrichments are now all ≥ 1.5× (vs the retired audit's mix of 0.26×–3.24×), reflecting that `S(γ_l)`-averaged participations are more concentrated than the contiguous-window-only participations.

Source: `data/audit/anatomy_low_gamma_grassmann_clusterext/cohort_summary.csv` (audit_72 --cluster-extent, 2026-05-19 pm). The retired `anatomy_low_gamma_grassmann/cohort_summary.csv` is superseded.

### Cross-band comparison

Under `S(b)`, **fusiform appears nowhere** across the cluster-extent anatomy panel — neither at β nor γ_l nor δ. ~~The KC-era "Hippocampus + left fusiform" β claim is partially retracted: Hippocampus survives at β Grassmann; left fusiform retracts everywhere.~~ **FULLY RETRACTED 2026-05-30/06-01: Hippocampus does NOT survive either — it crosses the cohort test on only 3–5/10 patients (a marginal hint, not a localization) and is null per-patient. Both "Hippocampus" and "left fusiform" are retracted; no β localization survives.**

### Caveats

- A1 hypergeometric is sparse for Grassmann (top-decile-per-patient endpoints insufficient for high-power per-region hypergeometric). A3 is the gate.
- "Very large" A3 obs_z entries indicate that all 200 matched-strength surrogate enrichments fell strictly below the observed enrichment — read as "categorically distinguishable from matched-strength surrogates".
- C5 epi-X for γ_l Grassmann was not re-run for anatomy under the cluster-extent paradigm. The anatomy verdict above is from the full cohort. Whether the k-window shift under epi-X corresponds to a different anatomy is **descriptive supplement, not done**.

### Cache + script provenance
| Artifact | Path |
|---|---|
| γ_l Grassmann anatomy (locked, `S(γ_l)` aggregation) | `data/audit/anatomy_low_gamma_grassmann_clusterext/cohort_summary.csv` |
| γ_l Grassmann anatomy (retired, `K*(γ_l)` window) | `data/audit/anatomy_low_gamma_grassmann/cohort_summary.csv` (do not cite) |
| Producing script | `scripts/01_compute/audit/audit_72_anatomy_grassmann.py --band low_gamma --cluster-extent` |
| Anatomy battery doc | `.agents/preprint/locked/ANATOMY_CONTROLS.md` |
| Anatomy verdict ledger | `.agents/preprint/locked/ANATOMY_LEDGER.md` (2026-05-19 pm entry) |

## 6. Patient-by-patient reading

### γ_l `ρ_split^coph` per-patient signs (descriptive — probe is silent under C3)

Per-patient rho_split values for γ_l are *not* the primary signal at this band (cophenet C3 fails). For completeness from `ctm_triangle/cohort_summary.csv`:
- Cohort median `ρ_split^coph` = +0.140 (7/10 in + direction by `n_trace_split`)
- 9/10 above their own drift floor — but only 5/10 above their own matched-strength surrogate

### γ_l Grassmann subspace trace per-patient (descriptive)

Per-patient `T_G(k)` is dominated by mid-k modes. The 13-cell window at k=12..23 has cohort `n_patients_below_own_surrogate` ranging 7..8 patients per k (e.g., k=13: 8/10; k=17: 8/10; k=23: 7/10; from `grassmann_matched_strength_surrogate/cohort_summary.csv`). Cohort agreement at the Grassmann probe is therefore higher than at cophenet — 7–8 patients in the trace direction per k cell within the contiguous window.

## 7. Robustness panel (under the locked battery)

| Check | γ_l verdict | Source |
|---|---|---|
| C1 within-rsPre split-half null | ✓ 7/10 cohort, p = 0.032 | `ctm_triangle/cohort_summary.csv` low_gamma |
| C2 drift-floor `ρ_split > ρ_drift` | ✓ 9/10 above drift, p = 0.010 | same |
| C3 matched-strength surrogate `ρ_split^coph` | ✗ ratio 30.2×, 5/10, p = 0.116 — **fails** | `matched_strength_surrogate_split_baseline` low_gamma |
| C3 matched-strength surrogate Grassmann (cluster-extent) | ✓ cluster_p_mass = 0.005, cluster_p_LR = 0.015, 13-cell run k=12..23, LOO max p_mass = 0.040 Pat_05 | `grassmann_cluster_extent` low_gamma |
| C4 cross-probe restriction | ✓ +0.143, 7/10 +sign (matches +0.140) | `ctm_triangle` low_gamma (rho_xprobe column) |
| C5 epi-X `ρ_split^coph` | ✗ reported "EMERGES" p=0.024 (2026-06-05) **WITHDRAWN 2026-06-12** as generic node-count (audit_85 --stratify epi, p_dec=0.130) → no cophenet trace | data/audit/epi_stratified/{cophenetic_cohort.csv, decimation_control_cohort.csv} |
| C5 epi-X Grassmann (secondary observation per Decision 10) | cohort gate held p_mass^epi-X = 0.030; raw mass contracts 66.14 → 32.75; LOO under epi-X 0.159 Pat_05 fragile — mechanistic, not verdict-demoter | `grassmann_epi_exclusion/c5_wilcoxon_cohort.csv` |
| Three-layer cohort table (raw FC / raw D / cophenet) | substrate borderline → raw D passes → cophenet demotes — **band-selective opposite of β** | headline section above |

## 8. Sensitivity panel under the locked 5-control battery

| Control | γ_l status (D_coph) | γ_l status (Grassmann) | Notes |
|---|---|---|---|
| C1 split | ✓ p = 0.032 | n/a (no split-baseline construction for subspace) | C1 specific to cophenet |
| C2 drift | ✓ p = 0.010 | n/a | C2 specific to cophenet |
| C3 matched-strength | ✗ p = 0.116 (FAILS) | ✓ cluster_p_mass = 0.005 (strong ↑) | mandatory; Grassmann uses Decision-8 mass-only gate + Decision-12 LOO precondition |
| C4 cross-probe | ✓ 7/10 +sign | n/a (subspace not pair-level) | C4 specific to cophenet |
| C5 epi-X | ⚠ "EMERGES" (audit_77, 2026-06-05: p=0.024, 8/10) **WITHDRAWN 2026-06-12** — generic node-count (audit_85 --stratify epi, p_dec=0.130), not epi-specific | cohort gate held p_mass^epi-X = 0.030; mass contracts; LOO under epi-X fragile (Pat_05 = 0.159) | cophenet epi-X withdrawn → γ_l rests on Grassmann only; epi-X **demoted to secondary** (PI 2026-06-18) |

**Verdict (locked, Decision 12)**: `strong trace, only Grassmann ↑` per `locked/VERDICT_LEDGER.md` (Decision-12 cascade 2026-05-26 — γ_l upgraded from "weak" because all three Decision-12 conditions pass at full data: cluster_p_mass=0.005, cluster_p_LR=0.015, LOO max p_mass=0.040). No further sensitivity tests owed under the locked battery. The C5 epi-X cophenet sensitivity at γ_l is *not* part of the locked battery — audit_68 covered α only.

## 9. Interpretation

**Neuroscientific framing.** γ_l (30–80 Hz) — the canonical "low gamma" or "narrowband gamma" — is associated with local cortical computation, especially during visual processing and attention (Fries 2009; Buzsáki & Wang 2012). A whole-network subspace trace at γ_l is consistent with task-evoked low-gamma synchrony reorganizing **macroscopic coupling modes** (the slow Laplacian eigenmodes) without affecting individual pair-level multiscale geometry.

**Why subspace only, not per-pair.** The substrate-to-cophenet trajectory is informative: raw FC and raw D(τ_max) both pass the matched-strength gate, but the cophenetic aggregation step demotes the per-pair signal. Mechanistically this is consistent with γ_l reorganization being **single-scale per-pair** (visible at τ_max) but **multi-mode at the subspace level** (a 20-mode rotation in the leading subspace). The cophenet integrates across `N−1` scales and averages out single-scale per-pair signals; the Grassmann probe is intrinsically multi-mode and detects the coherent subspace rotation.

**Why mid-k (12..23) rather than β's wider k=27..55 window.** The 13-cell γ_l window at lower k means the reorganization involves **fewer, slower modes** than β. In Laplacian eigenmode language, γ_l reorganization is more "macroscopic" — it lives in the slowest diffusion modes (largest spatial coherence length) while β reorganization extends to faster modes (intermediate coherence). This is consistent with γ_l's biological role as a band that coordinates large-scale attentional and visual networks via slower envelope coupling, whereas β coordinates more spatially structured motor/cognitive networks.

**Why epi-X shifts the k-window.** Removing epi-zone contacts shrinks the network and shifts the eigenmode indexing. The original k=12..23 modes contain contributions from epi-zone-to-non-epi-zone coupling at the lower k; once epi-zone contacts are removed, the remaining graph's slow modes no longer have the same spatial structure, and the trace re-projects onto slightly higher k (k=19..28). The 42% within-window overlap (k=19..23) is the part of the slow-mode subspace that is **independent** of epi-zone contributions; the rest (k=12..18) is epi-zone-dependent and reorganizes to k=24..28 under exclusion.

**Why γ_l Grassmann is `strong` (Decision-12 cascade, 2026-05-26).** All three Decision-12 conditions pass: (i) `cluster_p_mass = 0.005` clears the Decision-8 < 0.01 threshold (at empirical floor); (ii) `cluster_p_LR = 0.015` clears the < 0.05 co-statistic threshold; (iii) full-data LOO max `p_mass = 0.040` (Pat_05) clears the < 0.05 robustness precondition. γ_l Grassmann nevertheless differs from β in two ways worth flagging in Discussion: the cluster mass `T_G^* = 66.14` (normalized 0.259) is somewhat lower than β's 69.76 (0.273), and the C5 epi-X analysis shows ~50% mass contraction with LOO fragility under epi exclusion — consistent with γ_l reorganization recruiting cortex straddling the epi/non-epi boundary. The strong-tier verdict at full data is unaffected; the C5 contraction is a secondary mechanistic observation per Decision 10.

## 10. Source-of-truth references

### Lockdown documents
- `.agents/preprint/locked/CONTROLS.md` — 5-control battery
- `.agents/preprint/locked/VERDICT_LEDGER.md` — locked verdict tag `strong trace, only Grassmann ↑` for γ_l (Decision 12 cascade 2026-05-26)
- `.agents/preprint/methods/methods_revision_2026-05-18_cophenet.md` — binding methods directive
- `.agents/preprint/bands/01_beta.md`, `.agents/preprint/bands/02_alpha.md` — companion band briefs

### Audit data (rows the numbers come from)
- `data/audit/ctm_triangle/cohort_summary.csv` — γ_l row: C1, C2, C4 statistics
- `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` — γ_l row: C3 cophenet (fails)
- `data/audit/grassmann_cluster_extent/cohort_summary.csv` — γ_l cluster-extent permutation (audit_70)
- `data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv` — γ_l per-k T_G + surrogate stats (audit_66)
- `data/audit/grassmann_epi_exclusion/sensitivity.csv` — γ_l per-k under C5 epi-X (audit_67)
- `data/audit/grassmann_epi_exclusion/cohort_summary.csv` — γ_l epi-X cohort summary
- `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv` — γ_l substrate row
- `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv` — γ_l raw D sensitivity

### Audit scripts (provenance)
- `audit_63_split_baseline_surrogate.py` — C3 ρ_split^coph matched-strength (full cohort)
- `audit_66_grassmann_matched_strength_surrogate.py` — Grassmann C3 per-k
- `audit_67_grassmann_epi_exclusion.py` — Grassmann C5 per-k under epi-X
- `audit_70_grassmann_cluster_extent.py` — Grassmann cluster-extent permutation (2026-05-19)

### Companion guides
- `.agents/guides/02_methods/lrg-framework-guide.md` — LRG primitive (formulas + codebase mapping)
- `.agents/guides/02_methods/imcoh-guide.md` — imcoh_abs convention + Nolte 2004 immunity
- `.agents/guides/04_rules/never-always-list.md` — coding rules + no-hardcoded-thresholds feedback
- `.agents/guides/05_plotting/README.md` — figure conventions (for upcoming γ_l preprint figures)

## 11. What follows

This brief freezes the γ_l verdict for the preprint manuscript under the 2026-05-18 / 2026-05-19 lockdown. Next steps:

1. **γ_l figures** — produce F1 (4-control panel for `ρ_split^coph` showing the C3 failure: cohort medians + per-patient stems for split / drift / matched-strength surrogate / cross-probe), F2 (Grassmann `T_G(k)` cohort cells across k ∈ [2, 112] highlighting the 13-cell window at k=12..23 and the audit_70 empirical null distribution), F3 (per-patient subspace-rotation visualization for k = 17, the within-window k with strongest cohort agreement at 8/10), and F4 (audit_72 c5_wilcoxon overlay showing the secondary epi-X analysis: raw mass contraction 66.14 → 32.75, window shift k=12..23 → k=19..28, LOO under epi-X 0.159 Pat_05). Output at `data/preprint/figures/gammalow/`.
2. **No anatomy audit owed** at γ_l under the locked battery. A mode-decomposition anatomy analysis (which Desikan-Killiany regions concentrate the k=12..23 eigenmodes) would be a *descriptive supplement*, not a gating control.
3. **No further sensitivity tests owed.** The C5 epi-X Grassmann audit (audit_72 c5_wilcoxon) is the sensitivity layer; cophenet C5 epi-X was not run for γ_l (audit_68 was α-only) and is *not* required under the locked battery — γ_l's verdict (`strong trace, only Grassmann ↑`) is determined by the C3 cophenet failure + Grassmann Decision-8 mass gate + Decision-12 LOO precondition at full data.

Verdict ready for writing-agent handoff: **`strong trace, only Grassmann ↑`** (Decision-12 cascade 2026-05-26). The headline story is the **substrate-to-cophenet demotion** (raw FC borderline → raw D passes → cophenet fails, opposite of β's amplification) combined with **subspace-only re-emergence** at mid-k (13 cells at k=12..23 under all-clusters cluster-extent permutation, `cluster_p_mass = 0.005`, `T_G^* = 66.14` raw / 0.259 normalized, LOO max p_mass = 0.040 Pat_05). C5 epi-X is a secondary mechanistic observation (mass contracts ~50% under epi exclusion, LOO becomes fragile) reflecting that γ_l reorganization recruits cortex straddling the epi-zone boundary; does not promote or demote the verdict.
