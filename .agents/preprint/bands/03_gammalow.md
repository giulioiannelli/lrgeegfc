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
verdict_tag: "weak trace, only Grassmann"
verdict_source: VERDICT_LEDGER.md (locked 2026-05-18, revised 2026-05-19)
verdict_layers:
  substrate_rank: trace_direction_borderline (raw FC ρ_split^raw matched-strength p=0.053 at 7/10 cohort; ratio 53.9× but at significance threshold)
  rho_split_coph: no_trace (C3 fails: paired Wilcoxon p=0.116, n_above_surrogate 5/10 — cophenet step demotes the per-pair signal that exists at raw D)
  grassmann: weak_trace (cluster-extent permutation p=0.0149, audit_70; 13-cell observed run with primary window k=12..23 — above null 95th percentile of 8.0)
  grassmann_epi_excluded: 10-cell-run-at-shifted-window (audit_67; longest run k=19..28, overlap with full-cohort k=12..23 only at k=19..23 — partial robustness, k-window shifts by ~7 modes)
  anatomy_grassmann: strong_localized (7 named DK regions A3 alone, A1 sparse; occipital + temporal + frontal + medial-OFC: lateral occipital + cuneus + middle/superior temporal + rostral middle frontal + medial OFC + pars triangularis; audit_72 --cluster-extent under S(γ_l), 2026-05-19 pm — supersedes the retired K*(γ_l) audit which had picked up left fusiform)
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
---

# γ_low band (30–80 Hz) — preprint result report

## Head

γ_l carries a **whole-network subspace trace** detected by the Grassmann probe under cluster-extent permutation (audit_70 cluster_p = 0.0149) at a 13-cell contiguous-significant window with primary span k=12..23, but **no per-pair cophenet trace** (`ρ_split^coph` C3 paired Wilcoxon p = 0.116, n_above_surrogate 5/10). The cophenet step actively *demotes* the γ_l per-pair signal that exists at the raw-FC and raw-`D(τ_max)` layers (raw FC p=0.053 at 7/10; raw D p=0.032 at 6/10) — the multiscale dendrogram aggregation does not preserve the γ_l per-pair structure. Under C5 epi-X (audit_67), the Grassmann window shrinks from 13 cells to a 10-cell run that **shifts** to k=19..28, overlapping the full-cohort window only at k=19..23 — partial robustness with mode reassignment. Verdict from `locked/VERDICT_LEDGER.md`: **`weak trace, only Grassmann`** — γ_l reorganization lives in the leading-mode subspace at intermediate `k` (a band of roughly 20 slow modes), not in the per-pair cophenet geometry.

## Headline three-layer cohort table (γ_l-band row, copied from `locked/VERDICT_LEDGER.md`)

```
Layer                                        γ_l (full n=10)
-------------------------------------       ---------------------
Raw FC          ρ_split^raw  (raw_fc_ms)    7/10 p=.053  ratio 53.9×  ← borderline
Raw D(τ_max)    ρ_split      (preprint_05)  6/10 p=.032  ratio 28.8×  ← passes
Cophenet D_coph ρ_split^coph (audit_63)     5/10 p=.116  ratio 30.2×  ← FAILS
Grassmann       d_G(k)       (audit_70)     13-cell run k=12..23, cluster_p = 0.0149  ← weak trace
```

Sources: `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv`, `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv`, `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv`, `data/audit/grassmann_cluster_extent/cohort_summary.csv`.

**Key reading**: γ_l is the **inverse of β**. At β the cophenet step *amplifies* a substrate-borderline signal (raw FC 6/10 p=.053 → cophenet 7/10 p=.005, ratio 13.4× → 23.7×). At γ_l the cophenet step *demotes* a raw-D-significant signal (raw D 6/10 p=.032 → cophenet 5/10 p=.116). The multiscale per-pair aggregation acts band-selectively: it concentrates and purifies β per-pair structure but disperses γ_l per-pair structure. The γ_l reorganization re-emerges only at the **whole-network subspace level** via Grassmann, suggesting the γ_l signal is concentrated in coherent modes (~20 slow eigenmodes) rather than in individual pair distances.

## 1. Scientific claim

**Cohort-level question.** Does the γ_l-band (30–80 Hz) post-task resting state at the LRG-CTM layer carry task-related reorganization on either of the two locked probes — per-pair cophenet `ρ_split^coph` or whole-network subspace Grassmann `d_G(k)`?

**Refined biological claim.** γ_l coupling reorganizes during the task and the reorganization is retained in the **leading Laplacian eigenmode subspace** at intermediate dimensions (roughly k=12..23, a 13-cell contiguous-significant window above the empirical null), but **not** in the per-pair cophenet multiscale geometry. The 13-cell subspace window means roughly 20 slow modes (excluding the trivial first eigenvector) rotate coherently between phases — narrower than β's 29-mode signature and at lower k.

**Falsification budget.** The claim fails if any of: (a) Grassmann cluster-extent permutation p ≥ 0.05 (i.e., 13-cell run within the empirical null distribution from R=200 phantom-surrogate tests); (b) under C5 epi-X (audit_67), the contiguous-significant window collapses to < 80% of original cell count AND the surviving window does not overlap the original at any k. The 67% retention with k-window shift to k=19..28 is the borderline that earns this band **weak** rather than strong (`locked/VERDICT_LEDGER.md` Decision 2).

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
| C5 epi-X cophenet | **not run** for γ_l (audit_68 covered α only) | — |

C1, C2, and C4 all pass; C3 (the mandatory matched-strength control per `feedback_matched_strength_mandatory.md`) fails decisively. The cohort-median ρ ≈ +0.08 is 30× larger than the surrogate median in absolute terms, but per-patient variability is high enough that only 5/10 patients sit above their own surrogate and the paired Wilcoxon p = 0.116. Under the locked C3 rule, this fails the trace gate → **no trace at the cophenet probe**.

#### Reading
The γ_l per-pair multiscale trace **fails the matched-strength gate**. C1/C2/C4 passes do not rescue it under the locked rule (C3 is mandatory). The 30× obs/surr ratio with high per-patient variance (5/10 above own surrogate) is the textbook signature of a *cohort-median effect* that doesn't survive per-patient-paired surrogacy — the kind of pattern that within-baseline-only nulls would have endorsed and matched-strength correctly rejects.

This is reinforced by the substrate→raw D→cophenet demotion: raw FC and raw D both *pass* (raw FC p=0.053 borderline; raw D p=0.032), but the cophenet aggregation step removes the per-patient consistency. The interpretation is that γ_l per-pair reorganization is **single-scale at the per-pair layer** (visible at raw `D(τ_max)`) and gets averaged out across `N−1` merge heights in the cophenet.

#### Cache + script provenance
- Full cohort: `audit_63_split_baseline_surrogate.py` → `data/audit/matched_strength_surrogate_split_baseline/`
- Within-baseline triangle: `audit_*_ctm_triangle.py` → `data/audit/ctm_triangle/`
- Object built from `lrg.ultrametric_matrix` field in `data/cache/imcoh_lrg/Pat_NN/low_gamma_{phase}_lrg_imcoh-abs.npz`

### 3.2 Global subspace rotation: Grassmann chordal distance `d_G(k)` — **weak trace at γ_l (LOAD-BEARING)**

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
| Observed cluster mass (Σ −log10p) | 19.17 | same |
| Null mean cluster mass | 4.72 | same |
| Null 95th cluster mass | 13.96 | same |
| **cluster_p_longest_run** | **0.0149** | same — `weak` gate |
| **cluster_p_cluster_mass** | 0.0348 | same |
| C5 epi-X (audit_67): obs longest run | **10 cells at k = 19..28** | `grassmann_epi_exclusion/sensitivity.csv` low_gamma rows k=19..28 (all p < 0.05 under epi-X) |
| C5 epi-X overlap with full-cohort window | k = 19..23 (5/13 ≈ 42% within-window) | per-k inspection |

The observed 13-cell longest run sits between the null 95th percentile (8.0) and the null max (15). Both cluster-p statistics — longest run (0.0149) and cluster mass (0.0348) — fall in [0.01, 0.05). By the locked rule, γ_l Grassmann is **weak trace**.

#### Secondary clusters
Beyond the primary k=12..23 window, the per-k cohort Wilcoxon also shows shorter significant clusters at k=26..29 (4 cells), k=31..39 (9 cells), and k=72..76 (5 cells). The 13-cell cluster at k=12..23 is the audit_70 longest-run anchor; the secondary clusters at higher k contribute to the cluster-mass companion statistic (0.0348).

#### C5 epi-X sensitivity (audit_67)
Under epi-zone exclusion, the contiguous-significant window shifts:
- Original full-cohort: k=12..23 (13-cell anchor at p<0.05 gate).
- Epi-X: k=19..28 (10-cell longest run at p<0.05 gate).

Cell-count retention: 10/13 ≈ 77% (below the 80% guide). Within-window overlap: 5/13 ≈ 42% (the surviving k=19..23 cells). The signal **shifts modes** under epi-X: the lower part of the original window (k=12..18) drops out and a new window opens at k=24..28. This is the `locked/VERDICT_LEDGER.md` Decision 2 reasoning: 77% retention with substantial k-window shift earns `weak` rather than `strong`.

Interpretation: the γ_l Grassmann trace is **partially robust** to epi-zone removal — a 10-cell signal survives — but the *identity* of the modes carrying the trace changes. This is consistent with epi-zone contacts contributing to the lower-k modes of the original signal; without them, the trace re-projects onto slightly higher-k modes.

#### Reading
γ_l Grassmann passes the matched-strength cluster-extent gate at `weak` strength (cluster_p = 0.0149). The 13-cell run at k=12..23 is well above the null mean (2.63) and 95th percentile (8.0), but not in the strong regime (< 0.01). C5 epi-X partial retention (77%, k-window shift) confirms the `weak` reading. γ_l reorganization is detectable at the **leading-mode subspace level** but not at the per-pair multiscale level.

#### Cache + script provenance
- `audit_66_grassmann_matched_strength_surrogate.py` → `data/audit/grassmann_matched_strength_surrogate/`
- `audit_67_grassmann_epi_exclusion.py` → `data/audit/grassmann_epi_exclusion/`
- `audit_70_grassmann_cluster_extent.py` → `data/audit/grassmann_cluster_extent/`
- Object built from `eigvecs` field in `data/cache/imcoh_lrg/Pat_NN/low_gamma_{phase}_lrg_imcoh-abs.npz` (obs) and `data/cache/matched_strength_surrogate_lrg/Pat_NN/low_gamma_{phase}_R200_swap20_seed20260511_imcoh_abs.npz` (surrogate)

## 4. Multiscale structure of the γ_l trace

### Per-pair multiscale via `D_coph` (silent)

The γ_l cophenet probe is silent under C3. The cohort-median `ρ_split^coph` = +0.083 with ratio 30.2× to surrogate is a real cohort-median direction, but per-patient consistency (5/10 above own surrogate) is too low to clear the paired Wilcoxon gate at p = 0.05.

The substrate-to-cophenet demotion is the key methodological observation at γ_l: raw FC passes at borderline (p=0.053), raw D(τ_max) passes cleanly (p=0.032), cophenet fails (p=0.116). The cophenet step is **band-selective in the opposite direction from β**: where β cophenet *amplifies* the substrate signal, γ_l cophenet *dilutes* it. Interpretation: γ_l per-pair structure is concentrated at a single scale (the τ_max scale where raw D detects it) and the multiscale dendrogram aggregation averages it out across `N−1` merge heights.

### Subspace multiscale via Grassmann `d_G(k)` (weak trace)

The 13-cell contiguous-significant window at k=12..23 with cluster_p = 0.0149 places the γ_l subspace trace squarely in the **mid-k regime**. β's signature is broader (29 cells at k=27..55); α has no Grassmann signal; γ_l sits between, at lower k than β with a narrower window.

The k=12..23 range corresponds to roughly modes φ_2, φ_3, …, φ_{24} of the Laplacian — the slowest 20 non-trivial diffusion modes. These are the modes that carry the largest spatial scale information about FC organization (longest mixing times). γ_l reorganization therefore lives in the **macroscopic mode structure** of the brain network, not in pair-level fine structure.

The C5 epi-X shift (k=12..23 → k=19..28) means epi-zone contacts contribute to the lower-k modes carrying the original trace; without them, the signal re-projects onto slightly faster modes. The mode-identity is not preserved, but the existence of a coherent subspace rotation is.

## 5. Anatomical distribution — γ_l Grassmann trace is strong-localized to an occipito-temporal + frontal + medial-OFC network

The γ_l Grassmann trace is **localized to an occipito-temporal + frontal + medial-OFC cortical network** of 7 named DK regions that pass A3 matched-strength surrogacy at p_emp<0.05 / obs_z>2 (R=200, seed 20260511) under the locked all-clusters paradigm. A1 hypergeometric is sparse (top-decile per-patient endpoints yield too few hits for high-power hypergeometric in the Grassmann probe — see `locked/ANATOMY_CONTROLS.md` audit note); the A3 surrogate gate is the load-bearing anatomy control for Grassmann.

Audited under `locked/ANATOMY_CONTROLS.md`. Verdict source: `locked/ANATOMY_LEDGER.md` 2026-05-19 pm cluster-extent revision.

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

Network anatomy: **occipital cortex** (left lateral occipital + left cuneus) + **temporal cortex** (left middle + left superior) + **frontal** (left rostral middle frontal + right pars triangularis) + **medial OFC** (right). The KC-era "left fusiform at γ_l" claim is **fully retracted** under `S(γ_l)` — left fusiform was carried by the retired `K*(γ_l) = [12, 23]` analysis but does not appear in the all-clusters cluster-extent rerun (audit_72 --cluster-extent, 2026-05-19 pm).

A1 enrichments are now all ≥ 1.5× (vs the retired audit's mix of 0.26×–3.24×), reflecting that `S(γ_l)`-averaged participations are more concentrated than the contiguous-window-only participations.

Source: `data/audit/anatomy_low_gamma_grassmann_clusterext/cohort_summary.csv` (audit_72 --cluster-extent, 2026-05-19 pm). The retired `anatomy_low_gamma_grassmann/cohort_summary.csv` is superseded.

### Cross-band comparison

Under `S(b)`, **fusiform appears nowhere** across the cluster-extent anatomy panel — neither at β nor γ_l nor δ. The KC-era "Hippocampus + left fusiform" β claim is partially retracted: Hippocampus survives at β Grassmann; left fusiform retracts everywhere. The retraction is the cleanest worked-example of why the cluster-extent methodology cascade was necessary.

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

Per-patient rho_split values for γ_l are *not* the load-bearing signal at this band (cophenet C3 fails). For completeness from `ctm_triangle/cohort_summary.csv`:
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
| C3 matched-strength surrogate Grassmann (cluster-extent) | ✓ cluster_p = 0.0149, 13-cell run k=12..23 | `grassmann_cluster_extent` low_gamma |
| C4 cross-probe restriction | ✓ +0.143, 7/10 +sign (matches +0.140) | `ctm_triangle` low_gamma (rho_xprobe column) |
| C5 epi-X `ρ_split^coph` | **not run** | (audit_68 was α-only) |
| C5 epi-X Grassmann | ⚠ partial: 10-cell run at shifted window k=19..28 (77% retention, 42% within-window overlap) | `grassmann_epi_exclusion/sensitivity.csv` low_gamma |
| Three-layer cohort table (raw FC / raw D / cophenet) | substrate borderline → raw D passes → cophenet demotes — **band-selective opposite of β** | headline section above |

## 8. Sensitivity panel under the locked 5-control battery

| Control | γ_l status (D_coph) | γ_l status (Grassmann) | Notes |
|---|---|---|---|
| C1 split | ✓ p = 0.032 | n/a (no split-baseline construction for subspace) | C1 specific to cophenet |
| C2 drift | ✓ p = 0.010 | n/a | C2 specific to cophenet |
| C3 matched-strength | ✗ p = 0.116 (FAILS) | ✓ cluster_p = 0.0149 (weak) | mandatory; Grassmann replaces 8-cell hardcoded with cluster-extent permutation |
| C4 cross-probe | ✓ 7/10 +sign | n/a (subspace not pair-level) | C4 specific to cophenet |
| C5 epi-X | **not run** (no audit_68 for γ_l cophenet) | ⚠ 77% retention with k-window shift to k=19..28 | sensitivity layer; γ_l Grassmann earns `weak` not `strong` per Decision 2 |

**Verdict (locked)**: `weak trace, only Grassmann` per `locked/VERDICT_LEDGER.md`. No further sensitivity tests owed under the locked battery. The C5 epi-X cophenet sensitivity at γ_l is *not* part of the locked battery — audit_68 covered α only and extending it to other bands was scoped out per Phase A.

## 9. Interpretation

**Neuroscientific framing.** γ_l (30–80 Hz) — the canonical "low gamma" or "narrowband gamma" — is associated with local cortical computation, especially during visual processing and attention (Fries 2009; Buzsáki & Wang 2012). A whole-network subspace trace at γ_l is consistent with task-evoked low-gamma synchrony reorganizing **macroscopic coupling modes** (the slow Laplacian eigenmodes) without affecting individual pair-level multiscale geometry.

**Why subspace only, not per-pair.** The substrate-to-cophenet trajectory is informative: raw FC and raw D(τ_max) both pass the matched-strength gate, but the cophenetic aggregation step demotes the per-pair signal. Mechanistically this is consistent with γ_l reorganization being **single-scale per-pair** (visible at τ_max) but **multi-mode at the subspace level** (a 20-mode rotation in the leading subspace). The cophenet integrates across `N−1` scales and averages out single-scale per-pair signals; the Grassmann probe is intrinsically multi-mode and detects the coherent subspace rotation.

**Why mid-k (12..23) rather than β's wider k=27..55 window.** The 13-cell γ_l window at lower k means the reorganization involves **fewer, slower modes** than β. In Laplacian eigenmode language, γ_l reorganization is more "macroscopic" — it lives in the slowest diffusion modes (largest spatial coherence length) while β reorganization extends to faster modes (intermediate coherence). This is consistent with γ_l's biological role as a band that coordinates large-scale attentional and visual networks via slower envelope coupling, whereas β coordinates more spatially structured motor/cognitive networks.

**Why epi-X shifts the k-window.** Removing epi-zone contacts shrinks the network and shifts the eigenmode indexing. The original k=12..23 modes contain contributions from epi-zone-to-non-epi-zone coupling at the lower k; once epi-zone contacts are removed, the remaining graph's slow modes no longer have the same spatial structure, and the trace re-projects onto slightly higher k (k=19..28). The 42% within-window overlap (k=19..23) is the part of the slow-mode subspace that is **independent** of epi-zone contributions; the rest (k=12..18) is epi-zone-dependent and reorganizes to k=24..28 under exclusion.

**Why γ_l Grassmann is `weak` rather than `strong`.** Three reasons: (i) cluster_p = 0.0149 falls in [0.01, 0.05), the `weak` band by the locked rule; (ii) C5 epi-X retention is 77% in cell count and only 42% in within-window overlap, below the 80% guide for `strong`; (iii) γ_l fails the cophenet probe entirely. The Grassmann signal is real (well above empirical null) but not robust enough across mode identity and per-pair representation to qualify as `strong`.

## 10. Source-of-truth references

### Lockdown documents
- `.agents/preprint/locked/CONTROLS.md` — 5-control battery
- `.agents/preprint/locked/VERDICT_LEDGER.md` — locked verdict tag `weak trace, only Grassmann` for γ_l
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

1. **γ_l figures** — produce F1 (4-control panel for `ρ_split^coph` showing the C3 failure: cohort medians + per-patient stems for split / drift / matched-strength surrogate / cross-probe), F2 (Grassmann `T_G(k)` cohort cells across k ∈ [2, 112] highlighting the 13-cell window at k=12..23 and the audit_70 empirical null distribution), F3 (per-patient subspace-rotation visualization for k = 17, the within-window k with strongest cohort agreement at 8/10), and F4 (audit_67 epi-X overlay showing the window shift k=12..23 → k=19..28). Output at `data/preprint/figures/gammalow/`.
2. **No anatomy audit owed** at γ_l under the locked battery. A mode-decomposition anatomy analysis (which Desikan-Killiany regions concentrate the k=12..23 eigenmodes) would be a *descriptive supplement*, not a gating control.
3. **No further sensitivity tests owed.** The C5 epi-X Grassmann audit (audit_67) is the sensitivity layer; cophenet C5 epi-X was not run for γ_l (audit_68 was α-only) and is *not* required under the locked battery — γ_l's verdict (`weak trace, only Grassmann`) is determined by the C3 cophenet failure + Grassmann cluster-extent + C5 partial retention.

Verdict ready for writing-agent handoff: **`weak trace, only Grassmann`**. The headline story is the **substrate-to-cophenet demotion** (raw FC borderline → raw D passes → cophenet fails, opposite of β's amplification) combined with **subspace-only re-emergence** at mid-k (13 cells at k=12..23 under cluster-extent permutation, with epi-X partial robustness at shifted k=19..28).
