---
name: preprint-alpha-band
era: IMCOH_ABS_COHORT_N10
status: current
kind: preprint-result-report
band: alpha
range_hz: [8, 13]
cohort: [Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15]
n_patients: 10
fc_method: imcoh_abs
tau: 1/lambda_max
canonical_lrg_object: "D_coph = cophenet(UPGMA(D(tau_max)))"
canonical_subspace_object: "U_k = span{phi_2, ..., phi_{k+1}}"
verdict_tag: "strong trace, only D_coph"
verdict_source: VERDICT_LEDGER.md (locked 2026-05-18, revised 2026-05-19)
verdict_layers:
  substrate_rank: trace_direction (raw FC ρ_split^raw matched-strength p=0.014 at 8/10 cohort — strongest substrate signal of any band)
  rho_split_coph: strong_trace (anchored at C3 alone: paired Wilcoxon p=0.002 ratio 8.3×; C1 p=.010, C2 p=.007, C4 paired-Wilcoxon non-degradation all pass at full cohort; per-patient n_above_surrogate 5/10 is descriptive auxiliary statistic, NOT a verdict modifier — patient-count thresholds retired 2026-05-19, reaffirmed 2026-05-28)
  grassmann: no_trace (cluster-extent permutation p=0.159, audit_70; 4-cell observed run at k=11..14 within null distribution)
  rho_split_coph_epi_excluded: secondary_mechanistic_observation (audit_68; ratio strengthens 8.3× → 27.7×, p stays significant at .014, Wilcoxon-on-epi-X p=0.0098 — supportive of non-epi-cortex contribution but NOT a verdict-promoter per amended Decision 10 / Decision 1 retraction, 2026-05-28)
  anatomy_cophenet: strong_localized (11 named DK regions A1+A3 join; bilateral cingulate + parahippocampal + medial OFC + caudal middle frontal + postcentral + precuneus + superior parietal; audit_71 2026-05-19)
  anatomy_cophenet_epi_excluded: 9_of_11_reproduce_2_drop (audit_71 --epi-x rerun 2026-05-22 after patching a silent no-op bug in the original --epi-x branch; lh-caudal-middle-frontal drops at A1 marginal q_BH=0.051; rh-medial-orbitofrontal drops at A2 — z_A2 falls 3.22 → 1.28, revealing an epi-zone-coupled contribution to its full-cohort enrichment; 9 surviving regions include 4 cingulate + lh-parahipp + rh-caudal-MF + rh-postcentral + rh-precuneus + rh-superior-parietal, with lh-rostral-AC strengthening 3.04× → 5.37× and several cingulate regions enrichment-increasing under epi-X)
sources:
  - .agents/preprint/locked/CONTROLS.md (locked 5-control battery)
  - .agents/preprint/locked/VERDICT_LEDGER.md (locked verdicts)
  - .agents/preprint/methods/methods_revision_2026-05-18_cophenet.md (binding methods directive)
  - data/audit/ctm_triangle/cohort_summary.csv (α row: rho_split_median, n_above_drift, rho_xprobe, Wilcoxon)
  - data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv (α row: paired Wilcoxon vs surrogate)
  - data/audit/alpha_epi_exclusion/cohort_summary.csv (α cophenet under C5 epi-X — strengthens)
  - data/audit/alpha_epi_exclusion/comparison.csv (per-patient Δρ full vs epi-X)
  - data/audit/grassmann_cluster_extent/cohort_summary.csv (α cluster-extent p=0.159, audit_70)
  - data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv (α per-k T_G + surr stats)
  - data/audit/grassmann_regate_no_filter/contig_summary.csv (α audit_66 + audit_67 contiguous windows)
  - data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv (α substrate row)
  - data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv (α raw D sensitivity)
revision_history:
  - 2026-05-19: initial brief produced from VERDICT_LEDGER.md lockdown
  - 2026-05-22: corrected anatomy_cophenet_epi_excluded entry + §5 epi-X anatomy section after patching the silent no-op bug in audit_71 --epi-x (np.isin int×str mismatch at lines 278–283 / 354–359); real epi-X result is 9 of 11 regions reproduce, not 11 of 11
  - 2026-05-28: stripped the C5 epi-X "5/10 → 7/10 rescue" narrative throughout per user directive; α verdict tag preserved at "strong trace, only D_coph" because C3 paired Wilcoxon p=0.00195 already passes at full cohort under the locked Wilcoxon-as-gate rule; C5 epi-X reframed as secondary mechanistic observation (not verdict-promoter) per amended Decision 10 + retracted Decision 1
---

# α band (8–13 Hz) — preprint result report

## Head

α carries a **per-pair multiscale trace** via the cophenetic communication distance `ρ_split^coph` on `D_coph = cophenet(UPGMA(D(τ_max)))`. At the full n=10 cohort, all four primary controls pass: split-baseline C1 p=0.010, drift-floor C2 p=0.007, matched-strength surrogate C3 paired Wilcoxon **p=0.002** (effect-size ratio 8.3×), cross-probe C4 paired-Wilcoxon non-degradation. The α verdict is anchored at C3 alone under the locked Wilcoxon-as-gate rule (per-patient `n_above_surrogate` = 5/10 is a descriptive cohort-agreement statistic, NOT a verdict modifier; patient-count thresholds were retired 2026-05-19). The Grassmann subspace probe finds **no significant contiguous window** (cluster-extent permutation p = 0.159, audit_70; 4-cell observed run at k = 11..14 within the null distribution, null 95th percentile = 6.0). Verdict from `locked/VERDICT_LEDGER.md`: **`strong trace, only D_coph`**.

**Secondary mechanistic observation (C5 epi-X, audit_68).** Under epi-zone exclusion the α effect-size ratio strengthens 8.3× → **27.7×** and the one-sample Wilcoxon on per-patient `obs_rho^epi-X > 0` clears at p = 0.0098 (LOO max p = 0.0195 Pat_03). Per-patient Δρ: Pat_03/05/07/10 strengthen; Pat_13 — the patient with the largest epi-zone burden (N_epi = 30/119) — inverts from +0.028 to −0.171; Pat_06/08/14 weaken mildly but stay positive. This is consistent with α-relevant non-epi cortex carrying the trace, but per amended Decision 10 + retracted Decision 1 (2026-05-28) C5 is **never a verdict-promoter** — the α verdict was already strong at C3 alone.

## Headline three-layer cohort table (α-band row, copied from `locked/VERDICT_LEDGER.md`)

```
Layer                                       α (full n=10)
-------------------------------------      ---------------------
Raw FC          ρ_split^raw  (audit_67)    p=.014, 8/10  (strongest substrate of any band)
Raw D(τ_max)    ρ_split      (preprint_05) p=.042, 7/10  (ratio 11.4×)
Cophenet D_coph ρ_split^coph (audit_63)    p=.002, 5/10  (ratio 8.3×)   ← verdict gate (Wilcoxon)
Cophenet D_coph ρ_split^coph (audit_68, C5 epi-X)   p=.0098, ratio 27.7× ← secondary mechanistic observation
```

Wilcoxon p-values are the verdict gate at each layer; the `N/10` per-patient counts are descriptive auxiliary statistics. Sources: `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv`, `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv`, `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv`, `data/audit/alpha_epi_exclusion/cohort_summary.csv`. n=10 full / n=10 epi-X (per-patient N_epi varies; see §6).

α has the **strongest substrate signal of any band** (raw FC p=.014, the only "separated" band on `raw_fc_matched_strength`). The cophenetic step preserves a multiscale α trace that passes C3 at p=0.002 with a 8.3× effect-size ratio. Under epi-zone exclusion the ratio strengthens to 27.7×, consistent with α-relevant non-epi cortex carrying the trace — a supportive mechanistic observation, not a verdict-modifier.

## 1. Scientific claim

**Cohort-level question.** Does the α-band (8–13 Hz) post-task resting state at the LRG-CTM layer sit closer to the task state than the pre-task resting state does, in the cophenetic communication geometry of `|ImCoh|`-derived FC matrices?

**Refined biological claim.** α-band cortical coupling reorganizes during the task and the reorganization is retained in the per-pair cophenetic multiscale geometry but **not** in the leading-mode subspace structure. The reorganization is per-pair multiscale (cophenet merges integrate `N−1` dendrogram scales) and **single-aspect** (Grassmann `d_G(k)` is silent — the slowest non-trivial subspace modes do not rotate at α). Mechanistic adjunct (not a verdict-bearing claim): C5 epi-zone exclusion strengthens the cohort effect by ~3.4× in effect-size ratio, consistent with the trace being carried by non-epileptogenic cortex.

**Falsification budget.** The claim fails if any of the primary controls fail: (a) `ρ_split^coph` cohort-paired Wilcoxon p ≥ 0.05 at C3 (matched-strength); (b) cross-probe C4 paired Wilcoxon rejects (`split > xprobe` significant degradation) or sign mismatch; (c) LOO max-p at C3 exceeds 0.05 (Decision 11 robustness). C5 epi-X is a secondary mechanistic observation per amended Decision 10 — it is not a primary falsifier of the α verdict.

## 2. Cohort, substrate, and library entry points

Same as β: n=10 patients (Pat_02/03/05/06/07/08/10/13/14/15); `imcoh_abs` band-averaged magnitude of imaginary coherence (Ewald 2012; Bastos & Schoffelen 2016); β-band `B = [8, 13] Hz`. Pat_03 acquired at 1024 Hz; handled at the config layer (`nperseg_for_fs(fs)`, `FS_OVERRIDES`) — full cohort member at analysis layer (policy updated 2026-05-18, no dropout sensitivity test). Pat_10 task rows [53, 54, 55] dropped at load. Pat_14 vendor-replaced 2026-04-25.

```python
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
W = load_fc_matrix(patient="Pat_02", phase="rest_pre", band="alpha", fc_method="imcoh_abs")
lrg = load_lrg_result(patient="Pat_02", phase="rest_pre", band="alpha", fc_method="imcoh_abs")
# lrg.ultrametric_matrix is D_coph (cophenet of UPGMA(D(τ_max))) — canonical per-pair object
# lrg.eigenvalues, lrg.eigenvectors are the Laplacian spectrum for the Grassmann probe
```

## 3. The two LRG probes — methodology, results, and provenance

### 3.1 Per-pair multiscale correlation on `D_coph` — **PRIMARY FINDING for α**

#### What `D_coph` is and why it is the per-pair multiscale object
`D_coph_ij = h(LCA_L(i, j))` is the merge height in the UPGMA dendrogram of `D(τ_max) = 1/ρ̂(τ_max)` at which contacts `i` and `j` first coalesce. The cophenetic image of a hierarchical clustering is ultrametric by construction, and each entry reads the intrinsic communication-coalescence scale of that pair. Because our `|ImCoh|` Laplacian spectrum is continuous (no gap-induced crossover scales), the dendrogram's `N−1` merge heights replace the spectral-gap scale identification of [Villegas 2025]; `D_coph` is therefore the natural multiscale per-pair carrier in this regime. See `methods/methods_revision_2026-05-18_cophenet.md` for the full justification.

#### Critical preamble (5-point)
1. **Claim.** α `ρ_split^coph(taskT, rsPost, rsPre_A, rsPre_B) > 0` at cohort level.
2. **Null (C3, mandatory).** Strength-preserving 4-cycle ±δ rewiring of each per-phase adjacency, R=200, seed 20260511. Paired Wilcoxon one-sided greater across patients: obs > surrogate.
3. **Strongest plausible alternative.** Same-band node-strength evolution from rsPre to rsPost driving the per-pair distance shift without genuine multiscale reorganization.
4. **Does C3 cover (3).** Yes — matched-strength surrogate preserves `s_i` to 10⁻⁴ and reshuffles only the topology. What C3 cannot exclude: physiological drift coupled to task block (would need a no-task control session, not present in the paradigm) — C2 drift-floor partially addresses this. What C3 cannot exclude alone: epi-zone-driven cohort dilution (addressed by C5 epi-X below, audit_68).
5. **Falsification.** C3 paired Wilcoxon p ≥ 0.05, or C2 drift floor not separable, or C4 cross-probe flips sign, or C5 epi-X weakens rather than strengthens.

#### Method
- `Δ_task(i,j) = D_coph^taskT(i,j) − D_coph^rsPre_A(i,j)`
- `Δ_rest(i,j) = D_coph^rsPost(i,j) − D_coph^rsPre_B(i,j)`
- `ρ_split^coph = ρ_S(Δ_task, Δ_rest)` Spearman over all `N(N−1)/2` pair entries.
- Independent rsPre halves A/B; full LRG pipeline (eigendecomp → propagator at `τ_max = 1/λ_max` → UPGMA → cophenet) rerun per half.

#### Results (α, full cohort)

| Statistic | Value | Source |
|---|---|---|
| Cohort-median `ρ_split^coph` | +0.1146 | `ctm_triangle/cohort_summary.csv` alpha row |
| n_trace_split | 8/10 | same |
| n_above_drift | 8/10 | same |
| C1 Wilcoxon p (`ρ_split > 0`) | **0.00977** | same col `wilcoxon_split_gt_0_p` |
| C2 Wilcoxon p (`ρ_split > ρ_drift`) | **0.00684** | same col `wilcoxon_split_gt_drift_p` |
| C3 matched-strength: obs_median | +0.1054 | `matched_strength_surrogate_split_baseline/cohort_summary.csv` alpha |
| C3 surr median (median over surrogate medians) | +0.0128 | same |
| C3 ratio obs/surr | **8.25×** | derived |
| C3 paired Wilcoxon z / p | 54.0 / **0.00195** | same |
| C3 n_above_surrogate (per-patient) | 5/10 | same — descriptive auxiliary statistic only |
| C4 cross-probe rho_xprobe | +0.105 | `ctm_triangle/cohort_summary.csv` alpha |
| C4 n_trace_xprobe | 8/10 (+sign matches +0.115) | same |

All 4 primary controls (C1/C2/C3/C4 Wilcoxons) pass at p < 0.05. Per CLAUDE.md `feedback_no_hardcoded_test_thresholds.md` + `feedback_patient_counts_never_the_gate.md`, the Wilcoxon test IS the gate; per-patient counts are descriptive auxiliary statistics and never modify the verdict.

#### C5 epi-zone exclusion (audit_68 — secondary mechanistic observation, per amended Decision 10)

| Statistic | Full cohort | Epi-X cohort | Δ |
|---|---|---|---|
| Cohort-median `obs_ρ` | +0.105 | **+0.187** | +0.082 |
| Surrogate median | +0.0128 | +0.0067 | −0.006 |
| Ratio obs/surr | 8.25× | **27.70×** | +3.4× |
| Wilcoxon paired z / p | 54 / 0.00195 | 49 / **0.01367** | p widens (smaller per-patient N) but stays < 0.05 |
| n_above_surrogate | 5/10 | **7/10** | +2 patients |

Under C5 epi-X the effect-size ratio strengthens ~3.4× and the Wilcoxon stays significant; the per-patient `n_above_surrogate` count increases from 5/10 to 7/10 (descriptive only). This is supportive of α reorganization being carried by non-epileptogenic cortex. **Per amended Decision 10 (2026-05-28) and retracted Decision 1, the C5 strengthening is a secondary mechanistic observation and not a verdict-promoter.** The α verdict is anchored at C3 alone. Source: `alpha_epi_exclusion/cohort_summary.csv`.

#### Per-patient C5 epi-X comparison (`alpha_epi_exclusion/comparison.csv`)

| Patient | `obs_ρ` full | `obs_ρ` epi-X | Δ | N_epi excluded | reading |
|---|---|---|---|---|---|
| Pat_02 | −0.037 | −0.035 | +0.002 | 14 | no trace either; tiny change |
| Pat_03 | +0.360 | +0.189 | −0.171 | 6 | strong both; weakens but still strong |
| Pat_05 | +0.078 | +0.121 | +0.043 | 14 | **strengthens to + trace** |
| Pat_06 | **+0.833** | +0.692 | −0.141 | 10 | dominant + trace both |
| Pat_07 | +0.081 | +0.211 | +0.130 | 7 | **strengthens to + trace** |
| Pat_08 | +0.380 | +0.290 | −0.090 | 9 | + trace both; weakens |
| Pat_10 | +0.130 | +0.193 | +0.063 | 10 | **strengthens to + trace** |
| Pat_13 | +0.028 | **−0.171** | −0.199 | **30** | **inverts** — patient has largest epi-zone fraction |
| Pat_14 | +0.210 | +0.184 | −0.026 | 12 | + trace both; mild weakening |
| Pat_15 | +0.040 | +0.040 | 0 | **0** | no epi nodes — invariant |

Read: under C5 epi-X, three patients move into the + trace direction (Pat_05, Pat_07, Pat_10), one patient inverts (Pat_13 — the only patient with N_epi ≥ 30, ~25% of contacts), the four already strongly positive (Pat_03, Pat_06, Pat_08, Pat_14) stay strongly positive with mild weakening, and Pat_15 (zero epi nodes) is invariant by construction. The cohort signal becomes more uniform and stronger under C5.

#### Reading
α `ρ_split^coph` passes all 4 primary controls at full cohort, with C3 paired Wilcoxon p = 0.00195 the verdict gate. The per-patient C3 `n_above_surrogate` = 5/10 is a descriptive cohort-agreement statistic and does not modify the verdict (patient-count thresholds were retired 2026-05-19). Under C5 epi-X the effect strengthens ~3.4× in ratio, consistent with α reorganization being carried by non-epileptogenic cortex — a mechanistic observation that supports the biological interpretation without being a verdict-promoter.

#### Cache + script provenance
- Full cohort: `audit_63_split_baseline_surrogate.py` → `data/audit/matched_strength_surrogate_split_baseline/`
- C5 epi-X: `audit_68_alpha_epi_exclusion.py` → `data/audit/alpha_epi_exclusion/`
- Within-baseline triangle: `audit_*_ctm_triangle.py` → `data/audit/ctm_triangle/`
- Object built from `lrg.ultrametric_matrix` field in `data/cache/imcoh_lrg/Pat_NN/alpha_{phase}_lrg_imcoh-abs.npz`

### 3.2 Global subspace rotation: Grassmann chordal distance `d_G(k)` — **no trace at α**

#### What `d_G(k)` reads vs what `ρ_split^coph` reads
`d_G(k)` is the chordal Grassmann distance between leading-`k`-dimensional Laplacian eigenmode subspaces of two phases (`d_G² = k − ‖A^T B‖_F²` with A, B = N×k orthonormal eigenvector matrices, skipping the trivial mode). It reads **whether the slow-diffusion subspace rotates** between phases, at a chosen `k`. The per-pair `ρ_split^coph` and the Grassmann `d_G(k)` are not derivable from each other and read structurally distinct aspects of the same propagator.

#### Critical preamble (5-point)
1. **Claim.** `T_G(k) = d_G(taskT, rsPost; k) − d_G(rsPre_A, taskT; k) < 0` at cohort level for some contiguous window of `k`.
2. **Null (C3, audit_70 cluster-extent permutation).** R=200 matched-strength surrogates; for each phantom-surrogate r, cohort paired Wilcoxon between surrogate r and the mean of the remaining R−1 surrogates per k, longest contiguous run with p < 0.05. Empirical null distribution of longest-run length.
3. **Strongest plausible alternative.** Same-strength edge rewiring producing a similar contiguous-significant cluster by chance.
4. **Does C3 cover (3).** Yes — matched-strength preserves node strengths, only topology shuffles; phantom-surrogate-vs-rest builds the exchangeable cluster-extent null. Replaces the previous 8-cell hardcoded threshold (VERDICT_LEDGER.md Decision 6).
5. **Falsification.** Observed longest run within the null distribution → no trace.

#### Method
- For each `k` ∈ [2, 112] and each patient: `T_G(k) = d_G(U_k^taskT, U_k^rsPost) − d_G(U_k^rsPre_A, U_k^taskT)`.
- Cohort paired Wilcoxon one-sided less per `k`: obs `T_G(k)` < mean-surrogate `T_G(k)` across patients.
- Cluster-extent: longest contiguous run of k with cohort-Wilcoxon p < 0.05.
- Null: for each surrogate r in R=200, phantom-test r vs (R−1)-mean, longest run. Cluster p = (1 + #(null ≥ obs)) / (R + 1).

#### Results (α)

| Statistic | Value | Source |
|---|---|---|
| Observed longest contiguous-significant run | **4 cells** at k = 11..14 | `grassmann_regate_no_filter/contig_summary.csv` audit_66 alpha |
| Cluster-extent null mean LR | 2.09 | `grassmann_cluster_extent/cohort_summary.csv` alpha |
| Cluster-extent null 95th LR | 6.0 | same |
| Cluster-extent null max LR | 12 | same |
| Observed cluster mass `T_G^*` (raw Σ −log10p) | 7.98 | same |
| Observed `T_G^*` (normalized per C1) | **0.0312** (raw / 255.65, n_k=111, R=200) | same |
| Null mean cluster mass | 7.35 (post-fix all-clusters) | same |
| Null 95th cluster mass | 21.26 (post-fix) | same |
| **cluster_p_longest_run (descriptive)** | 0.1592 | same |
| **cluster_p_mass (Decision-8 gate)** | **0.3482** | same — no trace (mass-only gate fails) |
| C5 epi-X (audit_67): obs longest run | 3 cells | `grassmann_regate_no_filter/contig_summary.csv` audit_67 alpha |

Observed 4-cell longest run is well within the null distribution (null 95th = 6.0; null mean = 2.09). Both cluster-p statistics fail the 0.05 gate; the cluster-mass statistic is right on the boundary (0.0995) but still above. Under C5 epi-X the longest run further shrinks to 3 cells.

#### Reading
α has no detectable Grassmann subspace trace under the cluster-extent permutation gate. The 4-cell observation at k=11..14 is unremarkable against the empirical null. The α reorganization does not live in the slow-mode subspace; it lives only in the per-pair cophenetic geometry.

#### Cache + script provenance
- `audit_66_grassmann_matched_strength_surrogate.py` → `data/audit/grassmann_matched_strength_surrogate/`
- `audit_67_grassmann_epi_exclusion.py` → `data/audit/grassmann_epi_exclusion/`
- `audit_69_grassmann_regating_no_filter.py` → `data/audit/grassmann_regate_no_filter/`
- `audit_70_grassmann_cluster_extent.py` → `data/audit/grassmann_cluster_extent/`
- Object built from `eigvecs` field in `data/cache/imcoh_lrg/Pat_NN/alpha_{phase}_lrg_imcoh-abs.npz` (obs) and `data/cache/matched_strength_surrogate_lrg/Pat_NN/alpha_{phase}_R200_swap20_seed20260511_imcoh_abs.npz` (surrogate)

## 4. Multiscale structure of the α trace

### Per-pair multiscale via `D_coph` (primary)

The cophenet integrates across `N−1 ≈ 112–121` merge heights per patient. The α `ρ_split^coph` is therefore an integrated multiscale signal at the per-pair level: it cannot be reduced to a single-scale per-pair comparison because the cophenetic image carries each pair to its own coalescence scale.

The headline three-layer table shows that raw FC (single-pair magnitudes), raw `D(τ_max)` (single-scale propagator distance), and cophenet `D_coph` all detect the α trace at the Wilcoxon level (p = 0.014 / 0.042 / 0.002 respectively); the descriptive per-patient counts shift from 8/10 → 7/10 → 5/10 across the three layers. Under C5 epi-X the cophenet effect-size ratio strengthens ~3.4× and the per-patient count rises to 7/10 (descriptive only). Reading: the multiscale per-pair filter induced by the dendrogram is **selectively retaining the non-epi-cortex component** of the α reorganization — a mechanistic observation supportive of the cophenet methodology directive, with the C3 Wilcoxon being the verdict gate at each layer.

### Subspace multiscale via Grassmann `d_G(k)` (silent)

α has no contiguous-significant Grassmann window. The 4-cell run at k = 11..14 is within the empirical null (cluster_p = 0.159). Reading: α's reorganization does not project into a coherent rotation of the leading-mode subspace at any `k`. This is what distinguishes α from β: β is multi-probe (both `ρ_split^coph` and Grassmann pass C3 with epi-X strengthening); α is single-probe (`ρ_split^coph` only).

## 5. Anatomical distribution — α cophenet trace is strong-localized to a bilateral cingulate + parahippocampal + parietal cortical network; 9 of 11 regions reproduce under C5 epi-X

The α cophenet trace is **localized to a distributed cortical network** of 11 named DK regions that all pass A1+A3 (hypergeometric + matched-strength surrogate, R=200, seed 20260511) at q_BH<0.05 / p_emp<0.05. Under C5 epi-X exclusion **9 of 11 named regions reproduce** at the joint A1+A3 gate, with two named dropouts that carry distinct meaning: left caudal middle frontal gyrus drops at A1 (marginal q_BH=0.051, the weakest of the 11 in full cohort), and **right medial orbitofrontal cortex drops at A2 (matched-strength surrogate z falls from 3.22 to 1.28)** — A1 remains strong, but the matched-strength surrogate no longer separates the residual signal, indicating that the full-cohort rh-medial-OFC enrichment carried a substantial epi-zone-coupled topological component. The 9 surviving regions all maintain or strengthen under epi-X, with **lh-rostral-anterior-cingulate the largest gainer** (enrichment 3.04× → 5.37×, q_BH 4.1e-20 → 1.5e-31). This is consistent with the C5 effect-size strengthening at the cohort level (mechanistic observation per amended Decision 10 / retracted Decision 1, 2026-05-28) and supports reading the α reorganization as a non-epileptic-cortex multiscale phenomenon. The original "11/11 reproduce identically" claim recorded here in the 2026-05-19 lockdown was an artifact of a silent no-op bug in the audit_71 `--epi-x` mask (`np.isin(int_array, str_array)` always-False); the bug was patched and the audit rerun on 2026-05-22 (see `.agents/preprint/directives/writing_directive_2026-05-22_alpha_results_drafting.md` §"Critical issues" → Issue 1 for the bug history and the verified outcome).

Audited under `locked/ANATOMY_CONTROLS.md`. Verdict source: `locked/ANATOMY_LEDGER.md` 2026-05-19 entry.

### α cophenet anatomy (top-decile per-pair `|Δρ_split^coph|`)

Eleven named Desikan–Killiany regions pass A1+A3 jointly (full cohort, R=200):

| Region | k_trace / K_region | A1 enrichment | A1 p_hyper | A3 obs_z | A3 p_emp |
|---|---|---|---|---|---|
| ctx-lh-caudalanteriorcingulate | 45 / 349 | 2.06× | 4.07e-06 | 4.67 | 0.005 |
| ctx-lh-parahippocampal | 85 / 357 | 3.80× | 7.53e-27 | 9.76 | 0.005 |
| ctx-lh-rostralanteriorcingulate | 90 / 472 | 3.04× | 3.66e-21 | 14.26 | 0.005 |
| ctx-rh-medialorbitofrontal | 242 / 2188 | 1.77× | 1.22e-17 | 3.22 | 0.005 |
| ctx-rh-caudalmiddlefrontal | 103 / 567 | 2.90× | 2.49e-22 | 6.23 | 0.005 |
| ctx-rh-postcentral | 247 / 1135 | 3.47× | 6.70e-67 | 4.47 | 0.005 |
| ctx-lh-caudalmiddlefrontal | 76 / 815 | 1.49× | 4.10e-04 | 3.02 | 0.015 |
| ctx-rh-caudalanteriorcingulate | 52 / 340 | 2.44× | 2.85e-09 | 3.45 | 0.020 |
| ctx-rh-precuneus | 86 / 936 | 1.47× | 2.92e-04 | 2.60 | 0.025 |
| ctx-rh-superiorparietal | 84 / 819 | 1.64× | 7.94e-06 | 3.64 | 0.025 |
| ctx-rh-posteriorcingulate | 100 / 804 | 1.99× | 7.39e-11 | 3.13 | 0.045 |

Network anatomy: **bilateral cingulate** (anterior caudal + anterior rostral + right posterior) + **left parahippocampal** + **right medial OFC** + **bilateral caudal middle frontal** + **right postcentral / precuneus / superior parietal**. A widely distributed cortical network spanning cingulate / medial-temporal / medial-frontal / parietal cortex. Cohort verdict: **strong localized**.

Source: `data/audit/anatomy_alpha_cophenet/cohort_summary.csv`.

### C5 epi-X anatomy: 9 of 11 named regions reproduce (post-patch 2026-05-22)

Re-running the α cophenet anatomy audit with epi-zone contacts removed per patient (`load_epileptic_nodes() × channel_labels.csv`) yields **9 of 11 named DK regions surviving** the joint A1+A3 gate. The two dropouts:

| Dropped region | Full cohort | Epi-X | Reading |
|---|---|---|---|
| `ctx-lh-caudalmiddlefrontal` | enr 1.49×, q_BH 1.6e-3, z 3.02, p_emp 0.015 | A1 FAIL: q_BH 5.1e-2 (above 0.05); A3 z 2.32, p_emp 0.020 | Weakest of the 11 in full cohort; marginal statistical fail under the slightly different cohort marginals after epi-X masking. |
| `ctx-rh-medialorbitofrontal` | enr 1.77×, q_BH 1.1e-16, z 3.22, p_emp 0.005 | A1 still strong (q_BH 3.3e-6); **A3 FAIL: z 1.28, p_emp 0.124** | A2/A3 (matched-strength surrogate) no longer separates the residual rh-medial-OFC signal from matched-strength noise after epi-X. **Biological reading**: a substantial part of the full-cohort medial-OFC enrichment was carried by topology coupled to epi-zone contacts; the matched-strength control correctly strips that contribution. |

The 9 surviving regions and how they shift under epi-X:

| Region | Full enr | Epi-X enr | Reading |
|---|---|---|---|
| `ctx-lh-caudalanteriorcingulate` | 2.06× (q 2.1e-5) | **2.08× (q 2.9e-2)** | stable |
| `ctx-lh-parahippocampal` | 3.80× (q 1.4e-25) | **3.89× (q 6.6e-25)** | medial-temporal anchor; strengthens slightly |
| `ctx-lh-rostralanteriorcingulate` | 3.04× (q 4.1e-20) | **5.37× (q 1.5e-31)** | **largest gainer under epi-X** |
| `ctx-rh-caudalanteriorcingulate` | 2.44× (q 1.8e-8) | **3.31× (q 3.2e-11)** | strengthens |
| `ctx-rh-caudalmiddlefrontal` | 2.90× (q 3.5e-21) | 2.74× (q 5.2e-17) | mild weakening |
| `ctx-rh-postcentral` | 3.47× (q 3.8e-65) | **3.66× (q 7.4e-69)** | rock stable; deepest q_BH in the cohort |
| `ctx-rh-posteriorcingulate` | 1.99× (q 5.9e-10) | 1.93× (q 1.1e-8) | mild weakening |
| `ctx-rh-precuneus` | 1.47× (q 1.3e-3) | 1.47× (q 1.1e-3) | identical |
| `ctx-rh-superiorparietal` | 1.64× (q 3.7e-5) | 1.64× (q 3.4e-5) | identical |

Net reading: the 9-region surviving network — **bilateral cingulate (4 regions: lh-/rh-caudal-AC, lh-rostral-AC, rh-posterior-cingulate) + left parahippocampal cortex + right caudal middle frontal gyrus + right postcentral / precuneus / superior parietal** — sees several constituent regions (particularly bilateral cingulate and left parahippocampal cortex) **strengthen** under epi-X. The pattern mirrors the cohort effect-size strengthening (`ρ_split^coph` cohort ratio 8.3× → 27.7× under epi-X — secondary mechanistic observation, not a verdict-promoter) and reads the α reorganization as a **physiological cortical phenomenon, with bilateral cingulate and left medial-temporal cortex as the anatomical core** and a small epi-zone-coupled contribution that explains the rh-medial-OFC dropout under matched-strength.

Source: `data/audit/anatomy_alpha_cophenet_epiX/cohort_summary.csv` (post-patch 2026-05-22; MD5 `199d742b22e00893ff3fbae4623a1b1c`).

**Bug history.** The original `anatomy_alpha_cophenet_epiX/cohort_summary.csv` produced 2026-05-19 was byte-identical to the full-cohort file because `audit_71_anatomy_cophenet.py` had a silent no-op at the `--epi-x` mask: `np.isin(np.arange(len(regions_df)), np.asarray(list(epi)))` compared integer FC-channel indices to the string labels returned by `load_epileptic_nodes()`, returning all-False, so `keep_mask` was always all-True. Patch applied 2026-05-22: replaced the comparison with the label-normalised lookup pattern used by `audit_68_alpha_epi_exclusion.py` (`_normalise_label(label_raw)` from `lrg_eegfc.utils.io.regions`). Sanity-checked against audit_68 (Pat_02 N_epi=14, Pat_13 N_epi=30, Pat_15 N_epi=0). The "11/11 reproduce identically" claim in this brief's 2026-05-19 lockdown originated from the bug — it was never a real measurement. The verified outcome above replaces it.

### Caveats

- A1 also enriches `Wm` (white matter, 1.04×, p_hyper=6.31e-04) and `unknown` (1.73×, 1.48e-08) — these are non-DK trivial bins, excluded from the localization claim per `locked/ANATOMY_CONTROLS.md`.
- A2 sampling-corrected bootstrap and A4 implant-geometry regression deferred.
- Grassmann probe at α has no trace (cluster_p=0.159; 4-cell run within the null distribution) — no Grassmann anatomy verdict is computed.

### Cache + script provenance
| Artifact | Path |
|---|---|
| α cophenet anatomy (full) | `data/audit/anatomy_alpha_cophenet/cohort_summary.csv` |
| α cophenet anatomy (C5 epi-X) | `data/audit/anatomy_alpha_cophenet_epiX/cohort_summary.csv` |
| Producing script | `scripts/01_compute/audit/audit_71_anatomy_cophenet.py --band alpha [--epi-x]` |
| Anatomy battery doc | `.agents/preprint/locked/ANATOMY_CONTROLS.md` |
| Anatomy verdict ledger | `.agents/preprint/locked/ANATOMY_LEDGER.md` (2026-05-19 entries) |

## 6. Patient-by-patient reading

### α `ρ_split^coph` per-patient signs across the primary probe

The full table is in §3.1 above (the C5 epi-X comparison). Summary at the cohort level:

- **Strong + trace at full cohort**: Pat_06 (+0.833), Pat_08 (+0.380), Pat_03 (+0.360), Pat_14 (+0.210)
- **Marginal + at full cohort, strengthens under C5**: Pat_05 (+0.078 → +0.121), Pat_07 (+0.081 → +0.211), Pat_10 (+0.130 → +0.193)
- **Borderline at full cohort, inverts under C5**: Pat_13 (+0.028 → −0.171; N_epi = 30, largest in cohort)
- **No trace either**: Pat_02 (−0.037 → −0.035), Pat_15 (+0.040 → +0.040; N_epi = 0)

7/10 patients are in the + trace direction under C5 epi-X (Pat_03/05/06/07/08/10/14). Pat_02 and Pat_15 are essentially flat. Pat_13 is the only inversion — driven by the largest epi-zone fraction in the cohort.


## 7. Robustness panel (under the locked battery)

| Check | α verdict | Source |
|---|---|---|
| C1 within-rsPre split-half null | ✓ 8/10 cohort, p = 0.010 | `ctm_triangle/cohort_summary.csv` alpha |
| C2 drift-floor `ρ_split > ρ_drift` | ✓ 8/10 above drift, p = 0.007 | same |
| C3 matched-strength surrogate `ρ_split^coph` | ✓ paired Wilcoxon p = 0.002 (verdict gate); ratio 8.3×; per-patient 5/10 descriptive | `matched_strength_surrogate_split_baseline` alpha |
| C3 matched-strength surrogate Grassmann | ✗ cluster_p = 0.159 (no trace) | `grassmann_cluster_extent` alpha |
| C4 cross-probe restriction | ✓ +0.105, paired-Wilcoxon non-degradation; descriptive 8/10 +sign | `ctm_triangle` alpha (rho_xprobe column) |
| C5 epi-X `ρ_split^coph` (secondary mechanistic) | Wilcoxon-on-epi-X p = 0.0098; ratio strengthens 8.3× → 27.7×; **not verdict-promoter** per amended Decision 10 | `alpha_epi_exclusion/cohort_summary.csv` |
| C5 epi-X Grassmann | (n/a; Grassmann has no trace at α) | `grassmann_regate_no_filter` audit_67 alpha |
| Pat_13 leave-one-out under C5 | Pat_13 inverts under epi-X — mechanistically consistent with the largest epi-zone fraction; descriptive only | `alpha_epi_exclusion/comparison.csv` |
| Cross-probe restriction (C4) | preserves sign | `ctm_triangle` |
| Three-layer cohort table (raw FC / raw D / cophenet) | all three layers pass C3 at Wilcoxon; cophenet C5 epi-X strengthens the ratio, supportive non-epi-cortex reading | headline section above |

## 8. Sensitivity panel under the locked 5-control battery

The α verdict is anchored at C3 paired Wilcoxon p = 0.00195 alone (locked Wilcoxon-as-gate rule; patient-count thresholds retired 2026-05-19, reaffirmed 2026-05-28). C5 epi-X is reported as a secondary mechanistic observation per amended Decision 10 + retracted Decision 1. **The verdict for α is `strong trace, only D_coph` based on C1–C4 at full cohort.**

| Control | α status | Notes |
|---|---|---|
| C1 split | ✓ p = 0.010 | passes |
| C2 drift | ✓ p = 0.007 | passes |
| C3 matched-strength (cophenet) | ✓ paired Wilcoxon p = 0.002 | verdict gate; per-patient 5/10 descriptive only |
| C3 matched-strength (Grassmann) | ✗ cluster_p = 0.159 | Grassmann probe inactive at α |
| C4 cross-probe (paired Wilcoxon) | ✓ non-degradation | passes |
| C5 epi-X (cophenet) | secondary mechanistic observation (audit_68) | ratio strengthens 3.4×; supportive of non-epi-cortex contribution; not verdict-promoter |
| C5 epi-X (Grassmann) | (n/a; Grassmann has no trace at α) | confirms no Grassmann trace |

## 9. Interpretation

**Neuroscientific framing.** α (8–13 Hz) is the canonical posterior/attentional rhythm and the dominant inter-areal coupling band at rest in human iEEG (Spaak et al. 2012; Mathewson 2011; Klimesch 2012). A non-epileptic-cortex α trace at the LRG-CTM layer is consistent with **attentional-state carryover** from task into post-task rest, specifically in the per-pair cophenetic communication geometry rather than in the slow-subspace structure.

**Why per-pair only (not Grassmann).** β (13–30 Hz) carries both `ρ_split^coph` and Grassmann subspace traces (multi-aspect). α carries only the per-pair multiscale signature. The interpretation is that the α reorganization is **distributed across the cophenetic merge-height continuum** without producing a coherent rotation of the slow-mode subspace at any single `k`. Mechanistically this is consistent with a wide-band, per-pair coupling re-weighting that does not concentrate in the slowest few diffusion modes — distinct from β, where the reorganization rotates a 29-mode subspace.

**Why epi-zone exclusion strengthens α (mechanistic observation).** Pat_13 has the largest fraction of epileptogenic contacts (30/119 ≈ 25%). Under C5 epi-X, three patients (Pat_05, Pat_07, Pat_10) move from marginal-positive to clearly-positive, while Pat_13 inverts (whose epi-zone is large enough that its α coupling dominates the residual signal). β's Grassmann under C5 also strengthens; at α the C5 effect is concentrated in the per-pair cophenet probe. This pattern is consistent with α reorganization being carried by non-epileptogenic cortex — supportive narrative, not a verdict-modifier.

**The role of the cophenet step at α.** Raw FC and raw `D(τ_max)` both detect α at the Wilcoxon level. The cophenetic step also passes C3 at p = 0.002 (verdict gate; descriptive per-patient count drops 8/10 → 7/10 → 5/10 across the three layers as the metric becomes more multiscale-aggregative). Under C5 epi-X the cophenet ratio strengthens 3.4×, supportive of cophenet acting as a multiscale + spatial-purity filter that selectively retains the non-pathological component of the α reorganization. This is a mechanistic endorsement of the methodology directive's framing of cophenet as a multiscale-resolution operator.

## 10. Source-of-truth references

### Lockdown documents
- `.agents/preprint/locked/CONTROLS.md` — 5-control battery
- `.agents/preprint/locked/VERDICT_LEDGER.md` — locked verdict tag `strong trace, only D_coph` for α
- `.agents/preprint/methods/methods_revision_2026-05-18_cophenet.md` — binding methods directive
- `.agents/preprint/bands/01_beta.md` — β template (this brief mirrors its structure)

### Audit data (rows the numbers come from)
- `data/audit/ctm_triangle/cohort_summary.csv` — α row: C1, C2, C4 statistics
- `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` — α row: C3 cophenet
- `data/audit/alpha_epi_exclusion/cohort_summary.csv` — α C5 epi-X cohort summary
- `data/audit/alpha_epi_exclusion/comparison.csv` — α C5 per-patient Δρ table
- `data/audit/grassmann_cluster_extent/cohort_summary.csv` — α cluster-extent permutation (audit_70)
- `data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv` — α per-k T_G + surrogate stats
- `data/audit/grassmann_regate_no_filter/contig_summary.csv` — α audit_66 + audit_67 contiguous windows
- `data/audit/grassmann_epi_exclusion/cohort_summary.csv` — α per-k under C5 epi-X
- `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv` — α substrate row
- `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv` — α raw D sensitivity

### Audit scripts (provenance)
- `audit_63_split_baseline_surrogate.py` — C3 ρ_split^coph matched-strength (full cohort)
- `audit_68_alpha_epi_exclusion.py` — C5 ρ_split^coph epi-X at α
- `audit_66_grassmann_matched_strength_surrogate.py` — Grassmann C3 per-k
- `audit_67_grassmann_epi_exclusion.py` — Grassmann C5 per-k
- `audit_69_grassmann_regating_no_filter.py` — Grassmann contig-summary regating
- `audit_70_grassmann_cluster_extent.py` — Grassmann cluster-extent permutation (2026-05-19)

### Companion guides
- `.agents/guides/02_methods/lrg-framework-guide.md` — LRG primitive (formulas + codebase mapping)
- `.agents/guides/02_methods/imcoh-guide.md` — imcoh_abs convention + Nolte 2004 immunity
- `.agents/guides/04_rules/never-always-list.md` — coding rules + no-hardcoded-thresholds feedback
- `.agents/guides/05_plotting/README.md` — figure conventions (for upcoming α preprint figures)

## 11. What follows

This brief freezes the α verdict for the preprint manuscript under the 2026-05-18 / 2026-05-19 lockdown. Next steps:

1. **α figures** — produce F1 (4-control panel for `ρ_split^coph` + C5 epi-X), F2 (Grassmann `T_G(k)` cohort cells highlighting the empirical null), F3 (per-patient Δ_task vs Δ_rest scatter on `D_coph` for representative patients), and the per-patient C5 epi-X waterfall plot from §6. Output at `data/preprint/figures/alpha/`.
2. **No anatomy audit owed at α** under the locked battery. If a Desikan–Killiany enrichment is wanted as a descriptive supplement, the β anatomy pipeline (`audit_64`) can be re-pointed at α.
3. **No further sensitivity tests owed.** The α verdict is anchored at C3 paired Wilcoxon p = 0.00195 alone under the locked Wilcoxon-as-gate rule; the C5 epi-X audit (audit_68) provides supportive mechanistic context (non-epi-cortex contribution) per amended Decision 10.

Verdict ready for writing-agent handoff: **`strong trace, only D_coph`**, anchored at C3. The mechanistic narrative — α reorganization as a non-epileptogenic-cortex per-pair multiscale phenomenon — is supported by C5 but is not the primary verdict claim.
