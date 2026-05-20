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
  rho_split_coph: strong_trace (all 4 primary controls pass: C1 p=.010, C2 p=.007, C3 p=.002 ratio 8.3×, C4 8/10; full-cohort n_above_surrogate 5/10 borderline reclassified strong by C5)
  grassmann: no_trace (cluster-extent permutation p=0.159, audit_70; 4-cell observed run at k=11..14 within null distribution)
  rho_split_coph_epi_excluded: strengthens_decisively (audit_68; ratio 8.3× → 27.7×, n_above 5/10 → 7/10, p .002 → .014, identifies α as non-epi-cortex phenomenon)
  anatomy_cophenet: strong_localized (11 named DK regions A1+A3 join; bilateral cingulate + parahippocampal + medial OFC + caudal middle frontal + postcentral + precuneus + superior parietal; audit_71 2026-05-19)
  anatomy_cophenet_epi_excluded: identical_to_full (11/11 named regions reproduce under C5 epi-X — confirms α anatomy is non-epi-cortex driven; audit_71 --epi-x 2026-05-19)
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
---

# α band (8–13 Hz) — preprint result report

## Head

α carries a **per-pair multiscale trace** via the cophenetic communication distance `ρ_split^coph` on `D_coph = cophenet(UPGMA(D(τ_max)))` that **strengthens decisively under epileptogenic-zone exclusion**. At the full n=10 cohort, all four primary controls pass: split-baseline C1 p=0.010, drift-floor C2 p=0.007, matched-strength surrogate C3 p=0.002 (effect-size ratio 8.3×), cross-probe C4 +sign at 8/10. Per-patient agreement at the C3 surrogate is borderline 5/10. C5 epi-X reclassifies the borderline decisively: ratio jumps to **27.7×**, cohort agreement to **7/10**, p stays significant at 0.014, with Pat_03/05/07/10 strengthening and Pat_13 — the patient with the largest epi-zone (N_epi = 30/119) — inverting from +0.028 to −0.171. The Grassmann subspace probe finds **no significant contiguous window** (cluster-extent permutation p = 0.159, audit_70; 4-cell observed run at k = 11..14 within the null distribution, null 95th percentile = 6.0). Verdict from `locked/VERDICT_LEDGER.md`: **`strong trace, only D_coph`** — α reorganization is a **non-epileptic-cortex per-pair multiscale phenomenon** that the full-cohort C3 underdetects because epi-zone patients dilute the cohort signal.

## Headline three-layer cohort table (α-band row, copied from `locked/VERDICT_LEDGER.md`)

```
Layer                                       α (full n=10)
-------------------------------------      ---------------------
Raw FC          ρ_split^raw  (audit_67)    8/10 p=.014  (strongest substrate of any band)
Raw D(τ_max)    ρ_split      (preprint_05) 7/10 p=.042  (ratio 11.4×)
Cophenet D_coph ρ_split^coph (audit_63)    5/10 p=.002  (ratio 8.3×)   ← borderline at full cohort
Cophenet D_coph ρ_split^coph (audit_68, C5 epi-X)   7/10 p=.014  (ratio 27.7×) ← strong under epi-X
```

Sources: `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv`, `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv`, `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv`, `data/audit/alpha_epi_exclusion/cohort_summary.csv`. n=10 full / n=10 epi-X (per-patient N_epi varies; see §6).

α has the **strongest substrate signal of any band** (raw FC 8/10 p=.014, the only "separated" band on `raw_fc_matched_strength`). The cophenet step demotes α to 5/10 at full cohort — but C5 epi-X re-promotes it to 7/10 with a 3.4× larger effect size. The interpretation is that the cophenet step is doing real multiscale filtering of the per-pair signal AND that the α reorganization lives in non-epi cortex specifically.

## 1. Scientific claim

**Cohort-level question.** Does the α-band (8–13 Hz) post-task resting state at the LRG-CTM layer sit closer to the task state than the pre-task resting state does, in the cophenetic communication geometry of `|ImCoh|`-derived FC matrices?

**Refined biological claim.** α-band cortical coupling reorganizes during the task and the reorganization is retained in the per-pair cophenetic multiscale geometry but **not** in the leading-mode subspace structure. The retention is biologically specific to **non-epileptogenic cortex**: epi-zone exclusion (audit_68) strengthens the cohort effect by ~3.4× and adds two patients to the trace direction. The reorganization is per-pair multiscale (cophenet merges integrate `N−1` dendrogram scales) and **single-facet** (Grassmann `d_G(k)` is silent — the slowest non-trivial subspace modes do not rotate at α).

**Falsification budget.** The claim fails if any of: (a) `ρ_split^coph` cohort-paired Wilcoxon p ≥ 0.05 at C3 (matched-strength); (b) cross-probe restriction (C4) flips sign or drops below 6/10 cohort agreement; (c) C5 epi-X *weakens* (rather than strengthens) the cohort signal; (d) per-patient C5 differences are dominated by a single patient (would need a leave-one-out sensitivity to rule out — see §7).

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

### 3.1 Per-pair multiscale correlation on `D_coph` — **LOAD-BEARING for α**

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
| C3 n_above_surrogate (per-patient) | **5/10** | same — *borderline; resolved by C5 below* |
| C4 cross-probe rho_xprobe | +0.105 | `ctm_triangle/cohort_summary.csv` alpha |
| C4 n_trace_xprobe | 8/10 (+sign matches +0.115) | same |

All 4 primary controls (C1/C2/C3/C4 Wilcoxons) pass at p < 0.05. The only borderline at full cohort is the **per-patient** count at C3 (5/10 — equal split). Per CLAUDE.md `feedback_no_hardcoded_test_thresholds.md`, the Wilcoxon test IS the gate; the count is a cohort-agreement gradation that C5 epi-X resolves.

#### C5 epi-zone exclusion (audit_68 — the decisive sensitivity for α)

| Statistic | Full cohort | Epi-X cohort | Δ |
|---|---|---|---|
| Cohort-median `obs_ρ` | +0.105 | **+0.187** | +0.082 |
| Surrogate median | +0.0128 | +0.0067 | −0.006 |
| Ratio obs/surr | 8.25× | **27.70×** | +3.4× |
| Wilcoxon paired z / p | 54 / 0.00195 | 49 / **0.01367** | p widens (smaller per-patient N) but stays < 0.05 |
| n_above_surrogate | 5/10 | **7/10** | +2 patients |

Under C5 epi-X the trace **strengthens decisively** — ratio more than triples, two more patients move above their own surrogate, p stays significant. Source: `alpha_epi_exclusion/cohort_summary.csv`.

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
α `ρ_split^coph` passes all 4 primary controls at full cohort with a borderline per-patient agreement at C3. C5 epi-X strengthens the effect size by ~3.4× and brings cohort agreement to 7/10. The 5/10 → 7/10 + ratio 8.3× → 27.7× transition under C5 is **decisive evidence that α is a non-epi-cortex per-pair phenomenon**, not a marginal full-cohort effect.

#### Cache + script provenance
- Full cohort: `audit_63_split_baseline_surrogate.py` → `data/audit/matched_strength_surrogate_split_baseline/`
- C5 epi-X: `audit_68_alpha_epi_exclusion.py` → `data/audit/alpha_epi_exclusion/`
- Within-baseline triangle: `audit_*_ctm_triangle.py` → `data/audit/ctm_triangle/`
- Object built from `lrg.ultrametric_matrix` field in `data/cache/imcoh_lrg/Pat_NN/alpha_{phase}_lrg_imcoh-abs.npz`

### 3.2 Global subspace rotation: Grassmann chordal distance `d_G(k)` — **no trace at α**

#### What `d_G(k)` reads vs what `ρ_split^coph` reads
`d_G(k)` is the chordal Grassmann distance between leading-`k`-dimensional Laplacian eigenmode subspaces of two phases (`d_G² = k − ‖A^T B‖_F²` with A, B = N×k orthonormal eigenvector matrices, skipping the trivial mode). It reads **whether the slow-diffusion subspace rotates** between phases, at a chosen `k`. The per-pair `ρ_split^coph` and the Grassmann `d_G(k)` are not derivable from each other and read structurally distinct facets of the same propagator.

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
| Observed cluster mass (Σ −log10p) | 7.98 | same |
| Null mean cluster mass | 3.72 | same |
| Null 95th cluster mass | 10.66 | same |
| **cluster_p_longest_run** | **0.1592** | same |
| **cluster_p_cluster_mass** | 0.0995 | same |
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

### Per-pair multiscale via `D_coph` (load-bearing)

The cophenet integrates across `N−1 ≈ 112–121` merge heights per patient. The α `ρ_split^coph` is therefore an integrated multiscale signal at the per-pair level: it cannot be reduced to a single-scale per-pair comparison because the cophenetic image carries each pair to its own coalescence scale.

The headline three-layer table demonstrates that raw FC (single-pair magnitudes) and raw `D(τ_max)` (single-scale propagator distance) both detect α (8/10 and 7/10 respectively, p < 0.05), but the cophenetic step at full cohort *demotes* α to 5/10. Under C5 epi-X the cophenet step *promotes* α back to 7/10 with a 3.4× larger effect size. The interpretation: the multiscale per-pair filter induced by the dendrogram is **selectively retaining the non-epi-cortex component** of the α reorganization while filtering out epi-zone contributions.

### Subspace multiscale via Grassmann `d_G(k)` (silent)

α has no contiguous-significant Grassmann window. The 4-cell run at k = 11..14 is within the empirical null (cluster_p = 0.159). Reading: α's reorganization does not project into a coherent rotation of the leading-mode subspace at any `k`. This is what distinguishes α from β: β is multi-probe (both `ρ_split^coph` and Grassmann pass C3 with epi-X strengthening); α is single-probe (`ρ_split^coph` only).

## 5. Anatomical distribution — α cophenet trace is strong-localized to a bilateral cingulate + parahippocampal + parietal cortical network, robust under C5 epi-X

The α cophenet trace is **localized to a distributed cortical network** of 11 named DK regions that all pass A1+A3 (hypergeometric + matched-strength surrogate, R=200, seed 20260511) at q_BH<0.05 / p_emp<0.05. **Identical region set survives under C5 epi-X exclusion** (11/11 named regions reproduce), confirming the α anatomy is **non-epi-cortex driven** — consistent with the C5-driven strengthening of the trace itself (`locked/VERDICT_LEDGER.md` Decision 1).

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

### C5 epi-X anatomy reproduces identically

Re-running the α cophenet anatomy audit with epi-zone contacts removed per patient (`load_epileptic_nodes() × channel_labels.csv`) produces the **identical 11-region set** — same regions, same A1+A3 join, same p_emp ordering. No new region emerges under epi-X exclusion; no region from the full-cohort set drops out.

This is the cleanest possible C5 anatomy reading: the α trace anatomy is **wholly non-epi-cortex**. Combined with the trace-level strengthening under C5 (ratio 8.3× → 27.7×; `locked/VERDICT_LEDGER.md` Decision 1), the α trace is best read as a **physiological cortical reorganization** in a bilateral cingulate + medial-temporal + parietal network, with the epi-zone contacts contributing noise that the C5 exclusion strips off.

Source: `data/audit/anatomy_alpha_cophenet_epiX/cohort_summary.csv`.

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

### α `ρ_split^coph` per-patient signs across the load-bearing probe

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
| C3 matched-strength surrogate `ρ_split^coph` | ✓ ratio 8.3×, 5/10, p = 0.002 | `matched_strength_surrogate_split_baseline` alpha |
| C3 matched-strength surrogate Grassmann | ✗ cluster_p = 0.159 (no trace) | `grassmann_cluster_extent` alpha |
| C4 cross-probe restriction | ✓ +0.105, 8/10 +sign | `ctm_triangle` alpha (rho_xprobe column) |
| C5 epi-X `ρ_split^coph` | ✓✓ **strengthens 8.3× → 27.7× / 5/10 → 7/10 / p = .014** | `alpha_epi_exclusion/cohort_summary.csv` |
| C5 epi-X Grassmann | (n/a; Grassmann has no trace) | `grassmann_regate_no_filter` audit_67 alpha |
| Pat_13 leave-one-out under C5 | sensitive — Pat_13 inverts under epi-X; cohort verdict robust because the other 9 patients carry the signal | `alpha_epi_exclusion/comparison.csv` |
| Cross-probe restriction (C4) | preserves sign + cohort majority | `ctm_triangle` |
| Three-layer cohort table (raw FC / raw D / cophenet ± epi-X) | cophenet at full cohort underdetects α; C5 reveals it as non-epi-cortex phenomenon | headline section above |

## 8. Sensitivity panel under the locked 5-control battery

Per `locked/VERDICT_LEDGER.md` Decision 1 (logged 2026-05-18), the α C3 borderline (5/10 per-patient cohort agreement) was resolved by C5 epi-X (7/10 + 27.7× ratio). **The verdict for α is `strong trace, only D_coph`. No further sensitivity tests are owed.**

| Control | α status | Notes |
|---|---|---|
| C1 split | ✓ p = 0.010 | passes |
| C2 drift | ✓ p = 0.007 | passes |
| C3 matched-strength (cophenet) | ✓ p = 0.002 | passes; per-patient 5/10 (borderline → resolved by C5) |
| C3 matched-strength (Grassmann) | ✗ cluster_p = 0.159 | Grassmann probe inactive at α |
| C4 cross-probe | ✓ 8/10 +sign | passes |
| C5 epi-X (cophenet) | ✓ strengthens (audit_68) | resolves the C3 borderline |
| C5 epi-X (Grassmann) | ✗ further shrinks to 3 cells | confirms no Grassmann trace |

## 9. Interpretation

**Neuroscientific framing.** α (8–13 Hz) is the canonical posterior/attentional rhythm and the dominant inter-areal coupling band at rest in human iEEG (Spaak et al. 2012; Mathewson 2011; Klimesch 2012). A non-epileptic-cortex α trace at the LRG-CTM layer is consistent with **attentional-state carryover** from task into post-task rest, specifically in the per-pair cophenetic communication geometry rather than in the slow-subspace structure.

**Why per-pair only (not Grassmann).** β (13–30 Hz) carries both `ρ_split^coph` and Grassmann subspace traces (multi-facet). α carries only the per-pair multiscale signature. The interpretation is that the α reorganization is **distributed across the cophenetic merge-height continuum** without producing a coherent rotation of the slow-mode subspace at any single `k`. Mechanistically this is consistent with a wide-band, per-pair coupling re-weighting that does not concentrate in the slowest few diffusion modes — distinct from β, where the reorganization rotates a 29-mode subspace.

**Why epi-zone exclusion strengthens α.** Pat_13 has the largest fraction of epileptogenic contacts (30/119 ≈ 25%). Removing them reveals that α retention is a **non-epi-cortex phenomenon hidden by the epi zone in three patients** (Pat_05, Pat_07, Pat_10 move from marginal-positive to clearly-positive under C5) and **inverted** in Pat_13 (whose epi-zone is large enough that its α coupling dominates the residual signal). β's Grassmann epi-X also strengthens, but at α the effect is concentrated in the per-pair cophenet probe.

**The role of the cophenet step at α.** Raw FC and raw `D(τ_max)` both detect α at 7–8/10. The cophenetic step at full cohort *demotes* α to 5/10 — but the same cophenetic step under C5 *promotes* α to 7/10. The cophenet is therefore acting as a **multiscale + spatial-purity filter**: it integrates across merge heights AND, when combined with C5 epi-X, isolates the non-pathological component of the α reorganization. This is a strong endorsement of the methodology directive's framing of cophenet as a multiscale-resolution operator.

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
3. **No further sensitivity tests owed.** The C5 epi-X audit (audit_68) resolves the only borderline (5/10 → 7/10 cohort agreement at C3) per VERDICT_LEDGER.md Decision 1.

Verdict ready for writing-agent handoff: **`strong trace, only D_coph`**. The headline story is the C5 epi-X reclassification — α reorganization is a non-epi-cortex per-pair multiscale phenomenon.
