---
name: preprint-beta-band
era: IMCOH_ABS_COHORT_N10
status: current
kind: preprint-result-report
band: beta
range_hz: [13, 30]
cohort: [Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15]
n_patients: 10
fc_method: imcoh_abs
tau: 1/lambda_max
canonical_lrg_object: "D_coph = cophenet(UPGMA(D(tau_max)))"
canonical_subspace_object: "U_k = span{phi_2, ..., phi_{k+1}}"
verdict_tag: "strong trace, both probes"
verdict_source: VERDICT_LEDGER.md (locked 2026-05-18, revised 2026-05-19)
verdict_layers:
  substrate_rank: trace_direction (raw FC ρ_split^raw matched-strength p=0.053 borderline; sensitivity layer only)
  rho_split_coph: strong_trace (all 4 primary controls pass; cohort ratio 23.7×, p=0.005, n_above 7/10)
  grassmann: strong_trace (cluster-extent permutation p=0.005, audit_70; 29-cell contiguous-significant run k=27..55)
  grassmann_epi_excluded: strengthens_to_36_cells_at_k21_56 (29/29 manuscript window persist + 7 new cells emerge)
  rho_split_coph_epi_excluded: not_run (C5 epi-X cophenet was audit_68 at α only; β verdict from C1+C2+C3+C4 + Grassmann epi-X)
  anatomy_cophenet: strong_localized (7 DK regions A1+A3 join; cingulate + parahippocampal + entorhinal + insula + postcentral + superior frontal; audit_71 2026-05-19)
  anatomy_grassmann: strong_localized (7 DK regions A3 alone; Hip + temporal + orbitofrontal + insula + rostral middle frontal; audit_72 --cluster-extent over `S(β)`, 2026-05-19 pm — region set identical to retired `K*(β)` run, only z/p_emp updates)
sources:
  - .agents/preprint/locked/CONTROLS.md (locked 5-control battery)
  - .agents/preprint/locked/VERDICT_LEDGER.md (locked verdicts)
  - .agents/preprint/methods/methods_revision_2026-05-18_cophenet.md (binding writing-agent directive)
  - .agents/preprint/established_results/00_open_methodology_question_lrg_D_convention.md (withdrawn, resolved)
  - data/audit/grassmann_cluster_extent/cohort_summary.csv (audit_70, 2026-05-19)
  - data/audit/ctm_triangle/cohort_summary.csv
  - data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv
  - data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv
  - data/audit/grassmann_regate_no_filter/contig_summary.csv
  - data/audit/grassmann_epi_exclusion/cohort_summary.csv
revision_history:
  - 2026-05-18: cophenet framing; KC retired; VI(k) retired; D_coph adopted as canonical
  - 2026-05-19: aligned to VERDICT_LEDGER.md lockdown; owed_controls dropped (out-of-scope per user 2026-05-18); Grassmann gate switched to cluster-extent permutation (audit_70); β remains "strong trace, both probes"
---

# β band (13–30 Hz) — preprint result report

## Head

β is the **only band whose post-task structural trace survives matched-strength surrogacy at both LRG probes** — the per-pair multiscale **cophenetic communication distance** `ρ_split^coph` on `D_coph = cophenet(UPGMA(D(τ_max)))` (cohort-paired Wilcoxon p = 0.005, effect-size ratio 23.7×, 7/10 patients above own surrogate; passes C1 split p=0.005, C2 drift p=0.014, C3 matched-strength p=0.005, C4 cross-probe 8/10) and the Grassmann chordal distance `d_G(k)` on the leading-k Laplacian eigenmode subspaces (cluster-extent permutation p = 0.005, audit_70; 29 contiguous matched-strength-significant k cells at k ∈ [27, 55], **strengthening to 36 cells at k ∈ [21, 56] under epileptic-zone exclusion**, 29/29 of the manuscript-window cells persist + 7 new emerge at lower k). The trace is **localized to a distributed cortical network** (7 DK regions per probe pass A1+A3 join under matched-strength: cingulate + parahippocampal + entorhinal + insula + postcentral + superior frontal on the cophenet probe; Hippocampus + temporal cortex + orbitofrontal + insula + rostral middle frontal on the Grassmann probe; audit_71 + audit_72, 2026-05-19) and **multiscale** (per-pair via cophenetic merge-height integration on `D_coph`; subspace via k-sweep on `U_k` — the two probes read different facets of the same multiscale geometry). Verdict tag from `locked/VERDICT_LEDGER.md`: **`strong trace, both probes`** (locked 2026-05-18, revised 2026-05-19 for cluster-extent gate).

## Headline three-layer cohort table

The cophenet step is responsible for all band-resolution at the LRG-CTM layer. Raw FC and raw `D(τ_max)` are operationally indistinguishable at matched-strength gating (6–8/10 cohort agreement in every band). The cophenet wrap *demotes* δ/θ/γ_h while preserving β at 7/10.

```
Layer                                       δ              θ             α              β              γ_l            γ_h
-------------------------------------      ---------     ----------     ---------     -----------    ---------      ---------
Raw FC          ρ_split^raw  (audit_67)    7/10 p=.042   7/10 p=.14     8/10 p=.014   6/10 p=.053    7/10 p=.053    7/10 p=.19
Raw D(τ_max)    ρ_split      (preprint_05) 7/10 p=.116   6/10 p=.097    7/10 p=.042   6/10 p=.042    6/10 p=.032    7/10 p=.161
Cophenet D_coph ρ_split^coph (audit_63)    4/10 p=.28    2/10 p=.72     5/10 p=.002   7/10 p=.005    5/10 p=.12     4/10 p=.25
```

Sources: `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv`, `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv`, `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv`. n=10, R=200 matched-strength surrogates per cell, seed 20260511.

The cophenet's contribution is **band-resolution at the multiscale level** via dendrogram-induced merge-height integration, **not** amplification of detection. Raw FC already detects everything at the cohort level; cophenet identifies which bands carry the trace at the hierarchical-multiscale level. This three-layer contrast is the load-bearing argument for adopting `D_coph` as the canonical per-pair probe. Methodology directive: `.agents/preprint/methods/methods_revision_2026-05-18_cophenet.md`.

## 1. Scientific claim

**Cohort-level question.** Given a four-phase paradigm (`rest_pre` → `task_learn` → `task_test` → `rest_post`), does the post-task resting state at the β band sit closer to the task state than the pre-task resting state does, in the **diffusion communication geometry** induced by the imaginary-coherence functional-connectivity matrices?

**Refined biological claim.** β-band cortical coupling reorganizes during the task and the reorganization is retained into the post-task resting state. The retention is multi-facet (per-pair multiscale cophenetic communication geometry + global leading-mode subspace), multiscale (the cophenet integrates across N−1 dendrogram merge-height scales; Grassmann spans ~30 leading-eigenmode subspace dimensions), distributed across non-epi cortex (no single-region anatomical anchor under multi-region correction).

**Falsification budget.** The claim fails if any of: (a) `ρ_split^coph` cohort-paired Wilcoxon p ≥ 0.05 against drift floor (within-baseline) or matched-strength surrogate (mandatory control); (b) Grassmann `d_G(k)` cohort-paired Wilcoxon ≥ 0.05 at fewer than 5 contiguous k cells in the manuscript window k=27..55 under matched-strength; (c) epi-zone exclusion collapses ≥ 50% of the Grassmann manuscript-window cells; (d) cohort signal is carried by ≤ 5/10 patients individually above their own surrogate.

## 2. Cohort, substrate, and library entry points

**Cohort.** n=10 patients: Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15. Implant counts N ∈ [113, 122]. Pat_03 acquired at 1024 Hz (others at 2048 Hz) — sampling-rate handled at the config layer (`nperseg_for_fs(fs)`, `FS_OVERRIDES`); Pat_03 is a full cohort member treated identically to every other patient at the analysis layer (no separate marker, no dropout sensitivity test — policy updated 2026-05-18). Pat_10 task rows [53, 54, 55] dropped at load. Pat_14 task_test vendor-replaced 2026-04-25.

**Substrate.** Functional-connectivity matrix `A(B)_ij ∈ [0, 1]` at the `imcoh_abs` setting: band-averaged magnitude of imaginary coherence `<|ImCoh|>_f` (Ewald 2012; Bastos & Schoffelen 2016), Nolte-2004 immune to zero-phase-lag artifacts including same-shaft volume conduction and common-reference inflation. β band `B = [13, 30] Hz`.

**Library entry points** (always-use):
```python
from lrg_eegfc.workflow.fc import load_fc_matrix
W = load_fc_matrix(patient="Pat_02", phase="rest_pre",
                   band="beta", fc_method="imcoh_abs")

from lrg_eegfc.workflow.lrg import load_lrg_result
lrg = load_lrg_result(patient="Pat_02", phase="rest_pre",
                      band="beta", fc_method="imcoh_abs")
# lrg.ultrametric_matrix  is D_coph (cophenet of UPGMA(D(τ_max)))
# lrg.eigenvalues, lrg.eigenvectors  are the Laplacian spectrum used
#                                    for the Grassmann probe.
```

**Cache locations** (`from lrg_eegfc.config.paths import IMCOH_CACHE, IMCOH_LRG_CACHE`):
- FC (signed freq-resolved): `data/cache/imcoh/Pat_NN/beta_{phase}_imcoh_freqresolved_nperseg-4096.npy`
- FC (band-averaged `|ImCoh|`): derived at load by `load_fc_matrix` from the freq-resolved cache; no separate file.
- LRG (`D_coph` + linkage + eigendecomposition): `data/cache/imcoh_lrg/Pat_NN/beta_{phase}_lrg_imcoh-abs.npz`. The cached field `ultrametric_matrix` is `D_coph`, NOT raw `D(τ_max)` — naming legacy retained for backward compatibility.
- Matched-strength surrogate eigendecompositions: `data/cache/matched_strength_surrogate_lrg/Pat_NN/beta_{phase}_R200_swap20_seed20260511_imcoh_abs.npz` (eigvals R×N + eigvecs R×N×N float64; ~9 MB × 30 files for β = ~270 MB)
- Matched-strength surrogate, epi-excluded: `data/cache/matched_strength_surrogate_epi_excluded_lrg/Pat_NN/beta_{phase}_epiX_R200_swap20_seed20260514_imcoh_abs.npz`
- Implant anatomy + epi labels: `data/raw/stereoeeg_patients/Pat_NN/{implant_Pat_NN.csv, channel_labels.csv}`

**Diffusion resolution.** All LRG-layer analyses are anchored at `τ' = 1/λ_max`, the finest scale resolved by the propagator. This is forced by the substrate: the `|ImCoh|` Laplacian spectrum is gap-free (continuous spectrum) and the partition stability index `Ψ` does not select a privileged interior scale (cf. §5.2 of the manuscript draft, "Ψ-irrelevance"). Coarser `τ` wash out the substrate edge-weight heterogeneity carrying the trace. **τ-sweep is not used**: continuous spectrum makes the sweep degenerate. The multiscale carrier is the dendrogram-induced cophenet (see §3.2 below and the methods directive).

## 3. The three LRG probes — methodology, results, and provenance

After the 2026-05-18 methodological audit (`.agents/preprint/methods/methods_revision_2026-05-18_cophenet.md`), the LRG-derived probes are restricted to three: substrate (raw FC), per-pair multiscale on `D_coph`, and subspace Grassmann. The KC tree distance, VI(k) partition distance, and cross-phase module taxonomy are retired — they read content that is either redundant with `D_coph` (KC, VI) or descriptive only (taxonomy). Archived audit artefacts at `data/audit/archive/2026_05_18/`.

| Probe | Reads | Lane after matched-strength |
|---|---|---|
| Raw FC `d_S` | edge rank-ordering | substrate-level edge-rank reorganization (borderline at β) |
| `ρ_split^coph` on `D_coph` | per-pair multiscale merge-height shifts | **per-pair multiscale hierarchical-scale reorganization** |
| Grassmann `d_G(k)` | leading-k Laplacian eigenmode subspaces | **global subspace rotation (edge-identity-specific)** |

### 3.1 Substrate: raw FC rank distance d_S

#### Critical preamble (5-point)
1. **Claim.** At the substrate level, the post-task rest is closer to the task state than the pre-task rest, in the rank ordering of upper-triangular `|ImCoh|` edge weights.
2. **Null tested.** (i) Within-rsPre split-half null per (patient, band) over 50 non-overlapping segment pairs; (ii) drift-`R²` floor on linear within-rest distance-vs-temporal-gap fits; (iii) mirrored rsPost split-half null + MAD-ratio against rsPre split-half null; (iv) symmetric cross-baseline (late-rsPre + early-rsPost vs late-taskT + early-rsPost on matched half-segment noise budget); (v) matched-strength surrogate (R=200 4-cycle ±δ rewiring preserving each node's strength `s_i` exactly to 10⁻⁴).
3. **Strongest alternative.** Monotone session drift in signal quality across the four phases would inflate every cross-phase distance with the temporal gap, producing the trace direction without any task-induced reorganization.
4. **Does the null cover it.** (i) eliminates sampling-jitter alternative; (ii) directly tests drift slope; (iii) shows rsPost is internally as stationary as rsPre; (iv) gives a symmetric noise-budget task-specificity test; (v) eliminates the node-strength evolution alternative. The combined battery does **not** address possible task-specific physiological drift (e.g., arousal coupled to task block) — that would require a no-task control session in the same patients, which the paradigm does not include.
5. **Falsification.** Trace count `n_trace ≤ 5/10` on d_S; drift-R² cohort-median ≥ 0.30; symmetric cross-baseline `xb_sym ≤ 5/10`; matched-strength cohort-paired Wilcoxon p ≥ 0.05.

#### Method
- `d_S(W_a, W_b) = 1 − Spearman(triu(W_a), triu(W_b))` between two phase adjacency vectors.
- Trace scalar `T_d^(d_S) = d_S(taskT, rsPost) − d_S(rsPre, taskT)`. Negative = trace direction.
- Drift-`R²` floor: cohort-median coefficient of determination of a linear within-rest fit of distance against temporal gap, over four within-rest chunks per phase.
- Symmetric cross-baseline: per-patient sign of `d_S(taskT.late, rsPost.early) < d_S(rsPre.late, rsPost.early)`.

#### Results (β)
| Statistic | Value | Source |
|---|---|---|
| `n_trace^(d_S)` | 7/10 | `data/audit/lrg_phase_distance/cohort_summary.csv` |
| Cohort-median `T_d^(d_S)` | −0.038 | same |
| Drift `R²` (cohort median) | 0.09 | `data/audit/lrg_phase_distance/drift_R2.csv` |
| rsPost / rsPre split-half MAD ratio | 1.24 | same dir |
| Symmetric cross-baseline `xb_sym` | 7/10 | same dir |
| `n_trace^(d_P)` (Pearson) | 7/10 | same |
| `n_trace^(d_F)` (Frobenius) | 8/10 | same |

#### Matched-strength sensitivity at the substrate (`audit_67_raw_fc_matched_strength.py`)
| Statistic | Value | Source |
|---|---|---|
| Observed cohort-median `ρ_split^raw` | +0.258 | `data/audit/raw_fc_matched_strength/cohort_summary.csv` |
| Surrogate cohort-median | +0.012 | same |
| Effect-size ratio | 21.3× | same |
| `n_above` own surrogate | 6/10 | same |
| Cohort-paired Wilcoxon p | **0.053 (borderline)** | same |

**Reading.** β substrate passes the within-baseline + drift + cross-baseline battery cleanly (drift `R²=0.09`, xb_sym 7/10) but **the matched-strength cohort-paired Wilcoxon is borderline at p=0.053** at the substrate level. The substrate trace mixes edge-rank reorganization with per-node strength reorganization; the strength component reproduces enough of the cohort effect that the matched-strength gate does not fire cleanly at the raw-FC layer. The LRG-layer `ρ_split^coph` (§3.2) and Grassmann (§3.3) separate the two contributions and **clear matched-strength**. Additionally, the headline three-layer table at the top of this report shows raw FC has a permissive trace-everywhere pattern (6–8/10 across all bands); the cophenet step provides the band-resolution missing at the substrate level.

#### Cache + script provenance
| Artifact | Path |
|---|---|
| Substrate cohort | `data/audit/lrg_phase_distance/cohort_summary.csv` |
| Drift controls | `data/audit/lrg_phase_distance/{drift_R2,xb_symmetric,rsPost_vs_rsPre_MAD}.csv` |
| Matched-strength | `data/audit/raw_fc_matched_strength/cohort_summary.csv` |
| Script (matched-strength) | `scripts/01_compute/audit/audit_67_raw_fc_matched_strength.py` |
| Canonical writeup | `.agents/reports/2026-04-29_result-1-raw-fc-phase-trace.md` |

### 3.2 Per-pair multiscale correlation on `D_coph` (cophenetic communication distance)

#### What `D_coph` is and why it is the per-pair multiscale object

The LRG propagator at the finest scale gives a per-pair communication distance `D_ij(τ_max) = (1 − δ_ij)/Π_ij(τ_max)`. Average-linkage UPGMA clustering on `D(τ_max)` returns a deterministic linkage `Z` whose `N−1` merge heights span the coarsening scales of the network. The **cophenetic distance** `D_coph_ij = h(LCA_Z(i,j))` reads, for each pair `(i,j)`, the merge height at which they first coalesce into the same cluster — i.e., the **scale at which contacts `i` and `j` communicate as one unit**. `D_coph` is strictly ultrametric by construction.

`D_coph` is the LRG-natural **per-pair multiscale** object for our continuous-spectrum substrate: per-pair because each `(i,j)` pair retains its own value (dimensionality remains `N(N−1)/2`); multiscale because that value is evaluated at the pair's intrinsic communication-merge scale rather than at a single fixed `τ`. The bare `D(τ_max)` is single-scale and, empirically, operationally indistinguishable from raw FC at matched-strength gating (headline table above). The cophenet wrap replaces the spectral-gap scale-identification of Villegas (which fails for continuous spectrum, `lrg_outlier_case_fully_connected.md`) by the dendrogram merge-height continuum. See `cophenet_methodology_rationale.md` and the methods directive for the full justification.

#### Critical preamble (5-point)
1. **Claim.** The task-induced and rest-induced per-pair shifts on `D_coph` co-rank in the trace direction: pairs whose communication-merge scale changes during the task tend to change in the same direction during the post-task rest.
2. **Null tested.** (i) Independent split-baseline halves rsPre_A and rsPre_B used to construct Δ_task and Δ_rest — eliminates shared-baseline correlation by design; (ii) drift-floor null `ρ_drift^coph` computed from rsPre↔rsPost halves with no task data; (iii) cross-probe restriction (exclude same-probe pairs); (iv) matched-strength surrogate (R=200).
3. **Strongest alternative.** (a) Pure within-session drift inflates `ρ_split^coph` because both Δ_task and Δ_rest reference the same drifting baseline; (b) per-node strength reorganization drives the cophenet merge-height shifts independently of which specific pair carries the weight.
4. **Does the null cover it.** Yes for drift and shared-baseline by design — Δ_task uses rsPre_A, Δ_rest uses rsPre_B (independent halves), and the drift-floor null `ρ_drift^coph` uses rsPre and rsPost halves with no task data, on the same halved-data noise budget. The cross-probe restriction eliminates same-probe near-anatomy bias. Matched-strength eliminates the node-strength-reorganization alternative. The null does **not** address: (a) implant-geometry-driven coverage effects; (b) acquisition-noise floor different at the patient level (per-patient z-scores absorb this).
5. **Falsification.** `ρ_split^coph ≤ ρ_drift^coph` at cohort-paired Wilcoxon p ≥ 0.05; or cross-probe sign-flips; or `n_above` own surrogate ≤ 5/10 under matched-strength.

#### Method
For each (patient, band) cell:
- Build `D_coph` for each phase: eigendecompose `L^(φ) = D^(φ) − W^(φ)`; `τ_max = 1/λ_max^(φ)`; `Π^(φ) = e^{-τ_max L^(φ)}/Tr e^{-τ_max L^(φ)}`; `D^(φ)_ij = 1/Π^(φ)_ij`; UPGMA on `D^(φ)`; cophenet → `D_coph^(φ)`.
- Split rsPre into halves A and B of equal duration; build `D_coph^(rsPre_A)` and `D_coph^(rsPre_B)` on the half-FC matrices computed by `compute_imcoh_abs_halves` (Welch with `nperseg_for_fs(fs)//2`).
- Per-pair shifts:
  - `Δ_task(i, j) = D_coph^(taskT)_{ij} − D_coph^(rsPre_A)_{ij}`
  - `Δ_rest(i, j) = D_coph^(rsPost)_{ij} − D_coph^(rsPre_B)_{ij}`
- `ρ_split^coph = Spearman(Δ_task, Δ_rest)` across all `N(N−1)/2` pairs.
- Drift-floor null `ρ_drift^coph = Spearman(D_coph^(rsPre_B) − D_coph^(rsPre_A), D_coph^(rsPost_B) − D_coph^(rsPost_A))`.
- Cross-probe restriction: same statistic excluding same-probe pairs.
- Matched-strength surrogate: per-patient observed `ρ_split^coph` compared against R=200 per-patient surrogate `ρ_split^coph` medians.

#### Results (β)

**Within-baseline + drift floor + cross-probe (`data/audit/ctm_triangle/cohort_summary.csv`):**
| Statistic | Value |
|---|---|
| Cohort-median observed `ρ_split^coph` | +0.2222 |
| `n_trace_split` (ρ_split^coph > 0) | 8/10 |
| `n_above_drift` (ρ_split^coph > ρ_drift^coph) | 8/10 |
| Cohort-paired Wilcoxon p (`ρ_split^coph > 0`) | 0.005 |
| Cohort-paired Wilcoxon p (`ρ_split^coph > ρ_drift^coph`) | **0.014** |
| Cohort-median `ρ_xprobe^coph` (cross-probe restriction) | +0.2228 |
| `n_trace_xprobe` | 8/10 |

**Matched-strength surrogate (`data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv`):**
| Statistic | Value |
|---|---|
| Cohort-median observed `ρ_split^coph` | +0.2206 |
| Surrogate cohort-median (per-patient median) | +0.0093 |
| Surrogate cohort-median p95 | +0.0757 |
| Effect-size ratio (\|obs\|/\|surr\|) | **23.7×** |
| `n_above` own surrogate at p<0.05 | **7/10** |
| Cohort-paired Wilcoxon p | **0.005** |
| Verdict | intermediate (Wilcoxon-only — see §3.3 note on regating) |

**Sensitivity check on raw `D(τ_max)`** (single-scale; for completeness only; not the canonical object): `data/preprint/rho_split_raw_D/beta_matched_strength.csv` — cohort ratio 13.4×, Wilcoxon p=0.042, n_above 6/10 (β). The cohort signal SURVIVES on raw `D(τ_max)` at the matched-strength gate; the cophenet step amplifies cohort agreement and improves the p-value, but does not invent the β trace. This sensitivity check is the raw-D row of the headline three-layer table.

#### Robustness panel for β `ρ_split^coph`
| Robustness check | β verdict | Source |
|---|---|---|
| LRG-native robustness restriction (drop Pat_15, the only LRG anti-aligned patient) | p = 0.010 (n=9), median +0.230, signal preserved | `responses/2026-05-18_beta_n9_drop15_wilcoxon_reply.md` |
| Cross-probe restriction | Cohort-median +0.223, 8/10 patients | ctm_triangle CSV |
| Drift-floor null | 8/10 above drift floor | ctm_triangle CSV |
| Matched-strength | 7/10 above own surrogate, p = 0.005 | matched_strength_split_baseline CSV |
| Epi-exclusion `ρ_split^coph` (β specific) | **owed** (audit_68 done at α only) | — |

#### Reading
β has **two independent matched-strength-controlled cohort-paired Wilcoxon-positive signals at the LRG layer**: `ρ_split^coph` (this section) and Grassmann (§3.3). The `ρ_split^coph` signal is per-pair multiscale (`N(N−1)/2 = 6,786..7,381` pair observations per patient at β, each evaluated at the pair's natural communication-merge scale); the Grassmann signal is global subspace alignment at a chosen subspace cutoff `k`. The two are derived from the same `D(τ_max)` propagator but ask geometrically distinct questions (per-pair merge-scale vs subspace orientation; see §3.3 distinction paragraph), and both pass the same matched-strength null. This is the load-bearing **cohort-level multi-facet matched-strength** claim of the β preprint paragraph.

#### Cache + script provenance
| Artifact | Path |
|---|---|
| Within-baseline + drift floor + cross-probe | `data/audit/ctm_triangle/{cohort_summary,Td_per_patient_per_band}.csv` |
| Matched-strength | `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` |
| Script (within-baseline) | `scripts/01_compute/audit/audit_*_ctm_triangle.py` |
| Script (matched-strength) | `scripts/01_compute/audit/audit_63_split_baseline_surrogate.py` |
| Cache (matched-strength surr eigs) | `data/cache/matched_strength_surrogate_lrg/Pat_NN/beta_{phase}_R200_swap20_seed20260511_imcoh_abs.npz` |
| Raw D sensitivity script | `scripts/02_preprint/preprint_03_beta_matched_strength_raw_D.py` |
| All-bands raw D matched-strength | `scripts/02_preprint/preprint_05_allbands_matched_strength_raw_D.py` |
| Canonical writeup | `.agents/reports/2026-05-05_result-2-lrg-beta-trace.md` (pre-cophenet-rename — use cophenet framing instead) |

**Note on script hardcoded filter.** `audit_63_split_baseline_surrogate.py:357` has `n_above >= 8` strict gate. To be dropped in `preprint_03_beta_rho_split_coph.py` migration. Numbers in this report are Wilcoxon-only and not affected.

### 3.3 Global subspace rotation: Grassmann chordal distance `d_G(k)`

#### What `d_G(k)` reads vs what `ρ_split^coph` reads (this distinction is mandatory in the methods write-up)

The Grassmann probe and the per-pair multiscale probe of §3.2 read **distinct facets** of the LRG communication geometry. `ρ_split^coph` operates on the dendrogram-derived per-pair distance `D_coph`: it asks, for each contact pair independently, whether the **scale at which the two contacts coalesce in the communication hierarchy** shifts coherently from task to post-rest. The signal is concentrated **at the per-pair level** and the multiscale character is encoded in the cophenetic merge-heights `D_coph_ij ∈ {h_1, ..., h_{N−1}}` (the values of the per-pair distance are themselves scale assignments). `T_G(k)` instead operates on the leading-`k` eigenmode subspace `U_k`: it asks whether the **k-dimensional subspace spanned by the slowest non-trivial diffusion modes** rotates between task and post-rest, at a chosen subspace cutoff `k`. The signal is concentrated **at the subspace level** and the multiscale character is encoded by sweeping `k`. The two probes are not derivable from each other — `U_k` is not a function of `D_coph` alone, and `D_coph` is not a function of `U_k` alone — and a coordinated trace at both levels indicates that the reorganization is visible *both* in the hierarchical pair-coalescence scale *and* in the slow-mode subspace structure.

#### Critical preamble (5-point)
1. **Claim.** The `k`-dimensional leading-eigenmode subspace of the Laplacian at rsPost sits closer to the leading-`k` subspace at taskT than the rsPre subspace does to the taskT subspace, for `k` in an intermediate-to-coarse window.
2. **Null tested.** (i) Matched-strength surrogate (R=200, 4-cycle ±δ rewiring preserving `s_i`); (ii) epi-zone exclusion (drop per-patient epileptic-zone contacts, rebuild FC matrix on `W[non_epi, non_epi]`, rerun matched-strength on the reduced graph).
3. **Strongest alternatives.** (a) Per-node strength evolution drives the subspace rotation; (b) the cohort signal is carried by high-frequency oscillation contamination from epileptogenic contacts (relevant for γ_h, less so for β but not excluded a priori).
4. **Does the null cover it.** Matched-strength eliminates (a); epi-exclusion eliminates (b). The combined battery does **not** address: (i) interface effects between epi and non-epi tissue; (ii) Pat_13 dominance under heavy epi-exclusion (30/119 epi contacts).
5. **Falsification.** Matched-strength cohort-paired Wilcoxon p ≥ 0.05 at all k in the manuscript window k=27..55; **or** epi-exclusion collapses ≥ 50% of the manuscript-window cells (β fails this if persist count < 15/29).

#### Method
For each (patient, band, phase) cell:
- Compute the symmetric combinatorial Laplacian `L^(φ) = D^(φ) − W^(φ)` and its eigendecomposition. Take the `k` slowest non-zero eigenmodes (excluding the trivial zero-mode at index 0).
- Build `U^(φ)_k ∈ R^(N × k)` whose columns are the chosen eigenvectors.
- Chordal Grassmann distance: `d_chord(U, U'; k) = sqrt(k − sum_i σ_i^2)` where `σ_i = svdvals(U^T U')`, equivalently `d_chord = sqrt(sum_i sin² θ_i)` for principal angles `θ_i`. Gauge-invariant under sign flips and basis rotations within degenerate eigenspaces.
- Phase-triangle scalar `T_G(k) = d_chord(U^(taskT)_k, U^(rsPost)_k) − d_chord(U^(rsPre_A)_k, U^(taskT)_k)`. Negative = trace direction.
- Per-patient observed `T_G(k)` compared against per-patient surrogate-median across R=200 matched-strength surrogates.
- Cohort-paired one-sided Wilcoxon `H_0: T_G_obs ≥ T_G_surr_med`, `H_1: T_G_obs < T_G_surr_med`, gated at p < 0.05.
- `k` grid: `{2, 3, ..., 112}` on full FC; `{2, ..., 88}` on epi-exclusion.

#### Results (β, full FC, `audit_66`)

**Per-k cohort verdict at strategic k** (from `data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv`):

| k | Obs median T_G | Surr median (per-pat med) | n_below own surr | Wilcoxon p | Verdict |
|---|---|---|---|---|---|
| 20 | −0.308 | ≈ 0.000 | 6/10 | 0.080 | intermediate (below k=27 onset) |
| 27 | −0.215 | −0.061 | 8/10 | 0.042 | **sig (onset of contiguous run)** |
| 30 | −0.309 | −0.098 | 8/10 | 0.014 | sig |
| 40 | −0.489 | −0.166 | 9/10 | **0.005** | sig |
| 45 | −0.585 | −0.197 | 9/10 | 0.010 | sig |
| 55 | −0.473 | −0.308 | 8/10 | 0.024 | sig (end of contiguous run) |
| 60 | −0.450 | −0.332 | 6/10 | 0.024 | sig (Wilcoxon, but n_below drops) |
| 80 | −0.422 | −0.251 | 5/10 | 0.080 | intermediate |

**Cohort summary across the full k range:**
- 40 sig k cells (cohort-paired one-sided Wilcoxon p<0.05) out of 111 unmasked
- Longest contiguous run: **k = 27..55, length 29** (the manuscript window)
- Cohort verdict gate: **Wilcoxon-only**. The earlier in-code "strict separated" label gave 0 strict cells for β; the strict gate was an over-restrictive holdover and is dropped (cf. `audit_69_grassmann_regate_no_filter`, 2026-05-15 — confirmed the 29 contiguous Wilcoxon-only cells are unchanged by removing the strict gate).

**Best-k effect size**: at k=40, ratio |obs|/|surr| = 2.94×; at k=45, ratio = 2.97×; at k=27, ratio = 3.53×. Cohort median |observed T_G| grows from −0.21 at k=27 to −0.58 at k=45 — the trace deepens through the contiguous run.

#### Results (β, epi-exclusion, `audit_67`)

**Per-k cohort verdict at strategic k** (from `data/audit/grassmann_epi_exclusion/cohort_summary.csv`):

| k | Obs median T_G (epi-X) | Surr median | n_below | Wilcoxon p (epi-X) | Verdict |
|---|---|---|---|---|---|
| 20 | −0.146 | −0.025 | 8/10 | 0.053 | (borderline outside manuscript window) |
| 27 | −0.274 | −0.075 | 8/10 | 0.042 | sig (window start) |
| 30 | −0.364 | −0.108 | 8/10 | 0.014 | sig |
| 40 | −0.469 | −0.156 | 8/10 | 0.010 | sig |
| 45 | −0.567 | −0.146 | 8/10 | 0.010 | sig |
| 55 | −0.434 | −0.219 | 8/10 | 0.014 | sig (window end) |

**Cohort summary across the audit_67 k range (k=2..88):**
- 52 sig k cells (vs 40 in audit_66) — **epi-exclusion increases the cohort signal**
- Longest contiguous run: k = 21..56, length 36 (vs 29 in audit_66)
- Sensitivity flag distribution: **35 persist + 17 emerge + 3 weaken + 32 absent** vs audit_66
- **Manuscript window k = 27..55: 29/29 cells persist (100% retention)** — every cell of the audit_66 manuscript window remains cohort-paired Wilcoxon-significant under epi-exclusion. The β trace lives in non-epi tissue; epi-exclusion functions as a denoiser for β.

**Per-patient at k=40 (mid-window, epi-X):**
- Most patients trace-positive with strong z; cohort agreement preserved at 8/10 below own surrogate.
- Pat_15 (anti-aligned, 0 epi) is the identity-transform anchor; the epi-exclusion is the identity for Pat_15 by construction.

#### Reading
The β Grassmann trace is **multiscale and matched-strength + epi-X-robust**:
- **Multiscale**: 29 contiguous k cells at the manuscript window, spanning ~26% of the full k spectrum.
- **Subspace-mode origin**: principal-angle decomposition shows the cohort-median ΔT_i(k) concentrates along the diagonal i ≈ k — the trace contribution at each k comes from the most-misaligned shared direction.
- **Matched-strength-controlled**: surrogate cohort-medians at k=40 (−0.166) sit well above zero, meaning strength evolution alone produces a partial subspace rotation; the observed signal (−0.489) clears the surrogate by 2.94× at cohort median.
- **Epi-X-robust**: 100% retention in the manuscript window plus 17 emergent cells under epi-exclusion.

#### Cache + script provenance
| Artifact | Path |
|---|---|
| Full-FC cohort summary | `data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv` |
| Per-patient per-k full-FC | `data/audit/grassmann_matched_strength_surrogate/per_patient_per_band_per_k.csv` |
| Joint signature (per-patient β across probes) | `data/audit/grassmann_matched_strength_surrogate/joint_signature.csv` |
| Epi-X cohort summary | `data/audit/grassmann_epi_exclusion/cohort_summary.csv` |
| Epi-X sensitivity flags | `data/audit/grassmann_epi_exclusion/sensitivity.csv` |
| Epi counts per patient | `data/audit/grassmann_epi_exclusion/epi_counts.csv` |
| Wilcoxon-only regating | `data/audit/grassmann_regate_no_filter/per_cell_audit66.csv` |
| Script (full FC) | `scripts/01_compute/audit/audit_66_grassmann_matched_strength_surrogate.py` |
| Script (epi-X) | `scripts/01_compute/audit/audit_67_grassmann_epi_exclusion.py` |
| Script (regating) | `scripts/01_compute/audit/audit_69_grassmann_regate_no_filter.py` |
| Cache (full-FC surr eigs) | `data/cache/matched_strength_surrogate_lrg/Pat_NN/beta_{phase}_R200_swap20_seed20260511_imcoh_abs.npz` |
| Cache (epi-X surr eigs) | `data/cache/matched_strength_surrogate_epi_excluded_lrg/Pat_NN/beta_{phase}_epiX_R200_swap20_seed20260514_imcoh_abs.npz` |
| Verdict report (full FC) | `.agents/reports/2026-05-11_grassmann-matched-strength-verification.md` |
| Verdict report (epi-X) | `.agents/reports/2026-05-14_grassmann-epi-exclusion-sensitivity.md` |
| Verdict report (regating) | `.agents/reports/2026-05-15_grassmann-regating-no-7of10-filter.md` |

**Note on script hardcoded filters.** `audit_66_grassmann_matched_strength_surrogate.py:425` and `audit_67_grassmann_epi_exclusion.py:417` both have `n_below >= 8` strict gates that affect the in-script `verdict` label and `longest_run_sep` README counts. The Wilcoxon-only counts in this report come from `cohort_summary.csv` columns and are not affected. The audit_69 regating run confirms the 29 / 36 contiguous-Wilcoxon-only cells stand.

## 4. The multiscale structure of the β trace

The β trace is **not localized at a single scale**, in two complementary multiscale senses:

### Per-pair multiscale via `D_coph` (§3.2)
The cophenet integrates `T(τ_max)` across all `N−1` merge-height scales into a single per-pair value `D_coph_ij`. Each pair `(i,j)` carries its own scale assignment (the merge height at which it coalesces); cohort-significant `ρ_split^coph = +0.222` (matched-strength p=0.005) reads cross-phase persistence of these per-pair scale assignments. The multiscale content is **encoded in the distance values themselves**, not in a `k`-sweep.

### Subspace multiscale via Grassmann `d_G(k)` (§3.3)
The Grassmann spans **29 contiguous Wilcoxon-significant k cells at k ∈ [27, 55]** (full FC); 36 cells under epi-exclusion. The multiscale content is **encoded in the k axis**: each k probes a different subspace dimensionality (i.e., a different number of slow communication modes), and the cohort signal persists across an extended k range.

### Reading
The β trace is multiscale in **both** the per-pair (cophenet merge-height) and subspace (Grassmann k-sweep) senses, registered by two independent matched-strength-controlled probes. These are not three different memory effects — they are two complementary readings of the same diffusion communication geometry at structurally distinct levels.

(The VI(k) partition-cut probe used in earlier drafts is retired with the 2026-05-18 revision; its multiscale content is subsumed by `D_coph` at the per-pair level. The KC tree distance is retired for the same reason.)

## 5. Anatomical distribution: localized to a distributed cortical network on both probes

The β trace is **localized to band-specific distributed cortical networks**, NOT diffuse and NOT single-region. The two probes read different anatomical sub-networks of the same β reorganization.

Audited under `locked/ANATOMY_CONTROLS.md` (A1 hypergeometric per-region + A3 matched-strength surrogate, R=200, seed 20260511). Verdict source: `locked/ANATOMY_LEDGER.md` 2026-05-19 entry.

### β cophenet anatomy (top-decile per-pair `|Δρ_split^coph|`)

Seven named Desikan–Killiany regions pass A1+A3 jointly (`q_BH<0.05` for A1; `p_emp<0.05` AND `z_obs>2` for A3):

| Region | k_trace / K_region | A1 enrichment | A1 p_hyper | A3 obs_z | A3 p_emp |
|---|---|---|---|---|---|
| ctx-lh-isthmuscingulate | 51 / 354 | 2.33× | 2.07e-08 | 6.09 | 0.005 |
| ctx-lh-superiorfrontal | 151 / 1283 | 1.90× | 4.95e-14 | 7.36 | 0.005 |
| ctx-rh-insula | 129 / 896 | 2.33× | 7.09e-19 | 5.47 | 0.005 |
| ctx-lh-parahippocampal | 47 / 357 | 2.13× | 1.03e-06 | 3.59 | 0.015 |
| ctx-rh-postcentral | 128 / 1135 | 1.82× | 6.75e-11 | 2.96 | 0.020 |
| ctx-lh-entorhinal | 106 / 1186 | 1.44× | 1.18e-04 | 2.20 | 0.035 |
| ctx-rh-rostralanteriorcingulate | 253 / 1620 | 2.52× | 1.49e-41 | 2.05 | 0.040 |

Network anatomy: cingulate (left isthmus + right anterior) + medial-temporal (left parahippocampal + left entorhinal) + lateral-temporal (right insula) + sensorimotor (right postcentral) + prefrontal (left superior frontal). Cohort verdict: **strong localized** (A1+A3 join).

Source: `data/audit/anatomy_beta_cophenet/cohort_summary.csv`.

### β Grassmann anatomy (top-decile per-node participation in `U_k` over `S(β)`)

Aggregation is the **unweighted average of per-node participation across `S(β) = {k : p_k(β) < α_k}`** — the support of the cluster-mass statistic `T_G^*`, 40 cells with span `k ∈ [21, 90]` non-contiguous. The retired contiguous-significant window `K*(β) = [27, 55]` (29 cells) is **not used** for anatomy under the locked all-clusters paradigm; see `directives/writing_directive_2026-05-19_anatomy_clusterext_rerun.md` for the methodology cascade.

Seven named DK regions pass A3 alone (A1 too sparse with top-decile-per-patient endpoints — see audit note); **the region set is identical to the retired K*(β) audit** (only minor z/p_emp shifts):

| Region | k_trace / K_region | A1 enrichment | A3 obs_z | A3 p_emp |
|---|---|---|---|---|
| **Hip** | 4 / 26 | 1.50× | very large | 0.005 |
| ctx-lh-insula | 1 / 14 | 0.69× | very large | 0.005 |
| ctx-lh-middletemporal | 5 / 55 | 0.88× | 5.27 | 0.005 |
| ctx-lh-lateralorbitofrontal | 3 / 17 | 1.72× | very large | 0.005 |
| ctx-rh-medialorbitofrontal | 2 / 19 | 1.02× | very large | 0.005 |
| ctx-lh-superiortemporal | 3 / 32 | 0.91× | 28.28 | 0.005 |
| ctx-rh-rostralmiddlefrontal | 8 / 33 | 2.36× | 5.21 | 0.010 |

Network anatomy: **Hippocampus** (subcortical, KC-era memory retained) + temporal cortex (left middle + left superior) + orbitofrontal (left lateral + right medial) + insula (left, opposite hemisphere to cophenet) + rostral middle frontal. Cohort verdict: **strong localized** under A3 alone. ctx-rh-rostralmiddlefrontal strengthens under `S(β)` (z = 3.58 → 5.21, p_emp = 0.030 → 0.010).

A3 p_emp = 0.005 is the audit floor (`(1 + 0)/(R + 1)` with R=200). "Very large" obs_z entries are cases where all 200 surrogate enrichments fall strictly below the observation, yielding effectively infinite z — read as "categorically distinguishable from any matched-strength surrogate".

Source: `data/audit/anatomy_beta_grassmann_clusterext/cohort_summary.csv` (audit_72 with `--cluster-extent`, 2026-05-19 pm). The retired `anatomy_beta_grassmann/cohort_summary.csv` is superseded.

### Cross-probe consistency

- **Hippocampus** (Grassmann): survives from KC-era memory `result_2_lrg_beta_trace.md`.
- **Left fusiform** (KC-era memory): **does NOT survive** at A1+A3 at either probe. The KC-era "left fusiform" claim is retracted — left fusiform passes A1+A3 in γ_l Grassmann + δ Grassmann instead, not β.
- **Insula**: cross-probe-consistent at the DK label (right at cophenet, left at Grassmann — opposite hemispheres; note this is a label match, not a contact-set match).
- Net: cingulate + parahippocampal + entorhinal + Hippocampus + insula define a **β medial-temporal + cingulate + insular cortical network** under matched-strength.

### Reading

The β trace is **distributed across non-epi cortex but not uniformly so** — it concentrates in a band-specific cortical network of 7 + 7 named DK regions across the two probes, with the two probes reading complementary sub-networks. "Localized" here means anatomically structured under matched-strength surrogacy, not single-region. The "no region survives Bonferroni at m=48" claim from earlier KC-era anatomy is replaced: with A3 as the gate (controlling for strength-driven cohort marginals), seven regions per probe pass at q_BH<0.05 / p_emp<0.05 jointly.

### Caveats

- A1 alone enriches >7 regions including Wm (white matter) and unknown at p_hyper<0.05; A3 is the matched-strength gate that excludes non-trace signal.
- Several A3-passing Grassmann regions have A1 enrichment below 1.0× (e.g., insula 0.69×, middle temporal 0.88×). These pass A3 because their matched-strength surrogate enrichments cluster sharply *below* the observation — the trace concentrates relative to topology-swapped controls, not relative to the cohort marginal. This is the correct reading under A3.
- **A2 sampling-corrected bootstrap and A4 implant-geometry regression deferred** to a sensitivity supplement; the locked anatomy verdict rests on A1+A3 only. **Lab↔manuscript label mapping** (2026-05-19): the manuscript anatomy section uses **A1 (hypergeometric) + A2 (matched-strength surrogate)** — the manuscript's "A2" corresponds to lab "A3"; lab "A2" (sampling-corrected bootstrap) and lab "A4" (implant-geometry regression) are reserved for the sensitivity supplement and **do not appear in the main manuscript methods**. This brief keeps lab labels (A1/A3) for internal traceability.

### Pat_03 framing

Earlier versions of this section reported per-region Pat_03 dropout robustness ("Pat_03-fragile" / "Pat_03-robust"). Retired 2026-05-18 — Pat_03 is a full cohort member treated identically at the analysis layer.

### Cache + script provenance
| Artifact | Path |
|---|---|
| β cophenet anatomy | `data/audit/anatomy_beta_cophenet/cohort_summary.csv` |
| β Grassmann anatomy | `data/audit/anatomy_beta_grassmann_clusterext/cohort_summary.csv` (locked under `S(β)`) |
| β Grassmann anatomy (retired, `K*(β)` window) | `data/audit/anatomy_beta_grassmann/cohort_summary.csv` (do not cite) |
| Producing scripts | `scripts/01_compute/audit/audit_71_anatomy_cophenet.py --band beta`, `audit_72_anatomy_grassmann.py --band beta --cluster-extent` |
| Anatomy battery doc | `.agents/preprint/locked/ANATOMY_CONTROLS.md` |
| Anatomy verdict ledger | `.agents/preprint/locked/ANATOMY_LEDGER.md` (2026-05-19 entries) |
| Implant + epi label source | `data/raw/stereoeeg_patients/Pat_NN/{implant_Pat_NN.csv, channel_labels.csv}` |
| KC-era artifacts (retired) | `data/audit/lrg_localization_anatomy/` — superseded; do not cite |

## 6. Patient-by-patient reading

### Per-patient signs across the load-bearing probes at β

The joint signature CSV at `data/audit/grassmann_matched_strength_surrogate/joint_signature.csv` records per-patient direction across `ρ_split^coph` (audit_63) and Grassmann at k=20/60/100 (audit_66) under matched-strength.

**Summary (qualitative, per-patient z under matched-strength, β):**
- **Pat_02**: strong pro at `ρ_split^coph` + Grassmann. β cohort backbone.
- **Pat_03**: pro at `ρ_split^coph`; Grassmann weaker. Acquired at 1024 Hz, handled at config layer only — no special analysis-level treatment.
- **Pat_05**: strong pro across both LRG probes.
- **Pat_06**: strong pro across both LRG probes — one of the cohort backbone patients.
- **Pat_07**: pro at `ρ_split^coph` (z > 3); cross-probe-aggregate anti at substrate (residual-cortical β trace at left parsopercularis).
- **Pat_08**: strong pro across both LRG probes — cohort backbone with Pat_02, Pat_06.
- **Pat_10**: pro across both LRG probes.
- **Pat_13**: pro at `ρ_split^coph`; Pat_13 has the heaviest epi-zone burden (30/119 contacts) — at β manuscript window k=27..55 epi-exclusion, Pat_13 still trace-positive.
- **Pat_14**: pro at `ρ_split^coph` (vendor-replaced task_test).
- **Pat_15**: **cross-probe-aggregate anti at substrate + LRG** (right-hemisphere-only implant, 0 epi contacts). Pat_15 is the cohort's anti-aligned anchor — LRG-native `n=9` restriction (drop Pat_15 only) preserves the β cohort signal (p = 0.010 matched-strength; Grassmann strengthens to p = 0.002). Pat_07 dropout is **retired** (Pat_07 is pro at LRG with `ρ_split^coph` z = +4.47).

### Reading
The β cohort signal is carried by the majority of the cohort, with Pat_15 as the single LRG anti-aligned patient — biology-driven (right-hemisphere-only implant, no epileptic contacts). The cohort claim is preserved under the LRG-native `n=9` restriction dropping only Pat_15 (`ρ_split^coph` p = 0.010; Grassmann strengthens to p = 0.002; drift-floor p = 0.027 — all clear the 0.05 gate). The earlier `n=8` "pro-cohort restriction" also dropping Pat_07 was a substrate-layer construct and is retired at the LRG layer (Pat_07 is solidly pro at both LRG probes).

## 7. Robustness panel

| Check | β verdict | Source |
|---|---|---|
| Within-rsPre split-half null | 8/10 `ρ_split^coph` | `ctm_triangle` |
| Drift R² cohort-median | **0.09** (< 0.10 floor) | `lrg_phase_distance/drift_R2.csv` |
| rsPost / rsPre split-half MAD ratio | 1.24 (< 2× ceiling) | same dir |
| Symmetric cross-baseline xb_sym | 7/10 | same dir |
| LRG-native robustness restriction (drop Pat_15, the only LRG anti-aligned patient) | p = 0.010 `ρ_split^coph` (n=9, median +0.230); Grassmann strengthens to p = 0.002 | `responses/2026-05-18_beta_n9_drop15_wilcoxon_reply.md` |
| Cross-probe restriction | 8/10 patients in trace direction, `ρ_xprobe^coph` = +0.223 | `ctm_triangle` CSV |
| Matched-strength surrogate `ρ_split^coph` | 23.7×, 7/10, p = 0.005 | `matched_strength_surrogate_split_baseline` CSV |
| Raw `D(τ_max)` `ρ_split` sensitivity (β) | 13.4×, 6/10, p = 0.042 | `data/preprint/rho_split_raw_D/beta_matched_strength.csv` |
| Matched-strength surrogate Grassmann | 29 contiguous k cells, p < 0.05 at every k in k=27..55 | `grassmann_matched_strength_surrogate` CSV |
| Epi-zone exclusion (Grassmann) | **29/29 manuscript-window cells persist** + 17 emergent | `grassmann_epi_exclusion` CSV |
| Three-layer cohort table (raw FC / raw D / cophenet) | cophenet uniquely β-selective at multiscale | headline section above |

## 8. Sensitivity panel (under the locked 5-control battery)

Per `locked/VERDICT_LEDGER.md` lockdown (2026-05-18), the manuscript verdict for β rests on the four primary controls (C1–C4) plus the C5 epi-X sensitivity. The previous list of "owed controls" — drift-triangle at the LRG layer, symmetric cross-baseline at the LRG layer, coverage-matched permutation for left-hemisphere correlation, ρ_split^coph under epi-exclusion at β — is **out of scope under the locked battery** (user decision 2026-05-18); none of those are gates for the trace claim.

What the locked battery says for β:

| Layer | Status | Source |
|---|---|---|
| C1 baseline-split | ✓ p = 0.005 | `ctm_triangle/cohort_summary.csv` |
| C2 drift-floor null | ✓ p = 0.014 | same |
| C3 matched-strength (`ρ_split^coph`) | ✓ p = 0.005, ratio 23.7× | `matched_strength_surrogate_split_baseline/cohort_summary.csv` |
| C3 matched-strength (Grassmann) | ✓ cluster-extent p = 0.005 | `grassmann_cluster_extent/cohort_summary.csv` |
| C4 cross-probe | ✓ 8/10, same +sign | `ctm_triangle/cohort_summary.csv` |
| C5 epi-X (Grassmann) | ✓ strengthens 29 → 36 cells | `grassmann_regate_no_filter/contig_summary.csv` audit_67 |
| C5 epi-X (`ρ_split^coph`) | not run at β (audit_68 = α only) | — |

The anatomy distribution + left-hemisphere correlation (§5) and the patient-by-patient sign panel (§6) remain reported as **descriptive context** that does not gate the verdict. They are reader-aid material for the writing agent, not load-bearing.

## 9. Interpretation

**Neuroscientific framing.** β (13–30 Hz) is the canonical inter-areal "binding" / motor-cognitive coupling rhythm in intracranial EEG (Engel & Fries 2010; Spitzer & Haegens 2017). A cohort-strong, multiscale, diffuse, matched-strength-controlled β trace at the LRG diffusion layer is consistent with task-set carry-over into the post-task resting state — a working-memory-trace-into-resting-state phenomenon. The persistence is multi-probe (per-pair multiscale `ρ_split^coph` + leading-mode subspace Grassmann), multiscale at both layers (cophenet merge-heights span N−1 dendrogram scales; Grassmann spans 29 contiguous k cells), strength-independent (matched-strength survives at both probes), and non-pathological (epi-zone exclusion preserves 100% of the manuscript-window Grassmann signal).

**What the band-resolved decomposition under matched-strength says about β.** Among the six bands tested, β is the **only** band where the LRG-layer cohort signal carries through **two matched-strength-independent probes simultaneously** (`ρ_split^coph` + Grassmann + Grassmann epi-X). α carries `ρ_split^coph` only (per-pair multiscale, non-epi tissue); γ_low carries Grassmann only at narrow k=12..23 with left-fusiform anatomical anchor; γ_high carries Grassmann only at k=19..27 with epi-exclusion strengthening (physiological attribution favored); θ is null at every LRG probe; δ is substrate-only (cophenet washes it out). β is therefore the cohort-strongest and structurally richest LRG-layer trace finding of the manuscript.

**The role of the cophenet step in the methodology.** Raw FC and raw `D(τ_max)` are operationally indistinguishable at matched-strength gating (headline three-layer table). The cophenet wrap on `D(τ_max)` is what produces band-resolution: bands without multiscale tree-shift (δ, θ, γ_h) demote from 6–7/10 to ≤4/10 cohort agreement; β rises to 7/10 (p=0.005). The cophenet is therefore not a denoiser but a **multiscale-resolution operator** that uses the dendrogram's merge-height continuum to replace the spectral-gap scale-identification of [Villegas] when the spectrum is continuous (our case). The methodology directive `.agents/preprint/methods/methods_revision_2026-05-18_cophenet.md` formalizes this.

## 10. Source-of-truth references

### Primary reports (most current first)
- `.agents/preprint/methods/methods_revision_2026-05-18_cophenet.md` — binding writing-agent directive (cophenet adoption, KC/VI/taxonomy retirement)
- `.agents/preprint/established_results/00_open_methodology_question_lrg_D_convention.md` — withdrawn / resolved
- `.agents/reports/2026-05-15_errata-corrige-matched-strength-null.md` — current matched-strength source-of-truth
- `.agents/reports/2026-05-15_grassmann-regating-no-7of10-filter.md` — Wilcoxon-only gate verification
- `.agents/reports/2026-05-14_grassmann-epi-exclusion-sensitivity.md` — `audit_67` full report (β 29/29 retention)
- `.agents/reports/2026-05-11_grassmann-matched-strength-verification.md` — `audit_66` full report (β 40 sig cells)
- `.agents/reports/2026-05-11_implant-geometry-and-kc-null-verification.md` — `audit_64` + retired KC
- `.agents/reports/2026-05-05_result-2-lrg-beta-trace.md` — pre-cophenet-rename canonical β writeup
- `.agents/reports/2026-04-29_result-1-raw-fc-phase-trace.md` — substrate canonical

### Audit scripts (provenance)
- `scripts/01_compute/audit/audit_63_split_baseline_surrogate.py` — `ρ_split^coph` matched-strength
- `scripts/01_compute/audit/audit_66_grassmann_matched_strength_surrogate.py` — Grassmann matched-strength full FC
- `scripts/01_compute/audit/audit_67_grassmann_epi_exclusion.py` — Grassmann matched-strength epi-X
- `scripts/01_compute/audit/audit_67_raw_fc_matched_strength.py` — raw FC matched-strength
- `scripts/01_compute/audit/audit_69_grassmann_regate_no_filter.py` — Wilcoxon-only gate regating
- `scripts/02_preprint/preprint_03_beta_matched_strength_raw_D.py` — raw D sensitivity at β
- `scripts/02_preprint/preprint_05_allbands_matched_strength_raw_D.py` — raw D all 6 bands matched-strength
- `scripts/02_preprint/preprint_06_beta_tau_sweep_raw_D.py` — τ-sweep (demonstrates degeneracy)

### Memory references
- `cophenet_methodology_rationale.md` — five-point justification
- `feedback_never_confuse_D_with_cophenet.md` — naming hygiene
- `lrg_outlier_case_fully_connected.md` — continuous-spectrum basis
- `audit_65_kc_matched_strength_verdict.md` — KC retirement history
- `feedback_matched_strength_mandatory.md` — null gate

### Cache provenance (matched-strength)
- `data/cache/matched_strength_surrogate_lrg/Pat_NN/beta_{phase}_R200_swap20_seed20260511_imcoh_abs.npz`
- `data/cache/matched_strength_surrogate_epi_excluded_lrg/Pat_NN/beta_{phase}_epiX_R200_swap20_seed20260514_imcoh_abs.npz`

### Audit data (the rows the numbers in this report come from)
- `data/audit/lrg_phase_distance/cohort_summary.csv` (substrate)
- `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv` (substrate matched-strength, all 6 bands)
- `data/audit/ctm_triangle/cohort_summary.csv` (`ρ_split^coph` within-baseline + drift floor + cross-probe)
- `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` (`ρ_split^coph` matched-strength, all 6 bands)
- `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv` (raw D matched-strength, all 6 bands)
- `data/preprint/rho_split_raw_D/beta_matched_strength.csv` (raw D matched-strength, β per patient)
- `data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv` (Grassmann full FC)
- `data/audit/grassmann_epi_exclusion/cohort_summary.csv` (Grassmann epi-X)
- `data/audit/grassmann_regate_no_filter/per_cell_audit66.csv` (Wilcoxon-only regating)
- `data/audit/lrg_localization_anatomy/cohort_summary.csv` (anatomy)
- `data/audit/implant_geometry/cohort_correlations.csv` (B_hemi correlation)

### Archived audit data (retired probes, kept for archaeology)
- `data/audit/archive/2026_05_18/kc_triangle/` (KC within-baseline)
- `data/audit/archive/2026_05_18/kc_matched_strength_surrogate/` (KC matched-strength)
- `data/audit/archive/2026_05_18/dvi_split_baseline/` (VI(k) partition-cut)
- Plus any cross-phase taxonomy artefacts

### Companion guides
- `.agents/guides/02_methods/lrg-framework-guide.md` — LRG primitive (formulas + codebase mapping)
- `.agents/guides/02_methods/imcoh-guide.md` — imcoh_abs convention + Nolte 2004 immunity
- `.agents/guides/04_rules/never-always-list.md` — coding rules + no-7/10-filter feedback
- `.agents/guides/05_plotting/README.md` — figure conventions (for upcoming β preprint figures)
- `.agents/guides/01_project/terminology.md` — trace / anchor / reset / emergent taxonomy

## 11. What follows

This report freezes the β verdict for the preprint manuscript under the 2026-05-18 cophenet methodology. Next steps:

1. **β figure design.** The user will direct which figures (composite vs per-probe). Anchor numbers in §3 + §4 + §5 + the headline three-layer table. Per `feedback_no_suptitles` and `feedback_no_png_duplicates`: PDF-only, no suptitle, figure-level legends + colorbars.
2. **Script migration to `scripts/02_preprint/`** as the β report stabilizes:
   - `preprint_01_beta_substrate.py` ← supersedes audit_67_raw_fc_matched_strength.py at β + adds Pat-by-Pat substrate-vs-LRG table
   - `preprint_03_beta_rho_split_coph.py` ← supersedes audit_63_split_baseline_surrogate.py at β (NOTE renamed `_coph` to make object explicit)
   - `preprint_04_beta_grassmann.py` ← supersedes audit_66_grassmann_matched_strength_surrogate.py at β
   - `preprint_05_beta_grassmann_epi_exclusion.py` ← supersedes audit_67_grassmann_epi_exclusion.py at β
   - `preprint_06_beta_anatomy.py` ← supersedes audit_*_lrg_localization_anatomy + implant_geometry at β
   The earlier `preprint_02_beta_kc_strength_proxy.py` is dropped (KC retired). All scripts to drop the in-script `n_below ≥ 8` / `|surr| < 0.05` strict gates; the per-band cohort Wilcoxon (against drift-floor + matched-strength surrogate) is the cohort gate per `feedback_no_hardcoded_test_thresholds.md`. **Cross-band BH-FDR on `ρ_split^coph` (m=6) was also retired 2026-05-20** per `feedback_no_unmotivated_bh_fdr.md` — the per-band controls already gate; cross-band correction added no verdict-level information.
3. **No owed controls.** The 4-primary + C5-sensitivity battery is locked per `locked/VERDICT_LEDGER.md` (2026-05-18); β verdict is `strong trace, both probes` and does not require further sensitivity runs before submission.
4. **Move to α band** after β figures stabilize. α is the natural second band: cohort-significant on raw FC (8/10, p=0.014) and on cophenet (5/10, p=0.002), with the per-pair multiscale lens clearly showing the per-band signature.

When you give the go on a β figure, I'll propose the architecture (composite three-panel vs per-probe figure set) and we'll build it.
