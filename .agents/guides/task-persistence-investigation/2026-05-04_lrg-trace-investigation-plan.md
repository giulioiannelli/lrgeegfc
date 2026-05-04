---
name: lrg-trace-investigation-plan
type: scope
era: COHORT_N10
status: draft
created: 2026-05-04
updated: 2026-05-04
pointers:
  - .agents/reports/2026-04-29_result-1-raw-fc-phase-trace.md
  - .agents/reports/2026-04-28_raw-fc-phase-distance-verdict.md
  - .agents/reports/2026-04-27_task-persistence-reconciliation.md
  - .agents/reports/2026-04-29_critical-post-mortem.md
  - .agents/reports/2026-04-29_measure-correctness-audit.md
  - .agents/guides/task-persistence-investigation/2026-04-29_decision-rules.md
  - .agents/guides/task-persistence-investigation/2026-04-29_cohort-coverage-matrix.md
  - .agents/guides/task-persistence-investigation/2026-04-25_cbr-investigation.md
  - .agents/guides/task-persistence-investigation/2026-04-26_multiscale-partition-coherence.md
  - data/audit/raw_fc_phase_distance/final_verdict_table_with_controls.csv
  - data/audit/raw_fc_phase_distance/Td_per_patient_per_band.csv
---

# Section 5 LRG-trace investigation plan

**Section 4 of `notes_imcoh.tex` established at the substrate level — on raw `imcoh_abs` matrices, after three drift controls — that the task leaves a band-specific, per-patient trace under the rank distance d_S, with β + low_γ passing all controls cleanly, α passing with a drift caveat, δ a soft trace, and θ + high_γ at or below noise. Section 5 lifts the same triangle question
T_LRG(p, b) = d_LRG(task_test, rest_post) − d_LRG(rest_pre, task_test)
to the LRG hierarchy at τ = 1/λ_max, on three primitives — the diffusion propagator ρ(τ), the communication ultrametric distance D(τ), and the dendrogram T(τ) — to ask whether the LRG layer agrees with the substrate, sharpens it, or surfaces structure (θ-resolution, high_γ recovery) that the substrate could not see. This document is the methods-review + empirical-pilot plan that selects the load-bearing global probe and the load-bearing localization probe Section 5 will cite, with bidirectional cross-validation between the two layers and explicit coverage of the negative-verdict probes (FTD, RSA, E1 Grassmann, Ψ_L) so the recommendation is grounded in the full phase-comparison repertoire that has been tried.**

## Reading guide

The document has five parts:

- **Part 1 — Methods review.** Comparative table over every repo-resident probe (CBR family, MSPC, CNP, CTM, VI(k), FTD, RSA, E1 Grassmann, Ψ_L, EH-discriminability, four-phase-τ diagnostics, cohort-coverage-matrix integrator) and every literature probe we considered (Frobenius / scaled-Frobenius on D, rank correlation on triu(D), cophenetic correlation, Baker's γ, Robinson–Foulds variants, KC distance with λ-blend, matching-cluster, weighted-RF, principal angles / Grassmann, von-Neumann entropy difference, NMI(k), ARI(k)). Selection rule for Part 2: at minimum one matrix-level probe on D(τ), one spectral probe on ρ(τ) or L̂, one partition-based probe on T(τ).

- **Part 2 — Empirical pilot of the global per-band trace measure.** Per-(patient, band) T_LRG distributions, per-band cohort `n_trace` counts, side-by-side comparison to the substrate d_S verdict from `final_verdict_table_with_controls.csv`. Includes VI(k) recast as a multiscale band × k heatmap reported as a Section-5 result per se (not collapsed over k).

- **Part 3 — Empirical pilot of the localization probe.** Cohesion-CBR re-anchored on rest_pre + triplet-consensus CBR; per-(patient, band) trace-subtree catalogs with leaf membership, sEEG probe assignment, dendrogram height, and brain-space layout.

- **Part 4 — Cross-validation between global and localization layers.** Bidirectional. Pairs flagged by the global probe must concentrate inside subtrees flagged by the localization probe, and vice versa.

- **Part 5 — Recommendation.** What Section 5 of the manuscript should cite as the load-bearing global probe and the load-bearing localization probe, with quantitative justification per band.

τ stays fixed at 1/λ_max throughout. No probe relies on a privileged dendrogram cut n*. Per-patient distributions are reported alongside every cohort scalar — the trace is fundamentally a per-patient quantity and cohort verdicts read as "signals toward a band-level trace" rather than universal claims. Inter-patient variability is itself a finding.

## Substrate baseline (Section 4 verdict, mirrored here for reference)

From `data/audit/raw_fc_phase_distance/final_verdict_table_with_controls.csv` and `.agents/reports/2026-04-29_result-1-raw-fc-phase-trace.md`. Cohort N=10 (Pat_14 vendor-replaced 2026-04-25). Trace direction = T_d < 0.

| band | n_trace d_S | drift R² | rsPost null ratio | xb n_trace | substrate verdict |
|:---:|:---:|:---:|:---:|:---:|:---|
| β | 7/10 | pass | pass | 7 | trace (controls pass) |
| α | 8/10 | 0.20 caveat | pass | 7 | trace (1 caveat) |
| low_γ | 7/10 | pass | pass | 7 | trace (controls pass) |
| δ | 6/10 | pass | pass | 6 | soft trace (xb=6 caveat) |
| high_γ | 5/10 | pass | borderline | 5 | borderline |
| θ | 3/10 | fails | fails | 3 | drift-only |

Section 5 must mirror this table at the LRG level for every Part-2 candidate probe and report the deltas.

---

## Part 1 — Methods review

**Most candidate LRG-level phase-comparison probes have already been tried in the past month, and most came back negative or shelved. This is ballast, not a setback — Section 5's recommendation is grounded in a documented repertoire rather than speculation, and the residual surviving probes (the L1 continuous-trace family, the L5(k) VI(k) multiscale heatmap, the L4/L6/L7 leaf-level CBR-CNP-MSPC family) are the natural Part-2 / Part-3 candidates. The probes that came back negative (FTD, RSA, E1 Grassmann, Ψ_L, MRL strict-J) appear in the table for honesty and feed Part 5's rationale by elimination, not Part 2 directly.**

### 1.1 Geometric-ladder framing

Section 5 inherits the cohort-coverage-matrix's seven-rung decomposition of the LRG geometric ladder (per `.agents/guides/task-persistence-investigation/2026-04-29_decision-rules.md`). L2 (entropy curve `S(τ)` / specific heat `C(τ)`) is permanently removed from the ladder for our continuous-spectrum outlier case; L0 (residual subspace) is closed by the 2026-04-28 RSA diagnostic. The remaining seven rungs each consume one of {ρ(τ), D(τ), T(τ), L̂ eigenstructure, partitions of T} and carry one or more candidate probes. Section 5's Part-2 picks at minimum one matrix-level probe (acting on D(τ)), one spectral probe (acting on ρ(τ) or L̂), and one partition-based probe (acting on T(τ) treated multiscale).

| rung | primitive | repo-resident probes | role in Section 5 |
|:---:|:---|:---|:---|
| L1 | D(τ), ρ(τ) (matrix-level) | CTM (continuous trace matrix); raw-FC d_S/d_P/d_F at substrate level (Section 4) | matrix-level Part-2 candidate; substrate analogue is Section 4's headline |
| L3 | T(τ) (tree-distance) | KC distance with λ-blend (`tree_distance.kc_distance`); FTD (`functional_tree_distance`, integrated over τ — shelved) | tree-distance triangle Part-2 candidate; FTD already negative |
| L4 | T(τ) (subtree-leafset) | Trace-Modules (audit_15, strict J=0.9, cohort-null); MRL (audit_07, superseded); Consensus-subtree (scope only) | enumeration probe; null at strict gate, secondary at relaxed gates |
| L5(k) | T(τ) (integer-k flat partition) | VI(k) full profile (`compute_imcoh_vi.py`); ΔVI(k) multiscale heatmap (`h2_partition_multiscale_raw.csv`) | partition-based Part-2 candidate as multiscale band × k heatmap |
| L5(h_rel) | T(τ) (fractional-depth partition) | ΔVI(h_rel) (`dvi_hrel_n10_imcoh_abs.csv`) | partition-based Part-2 secondary; partition-resolution-locked for δ |
| L6 | T(τ) (per-leaf cophenetic vector) | Cohesion-CBR (audit_12), CNP (audit_13) | localization Part-3 primary |
| L7 | T(τ) (per-leaf, per-scale flat partition) | MSPC (audit_14) | localization Part-3 secondary; multiscale strip |
| (closed) | residual subspace of D | RSA (α_k, β_M; killed 2026-04-28) | negative; cited in Part 5 |
| (closed) | leading-k eigenspace of L̂ | E1 Grassmann (`audit_27`; negative 2026-04-29) | negative; triangle form not yet computed (Part-2 fresh candidate) |
| (out of ladder) | Ψ(n; τ) | partition selector — DEAD (bimodal) | constraint — confirms no privileged n* |
| (out of ladder) | C(τ) susceptibility | continuous-spectrum case — non-informative | permanently removed, never used |

### 1.2 Repo-resident probe inventory

Each row carries: scope-report file (under `.agents/guides/task-persistence-investigation/`), implementation script (under `scripts/`), output cache (under `data/`), what it consumes, the per-(p, b) yield, the per-band cohort scalar, the multiscale axis (if any), and the current verdict status.

| probe | scope | implementation | output | consumes | yield per (p, b) | cohort scalar | multiscale | current status |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| **legacy CBR** (audit_08) | none | `audit_08_per_patient_hierarchy.py` ✓ | `data/audit/per_patient_hierarchy/` | D(τ) ultrametric, Jaccard | per-leaf 4-state classification | n_T-leaves(p, b) at hard τ | size sweep 5..60 | **Superseded** by Cohesion-CBR; hard-threshold null in 4/9 patients |
| **task-anchored CBR** (TA) | `2026-04-25_task-anchored-cbr.md` | none | none | D, Jaccard | as legacy, anchor=task_test | as legacy | size sweep | **Abandoned** (superseded by Cohesion before implementation) |
| **containment-CBR** (B) | `2026-04-25_containment-cbr.md` | none | none | D, directional containment | per-leaf 4-state | n_T-leaves at hard τ | size sweep | **Abandoned** (superseded by Cohesion) |
| **Cohesion-CBR** (E) | `2026-04-25_cohesion-cbr.md` | `audit_12_cohesion_cbr.py` ✓ | `data/audit/per_patient_hierarchy_cohesion/{anchor_classification.csv, *.pdf}` | D, Jaccard + tightness, soft 4-corner affinity | per-leaf soft affinity α_p(ℓ), confidence | leaf-fraction TRACE per (p, b) | size sweep | **Active** — production probe; rest_post anchored; Section 5 will re-anchor on rest_pre |
| **consensus-subtree** (D) | `2026-04-25_consensus-subtree.md` | none | none | T leafset families | per-leaf 4-state via exact equality | n_T at exact-match | size sweep | **Abandoned** (strictest variant; absence-test paused) |
| **MRL** | `2026-04-25_module-retention-landscape.md` | embedded in early MRL pipeline | `data/audit/per_patient_hierarchy/` (legacy coloured dendrograms) | D, Jaccard via match_τ | scalar field M̄(b, h_rel) | cohort-null | h_rel sweep | **Superseded** 2026-04-25 (replaced by Cohesion-CBR + CTM) |
| **MRL ↔ CBR reconciliation** | `2026-04-25_mrl-vs-cbr-reconciliation.md` | analysis only | `data/reports/imcoh_mrl/` | MRL + Cohesion-CBR comparative | qualitative | — | — | **Reference** (post-mortem) |
| **CNP** (cophenetic neighbourhood) | `2026-04-25_cophenetic-neighbourhood.md` | `audit_13_cnp.py` ✓ | `data/audit/per_patient_hierarchy_cnp/{leaf_assignment.csv, Pat_XX/*.pdf}` | per-leaf cophenetic distance vector | per-leaf soft affinity on (ρ_pre_post, ρ_task_post) | leaf-fraction TRACE per (p, b) | scale-collapsed via Spearman | **Active**, uncontrolled cohort-wide |
| **Continuous Trace Matrix (CTM)** | `2026-04-26_continuous-trace-matrix.md` | `continuous_trace_matrix.py` (hypothesis_tests) ✓ | `data/reports/imcoh_continuous_trace/` | D upper-triangle | per-(p, b) Spearman ρ(Δ_task, Δ_rest); per-pair sign σ ∈ {-1, 0, +1} | n_trace(b) of patients with ρ > drift floor | matrix-level | **Active controlled** for α/β/low_γ; **fails split-baseline** for δ/high_γ |
| **MSPC** | `2026-04-26_multiscale-partition-coherence.md` | `audit_14_mspc.py` ✓ | `data/audit/per_patient_hierarchy_mspc/{multiscale_assignment.csv (407k rows), leaf_assignment.csv, *.pdf}` | flat partitions c_k(·) at every k ∈ [2, ⌊N/2⌋]; cluster-mate Jaccard | per-(leaf, scale) soft affinity α_p(ℓ, k) | per-band cohort fraction TRACE cells | k-resolved | **Active**, uncontrolled cohort-wide |
| **Trace-Modules** (audit_15) | `2026-04-25_trace-modules.md` | `audit_15_trace_modules.py` ✓ | `data/audit/trace_modules/trace_subtrees_n10_imcoh_abs.csv` | D, Jaccard at J_min=0.9 | per-(p, b, k, subtree) row | n_T(b) ≥ 1/10 | integer-k cuts | **Active** — strict-Jaccard cohort-null (1 row total); confirms strict-identity null |
| **Functional Tree Distance (FTD)** | `2026-04-28_functional-tree-distance.md` | `diag_functional_tree_distance.py` ✓ | `data/cache/functional_tree_distance/summary_n10_imcoh_abs.csv` | ρ(τ) and D(τ) over τ-grid [1/λ_max, τ*]; split-half null | δ_S, δ_P integrated over log-τ | n_trace per band | τ-integrated | **Shelved** 2026-04-28 — 5/60 cohort cells passed rank gate, 0/6 bands ≥3/10 |
| **Residual Subspace Alignment (RSA)** | `2026-04-28_residual-subspace-trace.md` | `diag_residual_subspace.py` ✓ | `data/audit/residual_subspace/{alpha_beta_grid.csv, diagnostic.pdf}` | residual eigenvectors of D^task − D^pre vs D^post − D^pre | α_1, α_3, α_5 canonical correlations; β_R rank fraction | — | k-grid {1,3,5} | **Killed** 2026-04-28 — α_1 below null, β_R median 0.88 → rank-1 trivializes α_1 |
| **E1 Grassmann** | `2026-04-29_e1-spectral-subspace-alignment.md` | `audit_27_e1_subspace_alignment.py` ✓ | `data/audit/spectral_subspace/e1_subspace_alignment_n10_imcoh_abs.csv` | top-k eigenvectors of L̂ at τ=1/λ_max | Δ_E1(k) = d(V_test, V_post) − d(V_pre, V_post) per (p, b, k) | frac_pos(b, k) | k-grid {2,3,5,8,13,21} | **Negative** 2026-04-29 — no band reaches 0.8 cohort fraction at any k; θ pushes eigenspace away from post |
| **Ψ(n; τ) partition selector** | — | `diag_psi_tau_scan.py` ✓ | `data/audit/psi_tau_scan/{psi_grid.csv, psi_diagnostic.pdf}` | dendrogram merge tree | argmax_n Ψ, Ψ_L | — | τ-sweep | **Dead** 2026-04-28 — bimodal landscape, no mesoscale gap → no privileged n* |
| **Edge-vs-Hierarchy discriminability (EH)** | `2026-04-28_edge-vs-hierarchy-discriminability.md` | `audit_2X_edge_vs_hierarchy.py` (TBD) | TBD | triu(A), triu(M), T | (d_E, d_H^S, d_H^KC) Z-scores per (p, b, phase-pair) | n_+ cohort | — | **Draft** — premise audit, not a primary trace measure |
| **Four-phase τ anchors** | — | `audit_17_taumin_four_phase.py`, `audit_18_taustar_four_phase.py`, `audit_19_lambda_idx_four_phase.py` ✓ | `data/audit/four_phase_{taumin, tau_star, lambda_idx1}/` | 4 phases LRG | per-band τ-flatness diagnostics | — | τ-grid | **Active diagnostic** — confirms τ=1/λ_max flat across phases |
| **VI(k) full profile** | `2026-04-24_h1-h4-vi-results.md` (canonical results); decision-rules L5(k) | `compute_imcoh_vi.py` (batch) ✓ | `data/reports/imcoh_vi/{vi_raw_profiles.csv, hypothesis_contrasts.csv, h2_partition_multiscale_raw.csv}` | dendrograms, fcluster at every k | VI per phase-pair per (p, b, k); ΔVI(k) | per-(b, k) cohort frac_pos | k-resolved | **Active** — only surviving cohort partition probe (α p=0.008, β p=0.014 under H2c continuous controls); recast as triangle in Part 2 |
| **VI(h_rel) full profile** | decision-rules L5(h_rel) | `compute_imcoh_vi.py` ✓ | `data/audit/dvi_split_baseline/dvi_hrel_n10_imcoh_abs.csv` | fcluster_at_h_rel | ΔVI(h_rel) per (p, b, h_rel) | per-(b, h_rel) cohort frac_pos | h_rel grid | **Active** — partition-resolution-locked vs L5(k) for δ (per decision-rules) |
| **Cohort coverage matrix** (integrator) | `2026-04-29_cohort-coverage-matrix.md` | `audit_18_cohort_coverage_matrix.py` (TBD per rebuild plan §10 step 3) | `data/audit/cohort_coverage_matrix/{cohort_coverage_matrix_n10_imcoh_abs.csv, triangulation_n10_imcoh_abs.csv}` | all rungs above | per-cell V(b, r), per-band T(b) | meta | multi-rung | **In progress** — v1.2 covers L1/L4/L4_aux/L5_k/L5_hrel/L6/L7; L3 deferred (KC λ-sensitivity TBD) |

### 1.3 Literature-only probes (not yet implemented)

These are candidates that the Part-1 review must consider but that no script in the repo currently runs as a phase-comparison triangle. The library helpers for the first three exist; the rest would need fresh code.

| literature probe | family | library coverage | per-(p, b) yield | multiscale-honest | traceability | viability for Part 2 |
|:---|:---|:---|:---|:---|:---|:---|
| Frobenius `‖D^a − D^b‖_F` | matrix on D | `numpy` only | scalar per phase-pair | no — single τ | none | weak (norm-driven, no rank) |
| Scaled Frobenius `‖D^a − D^b‖_F / ‖D^a‖_F` | matrix on D | `numpy` only | scalar per phase-pair | no | none | weak |
| **Spearman rank on triu(D)** | matrix on D | `scipy.stats.spearmanr`; treated as 1 − ρ | scalar per phase-pair | no | per-pair contributions yes | **strong** — direct LRG analogue of substrate d_S; Part-2 candidate |
| Cophenetic correlation between T^a, T^b | tree (height) | `scipy.cluster.hierarchy.cophenet` | scalar per phase-pair | no — collapses tree | per-pair | medium — treats tree as cophenetic distance vector |
| Baker's γ | tree (height) | one-line: Spearman of cophenetic vectors | scalar per phase-pair | no | per-pair | medium |
| Robinson-Foulds (unweighted) | tree (topology) | not in repo | scalar per phase-pair | no | per-internal-edge | weak — too coarse on ultrametric trees |
| Weighted RF | tree (topology + heights) | `tree_distance.weighted_rf_distance` ✓ | scalar per phase-pair | no | per-internal-edge | medium — coded but not yet cast as triangle |
| KC distance (Kendall-Colijn, λ-blend) | tree (topology + heights) | `tree_distance.kc_distance` ✓ | scalar per phase-pair per λ | λ-sweep is multiscale-adjacent | per-leaf-pair via m, M vectors | **strong** — wired into cohort matrix L1; not yet cast as triangle for Section 5; Part-2 candidate |
| Matching-cluster distance (Bogdanowicz-Giaro) | tree (cluster set) | `tree_distance.matching_cluster_distance` ✓ | scalar per phase-pair | no | per-internal-cluster | medium — coded, not yet used |
| **Principal angles / Grassmann** on top-k eigenvectors of L̂ | spectral on ρ(τ) | `scipy.linalg.subspace_angles` + custom; E1 audit covers pairwise; **triangle form not yet computed** | per (p, b, k) | k-grid is multiscale-adjacent | per-eigenmode | **strong** — Part-2 fresh candidate (triangle of pairwise E1) |
| von Neumann entropy difference | spectral on ρ(τ) | `lrg_eegfc.utils.lrg.*` (entropy primitives) | scalar per (p, b, phase) | no — single τ | none | weak (continuous spectrum case; would need re-derivation) |
| NMI(k), ARI(k) | partition on T | `sklearn.metrics`; existing partition loops | per-(p, b, k) | k-resolved | per-leaf | medium — equivalent information to VI(k); use VI as canonical |

### 1.4 Selection for Part 2

**Part-2 candidate slate (5 probes; one per family + the partition heatmap as a result per se):**

| family | probe | rationale |
|:---:|:---|:---|
| matrix on D(τ) | **Spearman rank on triu(D) recast as triangle** T_D-rank(p, b) = (1 − ρ(triu D^TT, triu D^RPost)) − (1 − ρ(triu D^RPre, triu D^TT)) | direct LRG analogue of substrate d_S; mirrors Section 4's compute scaffold |
| matrix on D(τ) | **CTM σ-aggregate** — n_trace_pairs / n_total_pairs | re-mine of `imcoh_continuous_trace/`; controlled at α/β/low_γ; fails δ/high_γ split-baseline (documented) |
| spectral on L̂ | **Top-k eigenspace overlap as a triangle** T_E1(p, b, k) = d_chord(V^TT, V^RPost; k) − d_chord(V^RPre, V^TT; k) at k ∈ {3, 5, 8, 13} | recast E1 audit's pairwise frac_pos into the triangle direction; the triangle has not been computed and may surface structure pairwise missed |
| tree-distance | **KC distance λ-blend triangle** T_KC(p, b, λ) = d_KC(λ; T^TT, T^RPost) − d_KC(λ; T^RPre, T^TT) at λ ∈ {0, 0.25, 0.5, 0.75, 1.0} | new; library exists; per the decision-rules L3 spec |
| partition (T) | **VI(k) multiscale band × k heatmap** T_VI(p, b, k) = VI(c^TT, c^RPost; k) − VI(c^RPre, c^TT; k) per (p, b, k); reported as per-patient heatmap + cohort `n_trace(b, k)` heatmap, no integration over k | the only surviving cohort partition probe at controlled n=10; framed as a Section-5 result per se; same-k cross-patient comparison flagged as a dirty operation but historically load-bearing |

**Out of Part 2 (negative or shelved):**
- FTD integrated form — already shelved with documented numbers; Part 1 row is sufficient; Part 5 cites the negative.
- RSA — killed; Part 1 row is sufficient.
- E1 Grassmann pairwise — negative; the triangle form (above) is the Part-2 carry-over.
- Ψ_L — dead; out of ladder.
- MSPC as a global per-band scalar — used as localization in Part 3 (per the user's prompt, MSPC is the secondary localization probe), not as a global probe.

**Out of Part 1 entirely (out of scope):**
- C(τ) entropy susceptibility (permanently removed for outlier case).
- Anything that requires picking n*.

## Part 2 — Empirical pilot of the global per-band trace measure

**Each of the five Part-1 candidate probes is run on the n=10 |ImCoh|_abs cohort at τ = 1/λ_max for all six bands and recast into the T_LRG triangle form. The per-(patient, band) distribution is reported alongside the per-band cohort scalar — no probe is collapsed to a cohort scalar without showing the underlying patient distribution. Each probe's per-band output is rendered side-by-side against the substrate d_S verdict. The probe that produces the cleanest cohort verdict on β + low_γ + α + δ AND either tightens α's drift caveat, raises δ above noise, or surfaces θ / high_γ structure the substrate missed becomes the headline; this selection is empirical and not pre-committed.**

### 2.1 Triangle definition

For every candidate probe with phase-pair distance `d_LRG(a, b)`, define
```
T_LRG(p, b) = d_LRG(task_test, rest_post) − d_LRG(rest_pre, task_test)
```
Negative T_LRG = trace direction. The triangle is asymmetric in the substrate sense (Section 4 audit_25's d_S) and inherits that asymmetry at the LRG level. Per-band cohort scalars:
- `n_trace(b) = |{p ∈ P : T_LRG(p, b) < 0}|`
- `T̃_LRG(b) = median_{p} T_LRG(p, b)`
- `n_trace_BH(b) = |{p ∈ P : T_LRG(p, b) < drift_floor(p, b)}|` for probes with a per-rung null

### 2.2 Probe-by-probe specification

#### 2.2.1 Spearman rank on triu(D) — D-rank triangle

**What it computes.** For each (p, b, φ), load the LRG ultrametric distance D^φ(τ=1/λ_max) = 1/ρ^φ(τ=1/λ_max) via `lrg_eegfc.workflow.lrg.load_lrg_result`. Stack triu(D^φ) into a per-phase pair-distance vector. Then
```
d_S^LRG(a, b) := 1 − Spearman(triu(D^a), triu(D^b))
```
Triangle: T_DS(p, b) = d_S^LRG(TT, RPost) − d_S^LRG(RPre, TT). Direct LRG analogue of Section 4's substrate d_S.

**Compute scaffold.** New audit script `scripts/01_compute/audit/audit_35_lrg_phase_distance.py`, modeled on `audit_25_raw_fc_phase_distance.py`. Reuses audit_25's V*(p, b) intersection logic (giant-component intersection across phases, after epileptic-node masking per `utils.io.patient.load_epileptic_nodes`) and Z-score machinery against a within-rest_pre split-half null (50 splits via Welch segments). (Numbering note: audit_28 through audit_31 are already taken by drift-triangle-null / lrg-inspection / rsPost-split-half / cross-baseline; this plan starts at audit_32.)

**Output.**
- `data/audit/lrg_phase_distance/{Pat_XX}/audit.csv` — per-patient per-(band, distance, phase-pair) aggregates.
- `data/audit/lrg_phase_distance/Td_per_patient_per_band.csv` — per-(p, b) T_LRG^DS scalar (10 × 6 × 3 = 180 rows; F/P/S as in substrate Td CSV).
- `data/audit/lrg_phase_distance/cohort_summary.csv` — per-band `n_trace`, `T̃_LRG`, drift_R², rsPost null ratio, xb counts mirroring substrate's column schema.
- `data/audit/lrg_phase_distance/final_verdict_table_with_controls.csv` — direct LRG analogue of substrate's final verdict table.

**Compute cost.** N=120 contacts × 119 / 2 ≈ 7140 pairs per phase × 4 phases × 10 patients × 6 bands × 51 distance-Z calculations = O(10⁷) Spearman calls. Trivial.

**Expected reading rule.** If T_DS LRG verdict matches substrate d_S verdict per band within ±1 patient, the LRG layer agrees with substrate (no sharpening, no contradiction). If T_DS sharpens substrate (β/low_γ moves from 7/10 → 9/10, α drift caveat resolves, δ moves above 7/10), the LRG layer is doing real work — the rank structure on D(τ) carries information the raw matrix did not.

#### 2.2.2 CTM σ-aggregate

**What it computes.** Re-mine of existing `data/reports/imcoh_continuous_trace/per_cell_summary_split.csv` outputs. CTM already computes per-pair sign agreement σ(i, j) ∈ {-1, 0, +1} on Δ_task = D^TT − D^RPre and Δ_rest = D^RPost − D^RPre. Per-(p, b):
```
n_trace_pairs(p, b) = |{(i, j) ∈ V*(p, b) : σ(i, j) = +1}|
n_total_pairs(p, b) = |V*(p, b)| × (|V*(p, b)| − 1) / 2
T_CTM(p, b) = − n_trace_pairs(p, b) / n_total_pairs(p, b)   (sign-flipped so trace = negative)
```
This is approximately equivalent to Spearman ρ on (Δ_task, Δ_rest) but at the per-pair sign level, not the rank level.

**Compute scaffold.** Re-mining script `scripts/01_compute/audit/audit_33_ctm_triangle.py` — reads CTM CSV outputs, casts σ-aggregates into the T_LRG triangle form, writes CSV. ~30 lines.

**Output.** `data/audit/ctm_triangle/{Td_per_patient_per_band.csv, cohort_summary.csv}`.

**Carry-over from prior verdict.** CTM's existing verdict (controlled α/β/low_γ pass split-baseline + cross-probe + drift-floor; δ/high_γ fail split-baseline shared-baseline). The triangle recast preserves this: bands that already failed CTM controls will produce a near-zero T_CTM cohort scalar; the existing controls are inherited.

#### 2.2.3 Top-k eigenspace overlap as triangle (Grassmann)

**What it computes.** For each (p, b, φ), load the leading-k eigenvectors V^φ_k of L̂(τ=1/λ_max) — already cached in the E1 audit's compute pipeline. Chordal Grassmann distance between two phases at fixed k:
```
d_chord(V^a_k, V^b_k) = sqrt(k − Σ_{i=1..k} σ_i²)
```
where σ_i are the singular values of (V^a_k)ᵀ V^b_k (canonical correlations). Triangle:
```
T_E1(p, b, k) = d_chord(V^TT_k, V^RPost_k) − d_chord(V^RPre_k, V^TT_k)
```
at k ∈ {3, 5, 8, 13}. Per-band cohort scalars per k; report all four panels, no integration over k.

**Compute scaffold.** New script `scripts/01_compute/audit/audit_37_e1_grassmann_triangle.py` — reuses audit_27's eigenvector loading, recasts as triangle. ~80 lines.

**Output.** `data/audit/e1_grassmann_triangle/{Td_per_patient_per_band_k.csv (10 × 6 × 4), cohort_summary_k.csv (6 × 4)}`.

**Carry-over from prior verdict.** E1 pairwise frac_pos was negative; the triangle form may differ because the asymmetry (TT-RPost vs RPre-TT) is not the same statistical question as pairwise frac_pos(b, k). The audit_30 verdict is genuinely unknown.

#### 2.2.4 KC distance λ-blend triangle

**What it computes.** For each (p, b, φ), build the KC vectors (m, M) per `lrg_eegfc.utils.metrics.tree_distance.kc_vectors(Z^φ)` from the linkage matrix Z^φ. KC distance with λ-blend:
```
d_KC(λ; T^a, T^b) = sqrt(Σ_pairs((1 − λ) (m^a − m^b)² + λ (M^a − M^b)²))
```
λ=0 = topology-only (m vector counts edges between leaf pairs); λ=1 = height-only (M vector carries the heights). Triangle:
```
T_KC(p, b, λ) = d_KC(λ; T^TT, T^RPost) − d_KC(λ; T^RPre, T^TT)
```
at λ ∈ {0, 0.25, 0.5, 0.75, 1.0}. Per-band cohort scalars per λ.

**Compute scaffold.** New script `scripts/01_compute/audit/audit_36_kc_triangle.py` — calls `kc_distance` directly. ~60 lines. **Note:** the cohort-coverage-matrix L3 spec deferred KC λ-sensitivity (per cohort-coverage scope §6 v1); audit_36 is the missing rung.

**Output.** `data/audit/kc_triangle/{Td_per_patient_per_band_lambda.csv (10 × 6 × 5), cohort_summary_lambda.csv (6 × 5)}`.

**Reading rule.** If T_KC at λ=0 (topology) and T_KC at λ=1 (heights) differ by band, the trace is in one but not the other — informative for diagnosing whether task induces topological vs metric reorganization in the dendrogram.

#### 2.2.5 VI(k) multiscale band × k heatmap (result per se)

**What it computes.** Re-mine of `data/reports/imcoh_vi/vi_raw_profiles.csv`. For each (p, b, k):
```
T_VI(p, b, k) = VI(c^TT_k, c^RPost_k) − VI(c^RPre_k, c^TT_k)
```
where c^φ_k = `fcluster(Z^φ, k, criterion='maxclust')`. **Do not collapse over k.** The result is a per-patient T_VI(b, k) heatmap and a cohort heatmap whose cell (b, k) carries
```
n_trace(b, k) = |{p : T_VI(p, b, k) < 0}|
```

**Framing.** Reported as a Section-5 result per se, not a step toward a scalar. The same-k cross-patient comparison is acknowledged in the body as a dirty operation — patients have different |L_p| (number of contacts) and slightly different dendrogram shapes, so cluster-id k=20 in Pat_02 is not literally the same partition as k=20 in Pat_06. But VI(k) profiles have historically guided the analysis (the H1-H4 VI results, the cohort-coverage L5(k) ridge at δ k=20–32, the Section-4 substrate d_S verdict), and any cohort-wide regularities (vertical bands at specific k where ≥7/10 patients trace; band-specific k-windows where the heatmap is consistently negative) deserve documentation regardless of the dirty-comparison caveat.

**Compute scaffold.** Re-mining script `scripts/01_compute/audit/audit_32_vi_triangle_heatmap.py` — reads `vi_raw_profiles.csv`, casts triangles per (p, b, k), writes CSV + heatmap PDF. (Earliest free audit number above the 28-31 block.)

**Output.**
- `data/audit/vi_triangle_heatmap/T_VI_per_patient_per_band_per_k.csv` (10 × 6 × K_max ≈ 6000 rows).
- `data/audit/vi_triangle_heatmap/cohort_n_trace_band_k.csv` (6 × K_max).
- `data/outputs/figures/section_5_lrg_trace/vi_triangle_heatmap_cohort.pdf` — single figure: 6 rows (bands), x-axis k ∈ [2, K_max], color = n_trace(b, k) on viridis.
- `data/outputs/figures/section_5_lrg_trace/vi_triangle_heatmap_per_patient/{Pat_XX}.pdf` — per-patient panel: 6 rows (bands), x-axis k, color = T_VI(p, b, k) on diverging RdBu (centered at 0; negative = trace).

**Reading rule.** Vertical band of red on the cohort heatmap at a specific k = many patients trace at the same scale = the manuscript's "k-locked trace" claim for that band. Diffuse red across most of the k-axis = scale-coherent trace (better story). Red-to-blue transition along k = scale-specific reorganization (band-band difference is informative).

### 2.3 Per-band side-by-side comparison

For each candidate probe, emit one row per band with columns:

| band | substrate d_S n_trace | T_LRG n_trace | substrate verdict | LRG verdict | delta |
|:---:|:---:|:---:|:---|:---|:---|
| β | 7/10 | _ | trace | _ | _ |
| α | 8/10 | _ | trace (caveat) | _ | _ |
| low_γ | 7/10 | _ | trace | _ | _ |
| δ | 6/10 | _ | soft trace | _ | _ |
| high_γ | 5/10 | _ | borderline | _ | _ |
| θ | 3/10 | _ | drift-only | _ | _ |

Across the five probes, this generates 5 × 6 = 30 cells of LRG verdict + delta. The pattern across the table feeds Part 5's selection rule.

### 2.4 Figures emitted by Part 2

- **Per-patient T_LRG box-plots per probe.** One PDF per probe. 6 panels (bands), y-axis T_LRG, points = patients, box = cohort. Substrate d_S T_d cohort overlaid as a dashed line for visual comparison.
- **Probe × band cohort scalar table.** One PDF table-figure with 5 rows (probes) × 6 columns (bands), cell = `n_trace / 10`, color-coded to substrate verdict alignment (green = matches or sharpens, yellow = matches with delta ±1, red = disagrees by ≥2).
- **VI(k) cohort heatmap** (above).
- **VI(k) per-patient heatmaps** (above).

All PDFs vector, no `set_rasterized`, no `fig.suptitle`, figure-level legends via `lrg_eegfc.visuals.layout.figure_legend`.

### 2.5 Selection rule for Part 5 headline (pre-registered)

The headline global probe is selected by the following pre-registered rule (so the choice is not post-hoc tunable):

1. **Substrate-trace bands (β, α, low_γ, δ) must be at least matched.** A probe must produce `n_trace(b) ≥ substrate_n_trace(b) − 1` for b ∈ {β, α, low_γ, δ} to be eligible.
2. **Among eligible probes, prefer the probe that adds the most controlled cells.** A controlled cell = `n_trace ≥ 7/10` AND passes a probe-matched within-baseline null. The probe with the highest count wins.
3. **Tie-breakers in order:** (a) does it tighten α's drift caveat (R² < 0.10)? (b) does it raise δ to 7/10 with controls? (c) does it surface θ or high_γ above noise (n_trace ≥ 7/10 with controls)? (d) cleanest visual story in the per-patient box-plots (lowest cohort variance).

The VI(k) multiscale heatmap is a parallel result-per-se, not a competitor for the global-probe headline. It appears in Section 5 alongside the headline probe regardless of which probe wins (1) and (2).

## Part 3 — Empirical pilot of the localization probe

**The localization layer asks per-patient, per-band: which subtrees of the dendrogram "moved closer to task in rest_post than in rest_pre"? Cohesion-CBR is the primary localization probe — already implemented (`audit_12_cohesion_cbr.py`), already cached at `data/audit/per_patient_hierarchy_cohesion/`, but anchored on rest_post. Section 5 needs it re-anchored on rest_pre to match the triangle-form question, and a new triplet-consensus variant that is bias-free relative to phase anchor. MSPC (`audit_14`, already cached) is the secondary localization probe — per-leaf trace fraction across scales.**

### 3.1 Cohesion-CBR re-anchored on rest_pre

**Compute scaffold.** Clone `scripts/01_compute/audit/audit_12_cohesion_cbr.py` into `scripts/01_compute/audit/audit_38_cohesion_cbr_rpre_anchor.py`. The only structural change is which dendrogram drives the anchor enumeration: instead of iterating over internal nodes of T^RPost, iterate over internal nodes of T^RPre. For each rest_pre subtree S:

```
for each internal node v of T^RPre with size_min ≤ |leaves(v)| ≤ size_max:
    L = leaves(v)
    h_rpre = native height of v in T^RPre
    J_rpre_tt = best-Jaccard match of L against any subtree of T^TT
    J_rpre_rpost = best-Jaccard match of L against any subtree of T^RPost
    J_tt_rpost = best-Jaccard match of (best-tt-match of L) against any subtree of T^RPost
    cohesion_pre_to_task = soft 4-corner affinity on (J_rpre_tt, J_rpre_rpost) corners (TRACE corner: low J_rpre_tt low J_rpre_rpost; PERSIST: high high; etc.)
    cohesion_task_to_post = soft 4-corner affinity on (J_rpre_tt, J_tt_rpost) corners
    classify subtree as TRACE if cohesion_pre_to_task is dominant TRACE corner AND cohesion_task_to_post is dominant PERSIST corner
```

The two-stage cohesion is the operationalization of "moved away from rpre during task AND stayed away in rpost". The existing audit_12 machinery handles the corner-distance affinities; what changes is the anchor and the predicate.

**Output.**
- `data/audit/per_patient_hierarchy_cohesion_rpre/{Pat_XX}/{band}_hierarchy-crossphase.pdf` — 4-phase dendrogram with rest_pre subtrees colored by classification.
- `data/audit/per_patient_hierarchy_cohesion_rpre/anchor_classification.csv` — per-(patient, band, anchor_size, h_rpre, leafset_id) classification with leaf membership.
- `data/audit/per_patient_hierarchy_cohesion_rpre/trace_subtrees_n10_imcoh_abs.csv` — filtered subset where classification = TRACE; per-row: patient, band, h_rpre, |L|, leaf_ids, sEEG probe assignment, trace strength (cohesion confidence).

**Per-patient deliverable per (p, b).** A list of 0..K trace subtrees with leafsets; the dendrogram with those subtrees highlighted; the brain-space view via `lrg_eegfc.visuals.spatial.view_brain_connectome` zoomed to the bbox of the trace subtrees' channels (per the never-always rule on glass-brain zoom).

### 3.2 Triplet-consensus CBR

**Compute scaffold.** New script `scripts/01_compute/audit/audit_39_triplet_consensus_cbr.py`. Bias-free relative to which phase is anchored. For each candidate leafset L (enumerated over the union of rest_pre, task_test, rest_post internal nodes):

```
for each L in union(internal-leafsets(T^RPre), internal-leafsets(T^TT), internal-leafsets(T^RPost)):
    if size_min ≤ |L| ≤ size_max:
        c_RPre_L = best-Cohesion-score of L in T^RPre
        c_TT_L = best-Cohesion-score of L in T^TT
        c_RPost_L = best-Cohesion-score of L in T^RPost
        TT_RPost_score = c_TT_L * c_RPost_L
        RPre_TT_score = c_RPre_L * c_TT_L
        RPre_RPost_score = c_RPre_L * c_RPost_L
        if TT_RPost_score > max(RPre_TT_score, RPre_RPost_score) AND TT_RPost_score > tau_consensus:
            classify L as TRIPLET-TRACE
```

The product `TT_RPost_score = c_TT_L * c_RPost_L` is high iff L is a coherent subtree in BOTH task and post. The "trace" predicate is then "more coherent in (task, post) than in (rpre, *)". `tau_consensus = 0.5` (pre-registered).

**Output.**
- `data/audit/per_patient_hierarchy_triplet_consensus/{Pat_XX}/{band}_hierarchy-crossphase.pdf` — 4-phase dendrogram with triplet-trace subtrees highlighted across all three phase-trees (the same leafset shows in all three trees if it survives).
- `data/audit/per_patient_hierarchy_triplet_consensus/triplet_trace_subtrees.csv` — per-row: patient, band, leafset_id, |L|, c_RPre, c_TT, c_RPost, TT_RPost_score, leaf_ids, sEEG probe assignment.

**Per-patient deliverable.** As §3.1 above, but the highlight is consistent across all three phase dendrograms (same leafset, same color). This is the figure that makes the "task-induced and post-task-retained module" claim visually undeniable when the underlying numbers support it.

### 3.3 MSPC — per-leaf trace fraction (secondary)

**Compute scaffold.** Re-mining script `scripts/01_compute/audit/audit_34_mspc_leaf_trace_fraction.py` — reads existing `data/audit/per_patient_hierarchy_mspc/multiscale_assignment.csv` (407k rows, already cached). Per-leaf:
```
f_trace(p, b, ℓ) = |{k ∈ k_grid : dominant(p, b, ℓ, k) = TRACE}| / |k_grid|
```
Per-band cohort summary: distribution of `f_trace` across leaves and patients. Per-leaf threshold `f_trace ≥ 0.5` flags leaf as a TRACE-leaf. Per-patient TRACE-leaf list with sEEG probe assignment.

**Output.**
- `data/audit/mspc_leaf_trace_fraction/leaf_trace_fraction.csv` — per-(p, b, ℓ) f_trace + flag.
- `data/audit/mspc_leaf_trace_fraction/cohort_summary.csv` — per-band cohort fraction of TRACE-leaves.

**Connection to §3.1/§3.2.** MSPC's TRACE-leaves should be a superset of the leaves contained in Cohesion-CBR's flagged subtrees (the per-leaf granularity is finer than the per-subtree granularity per the MSPC ↔ CBR dichotomy in §10 of the MSPC scope). If they disagree (a TRACE-leaf in MSPC is not inside any flagged Cohesion-CBR subtree), that leaf is a "distributed-reorganization" leaf — task changed its neighbours but the new neighbour set is not a coherent subtree. Report the distribution per band.

### 3.4 Brain-space localization figures

For each (patient, band) cell with at least one flagged trace subtree (from §3.1 or §3.2), emit a 2-panel figure:

- **Left panel:** the dendrogram of the anchor phase (rest_pre for §3.1, the consensus subtree's natural phase for §3.2) with the trace subtree highlighted in `#d62728` (the user-pinned trace color). Codebar at top per the audit_08 layout.
- **Right panel:** nilearn glass brain (multi-view: sagittal + coronal + axial) zoomed to the electrode bbox of the trace subtree's channels. Edge thresholding reuses existing `view_brain_connectome` defaults; channels colored by trace-strength on a viridis gradient.

Output per `(p, b)` as `data/outputs/figures/section_5_lrg_trace/localization_per_patient/{Pat_XX}_{band}.pdf`.

### 3.5 Cross-cohort regularity panel

A second figure aggregating across patients per band: which sEEG probes (anatomical regions, via the implant CSV `data/raw/stereoeeg_patients/Pat_NN/implant_pat_NN.csv` Desikan-Killany column — never the channel-letter prefix per never-always rule) most often participate in trace subtrees? Per band, a horizontal bar chart of "n_patients with at least one trace subtree containing a contact in this region" — top 10 regions per band. Output: `data/outputs/figures/section_5_lrg_trace/cross_cohort_region_regularity.pdf` (6 facets, one per band).

This panel is the cohort-level interpretive layer that turns the per-patient subtree catalog into a band-specific regional regularity statement (within the framing constraint that variability is itself a finding).

## Part 4 — Cross-validation between global and localization layers

**A pair `(i, j)` flagged as trace-direction by the global probe (Part 2) should sit, with elevated probability, inside a subtree flagged by the localization probe (Part 3). Conversely, the bottom decile of pair-level T_LRG values (the most trace-leaning pairs) should concentrate inside flagged subtrees. Both directions must agree per (patient, band) and cohort-aggregated. Disagreements are diagnostic, not hidden — they tell us which probe is leaking. The Part-4 cross-validation is internal to the LRG layer (LRG global probe ↔ LRG localization probe); per the user's directive on per-pair scope, no per-pair substrate distance recomputation is in scope.**

### 4.1 Per-pair contributions from the global probe

For the headline global probe selected by Part 2's pre-registered rule (§2.5), compute per-pair T_LRG_pair(p, b, i, j) at τ=1/λ_max. The exact form depends on which probe is the headline:

- **If headline = D-rank triangle (§2.2.1):** T_DS_pair(p, b, i, j) = rank(D^TT(i, j) − D^RPost(i, j)) − rank(D^RPre(i, j) − D^TT(i, j)) computed against the per-band rank distribution. Negative = the (i, j) pair tracks the trace direction. Equivalently, the Spearman ρ on triu(D) decomposes into per-pair rank-difference contributions; that decomposition is the per-pair signal.
- **If headline = CTM σ-aggregate:** T_CTM_pair = − σ(p, b, i, j) ∈ {−1, 0, +1}; trivially per-pair (already cached in `imcoh_continuous_trace/`).
- **If headline = E1 Grassmann triangle:** the per-pair quantity is undefined (Grassmann distance is a global property of an eigenspace pair; it cannot be decomposed into per-pair contributions). In this case, fall back to the second-place probe in Part 2's selection rule for the per-pair role. (This is a documented limitation of the spectral family.)
- **If headline = KC λ-blend triangle:** at λ=0, per-pair contribution is the per-pair difference of m vectors (path-length between the two leaves through the tree); at λ=1, per-pair difference of M vectors. Pick the λ that produced the cleanest cohort verdict in §2.2.4.
- **If headline = VI(k) heatmap:** VI is a partition-level statistic; per-pair decomposition is via mutual-information per pair-ish, but the natural per-pair object here is the same-cluster indicator `1[c^φ(i) = c^φ(j)]`. The Hamming difference of these indicators across phase pairs is a per-pair triangle. (Used only if VI(k) becomes the headline — which is unlikely given the heatmap is a result-per-se, not a competitor.)

**Output.** Per-(p, b) `pair_contributions_lrg.npz` keyed by `(i, j)` carrying `T_LRG_pair`. Roughly 10 patients × 6 bands × O(N²)≈O(10⁴) pairs per cell ≈ 600k pair-rows total. Compute trivial.

### 4.2 Forward direction: enrichment inside flagged subtrees

For each (p, b) with at least one flagged trace subtree S_local from Part 3 (Cohesion-CBR rest_pre anchor or triplet-consensus, whichever is the headline localization probe per Part 5's rule):

```
P_local(p, b) = {(i, j) : i, j ∈ ∪ S_local AND i ≠ j}
P_complement(p, b) = {(i, j) : i, j ∈ V*(p, b) AND (i, j) ∉ P_local}
enrich(p, b) = mean(T_LRG_pair | P_local) − mean(T_LRG_pair | P_complement)
```

Permutation null: shuffle subtree assignment 1000× by randomly permuting which leaves belong to S_local while preserving |∪ S_local|. Compute permutation distribution of `enrich`. Report `p_perm(p, b) = fraction of permutations with enrich_perm ≤ enrich_obs`.

Cohort aggregation per band:
```
n_enriched(b) = |{p : enrich(p, b) < 0 AND p_perm(p, b) ≤ 0.05}|
```
where `enrich(p, b) < 0` means trace-direction pairs are concentrated inside flagged subtrees.

### 4.3 Backward direction: concentration of bottom-decile pairs

For each (p, b), define
```
P_bottom(p, b) = {(i, j) : T_LRG_pair(p, b, i, j) is in bottom decile across all pairs in V*(p, b)}
concentration(p, b) = |P_bottom ∩ P_local| / |P_bottom|
```

Null expectation: if pairs were uniformly distributed across the network, `concentration ≈ |P_local| / |V*² / 2|`. Compute the same permutation null as in §4.2 and report `p_perm` for `concentration` exceeding null.

Cohort aggregation per band:
```
n_concentrated(b) = |{p : concentration(p, b) > null_expectation(p, b) AND p_perm(p, b) ≤ 0.05}|
```

### 4.4 Cross-direction agreement

Per (p, b), the two cells are:
- enriched? (forward direction, §4.2)
- concentrated? (backward direction, §4.3)

Both are diagnostic of "global probe and localization probe agree". The agreement matrix per band:

| | concentrated | not concentrated |
|:---:|:---:|:---:|
| **enriched** | strong cross-validation | concentration leaks (forward but not backward — localization may be too narrow) |
| **not enriched** | enrichment leaks (backward but not forward — global probe finds pairs that the localization missed) | no localization-global agreement (both probes are silent at this cell, or one of them is misclassifying) |

Per-band cohort scalar: `n_agree(b) = |{p : enriched AND concentrated}|`. For a band with a real per-patient trace, expect `n_agree(b) ≥ n_trace(b) − 2` (allowing 2 patients of slack for the localization probe's own miss-rate).

### 4.5 Compute scaffold and output

**Compute scaffold.** New script `scripts/01_compute/audit/audit_40_localization_global_crossvalidation.py`. Reads:
- `data/audit/lrg_phase_distance/Td_per_patient_per_band.csv` (or whichever Part-2 probe is headline) for `T_LRG_pair`.
- `data/audit/per_patient_hierarchy_cohesion_rpre/trace_subtrees_n10_imcoh_abs.csv` (or `triplet_trace_subtrees.csv`) for `S_local`.
- Patient layout from `data/audit/cohort_metadata.csv` and per-patient channel masks via `lrg_eegfc.utils.io.patient.load_epileptic_nodes` + `V*(p, b)` intersection logic from audit_25.

**Output.**
- `data/audit/localization_global_crossvalidation/per_patient_per_band_crossvalidation.csv` — per-(p, b): `enrich`, `p_perm_enrich`, `concentration`, `p_perm_concentration`, `enriched`, `concentrated`, `agree`.
- `data/audit/localization_global_crossvalidation/cohort_summary.csv` — per-band: `n_enriched`, `n_concentrated`, `n_agree`, `expected_n_agree (= n_trace − 2)`, `agreement_pass (boolean)`.
- `data/outputs/figures/section_5_lrg_trace/crossvalidation_panel.pdf` — 6-band panel: per band, scatter of (enrich, concentration) per patient, colored by agree/disagree, with permutation-null contour.

### 4.6 Handling disagreements

Disagreements are reported, not hidden. The interpretation rules:

- **Forward-only (enriched but not concentrated).** Localization probe is too narrow — the flagged subtrees miss most of the trace pairs. Action: report alongside the cohort scalar; flag the band as "localization-narrow"; the manuscript carries this as a caveat on the localization claim for that band.
- **Backward-only (concentrated but not enriched).** Global probe finds trace pairs scattered across the network — most are inside flagged subtrees but the within-subtree mean is washed out by the rest of the subtree. This is the rare regime where the localization is correct but the per-pair trace signal is sparse within the flagged region. Action: report; flag as "sparse-within-subtree"; the manuscript carries this as evidence the trace is module-localizing without being module-uniform.
- **Both null.** The localization probe found no subtrees, or the global probe is silent at this band. The cross-validation contributes nothing here; Part 5 falls back to the global probe's verdict alone.
- **Both positive at cohort level.** Strong cross-validation. The manuscript can claim "the trace is both visible at the network level AND localizable to specific subtrees" for this band.

## Part 5 — Recommendation

**The recommendation Section 5 will cite is selected by the pre-registered rules in §2.5 (global probe) and §3.1/§3.2/§4.4 (localization probe), grounded in the per-band tables produced by Parts 2–4. The recommendation respects the per-patient framing — cohort verdicts are signals toward a band-level trace, not universal claims — and acknowledges that the trace is fundamentally a per-patient quantity. Inter-patient variability (likely tied to implant geometry, see Pat_03 1024 Hz outlier flag) is itself a finding of interest, not noise to be aggregated away.**

### 5.1 Decision framework (filled in once compute lands)

The actual recommendation requires the per-band tables from Parts 2–4. The decision framework is pre-registered here:

#### Global probe selection
1. Apply §2.5 selection rule: among the five candidate probes from Part 1, pick the one that (a) matches substrate d_S verdict on β/α/low_γ/δ within ±1 patient, (b) controls the most cells against probe-matched within-baseline null, (c) tightens at least one of α drift / δ noise-floor / θ-or-high_γ recovery.
2. If two probes tie on (a)+(b), break ties by (c) priority order: α drift tighten > δ raise to controlled > θ recovery > high_γ recovery > cleanest visual story.
3. The headline global probe must be cited with: per-band T_LRG cohort scalar table, per-patient T_LRG box-plot figure, side-by-side comparison to substrate d_S, and an explicit statement of how it agrees with / sharpens / contradicts the substrate verdict per band.

#### Localization probe selection
1. **Cohesion-CBR re-anchored on rest_pre** is the prior frontrunner (the existing audit_12 machinery is mature and the rest_pre-anchor question matches the triangle direction).
2. **Triplet-consensus** is the supplement when the rest_pre anchor produces strong per-band candidates — the consensus subtree is the visually compelling object for the manuscript figure.
3. **MSPC f_trace per-leaf** is the cohort-level summary that complements the per-subtree catalog. It does not name modules but it provides the per-leaf scale-resolved field that lets the reader see which leaves trace at which scales.
4. The headline localization probe must be cited with: per-(patient, band) trace-subtree catalog, the per-patient brain-space figures, and the cross-cohort regional regularity panel from §3.5.

#### Cross-validation requirement
1. The headline global probe and the headline localization probe must produce `n_agree(b) ≥ n_trace(b) − 2` per band on the substrate-trace bands (β, α, low_γ, δ). Bands failing this requirement are flagged as "global-only" or "localization-only" in the manuscript text, never silently merged into the headline claim.

### 5.2 Documented negatives — why other paths were not chosen

The following probes appeared in Part 1's table and did not enter Part 2 / 5:

- **Functional Tree Distance (FTD), τ-integrated rank `δ_S` and height `δ_P`.** Shelved 2026-04-28 (`.agents/reports/2026-04-28_functional-tree-distance-verdict.md`). 5/60 cohort cells passed the rank gate at strict threshold; 0/6 bands ≥ 3/10. The integrated form does not surface a cohort signal at n=10. Section 5 does not pursue τ-integrated tree distances.
- **Residual Subspace Alignment (RSA), α_k canonical correlations on (D^task − D^pre, D^post − D^pre) eigenvectors.** Killed 2026-04-28 (`.agents/reports/2026-04-28_residual-subspace-diagnostic.md`). Mean α_1 sits below permutation null by −0.066 across cohort; β_R median 0.88 puts the residual matrix in the rank-1 regime where high α_1 is trivial. The probe does not separate trace from drift.
- **E1 Grassmann pairwise frac_pos.** Negative 2026-04-29 (`.agents/reports/2026-04-29_e1-cohort-verdict.md`). No band reaches 0.8 cohort fraction at any k ∈ {2, 3, 5, 8, 13, 21}; θ pushes the eigenspace away from post (anti-aligned). The triangle recast (§2.2.3) is in Part 2 because the asymmetry (TT-RPost vs RPre-TT) is a different statistical question than pairwise frac_pos, but if the triangle also returns negative the spectral family is closed.
- **Ψ(n; τ) partition selector.** Dead 2026-04-28. Bimodal landscape, no mesoscale gap — confirms that |ImCoh|_abs networks have no privileged dendrogram cut, which is why the entire Part-2/3 framework is constrained to multiscale-honest probes.
- **MRL strict-Jaccard module-retention scalar field.** Superseded 2026-04-25 by Cohesion-CBR + Continuous Trace Matrix. The (J_pre, J_post) cloud sits on the diagonal at high Jaccard; the trace zone is empty. The MRL ↔ CBR reconciliation post-mortem (`2026-04-25_mrl-vs-cbr-reconciliation.md`) is the canonical write-up.

These negatives are not an embarrassment — they constitute the audit trail that makes the Part 5 recommendation defensible. A reviewer asking "why not RSA / FTD / Ψ?" gets a per-probe verdict file with cohort numbers. The selection rule selects from a complete repertoire, not from a curated subset.

### 5.3 Framing constraints carried into the manuscript

The manuscript's Section 5 inherits the following framing constraints from the Section 4 substrate result and from the per-patient nature of the data:

1. **Per-patient T_LRG distributions are reported alongside every cohort scalar.** Box-plots are non-negotiable; cohort `n_trace(b)` without the underlying per-patient distribution is misleading.
2. **Cohort verdicts read as "signals toward a band-level trace", not universal claims.** Language like "the cohort shows" or "at the cohort level, β shows trace in 7/10 patients" is correct; "task induces a β-band trace" is over-claim.
3. **Inter-patient variability is itself a finding.** The bands where inter-patient variance is largest are flagged in the manuscript text and tied (where possible) to implant geometry via the cross-cohort regional regularity panel from §3.5.
4. **τ = 1/λ_max is the canonical choice.** Any τ-sweep that surfaces in supplementary material is framed as a robustness check, not as part of the headline.
5. **The trace / anchor / reset / emergent taxonomy** (per `.agents/guides/01_project/terminology.md`) is used throughout — `n_trace` not `n_persist`, "trace zone" not "persistence zone", etc.
6. **The β + low_γ + α + δ four-band reading is the load-bearing claim from Section 4** — Section 5 either reproduces or sharpens this list. Bands moving in or out of the trace list relative to substrate must be explicitly flagged with the LRG-vs-substrate delta.

### 5.4 What Section 5 does NOT claim

To be explicit about scope:

- Section 5 does NOT claim a band-specific neural mechanism — it characterizes geometry, not biology.
- Section 5 does NOT claim that the trace is universal across patients — the per-patient framing is non-negotiable.
- Section 5 does NOT claim that the LRG layer is necessary in addition to the substrate layer — the right reading depends on the empirical Part-2 deltas. If LRG only reproduces substrate, the claim is "LRG confirms substrate"; if LRG sharpens substrate, the claim is "LRG adds resolution"; if LRG contradicts substrate, the claim is "the layers carry different aspects of the same task-induced reorganization, both reported".
- Section 5 does NOT claim a τ-dependent multiscale story beyond τ=1/λ_max — the manuscript holds τ fixed and acknowledges in the methods note that the choice is the canonical Villegas-2025 default for our outlier (continuous-spectrum, fully-connected, weight-heterogeneous) FC case.

### 5.5 Implementation handoff

When this plan transitions from `status: draft` to `status: current`, the following must be in place:

- All five Part-2 audit scripts (`audit_32` VI heatmap, `audit_33` CTM triangle, `audit_35` D-rank, `audit_36` KC λ-blend, `audit_37` Grassmann) implemented and producing CSV outputs at the paths listed.
- The two Part-3 audit scripts (`audit_38` Cohesion-CBR rest_pre, `audit_39` triplet-consensus) plus the MSPC re-mining script (`audit_34`) implemented.
- The Part-4 cross-validation script (`audit_40`) implemented and producing the per-band agreement matrix.
- The figures listed in §2.4, §3.4, §3.5, §4.5 produced as PDFs at the paths listed.
- The plan's Part 5 §5.1 actual-recommendation cells filled in with the empirical numbers.
- The manuscript's Section 5 draft cross-references this plan's CSV outputs and figure paths verbatim.

The implementation effort estimate (from the plan-mode plan) is ~4 days of focused engineering work; compute is uncapped and is not the bottleneck.

---

## Working notes (will be removed before this becomes `status: current`)

### Structure deviation from the README scope-report template

This document is a methods-review meta-document, not a single-measure scope. It does not carry the per-measure `notation → definitions → properties → caveats → pseudocode → visualization → connection-to-prior-tools → implementation-plan → open-questions` template — that template applies to each row of Part 1's table, where the row points back to the originating per-measure scope report under this same folder.

### Pointers to the canonical compute scaffolds we will mirror

- Substrate analogue: `scripts/01_compute/audit/audit_25_raw_fc_phase_distance.py`. Section 5's LRG-side analogue should sit at `scripts/01_compute/audit/audit_2X_lrg_phase_distance.py` and reuse audit_25's V*(p, b) intersection logic + Z-score machinery as much as possible.
- Localization re-anchoring: `scripts/01_compute/audit/audit_12_cohesion_cbr.py`. The rest_pre-anchored variant clones the anchor-enumeration loop with a swapped anchor dendrogram.

### Library entry points (use; do not reinvent)

- `lrg_eegfc.workflow.lrg.load_lrg_result` — load ρ(τ), eigenstructure, dendrogram per (p, b, φ).
- `lrg_eegfc.utils.metrics.tree_distance` — `kc_distance` (λ-blend), `matching_cluster_distance`, `weighted_rf_distance`. KC is wired into cohort-coverage-matrix L1; matching-cluster and weighted-RF are coded but not yet cast as a T_LRG triangle.
- `lrg_eegfc.utils.metrics.compare_phases`, `compare_ultrametric_matrices` — phase-comparison helpers.
- `lrg_eegfc.utils.metrics.functional_tree_distance` — FTD primitives (negative verdict integrated form; reusable for single-τ triangle).
- `lrg_eegfc.utils.metrics.hypothesis` — `wilcoxon_z`, `bh_fdr`, `cluster_stats`, `boot_ci_mean`, `rank_biserial`.
