---
name: 2026-04-29_eigenvector-direct-pivot-plan
type: plan
era: COHORT_N10
status: active
created: 2026-04-29
updated: 2026-04-29
pointers:
  - .agents/reports/2026-04-29_critical-post-mortem.md
  - .agents/guides/02_methods/lrg-framework-guide.md
  - .agents/guides/task-persistence-investigation/2026-04-29_decision-rules.md
  - .agents/guides/task-persistence-investigation/2026-04-29_cohort-coverage-matrix.md
  - .agents/reports/2026-04-28_residual-subspace-diagnostic.md
---

# Eigenvector-direct pivot — plan for E1, E2, E3 rungs

## Renormalization head

**Add three new rungs (E1, E2, E3) to the cohort-coverage matrix that
probe diffusion-pattern reorganisation directly via the eigenvector
basis of the Laplacian, replacing the failed scalar-spectrum probes
(entropy curve, specific heat, Ψ-as-selector, rank-1 residual
subspace) that all collapse on our continuous spectrum. Each E-rung
uses the FULL top-k eigenvector basis at fixed τ = 1/λ_max (not
residuals, not scalar functions of eigenvalues), compares phases
through standard linear-algebra distances (principal angles for
subspace alignment, Euclidean distance for diffusion-map embedding,
cohort distribution shift for inverse participation ratio), and ties
into the within-baseline halves cache for nulls. Critical
methodological distinction from the dead residual-subspace probe:
that probe used rank-1 leading eigenvectors of *residuals* and was
killed by near-rank-1 dominance making alignment trivial; the
E-series uses rank-k eigenvectors of the *Laplacian itself* — different
math, different failure modes. Order of execution: cache the
eigenvectors first (one-line LRG cache extension), then E1
(cleanest), then E2 + E3 in parallel. Acceptance criterion: at least
one E-rung surfaces cohort-wide signal in α / β / γ_h that the current
dendrogram-based ladder declared ergodic — if all three E-rungs are
silent for those bands too, that is a much stronger negative result
than the current matrix gives.**

---

## TL;DR

1. **Three new spectral rungs replace the dead scalar probes.** E1 (subspace alignment), E2 (diffusion-map embedding drift), E3 (inverse-participation-ratio shifts).
2. **All three operate at fixed τ = 1/λ_max on the full eigenvector basis.** Not residuals, not scalar functions of eigenvalues. Different math from the dead probes.
3. **Step 0 is a one-line cache change.** Add eigenvectors to the LRG result NPZ. Re-cache n=10 × 6 bands × 4 phases × 2 halves ≈ 480 small files.
4. **Order: E1 first; if E1 yields nothing, E2 and E3 are exploratory only.** E1 is the cleanest direct test of "do diffusion modes rotate".
5. **Each rung is one scope report + one producer + one within-baseline null + one matrix-v1.2 integration step.** Six PRs total per rung; 18 PRs total but each is small.
6. **Pre-registered cohort threshold: ≥ 8 / 10 patients positive at sufficient k or τ-tilde and ridge ≥ 3 along the scale axis where applicable.**
7. **Acceptance gate: post-mortem decision at the v1.2 matrix.** If α / β / γ_h surface as positive at ≥ 1 E-rung, the eigenvector-direct hypothesis is supported and the paper headline pivots to the E-series. If silent, the v1.2 result strengthens the negative claim and the paper pivots to "we tested both tree-based and eigenvector-direct probes, neither gives cohort-wide multi-band trace beyond δ-partition-resolution-locked".

---

## 1. Notation and primitives

**Abbreviations.** *L̂* = Laplacian = D̂ − Â (D̂ degree matrix, Â functional-connectivity weight matrix; weighted, undirected, fully connected, N nodes per patient, 100 ≤ N ≤ 122 across cohort). *λᵢ* = i-th eigenvalue of L̂, ordered 0 = λ₁ < λ₂ ≤ … ≤ λ_N. *vᵢ* = corresponding eigenvector (column of the spectral basis V, N × N). *τ* = diffusion time, fixed at τ = 1/λ_max throughout (per the framework guide). *K̂(τ)* = exp(−τL̂) = V diag(exp(−τλᵢ)) V^T (heat kernel propagator). *ρ̂(τ)* = K̂(τ) / Tr(K̂(τ)) (density operator). *D* = communication distance, Dᵢⱼ = 1 / ρ̂ᵢⱼ. *Top-k eigenvectors* = the k eigenvectors corresponding to the k smallest non-zero eigenvalues (i.e. v₂, v₃, …, v_{k+1}; the trivial v₁ is dropped).

For each (patient, band, phase) and each (patient, band, half_phase) we already cache `LRGResult` containing the spectrum and the dendrogram. The eigenvectors are computed during `compute_lrg_analysis` but **not currently persisted in the NPZ cache**.

## 2. Step 0 — extend the LRG cache to persist eigenvectors

**Why this is needed.** All three E-rungs require the full top-k eigenvector basis. Computing it on-demand from the cached FC matrix is cheap (eigh on a ~120×120 matrix takes < 100 ms), but doing it 480 times per producer run and across multiple producers is wasteful. Cache once, reuse everywhere.

**One-line change in `src/lrg_eegfc/workflow/lrg.py:compute_lrg_analysis`.**
Add `eigenvectors=V` (shape `(N, N)`, float64, ~120 KB per file) to the saved NPZ. Update `LRGResult` schema (`workflow/lrg.py:61-102`) to include the field. Backwards-compatible: if loading an old NPZ without `eigenvectors`, fall back to `None` and log a warning.

**Re-cache cost.** ~480 LRG NPZ files (n=10 × 6 bands × 4 phases + n=10 × 6 bands × 4 phases × 2 halves). Each LRG computation takes ≈ 0.5 s at the current cache density. Total < 5 minutes single-threaded; CLI command `lrg-eegfc compute lrg --patients all --fc-method imcoh_abs --include-halves --refresh`.

**Verification.** After re-cache, `LRGResult.eigenvectors @ np.diag(LRGResult.spectrum) @ LRGResult.eigenvectors.T` must reconstruct L̂ to float tolerance (orthonormality + spectral theorem). Add this as a smoke check in `audit_25_eigenvector_cache_check.py` (one row per patient/band/phase to a CSV, all `pass=True` expected).

**PR sequence (step 0).**
- 0a. Edit `compute_lrg_analysis` to persist eigenvectors. *Code change.*
- 0b. Re-cache n=10. *Compute, no code change.*
- 0c. `audit_25_eigenvector_cache_check.py` smoke test. *Audit code.*

---

## 3. E1 — Spectral subspace alignment

**The cleanest E-rung. Does the diffusion mode subspace rotate across phases?**

### 3.1 Definition

For each (patient, band, phase φ ∈ {pre, test, post}), let `V^φ_k = [v_2^φ, v_3^φ, …, v_{k+1}^φ]` be the N × k matrix of top-k non-trivial eigenvectors of L̂^φ. The **subspace** spanned by `V^φ_k` is a point on the Grassmann manifold of k-dimensional subspaces of ℝ^N.

The standard Grassmann distance between two k-subspaces V_a, V_b is computed via principal angles. Numerically: form `M = V_a^T V_b` (k × k); singular values of `M` are the cosines of principal angles `σᵢ = cos(θᵢ)` for i=1…k. Then:

- **Chordal distance:** `d_chord(V_a, V_b) = sqrt( k − Σᵢ σᵢ^2 )`. Range `[0, sqrt(k)]`. Zero when subspaces coincide.
- **Geodesic distance:** `d_geo(V_a, V_b) = sqrt( Σᵢ θᵢ^2 )` where `θᵢ = arccos(min(σᵢ, 1))`. Range `[0, k·π/2]`. Zero when subspaces coincide.

We use **chordal** distance throughout (numerically stable, no arccos near 1). The asymmetry in trace direction is handled by per-patient differences:

`Δ^E1_k(p, b) = d_chord(V^test_k, V^post_k) − d_chord(V^pre_k, V^post_k)`

**Trace direction:** `Δ^E1_k < 0` means test-subspace is closer to post-subspace than pre-subspace is. Equivalently, **negative Δ = trace** (the diffusion modes during task align more with post than the pre baseline did).

### 3.2 Predicate

Per (patient, band, k) cell positivity: `Δ^E1_k(p, b) < 0 AND |Δ^E1_k(p, b)| > Δ^E1_null_drift_k(p, b)`.

### 3.3 Null

Half-baseline drift: `Δ^E1_null_drift_k = d_chord(V^rpre_B_k, V^rpost_A_k) − d_chord(V^rpre_A_k, V^rpost_A_k)` from the halves cache. Same Δ-formula but with `(rpre_A, rpre_B, rpost_A)` substituted for `(pre, test, post)`. Captures within-rest subspace drift baseline.

### 3.4 Scale axis

k ∈ {2, 3, 5, 8, 13, 21}. (Sparse log-ish grid; matches the "geometric" feel without being k=2…N. Six grid points, ridge condition becomes ridge ≥ 2 contiguous.) Eligibility: `k ≤ N − 1` for every patient (always true since N ≥ 100).

### 3.5 Pre-registered cohort decision

- Per-(p, b, k) cell: positive iff `Δ^E1_k < 0` AND `|Δ^E1_k| > Δ^E1_null_k`.
- Per-(b, k): cohort `frac_pos ≥ 0.8` AND Wilcoxon q (BH-FDR within rung, m=6 bands) ≤ 0.05.
- Per-band V verdict: positive iff `frac_pos ≥ 0.8` at any k AND ridge ≥ 2 contiguous along the k-axis.

Same structure as L5_k in the existing matrix; ridge length tightened to 2 because the k-grid has only 6 points.

### 3.6 Sanity-test gate (must pass before producer is run on real data)

Construct a stochastic block model graph with 4 blocks of 25 nodes each, intra-block weight 0.8, inter-block weight 0.1 + Gaussian noise 0.05. Compute L̂; eigenvectors v₂, v₃, v₄ are the block-membership indicators. **Verify:** chordal distance between two independent realisations of the same SBM is < 0.3 at k=4; chordal distance between SBM and shuffled-block SBM is > 0.8 at k=4. If gate fails, fix the implementation before running on patient data.

### 3.7 Implementation

**New script:** `scripts/01_compute/audit/audit_26_e1_subspace_alignment.py`. Reads `LRGResult.eigenvectors` for each (patient, band, phase, k); computes chordal distance pairs; emits

`data/audit/spectral_subspace/e1_subspace_alignment_n10_imcoh_abs.csv`

with columns: `(patient, band, k, d_pre_post, d_test_post, delta_e1, d_null_pre_post, d_null_rpre_rpost, delta_null, passes_null)`.

**Library entrypoints:** `np.linalg.svd`, `np.linalg.eigh`. No new helpers needed unless `chordal_distance(V_a, V_b)` is called from a second script later → promote to `lrg_eegfc.utils.metrics.spectral.chordal_distance` at second use.

**Cost.** Per (patient, band, k): one SVD on a k × k matrix, plus reading two NPZs. Total cohort run < 30 s.

**PR sequence (E1).**
- E1a. Scope report: `.agents/guides/task-persistence-investigation/2026-04-29_e1-spectral-subspace-alignment.md` per the README scope-report structure.
- E1b. Implementation + sanity-test gate.
- E1c. Cohort run + integration into matrix v1.2 (one new column).

---

## 4. E2 — Diffusion-map embedding drift

**Per-node embedding distance across phases. Tests local rather than global reorganisation.**

### 4.1 Definition

For each (patient, band, phase), define the Coifman–Lafon diffusion-map embedding of node x:

`Ψ^φ_k(x) = ( exp(−τλ_2^φ/2) v_2^φ(x), exp(−τλ_3^φ/2) v_3^φ(x), …, exp(−τλ_{k+1}^φ/2) v_{k+1}^φ(x) )`

This is a k-dimensional point in ℝ^k. Distances in this space approximate diffusion distances: nodes that are "close in the diffusion sense" embed close together.

**Per-node trace measure:**
`δ^E2_k(p, b, x) = ‖Ψ^test_k(x) − Ψ^post_k(x)‖_2 − ‖Ψ^pre_k(x) − Ψ^post_k(x)‖_2`

**Trace direction:** `δ^E2_k < 0` for node x means x's diffusion-position during task is closer to its post-position than pre-position is. Negative = trace.

**Cohort statistic per (p, b, k):** `n_T-nodes(p, b, k) = #{x : δ^E2_k(p, b, x) < 0 AND |δ^E2_k| > null_drift}`. Per-cell positivity threshold: `n_T-nodes ≥ τ_n(b, k)` where `τ_n(b, k)` is the 95th percentile of the null distribution.

### 4.2 Null

Half-baseline per-node drift: `δ^E2_null_k(p, b, x) = ‖Ψ^rpre_B_k(x) − Ψ^rpost_A_k(x)‖_2 − ‖Ψ^rpre_A_k(x) − Ψ^rpost_A_k(x)‖_2`. 95th percentile of `n_T-nodes_null` cohort distribution gives `τ_n(b, k)`.

### 4.3 Scale axis

Same as E1: k ∈ {2, 3, 5, 8, 13, 21}.

### 4.4 Critical caveat — eigenvector sign and ordering ambiguity

Eigenvectors of L̂_pre and L̂_post are computed independently; their *signs* and *orderings* are not canonical. If we compute `Ψ^pre_k(x) − Ψ^post_k(x)` naively, an eigenvector that flips sign between phases gives a meaningless 2× scale-up of the difference.

**Mitigation:** before computing Ψ, **align the eigenvectors of L̂^post to those of L̂^pre** via signed best-match (per-eigenvector: pick the sign of `vᵢ^post` that maximises `<vᵢ^post, vⱼ^pre>` over j ∈ {1, …, k}; if best-match j ≠ i, swap eigenvectors). This is the "eigenvector continuation" / Procrustes-on-spectral-basis step. Document carefully in scope.

Without this alignment E2 cannot give meaningful per-node distances; with it, E2 measures "how much each node's diffusion position drifts in the aligned spectral basis".

### 4.5 Sanity-test gate

Same SBM as E1. **Verify:** within-realisation per-node embedding-distance < 0.1 (median); shuffled-block embedding-distance > 0.5 (median). Eigenvector sign-and-order alignment must reduce diff to within-realisation noise.

### 4.6 Implementation

**New script:** `scripts/01_compute/audit/audit_27_e2_diffusion_map_drift.py`. Output:

`data/audit/spectral_subspace/e2_diffusion_map_drift_n10_imcoh_abs.csv` columns `(patient, band, k, x, delta_e2, delta_null, passes_null)`. Emits per-(patient, band, k) summary `n_T-nodes` to `e2_summary_n10_imcoh_abs.csv`.

**Library promotion:** `lrg_eegfc.utils.metrics.spectral.align_eigenvectors(V_target, V_source, k)` is a new helper used by E2 and E3; promote on first use here. Keep `chordal_distance` and `align_eigenvectors` together in the same module.

**Cost.** Per (patient, band, k): one alignment (k SVDs), N node-distances. Total cohort run < 1 min.

**PR sequence (E2).** Same shape as E1: scope, implementation+sanity, cohort run+matrix integration.

---

## 5. E3 — Inverse-participation-ratio shifts

**Eigenvector localisation. Are diffusion modes more concentrated on specific node clusters during/after task?**

### 5.1 Definition

For eigenvector `vᵢ^φ` (length N, normalized to ‖vᵢ^φ‖₂ = 1), the **inverse participation ratio** is:

`IPR(vᵢ^φ) = Σ_x vᵢ^φ(x)^4 ∈ [1/N, 1]`

`IPR = 1/N` when the eigenvector is fully delocalised (uniform on all nodes); `IPR = 1` when fully localised on one node. The IPR distribution `{IPR(v_2^φ), IPR(v_3^φ), …, IPR(v_{k+1}^φ)}` characterises how "modular" the diffusion basis is.

**Per-(p, b) trace statistic:** Wasserstein-1 distance between IPR distributions across phases.

`W^E3(p, b, φ_a, φ_b) = W₁({IPR(vᵢ^{φ_a})}_{i=2..k+1}, {IPR(vᵢ^{φ_b})}_{i=2..k+1})`

**Trace direction:** `W^E3(test, post) < W^E3(pre, post)` for a band means test localisation distribution is closer to post than pre is.

`Δ^E3(p, b) = W^E3(test, post) − W^E3(pre, post)`

Negative Δ = trace.

### 5.2 Null

Halves drift: `Δ^E3_null(p, b) = W^E3(rpre_B, rpost_A) − W^E3(rpre_A, rpost_A)`. Per-cell: `Δ^E3 < 0 AND |Δ^E3| > Δ^E3_null`.

### 5.3 Scale axis

None — E3 has no scale axis. The "k" of the IPR distribution is fixed at k=20 (top 20 non-trivial eigenvectors). Cohort verdict per band requires only `frac_pos ≥ 0.8` and Wilcoxon q ≤ 0.05 (no ridge condition).

### 5.4 Eigenvector ordering caveat

IPR is sign-invariant (fourth power) so signs don't matter. BUT eigenvector *ordering* across phases matters — comparing `IPR(v_5^pre)` to `IPR(v_5^post)` is meaningless if they correspond to different physical modes. Wasserstein on the *distribution* sidesteps this (compares the sorted IPR values not paired ones). That's the correct framing for E3.

### 5.5 Sanity-test gate

SBM example: single block with weak inter-block coupling has eigenvectors with IPR ≈ 1/(block size) (delocalised within block); two-block SBM has more localised modes. **Verify:** IPR distribution Wasserstein distance between SBM and shuffled-block SBM > 0.05; within-realisation distance < 0.01.

### 5.6 Implementation

**New script:** `scripts/01_compute/audit/audit_28_e3_ipr_shifts.py`. Output:

`data/audit/spectral_subspace/e3_ipr_shifts_n10_imcoh_abs.csv` with `(patient, band, w_pre_post, w_test_post, delta_e3, delta_null, passes_null)`.

**Library:** `scipy.stats.wasserstein_distance` directly. No new helpers.

**Cost.** Per (patient, band): one Wasserstein computation on length-20 distributions. Total cohort run < 5 s.

**PR sequence (E3).** Same shape as E1, E2.

---

## 6. Matrix v1.2 integration

After steps 0, E1, E2, E3 land, extend the cohort-coverage matrix:

- New rungs: `E1_subspace`, `E2_diffmap`, `E3_ipr`. RUNG_ORDER becomes 11.
- New consumer entries in `audit_24_cohort_coverage_matrix.py` CONSUMERS dict.
- Triangulation predicate update: add `eigenvector-confirmed` label when `V_E1 = positive AND V_E5_k = positive` (subspace AND modular agreement).
- The retired headline-triangulated predicate in v1 (which required L1 + L5_k + (L3 OR L7) and is structurally unreachable in v1) gets revisited: if E1 is positive cohort-wide, define `headline-triangulated-spectral` requiring V_L5_k AND V_E1 AND (V_E2 OR V_E3) positive.

This is **matrix v1.2**, not v2 — v2 was the dendrogram-based ladder with §7 controls landed; that path is paused per the post-mortem decision.

## 7. Pre-registered acceptance gate

Before running E1 on real data, the user commits to the following decisions for the v1.2 outcome:

**Outcome A (E1 surfaces α / β / γ_h cohort-positive that current matrix declared ergodic).** Eigenvector-direct hypothesis is supported. Paper headline pivots to "diffusion-mode reorganisation as the load-bearing trace measure"; dendrogram-based ladder relegated to supplementary; E2, E3 land as triangulation. δ stays as partition-resolution-locked. New scientific story is more positive.

**Outcome B (E1 confirms ergodic for α / β / γ_h).** Negative result is much stronger now (it covers both tree-based and eigenvector-direct probes). Paper headline is "we tested both tree-based and eigenvector-direct LRG probes; only δ partition-resolution-locked survives controls." E2 and E3 still run for completeness but become diagnostic figures, not headline-bearing.

**Outcome C (E1 mixed — some bands pass, some fail).** Paper reports the rung-by-rung and band-by-band gradient honestly, with E1 as the headline graded-evidence figure and the binary triangulation downgraded to supplementary.

In all three outcomes, the post-mortem's K1 (δ partition-resolution-locked is real), K2 (strict-J cohort-null), K4 (spectrum-scalar dead) stay. The pivot does not invalidate prior results — it adds a complementary test of the same claim using a different geometric reduction.

## 8. PR sequence (full)

| # | Step | Type | Owner |
|---|---|---|---|
| 0a | Persist eigenvectors in `LRGResult` cache | code change | one PR |
| 0b | Re-cache n=10 LRG (full + halves) | compute | one PR |
| 0c | `audit_25_eigenvector_cache_check.py` smoke test | audit code | one PR |
| 1a | E1 scope report `2026-04-29_e1-spectral-subspace-alignment.md` | scope | one PR |
| 1b | `audit_26_e1_subspace_alignment.py` + SBM sanity gate | audit code | one PR |
| 1c | E1 cohort run + matrix v1.2 integration | audit + matrix | one PR |
| 2a | E2 scope report `2026-04-29_e2-diffusion-map-drift.md` | scope | one PR |
| 2b | `audit_27_e2_diffusion_map_drift.py` + SBM sanity + `align_eigenvectors` library helper | audit code | one PR |
| 2c | E2 cohort run + matrix v1.2 integration | audit + matrix | one PR |
| 3a | E3 scope report `2026-04-29_e3-ipr-shifts.md` | scope | one PR |
| 3b | `audit_28_e3_ipr_shifts.py` + SBM sanity | audit code | one PR |
| 3c | E3 cohort run + matrix v1.2 integration | audit + matrix | one PR |
| 4 | Graded companion verdict added to matrix figure (alongside binary V) | figure code | one PR |
| 5 | v1.2 narrative report `2026-04-29_eigenvector-direct-verdict.md` | scope | one PR |

**Critical path:** 0a → 0b → 0c → 1a → 1b → 1c. After 1c, decide based on outcome A / B / C whether to continue with E2, E3 or stop. Steps 2a-c, 3a-c can run in parallel after 1c if the outcome warrants.

## 9. Verification gates

**Gate-0** (cache): `LRGResult.eigenvectors @ diag(spectrum) @ eigenvectors.T ≈ L` to float tolerance for every cached file.

**Gate-1** (E1): SBM sanity test passes (within-realisation chordal < 0.3, shuffled chordal > 0.8 at k=4). Cohort run reproduces bit-identical on re-run with no RNG.

**Gate-2** (E2): SBM sanity test passes; `align_eigenvectors` reduces within-realisation embedding distance to within-realisation noise level. Wasserstein on per-node distribution is rotation-invariant.

**Gate-3** (E3): SBM sanity test passes; `wasserstein_distance` matches scipy reference; IPR sums to 1 per eigenvector.

**Gate-Cohort** (every E-rung): Phase 0 audit verdict on the new measure is `current` before integration into matrix v1.2.

**Gate-Triangulation** (matrix v1.2): pre-registered triangulation predicates do not change post-hoc based on outcome. The `eigenvector-confirmed` label is committed to the decision-rules scope BEFORE the cohort run.

## 10. Risks and known weak points

- **Eigenvector sign / ordering ambiguity (E2 critical).** The `align_eigenvectors` helper is the most error-prone piece of the plan. SBM sanity gate must catch alignment bugs before cohort run.
- **k-grid arbitrariness.** k ∈ {2, 3, 5, 8, 13, 21} is hand-picked. A sensitivity sweep over k is owed before the matrix v1.2 verdict is published.
- **τ-sensitivity.** All three E-rungs operate at τ = 1/λ_max. Whether E1/E2/E3 verdicts survive at coarser τ (e.g. 1/λ_gap) is an open question, same as for the existing matrix.
- **SBM is artificial.** The sanity gates use stochastic block model graphs, which have discrete topology. They confirm the implementations are correct but not that the measures are informative on continuous-spectrum graphs. The cohort run is the actual test.
- **n=10 is small for cohort thresholds.** Same caveat as before.
- **Cohort-frame vs implant-heterogeneity.** Open question O3 from the post-mortem applies equally here: a cohort-positive E1 verdict is still a topology-only claim.

## 11. What this plan does not address

- The biological-interpretation question (sham-task null) — separate scope, parked.
- Per-patient case-study framing — separate paper-framing decision, parked.
- The KC λ-sweep (existing §7.2 scope) — paused indefinitely; KC is tree-based and may inherit the dendrogram-compression issue.
- L4 / L5_hrel / L6 / L7 within-baseline nulls (existing §7.3, §7.4, §7.5, §7.6 scopes) — paused; the dendrogram-based ladder will not be extended further.

## 12. Pointer to next conversation

After compaction, resume with: "Step 0a — extend `LRGResult` to persist eigenvectors per the eigenvector-direct pivot plan."

Reading order on resume:
1. `.agents/reports/2026-04-29_critical-post-mortem.md` (the why).
2. This document (the what).
3. `.agents/guides/02_methods/lrg-framework-guide.md` §6 (the constraint that motivates the pivot).
4. `.agents/reports/2026-04-28_residual-subspace-diagnostic.md` (the dead probe to NOT repeat).
