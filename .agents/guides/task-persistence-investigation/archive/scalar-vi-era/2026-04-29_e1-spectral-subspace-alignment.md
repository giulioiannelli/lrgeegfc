---
name: 2026-04-29_e1-spectral-subspace-alignment
type: scope
era: COHORT_N10
status: draft
created: 2026-04-29
updated: 2026-04-29
pointers:
  - .agents/plans/active/2026-04-29_eigenvector-direct-pivot-plan.md
  - .agents/reports/2026-04-29_critical-post-mortem.md
  - .agents/guides/02_methods/lrg-framework-guide.md
  - .agents/guides/task-persistence-investigation/2026-04-29_decision-rules.md
  - .agents/reports/2026-04-28_residual-subspace-diagnostic.md
---

# E1 — Spectral subspace alignment

## Renormalization head

**E1 measures whether the diffusion-mode subspace of the Laplacian rotates between phases by computing the Grassmann chordal distance between the top-k non-trivial eigenspaces of `L̂^pre`, `L̂^test`, `L̂^post`. Trace direction is the per-patient difference `Δ^E1_k = d_chord(V^test_k, V^post_k) − d_chord(V^pre_k, V^post_k)`; negative means task-subspace aligned with post-subspace more than pre did. The null is the within-rest analogue computed from the halves cache (`rpre_A, rpre_B, rpost_A`). E1 sees what dendrogram-based rungs cannot see: a rotation of the diffusion basis that may not show up as either pair-distance shift (L1) or modular swap (L5_k).**

---

## TOC

1. **The probe is the angle between two k-dimensional eigenspaces of the Laplacian, not a scalar function of eigenvalues.** Subspace, not spectrum.
2. **Top-k eigenvectors v_2…v_{k+1} are non-trivial; v_1 (constant mode) is dropped.**
3. **Chordal distance is `sqrt(k − Σ σᵢ²)` where σᵢ are singular values of `V_a^T V_b`.** Standard Grassmann manifold metric.
4. **The trace statistic Δ^E1_k is sign-aligned with the dendrogram rungs: negative = trace.**
5. **Within-baseline null is computed from the halves cache, not a permutation null.**
6. **Six grid points `k ∈ {2, 3, 5, 8, 13, 21}`; ridge condition tightened to ≥ 2.**
7. **SBM (stochastic block model) sanity gate must pass before patient run.**

---

## 1. Notation

- *L̂^φ* = Laplacian D̂^φ − Â^φ for phase `φ ∈ {pre, test, post}`. Symmetric positive-semidefinite, fully connected, N nodes per patient (100 ≤ N ≤ 122 across cohort, see `data/audit/cohort_metadata.csv`).
- *λᵢ^φ* = i-th eigenvalue of L̂^φ, ascending: 0 = λ_1^φ < λ_2^φ ≤ … ≤ λ_N^φ.
- *vᵢ^φ* = i-th unit-norm eigenvector of L̂^φ (column of the spectral basis V^φ ∈ ℝ^{N × N}).
- *V^φ_k* = N × k matrix `[v_2^φ, v_3^φ, …, v_{k+1}^φ]`. The k non-trivial eigenvectors at phase φ.
- *𝒮(V^φ_k)* = the k-dimensional linear subspace of ℝ^N spanned by the columns of V^φ_k. A point on the Grassmann manifold `Gr(k, N)`.
- *τ* = diffusion time, fixed at τ = 1/λ_max throughout (per `lrg-framework-guide.md` §3 + framework decision 2026-04-29). E1 does not directly use τ but the eigenvectors are the same V used in the heat kernel `K̂(τ) = V diag(exp(−τλ)) V^T`.
- *halves cache* = `data/cache/imcoh_lrg_halves/<pat>/<band>_<phase>_<A|B>_lrg_imcoh-abs.npz` containing one `LRGResult` per (patient, band, phase ∈ {rest_pre, rest_post}, half ∈ {A, B}). Eigenvectors persisted 2026-04-29 (Phase 0, step 0a of pivot plan).
- *Cohort* = n=10, IMCOH_ABS, `Pat_{02, 03, 05, 06, 07, 08, 10, 13, 14, 15}`.

## 2. Definitions

### 2.1 Principal angles and chordal distance

Given two k-dimensional subspaces with orthonormal bases `V_a, V_b ∈ ℝ^{N × k}` (columns are unit-norm and pairwise orthogonal — the eigenvectors of a symmetric matrix automatically satisfy this), form the cross-Gram matrix `M = V_a^T V_b ∈ ℝ^{k × k}`. The singular values `σ_1 ≥ σ_2 ≥ … ≥ σ_k ≥ 0` of M are the cosines of the principal angles `θ_i = arccos(σ_i) ∈ [0, π/2]` between the two subspaces (Björck & Golub 1973; standard linear-algebra reference).

**Chordal distance:**

```
d_chord(V_a, V_b) = sqrt( k − Σ_{i=1}^{k} σ_i² )    ∈ [0, sqrt(k)]
```

`d_chord = 0` iff the subspaces coincide (all σ_i = 1, all θ_i = 0). `d_chord = sqrt(k)` iff the subspaces are orthogonal (all σ_i = 0).

**Why chordal not geodesic.** Both are valid metrics on `Gr(k, N)`. Chordal is numerically stable (no `arccos` near 1 where the derivative blows up); geodesic `sqrt(Σ θ_i²)` requires `arccos(min(σ, 1))` clipping and is harder to bootstrap. Chordal is the standard choice in subspace-comparison literature on noisy spectra (Edelman, Arias & Smith 1998).

### 2.2 The trace statistic Δ^E1_k(p, b)

Per (patient p, band b, scale k):

```
Δ^E1_k(p, b) = d_chord(V^test_k, V^post_k) − d_chord(V^pre_k, V^post_k)
```

**Sign convention:** Δ^E1_k < 0 means test-subspace is closer to post-subspace than pre-subspace was → task-induced rotation persists into rest_post → trace. Aligns with the dendrogram-based rungs L1/L5_k where `negative = trace`.

### 2.3 The within-baseline null Δ^E1_null,k(p, b)

Substitute the halves triple `(rpre_A, rpre_B, rpost_A)` for the cross-phase triple `(pre, test, post)`:

```
Δ^E1_null,k(p, b) = d_chord(V^rpre_B_k, V^rpost_A_k) − d_chord(V^rpre_A_k, V^rpost_A_k)
```

This captures the within-rest split-half subspace drift baseline. If `|Δ^E1_k|` does not exceed `|Δ^E1_null,k|`, the apparent task-induced rotation is no larger than within-rest noise.

### 2.4 Per-cell positivity predicate

```
positive(p, b, k) ⇔ ( Δ^E1_k(p, b) < 0 )  AND  ( |Δ^E1_k(p, b)| > |Δ^E1_null,k(p, b)| )
```

### 2.5 Per-band cohort verdict

Mirrors L5_k structure:

```
frac_pos(b, k)  = (1 / |Cohort|) Σ_{p ∈ Cohort} 1{positive(p, b, k)}
```

The Wilcoxon signed-rank test on `(Δ^E1_k − Δ^E1_null,k)` across patients gives a p-value per (b, k); BH-FDR within rung over m = 6 bands gives q.

```
V_E1(b) = positive   if  ∃ k* such that  frac_pos(b, k*) ≥ 0.8
                                       AND q(b, k*) ≤ 0.05
                                       AND ∃ contiguous k-window of length ≥ 2
                                           where frac_pos ≥ 0.8 AND q ≤ 0.05
        = silent     if  no k passes the cohort threshold AND null is present
        = ineligible if  >50% of cells missing (eigenvector cache absent for too many)
```

Ridge condition tightened from L5_k's ≥ 3 to ≥ 2 because the k-grid has only 6 points (`{2, 3, 5, 8, 13, 21}`); a ridge of 3 contiguous points covers 50% of the axis, which is too permissive. Length 2 covers 33% — calibrated to the same proportion as L5_k's "≥ 3 contiguous in 28-step k=2…29 axis ≈ 11%" but recognising the coarser grid.

## 3. Properties

- **Range.** `d_chord ∈ [0, sqrt(k)]`. So `Δ^E1_k ∈ [−sqrt(k), sqrt(k)]`.
- **Invariance.** Chordal distance depends only on the *subspace*, not on the basis. Eigenvector sign flips and intra-block rotations of the eigenbasis (when there are repeated eigenvalues) do not affect d_chord. **This is why E1 is the cleanest E-rung — no sign-alignment helper needed.**
- **Compute cost.** Per (p, b, k): one SVD on a k × k matrix. With k ≤ 21 and 10 × 6 = 60 patient × band cells × 6 k-values × 3 phase pairs = 1080 SVDs, total runtime under 5 s.
- **Trivial-mode handling.** v_1 is the constant mode (eigenvalue 0, fully delocalised). Including it makes d_chord(V^φ_1, V^ψ_1) = 0 always (constant vector matches constant vector), shifting d_chord uniformly. Excluding it focuses the test on the diffusion modes that actually carry phase information.
- **What E1 cannot detect.** Eigenvalue *magnitude* changes that don't rotate the subspace (homothety of the spectrum); changes in eigenvector localisation that preserve the subspace span (rotation within the subspace itself, e.g. v_2 ↔ v_3 swap inside the (v_2, v_3) plane). These would show up in E2 / E3 but not E1.
- **What E1 detects that L1/L5_k miss.** If the diffusion modes rotate (different functional sub-networks become "principal") but pair-distances on D average out (no ultrametric shift cohort-wide) and partition labels at integer-k stay the same (no modular swap), E1 still flags it. The hypothesis is that this regime exists for α / β / γ_h.

## 4. Caveats and failure modes

- **Eigenvalue degeneracy.** Repeated eigenvalues mean the corresponding eigenvectors are only defined up to a rotation within the eigenvalue's eigenspace. For our continuous-spectrum FC graphs with weight heterogeneity, exact degeneracies are measure-zero. Numerically, near-degenerate eigenvalues at the boundary k vs k+1 could include / exclude an eigenvector by chance. **Mitigation:** the sparse k-grid `{2, 3, 5, 8, 13, 21}` reduces the chance that a near-degeneracy lands exactly at the cutoff; a sensitivity check using k → k+1 should not flip the verdict.
- **Anatomy-dominated eigenvectors at low k.** The first non-trivial mode v_2 (Fiedler vector) is often dominated by anatomy/probe-bias for weighted FC graphs. Trace at k=2 may be anatomy alignment, not task. **Mitigation:** report both `frac_pos(b, k)` curves so the reader can see whether trace is k-localised at low k (anatomy-suspect) vs distributed (cleaner subspace effect).
- **Half-baseline dilution.** The halves are computed with halved Welch nperseg; the resulting eigenvectors are noisier than full-phase. This *inflates* the null and makes E1 conservative — Δ^E1_null,k > 0 noise floor. False-positive rate is suppressed but false-negative rate may be elevated. **Mitigation:** report cohort `|Δ^E1_null,k|` distribution alongside `|Δ^E1_k|` so the reader can see the noise floor.
- **N varies across patients.** Pat_03 has N ≈ 100, Pat_07 has N = 122. Chordal distance is normalised by k (range `[0, sqrt(k)]` independent of N) so the per-cell statistic is comparable across patients. But k=21 is a larger fraction of N for Pat_03 than Pat_07. **Mitigation:** report cohort distribution of `Δ^E1` per k; if Pat_03 systematically separates, flag as anatomy-suspect.
- **Trivial mode v_1.** L̂ has λ_1 = 0 with eigenvector v_1 ∝ 1 (constant). Always exclude it from V^φ_k. The persisted `eigenvalues, eigenvectors` arrays are sorted ascending, so `V^φ_k = eigenvectors[:, 1:k+1]`. **Implementation gate:** verify `|λ_1| < 1e-8` for every cell (smoke check in audit_26).
- **k > N − 1.** Edge case. Always satisfied for our cohort since N ≥ 100 > 21. Producer should still assert.

## 5. Pseudocode

```
function chordal_distance(V_a, V_b):
    # V_a, V_b are N × k orthonormal bases
    M = V_a.T @ V_b                   # k × k cross-Gram
    sigma = singular_values(M)        # length-k array, σ_i ∈ [0, 1]
    return sqrt( k − sum(sigma**2) )

function delta_e1_cell(pat, band, k, eigvecs_pre, eigvecs_test, eigvecs_post):
    # Drop trivial mode (column 0); take columns 1..k.
    Vpre  = eigvecs_pre[:,  1 : k+1]
    Vtest = eigvecs_test[:, 1 : k+1]
    Vpost = eigvecs_post[:, 1 : k+1]
    d_pre_post  = chordal_distance(Vpre,  Vpost)
    d_test_post = chordal_distance(Vtest, Vpost)
    return d_test_post − d_pre_post

function delta_e1_null_cell(pat, band, k, eigvecs_rpreA, eigvecs_rpreB, eigvecs_rpostA):
    VA = eigvecs_rpreA[:,  1 : k+1]
    VB = eigvecs_rpreB[:,  1 : k+1]
    VC = eigvecs_rpostA[:, 1 : k+1]
    d_A_C = chordal_distance(VA, VC)
    d_B_C = chordal_distance(VB, VC)
    return d_B_C − d_A_C

for each patient p in Cohort:
    for each band b in BANDS:
        load LRG eigvecs for (p, b, pre|test|post) from full cache
        load LRG eigvecs for (p, b, rpre_A|rpre_B|rpost_A) from halves cache
        for each k in K_GRID:
            cells.append(delta_e1_cell(...), delta_e1_null_cell(...))

aggregate:
    per (b, k): frac_pos = mean(passes_null AND delta < 0)
                wilcoxon p-value on (delta - delta_null)
                BH-FDR q
    per b:      verdict per §2.5
```

## 6. Visualisation spec

Two-panel figure per band (faceted 6-band grid in a single PDF):

- **Panel A — frac_pos(k) curve.** x = k, y = cohort fraction positive ∈ [0, 1]. Solid line = `frac_pos(b, k)` for E1; dashed reference = 0.8 cohort threshold; markers `×`/`+` at k where Pat_02 / Pat_03 contribute to `frac_pos`. Title states the v_E1 verdict and ridge length.
- **Panel B — per-patient |Δ^E1_k| vs |Δ^E1_null,k| at the verdict-k.** Scatter, one point per patient, x = |null|, y = |signal|. Diagonal y=x = "no separation"; points above diagonal = signal exceeds null. Shape: square = positive cell, open circle = silent. Pat_02 marked with `×`, Pat_03 with `+`.

Reader rules: `frac_pos ≥ 0.8` line crossing at any k is necessary; ridge ≥ 2 contiguous = sufficient. Scatter above diagonal at the verdict-k confirms cohort-wide separation, not just one outlier.

No `fig.suptitle`. PDF only, full vector (do not call `set_rasterized`).

## 7. Connection to prior tools

| Prior tool | Rung | Substrate | What E1 adds |
|---|---|---|---|
| L1 (H2c continuous trace, Spearman ρ on D upper triangle) | L1 | pair distances on D | E1 sees subspace rotation that pair-averages cancel |
| L5_k (Δ_VI(k) on integer-k partitions) | L5_k | partition labels | E1 sees rotation that doesn't change partition labels (label-invariant subspace rotation) |
| L5_hrel (Δ_VI at fixed h_rel) | L5_hrel | partition labels at shared resolution | same as L5_k caveat |
| Residual-subspace `α_1` (closed 2026-04-28) | (dead) | rank-1 leading eigenvector of the *residual* `D^test − D^pre` | E1 uses rank-k eigenvectors of the *Laplacian itself*, not residuals — different math |

E1 complements the dendrogram-based rungs by probing a different geometric reduction of the same primitive `L̂` that builds them. **It does not replace L5_k**: a band that is partition-resolution-locked at L5_k (δ) is allowed to also be subspace-positive at E1 (consistent), or subspace-silent (the partition swap was modular but not basis-rotating). The matrix triangulation rules (`.agents/guides/task-persistence-investigation/2026-04-29_decision-rules.md`) define the cross-rung verdicts.

E1 explicitly does not subsume the failed residual-subspace probe (`.agents/reports/2026-04-28_residual-subspace-diagnostic.md`). That probe was killed by near-rank-1 dominance of the residual matrix — a property of `D^test − D^pre`, not of `L̂`. E1 operates on the Laplacian's eigenvectors, which are not rank-1-dominant for our continuous-spectrum FC graphs.

## 8. Implementation plan

**Producer:** `scripts/01_compute/audit/audit_27_e1_subspace_alignment.py` (audit_26 is the cache-check; audit_27 is the first E-series script).

**Output:** `data/audit/spectral_subspace/e1_subspace_alignment_n10_imcoh_abs.csv` with columns:

```
patient, band, k,
d_pre_post, d_test_post, delta_e1,
d_null_rpreA_rpostA, d_null_rpreB_rpostA, delta_e1_null,
passes_null, sign_negative, positive
```

**Library use only — no new helpers in this PR:**
- `lrg_eegfc.workflow.lrg.load_lrg_result` for both full-phase and halves caches (cache_root override for halves).
- `np.linalg.svd` for principal angles.
- `lrg_eegfc.utils.metrics.hypothesis.wilcoxon_z, bh_fdr` for cohort stats.
- `lrg_eegfc.config.const.PATIENTS_LIST, BRAIN_BANDS_NAMES`.

**Promotion to library deferred:** if a second E-rung (E2 / E3) reuses `chordal_distance` or `principal_angles`, promote on second use to `lrg_eegfc.utils.metrics.spectral`. E2 will use a different distance (Euclidean on diffusion-map embedding), so promotion is unlikely to be needed at E1 alone.

**SBM sanity gate (BLOCKING).** Before the cohort run, the producer must pass:
- 4-block SBM, 25 nodes per block, intra-block weight 0.8, inter-block weight 0.1 + Gaussian noise σ=0.05.
- Compute L̂; the natural test scale is `k = n_blocks − 1 = 3` because v_1 is the trivial constant mode (always shared) and v_2, v_3, v_4 span the block-indicator subspace. Going to k=4 crosses the spectral gap into intra-block noise modes where two independent realisations genuinely differ — a correct chordal_distance WILL flag that, so we calibrate at k=3.
- **Pass criteria** (all three required):
  - chordal distance between two independent SBM realisations at k=3 < 0.3
  - chordal distance between SBM and shuffled-block SBM at k=3 > 0.8
  - reflexive: d_chord(V, V) < 1e-5  (sqrt of float64 ε at small k)

If gate fails, fix the implementation before patient run. CLI flag `--sanity-only` runs just the gate.

**CLI:**
```
python scripts/01_compute/audit/audit_27_e1_subspace_alignment.py [--sanity-only] [-v]
```

**Cohort run prerequisites:**
- Phase 0 step 0a (LRGResult eigvecs persisted) — DONE 2026-04-29.
- Phase 0 step 0b (full-phase + halves re-cache) — running 2026-04-29.
- Phase 0 step 0c (audit_26 spectral identity check) — required before audit_27 runs.

**Cost.** Cohort run < 30 s.

## 9. Decision rules and matrix integration

E1 enters the cohort-coverage matrix as a new rung `E1_subspace` consumed by `audit_24_cohort_coverage_matrix.py`. The CONSUMERS dict gains:

```python
"E1_subspace": {
    "csv": "data/audit/spectral_subspace/e1_subspace_alignment_n10_imcoh_abs.csv",
    "scale_axis": "k",
    "scale_grid": [2, 3, 5, 8, 13, 21],
    "value_col": "delta_e1",
    "null_col": "delta_e1_null",
    "ridge_min": 2,
}
```

Triangulation predicate update (decision-rules scope):

```
T_eigenvector_confirmed(b) ⇔ V_E1(b) = positive  AND  V_L5_k(b) = positive
```

Pre-registered before the cohort run lands per Gate-Triangulation in pivot plan §9.

## 10. Open questions

- **Should we also report d_geodesic alongside d_chordal?** Geodesic emphasises near-orthogonal subspaces (heavier penalty on θ_i near π/2). Chordal is robust; geodesic might surface band-specific behaviour at large rotations. Deferred — chordal is the registered metric for E1; geodesic is supplementary at most.
- **Sensitivity to k-grid choice.** `{2, 3, 5, 8, 13, 21}` is hand-picked Fibonacci-ish. A finer grid might reveal a narrow ridge at intermediate k. Deferred to post-v1.2 sensitivity analysis.
- **Should the null also include `(rpost_A, rpost_B, rpre_A)` symmetric variant?** E1 currently uses one halves triple. A symmetric null averaging both directions would be more conservative. Deferred.
- **τ-sensitivity.** All E-rungs operate at τ = 1/λ_max (the codebase default). Whether E1 verdicts survive at coarser τ (e.g. τ such that a fixed h_rel is achieved) is an open question, parallel to L5_k vs L5_hrel. Deferred.
- **What if E1 is silent for α / β / γ_h?** Then the eigenvector-direct hypothesis is not supported by the cleanest probe. E2 / E3 still run for completeness, but the paper headline pivots to the negative result per pivot plan §7 Outcome B.
