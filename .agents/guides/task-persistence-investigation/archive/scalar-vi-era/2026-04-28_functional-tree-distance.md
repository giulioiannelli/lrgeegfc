---
name: functional-tree-distance
type: scope
era: COHORT_N10
status: superseded
created: 2026-04-28
updated: 2026-04-28
superseded_by: .agents/reports/2026-04-28_functional-tree-distance-verdict.md
pointers:
  - .agents/reports/2026-04-28_functional-tree-distance-verdict.md
  - .agents/guides/02_methods/lrg-framework-guide.md
  - .agents/guides/task-persistence-investigation/2026-04-26_continuous-trace-matrix.md
  - .agents/guides/task-persistence-investigation/2026-04-25_task-trace-canonical.md
  - lrgsglib/src/lrgsglib/utils/lrg/spectral.py
  - lrgsglib/src/lrgsglib/utils/lrg/infocomm.py
  - src/lrg_eegfc/workflow/lrg.py
---

> **Status: SUPERSEDED 2026-04-28.** Implemented as scoped, run on the
> n=10 cohort with within-baseline split-half null; **no cohort signal
> survived the noise gate** (3 trace cells at strict rank gate, 10 at
> permissive height gate; best band 3/10 — far below the ≥ 8/10
> cohort threshold). Verdict and full numbers in
> [`2026-04-28_functional-tree-distance-verdict.md`](../../reports/2026-04-28_functional-tree-distance-verdict.md).
> Keep this scope as the definitional record of what was tried; do
> not re-implement variants without first reading the verdict.

# Functional tree distance — comparing two LRG flows as elements of `L²(I, ℝ^{N×N})`

**Each LRG output `(p, b, φ)` is treated as a τ-indexed family of matrices —
either the kernel density `ρ(τ)` or the communication distance
`D(τ) = 1/ρ(τ)`. Phase-pair similarity is the per-τ matrix-correlation
curve `δ(τ; A, B)`, evaluated as both a snapshot at the finest meaningful
scale `τ = 1/λ_max` and a `d(log τ)`-integral over the meaningful range
`I = [1/λ_max, τ\*]` where `τ\* = argmax_τ C(τ)` is the LRG characteristic
scale (a single peak in our fully-connected case). Two parallel views
run in lock-step: rank-based `δ_S` on `D = 1/ρ` (dendrogram-topology
agreement, robust to `1/ρ` divergences) and value-based `δ_P` on `ρ`
itself (height/scale agreement, well-behaved on the bounded kernel).
Together they distinguish phase-pairs that share *clustering ordering*,
share *absolute scale*, or share both.**

---

## 1. Notation

Indices and sets:

| Symbol | Domain | Meaning |
|:---|:---|:---|
| `p ∈ P` | primitive set `{Pat_02, Pat_03, Pat_06, Pat_13}` (`|P| = 4`); cohort `|P| = 10` after primitive clears | patient |
| `b ∈ B` | `{delta, theta, alpha, beta, gamma_low}` etc. (`|B| = 5–6` per `BRAIN_BANDS`) | frequency band |
| `φ ∈ Φ` | `Φ = (rest_pre, task_learn, task_test, rest_post)` per `config.const.PHASE_LABELS` | recording phase |
| `(A, B) ∈ Φ × Φ`, `A ≠ B` | 6 unordered pairs per `ALL_PHASE_PAIRS` | phase-pair |

Per `(p, b, φ)`:

| Symbol | Domain | Meaning |
|:---|:---|:---|
| `Â^{p,b,φ}` | `ℝ_{≥0}^{N×N}` | adjacency from `load_fc_matrix(p, φ, b, "imcoh_abs")` |
| `L^{p,b,φ} = D̂ − Â^{p,b,φ}` | symmetric PSD | fluid Laplacian |
| `λ_i^{p,b,φ}`, `i = 0,…,N−1` | `0 = λ_0 ≤ λ_1 ≤ ⋯ ≤ λ_{N−1}` | Laplacian eigenvalues |
| `λ_max^{p,b,φ} := λ_{N−1}^{p,b,φ}` | | largest eigenvalue |
| `τ ∈ ℝ_{>0}` | | diffusion time |
| `K^{p,b,φ}(τ) := exp(−τ L^{p,b,φ})` | element-wise positive | network propagator |
| `Z^{p,b,φ}(τ) := Tr K^{p,b,φ}(τ)` | | partition function |
| `ρ^{p,b,φ}(τ) := K^{p,b,φ}(τ) / Z^{p,b,φ}(τ)` | trace 1, PSD | density operator (similarity, NOT a metric) |
| `D^{p,b,φ}(τ)` | `N×N`, ultrametric (Villegas 2025 Eq. 1) | communication distance, `D_{ij} = (1−δ_{ij})/ρ_{ij}` symmetrised, zero-diagonal |
| `S^{p,b,φ}(τ)` | `[0, 1]` | normalised von Neumann entropy |
| `C^{p,b,φ}(τ) := −dS/d(log τ)` | | spectral complexity |
| `τ\*^{p,b,φ} := argmax_τ C^{p,b,φ}(τ)` | scalar | LRG characteristic scale (single peak in the fully-connected case) |

Triangular vectorisation: `triu(M) ∈ ℝ^{N(N−1)/2}` is the upper-triangle
(`i < j`) flattened to a vector, excluding the diagonal.

For phase-pair `(A, B)` over `(p, b)`:

```
τ_min^{p,b}(A,B)  := max( 1/λ_max^{p,b,A}, 1/λ_max^{p,b,B} )      (finest τ at which both phases are at-or-coarser than their fastest mode)
τ_max^{p,b}(A,B)  := min( τ\*^{p,b,A}, τ\*^{p,b,B} )                (coarsest τ at which both phases are still at-or-finer than their characteristic scale)
I^{p,b}(A,B)      := [τ_min^{p,b}(A,B), τ_max^{p,b}(A,B)]          (meaningful interval)
```

Validity: `I` is non-empty iff `τ_min < τ_max`; otherwise the phase-pair
is *interval-incompatible* and the integral form is undefined for that
`(p, b, A, B)`. Snapshot is still defined in this case.

---

## 2. Definitions

### 2.1 Per-τ topology similarity (Spearman, on `D = 1/ρ`)

```
δ_S(τ; A, B; p, b) := 1 − ρ_Spearman( triu(D^{p,b,A}(τ)),  triu(D^{p,b,B}(τ)) )    ∈ [0, 2]
```

Measures rank-agreement on the upper triangle of communication distances.
`δ_S = 0` ⇔ identical ordering of pair-distances ⇔ same dendrogram
topology under UPGMA. Invariant to any monotone transformation of the
matrix entries within either phase.

### 2.2 Per-τ height similarity (Pearson, on `ρ`)

```
δ_P(τ; A, B; p, b) := 1 − ρ_Pearson( triu(ρ^{p,b,A}(τ)),  triu(ρ^{p,b,B}(τ)) )    ∈ [0, 2]
```

Measures linear (value) agreement on the upper triangle of the kernel
density. `δ_P = 0` ⇔ kernels agree up to an affine transform. Sensitive
to *both* ordering AND magnitude. Note the kernel is computed directly
on `ρ`, NOT `1/ρ`, so the comparison runs on the bounded similarity
(`ρ_{ij} ∈ [0, 1/N]` since `Tr ρ = 1`, `ρ_{ii} > 0`, `ρ_{ij} ≥ 0`)
without divergences.

### 2.3 Snapshot scalars (at `τ = τ_min`)

```
δ̂_S^{p,b}(A, B) := δ_S(τ_min^{p,b}(A,B); A, B; p, b)
δ̂_P^{p,b}(A, B) := δ_P(τ_min^{p,b}(A,B); A, B; p, b)
```

The "single-τ comparison at `τ = 1/λ_max`" of the user's request, with
`1/λ_max` extended to the phase-pair-symmetric `max(1/λ_max^A, 1/λ_max^B)`
so both phases are at or below their canonical resolution. This recovers
the paper-prescribed `1/λ_max` view as a special case when
`λ_max^A = λ_max^B`.

### 2.4 Functional integral scalars (over `I`)

```
Δ_S^{p,b}(A, B) := (log τ_max − log τ_min)^{−1}  ·  ∫_{τ_min}^{τ_max} δ_S(τ; A, B; p, b) d(log τ)
Δ_P^{p,b}(A, B) := (log τ_max − log τ_min)^{−1}  ·  ∫_{τ_min}^{τ_max} δ_P(τ; A, B; p, b) d(log τ)
```

Mean dissimilarity over `I` in log-τ measure (LRG-natural). Normalising
by interval log-length keeps `Δ_S, Δ_P ∈ [0, 2]` and dimensionless,
making them comparable across `(p, b)` cells with different intervals.
The `d(log τ)` choice — not `dτ` — prevents the upper edge of a
log-spaced interval from artificially up-weighting the integrand.

### 2.5 Reading the dual `(δ_S, δ_P)`

| pattern at τ | interpretation |
|:---|:---|
| `δ_S ≈ 0`, `δ_P ≈ 0` | trees agree in topology AND scale at τ |
| `δ_S ≈ 0`, `δ_P > 0` | same dendrogram topology, different absolute kernel scale at τ — structure shifted in scale only |
| `δ_S > 0`, `δ_P ≈ 0` | kernels are linearly correlated but rank-ordering differs — rare; usually a sign of tied values or near-degenerate ρ pattern |
| `δ_S > 0`, `δ_P > 0` | trees genuinely re-arranged — both topology and scale differ |

---

## 3. Properties

| # | Statement |
|:--|:---|
| 3.1 | `δ_S(τ; A, B), δ_P(τ; A, B) ∈ [0, 2]`. `δ_S(τ; A, A) = δ_P(τ; A, A) = 0`. Symmetric: `δ_X(τ; A, B) = δ_X(τ; B, A)`. |
| 3.2 | `δ̂_S, δ̂_P, Δ_S, Δ_P ∈ [0, 2]` by inheritance + the `1/(log τ_max − log τ_min)` normalisation. |
| 3.3 | **Rank invariance of `δ_S`.** If `D^{p,b,A}(τ) = f(D^{p,b,B}(τ))` element-wise with `f` strictly monotone, then `δ_S(τ; A, B) = 0`. So `δ_S` cannot detect a uniform monotone rescaling of the communication distance. By design — that's how it dodges `1/ρ` divergences. |
| 3.4 | **Affine sensitivity of `δ_P`.** Pearson is invariant under independent affine rescaling per phase; `δ_P = 0` iff `ρ^A = α ρ^B + β` for some `α > 0, β`. So `δ_P` detects ANY non-affine reorganisation, including rank flips that `δ_S` would miss. |
| 3.5 | **Why `δ_S` on `D = 1/ρ` and `δ_P` on `ρ` directly.** `D = 1/ρ` blows up for weakly-coupled pairs at small τ (when `ρ_{ij} → 0^+`); rank-based comparison absorbs this. Pearson on `D` would be dominated by those tail divergences, so we run it on the bounded similarity `ρ` instead. The two choices are complementary, not arbitrary: each metric is paired to the kernel where it is well-behaved. |
| 3.6 | **For fully-connected weighted graphs (our case),** "graph topology" is vacuous (every edge present); the non-trivial structure lives in the *weights*. `δ_S` reads the dendrogram-topology induced by those weights — non-trivial. `δ_P` reads the same plus the absolute weight magnitudes. Neither is reading "graph topology" in the usual sparse-graph sense. |
| 3.7 | **Cross-patient comparability.** `Δ_S, Δ_P` are dimensionless (correlation distances integrated against a dimensionless `d(log τ)`) and bounded `[0, 2]` — directly comparable across `(p, b)`. The interval `I^{p,b}(A,B)` is patient/band-specific but enters only as a normalisation factor; the integrand at any `τ` is a dimensionless similarity. |
| 3.8 | **Scale-pair vs. trace-pair semantics.** The "task-trace persists" claim is operationally `δ(test, post) < δ(pre, post)` at some τ window — post sits closer to test than to pre on the τ-axis. This claim is detected at the `(p, b)` level and does NOT require pooling into a cohort scalar. |
| 3.9 | **What this measure cannot detect:** |
|     | (a) Per-leaf changes — single-electrode neighbourhood reorganization. Use cophenetic-neighbourhood (`2026-04-25_cophenetic-neighbourhood.md`). |
|     | (b) Subtree-identity persistence — "did module X re-form in post". Use Cohesion-CBR (`2026-04-25_cohesion-cbr.md`). |
|     | (c) Direction of pair-wise distance shift (closer vs. farther). Use continuous-trace matrix (`2026-04-26_continuous-trace-matrix.md`). |
|     | (d) Within-phase variance / split-half null. Use H2e (`scripts/01_compute/hypothesis_tests/h2e_split_half.py`). |
| 3.10 | **Time complexity per `(p, b, A, B)`:** `O(n_τ · N^3)` from the matrix exponential at each τ; `n_τ = 100`, `N ≈ 100–120`. Empirical: ~5–10 s per phase-pair on dense `expm`. Cohort-of-4 cost: 4·5·6 = 120 phase-pair evaluations ≈ 10–20 minutes; acceptable for a primitive look. |

---

## 4. Caveats and failure modes

| # | Caveat | Mitigation |
|:--|:---|:---|
| 4.1 | `D = 1/ρ` blows up at small τ when `ρ_{ij}` is near-zero (weakly-coupled pairs before any diffusion). | Spearman is rank-based and immune. For sanity, log a numerical floor `ρ_{ij} ≥ ε = 1e−12` and assert `δ_S` is unchanged across `ε ∈ {1e−12, 1e−10, 1e−8}`. |
| 4.2 | `τ\*^{p,b,φ}` per phase may differ enough that `τ_min ≥ τ_max` and `I` is empty. | Report per-(p, b) interval lengths; flag empty intervals as "interval-incompatible — snapshot only". Do not paper-over with arbitrary τ. |
| 4.3 | `C(τ)` peak detection: in our continuous-spectrum case the peak is broad. Argmax on the cached `entropy_C` (length 399 over `τ ∈ [10⁻³, 10⁵]`) may be sensitive to numerical noise at the peak. | Smooth `C(τ)` with a Savitzky-Golay filter (window 5, polyorder 2) before argmax. Verify the peak is interior (`0 < argmax < len(C) − 1`) per phase; flag boundary peaks. |
| 4.4 | Same-probe bias: contacts on the same sEEG probe carry trivially high coupling, dominating `ρ` at all τ (per `probe-bias-guide.md`). | Run a cross-probe-only re-evaluation (zero out same-probe entries before computing ρ) as a second-pass diagnostic. Document delta. |
| 4.5 | Pat_03 outlier: 1024 Hz (others 2048 Hz). MSC-era values were 3× higher; ImCoh-era values renormalised but `λ_max` and `τ\*` may still sit on a different scale. | Include Pat_03 in the primitive set explicitly (as the user requested) and flag it distinctly in the figure. |
| 4.6 | `λ_max^{p,b,φ}` per phase requires the giant component, which may differ across phases of `(p, b)`. | Restrict to the *intersection* of giant-component electrode sets across all four phases; document drop-rate per `(p, b)`. Failure on > 5 % drop is a flag to widen the FC threshold or skip the cell. |
| 4.7 | `ρ_{ij} > 0` strictly only on the giant component of the underlying graph. Disconnected components produce `ρ_{ij} = 0` and `D_{ij} = ∞`. | Operate on the giant component only (already enforced by `compute_lrg_analysis` via `get_giant_component`). |
| 4.8 | Pearson on `triu(ρ)` is dominated by the largest-magnitude entries; a few high-coupling pairs can drive the correlation. | Report Pearson computed on `triu(ρ)` and on `triu(log ρ)` (variance-stabilised) as a sensitivity check. The `log ρ` version IS the "−log K_{ij}" framework-guide §10 #2 candidate — but only as a sanity check inside this scope, not as the canonical metric (per user 2026-04-28: study −log behaviour separately before using it as a canonical info-distance). |

---

## 5. Pseudocode

```
function FUNCTIONAL_TREE_DISTANCE(p, b, A, B, n_tau = 100, savgol_window = 5):

    # 1) Load FC and giant-component intersection across the four phases
    Â_phases := { φ : load_fc_matrix(p, φ, b, "imcoh_abs") for φ ∈ Φ }
    nodes    := ⋂_{φ ∈ Φ} giant_component_nodes( graph_from(Â_phases[φ]) )
    L_X(φ)   := laplacian( restrict_to(Â_phases[φ], nodes) )
    spec(φ)  := eigenvalues_sorted_ascending( L_X(φ) )

    # 2) Per-phase τ\* via argmax of smoothed C(τ)
    for φ ∈ {A, B}:
        load entropy_tau, entropy_C from cached LRGResult(p, φ, b, "imcoh_abs")
        C_smooth := savgol_filter( entropy_C, window = savgol_window, polyorder = 2 )
        i_star   := argmax( C_smooth )
        assert 0 < i_star < length(C_smooth) − 1   else flag boundary
        τ_star(φ) := entropy_tau[ i_star ]   # NB: entropy_C is len(τ)−1, align by midpoint

    # 3) Meaningful interval
    τ_min := max( 1 / max(spec(A)),  1 / max(spec(B)) )
    τ_max := min( τ_star(A), τ_star(B) )
    if τ_min ≥ τ_max:
        emit "interval-incompatible"; integral undefined; continue with snapshot only

    # 4) τ-grid, log-spaced
    τ_grid := logspace( log10(τ_min), log10(τ_max), n_tau )

    δ_S := array(n_tau);  δ_P := array(n_tau)

    # 5) Per-τ correlation-distance evaluation
    for i, τ in enumerate(τ_grid):
        K_A := matrix_exp( -τ · L_X(A) );   ρ_A := K_A / trace(K_A)
        K_B := matrix_exp( -τ · L_X(B) );   ρ_B := K_B / trace(K_B)

        # Communication distance (1/ρ), symmetrised, zero diagonal
        D_A := elementwise_inverse(ρ_A); symmetrise_max(D_A); zero_diagonal(D_A)
        D_B := elementwise_inverse(ρ_B); symmetrise_max(D_B); zero_diagonal(D_B)

        δ_S[i] := 1 − spearman_rho( upper_triangle(D_A), upper_triangle(D_B) )
        δ_P[i] := 1 − pearson_r ( upper_triangle(ρ_A), upper_triangle(ρ_B) )

    # 6) Snapshot scalars at τ = τ_min
    δ̂_S := δ_S[0];   δ̂_P := δ_P[0]

    # 7) Normalised log-τ integrals
    log_τ        := log( τ_grid )
    L            := log_τ[-1] − log_τ[0]
    Δ_S          := trapezoid( δ_S, log_τ ) / L
    Δ_P          := trapezoid( δ_P, log_τ ) / L

    return record {
        τ_grid, δ_S, δ_P,
        δ̂_S, δ̂_P, Δ_S, Δ_P,
        τ_min, τ_max, τ_star_A, τ_star_B,
        λ_max_A := max(spec(A)),  λ_max_B := max(spec(B)),
        n_nodes := |nodes|
    }
```

Loop bounds explicit: `i ∈ [0, n_τ)`. No hidden iteration. All linear
algebra is dense `O(N³)` per τ via `scipy.linalg.expm` (acceptable for
N ≈ 100–120).

---

## 6. Visualization spec

### 6.1 Primary figure: per-(patient, band) `δ(τ)` curves

- **Layout.** 4 patients (Pat_02, Pat_03, Pat_06, Pat_13) as columns × 5 bands (delta, theta, alpha, beta, gamma_low) as rows. 20 panels total. Each panel ≈ 3"×2".
- **Axes.** x = `log10 τ`. y = `δ` (`[0, 2]`, but typical range `[0, 1]`). Log-x is the LRG-natural axis; share x within each row, share y within each column.
- **Curves per panel.** 6 phase-pairs × 2 metrics = 12 curves. Phase-pair distinguished by hue (categorical, 6-class palette ordered by canonical comparison: `(pre, learn) (pre, test) (pre, post) (learn, test) (learn, post) (test, post)`). Metric distinguished by linestyle: solid = `δ_S` (rank), dashed = `δ_P` (height).
- **Annotations.**
  - Vertical line at `τ_min` (snapshot τ) — black, 0.5 pt.
  - Shaded band over `[τ_min, τ_max]` — light grey, 30 % alpha.
  - Title: `Pat_NN · band · n_nodes = M` (no `suptitle`).
  - Pat_03 title in italic / different colour to flag the 1024 Hz outlier.
- **Reading rules.**
  - **Task-trace persistence at scale τ.** `δ(test, post)(τ) < δ(pre, post)(τ)` ⇒ at this τ, post-rest sits closer to test than to pre.
  - **Topology vs height divergence.** Solid `δ_S` and dashed `δ_P` track each other ⇒ both ordering and magnitude agree. Solid stays low while dashed rises ⇒ topology preserved, scale shifted. Solid rises while dashed stays low ⇒ rank-flip without affine breakdown (rare; flag).
  - **Where the structure lives.** Curves with low `δ` over a wide τ window ⇒ structurally similar across that whole scale range. Sharp dips at narrow τ ⇒ scale-localised similarity.
- **No PDF + PNG.** PDF only per never-rule. Rasterise the curves at `dpi=200` if they get heavy.

### 6.2 Companion figure: scalar heatmaps

- **Layout.** 4 panels (one per patient), rows = bands, cols = 6 phase-pairs.
- **Cell content.** Each cell split horizontally: top half `Δ_S` (rank-integral), bottom half `Δ_P` (height-integral). Diverging colormap anchored at the cohort-median `Δ`, range `[0, 1]`.
- **Annotations.** Print `δ̂_S / δ̂_P` (snapshot) as a small text in each cell. Highlight cells where `Δ(test, post) < Δ(pre, post)` with a black border (the task-trace claim at the integral level).

### 6.3 Sanity figure: τ-interval diagnostic

- **Layout.** 4 panels (one per patient), x = band index, y = `log10 τ`.
- **Per band, four vertical line segments**, one per phase, from `1/λ_max^φ` (bottom) to `τ\*^φ` (top), color-coded by phase (rest_pre = blue, task_learn = green, task_test = orange, rest_post = red).
- **Read.** Overlap structure across phases tells us whether `I^{p,b}(A,B)` is empty for any phase-pair. Lines with no overlap = "interval-incompatible" cells.

---

## 7. Connection to prior tools

| Prior tool | What it measures | Relation to this scope |
|:---|:---|:---|
| `h2c_ultrametric_drift` (Spearman ρ on `Δ_task` vs `Δ_rest`) | Per-patient correlation of pair-wise distance shifts at `τ = 1/λ_max` | This scope generalises to a τ-curve and adds the height-aware `δ_P` view. H2c is fixed-τ + cophenetic; here we work on `D(τ)` and `ρ(τ)` directly. |
| `2026-04-26_continuous-trace-matrix` | Per-pair signed shifts `Δ_task, Δ_rest` and sign-agreement maps | Both look at the raw distance matrix pre-dendrogram. Continuous-trace measures pre→test→post co-shifts (per-pair signed); functional-tree measures phase-pair similarity over τ (per-(phase-pair, τ) magnitude). Complementary. |
| `tree_distance.py` (KC, MC, weighted RF) | Single-τ post-UPGMA tree-structure metrics | Single-τ, post-projection. Functional-tree avoids the UPGMA projection by working on `D(τ)` and `ρ(τ)` directly across τ. |
| `2026-04-25_module-retention-landscape` | Discrete subtree retention via Jaccard | Hard-thresholded; superseded for cohort census. Functional-tree is the continuous, no-threshold counterpart at the *whole-tree* level (vs. per-subtree level for MRL). |
| `2026-04-25_cohesion-cbr` | Soft per-leaf affinities | Per-leaf classification; functional-tree is per-(phase-pair, scale) global similarity. Complementary, not redundant. |
| `partition_multiscale` (Δ_VI, Δ_H, Δ_NMI per k) | Cut-based partition similarity per cluster count k | Cut-based — has the k-arbitrariness problem (per user 2026-04-28). Functional-tree replaces the k axis with the τ axis (physical, not combinatorial). |

Subsumes (in spirit, not in implementation): cut-based Δ_VI/Δ_H/Δ_NMI over k. Complements: continuous-trace matrix, CBR family. Does NOT replace: per-leaf or per-subtree localisations.

---

## 8. Implementation plan

### 8.1 Library

New module `src/lrg_eegfc/utils/metrics/functional_tree_distance.py` exporting:

- `tau_interval(L_A, L_B, entropy_C_A, entropy_tau_A, entropy_C_B, entropy_tau_B, savgol_window=5) -> (τ_min, τ_max, τ\*_A, τ\*_B, λ_max_A, λ_max_B)`
- `delta_curves(L_A, L_B, τ_grid) -> (δ_S, δ_P)`
- `compute_functional_tree_distance(p, b, A, B, fc_method="imcoh_abs", n_tau=100, savgol_window=5, save_cache=True) -> FunctionalTreeDistanceResult` (dataclass with all fields from the pseudocode return).

**Reuses (no new primitives):**
- `compute_laplacian_properties(G, tau)` from `lrgsglib.utils.lrg.spectral` — for ρ(τ), `1/ρ(τ)`, `λ_max`.
- `entropy(G, ...)` from `lrgsglib.utils.lrg.infocomm` — for `S(τ)`, `C(τ)`, `entropy_tau`.
- `load_fc_matrix(p, φ, b, "imcoh_abs")` from `lrg_eegfc.workflow.fc`.
- `load_lrg_result(p, φ, b, "imcoh_abs")` from `lrg_eegfc.workflow.lrg` — to reuse the cached `entropy_tau, entropy_C`.
- `scipy.stats.spearmanr`, `scipy.stats.pearsonr`.
- `scipy.linalg.expm`.
- `scipy.signal.savgol_filter`.
- `scipy.integrate.trapezoid`.
- `np.fill_diagonal`, `np.maximum`, etc.

No private re-implementation of correlation, exponentiation, or smoothing. Per `coding-rules.md`.

### 8.2 Compute script

`scripts/01_compute/diagnostics/diag_functional_tree_distance.py` driving:

```
patients := [Pat_02, Pat_03, Pat_06, Pat_13]
bands    := list(BRAIN_BANDS.keys())
pairs    := ALL_PHASE_PAIRS  # 6 unordered phase-pairs
for (p, b, (A, B)) ∈ patients × bands × pairs:
    result := compute_functional_tree_distance(p, b, A, B)
    save result to data/cache/functional_tree_distance/<p>/<b>_<A>-<B>.npz
```

### 8.3 Visualization script

`scripts/01_compute/figures_embedded/fig_functional_tree_distance_primitive.py` producing the 20-panel `δ(τ)` figure described in §6.1; output `data/outputs/figures/2026-04-28_functional_tree_distance_primitive.pdf`. Companion heatmap and τ-interval diagnostic as additional pages or sibling PDFs.

### 8.4 Cache layout

```
data/cache/functional_tree_distance/
├── Pat_02/
│   ├── delta_rest_pre-task_learn.npz
│   ├── delta_rest_pre-task_test.npz
│   ├── delta_rest_pre-rest_post.npz
│   ├── delta_task_learn-task_test.npz
│   ├── delta_task_learn-rest_post.npz
│   ├── delta_task_test-rest_post.npz
│   ├── theta_rest_pre-task_learn.npz
│   └── …
├── Pat_03/  …
└── …
```

`.npz` keys: `tau_grid, delta_S, delta_P, delta_hat_S, delta_hat_P, Delta_S, Delta_P, tau_min, tau_max, tau_star_A, tau_star_B, lambda_max_A, lambda_max_B, n_nodes, patient, band, phase_A, phase_B, fc_method, savgol_window, n_tau`.

---

## 9. Open questions

1. **Smoothing window for `τ\*` detection.** Default `savgol_window = 5`; sensitivity to window size and polyorder not yet quantified. Defer until first compute lands; revisit if argmax flips with `window ∈ {3, 5, 7, 11}`.
2. **Common giant-component restriction across four phases.** Empirical drop-rate per `(p, b)` not yet known. Cohort sanity table to be produced as part of the first compute pass; flag any `(p, b)` cell with `drop > 5 %`.
3. **τ\* convention for phase-pair interval.** `τ_max = min(τ\*_A, τ\*_B)` (intersection, default) vs. `max` (union) vs. `τ\*_{rest_pre}` (baseline-anchored). Sensitivity to be checked on Pat_02 first by re-running with the alternative rule and comparing curves.
4. **Probe-bias control.** Re-evaluate on cross-probe-only ρ as a second-pass diagnostic per `probe-bias-guide.md`. Deferred to a follow-up scope (this scope = primary, no probe control yet).
5. **Cohort generalisation.** This scope is primitive (`|P| = 4`); generalise to `|P| = 10` only after the primitive figure is interpretable and the open questions above are resolved.
6. **`log ρ` Pearson sanity.** Run alongside the canonical `δ_P(ρ)` as a variance-stabilised companion — does NOT replace it as the canonical metric. Per user 2026-04-28: study `−log K` behaviour separately before adopting it as a canonical info-distance.
7. **Within-phase null.** What does `δ` look like between two split-halves of a single phase's data? H2e provides this for the current per-pair drift; extending to functional-tree-distance is a separate scope.

---

## 10. Cross-references

- LRG framework (single source of truth): `.agents/guides/02_methods/lrg-framework-guide.md` §2 (formulas), §3 (dendrogram), §5 (τ choices), §10 (open questions).
- ImCoh derivation: `.agents/guides/02_methods/imcoh-guide.md`.
- Probe-bias controls: `.agents/guides/02_methods/probe-bias-guide.md`.
- Continuous-trace matrix (the per-pair sibling): `.agents/guides/task-persistence-investigation/2026-04-26_continuous-trace-matrix.md`.
- Task-trace canonical (P/T/R/RA): `.agents/guides/task-persistence-investigation/2026-04-25_task-trace-canonical.md`.
- Decision rules: `.agents/guides/task-persistence-investigation/2026-04-29_decision-rules.md`.
- Coding rules (library-first): `.agents/guides/04_rules/coding-rules.md`.
- Never-always list: `.agents/guides/04_rules/never-always-list.md`.
