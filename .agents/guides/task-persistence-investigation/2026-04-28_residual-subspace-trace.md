---
name: 2026-04-28_residual-subspace-trace
type: scope
era: IMCOH_ABS × COHORT_N10
status: superseded
created: 2026-04-28
updated: 2026-04-28
pointers:
  - 2026-04-27_perspective-multidirectional-multiscale.md
  - 2026-04-26_continuous-trace-matrix.md
  - .agents/reports/2026-04-28_residual-subspace-diagnostic.md
---

> **2026-04-28: Superseded by failed diagnostic.** All three pre-registered
> kill conditions hit on the 3-patient × 6-band × 12-τ̃ scan. Path closed.
> See [`2026-04-28_residual-subspace-diagnostic.md`](../../reports/2026-04-28_residual-subspace-diagnostic.md)
> for the verdict. Scope kept for the record; no successor.

# Residual-subspace alignment — group-wise multiscale task trace

**Renormalization head.** A group-wise task trace is the leading principal
mode of the cross-phase residual ultrametric matrix at scale `τ̃`: nodes with
large loading on `u₁` of `D^post − D^pre` are the ones whose pairwise
ultrametric distances reorganised together. The trace claim becomes "the
top-`k` eigen-subspaces of `S(τ̃) = D^task − D^pre` and `R(τ̃) = D^post − D^pre`
align". This avoids the cluster-machinery (Ψ-as-selector, Ψ_L) which the
2026-04-28 Ψ τ-scan diagnostic showed has no purchase on imcoh_abs trees, and
recovers a group-wise narrative that pair-level residual tests dilute.

## Notation

- Patients `p ∈ COHORT_N10`; bands `b ∈ {δ, θ, α, β, γ_l, γ_h}`; phases
  `φ ∈ {pre = rest_pre, task = task_test, post = rest_post}`.
- For each `(p, b, φ)`: weighted FC matrix `A^φ` (giant component on `|ImCoh|`),
  Laplacian `L^φ = D^φ − A^φ`, eigendecomposition `L^φ = U^φ diag(λ^φ) (U^φ)ᵀ`
  with `λ^φ_max = max(λ^φ)`, `λ^φ_gap = min(λ^φ : λ^φ > λ^φ_max·10⁻¹⁰)` (Fiedler
  on the giant component, ignoring the trivial constant mode).
- Dimensionless resolution `τ̃ ∈ [1, min_{φ}(λ^φ_max / λ^φ_gap)]`. Each phase
  evaluates at `τ_φ = τ̃ / λ^φ_max` so the same `τ̃` represents the same
  fraction of the network-intrinsic spectral range in every phase.
- Heat-kernel density `ρ̂^φ(τ̃) = e^{−τ_φ L^φ} / Tr(e^{−τ_φ L^φ})`
  (Villegas 2025 §I, eq. 1).
- Communication / ultrametric distance `D^φ_ij(τ̃) = (1 − δ_ij) / ρ̂^φ_ij(τ̃)`,
  proven ultrametric on connected networks (Villegas 2025 §I).
- Common node-subset `V*(p, b)` = nodes present in the giant component of all
  three phases. All matrices are restricted to `V* × V*`.

## Definitions

**Residual matrices** (per `(p, b, τ̃)`, on `V* × V*`):

```
S(τ̃) = D^task(τ̃) − D^pre(τ̃)
R(τ̃) = D^post(τ̃) − D^pre(τ̃)
```

Both real, symmetric. Eigendecompose:

```
S = U^S Λ^S (U^S)ᵀ      R = U^R Λ^R (U^R)ᵀ
```

with eigenvalues sorted by `|λ^S_i| ≥ |λ^S_{i+1}|`, similarly for `R`.

**Subspace alignment**: principal-angle / mean canonical correlation between
the top-`k` eigenvector subspaces:

```
α_k(τ̃) = (1/k) · Σ_{i=1}^{k} σ_i( (U^S_{1:k})ᵀ U^R_{1:k} ) ∈ [0, 1]
```

where `σ_i` are the singular values of the `k × k` matrix product. `α_k = 1`
iff `span(U^S_{1:k}) = span(U^R_{1:k})`. Computed for `k ∈ {1, 3, 5}`.

**Energy concentration** (per residual `M ∈ {S, R}`):

```
β_M(τ̃) = (λ^M_1)² / Σ_i (λ^M_i)² ∈ [0, 1]
```

`β_M ≈ 1` ⇒ residual dominated by a single mode (a coherent group). Small
`β_M` ⇒ diffuse reorganisation, no single group exists.

**Phase-label permutation null**: for `n_perm` permutations `π` of
`{pre, task, post}`, recompute `S^π, R^π` with re-labelled matrices and
compute `α_1^π(τ̃)`. Compare observed `α_1` to median of `{α_1^π}`.

## Properties

- **Range**: `α_k ∈ [0, 1]`, `β_M ∈ [0, 1]`. Both are scalar per `(p, b, τ̃)`.
- **Cohort-aggregable**: scalars in `[0, 1]` aggregate cleanly across patients
  without leafset matching, anatomy mapping, or partition selection.
- **Group-wise**: `u_1^R` defines a continuous loading per node — the "trace
  cluster" is the heavy-tail of `|u_1^R(i)|`, no hard k-cut required.
- **What it cannot detect**: changes that have zero net effect on pairwise
  ultrametric distances (e.g. consistent permutation of identical heights);
  group-wise reorganisation orthogonal to the leading modes (caught by `α_3`,
  `α_5`); non-linear / multi-modal restructuring that splits the residual
  energy across many comparable eigenvalues (low `β`, `α` interpretation breaks).
- **Permutation invariance**: invariant to relabelling of nodes (since both
  `S` and `R` are computed on the same `V*` indexing). Anatomy must be
  mapped *post-hoc* via electrode coordinates — not a bug, since the analysis
  is anatomy-agnostic by design until the brain-figure step.

## Caveats and failure modes

| Caveat | Mitigation |
|---|---|
| Eigenvector signs are arbitrary | Use `\|u\|` for loadings; `α` uses subspace SVD which is sign-blind. |
| Multiple comparable eigenvalues collapse "the group" notion | Report `β`. If `β < 0.1` (ten or more comparable modes), don't call it a group. |
| Different node sets across phases | Restrict to `V*(p, b)` from the start. Drop `(p, b)` if `\|V*\| < 10`. |
| Surrogate-null too coarse with 5 perms | Acceptable for diagnostic; widen to circular-shift on time-series + recompute LRG before any cohort claim. |
| `τ̃`-window dispersion across patients | Common τ̃-grid intersected per patient via `min_φ(λ_max/λ_gap)`. Report which patients hit which `τ̃`-window. |
| Group-wise but anatomically scattered | Brain-figure check: `u_1^R` loadings must concentrate in coherent Desikan-Killiany regions to support a "group" claim. |

## Pseudocode

```
INPUT: patients P, bands B, phases Φ = {pre, task, post}, n_tau, n_perm
FOR each (p, b) IN P × B:
    FOR each φ IN Φ:
        A_φ ← load_fc_matrix(p, φ, b, "imcoh_abs")
        G_φ ← giant_component(graph_from(|A_φ|))
        (eigvals_φ, eigvecs_φ) ← eigendecompose(laplacian(G_φ))
        λmax_φ ← max(eigvals_φ); λgap_φ ← second-smallest non-trivial eigenvalue
    V* ← intersection over φ of nodes_kept_φ
    IF |V*| < 10: SKIP (p, b)
    τ̃_max ← min_φ(λmax_φ / λgap_φ)
    τ̃_grid ← geomspace(1.0, τ̃_max, n_tau)
    FOR each τ̃ IN τ̃_grid:
        FOR each φ IN Φ:
            ρ_φ ← rho_tau(eigvals_φ, eigvecs_φ, τ̃ / λmax_φ)
            D_φ ← (1 - I) / ρ_φ           # ultrametric distance
            D_φ ← restrict_to(V*) (D_φ)
        S ← D_task − D_pre;  R ← D_post − D_pre
        S ← (S + Sᵀ)/2;       R ← (R + Rᵀ)/2
        (Λ^S, U^S) ← eigendecompose(S); sort by |Λ^S| descending
        (Λ^R, U^R) ← eigendecompose(R); sort by |Λ^R| descending
        FOR k IN {1, 3, 5}:
            σ ← singular_values((U^S[:, :k])ᵀ · U^R[:, :k])
            α_k ← mean(σ)
        β_S ← (Λ^S_1)² / sum((Λ^S_i)²)
        β_R ← (Λ^R_1)² / sum((Λ^R_i)²)
        # Permutation null
        FOR i IN 1..n_perm:
            π ← random_permutation(Φ)
            (recompute S^π, R^π using re-labelled D matrices)
            α_1^π_i ← subspace_align(top-1 of S^π, top-1 of R^π)
        null_α_1 ← median(α_1^π)
        EMIT (p, b, τ̃, α_1, α_3, α_5, β_S, β_R, null_α_1, u^R_1, u^S_1)
```

## Visualization spec

**Page 1 — α_1(τ̃) vs τ̃, observed vs null.** Grid: rows = patients
(Pat_02, Pat_06, Pat_03), cols = 6 bands. Each panel: solid line = observed
`α_1(τ̃)`, dashed line = `null_α_1(τ̃)` median over phase-label permutations.
Reading rule: solid above dashed across a contiguous τ̃-window = candidate
trace at that scale. Solid below dashed = anti-trace (alignment by chance
or worse).

**Page 2 — β_S(τ̃), β_R(τ̃) vs τ̃.** Same grid. Two lines per panel
(`β_S` solid, `β_R` dashed). Reading rule: `β_R > 0.3` ⇒ residual has a real
group; `β_R < 0.05` ⇒ diffuse, no group, kill.

**Page 3 — α_3, α_5 vs τ̃.** Same grid. Two lines per panel (`α_3` and `α_5`).
Reading rule: if `α_3 ≫ α_1`, the trace lives in a multi-mode subspace, not a
single mode — group-wise but multi-component.

**Page 4 — `u_1^R` loadings sorted by magnitude**, one panel per `(patient, band)`
at the τ̃-bin where observed `α_1` peaks above null. X-axis: rank, y-axis:
`|u_1^R|`. Heavy tail = trace group present; flat distribution = no group.
Top-loading nodes annotated with their channel labels for anatomical inspection.

## Connection to prior tools

| This measure | Complements | Subsumes / replaces |
|---|---|---|
| α_k(τ̃), β_M(τ̃) on D residuals | Continuous-trace matrix (2026-04-26) — same observable, different summary. CTM = per-pair sign agreement; here = subspace-level agreement. | Pair-level matrix tests (Frobenius, Mantel) when the goal is group-wise narrative. |
| `u_1^R` loadings | Cohesion-CBR per-patient subtree census | Ψ_L-based cluster filtering (which the Ψ τ-scan empirically killed). |
| Phase-label permutation null | Circular-shift surrogate phase null on H2c | Provisional only; not the cohort-grade null. |

Does **not** replace MRL / CBR / consensus subtree — those are
*subtree-identity* tests with anatomical leafsets. Subspace alignment is a
*subspace-of-distances* test; the two answer different questions.

## Implementation plan

- Diagnostic script: `scripts/01_compute/diagnostics/diag_residual_subspace.py`
  (3 patients × 6 bands × 12 τ̃-points, `n_perm = 5`).
- Library entry points reused: `lrg_eegfc.workflow.fc.load_fc_matrix`,
  `lrgsglib.nx_patches.funcs.get_giant_component`. No new helpers promoted
  unless reused outside this script.
- Outputs: `data/audit/residual_subspace/{alpha_beta_grid.csv,
  loadings.npz, diagnostic.pdf}`.
- Verdict: `.agents/reports/2026-04-28_residual-subspace-diagnostic.md`.

**Kill conditions, encoded in the verdict:**
- `β_R < 0.1` everywhere → diffuse residual, no group → kill.
- Observed `α_1 ≤ null_α_1` cohort-wide → no task-specific subspace
  alignment → kill.
- Top-`u_1^R` loadings scatter across non-coherent anatomy → "group"
  is statistical artefact → don't headline.

## Open questions

- **Surrogate null grade.** Phase-label permutation here is a coarse null;
  the cohort-grade null is circular-shift on the underlying time-series with
  full LRG re-derivation. To be added if the diagnostic survives.
- **`τ̃`-grid density.** `n_tau = 12` is for diagnostic visibility; cohort
  run may need finer (≥ 25) to see a τ̃-window of width < 0.5 decades.
- **k choice for `α_k`.** `k = 1` is the cleanest narrative; `k = 3, 5`
  reported as robustness checks. If the trace is multi-modal,
  `k > 1` may dominate the picture and reframe the claim.
- **Anatomical aggregation.** Once `u_1^R` is in hand, mapping to
  Desikan-Killiany regions per patient and aggregating to cohort-region
  loadings is a follow-up; not in this diagnostic.
