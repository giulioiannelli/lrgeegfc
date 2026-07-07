---
name: diffusion-distance-trace
type: scope
era: IMCOH_ABS
status: draft
created: 2026-07-06
updated: 2026-07-06
pointers:
  - .agents/guides/task-persistence-investigation/2026-04-26_continuous-trace-matrix.md
  - .agents/guides/task-persistence-investigation/2026-06-22_tau-sensitivity-cophenetic-trace.md
  - .agents/guides/task-persistence-investigation/2026-06-12_cophenetic-consolidation-arc.md
  - .agents/guides/task-persistence-investigation/2026-04-25_task-trace-canonical.md
---

# Diffusion-distance trace — ρ^diff(τ)

**ρ^diff(τ) is the τ-resolved diffusion-distance version of ρ_split. Same
Spearman-on-Δ trace functional; the per-pair distance is the LRG propagator's
own multiscale geometry `D_τ(i,j) = √(Σ_{a≥2} e^(−2τλ_a)(u_a(i)−u_a(j))²)`
instead of a single-τ cophenetic ultrametric. The claim it tests: a
band-specific, learning-vs-test-dissociating task trace lives at *intermediate*
diffusion scales — provably invisible to the raw adjacency matrix (whose content
is the small-τ strength+edge regime) and to a Grassmann subspace probe
(eigenvalue-blind). Falsifier stated up front: if ρ^diff(τ) only fires as τ→0,
it adds nothing over the cophenetic trace (audit_121) and we say so.**

> Status: `draft`. Exploration only — this measure does **not** touch ρ_sym, the
> preprint, or any core artifact. It is not registered in the README Live-measures
> table and will not be until (if ever) it earns a headline. Scope-before-code
> per the folder rule; the reuse-map (§10) and prior-attempts (§9) sections are
> reconciled against the 2026-07-06 grounding sweep.

---

## 0. Grounding reconciliation (2026-07-06 archive + code sweep)

**This measure is the never-run E2, resurrected and upgraded.** The 2026-04-29
"eigenvector-direct pivot" scoped exactly this object — E2 "diffusion-map
embedding drift": Coifman–Lafon embedding via eigenvectors weighted by
`exp(−τλ/2)`, per-node cross-phase distance (`audit_27_e2_diffusion_map_drift.py`,
scoped, **never built**). It was frozen when its sibling **E1 (Grassmann
subspace, node-resolved, full top-k) went cohort-negative in every band**
(`.agents/reports/2026-04-29_e1-cohort-verdict.md`). So "use the eigenvectors
directly" is *necessary but not sufficient* — it was tried in its hard-cut form
and failed. ρ^diff differs from E1 by being **eigenvalue-weighted** (soft
`e^{−2τλ}`, not a hard k-cut) and from E2 by being fed through the now-canonical
cross-phase **ρ_split/ρ_sym** functional and gated by the mandatory
matched-strength null. A current methods doc
(`.agents/guides/02_methods/multiscale-notions-propagator-vs-grassmann.md:146`,
2026-07-04) names precisely this niche — "an eigenvalue-weighted subspace metric
weighting principal angles by `e^{−τλ}` … **Not run**" — as the documented open
interpolator between the surviving cophenetic probe (soft, full-spectrum) and the
dead Grassmann probe (hard, eigenvalue-blind).

**Two hard guardrails from prior τ / diffusion work — non-negotiable:**
1. **Collapse placebo** (audit_121, `.agents/reports/2026-06-22_tau-sensitivity-of-the-trace.md`).
   Past the Fiedler time the kernel goes low-rank and a *trace-free placebo*
   `ρ_indep` climbs ≈0.05 → ≈0.23 to **match** `ρ_split` — a false "coarse-τ
   strengthens the trace." Every τ-sweep here ships a collapse/placebo monitor
   (§6 C1); coarse-τ "gains" are assumed artifact until the placebo clears them.
   The genuine cophenetic trace is fine-scale and τ-flat.
2. **Magnitude → strength confound** (audit_73; the KC collapse, audit_65).
   Magnitude-weighted distances reintroduce strength structure and fail
   matched-strength (KC's β 10/10 headline died this way — ≈50% of the shift was
   pure strength through average-linkage). ρ^diff inherits ρ_split's **rank-based
   (Spearman) magnitude-blindness**, so the confound enters only through *which
   pairs rank where* — strength-dominated at small τ (§4), killed by
   matched-strength there. Intermediate τ is the only strength-independent regime.

**Reuse — the observable is nearly free.** Combinatorial-L̂ eigenpairs are
**already cached**: `load_lrg_result(p,φ,b,"imcoh_abs").eigenvalues/.eigenvectors`
(ascending λ, unit-norm columns, col 0 = trivial mode; `None` only for
pre-2026-04-29 legacy caches → fall back to `eigh(diag(W.sum(1))−W)`).
Matched-strength surrogate eigenpairs are **also cached**
(`matched_strength_surrogate_lrg/Pat_NN/{band}_{phase}_R200_swap20_seed20260511_imcoh_abs.npz`,
keys `eigvals (R,N)`, `eigvecs (R,N,N)`). Both the observable and its MS null
reduce to: load eigenpairs → form `D_τ` (the one new primitive) → the existing
`rho_split`/`ρ_sym` + `surrogate_p_value` + `wilcoxon_z` pipeline.

**Canonical functional = split-baseline ρ_sym.** The locked cross-phase functional
(`audit_67`, `audit_156`) splits rest_pre into independent halves A/B
(`Δ_task = X − preA`, `Δ_rest = R⁺ − preB`) and symmetrizes
`ρ_sym = ½[ρ(xA,yB) + ρ(xB,yA)]` to remove shared-baseline inflation.
- **Pass 1 (shape reconnaissance, THIS build):** shared full rest_pre baseline on
  cached eigenpairs; read **contrasts only** — band heterogeneity, learn-vs-test
  separation, hump location — which are robust to the common additive inflation.
  Not a claim; absolute ρ is not interpreted.
- **Pass 2:** upgrade to split-baseline ρ_sym + the collapse placebo +
  matched-strength before any cohort claim.

---

## 1. Motivation (why supersede the cophenetic snapshot)

The paper's thesis is that a **multiscale / higher-order** structural trace of
the task survives into `rest_post` that the raw adjacency matrix cannot resolve.
The current instrument — `ρ_sym`, a symmetric Spearman on the cophenetic
ultrametric `D_coph = cophenet(UPGMA(D(τ_max)))` — instantiates "multiscale" in a
weak way: it evaluates the propagator at a **single** scale `τ_max = 1/λ_max`,
then imposes a hierarchy through UPGMA agglomeration and collapses it to an
ultrametric. The multiscale ladder is manufactured by the linkage, not read from
the propagator's own scale-flow. Empirically it is a thin instrument (cohort
ρ ≈ 0.20, defended only by a 16.6× matched-strength margin, band gate marginal,
"compresses vs ρ_split").

This measure returns to the object the thesis actually names: the propagator
`e^(−τL)` as a function of τ, read through the diffusion geometry it induces on
the contacts, swept across scales, with **no** UPGMA and **no** ultrametric
collapse.

---

## 2. Notation

- Cohort `P = {Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15}`, `|P| = 10`, index `p`.
- Bands `B = {δ, θ, α, β, γ_l, γ_h}`, index `b`. (γ_h = 80–300 Hz per project
  definition, not textbook.)
- Phases `φ ∈ {R⁻ = rest_pre, L = task_learn, T = task_test, R⁺ = rest_post}`.
- FC matrix `W = W^{p,b,φ} ∈ ℝ^{N×N}`, symmetric, entrywise ≥ 0 (`imcoh_abs`),
  zero diagonal, `N = N_p` contacts. Pat_10 task rows `[53,54,55]` dropped at load.
- Strength (weighted degree) `s_i = Σ_j W_ij`; `D = diag(s)`.
- **Combinatorial Laplacian** `L = D − W` (the Villegas "fluid" Laplacian; the same
  operator the pipeline and the Grassmann probe diagonalise).
- Eigendecomposition `L = Σ_{a=1}^{N} λ_a u_a u_aᵀ`,
  `0 = λ_1 ≤ λ_2 ≤ … ≤ λ_N`, `{u_a}` orthonormal, `u_1 = 1/√N` (constant mode).
  `λ_2` = algebraic connectivity (Fiedler value); `λ_N = λ_max`.
- Upper-triangle pair index `m ↔ (i<j)`, `M = N(N−1)/2` pairs.
- τ-ladder `τ[1..K]`, log-spaced (see §5, §11); dimensionless scale `s = τ·λ_2`.

The operator is the combinatorial Villegas *fluid* Laplacian `L̂ = D̂ − W`
throughout — the same operator the LRG pipeline and the Grassmann probe
diagonalise. Strength-normalised / random-walk operators (`I − D^{-1}W`,
`I − D^{-1/2}WD^{-1/2}`) are a **different framework** (not the LRG propagator)
and are **out of scope**. Locked 2026-07-07.

---

## 3. Definitions

### 3.1 Heat kernel and diffusion coordinates (Coifman–Lafon 2006)

Heat kernel at scale τ:
```
K_τ = e^(−τL) = Σ_{a=1}^N e^(−τλ_a) u_a u_aᵀ.
```
Diffusion coordinate map (drop the constant mode a=1, which cancels in all
distances):
```
Ψ_τ(i) = ( e^(−τλ_a) u_a(i) )_{a=2..N}  ∈ ℝ^{N−1}.
```

### 3.2 Diffusion distance

```
D_τ(i,j)² = Σ_{a=2}^N e^(−2τλ_a) ( u_a(i) − u_a(j) )²
          = ‖ Ψ_τ(i) − Ψ_τ(j) ‖²
          = ‖ K_τ(i,·) − K_τ(j,·) ‖²      (rows of the heat kernel).
```
Domain: `i,j ∈ 1..N`. Range: `D_τ ≥ 0`, a genuine Euclidean metric on contacts
(triangle inequality holds). The a=1 term vanishes because `u_1` is constant.
Note the `e^(−2τλ_a)` weight (from squaring the kernel): the effective decay
scale is `2τ`.

### 3.3 Trace functional (drop-in ρ_split generalization)

For an arc through task reference `X ∈ {L, T}` (i.e. `R⁻ → X → R⁺`), per band,
per patient, per scale τ, on the pair vector over `m = 1..M`:
```
d^{φ}_m(τ) = D_τ^{p,b,φ}(i(m), j(m))                        # phase-φ distance vector
Δtask_m(τ) = d^{X}_m(τ) − d^{R⁻}_m(τ)                       # task-induced change
Δrest_m(τ) = d^{R⁺}_m(τ) − d^{R⁻}_m(τ)                      # rest-induced change

ρ^diff_{p,b,X}(τ) = Spearman_m( Δtask(τ), Δrest(τ) )        # tie-corrected
```
Sign convention (locked, `feedback_td_sign_convention`): **ρ^diff > 0 = TRACE**
(the pairs the task moved in diffusion geometry are the pairs that stayed moved
in `rest_post`); `< 0` = anti-trace; `≈ 0` = no trace. Identical reading to
ρ_split — only the underlying distance changed.

### 3.4 Cohort aggregation (no consensus scalar)

At each (b, X, τ):
```
median_p ρ^diff_{p,b,X}(τ),   IQR_p ρ^diff_{p,b,X}(τ)      # keep the spread
w, p_W = wilcoxon( {ρ^diff_{p,b,X}(τ)}_p , alternative='greater' )
```
The deliverable is the **τ-profile** (a curve per band per arc), never a pooled
scalar (pooling into a consensus scalar is project-forbidden). Per-patient
fluctuations are shown, not averaged away (`feedback_fluctuations_are_signal`).

### 3.5 Resolution axis (effective number of active modes)

```
n_eff(τ) = ( Σ_{a≥2} e^(−2τλ_a) )² / Σ_{a≥2} e^(−4τλ_a)     # participation ratio
```
`n_eff` runs from `N−1` (τ→0, all modes) to `→1` (τ→∞, Fiedler only). It gives a
per-patient-commensurable "how many modes are resolved" axis and defines the
collapse guard (§6 C1: interpret only where `median_p n_eff ≥ 2`).

---

## 4. The wedge — why this is irreducible to adjacency or Grassmann

Three regimes of τ, derived exactly:

**τ→0 is degenerate, not "the adjacency limit."** At τ = 0 every off-diagonal
pair has `D_0(i,j)² = 2` (using completeness `Σ_{a≥2} u_a(i)u_a(j) = δ_ij − 1/N`).
The geometry carries no pair information at τ = 0.

**Small τ is the strength + adjacency regime.** First order in τ:
```
D_τ(i,j)² ≈ 2 − 2τ ( s_i + s_j + 2 W_ij ) + O(τ²).
```
(from `Σ_{a} λ_a(u_a(i)−u_a(j))² = L_ii + L_jj − 2L_ij = s_i + s_j + 2W_ij`).
So at small τ the **ranking of pairs is dominated by summed endpoint strength**
plus the edge weight — precisely the quantity matched-strength preserves and the
raw adjacency matrix already carries. ⇒ small-τ ρ^diff is redundant with raw-FC
ρ_split and is the most matched-strength-vulnerable regime. This is the exact,
non-hand-wavy statement of "the adjacency matrix already sees this."

**Large τ is the Fiedler-collapse regime.** As τ grows, `K_τ → (1/N)11ᵀ +
e^(−2τλ_2) u_2 u_2ᵀ`, so `D_τ(i,j) → e^(−τλ_2)|u_2(i) − u_2(j)|` (rank-1). The
geometry degenerates to the single Fiedler coordinate — exactly what a
1-dimensional spectral cut / low-k Grassmann sees, and the source of the
audit_121 coarse-τ "collapse artifact."

**Intermediate τ is the genuinely higher-order regime.** Here `e^(−2τλ_a)`
weights many eigenvectors by a smooth, τ-tunable function of their eigenvalues.
This is the only regime in which the object is neither strength+adjacency nor a
Fiedler cut. **The paper's higher-order claim lives or dies at intermediate τ.**
- The **adjacency matrix** carries the small-τ content only (`W`, `s`).
- **Grassmann** `T_G*` compares `span{u_1..u_k}` by principal angles, swept over
  k — eigenVECTORS only, hard cut, eigenVALUES discarded, no scale continuum.
- **ρ^diff(τ)** is the eigenvalue-WEIGHTED, scale-continuous object. A task change
  that reweights eigenvalue spacing while fixing the dominant subspace shifts
  ρ^diff at intermediate τ while leaving Grassmann's subspace angles unchanged;
  a diffuse multi-mode change is integrated by the kernel while entrywise `W`
  comparison misses the coordination. That is the wedge.

---

## 5. Properties

- **Range / metric.** `ρ^diff ∈ [−1,1]`; `D_τ` is a metric (0 on diagonal,
  symmetric, triangle inequality). `n_eff ∈ [1, N−1]`.
- **Limits.** τ→0: uniform (degenerate, uninformative). small τ: strength+edge
  (§4). large τ: Fiedler-only (collapse). intermediate τ: multi-mode.
- **Invariance.** `D_τ` is invariant to orthogonal rotation within any degenerate
  eigenspace — the distance is a smooth sum over modes, so it is **stable even
  when individual eigenvectors are not** (near-degenerate λ). This is a
  robustness *advantage* over Grassmann-at-fixed-k, which is unstable exactly
  there. It is not invariant to strength rescaling (a feature: strength-sensitivity
  is confined to small τ and handled by matched-strength).
- **Complexity.** One `eigh` per (p,b,φ): `O(N³)`, `N ≈ 100–130` → trivial. Per τ,
  `G_τ = U diag(e^(−2τλ)) Uᵀ` then `D_τ² = diag(G)1ᵀ + 1 diag(G)ᵀ − 2G`: `O(N³)`
  but a single vectorised einsum; K τ-values → `O(K N³)`. Full observable
  (10 pat × 6 band × 4 phase × K≈30) is seconds. MS null multiplies by R (defer,
  reuse cached surrogate eigendecompositions).
- **Negative properties — what it CANNOT detect.** (a) Pure strength changes that
  leave the normalised geometry intact (invisible at intermediate τ — by design).
  (b) Anything finer than the mode structure at the collapse tail. (c) It is a
  *pairwise* summary — it does not localise to nodes on its own (node-resolved
  second pass deferred, §11). (d) Being Spearman-on-Δ it measures **co-directional
  reorganization, not magnitude** — same limitation as ρ_split; magnitude needs
  the `×null / %held-of-movers` reading (`rho_reading_key`).

---

## 6. Caveats & failure modes

- **C1 — Collapse tail (large τ).** `K_τ → rank-1`, ρ^diff → trivial Fiedler
  correlation. *Mitigation:* cap the ladder / interpretation at `median_p n_eff ≥ 2`;
  flag any signal that appears **only** in the collapsed tail as artifact — this is
  the audit_121 lesson made structural.
- **C2 — Strength confound (small τ).** `D_τ² ≈ 2 − 2τ(s_i+s_j+2W_ij)` (§4).
  *Mitigation:* matched-strength null; report intermediate-τ, not small-τ, as the
  claim; small-τ ρ^diff should track raw-FC ρ_split (a positive consistency check).
- **C3 — Lattice ties in Spearman.** Fully-connected FC ⇒ heavily tied pair
  vectors (cf. `perpair_cophenetic_trace_visualization`: ~98% tied for the
  cophenetic version). *Mitigation:* tie-corrected Spearman; visualize with the
  mover-density / t-map views, never a raw pair scatter.
- **C4 — Eigenvector degeneracy.** Near-degenerate λ ⇒ unstable individual `u_a`.
  *Mitigation:* none needed — `D_τ` is invariant to intra-eigenspace rotation
  (§5). Documented as an advantage, but verify numerically (perturb & recompute).
- **C5 — Same-probe bias.** Same-shaft contacts have trivially high FC and
  dominate the small-τ geometry (`probe-bias-guide`). *Mitigation:* run a
  zero-same-probe-edge sensitivity pass; intersect with `channel_labels.csv`.
- **C6 — τ-commensurability across patients.** Different spectra ⇒ raw τ not
  comparable. *Mitigation:* dimensionless `s = τλ_2` and/or the `n_eff` resolution
  axis; decide the primary axis after seeing the profiles (§11).
- **C7 — Operator is fixed to the combinatorial Laplacian.** The measure is
  defined on the Villegas fluid Laplacian `L̂ = D̂ − W` only — the LRG operator.
  Its `D_τ` is strength-scaled at fine τ (§4); that is handled by the
  matched-strength null, **not** by switching operators. Symmetric-normalised /
  random-walk Laplacians are a different framework (not the LRG propagator) and
  are out of scope. Locked 2026-07-07.
- **C8 — Not a new null.** ρ^diff reuses the **existing** matched-strength null,
  unchanged. No new surrogate construction is introduced (avoids the
  opaque-null trap, `feedback_no_opaque_matrix_nulls`).

---

## 7. Pseudocode

```
INPUT: cohort P, bands B, phases {R-,L,T,R+}, τ-ladder τ[1..K], fc_method=imcoh_abs

# ---- observable ----
for p in P:
  for b in B:
    for φ in {R-,L,T,R+}:
      W        ← load_fc_matrix(p, φ, b, imcoh_abs)     # N×N ≥0, zero-diag
      L        ← diag(rowsum(W)) − W                     # combinatorial Laplacian
      (λ[1..N], U[·,1..N]) ← eigh(L)                     # ascending, orthonormal cols
      for k in 1..K:
        E      ← exp(−2·τ[k]·λ[2..N])                    # drop constant mode a=1
        G      ← U[:,2..N] · diag(E) · U[:,2..N]ᵀ        # N×N Gram of Ψ_τ
        g      ← diag(G)
        d2     ← g·1ᵀ + 1·gᵀ − 2·G                       # squared diffusion distance
        Dvec[p,b,φ,k] ← sqrt(clip(d2,0,·))[upper-triangle]   # length-M vector
        neff[p,b,k]   ← (Σ E)² / Σ E²                     # resolution axis
    for X in {L,T}:                                       # the two arcs
      for k in 1..K:
        Δtask ← Dvec[p,b,X,k] − Dvec[p,b,R-,k]
        Δrest ← Dvec[p,b,R+,k] − Dvec[p,b,R-,k]
        rho[p,b,X,k] ← spearman(Δtask, Δrest)             # tie-corrected

# ---- cohort profiles (NO consensus scalar) ----
for b in B, X in {L,T}, k in 1..K:
  med[b,X,k], iqr[b,X,k] ← median_p, IQR_p  rho[·,b,X,k]
  wstat, pW[b,X,k]       ← wilcoxon(rho[·,b,X,k], alternative='greater')
collapse_ok[b,k] ← ( median_p neff[·,b,k] ≥ 2 )           # interpret only where true

# ---- matched-strength null (DEFERRED to after observable is inspected) ----
for p,b,φ: load R surrogates W_r → L_r → eigh → Dvec_r     # reuse cached eigendecomps
for b,X,k in intermediate window:
  rho_null ← { spearman(Δtask_r, Δrest_r) }_r
  pMS[b,X,k] ← mean( rho_null ≥ med/obs rho )              # upper tail
```

---

## 8. Visualization spec

- **Primary — per-band τ-profile.** x = `s = τλ_2` (log) or the `n_eff` resolution
  axis (decreasing modes → coarser, read right-to-left); y = ρ^diff. One line per
  band (cohort median), shaded IQR; two arcs (L vs T) as two panels or two
  linestyles. **Reading rules:** an intermediate-s *hump* that is band-specific =
  the higher-order trace; a monotone rise out of small s = fine-scale redundancy
  (the falsifier — no gain over cophenetic); horizontal contrast across bands =
  per-band heterogeneity; L-vs-T curve separation = learning-vs-test dissociation.
- **Significance marking.** Mark every matched-strength-significant τ cell
  *uniformly* — no `axhline(n/10)`, no pre-shaded windows, no threshold lines
  (`feedback_no_visual_threshold_lines`). Shade the `n_eff < 2` collapse tail as
  "do not interpret."
- **Secondary — head-to-head.** Per band, overlay the three probes on a common
  τ / resolution axis: ρ^diff(τ) curve; raw-FC ρ_split (a horizontal reference,
  it is scale-free); Grassmann `T_G*` extent (its k-window mapped to a scale).
  The figure must *show* where the propagator resolves structure the other two
  miss (intermediate τ) and where it does not (small τ, collapse tail).
- **Spread shown**, not just the median (`feedback_fluctuations_are_signal`).
- PDF only, `use_lrg_style()`, no `suptitle`, no watermark, transparent, vector.

---

## 9. Connection to prior tools

| Prior tool | Relation to ρ^diff(τ) |
|---|---|
| **Continuous-trace matrix (2026-04-26) / ρ_split / ρ_sym** | *Same* Δ_task/Δ_rest Spearman functional. They use a single-τ distance (`D(τ_max)=1/ρ̂`, or its cophenetic ultrametric). ρ^diff **swaps the distance** to the τ-resolved diffusion metric and **sweeps** it. Subsumes ρ_split at small τ (strength+edge); adds the intermediate-τ eigenvalue-weighted geometry ρ_split cannot see. |
| **τ-sensitivity of the cophenetic trace (2026-06-22, audit_121)** | Swept τ *inside the cophenetic pipeline*; found trace fine-scale + τ-robust, coarse-τ = collapse. ρ^diff re-opens the τ-sweep on the diffusion geometry (no UPGMA/cophenetic distortion). audit_121's collapse = my C1; its fine-scale finding = my falsifier. If ρ^diff *also* only fires at small τ, it confirms audit_121 and I report no gain — this reconciliation is mandatory. |
| **Grassmann `T_G*`** | Eigenvector-subspace, eigenvalue-blind, hard-k. The head-to-head comparator (§4). The irreducible content of ρ^diff = whatever survives at intermediate τ when Grassmann is null. |
| **Cophenetic consolidation arc (2026-06-12, N2) / inference-mark** | The learn/test + encoding/inference structure I mirror with the two arcs. A propagator-native L-vs-T dissociation at intermediate τ would be the diffusion analog of the N2 arc. |
| **Archived spectral / subspace attempts** (`eigenmode_embedding`, `e1-spectral-subspace-alignment`, `perspective-multidirectional-multiscale`, `intrinsic-scale-reorganization`, the removed S(τ)/C(τ) "L2 rung") | ⏳ reconcile with 2026-07-06 grounding sweep. Working hypothesis: those were eigenvalue-*only* (spectral-function) or subspace-*only* probes that failed on the continuous spectrum's lack of scale-ID; a node-resolved, eigenvalue-*weighted* diffusion distance needs no scale-ID and may sidestep that specific failure. **Must confirm no prior diffusion-distance attempt exists**; if one does, this section states exactly how ρ^diff differs. |

What it does **not** replace: the localization atlas (β→OFC), the per-node trace
decomposition, or ρ_sym itself unless/until it demonstrably beats it under
matched-strength and out-resolves both baselines.

---

## 10. Implementation plan

*Reuse-first; ⏳ paths reconciled against the 2026-07-06 grounding sweep.*

- **New library primitive** (general graph concept ⇒ belongs in the library, not a
  script; general names per `feedback_library_names_general`):
  `lrg_eegfc.utils.lrg.diffusion` (or `utils/metrics/diffusion.py`) —
  `heat_kernel_gram(eigvals, eigvecs, tau)`, `diffusion_distance(...)`,
  `n_eff(eigvals, tau)`. Justified promotion: ≥2 callers (observable + MS null +
  future node-resolved pass).
- **Reuse (confirmed 2026-07-06, do not reinvent):**
  `workflow.fc.load_fc_matrix(p, φ, b, "imcoh_abs")`;
  cached eigenpairs `workflow.lrg.load_lrg_result(...).eigenvalues/.eigenvectors`
  (fallback `eigh(diag(W.sum(1))−W)`; the same combinatorial L̂ as
  `metrics/functional_tree_distance.py:284 _laplacian` and
  `surrogate/matched_strength.py:155 _laplacian_eig`);
  the ρ_split template `audit_67_raw_fc_matched_strength.py:199 rho_split_from_phases`
  and the ρ_sym rows `audit_156...:125 _rho_sym_rows` (audit-only — promote to
  `utils/metrics/cross_phase.py` on second use);
  MS surrogates `surrogate/matched_strength.py:220 load_or_compute_surrogate_eigs`
  / `:352 coupled_surrogate_cophenet`; upper tail `hypothesis.py:144 surrogate_p_value`;
  `hypothesis.py:{39 wilcoxon_z, 100 bh_fdr, 78 boot_ci_mean, 197 loo_sensitivity}`;
  tie-corrected `scipy.stats.spearmanr`. `config.const.{PATIENTS_4PHASE,
  PHASE_LABELS, BRAIN_BANDS_NAMES}` (Pat_10 [53,54,55] drop retired — vendor
  re-supplied uniform 113-ch; `load_fc_matrix` returns consistent N, assert it).
- **Script (exploration, does not touch core):**
  `scripts/01_compute/diagnostics/diag_diffusion_distance_trace.py`.
  Outputs to non-preprint locations only: CSV `data/reports/diffusion_distance_trace/`,
  figures `data/outputs/figures/diffusion_trace/`. Nothing under `.agents/preprint/`.
- **Compute discipline** (`feedback_optimize_time_and_surface_progress`): eigh once
  per cell; vectorised einsum for `G_τ`/`D_τ²`; time 1–2 cells and extrapolate the
  full runtime before the full launch; live `[i/N] label elapsed ETA` with
  `flush=True`. Observable first (no claim); MS null only on intermediate-τ cells
  that show structure.

---

## 11. Open questions (load-bearing, deferred)

1. **Primary scale axis:** `s = τλ_2` vs the `n_eff` resolution axis vs `τλ_max`.
   `n_eff` is the most interpretable ("modes resolved"); decide after profiles.
2. **Laplacian:** fixed — the combinatorial Villegas operator (not an open
   question). The strength-sensitivity of `D_τ` (§4) is addressed by the
   matched-strength null, not by normalising the operator (which would leave the
   LRG framework). Locked 2026-07-07.
3. **Aggregation:** cohort median+IQR + per-τ Wilcoxon (keep the profile) — no
   τ-integrated consensus scalar (forbidden). Confirm this is the reporting unit.
4. **Localization (pass 2):** node-resolved diffusion-coordinate displacement
   `‖Ψ_τ^{R⁺}(i) − Ψ_τ^{R⁻}(i)‖` vs task, for the anatomy story. Deferred until the
   pairwise object earns its keep.
5. **Optional synthetic sanity check:** the user chose empirical-only for
   irreducibility; a minimal synthetic construction (eigenvalue-spacing rescale at
   fixed subspace → propagator sees it, Grassmann/adjacency do not) would cheaply
   de-risk the wedge. Flagged, not committed.
