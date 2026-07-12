---
name: 2026-07-12_drift-orthogonalized-persistence
type: guide
era: IMCOH_ABS × COHORT_N10 (post-sparsified-arc reframe)
status: current
created: 2026-07-12
updated: 2026-07-12
pointers:
  - .agents/reports/2026-07-12_brutal-review-failure-directions.md
  - .agents/guides/task-persistence-investigation/README.md
---

# Drift-orthogonalised persistence — a drift-immune trace that keeps the monotonic signal

**Head.** The cross-phase design is `rest_pre → task → rest_post`, intrinsically
time-ordered, so a *trace* (task changes the network AND it persists) is **monotonic
by construction** — and so is spontaneous drift. Any test that removes monotonic-in-time
structure (the linear detrend, `09_drift_detrend.py`) therefore removes the signal by
construction and is **conceptually void** (null = alternative). This measure fixes that:
instead of removing the *direction of change*, it removes the **spontaneous resting
repertoire** (the directions the brain explores on its own at rest, which *contain*
drift), and asks whether the task pushed the network **outside** that repertoire and
whether `rest_post` **stayed** outside. Monotonic change is preserved; drift is not.

## 1. Five-point critical preamble
1. **Claim.** The task drives a reorganisation of the FC/cophenetic structure that lies
   **outside the brain's spontaneous resting repertoire**, and that novel component
   **persists** into `rest_post`.
2. **Null.** Post's alignment with the task-specific (novel) direction is no greater than
   a **task-naïve resting segment's** alignment with it (held-out `rest_pre` stand-in for
   `post`) — i.e. the "persistence" is just spontaneous repertoire re-expression /
   under-sampling. Secondary null: matched-strength.
3. **Strongest alternative the null must control.** (a) **Spontaneous drift** — the
   network was sliding anyway; (b) **under-sampled repertoire** — short `rest_pre` makes
   ordinary rest look "novel"; (c) node **strength**.
4. **Does it control them — by mechanism.** (a) Drift lives in the resting repertoire `S`
   (the network drifts at rest), so `P_⊥` removes it — *without* removing the direction of
   change, so a monotonic trace survives (the detrend's fatal flaw is avoided). (b) The
   held-out-rest null re-does the identical projection with a task-naïve segment; if
   under-sampling manufactured novelty, the null shows the same `T_⊥`. (c) matched-strength.
   **Cannot reject:** a task effect that lies *within* the spontaneous repertoire (the task
   merely amplifies resting modes) → reported as a null, not a false negative; drift that
   is strongly **non-stationary** (rotates between `rest_pre` and task/post).
5. **Falsify.** Cohort `T_⊥ ≤ 0`, or `T_⊥` not exceeding the held-out-rest null, or
   `‖v_⊥‖ ≈ 0` (the task never left the repertoire) → no drift-immune persistent trace.

## 2. Notation
- Node-pairs `m = N(N−1)/2`. State vector `x ∈ ℝ^m` per phase — **primary space =
  cophenetic distance vector** `D(τ)` (trace's native space); FC-edge space is the
  alternative (Open Q1).
- `rest_pre` windows `k = 1..K` → states `x_pre^(k)`; mean `x̄_pre`.
- Centered rest matrix `M = [x_pre^(1) − x̄_pre, …, x_pre^(K) − x̄_pre]ᵀ ∈ ℝ^{K×m}`.
- Task state `x_task` (= `task_test`; optionally shared with `task_learn`, §6).
- Post state `x_post`.

## 3. Predicates & formulas
**Resting repertoire (drift + spontaneous fluctuation).**
```
S = span of top-r right singular vectors of M            (PCA of the rest_pre trajectory)
P_S = V_r V_rᵀ ,   P_⊥ = I − P_S                          (V_r = m×r basis of S)
```
`r` chosen by cross-validated variance (Open Q2); `S` captures drift because the network
drifts within rest.

**Task-specific (novel) direction and its persistence.**
```
v    = x_task − x̄_pre                         # task change
v_⊥  = P_⊥ v ,   v̂_⊥ = v_⊥ / ‖v_⊥‖           # part outside the resting repertoire
d    = x_post − x̄_pre                         # post change
T_⊥  = ⟨ d , v̂_⊥ ⟩                            # does post retain the novel task direction?
```
Report **`T_⊥`** (signed persistence, sign convention: `>0` = trace, per locked T_d rule)
and the **novelty magnitude** `‖v_⊥‖ / ‖v‖` (how much of the task change is non-spontaneous).

**Held-out-rest null (primary).** For each held-out window `k*`: rebuild `S` and `v̂_⊥`
from the remaining data, then set `x_post ← x_pre^(k*)` (task-naïve) and compute
`T_⊥^null(k*)`. The null distribution `{T_⊥^null(k*)}` is "how much a resting segment aligns
with the task direction by repertoire/under-sampling alone." Cohort gate:
`Wilcoxon(T_⊥^obs − median_k* T_⊥^null, greater)`.

## 4. Properties
- **Drift-immune, monotonicity-preserving.** `drift ⊂ S ⇒ P_⊥` kills it; the direction of
  change is *not* removed ⇒ a monotonic trace survives. (The exact repair of the detrend.)
- **Task-within-repertoire ⇒ `‖v_⊥‖ → 0`**, reported as "no novel structure," not a false
  negative.
- **Reference frame is the rest itself** — no external drift model assumed; "drift" is
  defined empirically as what the resting brain does.
- **Reduces to the naïve trace** as `r → 0` (`P_⊥ = I`): `T_⊥ = ⟨x_post−x̄_pre, v̂⟩`, a
  drift-blind cosine trace. `r` is the drift-immunity knob.

## 5. Caveats (must live in the writeup's first paragraph)
- **C-a Under-sampling.** Short `rest_pre` ⇒ small/biased `S` ⇒ inflated novelty. *Handled
  by the held-out-rest null*, which shares the bias. Do NOT interpret `T_⊥>0` without it.
- **C-b `r` sensitivity.** Too small → drift leaks in (false +); too big → `S` eats `v` (false −).
  Report `T_⊥(r)` as a curve; pick `r` by CV, not by the answer.
- **C-c Non-stationary drift.** `S` from early rest may miss late drift ⊥ to it. Diagnostic:
  re-estimate a late-rest drift direction (from `rest_post` windows) and check `v̂_⊥`'s overlap.
- **C-d Nonlinearity.** In cophenetic-distance space the linear subspace is a linearisation;
  a linear FC drift is not linear in `D`. Document; cross-check in FC-edge space (Open Q1).
- **C-e Scope.** Proves **non-spontaneous persistence**, NOT task-content specificity
  (replay). Headline accordingly.

## 6. Using task_learn (optional sharpening — NOT a trajectory point)
`learn`=encoding, `test`=inference are *different computations*; do not treat
`pre→learn→test→post` as one shape (that is the ill-posed "plateau" test). Instead:
- **Robust task direction:** `v = ½[(x_learn − x̄_pre) + (x_test − x̄_pre)]` — the reorg
  common to both task episodes is less likely to be a one-off.
- **Inference-specific variant:** `v_inf = x_test − x_learn`, orthogonalised against `S`,
  persistence into post = inference-specific drift-immune trace (ties to the enc-vs-inf thread).

## 7. Pseudocode
```
for (patient, band):
    x_pre_k   = cophenetic_vec(rest_pre window k)   for k=1..K
    xbar_pre  = mean_k x_pre_k
    M         = stack_k (x_pre_k - xbar_pre)
    V_r       = top-r right singular vectors of M          # repertoire basis
    Pperp     = I - V_r V_rᵀ
    v         = cophenetic_vec(task) - xbar_pre
    vperp     = Pperp @ v ;  vhat = vperp / norm(vperp)
    d         = cophenetic_vec(post) - xbar_pre
    T_perp    = dot(d, vhat)
    novelty   = norm(vperp) / norm(v)
    # held-out-rest null
    for k* in 1..K:
        rebuild V_r, vhat from data\{k*}
        T_null[k*] = dot(x_pre_{k*} - xbar_pre(\k*), vhat)
cohort: Wilcoxon(T_perp_obs - median(T_null), greater) ; also sweep r ; matched-strength 2ndary
```
Cost: cophenetic per rest-window + task + post (reuse the 1-Welch-bank-per-window /
all-bands trick from `08_drift_null_fair.py`). Cheap; deterministic except MS.

## 8. Connection to prior tools
- **ρ_sym** (`heat_multiscale.rho_sym`): the naïve drift-blind trace; this is its
  `r→0`-reduces-to, drift-immune generalisation. Same cophenetic substrate.
- **Detrend** (`09_drift_detrend.py`): the void test this replaces — it removed the change
  direction; this removes only the resting repertoire.
- **State-space / attractor** (audit_164) and **reinstatement** (N1/N4): the "resting
  repertoire" is the resting manifold; `v_⊥` = leaving it. Complementary framings.
- **Brutal review** `2026-07-12_brutal-review-failure-directions.md`: this addresses C6
  (nulls) for the trace; it does NOT address C3/C8 (is the multiscale/backbone even
  non-trivial) — run T1 (trivial-multiscale) independently.

## 9. Open questions (decide before coding)
1. **Space:** cophenetic-distance (native, nonlinear) vs FC-edge (drift more linear).
   Recommend: run BOTH; FC-edge as the drift-clean check, cophenetic as the trace claim.
2. **r selection:** variance-explained threshold? CV on held-out rest reconstruction?
   Report the full `T_⊥(r)` curve regardless — never pick `r` by the p-value.
3. **Global vs per-node:** start global (cohort trace), then per-node `v_⊥` for localisation
   (which nodes carry the non-spontaneous persistent structure).
4. **Which null is primary** — held-out-rest (controls under-sampling) vs matched-strength
   (controls strength). They control different things; likely need both.
