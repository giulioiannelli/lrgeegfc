---
name: tau-sweep-consolidation-arc
type: scope
era: IMCOH_ABS × COHORT_N10
status: executed_2026-06-22
created: 2026-06-22
result: .agents/reports/2026-06-22_tau-sweep-arc-result.md
extends: .agents/guides/task-persistence-investigation/2026-06-22_tau-sensitivity-cophenetic-trace.md
pointers:
  - scripts/01_compute/audit/audit_103_cophenetic_consolidation_arc.py
  - scripts/01_compute/audit/audit_110_inference_mark_localization.py
  - .agents/preprint/headlines/02_encoding_vs_inference.md
  - .agents/guides/02_methods/lrg-framework-guide.md
---

> **Head.** The whole encoding-vs-inference decomposition (N2) — the encoding echo
> `T_learn`, the inference-specific component `T_infspec·e`, the α/β dissociation,
> the encoding→OFC localization — was computed at the **single** diffusion time
> `τ = 1/λ_max` (the finest/fastest scale). At that scale the heat kernel has barely
> diffused, so the cophenetic distances are close to a reweighting of the raw FC; the
> multi-step / mesoscale structure that is LRG's actual value-add only develops at
> **larger τ**. This scope asks whether the *content decomposition* is robust across
> scale, and — more interestingly — whether **inference-specific consolidation is a
> fine-scale or a mesoscale phenomenon**. It is the arc-specific complement to the
> sibling τ-sensitivity scope (which covers the N1 *trace*); read that first for the
> general τ machinery and the "1/λ_max is borrowed from the sparse-spectrum regime"
> argument — not repeated here.

---

## 0. Five-point critical preamble (mandatory; written before any code)

**(1) Claim under interrogation.** Not a new positive result. The object is the
**scale-robustness of the N2 decomposition**: "the encoding echo (α & β), the
β-specific inference-specific component, and the encoding→OFC localization are
properties of the LRG geometry, not artifacts of evaluating it at the single finest
scale `τ = 1/λ_max`." A *secondary, exploratory* question rides along: does the
inference-specific component **strengthen at slower diffusion** (the mesoscale
regime), which would sharpen the cognitive interpretation (inference = integration
over multi-step relational paths)?

**(2) Reference / null.** The reference is each functional's value **at
`τ = 1/λ_max`** (current canonical; the diagnostic must reproduce it bit-exactly —
correctness anchor against `audit_103`). The sweep measures the **observed**
functionals `T_learn(τ)`, `T_infspec·e(τ)`, and the encoding/inference localization
as functions of τ over the meaningful spectral window. First-order question is
binary: **flat** (τ-irrelevant → N2 verdicts robust) or **structured** (τ matters →
N2 is scale-specific). **No matched-strength surrogate is run in the sweep itself**
— it is an observed-statistic diagnostic. A surrogate null is **mandatory** only if
the sweep surfaces a *new* positive claim at some τ ≠ 1/λ_max (e.g. "inference peaks
at τ ≈ 4/λ_max"), which would trigger a dedicated follow-up at that τ.

**(3) Strongest alternative the diagnostic must control for.** That any apparent
τ-structure is a **numerical artifact** of the cophenetic construction near the
spectral extremes (e.g. `1/ρ` blow-ups as the density matrix approaches uniformity
at large τ, or degenerate ties at tiny τ) rather than a real change in the
consolidation signal. Controlled by: (a) restricting τ to the validated window from
the sibling scope; (b) the bit-exact reproduction at `1/λ_max`; (c) reporting the
**raw cophenetic stability** (fraction of finite, non-tied entries) alongside each
τ so a degenerate regime is visible, not hidden.

**(4) Does it address the duration confound?** Orthogonal — τ is a graph-geometry
axis, the duration confound is a data-amount axis (already resolved: N2.5, β
duration-robust). The sweep neither re-opens nor depends on that. *But* one useful
cross-check: the encoding functional `e = D_taskLearn − D_preA` is duration-immune at
every τ, so `T_learn(τ)` is a clean reference curve uncontaminated by the
test/learn length asymmetry at all scales.

**(5) Falsification / limits.** If `T_infspec·e(τ)` is **flat and β-only across the
whole window**, N2.2 graduates from "robust at one scale" to "scale-robust" — a
strengthening, no new claim. If it is **structured** (e.g. inference emerges only at
larger τ, or the α/β dissociation closes/opens with scale), N2's framing must say so
explicitly and the canonical `1/λ_max` is re-justified or replaced. **Cannot do:**
prove a *causal* multi-step-integration interpretation of any τ-peak (that is a
correlational scale signature, not a mechanism); rescue the retired whole-brain
truncation null (a τ-sweep does not fix the truncation artifact).

---

## 1. Notation

Per phase x and diffusion time τ: `ρ̂_x(τ) = e^{−τ L̂_x} / Z`, ultrametric distance
`T_x(τ) = 1/ρ̂_x(τ)`, average-linkage cophenetic condensed vector `D_x(τ)`. The
four-phase arc functionals, now τ-resolved (signs per the locked convention,
positive = trace):

- encoding `e(τ) = D_taskLearn(τ) − D_preA(τ)`
- inference-specific `f(τ) = D_taskTest(τ) − D_taskLearn(τ)`
- persistence `p(τ) = D_restPost(τ) − D_preB(τ)`
- `T_learn(τ) = ρ_S(e(τ), p(τ))`  (encoding echo)
- `T_infspec·e(τ) = partial ρ_S(f(τ), p(τ) | e(τ))`  (inference, controlling encoding)

τ-grid: the **log-spaced spectral window** from the sibling scope (anchored so one
node is exactly `1/λ_max`), spanning `≈ [1/λ_max, ~10/λ_max]` — the SOZ-marker work
showed the informative slow-diffusion regime lives around `4–10/λ_max` in this exact
dataset.

## 2. Predicates (what we read off the curves)

- **P1 robustness:** `T_infspec·e(τ)` significant-sign-stable and **β-only** across
  the window ⇒ N2.2 scale-robust.
- **P2 scale-of-inference:** location of `argmax_τ T_infspec·e(τ)`. Fine
  (`≈1/λ_max`) vs mesoscale (`≫1/λ_max`) is the cognitively interesting bit.
- **P3 dissociation stability:** does `T_learn(τ)` stay α-and-β while
  `T_infspec·e(τ)` stays β-only across τ, or do they converge/cross?
- **P4 localization drift:** recompute the encoding→OFC / inference→(distributed)
  localization at 2–3 representative τ; does OFC sharpen, hold, or dissolve as the
  kernel coarse-grains? (Localization needs the matched-strength null → only at the
  handful of τ that P1–P3 flag as interesting, not the whole grid.)

## 3. Properties / correctness anchors

- **C-exact:** at the `1/λ_max` node, `T_learn`, `T_infspec·e` reproduce
  `audit_103` per-patient to ~1e-12 (same anchor audit_103 already uses vs audit_83).
- **C-monotone-floor:** report finite/non-tied fraction of each `D_x(τ)` per τ;
  flag any τ where it drops (degenerate regime → exclude from interpretation).
- **C-duration-immune:** `T_learn(τ)` is independent of `task_test` length at every
  τ (no `task_test` term) — a built-in clean reference.

## 4. Pseudocode

```
for pat in COHORT:
  Ws = {ph: load_phase_fc(pat, ph, band) for ph in ARC_PHASES}   # reuse audit_103
  for tau in TAU_GRID:                      # tau as multiple of 1/lambda_max
    D = {ph: cophenetic_condensed_at_tau(Ws[ph], tau) for ph in ARC_PHASES}
    record T_learn(tau), T_infspec_pe(tau)  # reuse audit_103 _obs_functionals
# cohort: per band, per tau -> median + one-sided sign test (observed only)
# localization (P4): only at flagged tau, run audit_110 with cophenet_at_tau
```

Only new primitive needed: `cophenetic_condensed_at_tau(W, tau)` = `audit_63`’s
`lrg_ultrametric_condensed` with the hardcoded `tau = 1/lam_max` lifted to a
parameter. Everything else (loaders, functionals, localizer) is reused verbatim.

## 5. Visualization

- **Per-band τ-response curves** `T_learn(τ)` and `T_infspec·e(τ)` (α & β
  overlaid), cohort median ± IQR, x = τ·λ_max (log), vertical line at the canonical
  `1`. The single most informative panel: is β-inference flat, or does it climb into
  the mesoscale?
- **Localization-vs-τ** small multiples at the flagged τ (encoding OFC weight).

## 6. Connection to prior tools

- Reuses `audit_103` (arc functionals) + `audit_110` (localizer) + `audit_63`
  ultrametric — only τ is lifted to a parameter.
- Complements the sibling **N1** τ-sensitivity scope (`2026-06-22_tau-sensitivity-
  cophenetic-trace.md`): that one sweeps the *single-phase cross-phase trace*; this
  one sweeps the *four-phase content decomposition*. Share the τ-grid and the
  `cophenetic_*_at_tau` primitive.
- The SOZ-marker τ work (audit_85/101/102) is the empirical precedent that slow
  diffusion is the informative regime here.

## 7. Open questions

- Is inference-specific consolidation **fine-scale or mesoscale**? (P2 — the
  scientifically loaded one; a mesoscale peak would support an
  "integration-over-relational-paths" reading.)
- Does the encoding→OFC localization **sharpen** at coarser τ (mesoscale
  communities) or is it a fine-scale effect?
- If `1/λ_max` turns out suboptimal for inference, does the **whole N2 headline**
  re-anchor to a different canonical τ — and does that change which band carries
  inference?
