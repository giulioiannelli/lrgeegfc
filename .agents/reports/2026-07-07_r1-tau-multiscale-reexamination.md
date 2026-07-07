---
name: r1-tau-multiscale-reexamination
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-07
decision: PARKED — crystallized, NOT in the paper
supersedes_claims: []
pointers:
  - scripts/01_compute/audit/audit_121_tau_sweep_cophenetic_trace.py
  - scripts/01_compute/audit/audit_157_arc_scale_null_rhosym.py
  - scripts/01_compute/audit/audit_150_rho_sym_gate.py
  - data/audit/tau_sweep_trace/cohort.csv
  - .agents/guides/task-persistence-investigation/2026-06-22_tau-sensitivity-cophenetic-trace.md
memory: tau_sensitivity_trace_2026_06_22
---

# R1 τ-multiscale re-examination (ρ_sym era) — verdict: narrow prize, parked out of the paper

## Head

R2 varies the LRG diffusion time τ to find its mesoscale scale-signature, so we
asked the same of R1's **bare** trace: does reading it across τ, per band, add a
multiscale lens worth putting in the paper? **Answer: no — the prize is real but
narrow, and it does not clear the bar for a headline.** The τ-axis was already
swept for the bare trace (`audit_121`, 2026-06-22) and the honest gain on top of
that is (i) one *marginal* candidate — an **α mesoscale peak** that would need a
ρ_sym re-estimate **and** a matched-strength null it has never had, and (ii) a
*reframe* of β from "τ-robust" to "scale-broad" that mostly serves to sharpen
**R2**, not R1. β stays fixed at τ = 1/λ_max; γ_low is **not** rescued at any τ.
Decision: crystallize here, keep τ = 1/λ_max in R1, do not run the follow-up
(`audit_159`) unless the α candidate is later wanted as a standalone note.

## 1. Why this came up again

This is the **second** time a τ-sweep of the R1 trace has been raised and closed.
The motivation is legitimate: R2's flagship uses τ as a *variation* axis and finds
the inference-specific β component is **mesoscale-favouring** (`audit_157`:
`T_infspec_pe` +0.091 → **+0.126** at τ≈2.6, matched-strength Wilcoxon p 0.0098 →
0.0049). If the *decomposed* component has a scale signature, does the *bare* trace
have one too, per band — enough to add a genuine multiscale lens to R1?

## 2. What audit_121 already established (real numbers, ρ_split, observed-only)

`data/audit/tau_sweep_trace/cohort.csv` — median per-pair split-baseline
concordance ρ_split^coph vs α = τ·λ_max (α=1 = current convention; Fiedler
collapse sets in past α ≈ 2–4 per band; the trace-free placebo ρ_indep is the
collapse guard):

| band | fine α≈0.5 | **α=1 (now)** | **meso α≈1.5–2** | past Fiedler | read |
|---|---|---|---|---|---|
| **β** | 0.33 | 0.22 | 0.29 | flat ~0.28 → collapse@α40 | **scale-broad, τ-robust** (α=1 is a slight local *min*; finest is strongest) |
| **α** | 0.10 | 0.115 | **0.22 (≈doubles)** | decays | **mesoscale peak**, fcoll=0, placebo ~0.07 (margin ~3×) |
| γ_low | 0.13 | 0.14 | 0.13 | flat → collapse | flat; fails the cohort gate at **every** τ |
| δ / θ / γ_high | ~0 / neg | ~0 / neg | ~0 | only rise *into* collapse (fcoll↑) | **coarse-τ artifacts** (placebo catches up) |

The 2026-06-22 conclusion stands: **τ = 1/λ_max is vindicated**, β is scale-broad,
the α≈2–3 lead is present but was judged *marginal* (PI declined follow-up). The
coarse-τ "gains" in δ/θ/γ_high are placebo-confirmed collapse artifacts.

## 3. Why audit_121 is not, by itself, the final word

Two gaps kept audit_121 from being able to *certify* any τ≠1 structure:

1. **Estimator.** audit_121 uses **ρ_split** with the arbitrary A/B split-half arm
   assignment — retired on 2026-07-06 in favour of **ρ_sym** (`audit_149`/`150`).
   The α mesoscale peak lives precisely in the split-baseline differences of
   near-floor patients, the cells ρ_sym moves most. (Reassurance: α's *τ=1* signal
   is ρ_sym-robust — 0.115 → 0.10, gate p=0.024 — so the anchor survives; whether
   the *doubling* survives is untested.)
2. **Null.** audit_121 carries only the placebo guard (ρ_indep), never the
   mandatory matched-strength null at any τ≠1. The placebo is a *weaker* referee
   than matched-strength.

## 4. What a proper τ-lens (ρ_sym + matched-strength across τ) would add — and its size

- **① α gains a *scale* though it has no *place* (the one real candidate).** R1 says
  α is "real but spatially scattered, no home." A certified α mesoscale peak would
  make that a positive fact: **α's trace has a characteristic scale (mesoscale)
  though no characteristic address.** Plausible — matched-strength barely dented the
  β *arc* at the mesoscale (surr med ~0.01) — but **not guaranteed**, and it is a
  band that already sits near the measurement floor for half the cohort. This is the
  only piece that could be *new*.
- **② β reframe → sharpens R2, not R1.** "β trace is scale-broad" is a cleaner
  statement than "τ-robust / not a cherry-pick," and it **dissociates** from R2's
  mesoscale inference component (*the trace persists at every scale; the part the
  brain computed concentrates at the mesoscale*). But this strengthens the **N1→N2
  arc**, i.e. it is really an R2 talking point; it adds no new R1 result.
- **③ Documented negatives.** γ_low is not rescued at any τ; δ/θ/γ_high coarse rises
  are collapse artifacts. Useful housekeeping (prevents a future false coarse-τ
  headline), not a result.

Net gain for **R1 as a paper section**: one marginal, uncertain candidate (α) plus
a reframe that mostly helps R2. The multiscale story the paper needs is **already
carried** — by R2's inference-mesoscale signature and by the fact that the
cophenetic distance is itself a full-hierarchy (multiscale) summary. That is why
the bar is not cleared.

## 5. Decision

- **Keep τ = 1/λ_max in R1.** N1.5's "τ-robust, fine-scale" stands as written; no
  change to the headline, no τ figure/supplement in the paper.
- **Do not run `audit_159`** (the bare-trace ρ_sym × τ × matched-strength follow-up)
  as part of the paper program. It is scoped below and cheap, but its only
  certifiable upside (α mesoscale) is marginal and its main narrative payoff (β
  scale-broad vs R2 inference-mesoscale) is an R2 point already made without it.
- **This record is the crystallization.** Nothing enters the preprint.

## 6. If ever revived — the exact, cheap recipe

The engine already exists; this is a small stitch, not new machinery:

- Feed `audit_150`'s **ρ_sym gate statistic** (the N1.1 symmetric split-baseline
  concordance) into `audit_157`'s **τ-null loop** (cached seed-20260511 ensemble,
  cophenetic reformed at each τ from cached eigenpairs, **no regeneration**), over a
  **finer α-grid through the α≈1.5–2 window**. Anchor bit-exact to `audit_150` at
  α=1.
- **Scope-report-first** (mandatory) in `.agents/guides/task-persistence-investigation/`
  with the 5-point critical preamble:
  1. *Claim* — bare β/α trace has a per-band scale signature (β broad, α mesoscale)
     under ρ_sym.
  2. *Null* — canonical matched-strength (4-cycle ±δ, seed 20260511, R=200),
     cophenetic reformed at each τ.
  3. *Strongest alternative* — a mesoscale bump is generic coarse-graining geometry
     of any strength-heterogeneous graph; the null bumps identically → the obs−surr
     gap stays flat. α-vs-neg-control-bands is the discriminator.
  4. *Neg-control discipline* — δ/θ/γ_high are dead-null at τ=1; if any "passes" at
     the mesoscale the surface is manufacturing significance → retire (the
     `audit_103b` lesson). Never read past collapse onset (fcoll>0, α≳10–15).
  5. *Falsification* — PASS = α clears matched-strength at α≈1.5–2 (LO-P15 robust)
     while controls stay null. FAIL = α doesn't clear, or a control does → α
     mesoscale is observed-only; headline stays at τ=1. β is locked at τ=1
     regardless.

## Provenance

- τ-sweep (bare trace, ρ_split, observed-only): `audit_121_tau_sweep_cophenetic_trace.py`
  → `data/audit/tau_sweep_trace/{cohort,per_patient}.csv` · 2026-06-22. Scope:
  `.agents/guides/task-persistence-investigation/2026-06-22_tau-sensitivity-cophenetic-trace.md`.
  Memory: `tau_sensitivity_trace_2026_06_22`.
- R2 τ-null engine reused as the template: `audit_157_arc_scale_null_rhosym.py`
  (ρ_sym × τ × matched-strength) → `data/audit/consolidation_arc_rhosym/arc_scale_cohort_verdict.csv`
  · 2026-07-07. Report: `.agents/reports/2026-07-07_r2-mesoscale-encoding-rhosym-migration.md`.
- ρ_sym gate statistic (α=1 anchor for any revival): `audit_150_rho_sym_gate.py`
  → `data/audit/rho_sym_gate/` · 2026-07-06.
