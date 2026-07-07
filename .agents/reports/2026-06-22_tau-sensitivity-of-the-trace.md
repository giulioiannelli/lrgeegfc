---
name: tau-sensitivity-of-the-trace
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-06-22
updated: 2026-06-22
pointers:
  - .agents/guides/task-persistence-investigation/2026-06-22_tau-sensitivity-cophenetic-trace.md
  - scripts/01_compute/audit/audit_121_tau_sweep_cophenetic_trace.py
  - scripts/01_compute/audit/audit_121b_tau_sweep_figure.py
  - data/audit/tau_sweep_trace/audit_121_tau_sweep_trace.pdf
  - .agents/guides/02_methods/lrg-framework-guide.md
---

> **Head.** We swept the LRG diffusion time τ — fixed at `1/λ_max` in every
> result so far — across the full meaningful spectral window (finest →
> Fiedler → collapse) and re-measured the cohort cross-phase trace at each τ.
> **Verdict: `τ = 1/λ_max` is vindicated, not refuted.** The genuine trace is a
> **fine-scale** phenomenon that is **τ-robust** (flat) from the finest scale out
> to the Fiedler time, then dissolves as the kernel collapses. **No coarser τ
> delivers a stronger, artifact-clean trace.** The SOZ precedent — where *slow*
> diffusion revealed distant seizure-onset structure — **does NOT transfer** to
> the cross-phase trace; the two objects respond to τ in opposite directions.
> The investigation paid for itself twice: it **discharges** the standing
> "sensitivity to τ not exhaustively tested" caveat that flagged every cohort
> claim, and it **caught and quarantined** a coarse-τ `ρ_split` inflation that
> could have become a false "coarse-τ strengthens the trace" headline (it is a
> collapse-regime artifact, placebo-confirmed). **One genuine lead survives:
> the α-band trace prefers slightly coarser τ (α ≈ 2–3), roughly doubling — the
> only cell where `1/λ_max` may leave signal on the table. It needs a
> matched-strength null to become a claim.** Observed-statistic only; no
> matched-strength null was run here.

---

## 1. Why this was worth doing

`τ = 1/λ_max` is the *finest/fastest* diffusion time, borrowed from the LRG
papers where it is "the right place to start" for **sparse topological networks
with discrete spectra**. Our graphs are **fully connected, weight-heterogeneous,
near-continuous-spectrum** — the regime the papers explicitly exclude
(guide §5.2, §6). At `1/λ_max` the heat kernel has barely diffused, so
`D_coph(1/λ_max)` ≈ a monotone reweighting of the raw FC; the multi-step /
mesoscale structure that is LRG's value-add only develops at **larger τ**. And
the in-house SOZ marker had *already shown* that τ matters and that the
informative regime there is **slow** diffusion (δ distant-SOZ at `τ ≈ 4–10/λ_max`;
audit_85/101/102). So the worry was concrete: have all the trace verdicts been
read off one — possibly suboptimal — corner of τ-space? guide open-question #1.

## 2. What we did (audit_121)

τ swept **per phase** as `τ_s = α / λ_max^s`, so `α = 1` reconstructs each
phase's own `1/λ_max` and reproduces the **cached cophenetic matrix bit-exactly**
(verified, max |Δ| = 2.8 × 10⁻¹⁴ on all 60 cells — the anchor). `D_coph(τ)`
reconstructs from cached Laplacian eigenpairs (`V e^{−τλ} Vᵀ/Z` → UPGMA →
cophenetic), so the full 10 × 6 × 29-α sweep is cheap. At each α we recomputed,
on the cophenetic geometry: the locked triangle `T_d = ρ^coph(task,post) −
ρ^coph(pre,task)`, the split-baseline `ρ_split = Spearman(d_task−d_preA,
d_post−d_preB)`, a **trace-free placebo** `ρ_indep = Spearman(d_preA−d_preB,
d_task−d_post)` (must stay ≈0 unless coarsening manufactures correlation),
self-similarity to the canonical geometry, and a collapse monitor.
α ∈ [0.5, 40] brackets finer-than-current → genuine LRG mesoscale → post-Fiedler
collapse.

## 3. Findings

**(F1) The trace is fine-scale and τ-robust; `1/λ_max` is near-optimal.**
- β triangle `T_d`: positive plateau over α ∈ [0.5, 1.8], peak +0.093 @ α=0.70 vs
  +0.088 @ α=1 (Δ = +0.004 — flat), decays to ≈0 by α≈2. **0 %** of coarse-τ
  beats α=1.
- β `ρ_split`: **flat at ≈ +0.28** across the whole clean window (α=1 sits at
  +0.222, ~0.04 below a tightly-clustered +0.24–0.30 plateau). Robust, not a
  knife-edge at `1/λ_max`.
- For `T_d`, α=1 is inside or above the clean-window IQR for every band → a
  sound choice.

**(F2) The "coarse-τ strengthens the trace" appearance is a collapse artifact —
placebo-confirmed.** As τ coarsens past Fiedler the kernel goes low-rank (the
slow anatomy/Fiedler mode dominates all phases), `D_coph` → near-constant
(self-similarity 1.0 → 0.5; `frac_collapsed` 0 → 0.9), and the **placebo
`ρ_indep` rises from ≈0.05 to ≈0.23 to MATCH `ρ_split`**. The bands that
"lit up" at coarse τ in a naïve argmax (δ +0.03→+0.23, θ, γ_hi swinging +0.7→−0.4)
are the **anatomy/epi channels**, and their coarse `ρ_split` sits at or below
their placebo (δ ratio 2.0, θ 4.6, γ_hi 7.3 → not genuine trace). Gating on the
placebo + collapse confines every trustworthy reading to α ≲ 10–15.

**(F3) The SOZ analogy does not transfer.** The δ *node marker* genuinely
improved at slow τ (vs a label-shuffle null) because slow global diffusion
reaches distant co-diffusing nodes. The δ *cross-phase trace* does the opposite —
it has no fine-scale signal and its coarse-τ appearance is artifact. Different
questions, opposite τ-response. Assuming the analogy would have been wrong.

**(F4) One real lead — the α-band prefers slightly coarser τ.** In the clean
window the α-band triangle `T_d` **doubles**: +0.066 @ α=1 → +0.131 @ α=3.26;
`ρ_split` +0.115 → +0.214 @ α=2.09. Both inside the no-collapse / placebo<0.12
window. Placebo margin is thinner here (ratio ≈0.85), so this is **a lead, not a
claim**: it is the one cell where `τ = 1/λ_max` plausibly under-reads, and it is
mechanistically sensible (α reorganization at a marginally coarser network scale
than β's). **Decisive test = matched-strength null rebuilt at α ≈ 2–3.**

## 4. Caveats (read first)

- **Observed-statistic only.** No matched-strength surrogate was run. F1–F3 are
  robustness / negative findings about the observed statistic's τ-profile and an
  internal placebo, which do not themselves require the matched-strength null;
  **F4 (a positive lead) does**, per the mandatory-null rule.
- The placebo is a *cheap* stand-in for "does coarsening manufacture
  correlation"; it is not the matched-strength null.
- This speaks to the **cohort scalar** trace only. β→OFC per-node localization
  τ-sensitivity is a separate, heavier follow-up.

## 5. Next steps (for PI decision — none run yet)

1. **Matched-strength null at α ≈ 2–3 for the α-band** (settle F4). Smallest,
   highest-value follow-up.
2. Optionally confirm β robustness with a matched-strength null at 2–3 α points
   across the clean plateau (turn F1 into a hardened "τ-robust" statement).
3. Re-localization (β→OFC) at α=1 vs the α-band's preferred τ — only if F4 holds.

Figure: `data/audit/tau_sweep_trace/audit_121_tau_sweep_trace.pdf`
(a = β headline; b/c = ρ_split / T_d all bands; d = self-similarity + collapse).
