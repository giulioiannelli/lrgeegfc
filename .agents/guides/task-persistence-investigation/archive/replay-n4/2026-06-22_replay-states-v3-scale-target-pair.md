---
name: replay-states-v3-scale-target-pair
era: IMCOH_ABS_COHORT_N10
status: CLOSED_negative_transient_sustained_is_headline
kind: scope-report
measure: replay states v3 — open the THREE axes prior tests collapsed (diffusion scale τ, replay target learn/test, per-pair/subtree)
headline: N4 (.agents/preprint/headlines/04_replay_states.md)
supersedes_framing: v1 transient-flash (audit_125/126/127 NEG) + v2 sustained/short-window/γ (audit_129/130/131 NEG-transient, POS-sustained)
prior_verdict: .agents/reports/2026-06-22_replay-states-verification.md
created: 2026-06-22
---

# Replay states v3 — scale / target / pair-resolved transient search

**Head.** Every prior replay test (audit_124–132) sampled ONE corner of a 4-D tensor
`T(scale τ, pair-set, target phase, time)`: finest scale τ=1/λmax, all-pairs pooled to a
global scalar, target=task_test, and asked "is the time series bursty." Clean cohort
NEGATIVE for transient states; strong POSITIVE for a sustained hold (=N1). v3 opens the
three COLLAPSED axes — **scale τ** (PI hypothesis: the transient may hide at a mesoscale
τ\*), **target** (learn vs test trajectory), **per-pair/subtree** (OFC) — staying strictly
in LRG ρ^coph on |ImCoh| windows. Either it surfaces a scale/target/pair-localized replay
state, or it makes the "sustained-not-transient" verdict scale/target/pair-resolved
(airtight). Constraint: NO new off-framework analysis; signed-Laplacian variants allowed.

## Notation
- Phases P = {R⁻ = rest_pre, L = task_learn, T = task_test, R⁺ = rest_post}.
- Windows: non-overlapping length ℓ = 20 s (the only |ImCoh|-feasible length for α/β,
  audit_128). Per window w: |ImCoh| adjacency `W_w` → `L_w = diag(W_w·1) − W_w` →
  eigendecomp `(λ_w, V_w)`. (Recompute + CACHE windowed eigenpairs once; the cophenetic
  cache from audit_125 stored only τ₀ vectors, not eigenpairs.)
- LRG propagator at diffusion time τ: `ρ̂_w(τ) = V_w exp(−τλ_w) V_wᵀ / Z`,
  `Δ_w(τ) = 1/ρ̂_w(τ)` (max-sym, zero-diag), cophenetic `Dᶜ_w(τ)` = average-linkage
  ultrametric of `Δ_w(τ)`.  **Reuse `audit_121.D_raw` / `D_coph` verbatim** (bit-exact
  anchor at τ=1/λmax). τ-grid: α/λmax, α ∈ [0.5, 40] (audit_121 grid; finest → past-Fiedler).
- Reference configs (full-phase, canonical imcoh_abs): `Dᶜ_X(τ)`, X ∈ {R⁻, L, T}.
- **Scale-resolved task-likeness contrast** (target Y∈{L,T}, baseline R⁻):
  `sᵞ_w(τ) = ρ_coph(Dᶜ_w(τ), Dᶜ_Y(τ)) − ρ_coph(Dᶜ_w(τ), Dᶜ_{R⁻}(τ))`.

## Critical 5-point preamble (per direction)
**R1 τ-resolved transient.** (1) Claim: ∃ scale τ\* where R⁺ is bursty toward T beyond
its stationary mean. (2) Null: time-shuffle of window order (AC1 / excess-flash vs
shuffled marginal) — the make-or-break, identical to audit_125/131. (3) Strongest
alternative: coarse-τ LRG collapse makes Dᶜ_w(τ) unstable window-to-window → spurious
"bursts." (4) Control: time-shuffle PRESERVES collapse-noise (same windows reordered) ⇒
robust to it (unlike the static τ-sweep, audit_121, which needed the ρ_indep placebo);
ADD a dynamic PLACEBO target (other-patient T, or R⁻-split independent diff) that must
NOT burst — if it bursts identically at coarse τ, the effect is collapse, not replay.
(5) Falsify: no τ with cohort burstiness > shuffle that also beats placebo.

**R2 dual-target trajectory.** (1) Claim: replay target drifts over R⁺ (early L/encoding
→ late T/inference) and/or bursts toward T-specific. (2) Null: window-index permutation
for the trend; time-shuffle for bursts. (3) Alternative: monotone arousal drift mimics a
trend in either target equally. (4) Control: test the DIFFERENCE `sᵀ_w − sᴸ_w` (a
common arousal drift cancels in the difference) for trend; placebo = R⁻ trajectory.
(5) Falsify: flat `sᵀ−sᴸ` trend + no T-specific burst. NB audit_127 already nulled the
POOLED sustained inference-preference; R2's NEW content is the TEMPORAL trajectory + the
τ/pair combination.

**R3 per-pair / OFC-subtree transient.** (1) Claim: a localized subtree (OFC, B1 carrier)
transiently snaps to task config though the global average is flat. (2) Null: time-shuffle
on the subtree-restricted task-likeness; per-pair bimodality vs shuffle. (3) Alternative:
N²/2 pairs → multiple-comparison false bursts. (4) Control: restrict to the a-priori OFC
subtree (NOT data-mined) to fix the comparison count; report effect not just p.
(5) Falsify: OFC-subtree task-likeness no more bursty than its time-shuffle.

## Predicates / statistics
- SUSTAINED (known, anchor): `mean_{R⁺} sᵀ_w(τ₀) > mean_{R⁻}` , τ₀=1/λmax → 10/10 p=0.001.
- R1: per τ, `shift(τ)`, `AC1_post(τ)` vs time-shuffle `p_time(τ)`, `excess_flash(τ)`,
  `spread(τ)`, `placebo_flash(τ)`. Cohort one-sided Wilcoxon per τ; report the τ-profile.
- R2: `trend = Spearman(window_index, sᵀ_w − sᴸ_w)` over R⁺ (and R⁻ control); cohort sign
  test. Plus T-specific burst = burstiness of `sᵀ_w − sᴸ_w`.
- R3: OFC-subtree-restricted `s` (endpoint-incidence keep, B1/audit_112 convention);
  burstiness + per-pair switching-rate vs time-shuffle.

## Nulls / controls (ladder)
1. **Time-shuffle** (window-order permutation) — burstiness null; robust to coarse-τ
   collapse. Primary make-or-break.
2. **Dynamic placebo target** — config the window cannot be replaying (other-patient T;
   or R⁻-split independent difference). Coarse-τ collapse manufactures structure for ANY
   target → placebo must stay flat. (Dynamic analog of audit_121 ρ_indep.)
3. **Matched-strength** — ONLY if a positive survives 1 & 2 (mandatory referee for a
   claim, moot for a null). Strength-preserving shuffle of the windowed W_w.

## Feasibility
- α/β feasible only at ℓ=20 s (audit_128) — v3 keeps ℓ=20 s and opens the τ axis, it does
  NOT reopen short windows (that's settled negative). γ at 2 s can ride along (R1 τ-sweep
  on the 2 s γ cube) but γ is not a trace band.
- Cost: per pat/band ~600 windows × eigh(N≈100) + (~15 τ) × linkage(N). ≈ minutes/band;
  cache windowed eigenpairs so R1/R2/R3 share. τ₀ row MUST reproduce cached burstiness
  (anchor).

## Pseudocode (R1)
```
for pat, band:
  refs: for X in {R⁻, L, T}: eig(full_phase W_X) → Dᶜ_X(τ) for τ in grid
  win:  for ph in {R⁻, R⁺}: for w: W_w=windowed_imcoh; (λ,V)=eigh; cache; Dᶜ_w(τ)
  for τ in grid:
     s_pre[τ]  = [ρcoph(Dᶜ_w(τ),Dᶜ_T(τ)) − ρcoph(Dᶜ_w(τ),Dᶜ_R⁻(τ)) for w in R⁻]
     s_post[τ] = [...                                                   for w in R⁺]
     shift, AC1_post, p_time(shuffle), excess_flash, spread
     placebo[τ] via Dᶜ_T → other-patient T
  → per-τ cohort Wilcoxon; flag any τ with burst>shuffle AND >placebo
```

## Visualization
- R1: τ-profile heatmap (band × τ) of burstiness p / excess_flash, placebo overlaid;
  the money panel = a single window's dendrogram at τ\* snapping to task next to baseline.
- R2: per-patient `sᵀ_w − sᴸ_w` trajectory over R⁺ windows (target-drift curve).
- R3: OFC-subtree task-likeness time-course, R⁺ vs R⁻, with bursts marked.

## Connection to prior tools (library-first; v3 is a COMPOSITION, not new analysis)
- τ reconstruction: `audit_121` `D_raw`, `D_coph` (cached-eig, bit-exact anchor).
- windowed |ImCoh|: `lrg_eegfc.utils.fc.coherence.windowed` (segment_ffts, imcoh_abs_cube,
  band_abs_average) + audit_125 windowing.
- burstiness battery: audit_125/130/131 (ac1, perm_p, excess_flash, shift-vs-states).
- OFC subtree: B1 `audit_132` (load_channel_regions, endpoint-incidence keep, epi_keep_mask).
- references/cophenetic: `cophenetic_condensed_from_adjacency`, `load_fc_matrix`.

## Open questions
- τ-grid density for the windowed pass (full 30-pt vs ~15-pt fine→Fiedler).
- Placebo target choice: other-patient T (anatomy-mismatched) vs within-patient R⁻-split.
- R3 "subtree config": OFC endpoint-incidence pairs (B1) vs the task's own tightest subtree.
- Whether to fold γ-2s into R1 (rides free; not a trace band).

## Execution order
1. **R1** (PI hypothesis, highest value) — audit_133. τ-profile of burstiness, time-shuffle
   + placebo, anchor τ₀==cached.
2. **R2** dual-target trajectory — audit_134 (reuses windowed eigen cache + L reference).
3. **R3** OFC-subtree transient — audit_135 (reuses cache + B1 keep-masks).
