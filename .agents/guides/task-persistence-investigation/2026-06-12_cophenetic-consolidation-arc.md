---
name: cophenetic-consolidation-arc
type: scope
era: IMCOH_ABS / COHORT_N10
status: active
created: 2026-06-12
updated: 2026-06-12
pointers:
  - scripts/01_compute/audit/audit_83_wm_stratified_cophenetic.py
  - scripts/01_compute/audit/audit_63_split_baseline_surrogate.py
  - .agents/reports/2026-06-12_trace-impact-and-leverage.md
---

# Four-phase cophenetic consolidation arc — what persists offline, the encoding or the inference?

**Head.** The established trace asks a binary question — *does `rest_post` echo
`task_test`?* (audit_83, positive = TRACE). The **arc** asks the sharper
*neurophysiological* question the four-phase design actually licenses: across
`rest_pre → task_learn → task_test → rest_post`, **when** is the persistent
β-OFC hierarchy reorganization born, and **what** does it consolidate — the
**encoding** of the premise pairs (`task_learn`) or the **inference-specific**
reorganization (`task_test` − `task_learn`, the part added when the brain must
infer the novel, never-seen relations)? This is the move from "a trace exists"
to "the trace is the offline residue of *relational inference*, not of mere
stimulus exposure." Pure reuse of verified primitives — no new connectivity
estimator, no behavioral data.

This is a **result-first** investigation: the LRG/cophenetic machinery is the
microscope; the deliverable is a statement about the brain. Scope companion to
`.agents/reports/2026-06-12_trace-impact-and-leverage.md` (Route-1-arc, the
behavior-free half of the consolidation route).

---

## Critical preamble (5 points, before code)

**(1) Claim.** The persistent offline reorganization `p = D_RP − D_pre` is not a
uniform echo of "the task." It is preferentially aligned with the
**inference-specific** reorganization `f = D_TT − D_TL` rather than with mere
**encoding** `e = D_TL − D_pre`, in the trace bands (β primary; α, low-γ
secondary). Strong form: `f` carries the offline trace *beyond* `e` (partialling
encoding out).

**(2) Null.** Matched-strength 4-cycle ±δ surrogate (R ≥ 200), generated per
phase, reusing the exact audit_63/audit_83 ensemble machinery
(`load_or_compute_eigs_at_path`, `cophenetic_condensed_from_eigs`). Per-pair
cophenetic quantities recomputed on surrogate vectors; upper-tail
p = mean(surr ≥ obs). Mandatory before any cohort claim (matched-strength rule).

**(3) Strongest plausible alternatives.**
  (i) **Encoding leakage** — `e` and `f` are not independent (`task_test` and
      `task_learn` share gross structure), so an apparent `f`-alignment could be
      `e` bleeding through.
  (ii) **Difference-of-differences noise** — `f = D_TT − D_TL` is noisier than a
      single phase contrast; a *weaker* `ρ(f, p)` could be noise, not absence of
      inference-specific persistence.
  (iii) **Phase-duration / SNR asymmetry** — if `task_learn` and `task_test`
      recordings differ in length, their cophenetic estimates differ in
      stability, biasing the `e`-vs-`f` comparison independent of cognition.
  (iv) **Strength drift, not hierarchy reorganization** — the whole arc could be
      per-node strength change rather than communication-hierarchy change.

**(4) Null's mechanical reach.** Matched-strength fixes the per-node strength
sequence exactly, so any surviving `f`-alignment is **not** explained by (iv). It
does **not** reach (i)–(iii). Those require, respectively: a **partial**
Spearman `ρ(f, p · e)` (does inference-specific persistence survive controlling
encoding); a **difference-noise floor** (compare `|f|` to the within-rest
`D_preA − D_preB` difference magnitude — if `f` is at the rest-split noise level,
abstain); and a **phase-length audit** (report `task_learn` vs `task_test`
sample counts per patient; subsample to equal length if they diverge).

**(5) Falsification + limitations.** "Inference-specific offline persistence"
requires, cohort-wide: `ρ(f, p) > 0` (Wilcoxon, one-sided) **AND** it survives
matched-strength **AND** the partial `ρ(f, p · e) > ρ(e, p · f)` (inference beats
encoding). It is **falsified** if `rest_post` aligns equally with `e` and `f`, or
only with `e` (pure encoding residue), or with neither beyond the strength null.
**Limitations:** no behavior → cannot link `f`-persistence to inference
*success* (that is the deferred behavioral test case TC1); exact `task_learn` /
`task_test` stimulus timing undocumented; n = 10; the `e`-vs-`f` contrast assumes
the two task recordings are comparable (audited under (4)/(iii)); cophenetic
vectors are heavily tied (≈ N−1 merge heights) so **Spearman throughout**
(rank-robust, identical choice to audit_83) — never Pearson on raw `D`.

---

## Notation

- `D_x ∈ ℝ^{N(N−1)/2}` — LRG cophenetic condensed vector for phase `x`, at
  `τ = 1/λ_max`, average-linkage (UPGMA), via
  `lrg_ultrametric_condensed` / `cophenetic_condensed_from_adjacency` (verified
  identical, `matched_strength.py:307`).
- Phases: `rest_pre` split into `preA`, `preB` (unbiased baseline, audit_63
  `ensure_half_fcs`); `task_learn` (TL); `task_test` (TT); `rest_post` (RP).
- Differential vectors (split baseline removes shared anatomy/strength):
  `e = D_TL − D_preA` (encoding), `g = D_TT − D_preA` (inference-online),
  `f = D_TT − D_TL` (inference-specific), `p = D_RP − D_preB` (persistent).
- Arc functionals (Spearman ρ; **positive = aligned = trace**):
  - `T_test = ρ(g, p)` — the established audit_83 full-graph trace (**reference /
    cross-check**: must reproduce `cophenetic_raw_per_patient.csv` config=full).
  - `T_learn = ρ(e, p)` — does RP echo the encoding phase.
  - `T_infspec = ρ(f, p)` — does the inference-specific reorganization persist.
  - `C_LT = ρ(e, g)` — online consistency (encoding maintained into inference).
  - Partials: `T_infspec·e = ρ(f, p · e)`, `T_learn·f = ρ(e, p · f)`
    (first-order partial Spearman on ranks).

## Properties / sanity contracts

- **Sign convention locked:** positive = trace (consistent with `T_d > 0` rule).
- **Cross-check gate:** `T_test` computed here MUST equal audit_83 full-graph
  `obs_stat` per (patient, band) to ~1e-10 (same inputs, same functions) — a
  built-in anti-hallucination check; the script asserts and reports the max
  deviation before any new number is trusted.
- **Quantization-safe:** Spearman only.

## Pseudocode

```
for pat in COHORT_N10:
  ensure_half_fcs(pat, BANDS)                       # audit_63 cache (one-time)
  for band in BANDS:
    D = { x: lrg_ultrametric_condensed(load_phase_fc(pat, x, band))
          for x in [preA, preB, TL, TT, RP] }
    e,g,f,p = D[TL]-D[preA], D[TT]-D[preA], D[TT]-D[TL], D[RP]-D[preB]
    row = { T_test=ρ(g,p), T_learn=ρ(e,p), T_infspec=ρ(f,p), C_LT=ρ(e,g),
            T_infspec_pe=partial(f,p|e), T_learn_pf=partial(e,p|f),
            n_TL=len(tl_samples), n_TT=len(tt_samples) }   # (iii) length audit
cohort: per band → median, k/10 positive, IQR of each functional
validate: assert max|T_test − cached audit_83 obs_stat| < 1e-6
```

## Visualization (next stage)

- Per band, cohort: bar of `T_test` vs `T_learn` vs `T_infspec` (median + per-
  patient dots) — the headline "what persists" panel.
- 5×5 cross-phase cophenetic-ρ heatmap (preA, preB, TL, TT, RP), cohort-median,
  per band — the trajectory geometry.
- Defer matplotlib until the observed pattern is read with the user (no
  pre-registered gate; compute + inspect first).

## Connection to prior tools

Extends audit_83's split-baseline full-graph trace (which collapses the task to
`task_test`) to the **full four-phase arc** and adds the **encoding vs
inference-specific decomposition**. Same primitives, same null, same sign
convention. Orthogonal to the WM/epi stratifications (those re-slice *nodes*;
this re-slices *phases*).

## Open questions

- Does `f`-persistence (inference-specific offline trace) predict per-patient
  inference performance / the symbolic-distance effect? → behavioral TC1
  (deferred; needs the obtainable RT/accuracy data).
- Is the OFC concentration of the *trace* also the OFC concentration of the
  *inference-specific* component `f`? (Re-run the audit_83 localization on `f`.)
- Does `task_learn` itself already show an OFC reorganization (encoding), or does
  OFC engage only at inference (`f`)?
