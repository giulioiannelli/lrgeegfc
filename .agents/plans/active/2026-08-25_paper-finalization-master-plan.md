---
name: 2026-08-25_paper-finalization-master-plan
kind: plan
era: IMCOH_ABS × COHORT_N10 → PAPER_FINALIZATION (fundamentals re-opened)
status: active
created: 2026-08-25
scope: Master plan for finalizing the paper. Re-opens the four fundamentals (FC transform, sparsification, τ/multiscale, nulls) as a pre-registered Wave 0, then runs the three science lanes (multiscale persistence / band-and-place dissociation / SOZ detection) in parallel worktrees, then synthesizes the Nature-Neuroscience-shaped narrative. Written after a fresh read of the repo state, not from prior summaries.
pointers:
  - .agents/reports/2026-07-23_old-vs-current-headlines-panoramic.md   # the honest July diff
  - .agents/reports/2026-07-23_raw-inclusion-scale-selection-and-global-methods.md
  - .agents/reports/2026-07-15_sparsifier-suite-verdict.md
  - .agents/guides/02_methods/sparsification-choice.md
  - .agents/preprint/established_results/2026-07-13_settled-three-results.md
  - scripts/01_compute/sparsified_arc/13_matched_strength_mst020.py    # the canonical trace-gate harness
  - src/lrg_eegfc/utils/fc/heat_multiscale.py                          # rho_sym estimator
---

# Paper finalization — master plan

## Head

The paper has three results and one operator, and the operator is fine. What is not fine is the **substrate** it runs on: every headline currently rests on three choices that were made by argument rather than by evidence (|ImCoh| over |ImCoh|²), by outcome-selection (backbone density `f = 0.20`, chosen because it is "the unique fraction where α and β are jointly clean"), and by a null injected at the very last stage of the pipeline (matched-strength, which shuffles the finished FC matrix and therefore cannot test the connectivity estimator, the band split, session nonstationarity, or drift). Wave 0 re-derives those three choices from scratch under pre-registered rules, and adds the null the project has never run — a **data-based, timeseries-level ladder** whose centrepiece is the circular-shift surrogate, which for a coherency-based measure preserves coupling *magnitude* exactly and destroys only the *lag* structure that |ImCoh| claims to measure. Only after that contract is locked do the three science lanes run. This ordering is the whole point: the reason the project "got lost multiple times" is that the substrate moved underneath results that had already been written up.

---

## 1. Where we actually stand (fresh read, 2026-08-25)

### 1.1 The pipeline, end to end

`timeseries → Welch CSD → complex coherency C(f) → band transform ⟨|Im C|⟩_f → dense W → backbone → L = D − W → eig → heat kernel ρ̂(τ) = e^{−τL}/Tr → T = 1/ρ̂ → UPGMA → cophenetic D(s), s ≡ τλ_max → cross-phase ρ_sym → cohort Wilcoxon vs null`.

The cross-phase estimator is `ρ_sym = ½[ Spearman(D_task − D_A, D_post − D_B) + Spearman(D_task − D_B, D_post − D_A) ]`, where A and B are split halves of `rest_pre`. It is a good estimator: split-half referenced, symmetric in the half assignment, invariant to global drift in the *level* of the distances. Its ceiling is the split-half reliability itself (≈ 0.5), not 1.

### 1.2 What July said — read as claims to re-derive, not as premises

Every row below was produced on the substrate this effort is rebuilding, so an era-specific *negative* carries no more authority than an era-specific *positive*. Nothing here licenses dropping a hypothesis from scope. See §7.

| result | current status | what it rests on |
|---|---|---|
| β trace held in `rest_post`, **scale-invariant** (clears 16/16 scales, Friedman p = 0.17) | strongest surviving claim | mst@0.20 + matched-strength only |
| α trace, **scale-tuned** (12/16, peak s ≈ 5) | holds, weaker | same |
| θ / low_γ **null** in the hierarchy — band-selectivity | holds, and is the *sharpest* multiscale contribution | same |
| Encoding persists (β 16/16); inference-specific persists — bands contested | **CLAIM TO RE-DERIVE** (see §7); no prior verdict on enc/inf is a premise | same |
| SOZ marker from seeded diffusion: δ 0.83 / low_γ 0.82 / β 0.745 (mst@0.20); fused 0.904 (TMFG) | holds | matched-strength fake-SOZ null |
| **Anatomical localization** (β→OFC, encoding→OFC, inference→cingulate) | **CLAIM TO RE-DERIVE** — the July sweep found 0 cells clearing BH, but on a substrate this effort is rebuilding, and it never tested intrinsic communities | 2026-07-13 sweep |
| **"Only the hierarchy sees the trace"** | **UNSUPPORTED** — `coph\|raw` residual never clears; `raw\|coph` clears strongly | killed 2026-07-23 |
| Interpatient spread of β explained by left-hemisphere implant coverage (ρ = 0.68) | holds, but is a **coverage confound**, not a biological result | — |

### 1.3 The six fractures the user named, diagnosed

1. **|ImCoh| vs |ImCoh|² was never decided by evidence.** The 2026-04-15 reset fixed a genuine order-of-operations bug (Jensen: `mean(|·|) ≠ |mean(·)|`) and picked `imcoh_abs` on the strength of a citation (Ewald 2012 / Bastos & Schoffelen 2016). No head-to-head comparison of reliability, null floor, or downstream behaviour was ever run. Both are defensible in the literature; we currently cannot say why ours is the right one.
2. **The backbone fraction is outcome-selected.** `.agents/guides/02_methods/sparsification-choice.md` §4a states it plainly: "α and β have opposite density needs (α peaks sparse 15/16 @ 0.05; β needs density 16/16 @ 0.20) — only 0.20 gives α strong + β maximal + θ/low_γ null." That is choosing the knob so the answer comes out right. The guide's own criterion 3 was **parameter-free**, and `mst@f` is listed in its own method table as "reject (free parameter f)". This is the single most reviewer-fatal item in the manuscript. The user's phrasing — "not trashly changing when we change thresholds" — is exactly correct as a requirement.
3. **The τ / multiscale claim is currently split between a dead version and a live one.** Dead: "the hierarchy detects what raw FC cannot" (script 21: `coph|raw` never clears for α or β; `raw|coph` clears strongly). Live: "the hierarchy *selects*" (script 28: raw fires in δ/α/β/low_γ non-selectively; the cophenetic keeps β at all scales, α at s ≥ 1, and rejects low_γ at **every** scale) plus "the hierarchy *characterizes*" (β scale-invariant vs α scale-tuned). The live version is real but has never been given a proper cross-patient significance test — the scale-invariance verdict is a Friedman on n = 10.
4. **Cross-patient validity is the load-bearing gate and it is under-specified.** The gate is a one-sided Wilcoxon on `obs_rho − surr_p50` across 10 patients, per scale. Known problems, from the project's own diary: per-patient `obs_rho` co-varies with its *own* null height (β ρ = +0.83), so raw `obs_rho` spreads partly reflect baseline non-identifiability, not effect size; multiplicity across 16 scales × 6 bands is handled inconsistently; the LOO fragility metric in use was degenerate (constant 1).
5. **Localization is dead at every imposed parcellation and the replacement is a coverage confound.** Nothing focal survives for α or β. β is "coarsely left-lateralized" — but β trace strength tracks left-hemisphere *contact fraction* at ρ = 0.68, which is a statement about where the electrodes are, not where the biology is. The one test never run is the intrinsic one: whether the trace concentrates in the diffusion's **own** communities at the band's own scale, rather than in imposed anatomical labels.
6. **The null is injected at the last possible stage.** Matched-strength shuffles the finished N×N FC matrix. It therefore cannot test: the coherency estimator, the choice of band transform, the frequency-band split, session nonstationarity, slow drift, artifact epochs, or the phase segmentation itself. Every claim in the paper is conditioned on "given this FC matrix". No timeseries-level null has ever been run for the trace.

### 1.4 The two that can sink the paper

Fracture **2** (knob-selected backbone) and fracture **6** (null too late in the pipeline). A referee who notices either one can reject the whole thing without engaging with the science. Wave 0 exists to close them.

---

## 2. The narrative arc we are aiming at

One operator, read three ways, answering three questions that no single-scale connectivity analysis can pose together:

> **A reasoning episode leaves a trace in the resting brain that is organized by scale, and the organizing scale is set by frequency band.** What is held is not a replay of the premises the patient was *shown* but the relational structure they had to *compute* — the five-phase design separates encoding (`task_learn`) from inference (`task_test`), and the two persistences are separable. Reading the same diffusion operator *across* phases recovers the trace (β holds it at every scale — no characteristic scale; α holds it only in a narrow band of scales); reading it *within* a phase recovers the epileptogenic zone as a co-diffusing community. The bands dissociate the two: cognition lives in α/β, disease in δ/low-γ, and β is the only band in both.

The claim that makes this Nature-Neuroscience-shaped rather than a methods paper is the middle one — **that the persistence of a cognitive state is a scale-structured object, not a scalar** — and it is the claim Wave 0 must either harden or force us to drop. The honest fallback, if the multiscale axis does not survive, is a strong band-dissociation + epilepsy-marker paper, which is a different and smaller paper. We decide that on evidence, in Wave 0, not in the writing.

Three things this arc requires and does not yet have: (i) a substrate that a referee cannot call tuned; (ii) a null that reaches back to the data; (iii) an account of *where* that is either genuinely intrinsic or honestly absent.

---

## 3. Wave 0 — fundamentals (three parallel worktrees, running now)

Rules binding on all three: pre-register the decision rule in the script docstring **before** looking at the numbers (5-point critical preamble); matched-strength stays as the incumbent null while the new ladder is built; report **per-scale, never best-scale**; recompute every load-bearing number fresh from cached `imcoh` freq-resolved arrays or from raw timeseries — never source a number from a prior report or CSV.

### W0-A — the substrate contract (FC transform × sparsification, jointly)

They interact: the transform changes the weight distribution, which changes what any density-ranked backbone keeps. Decide them together or not at all.

- **A1 — transform.** `|ImCoh|` vs `|ImCoh|²` vs the signed-derived alternatives, judged on properties that are *independent of our hypothesis*: split-half reliability of the FC matrix itself; the null floor under circular shift; dynamic range and weight-distribution shape; the effect on backbone composition. Downstream trace/AUC numbers are recorded but **must not be the selection criterion**.
- **A2 — sparsification stability map.** Sweep the density knob continuously (`f` ≈ 0.02 → 1.0, ~20 values) × {mst-union, disparity, TMFG, PMFG, plain threshold, percolation}. The deliverable is **not** a best `f`. It is the surface: over which region of the knob is the verdict invariant? A plateau is a result; the absence of a plateau is also a result and must be reported as knob-dependence.
- **A3 — kill the free parameter if possible.** Two routes, both tested: (i) a genuinely parameter-free filter whose density is fixed by construction, evaluated against the plateau found in A2; (ii) a **knob-integrated readout** — integrate or median the statistic over the invariance plateau so no single `f` is ever reported. Prefer (i) if a parameter-free filter sits inside the plateau; otherwise ship (ii) and report the plateau width.
- **Deliverable:** `.agents/preprint/locked/PIPELINE_CONTRACT.md` (the locked substrate, with the pre-registered rule and the evidence) + one library entry point `canonical_graph(patient, phase, band)` that every downstream lane calls, so the three lanes cannot silently diverge.

### W0-B — the null ladder (data-based, timeseries-level)

The centrepiece insight, which makes this both cheap and sharp: for a coherency-based measure a **circular time shift is a pure phase rotation in the frequency domain**. Shifting channel `c` by Δ maps `S_ij(f) → S_ij(f)·e^{−2πifΔ/L}` — coherence *magnitude* is preserved exactly, only the phase (the lag) is randomized. So the circular-shift surrogate preserves every bit of the coupling magnitude structure, including all volume conduction, and destroys precisely the time-lagged interaction that |ImCoh| claims to measure. That is the fairest and most damaging null available for this project, and it has never been run.

It is also nearly free to compute: cache the **complex** band coherency `C(f)` once per (patient, phase) and each surrogate is an elementwise phase-ramp multiply. `complex_coherency_band()` already exists in `utils/surrogate/coherency_surrogate.py`.

- **N1 — circular shift** (per channel, whole recording). Null: "the trace is carried by coherence magnitude, not by lag." *The key null.*
- **N2 — phase randomization** (independent per-channel phases). Stricter floor; destroys both lag and cross-spectral consistency.
- **N3 — block phase-label permutation** ★. Concatenate the session, cut into contiguous blocks ≫ `nperseg`, reassign blocks to pseudo-{A, B, task, rest_post} preserving the real per-phase durations, recompute FC and ρ_sym. Null: "any partition of this session produces this geometry." This is the only null that reaches session nonstationarity, drift, and artifact epochs. Run both a free version and an **order-preserving** version (move the cut points, keep temporal order) — the latter is the honest drift-aware test that the retired drift null was trying and failing to be.
- **N4 — within-rest placebo.** All four pseudo-phases drawn from `rest_pre` alone. Null: "a session with no task produces this ρ_sym."
- **Deliverable:** `src/lrg_eegfc/utils/surrogate/timeseries_nulls.py` + a calibration report giving, per band and per scale, the observed ρ_sym against each rung of the ladder, and an explicit verdict on which rungs the α and β traces survive. If β survives N1 and N3, the paper's central claim is on ground it has never had before. If it does not, we need to know that now.

### W0-C — the cohort contract and the τ question

- **C1 — lock the cross-patient statistic.** One canonical `cohort_gate()` that every lane calls: margin-based (`obs − surr_p50`, never raw `obs_rho`), per-scale, whole-grid BH across the (band × scale) family, leave-one-patient-out of the *verdict* (not of the statistic), bootstrap CI on the cohort effect, and per-patient reporting that is descriptive only — patient counts are never the gate.
- **C2 — the τ question, decomposed into three separable claims** so we stop conflating them:
  - **(a) Detection** — does any scale see what `s → 0` (raw FC) cannot? Current answer: no. Re-verify across all scales and with a residualization that is not rank-linear, then close it.
  - **(b) Characterization** — is the *shape* of ρ_sym(s) band-discriminative, and is that discrimination significant *across patients*? Not a Friedman on the cohort median: a per-patient scale-profile test.
  - **(c) Selection** — the hierarchy rejects low_γ at every scale while raw FC accepts it. Is that rejection principled (and can we show *why*, e.g. the edge-locality hypothesis: low_γ's trace sits on a few strong edges that clustering absorbs, β's is distributed and survives coarse-graining)?
- **Deliverable:** the honest verdict on what τ buys, whatever it is, plus the locked cohort gate as library code.

**Wave 0 exit criterion:** `PIPELINE_CONTRACT.md` locked, null ladder run on α/β/θ/low_γ, cohort gate in the library. Nothing in Wave 1 starts before that.

---

## 4. Wave 1 — the three science lanes (parallel worktrees, after the contract)

Each lane inherits the locked contract, uses the null ladder, and passes the cohort gate. Each owns a distinct output namespace and a distinct set of files.

- **L1 — persistence as a multiscale feature.** The trace across scales, under the full null ladder; β scale-invariance vs α scale-tuning as a per-patient, cross-patient-tested claim; the selection result (what the hierarchy rejects and why); the physical meaning of `s` (`ℓ(s)`, cluster counts) so no scale is ever labelled micro/meso/macro without a number.
- **L2 — learning vs applying, and where/when each persists.** The lane's core is the **encoding vs inference dissociation**, re-derived from scratch: does the resting brain hold the premises it was *shown* (`task_learn`) or the ordering it had to *infer* (`task_test`), and are those two persistences separable in band, in scale, and in space? Band dissociation and localization are read off the *same* five-phase functionals rather than as separate analyses — localization redone as an **intrinsic** question (the diffusion's own communities at each band's own scale) rather than imposed parcels, with the implant-coverage confound handled explicitly rather than reported as laterality.
- **L3 — from persistence to SOZ detection.** The within-phase read of the same operator; backbone posture settled by the Wave 0 contract (this is where the two-backbone question is resolved); honest separation of what τ contributes (ranking) from what band fusion contributes (precision); leakage-free leave-one-patient-out evaluation.

## 5. Wave 2 — synthesis

Reconcile the three lanes into one Results narrative, rewrite `results_sec_1/2/3.tex` and `methods.tex` against the locked contract, retire every stale directive (the three 2026-07-13 conceptual-diff directives and the OVERVIEW §5 are stale on localization and must not be written from), and produce the figure set.

## 6. Ledger of what must be retired on contact

- `.agents/preprint/directives/writing_directive_2026-07-13_*` (§5 and the "What STAYS" localization bullets) — stale, assert β→OFC / encoding→OFC / inference→cingulate.
- Any claim of the form "only the multiscale hierarchy detects the trace" — unsupported.
- Any scale labelled micro / meso / macro without `ℓ(s)` or a cluster count.
- `mst@0.20` as a bare stated choice — either justified by a plateau or replaced.

---

## 7. Amendment 2026-08-25 — encoding vs inference is first-class, and no old verdict is a premise

**User correction, and it corrects §1.2 of this plan.** The July retirements of the encoding/inference layer ("inference in β *alone* is dead", "inference→cingulate is dead", "the δ inference-specific effect is a partial-correlation artifact") were carried into §1.2 above as *settled*. They are not settled and must not be treated as premises. **A prior verdict is a claim to re-derive, never a reason to drop a hypothesis from scope** — those verdicts were produced on a substrate (transform, backbone, null, statistic) that this very effort is rebuilding, so an era-specific negative carries no more authority than an era-specific positive.

**The learn-versus-apply contrast is fundamental to the paper.** The paradigm is transitive inference: `task_learn` presents ordered premises (A>B, B>C, …); `task_test` requires judging novel non-adjacent pairs, answerable only by inferring from the learned ordering. Comparing the persistence of *learning a structure* against the persistence of *applying* it is what makes this a cognition paper rather than a generic task-effect paper, and it is the axis along which the offline-consolidation interpretation lives.

**The construct to preserve** (`scripts/01_compute/sparsified_arc/05_enc_inf_arc.py`, five phases A, B, `task_learn`, `task_test`, `rest_post`; all functionals symmetrized over the A/B arms):

| term | definition | reads |
|---|---|---|
| `e` | `D_learn − D_A` | encoding change |
| `f` | `D_test − D_learn` | inference-specific change — **arm-invariant**, never references `rest_pre` |
| `g` | `D_test − D_A` | standard task change |
| `p` | `D_post − D_B` | persistence |
| `T_test` | `ρ(g,p)` | standard trace |
| `T_learn` | `ρ(e,p)` | does what was *encoded* persist? |
| `T_infspec` | `ρ(f,p)` | does what was *inferred* persist? |
| `T_infspec_pe` | `partial ρ(f,p \| e)` | does inference persist *beyond* encoding — the load-bearing cognitive claim |

`f` never touches `rest_pre`, so it is structurally immune to baseline-half artifacts. That is a real strength of the construct and should be stated as such.

**Standing caveat that travels with the conditional statistic.** `T_infspec_pe` is a partial correlation and can manufacture a positive value from null input: a sham arc built entirely from pre-task `rest_pre`, where no consolidation signal can exist, once returned a *significantly positive* conditional trace (p = 0.007), exceeding the real value (`.agents/preprint/supplementary/S2_drift_controls.md`). **Every null must therefore be calibrated on a no-signal arc before any p-value for a conditional functional is reported**, and a null that fails calibration for that functional is invalid for it and must be declared so rather than quoted.

**Consequences already dispatched to Wave 0:** the phase set is data, never a four-phase hardcode (W0-A); the substrate stability plateau is validated on `T_learn` and `T_infspec_pe` as well as `T_test`, and a disagreement between plateaus is a major finding, not something to resolve silently in favour of `T_test` (W0-A); the null ladder runs on the enc/inf functionals with `T_infspec_pe` first, and the calibration table (functional × rung × value-on-no-signal-input) is a first-class deliverable (W0-B); `cohort_gate()` is functional-agnostic and refuses to report a p-value for a statistic/null pair that fails calibration (W0-C).

**One new question added, and it may be the most valuable in Wave 0.** `T_learn(s)` and `T_infspec_pe(s)` are two curves over the same scale axis, for the same patient and band. **Do encoding and inference persist at different scales?** A scale dissociation between learning and applying would be a genuinely multiscale cognitive result that no single-scale method can produce, and it would tie the τ axis (§3 W0-C) directly to the cognitive dissociation instead of leaving them as separate stories. Assigned to W0-C; to be reported honestly including a negative.
