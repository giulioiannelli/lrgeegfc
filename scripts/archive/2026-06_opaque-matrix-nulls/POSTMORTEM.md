---
name: postmortem-opaque-matrix-nulls
type: postmortem
era: COHORT_N10
status: archived
created: 2026-06-12
---

# Post-mortem — the post-matched-strength null work (CRC + coordinated SB-CRC)

**Head.** Archived 2026-06-12 at user direction. After the (accepted)
matched-strength null, two further nulls were built — **CRC** (on-manifold
coherency surrogate, `2026-06-08`) and the **coordinated cross-phase null /
SB-CRC** (`2026-06-11`). Both operate on the *coherency matrices* via abstract
Haar rotations. They are opaque (not explainable in one sentence),
construction-dependent (verdicts moved unpredictably with backbone choice), and
were run onto real-data verdicts before being validated on synthetic
ground-truth. They are superseded by the honest recognition that the question
they tried to answer has no clean null without a control session.

## What was built
- **CRC** (`audit_96` pilot, `audit_97` cohort): rotate the complex coherency by a
  Haar real-orthogonal congruence `C̃ = O C Oᵀ`, re-derive `|ImCoh|`, push through
  the LRG cophenetic pipeline. Intended as an on-manifold *bracket* to
  matched-strength. Finding: it reproduced the matched-strength clear/no-clear
  decision 30/30 cells — i.e. matched-strength was adequate on the realizability
  axis. (That one reassurance is worth keeping as historical context.)
- **SB-CRC** (`audit_99`): hold a cross-phase shared backbone `B` fixed across all
  surrogate phases, rotate only each phase's deviation `Δ_φ`. Intended to separate
  *task-specific* persistence from a *task-independent stable fingerprint*.

## Why archived
1. **Opaque.** Matrix-level Haar-rotation surrogates cannot be explained in a
   sentence. The accepted matched-strength null can ("keep each node's total
   connectivity, shuffle who connects to whom"); these cannot.
2. **Construction-dependent.** The SB-CRC "fingerprint floor" moved — sometimes a
   lot, and in directions opposite to prediction — with the backbone weighting
   (equal-condition vs baseline). Precise decompositions ("X% task-specific") are
   not trustworthy when an unconstrained modelling knob swings them.
3. **Process error.** Ran on real-data verdicts before validating on synthetic
   ground-truth. The single-patient pilot (Pat_05) looked clean and hid a cohort
   calibration failure (non-PSD recombination → magnitude blow-up 1e2–1e8× for
   Pat_06/03/14), which only the full cohort exposed. The PSD-projection fix then
   raised floors everywhere and weakened the verdict — the off-manifold version had
   been anti-conservative all along.
4. **The deeper reason (the real one).** Separating *task-specific trace* from
   *stable trait fingerprint* has **no clean null — physical or matrix — without a
   task-free control session.** The trace lives in the cross-phase relationship of
   the very marginals any null would preserve. The machinery was papering over a
   gap in the experimental design.

## Phase-randomized (physical) surrogates — the clarification that closed this out
A physical, time-series surrogate (the principled, citable instinct) does **not**
rescue an `|ImCoh|`-based pipeline, and is **not** equivalent to matched-strength:
- **Multivariate phase randomization** (Prichard–Theiler: same random phase on all
  channels) applies one all-pass filter `H(f)=e^{iφ(f)}` to every channel. In the
  cross-spectrum `S_jk → H H* S_jk = |H|² S_jk = S_jk` — **ImCoh is preserved
  exactly.** It is a **no-op** for coherency-based FC (tests only nonlinearity,
  which `|ImCoh|` already ignores).
- **Univariate phase randomization** (independent phase per channel) destroys all
  cross-channel structure → `ImCoh → 0`. It tests "is there any connectivity at
  all," not "beyond node strength."
- So for a linear, second-order measure like `|ImCoh|`, the standard physical
  toolbox gives **either nothing (no-op) or total destruction** — no interpretable
  intermediate. That is precisely the niche a *matrix-level* null (matched-strength:
  preserve a coarse property, shuffle the fine pattern) fills, and why the field
  uses them. A physical null is a different test, not a prettier matched-strength —
  but for this pipeline it has no usable middle ground.

## What to keep / not redo
- **Keep matched-strength** as the principled null. It is enough for what can be
  claimed.
- **Do not rebuild** matrix-level coherency-rotation surrogates for the cross-phase
  trace.
- **Manuscript move:** state the limitation in one honest sentence — *"we cannot
  fully separate task-induced persistence from a stable trait-level connectivity
  backbone without a task-free control session"* — instead of a surrogate.
- Library helpers (`utils/surrogate/coherency_surrogate.py`,
  `complex_coherency_bands`) are left in place (general, harmless) but are orphaned;
  remove if a future cleanup wants them gone.

Archived files: `audit_96/97/99` (here); scope reports under
`.agents/guides/task-persistence-investigation/archive/2026-06/`; data under
`data/audit/archive/2026-06_opaque-matrix-nulls/`.
