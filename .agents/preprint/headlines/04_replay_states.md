---
name: headline-n4-replay-states
era: IMCOH_ABS_COHORT_N10
status: scoped_not_executed
kind: headline
scope: N4 (★ the paper's central idea) — the task-state RE-ENTERS rest as transient, recurring multiscale-multiband reorganization "states" (replay/reinstatement). Identifiable specifically because the LRG characterizes them as multiscale + multiband patterns. Dynamic complement to N1's static trace. Infra exists; NOT yet executed.
owner_agent: dynamic-FC / replay agent (new) + null-model chat
canonical_scope: .agents/guides/task-persistence-investigation/2026-06-22_replay-states-multiscale-reinstatement.md
updated: 2026-06-22
---

# N4 ★ — Replay states: the task configuration transiently re-enters rest as recurring multiscale-multiband reorganization

> **This is the paper's central idea** (PI, 2026-06-22). N1 establishes the
> *static, phase-averaged* task signature; N4 asks whether that signature recurs
> as discrete *states* in time-resolved rest. **STATUS: scoped, infra exists, NOT
> yet executed.** No results below are claimed — this file crystallizes the
> headline and points to the execution+verification brief (canonical scope, see
> frontmatter). Everything here is hypothesis + plan.

## §A — The idea, in plain language

N1 says: *on average*, the post-task resting network sits closer to the task
configuration than the pre-task resting network did. That is a **static** claim —
it averages over the whole rest period. The deeper, more biological question:
**is that average shift carried by discrete moments?** — brief windows in which the
brain *snaps back* into the task's connectivity configuration, then relaxes, then
snaps back again. If so, those windows are **replay / reinstatement states**: the
network transiently *re-enters* the task pattern offline, the way memory-
consolidation theory says a learned structure is rehearsed during rest.

The selling point — and why this needs *our* method — is **how the states are
identified**. A replay state here is not a single edge lighting up or one global
mode wobbling; it is a **multiscale, multiband reorganization** matching the task
signature: the same per-pair hierarchical geometry (ρ^coph) and the same
across-band pattern that defines the N1 trace, appearing transiently in a short
time window. Only a multiscale + multiband characterization can *recognize* such a
state; edge-wise or single-band or global-spectral views cannot (this is the CORE
payoff, made dynamic). So N4 turns the trace from a *static average* into a
*process*, and it does so with the method's signature capability.

**Why it matters:** it would be, to our knowledge, a human intracranial signature
of **offline reinstatement of a learned relational structure**, identified as a
recurring multiscale-multiband network state — the dynamic heart of the paper.

## §B — Technical statement (the plan; subheadlines)

> Full notation, predicates, pseudocode, nulls, and falsification are in the
> canonical scope report (frontmatter `canonical_scope`). Summary here.

**N4.1 — Time-resolved signature.** Slide a window over each phase; per window
build the |ImCoh| graph → LRG → per-pair `ρ^coph` (and the band-multiplex
signature). Define a per-window **task-likeness** score = similarity of the
window's multiscale-multiband configuration to the task-state configuration,
referenced to the `rest_pre` baseline (the same geometry that defines the N1
trace, evaluated per window instead of phase-averaged).

**N4.2 — Replay states = bursts beyond a shuffled-time null.** A replay state is a
window (or run of windows) whose task-likeness exceeds what a **time-shuffled /
phase-randomized** null and a **matched-strength per-window** null produce. The
claim is not "elevated mean" (that is just N1) but **transient, bursty structure**
— dwell times / event clustering above a stationary null.

**N4.3 — Identifiable only via multiscale + multiband.** Show the states are *not*
recoverable from a single edge, a single band, or the global leading-mode
(Grassmann) view alone — the recognition requires the per-pair multiscale +
multiband signature. This is N4's tie to CORE and its strongest methodological
selling point.

**N4.4 — rest_post vs rest_pre asymmetry.** Consolidation direction: replay states
should be **more frequent / more intense in `rest_post`** than `rest_pre`. The
`rest_pre` distribution is the within-subject control for "states that occur
anyway." (And the four-phase design lets `task_learn` vs `task_test`
configurations be tested as separate replay targets — see N2 / `task_learn`.)

**N4.5 — STATUS.** Scoped, not executed. Sliding-window infra exists
(`compute_time_windows.py`, `visualize_time_windows.py`). The **binding
precondition** is per-window |ImCoh| estimator stability (short windows × low-band
frequency resolution) — must be cleared *before* any replay claim.

## §C — Critical issues & powerful strengths

**Strengths (if it lands):** dynamic, mechanistic, theory-connected (replay /
consolidation); uses the method's unique multiscale-multiband recognition; the
rest_pre within-subject control is clean; it is the conceptual capstone that
unifies CORE + N1 + N2.

**Lead-with risks (all live, because nothing is run yet):**
- **|ImCoh| at short windows** — imaginary coherence needs enough cycles; low
  bands (δ/θ) may be unestimable at replay-relevant window lengths. *Precondition,
  not a footnote.*
- **"States" vs "smeared mean"** — the whole result dies if elevated task-likeness
  is stationary (no bursts). The burstiness/dwell-time test against a shuffled-time
  null is the make-or-break analysis, not the mean.
- **Specificity** — replay must match the *task* configuration, not any non-rest
  state; needs a random-alternative-config control.
- **n=10, sEEG sparsity** — states are within-subject; cohort aggregation must
  respect heterogeneous implants.
- **Outliers** — the N1 anti patients (Pat_15/10) may lack a task signature to
  replay; report per-patient.

## §D — To-dos & verifiables (owner: dynamic-FC/replay agent + null-model chat)

- [ ] **Execute the canonical scope** (frontmatter link) — it carries the 5-point
  preamble, predicates, nulls, pseudocode. Do the precondition first.
- [ ] **Precondition:** per-window |ImCoh| stability sweep (window length vs band);
  pick the shortest window with acceptable estimator variance per band.
- [ ] **Burstiness test:** dwell-time / event-rate of high-task-likeness windows vs
  a time-shuffled null (the core analysis).
- [ ] **Multiscale-multiband necessity:** ablate to edge-only / single-band /
  Grassmann-only and show the states vanish (the CORE-payoff demonstration).
- [ ] **rest_post vs rest_pre** asymmetry + (optional) learn-target vs test-target
  replay.
- [ ] Decide paper position: lead headline vs N1→N4 build (recommend stating it as
  the central idea, presented after N1 sets up the signature).

## §E — Figure / representation ideas

- **Task-likeness time-course** per phase with replay states highlighted (the
  money figure) — `rest_post` vs `rest_pre` side by side.
- **State raster** across patients (event times), + dwell-time distributions vs
  null.
- A **single replay-window dendrogram** snapping to the task configuration next to
  a baseline window (shows the multiscale recognition concretely).
- Ablation panel: detectability under edge-only / single-band / Grassmann vs full
  multiscale-multiband (CORE payoff).
- Keep figure rules: X-epi variant where relevant; no C4.

## §F — Provenance (infra + scope; NO results yet)

- Canonical scope / execution+verification brief →
  `.agents/guides/task-persistence-investigation/2026-06-22_replay-states-multiscale-reinstatement.md`.
- Sliding-window infra → `scripts/01_compute/batch/compute_time_windows.py`,
  `scripts/02_visualize/visualize_time_windows.py`, CLI `compute` time-window path.
- Static signature it reinstates → N1 (`01_trace.md`) + `../locked/VERDICT_LEDGER.md`.
- Idea origin → `.agents/reports/2026-06-18_consolidation-arc-handoff.md` ("replay
  events" next-step) + `2026-06-12_trace-impact-and-leverage.md`.
- **No `data/audit/` outputs exist for N4 yet** — by design; this headline is
  scoped, not executed.

## §G — Missing parts / open

- Everything empirical (it is unrun). The first deliverable is the precondition +
  burstiness test.
- Whether replay targets the `task_learn` (encoding) or `task_test` (inference)
  configuration preferentially — ties N4 to N2 and `task_learn`.
- Whether replay-state rate relates to the N1 trace magnitude per patient
  (within-subject, no behavior needed).
