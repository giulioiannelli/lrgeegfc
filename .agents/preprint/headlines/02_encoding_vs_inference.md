---
name: headline-n2-encoding-vs-inference
era: IMCOH_ABS_COHORT_N10
status: current
kind: headline
scope: N2 — what the persistence is ABOUT. The four-phase design (rest→LEARN→TEST→rest) and task_learn as the encoding reference let the trace be decomposed into memory encoding vs relational inference, separating by frequency and cortical territory. Duration control is live (audit_103b).
owner_agent: inference (consolidation arc) chat
updated: 2026-06-22
---

# N2 — The persistent trace decomposes into cognitive content: memory *encoding* vs relational *inference*, separating by frequency and cortical territory

> Emerges from N1: once we know a trace persists, we ask *what it is a trace of*.
> The four-phase design — and `task_learn`, which N1 ignored — is what makes the
> decomposition possible. Numbers in cached CSVs (§F).

## §A — Result in plain language

The task has four recordings: **rest → learn → test → rest**. In **learning** the
patient is *shown* premise pairs (this is **encoding**). In **testing** the
patient is given novel pairs and must *work out* the answer (this is
**inference**). N1 measured persistence using the test phase only; N2 uses
**both** task phases, and that is the whole point — `task_learn` is the reference
that lets us subtract "what was merely seen" from "what was figured out."

What the brain keeps offline splits cleanly:

- **The α rhythm keeps the *memorised* material only.**
- **The β rhythm keeps the memorised material *and* the *inferred* relations.**

So the inference-specific residue — a trace of *what the patient worked out*, not
just what they saw — lives specifically in **β**. The β-specificity is itself the
defense against the boring explanation: a recording-length difference would
affect all rhythms, not one. **This defense is now empirically confirmed**
(§C, duration control resolved 2026-06-22): the test phase is 1.35–2.53× longer
than the learn phase, but the β inference residue does **not** scale with that
length difference across patients, whereas **α** — the band that would have to
inflate to break β-specificity — does show the duration-tracking trend. The
confound that could fake the effect rides α, not β.

Where the two kinds of content sit also differs:

- **Encoding** concentrates in **orbitofrontal cortex** (the same OFC locus as the
  main β trace) — this is **robust**.
- **Inference** leans toward the **cingulate** — but this is currently a
  **directional hint, not an established localization** (see §C: the duration
  control demoted it).
- The cingulate carries *encoding* in α and low-γ but flips to *inference* in β —
  a content-by-frequency double dissociation derived from anatomy alone.

**Why it matters:** this turns "a trace persists" (N1) into "the trace carries the
*cognitive product* of the task." **It holds up under the duration control** (§C,
resolved 2026-06-22 — β is duration-robust). So this is an offline neural signature
of **consolidating an inferred relational structure** — the cognitive climax of the
paper. *Honest ceiling:* the inference component is distributed with **no
established anatomical hotspot** (cingulate is a directional hint) and **no
behavioral validation** (performance data unavailable); sell it as a signature of
inference-specific reorganization, never as a localized or behaviorally-tied effect.

## §B — Technical statement (per subheadline)

**N2.1 — The four-phase decomposition (task_learn's role).** Per region-pair
cophenetic distances: encoding `e = D_task_learn − D_rest_pre`, inference-specific
`f = D_task_test − D_task_learn`, persistence `p = D_rest_post − D_rest_pre`. The
consolidation arc (`audit_103`) tests whether task-change persists:
`T_learn = corr(e,p)` (encoding echo) and `T_infspec·e = partial corr(f,p | e)`
(inference, controlling encoding). Validated against the locked N1 trace exactly
(reproduces `audit_83` per-patient to ~1e-16). **`task_learn` is the pivot:**
`f` and `T_infspec·e` do not exist without it.

**N2.2 — Encoding echo persists (α & β); inference-specific is β-only.** Under the
mandatory matched-strength null: the encoding echo `T_learn` persists in α and β;
the inference-specific component `T_infspec·e` persists in **β only**
(p≈0.007, LOO-robust, all other bands null; `audit_103 --null`). **Duration-
controlled (resolved 2026-06-22, §N2.5):** the β component does not scale with the
test/learn length ratio across patients (Spearman ρ≈+0.25, n.s.), while α — the
only other candidate — carries the duration-tracking trend (ρ≈+0.53). So β-only is
the *duration-clean* verdict, not an artifact of the longer test recording.

**N2.3 — Content localizes differently.** Feeding the localizer the encoding vs
inference per-pair vectors (`audit_110`, β at R=1000): **standard/encoding →
OFC** (robust, BH-clears); **inference → cingulate** (BH-borderline at full
length). Double dissociation: cingulate = encoding in α/low-γ, inference in β.

**N2.4 — A focal low-γ encoding trace hidden by whole-brain averaging.** There is
*no* whole-brain low-γ trace, yet an **absolute** within-region test
(`audit_112`, no demeaning) finds a genuine **focal cingulate encoding trace** in
low-γ (8/8 patients; survives SOZ-contact removal). A real trace that cohort-wide
averaging conceals — and an argument for the per-region read-out. (Power floor:
systems sampled by ≤5 patients can't reach significance on this absolute test, so
OFC's evidence stays the demeaned localization of N2.3, not this.)

**N2.5 — The duration confound (RESOLVED 2026-06-22).** `task_test` is **1.35–2.53×
longer** than `task_learn` in every patient, and `f = D_test − D_learn` inherits
that asymmetry. The verdict, across two levels:

- *Localization* (`audit_113b`, demeaned per-system null): **splits.** *Encoding →
  OFC* is duration-invariant and **survives**; *inference → cingulate* does **not**
  survive matched length → **downgraded to a directional hint** (consistent sign,
  but significance was duration-assisted).
- *Whole-brain arc* (N2.2): the naive instrument — a length-matched strength null
  rebuilt on a **truncated** `task_test` (`audit_103b`) — is **INVALID and retired
  from the evidence.** Its negative control fails: **δ, dead-null at full length
  (p≈0.54), false-positives at p≈0.014 at both the head and center windows**, and
  at center even out-clears β. A contiguous truncation injects a common-mode
  structural shift into `D_test` that the *un-demeaned* whole-brain concordance
  reads as fake signal (worst in slow bands) — the same artifact the *demeaned*
  localization removes, which is why that control stayed clean. **Methodological
  takeaway worth one sentence in the paper:** whole-brain concordance is
  corruptible by contiguous-window truncation; the per-system demeaned read-out is
  not.
- *Clean replacement control:* a per-patient regression of the full-length
  `T_infspec·e` against the test/learn length ratio (artifact-free). **β does not
  scale with duration** (ρ≈+0.25, p≈0.49); **α carries a duration-tracking trend**
  (ρ≈+0.53). Combined with the full-length matched-strength null (the mandatory,
  validated referee, which β already clears at p≈0.007), this establishes the
  **β inference-specific arc as duration-ROBUST and α as duration-suspect** — so
  N2.2 ships **duration-controlled, not provisional.**

## §C — Critical issues & powerful strengths

**Strengths:** the encoding side (encoding echo + encoding→OFC) is robust and
duration-invariant; the four-phase decomposition is exact (reproduces N1);
β-specificity is a strong intrinsic argument against generic confounds; the
low-γ focal encoding trace is a clean "averaging hides it" demonstration.

**Lead-with weaknesses:**
- **The inference story is now duration-controlled but anatomically/behaviorally
  thin.** The whole-brain inference-specific β arc is **duration-robust** (clean
  ratio regression + full-length matched-strength; the naive truncation control was
  retired as invalid — N2.5), but it remains **distributed with no established
  hotspot** (inference→cingulate is a directional hint, downgraded by the
  matched-length localization null) and has **no behavioral anchor**. Sell it as a
  signature of inference-specific reorganization — not as localized or
  behaviorally-validated.
- **No behavioral/performance data — and none is obtainable** (PI 2026-06-22).
  The "inference mark" is therefore the neural residue of inference-specific
  *reorganization* and can **never** be validated against inference *success*;
  the claim's ceiling is the task-phase decomposition + anatomy + literature.
  Do not promise a brain–behavior correlation.
- **`task_learn` is short** (the duration asymmetry) and was historically
  neglected — its recording quality / length per patient should be audited.
- **What about α?** α = **encoding-only** and its localization is **diffuse** (no
  FDR-surviving concentration) — i.e. α is the broad "memory" band, β the focal
  "inference-carrying" band. State this contrast explicitly; it is a result, not
  a gap.

## §D — To-dos & verifiables (owner: inference chat)

- [x] **Stage-2-arc (DONE 2026-06-22):** the naive `audit_103b` truncation null is
  **invalid** (δ negative-control false-positives, N2.5) and retired; the duration
  verdict is carried by the clean ratio regression + full-length matched-strength →
  **β inference arc is duration-robust; α duration-suspect.** N2.2 ships
  duration-controlled.
- [x] **Split decision (DONE 2026-06-22): N2 stays unified — do NOT spin out a 5th
  headline.** Inference is duration-controlled enough to be N2's confident climax,
  but too thin anatomically (no hotspot) and behaviorally (no data) to carry its own
  headline next to N1. Headline count stays **4**.
- [ ] ~~Behavioral TC1~~ **REMOVED** — TI performance data is unavailable
  (PI 2026-06-22). There is no brain–behavior test; do not re-scope one.
- [ ] **`task_learn` audit:** lengths, quality, and whether a `rest_pre →
  task_learn → rest_post` encoding-trace localizes to OFC (coordinate with N1 §G).
- [ ] Re-verify all cited inference CSVs against the locked N1 trace (the exact
  reproduction is a selling point — keep it true).

## §E — Figure / representation ideas

- **Four-phase schematic** making `task_learn` vs `task_test` (encoding vs
  inference) explicit — fixes the neglected-phase problem visually.
- Encoding→OFC vs inference→cingulate **double-dissociation brain** (with the
  honest "inference duration-downgraded" caption), full-length labelled; keep an
  **X-epi** variant, no C4.
- Band × content matrix (α encoding-only; β encoding+inference) driven from CSV.
- low-γ focal cingulate encoding within-region plot (the "averaging hides it"
  panel).

## §F — Provenance (CSV · script · timestamp)

- Arc (encoding echo + inference-specific) → `data/audit/consolidation_arc/{arc_per_patient.csv,
  arc_null_per_patient.csv}` · `audit_103_cophenetic_consolidation_arc.py` ·
  2026-06-12 / 2026-06-18. Scope+verdict
  `.agents/guides/task-persistence-investigation/2026-06-12_cophenetic-consolidation-arc.md`;
  memory `arc_inference_consolidation_2026_06_18`.
- **Stage-2-arc (RESOLVED 2026-06-22)** → `data/audit/consolidation_arc/arc_lenmatched_null_R200_{band}{,_center,_tail}.csv`
  · `audit_103b_arc_lenmatched_null.py`. **Verdict: truncation null INVALID** (δ
  negative-control false-positives p≈0.014 head+center, full p≈0.54) — retired.
  Duration verdict carried by: (a) per-patient length-ratio regression of
  `arc_per_patient.csv` `T_infspec_pe` (β ρ≈+0.25 p≈0.49; α ρ≈+0.53) — β
  duration-robust, α duration-suspect; (b) full-length matched-strength
  `arc_null_per_patient.csv` (β p≈0.007). Localization duration split:
  `lenmatched_null_R200_{include,exclude}.csv` (`audit_113b`).
- Inference localization → `data/audit/inference_localization/inference_mark{,_R1000,
  _shaftcollapsed}_{include,exclude}.csv` · `audit_110_inference_mark_localization.py`
  (+ `audit_111_gen_task_learn_r1000_surrogates.py`) · 2026-06-18.
- Absolute within-system (low-γ focal) → `data/audit/inference_localization/within_system_trace_{include,exclude}.csv`
  · `audit_112_within_system_trace.py` · 2026-06-18.
- Duration control (localization) → `data/audit/inference_localization/{length_control.csv,
  lenmatched_null_R200_{include,exclude}.csv}` · `audit_113_inference_length_control.py`
  + `audit_113b_inference_lenmatched_null.py` · 2026-06-19. Memory
  `inference_mark_localization_2026_06_18`. Handoff
  `.agents/reports/2026-06-18_inference-mark-handoff.md`.

## §G — Missing parts / open

- ~~Stage-2-arc result~~ **RESOLVED 2026-06-22** — β inference arc is
  duration-robust (N2.5); the gate is closed.
- ~~Whether inference deserves its own headline~~ **DECIDED 2026-06-22** — no; N2
  stays unified, count stays 4 (§D).
- **No behavioral anchor** — performance data unavailable (PI 2026-06-22); the
  inference claim rests on neural decomposition + anatomy + literature.
- **`task_learn`'s own trace** and length/quality audit (the duration asymmetry
  makes `task_learn`'s recording length per patient worth auditing).
- **No anatomical home for inference** — cingulate is a directional hint; whether a
  cleaner localizer (or a τ-sweep — diffusion timescale is unexplored for the arc)
  recovers one is open.
