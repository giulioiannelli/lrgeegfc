---
name: consolidation-arc-handoff
type: handoff
era: IMCOH_ABS / COHORT_N10
status: active
created: 2026-06-18
updated: 2026-06-18
pointers:
  - .agents/guides/task-persistence-investigation/2026-06-12_cophenetic-consolidation-arc.md
  - .agents/reports/2026-06-12_trace-impact-and-leverage.md
  - scripts/01_compute/audit/audit_103_cophenetic_consolidation_arc.py
---

# RESUME-HERE handoff — "what does the brain keep after the task: memory or inference?"

**To restart after a compact: read this file top-to-bottom, then continue from
"NEXT STEP". Everything below is in plain words first; exact file paths and the
one command to run are at the bottom.**

---

## In plain words — what we did

- The patients (epilepsy, electrodes in the brain) did a reasoning task. First
  they **learn an order** by seeing pairs ("A beats B", "B beats C", …). Then
  they are **tested on pairs they never saw** ("does A beat C?") and have to
  **work it out**. This is called transitive inference.
- We already knew the task leaves a **lasting mark** in how brain regions are
  organised — it is still there during the rest period *after* the task. It is
  strongest in the **beta brain-rhythm** and sits mainly in the **orbitofrontal
  cortex** (a region known for building "mental maps").
- The complaint we were solving: that lasting mark felt like a clever
  measurement, not a real brain discovery. So we asked a sharper question.
- There are four recordings per patient: **rest → learning → testing → rest**.
  We used all four and asked: when the brain keeps a lasting mark afterwards, is
  it holding onto **the pairs it was shown** (memorising) or **the new
  relationships it figured out** (inferring)?

## In plain words — what we discovered

- The brain keeps a lasting mark of **both**, but **split by rhythm**:
  - the **alpha** rhythm keeps only the **memorised pairs**;
  - the **beta** rhythm keeps the memorised pairs **and** the **figured-out
    relationships**.
- Said simply: **the lasting trace of "what the patient worked out" (not just
  what they saw) lives specifically in the beta rhythm.**
- This passed the strict statistics test that is mandatory here (the one that
  has killed earlier claims): beta result p ≈ 0.007, still holds when we drop
  the one odd patient, and shows up **only in beta** — which rules out the main
  boring explanation (one recording being longer than another would affect all
  rhythms, not just beta).
- We proved the code is trustworthy: it reproduces the previous locked result
  **exactly** (to 15 decimal places, and the cohort numbers match for all six
  rhythms).

## Honest limits (do not oversell)

- The **memorising** part is the bigger effect; the **inferring** part is
  smaller but real and beta-specific.
- We have **not yet checked where** the beta "inference" mark sits — is it in the
  orbitofrontal "map" region too? (that is the next step).
- We do **not yet have the patients' task scores**, so we cannot yet say
  "stronger mark = better reasoning". That is the big future test.

---

## NEXT STEP (recommended, runs on existing data — no behaviour needed)

**Find out WHERE the beta "inference" mark sits.** Re-run the existing
localisation analysis on the *inference-specific* signal instead of the whole
trace. If it concentrates in the **orbitofrontal cortex** like the main trace
does, the story fuses into one headline: *the orbitofrontal cortex consolidates
the inferred structure offline.*

- Reuse `audit_83`/`audit_92` localisation machinery, but feed it the
  **inference-specific per-pair vector** `f = D_task_test − D_task_learn`
  (cophenetic distance of test minus learn), instead of the usual trace vector.
- First verify the exact plug-in point in `audit_83`'s localisation path before
  writing code (do NOT assume the API — we have been burned by that).

### Other queued items (all no-behaviour, do after the localisation)
1. **Length control:** the testing recording is 1.4–2.5× longer than learning.
   Cut testing down to the learning length, recompute, re-test the beta result.
   (Band-specificity already argues the confound is not the cause, but this
   closes it cleanly.)
2. **"Replay" moments:** look for brief instants during the after-rest where the
   brain momentarily snaps back into the task pattern (a separate, bigger idea —
   the sliding-window connectivity infra exists, but per-window imaginary-
   coherence stability must be checked first).
3. **Behaviour (future, needs scores):** does a stronger beta inference-mark
   predict better reasoning performance? This is the headline brain-behaviour
   result; needs the obtainable accuracy/reaction-time data.

---

## Exact state — what is DONE vs PENDING

DONE & validated:
- `scripts/01_compute/audit/audit_103_cophenetic_consolidation_arc.py`
  - observed pass: `python …/audit_103_cophenetic_consolidation_arc.py`
  - matched-strength null: `… audit_103_cophenetic_consolidation_arc.py --null`
  - (run with `/home/giulio/Documents/miniconda3/envs/lapbrain/bin/python`;
    `conda run -n lapbrain` is broken by a cross-compiler activation hook.)
- Outputs: `data/audit/consolidation_arc/arc_per_patient.csv`,
  `arc_null_per_patient.csv`, `null_run.log`.
- Scope + 5-point preamble + VERDICT:
  `.agents/guides/task-persistence-investigation/2026-06-12_cophenetic-consolidation-arc.md`.
- Result memory: `memory/arc_inference_consolidation_2026_06_18.md`.
- Deep-research literature + impact + leverage (clinical / methods /
  consolidation), with citations: `.agents/reports/2026-06-12_trace-impact-and-leverage.md`.

KEY NUMBERS (matched-strength null, n=10, one-sided paired Wilcoxon; LO-Pat_15 in parens):
- `T_learn` (encoding echo): α p=0.014 (0.010) ✓ ; β p=0.014 (0.014) ✓.
- `T_infspec·e` (inference-specific, controlling encoding): **β p=0.0068 (0.014), 7/10 ✓**;
  α p=0.080 ✗; all other bands null.
- Validation: observed `T_test` = audit_83 per-patient to 1.1e-16; null `T_test`
  cohort = audit_83 cached exactly, all 6 bands.

PENDING:
- localisation of `f` (NEXT STEP); length-equalised control; replay-events route;
  behavioural TC1 (needs scores).

DEFINITIONS (so the code reads cleanly):
- `D_x` = LRG cophenetic distance vector of phase x (one number per region-pair).
- encoding `e = D_task_learn − D_rest_pre_A`; inference-specific
  `f = D_task_test − D_task_learn`; persistent `p = D_rest_post − D_rest_pre_B`.
- `T_learn = corr(e,p)`, `T_infspec·e = partial corr(f,p controlling e)`,
  positive = the task-change stuck (a "trace").
