---
name: inference-mark-handoff
type: handoff
era: IMCOH_ABS / COHORT_N10
status: active
created: 2026-06-18
updated: 2026-06-18
pointers:
  - .agents/reports/2026-06-18_inference-mark-localization.md
  - .agents/guides/task-persistence-investigation/2026-06-18_inference-mark-localization.md
  - scripts/01_compute/audit/audit_110_inference_mark_localization.py
  - scripts/01_compute/audit/audit_112_within_system_trace.py
---

# RESUME-HERE — "where does the offline mark sit, and is α/low-γ locally strong?"

**To restart after a compact: read this top-to-bottom, then continue from
"NEXT STEP". Everything below is plain words first; exact commands at the bottom.**

---

## In plain words — what we did and found (this is DONE + saved)

The task is transitive inference: patients are **shown** premise pairs (learn),
then **tested** on novel pairs they must **infer**. Four recordings: rest →
learn → test → rest. We had already shown the β rhythm keeps an offline mark of
the *inference* part beyond the *memorised* part (the consolidation arc,
audit_103). This session answered **where that mark sits** and **whether the
weaker bands are locally strong**.

- **β splits anatomically by content.** The **memorised** material (encoding) and
  the standard trace over-accumulate in **orbitofrontal cortex (OFC)**; the
  **inferred** relations (inference, with memorising statistically removed)
  over-accumulate in the **CINGULATE**. At R=1000 both clear the strict
  multi-region correction: OFC q=0.010, cingulate q=0.050 (borderline q=0.060
  without epileptic contacts). *(An earlier R=200 run called the cingulate
  "suggestive q≈0.15" — that was a resolution-floor artifact; R=1000 fixed it.)*
- **Double dissociation (the strong evidence).** The **cingulate** carries
  *memorising* in α and low-γ, but flips to *inference* in β — independently
  re-deriving the arc's band split from anatomy alone.
- **α / low-γ question (user's hypothesis): confirmed for low-γ.** Whole-brain
  trace is significant only in α and β; **low-γ has no whole-brain trace** — yet
  an **absolute** within-region test (audit_112, no demeaning) finds a **genuine
  focal cingulate *memorising* trace in low-γ** (ρ=+0.40, 8/8 patients, q=0.035;
  survives removing epileptic contacts: ρ=+0.51, 7/8). A real trace that
  whole-brain averaging hides. α is *clear* overall and concentrated in
  temporal-limbic cortex — not murky.
- The recurring hub across all of it is the **cingulate**.
- **⚠ Duration control SPLITS the verdict (audit_113 + audit_113b, 2026-06-19).**
  `task_test` is 1.35–2.53× longer than `task_learn` in every patient (confound real,
  NOT moot). The matched-length strength null (audit_113b) splits the two sides:
  **memorizing → OFC SURVIVES** (standard OFC q=0.025/0.050; encoding ⊥`task_test`),
  but **inference → cingulate does NOT** — at matched length it sits inside its
  strength null (p=0.124 incl / 0.179 excl, q≥0.42) and isn't even the top inference
  system. The cingulate stays the *directional* home (stage-1: 8/8 sign, ρ=0.74) but
  its significance was duration-assisted → **inference → cingulate is downgraded to a
  directional hint, not an established localization.** ⚠ The consolidation ARC
  (audit_103) uses the same `f = D_TT − D_TL` and is NOT yet duration-tested.

## Honest limits (carry these forward)

- β cingulate **inference** clears BH only at the boundary (q=0.050/0.060); the
  corroboration is the 4-condition consistency + the double dissociation.
- The **absolute** within-system test (audit_112) is **underpowered for systems
  sampled by ≤5 patients** — OFC has a large trace (ρ=+0.56) it literally cannot
  prove at K=5. So "only cingulate clears the absolute test" is partly a sampling
  fact (cingulate K=8 best-sampled). **OFC's evidence is the demeaned
  localization (audit_110), not audit_112.** A demeaned *concentration* ≠ an
  absolute trace — report both.
- **No behavior yet** → "inference mark" = the offline residue of the
  inference-specific *reorganization*, not inference *success*.
- FRAMING (locked): a localization is an **overexpression on a distributed
  trace**, never a container ([[feedback-localization-is-overexpression-not-container]]).

---

## NEXT STEP (pick one; all run on existing data except #1)

1. **Behavioral TC1 (highest value, BLOCKED on data).** Does per-patient β
   cingulate inference-mark strength predict inference performance (symbolic-
   distance effect)? Needs the obtainable accuracy/RT scores — not in repo yet.
   When scores arrive: take per-patient β `inference_pe` cingulate value
   (recompute via audit_110 obs, no surrogate) and correlate with performance.
2. **Low-γ focal-trace standalone writeup (can do now).** The low-γ cingulate
   *encoding* focal trace (no whole-brain trace) is a result in its own right.
   Data ready in `within_system_trace_include.csv` / `_exclude.csv`. Would need a
   dedicated within-region figure + short report; decide framing with PI first.
3. **✅ DONE — duration control, BOTH stages (audit_113 + audit_113b, 2026-06-19).**
   `task_test` 1.35–2.53× LONGER than `task_learn` (NOT moot). Stage 1 (audit_113,
   observed, 3 placements): cingulate inference DIRECTION robust (8/8 sign, ρ=0.74,
   median 0.61× full; matched-vs-full paired p≈0.95 but underpowered K=8 = "didn't
   drop" ≠ "beats null"). Stage 2 (audit_113b, matched-strength null regen ON
   truncated `task_test`, R=200 head) = DECISIVE: verdict SPLITS — **standard→OFC
   survives** (p=0.005, q=0.025/0.050), **encoding→OFC** duration-invariant
   (⊥`task_test`, full q=0.010 stands), **inference→CINGULATE does NOT survive**
   (p=0.124/0.179, q=0.42/0.76; ranks 3rd/2nd, insula+MTL above; no system clears BH
   for inference). ⇒ inference→cingulate DOWNGRADED to directional hint; memorizing→OFC
   stays robust. Full-recompute==cache ρ=1.0. CSVs `length_control.csv`,
   `lenmatched_null_R200_{include,exclude}.csv`; cascaded into the report (head, title,
   Duration §, honest-limits, next-steps), memory, MEMORY.md, arc memory.
4. **⚠ NEW TOP CONTROL — stage-2-arc.** The arc (audit_103) uses the same
   `f = D_TT − D_TL`; duration-test the whole-brain `T_infspec·e` (β-vs-α dissociation)
   with a length-matched null — reuse audit_113b's truncated-`task_test` surrogate.
   Until then the arc headline is duration-uncontrolled.
5. Optional: R=1000 audit_113b (sharpens, cannot change the negative). Figures still
   show the full-length inference→cingulate BH ring — regenerate/annotate "full length".

---

## Exact state — DONE vs files

DONE & saved:
- `audit_110_inference_mark_localization.py` — demeaned localization (reuses the
  audit_83 OFC localizer verbatim; targets standard/encoding/inference_raw/
  inference_pe). `audit_110b_*` figures. `audit_111_*` generated the β
  task_learn R=1000 surrogate. `audit_112_within_system_trace.py` — absolute
  within-system trace.
- Outputs `data/audit/inference_localization/`:
  `inference_mark{,_R1000,_shaftcollapsed}_{include,exclude}.csv`,
  `within_system_trace_{include,exclude}.csv`,
  `figures/{brain_inference_dissociation_beta,systems_inference_dissociation}.pdf`.
- Report `.agents/reports/2026-06-18_inference-mark-localization.md`;
  scope `.agents/guides/task-persistence-investigation/2026-06-18_inference-mark-localization.md`;
  memory [[inference_mark_localization_2026_06_18]].

GOTCHAS (cost time last session):
- Run with `/home/giulio/Documents/miniconda3/envs/lapbrain/bin/python` —
  `conda run -n lapbrain` is broken (cross-compiler hook).
- **Two `audit_83` files (name collision):** the localizer is
  `audit_83_localization_matched_strength.py` (has `obs_trace`, `concordance`,
  `unit_means_from_s`, `--shaft-collapse`, `--R`); the WM one is
  `audit_83_wm_stratified_cophenetic.py`. Don't confuse them.
- **R=1000 is β-only** (audit_92 + audit_111 generated β only). α/low-γ
  localization is R=200.
- Background logs block-buffer when redirected to a file — check cache files /
  process state for progress, not the empty log.
- The localizer colours by **system** effect, not per-contact (single-contact
  cophenetic is too noisy — audit_76/108/109 lesson).

## Run commands

```
PY=/home/giulio/Documents/miniconda3/envs/lapbrain/bin/python
# demeaned localization (β at R=1000; α/low-γ at R=200; add --shaft-collapse / --epi-mode)
$PY scripts/01_compute/audit/audit_110_inference_mark_localization.py --band beta --R 1000
$PY scripts/01_compute/audit/audit_110_inference_mark_localization.py --R 200   # all bands, both epi
# absolute within-system trace (the audit_112 confirmation)
$PY scripts/01_compute/audit/audit_112_within_system_trace.py --R 200
# figures
$PY scripts/01_compute/audit/audit_110b_inference_localization_figure.py
```

## Key numbers (final)

- β demeaned (R=1000, BH/9 systems): standard→OFC q=0.010; encoding→OFC q=0.010;
  inference_pe→cingulate q=0.050 incl / 0.060 excl. Cingulate K=8, str_dev +0.04,
  5/8 pos, LOO-robust, shaft-robust.
- Whole-brain cohort trace (arc null): α std p=0.005 / enc p=0.014; β all ≤0.014;
  **low-γ null** (0.12 / 0.12 / anti).
- Within-system absolute (audit_112): low-γ encoding cingulate ρ=+0.40 8/8
  q=0.035 (incl) / ρ=+0.51 7/8 q=0.070 (excl). Power floor: K≤5 can't reach BH.
- Validation: standard OFC M_obs=+6.133e5 == audit_83 exactly; inference_pe
  partial == arc formula partial to 5.6e-17.
