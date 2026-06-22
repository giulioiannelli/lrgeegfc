---
name: verify-stage2-arc-duration
era: IMCOH_ABS_COHORT_N10
status: RESOLVED_2026-06-22
kind: verification-brief
headline: N2.2/N2.5 (.agents/preprint/headlines/02_encoding_vs_inference.md)
owner_agent: inference chat
updated: 2026-06-22
---

# Verify — does the whole-brain inference-specific arc survive a length-matched null? (stage-2-arc)

> **✅ RESOLVED 2026-06-22 — β inference arc is DURATION-ROBUST.** The instrument
> below (`audit_103b` whole-brain truncation null) was run and **found INVALID**:
> its negative control fails — δ, dead-null at full length (p≈0.54), false-positives
> at p≈0.014 at both head and center windows; a contiguous truncation injects a
> common-mode shift into `D_test` that the un-demeaned whole-brain concordance reads
> as fake signal. **Retired.** The duration verdict is instead carried by (a) a
> per-patient length-ratio regression (β ρ≈+0.25 n.s. = no duration scaling; α
> ρ≈+0.53 = duration-suspect) and (b) the full-length matched-strength referee
> (β p≈0.007) → β duration-robust, α duration-suspect. See N2.5 in the headline.
> The 5-point design below is kept for the record.

**Head.** This brief **gated the inference half of N2**. The inference-specific
component `T_infspec·e` (β-only, `audit_103`) is built on `f = D_test − D_learn`,
and `task_test` is **1.35–2.53× longer** than `task_learn`. The duration control
already **downgraded inference→cingulate localization**; the *whole-brain* arc uses
the same `f`. `audit_103b` ran the length-matched nulls (2026-06-22) — see the
resolution banner above.

## 5-point critical preamble

1. **Claim.** `T_infspec·e` (partial corr of `f` with persistence `p`, controlling
   encoding `e`) persists in **β** even when `task_test` is truncated to
   `task_learn` length — i.e. the inference component is **not** a duration
   artifact.
2. **Null.** At matched length, `T_infspec·e` ≤ a **matched-strength surrogate
   regenerated on the truncated `task_test`** (R≥200; the `audit_103b` null).
3. **Strongest alternative.** **Duration confound** — longer `task_test` gives `f`
   more samples / more drift, inflating the partial correlation independent of
   inference content. (Band-specificity argues against a *generic* duration effect,
   but does not by itself close it.)
4. **Does the null control it.** Truncate `task_test` to per-patient `task_learn`
   length, **recompute `f` and the matched-strength null on the truncated data**,
   re-test `T_infspec·e`. This matches the sample/length budget of the two phases
   and asks whether the inference excess survives. Cross-check vs the localization
   analogue (`audit_113b`) which already split the verdict.
5. **Falsification + limits.** Falsified (downgraded) if β `T_infspec·e` drops
   inside its matched-length null (as inference→cingulate did, p≈0.12–0.18).
   Limits: truncation **reduces power** (report both full and matched); n=10;
   `task_learn` length/quality is the binding asymmetry — audit it.

## Steps

1. Per patient, truncate `task_test` to `task_learn` duration (matched window).
2. Recompute `f = D_test_trunc − D_learn` and persistence `p`.
3. Regenerate the **matched-strength null on truncated `task_test`** (R≥200; head
   already at R=200, optionally R=1000).
4. Re-test `T_infspec·e` cohort (one-sided Wilcoxon), LOO-max, **all bands**
   (must stay β-specific to support the claim).
5. Report full-length vs matched-length side by side.

## Data / scripts

- `data/audit/consolidation_arc/arc_lenmatched_null_R200_{alpha,beta,delta,theta,low_gamma}.csv`
  · `audit_103b_arc_lenmatched_null.py` · 2026-06-22 (**running** — check
  completion across bands before adjudicating).
- Baseline arc: `data/audit/consolidation_arc/{arc_per_patient.csv,arc_null_per_patient.csv}`
  · `audit_103_cophenetic_consolidation_arc.py`.
- Localization analogue (already split): `audit_113b_inference_lenmatched_null.py`;
  memory `inference_mark_localization_2026_06_18`.

## Pass / fail

- **Pass (β survives matched-length, stays β-specific):** N2.2 inference component
  upgrades from *provisional* → established; N2 can credibly split into an
  inference headline (3→4→5 path).
- **Fail (drops inside null):** downgrade the inference arc to a directional hint
  (mirroring inference→cingulate); N2 stays an encoding-led headline with inference
  as a hint.

## Caveats

- Truncation power loss — never read a null result as "disproven" without the
  power caveat; report effect sizes both ways.
- This is the single highest-leverage open control for N2.

## What to return

Verdict (survives / downgraded), full-vs-matched numbers (from CSV), band-
specificity check, LOO-max, and the updated N2.2/N2.5 status line. Fill below.

> **Result (agent fills):** …
