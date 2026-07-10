---
name: supp-drift-controls
era: IMCOH_ABS_COHORT_N10
status: current
kind: supplementary
scope: Drift controls for the cross-phase trace. The locked C2 drift-floor covers the WHOLE-TASK trace (α/β/γ_low pass; β→OFC drift-clean) and raw-FC whole-task traces are entirely drift. IMPORTANT METHODOLOGICAL RECORD: a windowed drift null is INVALID for the four-phase CONDITIONAL functionals (T_infspec·e, T_learn) — it returns spurious positives on a no-signal pre-task arc; a brief 2026-07-10 "inference-specific is drift-confounded" verdict is WITHDRAWN. The inference-specific trace is validated by matched-strength, not a drift null.
updated: 2026-07-10
---

# S2 — Drift controls for the cross-phase trace

**Head.** Drift (slow within-recording nonstationarity making any later-vs-earlier
contrast correlate) is an orthogonal confound to matched-strength and needs its own
null **for the whole-task (simple-correlation) trace** — where the locked **C2
drift-floor** provides it and **passes (α 0.0068, β 0.0137, γ_low 0.0098)**; the raw-FC
whole-task correlation, by contrast, is entirely drift, which is exactly why the LRG
cophenetic trace (not raw FC) is the load-bearing probe. **A separate, cautionary
result:** we briefly tried to extend a *windowed* drift null (`audit_167`/`audit_168`)
to the four-phase **conditional** functionals (`T_infspec·e` inference-specific,
`T_learn` encoding) and read the output as "the inference-specific decomposition is
drift-confounded." **That reading was WRONG and is withdrawn** — a windowed drift arc
is not a valid null for a *partial* correlation. The inference-specific trace stands on
**matched-strength** (β-only p=0.0098) + duration-robustness; β→OFC whole-task stays
drift-clean via C2.

## Two orthogonal confounds; drift ≠ matched-strength
Matched-strength (C3) preserves node strength and rewires edges but leaves temporal
structure untouched → it does **not** control session **drift**. Drift needs its own
null **for a simple later-vs-earlier correlation** (the whole-task trace). This is what
C2 does.

## C2 (locked) — whole-task drift floor [VALID, authoritative]
`ρ_drift = ρ_S(D_coph^{preB} − D_coph^{preA}, D_coph^{postB} − D_coph^{postA})`;
gate `wilcoxon_split_gt_drift_p < 0.05`. Full-duration halves ⇒ well-powered.
Source `data/audit/ctm_triangle/cohort_summary.csv`.

| band | C2 `split>drift` p | verdict |
|---|---|---|
| alpha | **0.0068** | PASS |
| beta | **0.0137** | PASS |
| low_gamma | **0.0098** | PASS |
| delta | 0.246 | fail |
| theta | 0.278 | fail |
| high_gamma | 0.423 | fail |

## audit_167 — whole-task second construction [VALID]
Window a single **rest_pre** recording (pre-task ⇒ any positive whole-task trace is
drift by construction) into 5 ordered segments; run the ρ_sym arc for raw + cophenetic.
For the **whole-task** functional (a simple correlation, drift is a genuine confound):
β cophenetic +0.198 vs drift +0.051, **p=0.042 → drift-robust** (agrees with C2);
α/γ_low pass C2 but fail here (construction-dependent — C2 is the locked authority);
**raw edges = entirely drift** (raw whole-task drift ≈ real). Caveat: 1/5-length
windows are noisier than C2's full-duration halves → where they disagree, C2 wins.

## Why the windowed drift null does NOT apply to the conditional functionals [the withdrawn overreach]
`T_infspec·e = partial_ρ(f, p | e)` with `f = D_TT − D_TL`, `p = D_RP − D_br`,
`e = D_TL − D_bt`, is a **partial** correlation that already **conditions on encoding
`e`**. Drift-removal is therefore entangled with the estimator itself — there is **no
valid external windowed-drift-null analog** for a conditional trace (unlike a simple
correlation, which C2 handles). Direct evidence the *null* is invalid, not the trace:

- `audit_167` (equal windows) and `audit_168` (windows sized ∝ real phase durations
  `[L_pre/2, L_pre/2, L_TL, L_TT, L_RP]`) both build the sham arc **entirely from
  pre-task `rest_pre`, where no consolidation signal can exist**.
- Yet the conditional sham returns a **significantly positive** trace *on its own* in
  every band (`drift>0` p ≤ 0.014; β 0.007) and even *exceeds* the real value
  (β real +0.091 < sham +0.243). **A valid null must sit near zero on a no-signal
  arc.** A positive, significant sham ⇒ the windowed construction manufactures spurious
  *partial* correlation; it is measuring an estimator/construction artifact, not drift.

So the earlier verdict ("inference-specific is drift-confounded; drop N2's
encoding-vs-inference split") is **WITHDRAWN**. The reference numbers are retained below
purely as the cautionary record.

| band (audit_168, coph) | conditional real | conditional sham | sham>0 p |
|---|---|---|---|
| **beta** | +0.091 | +0.243 | **0.007** |
| alpha | −0.005 | +0.074 | 0.007 |
| low_gamma | −0.010 | +0.156 | 0.005 |
| delta | +0.036 | +0.103 | 0.001 |
| high_gamma | +0.023 | +0.159 | 0.014 |
| theta | −0.045 | +0.089 | 0.005 |

(The sham being significantly positive *everywhere* — including on a pre-task-only arc —
is the signature of an invalid null.) Source `data/audit/duration_matched_drift_infspec/`.

## Reading
- **Whole-task trace = drift-controlled.** β→OFC survives C2 (0.0137) + audit_167
  (0.042); raw whole-task FC is drift, the LRG cophenetic is not — the core reason to
  work in the cophenetic space. α/γ_low whole-task drift-robustness is
  construction-dependent (flag, not overturned).
- **Conditional (inference-specific / encoding) trace: NOT drift-testable this way.**
  A windowed drift arc is an invalid null for a partial correlation. Validate the
  conditional trace with **matched-strength** (β-only `T_infspec·e` p=0.0098, which
  *does* return ≈0 on no-signal), not a drift window. N2's inference-specific claim
  **stands**.
- **General rule (locked in memory):** whole-task simple correlation → windowed drift
  null is OK; conditional/partial trace → it is NOT — use matched-strength.

## Status
- **β whole-task → OFC:** CONFIRMED (matched-strength + C2 + audit_167).
- **N2 inference-specific decomposition:** stands on matched-strength + duration; the
  drift-null attempt is closed as a methodological cautionary record (not a verdict).
- **α whole-task:** drift-robustness construction-dependent — flag, not overturned.

Build: `audit_167_drift_null_ladder.py` (equal windows) +
`audit_168_duration_matched_drift_infspec.py` (duration-matched) — retained as the
cautionary record; data `data/audit/drift_null_ladder/`,
`data/audit/duration_matched_drift_infspec/`; C2 `data/audit/ctm_triangle/`.
