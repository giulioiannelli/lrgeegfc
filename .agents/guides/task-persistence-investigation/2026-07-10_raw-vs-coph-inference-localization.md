---
name: raw-vs-coph-inference-localization
era: IMCOH_ABS_COHORT_N10
status: result
kind: task-persistence-investigation
scope: The drift-free separator on the inference-specific cell. Where raw edges TIE cophenetic on matched-strength detection (T_infspec_pe β = 0.0098 both) and drift is invalid (conditional functional), does the hierarchy CONCENTRATE the trace anatomically more than raw? Answer: yes — coph→cingulate BH-clears (q 0.03/0.04), LOO-robust; raw is placeless (no BH-significant system). Localization, not detection, is the multiscale value-add on this cell. audit_171.
created: 2026-07-10
updated: 2026-07-10
---

# Raw-vs-cophenetic localization of the β inference-specific trace

**Head.** On the one cell where the pairwise/multiscale convergence is genuinely
unresolved — the β inference-specific functional `T_infspec_pe`, where raw edges
tie cophenetic under matched-strength (both `p = 0.0098`) and the drift null is
*invalid* (conditional/partial correlation) — the cophenetic representation
**concentrates the trace on cingulate and BH-clears** (`q = 0.03` include /
`0.04` exclude, LOO-robust), while the **raw-edge representation is placeless**
(no anatomical system survives BH; its uncorrected lead is a *different* system,
sensorimotor). Detection is tied *by construction* (coph is a lossy deterministic
function of the same `W`); **anatomical concentration is not** — and that is the
drift-free axis on which the hierarchy beats raw edges here.

## Why this is the test (motivation)
- Matched-strength: raw *ties or beats* coph on β detection (`T_infspec_pe` raw
  0.0098 = coph 0.0098; eigcent 0.0020). Verified by fresh recompute (bit-identical
  to `cohort_summary.csv`, max|Δ|=0 over 60 cells). So no *detection* test can
  separate them — structurally, `ρ_coph = f(W)` carries ≤ the information in `W`.
- Drift: **closed** for the conditional functional — a windowed drift null returns
  a spurious positive on a no-signal pre-task arc ([[feedback_drift_null_mandatory]]).
- ⇒ The only remaining drift-free axis is a *non-detection* property. Localization
  (anatomical concentration) is that axis.

## Method (audit_171)
Identical to the audit_160/audit_110 inference-mark localizer, run for **two
representations** with everything else held fixed — only `φ(W)` changes:
- `coph`: `D[ph] = ` canonical cophenetic of `W[ph]`; surrogate = `cophenetic_from_eigs`.
- `raw` : `D[ph] = ` raw upper-tri edges of `W[ph]`; surrogate =
  `raw_edges(adjacency_from_laplacian_eigs(·))`.
Per-pair `concordance_partial(f, p | e)` (rho_sym arm-symmetrized) → demeaned
per-system `unit_means` → **R=1000 matched-strength** cached ensemble
(seed 20260511, the SAME null draws for both reps) → per-system upper-tail
`MS_p` → BH over systems → LOO. Anti-hallucination anchor: the coph side
reproduces audit_160's cached `inference_pe` localizer **to the decimal**
(cingulate `M_obs` 430046 / `MS_p` 0.002997 include; 478816 / 0.003996 exclude).

## Result (β inference_pe, R=1000)
| repr | epi | leading system | MS_p | LOO max | BH_q |
|---|---|---|---|---|---|
| **coph** | include | **cingulate** | 0.003 | 0.049 | **0.030** ✓ |
| **coph** | exclude | **cingulate** | 0.004 | 0.032 | **0.040** ✓ |
| raw | include | sensorimotor | 0.032 | 0.118 | 0.320 ✗ |
| raw | exclude | sensorimotor | 0.022 | 0.027 | 0.220 ✗ |

- coph → cingulate clears BH and holds under LOO, in **both** epi conditions.
- raw → no BH-significant system; its uncorrected lead (sensorimotor) dies under
  BH-over-systems and is a *different, non-consolidation* location. **They do not
  co-localize** — the key point (not "coph sharper at the same spot" but "coph
  finds a coherent spot, raw finds noise").

## Honest scope / caveats
1. **Not a detection claim.** Raw ties coph on the p-value and *must* (coph = f(W)).
   The claim is narrow: *coph localizes the inference trace, raw does not.*
2. **Raw is not literally flat.** Its sensorimotor lead is uncorrected-significant
   (0.022–0.032); the correct statement is "no BH-significant coherent hotspot,"
   not "sees nothing."
3. **Duration/length control PENDING.** This is matched-strength only. Coph's
   inference→cingulate was previously *demoted to a hint* (~0.124) by the
   length-matched control (audit_113b/161). Whether the coph-vs-raw *contrast*
   survives that control is untested — and, per the drift lesson, whether a
   windowed length-resampling null is even *valid* for a conditional functional is
   an open question the user must rule on before it is run. Until then, the
   separator stands **under matched-strength**.
4. **cingulate ≠ OFC.** This is the inference-mark (R2.5/R2.6), distinct from the
   whole-task β→OFC flagship.

## Connection to prior tools
- Sibling of [[2026-07-09_pairwise-descriptor-ladder]] (audit_166): the ladder shows
  β detection converges; this shows localization does *not* converge.
- Reuses audit_110/audit_160 inference-mark localizer verbatim; audit_83 primitives.

## Open questions
- Does the contrast survive the length-matched null (if that null is valid here)?
- Is the same localization separation present for `encoding` / `standard` targets
  (where drift already separates raw from coph, so it is corroborative not load-bearing)?

Build: `scripts/01_compute/audit/audit_171_raw_vs_coph_inference_localization.py`.
Data: `data/audit/raw_vs_coph_inference_localization.csv`.
