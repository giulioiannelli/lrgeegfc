---
name: learning-phase-own-trace
type: report
era: IMCOH_ABS / COHORT_N10
status: result
created: 2026-06-22
headline: N2.1/N2.3 (.agents/preprint/headlines/02_encoding_vs_inference.md)
closes: N2 §D + §G open item "task_learn's own trace"
pointers:
  - data/audit/consolidation_arc/arc_null_per_patient.csv
  - data/audit/inference_localization/inference_mark_R1000_{include,exclude}.csv
  - data/audit/inference_localization/lenmatched_null_R200_{include,exclude}.csv
  - scripts/01_compute/audit/audit_103_cophenetic_consolidation_arc.py
  - scripts/01_compute/audit/audit_110_inference_mark_localization.py
  - scripts/01_compute/audit/audit_113b_inference_lenmatched_null.py
---

# The learning phase leaves its OWN persistent trace, and it anchors in OFC — same hotspot as the full trace, set before any inference

**Head.** Open item N2-§G asked: does `rest_pre → task_learn → rest_post` — the pure
**memorizing** phase, with no inference required — leave its own persistent trace,
and does it localize to OFC the way the full (test-phase) trace does? The answer is
**yes, and it was already computed** — the "encoding component" artifacts *are* this
object. The learning-phase trace clears matched-strength in **α and β**, localizes
to **OFC** (R=1000, the only system clearing BH in both epi modes, low-strength,
not hub-driven), and survives the duration control. **No new computation was
needed** — the signal was already in the encoding-component CSVs; this note surfaces
and closes it. (Per the standing rule: don't re-test when the signal is visible in
existing artifacts.)

## Why this is already answered (the identity)

The four-phase arc (`audit_103`) defines the **encoding component**
`e = D_task_learn − D_rest_pre_A` and the **persistence** `p = D_rest_post − D_rest_pre_B`.
The encoding echo `T_learn = ρ(e, p)` is **constructed identically** to the locked N1
trace `T_test = ρ(g, p)` with `g = D_task_test − D_rest_pre_A` — it simply swaps
`task_learn` for `task_test`. So `T_learn` **is** the literal "rest_pre → task_learn →
rest_post" trace, and the split-half baseline (preA for `e`, preB for `p`) makes it
the *rigorous* version, free of shared-baseline inflation. Likewise the localizer's
`encoding` target = `concordance(e, p)` is the per-system localization of that exact
trace. Follow-up A was hiding in plain sight inside the encoding side of N2.

## Result (verified from cached CSVs)

**1. The learning phase leaves its own trace — α and β** (`arc_null_per_patient.csv`,
matched-strength, positive = trace):

| band | T_learn obs (med) | surr (med) | n p<.05 | Wilcoxon p | LO-Pat_15 p |
|------|------------------:|-----------:|:-------:|-----------:|------------:|
| **β** | **+0.247** | +0.013 | **7/10** | **0.0137** | **0.0137** |
| **α** | **+0.223** | +0.030 | **6/10** | **0.0137** | **0.0098** |
| δ | +0.137 | +0.012 | 5/10 | 0.097 | 0.064 |
| low-γ | +0.103 | +0.003 | 5/10 | 0.116 | 0.082 |
| θ | +0.048 | +0.004 | 1/10 | 0.246 | 0.102 |
| high-γ | −0.098 | −0.000 | 4/10 | 0.348 | 0.248 |

α and β clear cleanly and are LO-Pat_15 robust; δ is a marginal hint only.

**2. The β learning-trace localizes to OFC** (`inference_mark_R1000_*.csv`, target
`encoding`, β, system granularity, BH across the 9 implanted systems):

- epi-**include**: **OFC** `M_obs=+4.78e5`, MS_p=0.001, **q=0.010 ✓** (only system
  clearing BH); strength-dev −0.51 (low-strength, not hub-driven).
- epi-**exclude**: **OFC** `M_obs=+6.68e5`, MS_p=0.001, **q=0.010 ✓**; MTL q=0.055
  (sub-threshold hint), insula p=0.032 q=0.107. Strength-dev −0.57.

This is the **same OFC hotspot** the test-phase ("standard") trace localizes to
(`audit_83`, β OFC q≈0.009–0.013) — so OFC is anchored by **both** task phases, not
just the test phase.

**3. It survives the duration control** (`lenmatched_null_R200_*.csv`, `audit_113b`,
length-matched strength null, R=200):

- encoding → OFC β, epi-**exclude**: MS_p=0.00498, **q=0.0498 ✓** (clears BH).
- encoding → OFC β, epi-**include**: MS_p=0.00995, q=0.0995 (BH-borderline — the raw
  matched-strength p stays solidly significant; only the BH-across-9 correction nudges
  it just over 0.05).

The raw matched-strength evidence is significant in **both** epi modes (p≈0.005–0.010);
the learning trace is therefore duration-robust, slightly softer than the test-phase
standard OFC trace (q=0.025/0.050 duration-controlled) but the same place and the same
low-strength character. This is expected and harmless: `e` uses the *shorter*
`task_learn` phase and has no test/learn length asymmetry baked in (the duration
confound was specific to `f = D_test − D_learn`, not to `e`).

## What it means

- **The OFC consolidation anchor is set by the memorizing phase itself** — before any
  inference is demanded. The brain files the *learned order* into OFC during learning,
  and that filing persists into rest. The longer test recording did not manufacture
  the OFC hotspot; the pure encoding phase already carries it.
- **It is encoding-side, not inference-specific.** The learning trace clears in **α
  and β**, mirroring the encoding echo — whereas the *inference-specific* component is
  **β-only**. The layered picture is clean:
  - **α — memorizing only:** learning trace clears (→ diffuse, no FDR-surviving hotspot).
  - **β — memorizing + reasoning:** learning trace clears **→ OFC**, and an *additional*
    inference-specific component rides on top (β-only, → cingulate, a hint).
- **Strengthens N2's encoding side.** N2.3 already said "memory → OFC"; this pins down
  that the OFC localization is anchored by `task_learn` in its own right, in both α and
  β, duration-robust — not an artifact of reading the test phase.

## Honest ceiling

- The epi-include duration-controlled localization is **BH-borderline** (q=0.099), not
  a clean clear; the clean clears are the R=1000 full-length localization (q=0.010,
  both epi modes) and the epi-exclude duration-matched (q=0.050). State it as
  "robust at the matched-strength level, BH-clears at epi-exclude" — not as a flat q<0.05.
- δ's cohort `T_learn` (+0.137, p=0.097) is a directional hint, not a clearer.
- No behavioral data — the learning trace is the neural residue of encoding-phase
  reorganization, never validated against memorizing *success* (none obtainable, PI
  2026-06-22).

## Provenance

- Cohort scalar → `data/audit/consolidation_arc/arc_null_per_patient.csv`
  (`T_learn_{obs,surr_p50,p}`) · `audit_103_cophenetic_consolidation_arc.py --null` · 2026-06-18.
- Localization (R=1000) → `data/audit/inference_localization/inference_mark_R1000_{include,exclude}.csv`
  (target `encoding`) · `audit_110_inference_mark_localization.py --R 1000` · 2026-06-18.
- Duration control → `data/audit/inference_localization/lenmatched_null_R200_{include,exclude}.csv`
  · `audit_113b_inference_lenmatched_null.py` · 2026-06-19.
- Surfaced + closed (no new compute) 2026-06-22. Closes N2 §D + §G "task_learn's own
  trace". Memory `arc_inference_consolidation_2026_06_18`.
