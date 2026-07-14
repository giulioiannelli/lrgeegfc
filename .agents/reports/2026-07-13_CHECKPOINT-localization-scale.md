---
name: checkpoint-2026-07-13-localization-scale
kind: verdict
era: IMCOH_ABS_COHORT_N10 (mst@0.20 recovery)
status: current
created: 2026-07-13
scope: THE resume anchor after compaction. Three results settled on mst@0.20 EXCEPT localization; the β trace is DELOCALIZED at every tested scale. New governing principle — everything is now tested against SCALE, localization included ("is there a scale where it holds / is spared?"). Open checklist + exact next steps to finalize.
pointers:
  - .agents/preprint/established_results/2026-07-13_settled-three-results.md   # the three results, locked
  - .agents/reports/2026-07-12_mst020-recovery-arc.md                          # recovery arc
  - scripts/01_compute/sparsified_arc/16_localization_mst020.py                # the localization port
---

# CHECKPOINT 2026-07-13 — localization + scale-everything

> ⛔ **SUPERSEDED (2026-07-13 PM) → resume from
> `.agents/reports/2026-07-13_CHECKPOINT-PM-narrative-and-cascade.md`.** All open items
> below (enc/inf localization, tissue, laterality, duration, Grassmann) are now DONE; the
> localization verdict here (β-trace-only, 3 scales) is extended by the full trace+enc+inf
> 16-scale sweep (all delocalized). Kept for history.
>
> 🔬 **AUDITED (2026-07-13 eve) → `.agents/reports/2026-07-13_localization-cohort-statistic-verdict.md`.**
> The "all delocalized" here was derived with a per-patient Wilcoxon that is DOA at K=5
> (broken cohort stat). A four-statistic bake-off (script 23) on the SAME backbone/null
> confirms **β genuinely delocalized on a VALID powered stat** (β→OFC stays retired), but
> shows the broken test had **hidden a real `low_gamma` concordance** — but low_γ is DIFFUSE
> (8/9 systems +, top = occipital K=3 artifact), NOT a focal cognitive home; the one specific
> low_γ result is **low_γ → SOZ** (disease). So: β delocalized ✓ (right reason now); low_γ
> diffuse + SOZ-only; NO band gives a focal cognitive localization. State PER BAND.

## Head

The three results are settled and verified on `mst@0.20` **except localization**. New
finding: **the β trace is DELOCALIZED at every scale** — "β → OFC" does **not** survive
the recovery (the whole-graph β→OFC was a degenerate-graph artifact, and the whole graph
is not a valid baseline). New governing rule, applied from now on: **every claim is tested
against scale (τ-swept)** — for each, ask *"at which scales does it hold, and is there a
scale where it is spared?"* Goal once localization resolves: **finalize everything.**

## Governing principle (NEW — apply to EVERY claim)

1. **The whole (dense) graph is NOT a baseline** — it is degenerate (audit_174). Test every
   claim independently on `mst@0.20`; never frame a result as "reproduces / contradicts the
   whole graph."
2. **Everything is τ-swept.** Report per-scale; ask *"which scales hold?"* This now includes
   **localization, tissue, SOZ-divergence, laterality, reinstatement** — not just the trace.
3. **No band is "false"** (the cophenetic result is not ground truth). Characterize bands by
   **scale-shape + anatomy**.
4. **Framing = information-DEPTH, not exclusive detection.** Simple pairwise FC = an
   info-poor scalar per band. Multiscale strengths that survive: ① **scale-signature**
   (β scale-invariant, α single-scale, δ fires-but-shallow), ③ **higher-order multi-step
   structure**. ② **localization is now WEAK** on the recovered scheme (see below).
5. **One null: matched-strength.** Read per-scale, never best-scale. Drift retired.

## Settled on mst@0.20 (verified) — full detail in `settled-three-results.md`

- **R1 trace:** β 16/16 scales (scale-invariant, p .014→.001); α 12/16 (mesoscale, .024→.007);
  θ/low-γ null; **δ weak+fragmented** (8/16 non-contiguous, fails fine+meso, rides coarse
  collapse — NOT clean null, NOT false); high-γ coarse-only. Pillar-6 laterality (β left,
  ρ.685 p.029) **only tested at s=5.6 → SCALE-TEST PENDING**.
- **R2 enc/inf:** encoding β 16/16 (scale-invariant) + α basic-metric-visible; inference
  (T_infspec_pe) δ/α/β **not β-only**; **δ-inference CUT** (partial-corr artifact).
  **LOCALIZATION of enc/inf components NOT yet tested on mst@0.20.**
- **R3 epi:** seed-affinity AUC δ.83/low-γ.82/β.745 (beats null 7-8/10); τ = ranking not
  precision (prec@5 .40 single, 60% is band-fusion); robust to sparsification. Solid.

## LOCALIZATION verdict (NEW 2026-07-13) — `16_localization_mst020.py`

TRACE localization (ρ_sym per-pair → per-system + SOZ, matched-strength R=200, BH/band),
scales s=1.0 / 2.83 / 5.65. Data: `data/sparsified_arc/localization_mst020_s{01.0,02.8,05.7}/`.

- **β = DELOCALIZED at every scale.** OFC +91k (p.50) @s1 → −322k (p.59) @meso; cingulate
  ~+150k (p.47) flat. **Nothing clears BH.** Strong but placeless — consistent with
  `node_ranksize_delocalisation` (β ~88% co-move).
- **δ:** diffuse, cingulate-leaning (+700k, **p.055** @meso), **AVOIDS the SOZ** (negative).
- **α:** weak insula (**p.031, 5/5** @s2.8 only; q.25) — scale-specific.
- **low-γ:** PFC across scales (p.02–.055, 6/7); SOZ @meso (p.053, 7/10).
- **NOTHING survives BH anywhere.** Port is validated (it *does* detect α→insula, low-γ→PFC
  trends → the β-null is a real delocalization, not a broken script).
- ⇒ **Localization IS weakly scale-dependent (a result), but no band has a BH-significant
  home on the recovered scheme.** "Strength ② (where it lives)" is weak; lean on ① + ③.

## OPEN CHECKLIST (all scale-aware) — before the "full picture"

**Localization (active thread):**
1. ⬜ **ENCODING (T_learn) component localization on mst@0.20, FULL scale sweep** — does
   enc→OFC survive at ANY scale? *Decisive for whether R2 keeps any anatomy.* (extend script 16)
2. ⬜ **INFERENCE (T_infspec_pe) component localization, scale-swept** — inf→cingulate at any scale?
3. ⬜ **Full 16-scale β localization sweep** — the user's explicit question: *"is there a scale
   the localization is spared?"* (I tested 3 scales; sweep all 16 to be certain.)

**Other whole-graph items → port to mst@0.20 + scale-test:**
4. ⬜ Tissue-class carrier (β gray–gray vs SOZ), scale-swept.
5. ⬜ SOZ-divergence (β spares / α recruits), scale-swept. [partial: δ avoids SOZ; low-γ→SOZ@meso]
6. ⬜ Reinstatement (held-not-replayed) — scale-test.
7. ⬜ Pillar-6 laterality — scale-test (currently s=5.6 only).
8. ⬜ Grassmann on mst@0.20, or drop from the story.

**Reframes (after localization resolves):**
9. ⬜ R1 §1.8 "β→OFC" → "delocalized" in settled doc + sec1 directive + talk slide 16.
10. ⬜ R2 enc/inf localization → per (1)(2) (OFC/cingulate may not survive).
11. ⬜ Finalize settled doc + 3 directives + slides once everything is scale-tested.

## How to run (env + scripts)

- Env python: `/home/giulio/Documents/miniconda3/envs/lapbrain/bin/python` (conda run fails).
- Pin `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`, ≤12 workers.
- Localization: `SA_SCALE=<s> $PY scripts/01_compute/sparsified_arc/16_localization_mst020.py`
  → `data/sparsified_arc/localization_mst020_s<NN.N>/`. ~35 s/scale. Currently tests the
  TRACE (T_test-like ρ_sym); needs extension to T_learn / T_infspec_pe for items 1–2.
- Trace gate: `13_matched_strength_mst020.py`. Enc/inf: `05_enc_inf_arc.py` (T_learn bug fixed).
  Controls ladder: `14_controls_ladder_mst020.py`. Epi: `06_epi_arc.py`. Pillar6: `15_patient_variation.py`.

## Git / state

- Last commit **73457fe** (directives + slides + T_learn fix + pillar6). **UNCOMMITTED:**
  `16_localization_mst020.py`, this checkpoint, settled-doc edits, memory. `data/` gitignored.
- Memory: `sparsified_recovery_mst020_2026_07_12.md` + MEMORY.md handoff (updated with loc verdict).

## Recommended next action

Run **item 1** (encoding component localization, full scale sweep) — it decides whether R2
keeps ANY anatomy. If enc→OFC and inf→cingulate are also delocalized at every scale, the
honest story loses its anatomy entirely and rests on **scale-signature + higher-order +
band-dissociation**; if they survive at some scale, localization becomes a *scale-resolved*
result. Either way, that run + items 3–7 complete the scale-tested picture → then finalize.
