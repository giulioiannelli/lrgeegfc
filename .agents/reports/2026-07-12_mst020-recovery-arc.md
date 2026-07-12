---
name: 2026-07-12_mst020-recovery-arc
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
status: current
created: 2026-07-12
updated: 2026-07-12
pointers:
  - .agents/reports/2026-07-12_brutal-review-failure-directions.md
  - memory/sparsified_recovery_mst020_2026_07_12.md
  - data/sparsified_arc/{sparsification_recovery,ms_mst020,enc_inf_arc_mst020,epi_arc_mst020,controls_ladder}
---

# mst@0.20 recovery arc — the multiscale LRG trace, rebuilt to survive

**Head.** After the drift cascade + brutal review broke the story, the fix was to
**pick the sparsifier, not chase p-values**. On a minimal-modification density backbone
(`mst_union_top_fraction(W, 0.20)` = MST ∪ strongest-20% edges) the old audit_150 α/β
cophenetic trace is **reproduced and sharpened** under matched-strength, now τ-resolved
and genuinely multiscale, with clean band-selectivity (θ/low_γ null; δ/γ_high are
scale-collapse artifacts, not claimed). **Drift is discarded** as a null (the task is
directional → a trace is monotonic by construction → drift ≡ the alternative).
Matched-strength is the only meaningful null. The propagator additionally serves as a
competitive SOZ *marker* (δ.83/γl.82/β.745) though not an improved *detector*. Encoding
vs inference dissociates, but the sharp "β-only inference" claim weakened; a
representation-ladder reframe is in progress.

## 1. Why mst@0.20 (minimal-modification sweep)
`12_sparsification_recovery.py`, `data/sparsified_arc/sparsification_recovery/`.
Swept `mst_union_top_fraction` over frac {1.0,…,0.05} + percolation + TMFG × 16 scales
`s ≡ τλmax ∈ [1,180]`. Findings:
- **dense @ s=1 == audit_150 exactly** (β .198 ≫ α .101 > rest) — the anchor.
- **Multiscale is a sparsification effect**: specific-heat C(τ) peak count = dense 1.0
  (degenerate, audit_174) → mst0.15–0.20 / percolation 1.4–1.5 (multiscale ladder turns
  on **below density ~0.2**).
- **β is sparsification-invariant** (ρ_sym 0.12–0.29 at every scale, dense→mst→perc→tmfg);
  **α is scale-tuned** (a mesoscale peak) and needs the **strong-edge** backbone —
  mst@0.20 keeps α (0.18) while percolation at the SAME density kills it (0.02). Edge
  selection matters, not just density. ⇒ mst@0.20 is the minimal departure that keeps
  both α and β AND switches multiscale on.

## 2. Trace verdict (matched-strength, τ-resolved)
`13_matched_strength_mst020.py`, R=200, `data/sparsified_arc/ms_mst020/`. Pipeline-
identical null: shuffle dense FC (audit_150 4-cycle) → mst@0.20 → LRG cophenetic at
scale s → ρ_sym. Read **per-scale** (never scale-max — the scale-max table falsely
flagged δ/γ_high CLEAR).

| band | scales clearing p<0.05 (of 16) | cohort p @ s=5.6 | verdict |
|---|---|---|---|
| **β** | **16/16**, s∈[1,180] | **0.001** | scale-invariant trace |
| **α** | **12/16**, s∈[1,90] | **0.007** | mesoscale trace (=old .024 @ s=1) |
| δ | 8/16 non-contiguous | 0.116 | patchy, 2–6 pt — not claimed |
| high_γ | 4/16, coarse only (s≥45) | 0.161 | null-collapse artifact — not claimed |
| θ | 0/16 | 0.577 | clean null |
| low_γ | 0/16 | 0.188 | clean null |

Per-patient @ s=5.6: **β all 10 positive, 7 beat own p95**; dissenters (Pat_13/14/15)
near-zero, none negative; **Pat_15** (β≈0, α<0) = known right-hemisphere-only implant;
**Pat_06** super-responder (δ/θ/α 0.78/0.87/0.90); Pat_13 globally low. Band-selectivity
is RESTORED — the percolation "all-bands-fire" (D2) was a percolation-specific pathology.

## 3. Encoding vs inference (`05_enc_inf_arc.py` SA_BACKBONE=mst020)
`data/sparsified_arc/enc_inf_arc_mst020/`. Per-scale MS gate on T_test/T_learn/T_infspec/
T_infspec_pe (formulas: e=D_L−D_A, f=D_T−D_L, g=D_T−D_A, p=D_P−D_B; T_x=ρ_sym of the pair).
- **T_test reproduces the trace** (β 16/16, α 12/16). Consistency ✓.
- **T_learn (encoding) is null in every band** (0/16) — obs high (0.4–0.6) but strength-
  matched → encoding reorganization is strength-explained.
- **T_infspec_pe** clears **δ/α/β** — NOT β-only (the full-graph β-only did not survive).
- **T_infspec (raw)** clears **β only**, fine scale s∈[1.4,4] — β is the sole standalone
  inference-specific persistence.
- Reframe in progress (`14_controls_ladder`): show T_learn is readable from a basic
  network metric (encoding = stored in simple connectivity) while T_infspec is readable
  only from the multiscale cophenetic (inference = multiscale).

## 4. Epilepsy marker/detector (`06_epi_arc.py` SA_BACKBONE=mst020)
`data/sparsified_arc/epi_arc_mst020/`. Seeded heat-kernel marker, leave-self-out,
strength-residualised, rest_post, strength-matched fake-SOZ null R=200.
- **MARKER (AUC)**: multiscale δ.83 / γl.82 / β.745 (best-scale .87/.85/.90), beats null
  **7–8/10 pt**; multiscale > single-scale in 7–8/10 for δ/β/γl. Matches/beats audit_132.
- **DETECTOR (precision)**: prec@5 ≈ 0.40 (δ/γl/γh; 4.4× the 0.09 prevalence) but
  **identical to single-scale** and **below the prior audit_117 0.60**. The τ-sweep helps
  ranking, **not** top-k precision. Occult candidates in `occult_candidates.csv`.
- **Band dissociation (talk hook)**: trace = α/β; epi = δ/γl/β. One operator, two
  band-specific readouts; β the only overlap.

## 5. Open / next
- **Controls ladder** (`14_controls_ladder_mst020.py`, running): raw-FC/strength/
  clustering/geodesic/**resistance(spectral)**/coph@τ_min/coph@meso × 4 functionals × MS.
  Claims: (A) band-selective trace UNIQUE to cophenetic; (B) encoding readable from a
  basic metric, inference only from multiscale cophenetic.
- **Patient variation ~ electrodes** (pillar 6): explain Pat_15/06/13 by implant geometry.
- Specific-heat C(τ) curves not yet plotted per-cell (only peak-count summarized).
- **Not migrated to preprint** — all under `data/sparsified_arc/`; talk-facing.

## Provenance
Env python `/home/giulio/Documents/miniconda3/envs/lapbrain/bin/python`; pin
OMP/OPENBLAS/MKL=1 + SA_WORKERS≤12 (LAPACK oversubscription). Scripts 12/13/14 +
05/06 (SA_BACKBONE=mst020) + fig_* under `scripts/01_compute/sparsified_arc/`.
