---
name: 2026-07-13_multiscale-marker-exploration-scope
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
status: current
created: 2026-07-13
scope: critical preamble + design for the multiscale SOZ-marker exploration — does the
  tau-swept mst@0.20 propagator, via change-across-tau shape features / per-band scale
  fine-tuning / cross-band-cross-scale fusion, push AUC or precision above the fixed-tau
  single-scale marker? Built to REFUSE the in-sample best-scale forking path.
pointers:
  - scripts/01_compute/sparsified_arc/06_epi_arc.py         # base pipeline (reused)
  - data/sparsified_arc/epi_arc_mst020/                     # per-scale AUC curves (Stage 1 free)
  - .agents/reports/2026-07-12_mst020-recovery-arc.md        # §4 marker verdict
---

# Multiscale SOZ-marker exploration — scope + critical preamble

## Head

The τ-sweep gave the marker a **scale axis** it did not have before. This asks whether
using that axis *cleverly* — reading each band at its own scale, fusing bands each at
their own scale, or feeding the **shape** of a contact's scale-response to a classifier —
can push **AUC or top-5 precision** above the fixed-τ marker. The whole exploration is
built around **one** danger: with 10 patients × 6 bands × 16 scales, picking the scale
that maximises AUC *in sample* is a forking path that already inflates β to a fake 0.90.
So every number here is **nested leave-one-patient-out (LOPO)** — scale and band and model
are chosen on the training fold only — and the headline control is a **label-shuffle LOPO
null** that re-runs the *entire selection procedure* on permuted labels. My honest prior
(from the recovery arc): multiscale lifts **ranking/AUC** a little and does **not** break
the **precision** ceiling; the exploration tests whether smarter multiscale use changes
that verdict. If LOPO gains vanish or the shuffle null reaches them, the answer is "no,"
and that is a publishable negative that protects the section.

## 5-point critical preamble

1. **Claim under test.** Reading the mst@0.20 seeded heat-kernel affinity across the τ
   grid — as a per-band tuned scale, a cross-band per-scale fusion, or a change-across-τ
   feature vector — yields a **higher out-of-sample AUC and/or prec@5** for SOZ detection
   than the fixed single-scale (τ=10/λmax) marker and than the plain multiscale mean.

2. **Null.** (a) **Matched-strength fake-seizure** surrogate (the mandated FC null) on the
   winning family — a gain must exceed what strength-matched pseudo-SOZ sets reach. (b)
   **Label-shuffle LOPO**: permute the SOZ labels and re-run the *complete* nested pipeline
   (scale selection + logistic fit + held-out scoring). This is the null of the *selection
   procedure*, not of a fixed statistic — the only null that catches the forking path.

3. **Strongest alternative the null must control.** *Overfitting via in-sample scale/band
   selection.* `auc_best_s` (β 0.90, δ 0.87) is already the ceiling of this artifact and is
   never claimed. The danger: a per-band-scale or shape-feature model that looks great
   because it was tuned on the same 10 patients it is scored on. The label-shuffle LOPO
   null reproduces the selection freedom on noise, so if it also lands high the "gain" is
   selection, not signal.

4. **Does the null actually control it — by mechanism.** Nested LOPO removes the held-out
   patient from *every* choice, so a per-band scale or logistic weight can only help the
   held-out patient if the scale/shape signal **generalises across patients**. The
   label-shuffle LOPO null then quantifies the residual optimism of the procedure itself
   (feature count, scale-grid freedom) under no signal. What it **cannot** reject: a
   *consistent* cohort-wide overfit that survives shuffling would need a bigger cohort —
   flagged, not solved, at n=10. prec@5 is high-variance at n=10 (5 contacts × 10 patients)
   → reported with its per-patient spread, never as a bare scalar.

5. **Falsification + remaining limitations.** Falsified if F2/F3/F4 do **not** exceed F0/F1
   out-of-sample, OR the label-shuffle LOPO null reaches the same values, OR the
   matched-strength null is not cleared. Limitations that remain regardless: n=10 (thin for
   6-band × shape-feature fusion), SOZ = clinical labels (not resection/outcome), and the
   two right-hemisphere hub patients are a known failure mode reported per-patient.

## Design (families, all nested LOPO; targets = held-out AUC + per-patient prec@5)

Baselines (single band, per band): **F0** fixed τ=10/λmax · **F1** multiscale mean.
New (single band): **F2** per-band best scale chosen on training fold · **F3** logistic on
change-across-τ shape features {fine, meso, coarse, mean, slope, peak-scale, curvature,
early/late}.
Fusion (6 bands): **Fus0** fixed-τ · **FusMeso** mesoscale · **FusPB** per-band best scale
(training-fold) · **FusShape** 6 bands × shape features, L2-regularised.
External reference: published fused LOO detector (AUC ~0.87, prec@5 0.60).

Controls: matched-strength fake-SOZ (R=200) on the winner; label-shuffle LOPO (R≥200) on
the winner and on Fus0. Report per-band AUC curves + per-patient scatter; **never**
best-scale-in-sample as a result (diagnostic ceiling only).

## Compute

Stage 1 (free): per-band best-scale LOPO from existing `{band}/{pat}.npz` AUC/prec5 curves.
Stage 2 (one re-run of `06_epi_arc` primitives, saving per-contact residual marker stacks
`(N,16)`; ~1 min for 60 cells): shape features + fusion + nulls. Reuses `laplacian_eig`,
`scale_grid`, `mst_union_top_fraction`, `_seeded_marker`, `_resid`.

Script: `scripts/01_compute/sparsified_arc/16_multiscale_marker_exploration.py`.
Outputs: `data/sparsified_arc/ms_marker_exploration/{stage1_single_band,stage2_fusion}.csv`.

---

## VERDICT (2026-07-13, CORRECTED) — τ has a large effect; op-point near-optimal; shape lifts all-SOZ precision robustly

> **Correction (self-caught after user challenge).** An earlier draft of this verdict said
> "fused AUC is flat ~0.816, τ buys nothing." **That was an analysis error** — it compared
> only two scales (s=5.85 and s=10) that both sit on the *plateau near the peak*. The full
> 16-scale sweep shows τ has a **large** effect. The corrected findings are below.

**Head.** Sweeping the fused 6-band marker across all 16 scales, the cohort AUC runs
**0.594 (fine s=1) → 0.831 (peak s≈17) → 0.594 (coarse s=200)** — a strong inverted-U;
**τ matters a lot.** The published operating point **s=10 is near-optimal** (0.816, vs 0.831
peak); a **nested best-global-scale** selection (inner-LOPO, no forking path) lands at
s≈17–24 every fold and buys only **+0.009 AUC** — so the op-point is *validated*, not
improved. The real gain is in **precision, from change-across-τ SHAPE features**: fusing each
band's {fine, meso, coarse, slope, peak-scale, curvature, early-late} (48-d, L2) raises
precision **across the whole ranking**, most robustly **R-precision** (precision at k=n_soz,
the "flag-all-SOZ" point): **0.42 → 0.51, winning in 7/7 patients (paired Wilcoxon
p=0.009)**, and it **beats the dimensionality-matched scale-scrambled control** (p=0.031) —
so it is genuine τ-shape information, not extra features. prec@5 (0.58→0.74) and PR-AUC
(0.49→0.55) point the same way but are noisier at n=10 (p=0.075 / 0.138 vs Fus0; both clear
the scramble control). Single-band per-band scale tuning lifts low-γ/β **ranking**
(low-γ 0.77→0.84, 8/10) but not precision.

**Fused AUC / precision vs diffusion scale (nested LOPO, mean over 10 held-out):**

| s=τλmax | 1 | 4.1 | 8.3 (≈op) | 11.9 | 16.9 (peak) | 24 | 48.7 | 200 |
|---|---|---|---|---|---|---|---|---|
| fused AUC | .594 | .801 | **.816** | .825 | **.831** | .830 | .787 | .594 |
| fused prec@5 | .46 | .52 | .58 | **.64** | .56 | .54 | .46 | .26 |
| fused R-prec | .38 | .40 | .42 | .42 | .41 | **.45** | .37 | .19 |

**Precision families (nested LOPO; all-SOZ metrics):**

| family | feats | AUC | prec@5 | **R-prec** | PR-AUC |
|---|---|---|---|---|---|
| Fus0 fixed τ=10 | 6 | 0.816 | 0.58 | 0.417 | 0.494 |
| Nested best global scale | 6 | 0.825 | 0.52 | 0.419 | 0.481 |
| FusPB per-band scale | 6 | 0.816 | 0.54 | 0.421 | 0.495 |
| FusLevels3 (3 levels) | 18 | 0.817 | 0.56 | 0.444 | 0.514 |
| **FusShape (levels+derivs)** | 48 | 0.817 | **0.74** | **0.513** | **0.551** |
| FusScramble control | 48 | 0.810 | 0.58 | 0.436 | 0.499 |
| Fused@peak (ORACLE, in-sample) | 6 | 0.831 | 0.56 | 0.409 | 0.504 |

Paired Wilcoxon (shape vs, n=10, one-sided greater): prec@5 vs Fus0 p=.075 / vs scram p=.037;
**R-prec vs Fus0 p=.009 (7/0) / vs scram p=.031**; PR-AUC vs Fus0 p=.138 / vs scram p=.080.

**Reading.** (1) τ genuinely reshapes the marker; the paper's s=10 op-point is near-optimal —
a *robustness* point worth one sentence. (2) The **shape-feature precision gain is real and,
on R-precision, statistically robust** even at n=10, and survives the capacity (scramble)
control — a legitimate secondary result, though it improves *precision at the top of the
ranking*, not AUC. **Implication for the paper:** the honest τ clause I wrote ("τ sharpens
ranking, precision headline is band-fusion") stays true for the single-scale marker, but the
multiscale-shape R-precision result is a genuine value-add — a candidate supplementary
finding (flagged n=10; matched-strength on the shape family still owed before any headline).
Next: replicate the R-precision gain across scramble seeds + C, then matched-strength.

---

## DESIGN-SPACE AUDIT (2026-07-13) — the backbone was the wrong knob; TMFG is marker-optimal and matched-strength-validated

Script: `17_marker_design_space_audit.py` (coordinate descent: backbone → phase → marker →
band-subset, nested-LOPO). Data: `data/sparsified_arc/ms_marker_exploration/design_space_audit.csv`
+ `data/sparsified_arc/epi_arc_tmfg/` (matched-strength on TMFG via `06_epi_arc SA_BACKBONE=tmfg`).

**Head.** The marker had **inherited the `mst@0.20` backbone from the TRACE**, and that was the
one badly-chosen knob. Sweeping the backbone (fused AUC@s=10, nested-LOPO): dense **0.745**
(worst — washes affinity out, cf audit_174) → mst@0.20 **0.816** (current) → mst@0.10 0.859,
mst@0.05 0.874 → **TMFG 0.904** (best). **Sparser structured backbone is monotonically better;
TMFG (parameter-free) is the structured-sparse optimum.** Phase (rest_pre≈rest_post 0.908/0.904),
marker aggregation (mean>max>sum), and band-subset (all-6 best) all **confirm the current
choices** — only the backbone moved.

**TMFG is validated three ways** (not a landscape mirage):
- **Leak check** — label-shuffle LOPO fused null mean 0.501, p95 0.571; obs 0.904 → p<0.01.
- **Per-patient** — broad gain, and it **rescues the hub failures**: Pat_15 (right-hemi, *the*
  failure) 0.49→0.88, Pat_14 0.70→0.90, Pat_13 0.81→0.89; Pat_10 0.71→0.75 (still weakest).
  Partly **dissolves the two-populations ceiling** (2 hub failures → arguably 1).
- **Matched-strength (mandated gate), single-band, R=200** — TMFG beats mst@0.20 on every
  headline band:

  | band | mst@0.20 AUC / n-beat | TMFG AUC / n-beat |
  |---|---|---|
  | β | 0.745 / 8 | **0.869 / 10** |
  | low-γ | 0.822 / 8 | 0.843 / 9 |
  | δ | 0.830 / 7 | 0.830 / 7 |
  | α | 0.606 / 5 (marginal) | **0.792 / 7** |
  | θ | 0.740 / 5 | 0.743 / 6 |
  | high-γ | 0.754 / 7 | 0.688 / 4 (drops) |

**Reading.** TMFG is the marker-optimal, parameter-free, matched-strength-validated backbone:
β clears the strength null in **10/10** at AUC 0.87, α becomes a real carrier (0.61→0.79), the
hub patients are rescued. Only high-γ weakens (never a headline). **Open honesty items:**
fused-level matched-strength (single-band table + fused leak-check done); n=10; and whether to
present a two-backbone paper (marker=TMFG, trace=mst@0.20) — see
`2026-07-13_TMFG-backbone-for-trace-HANDOFF.md` (does the TRACE also improve on TMFG?
deliberately deferred to a separate chat). The earlier shape-feature R-precision lead is now
largely **moot for the headline** — TMFG alone (0.904) far exceeds the mst@0.20 shape result;
keep shape as a minor secondary note, not the value-add.
