---
name: w0a_00_preregistration
kind: pre-registration
era: PAPER_FINALIZATION (Wave 0, lane W0-A)
status: locked
created: 2026-08-25
scope: The decision rule for the two substrate choices (band transform, backbone), written and frozen BEFORE any number was computed. Selection is on hypothesis-independent properties only. Downstream trace / AUC numbers are recorded but are explicitly NOT admissible as selection criteria.
pointers:
  - .agents/plans/active/2026-08-25_paper-finalization-master-plan.md
  - .agents/guides/02_methods/sparsification-choice.md
  - scripts/01_compute/paper_final/w0a_01_transform_headtohead.py
  - scripts/01_compute/paper_final/w0a_02_sparsification_stability.py
---

# W0-A pre-registration — the substrate decision rule

## Head

Two knobs (`imcoh_abs` vs `imcoh_sq`; backbone family × density `f`) are currently set by argument and by outcome-selection respectively. This document freezes, before any number is seen, the rule that will set them. The rule reads only properties that are independent of the task-trace hypothesis: estimator reliability, distance from a lag-destroying null floor, and invariance of the verdict under the knob. Trace p-values and SOZ AUCs are recorded for the record and are **not** inputs to the rule.

## 5-point critical preamble

**(1) The claim.** That the substrate of the paper — the FC transform `⟨|Im C|⟩_f` and the sparsified backbone — can be fixed by criteria that do not reference the cross-phase trace, and that the resulting substrate sits inside a region of the density knob over which the cross-phase verdict does not change.

**(2) The null.** H0 for the substrate contract: the cross-phase verdict is a function of the knob — i.e. there is no region of `f` over which the per-band, per-scale verdict is invariant, and any stated verdict is an artifact of the chosen `f`. For the transform arm the null is: the two transforms are rank-equivalent on the edge weights, so nothing downstream can distinguish them and the choice is vacuous.

**(3) The strongest plausible alternative the null should control for.** That a "plateau" in the verdict is manufactured by the statistic rather than by the graph — specifically, that adjacent `f` values give correlated verdicts simply because they share almost all their edges and almost all their surrogate draws, so the surface is smooth by construction and a flat region carries no information. A second alternative: that the transform comparison is decided by a monotone rescaling artifact (Spearman is invariant to monotone maps, Pearson is not), so "reliability" differences are units, not information.

**(4) Whether the null controls for it, by mechanism, and what it cannot reject.** The knob sweep uses *independent* surrogate draws per cell but the **same** observed graph family, so correlated verdicts across adjacent `f` are expected and are NOT by themselves evidence of a plateau. The control is therefore structural, not statistical: the plateau claim is only made where the verdict is invariant **across a change in the density-selection mechanism** (mst-union vs plain global threshold vs disparity), not merely across neighbouring `f` of one mechanism, and it is accompanied by the structural covariates (density, mean degree, clustering, algebraic connectivity, spectral gap, retained weight fraction) that show the graph itself changed materially over the plateau. It **cannot** reject the possibility that all mechanisms in the panel share a common bias (all are weight-ranked filters of the same dense matrix); a genuinely different construction — e.g. an inference-based network model — is outside this lane. For the transform arm, the monotone-rescaling alternative is controlled by reporting **Spearman** reliability as the primary (monotone-invariant) figure, with Pearson secondary; if the two transforms were exactly rank-equivalent, Spearman reliability would be numerically identical by construction, so any Spearman difference is real information about the band-averaging order and not about units.

**(5) What would falsify the claim, and remaining limitations.** The substrate contract is falsified if (a) no contiguous region of `f` of width ≥ one octave shows an invariant per-band verdict, or (b) the verdicts of mst-union and plain global threshold disagree at matched density, or (c) the transform arm shows the two transforms give materially different reliability *and* the higher-reliability transform is the one that kills the trace — in which case the honest report is that the substrate cannot be fixed hypothesis-independently and the paper must state knob-dependence explicitly. Remaining limitations: the null used inside the sweep is matched-strength, which is injected at the FC-matrix stage and therefore cannot test the coherency estimator, the band split, or session nonstationarity (lane W0-B); and the split-half reliability is a within-session test-retest, which cannot separate estimator noise from genuine within-session nonstationarity.

## The rule (frozen 2026-08-25, before any number)

### R1 — transform

Choose between `imcoh_abs = ⟨|Im C|⟩_f` and `imcoh_sq = ⟨(Im C)²⟩_f` on, in strict lexicographic order:

1. **Rank-equivalence gate.** If the two are rank-equivalent on the upper-triangle edge vector to within numerical tolerance (Spearman ≥ 0.9999 and top-`f` edge-set Jaccard = 1 at `f ∈ {0.05, 0.10, 0.20}`), then any rank-based backbone is invariant to the choice, the choice is downstream-vacuous, and we keep the incumbent `imcoh_abs` and **report the invariance as the justification**. Rule stops here.
2. **Split-half reliability.** Otherwise, prefer the transform with the higher cohort-median **Spearman** reliability between the two rest_pre halves' edge vectors (test-retest of the estimator itself), evaluated per band and aggregated over the n = 10 cohort. A difference is decisive only if the transform wins in ≥ 4 of the 6 bands **and** the cohort-median difference has a sign-consistent Wilcoxon at p < 0.05 across the 60 (patient, band) cells.
3. **Separation from the lag-destroyed floor.** Tie-break: prefer the transform whose observed edge weights sit further above their own circular-shift floor, measured as the cohort-median ratio `median(W_obs) / median(W_floor)` on the upper triangle, and as the rank correlation between observed and floor (lower = more of the observed structure is genuinely lag-carried).
4. **Dynamic range / heterogeneity.** Final tie-break: prefer the transform with the *lower* weight heterogeneity (Gini / CV of the edge weights), because Villegas-regime degeneracy is driven by weight heterogeneity on a dense graph and a less heterogeneous weight field leaves the propagator further from the single-peak collapse.

Trace gates and marker AUCs under each transform are computed and reported but are **not** admissible at any step of R1.

### R2 — backbone

1. **Parameter-free first.** A parameter-free filter (TMFG, PMFG, percolation-at-the-connectivity-bottleneck, disparity at a principled α) is adopted **iff** its (density, verdict) point lies inside the invariance plateau found in the sweep — i.e. its per-band, per-scale verdict is the same as the plateau verdict, and its density falls inside the plateau's density interval. Multiple qualifying filters: prefer the one whose density is closest to the plateau centre in log-density; break remaining ties by literature standing (PMFG > TMFG > percolation > disparity).
2. **Plateau definition, frozen.** A plateau for band `b` is the maximal contiguous set of `f` on the swept grid over which the per-scale cohort gate verdict (cleared / not cleared at α = 0.05, per scale, no best-scale collapse) is constant, requiring width ≥ one octave in `f` (i.e. `f_max / f_min ≥ 2`) and requiring the verdict to be reproduced by **at least two different density mechanisms** at matched density.
3. **Fallback: knob-integrated readout.** If no parameter-free filter qualifies, the reported statistic is the **median over the plateau** of the per-scale cohort statistic, with the plateau width stated, and no single `f` is ever reported as the pipeline setting. Stability of the knob-integrated readout is demonstrated by leave-one-`f`-out over the plateau.
4. **Failure mode, to be reported honestly.** If no plateau of width ≥ one octave exists for a band, that band's verdict is reported as **knob-dependent**, with the surface shown. It is not rescued by picking an `f`.

### R3 — the two-backbone question

A single substrate is adopted iff the trace plateau and the marker's parameter-free choice coincide (same filter, or the marker's filter lies inside the trace plateau). Otherwise the two-backbone posture stands and must be justified in Methods by the density argument, not by performance.

## Addendum R1-bis — added 2026-08-25 AFTER A1 ran, and BEFORE any A2 trace number was computed

This addendum is dated and firewalled: A1's four steps had been evaluated, but no cross-phase trace statistic under either transform had been computed when it was written. It does not relax R1 — R1's selection stands unless invariance makes it moot.

R1 selects `imcoh_sq`, at step 2, by a cohort-median Spearman split-half reliability margin of 0.0078 (≈ 1.2 % relative), winning 5 of 6 bands, Wilcoxon p = 5.0e-6 over the 60 (patient, band) cells. That margin is systematic but small, while the consequences are not: the two transforms disagree about 10–22 % of the strongest-edge backbone (top-5 % edge-set Jaccard 0.78–0.91) and `imcoh_sq` nearly doubles the edge-weight Gini (β: 0.27 → 0.51; participation ratio 0.76 → 0.31). A 1 % reliability margin is not a sufficient basis on which to regenerate every artifact in the project.

Therefore: the A2 stability sweep is run under **both** transforms. If the per-band, per-scale verdict is invariant across the transform inside the plateau, the contract records the transform as **immaterial within the plateau** and retains `imcoh_abs`, citing the demonstrated invariance — not the reliability comparison — as the justification. If the verdict is **not** transform-invariant, `imcoh_sq` is adopted exactly as R1 dictates, and every downstream artifact is regenerated on it.

## What is explicitly out of scope for this lane

The timeseries-level null ladder (W0-B), the cohort gate specification (W0-C), the SOZ side beyond a cheap check (L3). Any claim here is conditional on matched-strength being the null, which it should not remain.
