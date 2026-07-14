---
name: 2026-07-13_TMFG-backbone-for-trace-HANDOFF
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
status: open-handoff
created: 2026-07-13
scope: FLAG-ONLY handoff for a SEPARATE chat — does the cross-phase α/β cophenetic TRACE
  improve on the TMFG backbone the way the epilepsy MARKER did? Not run here on purpose.
pointers:
  - .agents/reports/2026-07-13_multiscale-marker-exploration-scope.md   # where TMFG was found (marker side)
  - scripts/01_compute/sparsified_arc/13_matched_strength_mst020.py     # the trace MS gate to re-point
  - scripts/01_compute/sparsified_arc/12_sparsification_recovery.py     # backbone sweep for the trace
  - src/lrg_eegfc/utils/fc/backbone.py                                  # tmfg_backbone(W)
---

# Does the cross-phase trace improve on TMFG? — open handoff (for another chat)

## Head

On the **marker** side (§3 epilepsy) we found the pipeline had inherited the backbone from
the trace (`mst@0.20`) and that this was **suboptimal for SOZ detection**: switching to the
**parameter-free TMFG** backbone lifted the fused marker AUC **0.816 → 0.904**, rescued the
right-hemisphere hub patients (Pat_15 0.49 → 0.88), and **passed matched-strength** more
strongly than mst@0.20 (β **10/10** at AUC 0.87, α 0.61→0.79, low-γ 9/10; only high-γ drops).
The open question this raises — **deliberately not run here** — is whether the **cross-phase
α/β cophenetic TRACE** (the paper's flagship, §1) *also* improves on TMFG, or whether
mst@0.20 remains the right backbone for the trace. These are different read-outs of the same
operator and need not share an optimum.

## Why this is plausible (and why it might not transfer)

- **Plausible yes:** TMFG keeps the **strongest edges + a triangulated (cycle-rich) skeleton**.
  The recovery arc established that α needs a **strong-edge** backbone (percolation at the same
  density killed α; mst@0.20 kept it) and that the multiscale ladder needs **cycles** (near-tree
  MST/TMFG-poor fails; cycle-rich d≈0.10–0.20 works — audit_175). TMFG is strong-edge AND
  triangulated, so it *could* keep both α and β while sharpening scale structure.
- **Plausible no / caution:** the trace is a **cross-phase** quantity (ρ_sym of cophenetic
  distances between phases) validated end-to-end on mst@0.20 under matched-strength
  (β 16/16, α 12/16). TMFG has a **fixed** edge budget (~3(N−2)), not a density knob, so the
  τ-sweep behaves differently; the trace could sharpen, wash out, or shift bands. The
  marker's win does **not** imply the trace's — verify, don't assume.

## What to run (separate chat)

1. **Re-point the trace MS gate to TMFG.** Add `SA_BACKBONE=tmfg` support to
   `13_matched_strength_mst020.py` (same one-line backbone branch already added to
   `06_epi_arc.py`: `elif BACKBONE=="tmfg": B = tmfg_backbone(W)`), run the τ-resolved
   matched-strength ρ_sym gate. Compare per-scale clears-p<0.05 (of 16) and cohort p to the
   mst@0.20 table (β 16/16 p .001; α 12/16 p .007).
2. **Read per-scale, never scale-max** (OVERVIEW rule 4.1) — TMFG's fixed budget may move the
   δ/high-γ collapse artifacts; keep them un-claimed.
3. **If the trace improves or holds on TMFG:** consider a **single-backbone paper** (marker +
   trace both TMFG) — cleaner than a two-backbone story. **If it degrades:** keep mst@0.20 for
   the trace and present the marker's TMFG as a *read-out-specific* backbone (defensible: the
   directive already says the marker "reads the propagator differently").
4. **Downstream if TMFG adopted for the trace:** localization (β→OFC, enc→OFC, inf→cingulate),
   Grassmann, reinstatement all sit on the cophenetic distances and would need a re-derivation
   pass. Large. Scope before touching.

## Do-not (carry the honesty rules)

Matched-strength is the only null; per-scale not best-scale; state what weakened up front; do
not re-attribute the Grassmann (dense/whole-task) to a sparsified backbone. The marker's TMFG
result is validated (leak-check + matched-strength); the trace's is **unknown until run**.

## Paste-ready kickoff for the other chat

> Investigate whether the cross-phase α/β cophenetic trace improves on the **TMFG** backbone
> vs the current `mst@0.20`. Context: on the epilepsy-marker side, TMFG lifted the fused SOZ
> marker 0.816→0.904 and passed matched-strength more strongly than mst@0.20 (β 10/10 @0.87,
> α 0.61→0.79). Read `.agents/reports/2026-07-13_TMFG-backbone-for-trace-HANDOFF.md` first,
> then re-point `13_matched_strength_mst020.py` to `SA_BACKBONE=tmfg` (backbone branch already
> in `06_epi_arc.py`), run the τ-resolved matched-strength ρ_sym gate, and report per-scale
> clears vs the mst@0.20 β 16/16 / α 12/16 baseline. Per-scale not best-scale. Verify, don't
> assume the marker's win transfers.
