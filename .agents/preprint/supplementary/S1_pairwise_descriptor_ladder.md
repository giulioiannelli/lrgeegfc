---
name: supp-pairwise-descriptor-ladder
era: IMCOH_ABS_COHORT_N10
status: current
kind: supplementary
scope: Classical low-level pairwise graph descriptors run through the identical ρ_sym + matched-strength pipeline as the cophenetic trace. Verdict — local scalars are blind (strength degenerate under MS by construction); band-selectivity is unique to cophenetic; the apparent β-inference exceptions (raw/eigcent/closeness under MS) are DRIFT-CONTAMINATED and not band-selective (audit_167/170), so not a genuine second witness.
updated: 2026-07-10
---

# S1 — Pairwise-descriptor ladder: do the usual low-level scalars reproduce the trace?

**Head.** A ladder of classical pairwise-network descriptors — node strength,
weighted clustering, eigenvector centrality, PageRank, closeness — run through the
**identical** cross-phase ρ_sym estimator and matched-strength null as the LRG
cophenetic trace, ordered by structural order (raw edges → node scalars →
cophenetic). Verdict: the *usual quantities* recover **none** of the trace (node
strength is degenerate under matched-strength **by construction**; weighted
clustering and PageRank are ≈strength-slaved and blind), and **band-selectivity —
recovering α while silent on the null band θ — is unique to the cophenetic
comparison**. The apparent exceptions — raw edges, eigenvector centrality and
closeness appearing to detect the β inference-specific trace under matched-strength
— are **drift-contaminated and not band-selective** (audit_167/170): matched-strength
does not control drift, and a within-rest drift arc reproduces their β value. They
are not a genuine pairwise witness.

## Design
- **Estimator held fixed** (ρ_sym cross-phase arc); only the per-phase
  representation φ(W) changes across rungs. raw edges / cophenetic are pair-level;
  the classical scalars are node-level. ρ_sym is a Spearman → representation-agnostic.
- **Null held fixed:** matched-strength, exact flagship ensemble (seed 20260511),
  so the cophenetic rung reproduces audit_150/152 bit-for-bit (anti-hallucination
  gate: passed).
- **Scored** on the cophenetic signatures: recovers α · silent on θ · survives MS ·
  β inference-specificity.

## Result (gate p; whole-task T_test | β inference T_infspec·e)

| rung | strength-slaved | α | β | θ | β inference |
|---|---|---|---|---|---|
| **cophenetic** | — | **.032** | **.032** | .784 (silent) | **.010** (β-only) |
| raw FC edges | — | .053 | .024 | .053 (not silent) | .010 → **drift-contam.** |
| node strength | 1.00 | 1.000 | 1.000 | 1.000 | 1.000 — **blind** |
| weighted clustering | 0.99 | .053 | .216 | .461 | .161 — blind |
| PageRank | 1.00 | .577 | .976 | .188 | .920 — blind |
| eigenvector centrality | 0.99 | .188 | .019 | .839 | .002 → **drift-contam.** |
| closeness | 0.90 | .216 | .097 | .216 | .024 → **drift-contam.** |

*"drift-contam." = clears matched-strength but its β-inference value is reproduced by
a within-rest drift arc (audit_167/170; drift floor 0.14–0.20 ≥ the real value; no
low-level rung's real value clears its own drift floor). Not band-selective either.*

## Reading
1. **The named "usual quantities" are blind.** Node strength is *degenerate under
   matched-strength by construction* — MS preserves it exactly, so the
   strength-trace is reproduced and can never clear the null (p = 1.000 all bands;
   the cleanest possible "strength is not causing our results"). Weighted
   clustering (0.99 strength-slaved) and PageRank (1.00) clear nothing.
2. **Band-selectivity is unique to cophenetic.** Only the cophenetic rung recovers
   α (.032) while staying silent on the null band θ (.784). Raw FC is broad
   (≈.05 every band, θ not silent); no pairwise scalar reproduces the α-recovery +
   θ-silence signature.
3. **The apparent exceptions are drift-contaminated, not a second witness.** Raw
   edges, eigenvector centrality and closeness *appear* to reproduce the β
   inference-specific trace under matched-strength — but MS does not control drift.
   The within-rest drift ladder (audit_167 raw+coph; audit_170 full menu) windows a
   single pre-task rest recording into 5 ordered segments — no task, no consolidation
   possible — and finds the conditional inference functional is **drift-biased
   positive for every low-level representation**: within-rest drift alone gives β
   `T_infspec·e` of raw 0.14 / strength 0.19 / clustering 0.17 / eigcent 0.17 /
   PageRank 0.17 / closeness 0.20 (drift-alone p=0.001), **as large as or larger than
   the real values**, and **no** low-level rung's real value significantly clears its
   own drift floor (raw sits *below* it). So their MS "detection" is drift-dominated.
   Cophenetic carries the **smallest** drift contamination (drift 0.073, ~half) and
   out-filters raw on the bias-robust differential (coph − raw of (real−drift): β
   p=0.032, β-only). **NB:** the windowed drift arc is an *invalid absolute gate* for a
   conditional functional (it is biased positive even on the no-signal arc, including a
   residual for cophenetic, `T_infspec·e` drift p=0.065) — so the inference-specific
   claim rests on **matched-strength (≈0 on no-signal, p=0.0098) + duration**, and the
   clean cross-descriptor discriminator is **band-selectivity** (only cophenetic
   recovers α while silent on θ), not the drift gate.

## Verdict
- The disaster ("node strength / clustering already reproduce our trace") is
  **averted** — the local scalars are blind.
- The cophenetic contribution is **band-selectivity** (the clean cross-descriptor
  discriminator) **+ least drift contamination + fine-grained** (Grassmann-null,
  audit_165), not mere detectability. The pairwise/scalar β-inference "convergence"
  is drift-contaminated and not band-selective — **not a genuine second witness**.
- **Strength-slaving diagnostic** (median |Spearman(φ, strength)| on dense FC):
  clustering 0.99, PageRank 1.00, eigcent 0.99, closeness 0.90 — near-total; the
  classical scalars are not independent of the strength the MS null already fixes.

Build: `scripts/01_compute/audit/audit_166_pairwise_descriptor_ladder.py` (MS ladder),
`audit_167_drift_null_ladder.py` (drift, raw+coph) + `audit_170_drift_ladder_descriptors.py`
(drift, full menu); library `src/lrg_eegfc/utils/metrics/graph_descriptors.py`;
data `data/audit/{pairwise_descriptor_ladder,drift_null_ladder,drift_ladder_descriptors}/`.
Full method + 5-point preamble: scope report
`.agents/guides/task-persistence-investigation/2026-07-09_pairwise-descriptor-ladder.md`.
