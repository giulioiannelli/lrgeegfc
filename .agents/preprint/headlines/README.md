---
name: headlines-readme
era: IMCOH_ABS_COHORT_N10
status: current
kind: headline
scope: purpose + routing for the headlines/ folder — one self-contained report per headline result of the multiscale-Laplacian programme
---

# `headlines/` — headline results of the multiscale-Laplacian programme

**Head.** One file per **research question** that becomes a part of the
paper on the *multiscale Laplacian density-matrix approach to reveal
interesting connectivity patterns that inform on higher-level brain
behaviour*. Each file is a **self-contained report**: the core observation,
the validated findings with quantities, the verification, the honest scope,
the discussion/positioning, the literature to cite (with **real DOIs —
never invented**), and the highest-value next steps. These are the
"juice-locked" headlines; the per-band briefs (`../bands/`) and locked
ledgers (`../locked/`) remain the operational sources of truth.

The unifying object across all headlines is the Laplacian density operator
ρ̂(τ) = e^{−τL̂}/Z (Villegas LRG): the same diffusion propagator whose
multiscale cophenetic structure carries the **β-band consolidation trace**
also, read as a node/pair diffusion affinity, exposes the **δ-band
epileptogenic diffusion community**. One description, multiple read-outs.

## Files

- [`epi-marker-analysis.md`](epi-marker-analysis.md)
  — the **complete** epileptic-marker investigation on ρ(τ)=e^{−τL}/Z (all bands,
  every audit): node-intrinsic markers fail (= hubness + electrode depth); the
  relational δ diffusion community is the surviving finding (distant SOZ, AUC 0.72,
  8/10, strength-orthogonal, phase-stable); includes the retraction ledger. Seed-based
  triage, not outcome-validated. Headline #1.
- (planned) the β-band cross-phase trace headline — to be migrated/condensed
  here from `../bands/01_beta.md` + `../locked/VERDICT_LEDGER.md` when the
  user asks.

## Rule for files in this folder

A headline file is **complete on its own** and **brutally honest**:
1. Core observation (one paragraph).
2. Findings + quantities (validated; cite the CSV/audit row).
3. Verification (the nulls and controls that make it not-an-artifact).
4. Honest scope + limitations — stated **first-class**, not buried.
5. Discussion / positioning within the multiscale-Laplacian programme.
6. Literature: only papers actually read; **copy-paste the DOI, never invent**.
7. Open questions / next steps.

Provenance (scripts, CSVs, figures, prior reports) goes at the foot of each file.
