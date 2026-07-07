---
name: 2026-07-06_raw-vs-cophenet-resolving-power
type: report
era: IMCOH_ABS_COHORT_N10
status: current
created: 2026-07-06
updated: 2026-07-06
pointers:
  - scripts/01_compute/audit/audit_156_raw_vs_cophenet_resolving_power.py
  - data/audit/raw_vs_multiscale/band_taxonomy_raw_vs_multiscale.csv
  - .agents/guides/02_methods/lrg-framework-guide.md
---

# Does the LRG cophenet resolve structure raw FC cannot? Yes — a double dissociation

**Head.** The worry "raw FC already shows the trace, so the multiscale machinery
is redundant" is answered by a constructive **double dissociation**: raw pairwise
FC and the LRG cophenet are provably **non-equivalent** representations. The
cophenet resolves persistence of the community **hierarchy** that raw FC loses;
raw FC responds to hierarchy-neutral edge **drift** that the cophenet filters.
Raw FC is not ill-defined — it is a *local* representation that is **blind to
community-level persistence** and **conflates it with drift**. Only the
multiscale representation separates the two.

## The test (audit_156, SBM ground truth, split-baseline trace matching ρ_sym)

`trace(X) = Spearman(X_task − X_preA, X_post − X_preB)`, `X ∈ {raw W, cophenet(W)}`.
We sweep a nuisance parameter and require the dissociation to be a **regime**, not
a lucky point (critical preamble in the script header).

**Regime H — hierarchy persists, edges noisy.** Task reorganizes the block
structure, rest_post keeps it, within-block edges are redrawn independently each
phase. As edge noise σ rises, the raw trace **collapses** while the cophenet
trace **holds**:

| σ (edge noise) | raw trace | cophenet trace | gap |
|---|---|---|---|
| 0.15 | +0.54 | +0.64 | +0.10 |
| 0.30 | +0.34 | +0.62 | **+0.29** |
| 0.40 | +0.23 | +0.55 | **+0.32** |
| 0.60 | +0.12 | +0.28 | +0.16 |

→ the community hierarchy persists; **only the multiscale method resolves it**
once edges are noisy. This is the β-like regime.

**Regime E — hierarchy-neutral drift.** Blocks static; a persistent within-block
zero-mean edge pattern. Raw sees it; the cophenet is blind until the drift grows
strong enough to become structural:

| within-block drift amp | raw trace | cophenet trace | gap |
|---|---|---|---|
| 0.15 | +0.10 | +0.07 | +0.03 |
| 0.22 | **+0.17** | +0.05 | **+0.12** |
| 0.40 | +0.30 | +0.28 | +0.02 |

→ raw FC flags non-hierarchical drift that the cophenet **correctly filters**.
This is the θ-like regime.

Figure: `data/reports/rho_sym_band_map/fig_raw_vs_cophenet_resolving_power.pdf`.

## The real data lives in both regimes (data/audit/raw_vs_multiscale)

- **β = Regime H.** Raw trace +0.258 (6/10) AND cophenet +0.221, cohort agreement
  **sharpened** 6→7/10. β's task reorganization is hierarchically coherent and
  persists at the community level.
- **θ = Regime E.** Raw trace +0.117 (7/10) but cophenet **−0.040 (2/10)**. θ's
  raw persistence is hierarchy-neutral drift; the cophenet filters it.
- **low-γ vs δ** (equal raw ≈0.11–0.13): cophenet keeps low-γ (0.083) and
  annihilates δ (0.008) — the filter discriminates on multiscale coherence, not
  raw magnitude.

So the discrimination raw FC **cannot** make — which bands carry hierarchical
persistence vs mere drift — is exactly the multiscale contribution, and it is
what carries the β→OFC system-scale localization (a hierarchy property, not an
edge property).

## Honest scope (what this proves and does not)

- **Proves:** cophenet ≠ raw (double dissociation); the cophenet resolves
  community-hierarchy persistence non-resolvable by raw FC (Regime H, robust across
  σ); raw FC is blind to it and conflates it with drift (Regime E).
- **Does NOT claim** the real β trace is raw-invisible — β has a raw trace too.
  The multiscale-specific content is (a) the **discrimination** (β vs θ, low-γ vs
  δ), (b) the **community-level persistence** phenomenon, (c) the **localization**.
- **Limitation:** SBM blocks are idealized; real FC is not exactly block-
  structured. The synthetic shows the phenomenon is *possible and non-resolvable*;
  the raw_vs_multiscale placement shows it is *realized* in the cohort.

## Recommended confirmatory follow-up (not yet run)

Empirical partial correlation on real β: per patient, test whether the cophenetic
cross-phase persistence survives regressing out the raw-edge cross-phase
persistence (pair-level partial ρ). A surviving hierarchy-specific component would
be the direct on-data version of Regime H. Needs W per phase + `canonical_cophenet`
+ the split baseline; ~1 audit.
