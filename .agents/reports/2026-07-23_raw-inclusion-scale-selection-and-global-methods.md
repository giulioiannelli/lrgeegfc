---
name: 2026-07-23_raw-inclusion-scale-selection-and-global-methods
kind: report
era: IMCOH_ABS_COHORT_N10 (mst@0.20 / tau-sweep)
status: current
created: 2026-07-23
scope: Fresh compute answering "does the raw-FC-vs-multiscale story make sense, and does it make sense of the difference with the global methods (Grassmann/geodesic)?" Runs the extended-scale cophenetic trace sweep (s=0.05..180) + dense-raw reference + scale-coherence under matched-strength (script 28), and assembles the global-method traces (ladder + probe table). VERDICT: the selection is REAL and the methods are consistent (not contradictory) — but the proposed selecting MECHANISM (scale-coherence) FAILED, and the hierarchy adds no new information (it is a selective filter, not a detector).
pointers:
  - scripts/01_compute/sparsified_arc/28_raw_inclusion_scale_selection.py    # CHECK 1 + 2
  - data/sparsified_arc/raw_inclusion_mst020/                                 # coph_scale/raw_ref/coherence/cohort_gate
  - scripts/01_compute/sparsified_arc/14_controls_ladder_mst020.py            # geodesic/resistance (CHECK 3)
  - scripts/01_compute/sparsified_arc/21_coph_beyond_raw.py                   # coph|raw ~ 0 (no new info)
  - .agents/reports/2026-07-23_old-vs-current-headlines-panoramic.md          # context
---

# Raw FC as the fine limit, the scale-selection, and the global methods

## Head

Three checks, fresh under matched-strength. **CHECK 1 works and is the good news:** dense raw FC registers the reorganization in four bands (δ, α, β, low_γ) — non-selective — while the multiscale cophenetic **selects**: it keeps β at every scale, α at s≥1, and **rejects low_γ at every single scale** (low_γ is raw-only), exposing δ as a coarse-scale artifact (dead through the fine/meso range, reappears only where the surrogate collapses). **CHECK 3 makes sense of the global methods:** the methods are not contradictory — they agree on β and dissociate the other bands by *which structure carries the trace* (β = every carrier; α = the nested hierarchy only; low_γ = global-modes/paths but NOT a hierarchy; δ = edge/coarse-artifact). **CHECK 2 is the honest bad news:** the selecting mechanism I proposed — scale-coherence of the task reorganization — **does not discriminate** (all six bands sit at 0.68–0.76; cross-band correlation with selection = +0.06, p=0.91; θ and low_γ are as "coherent" as β). So the selection is real and structured, but its *cause* is not scale-coherence — it is resolution-robustness (survival of the persistence under hierarchical coarse-graining), and that cause is not yet pinned. Two hard caveats stand: the hierarchy adds **no new information** beyond raw (script 21: coph|raw ≈ 0 — it is a selective filter, not a detector), and the "inclusion" is clean only for α/β (the UPGMA clustering already selects at s=0.05, so it is not a smooth τ-morph from raw's 4-band pattern).

## CHECK 1 — raw is the fine limit; the hierarchy selects (cohort matched-strength gate_p; * = clears)

Dense-raw reference (the s→0 single-edge trace) reproduces script 22 exactly (δ .024, θ .053, α .019, β .024, low_γ .032, high_γ .065 — non-selective, 4 bands clear). Then the cophenetic sweep, key columns:

| band | raw (dense) | s=0.05 | s≈0.4 | s≈1 | s≈5 (meso) | s≈15 | s≈50 | s=180 | pattern |
|---|---|---|---|---|---|---|---|---|---|
| **β** | *.024 | *.005 | *.003 | *.010 | *.002 | *.007 | *.005 | *.005 | **clears everywhere** (raw + all 28 scales) |
| **α** | *.019 | *.024 | .116 | *.014 | *.007 | *.010 | *.019 | .246 | raw + s≥1 (dips in the fine 0.07–0.42 zone; fades at the coarsest) |
| **low_γ** | *.032 | .539 | .539 | .423 | .246 | .246 | .278 | .278 | **raw ONLY — killed at every cophenetic scale** |
| **δ** | *.024 | .138 | .539 | .161 | .097 | *.014 | *.019 | .065 | raw + **coarse-only** (dead through fine/meso → artifact) |
| **θ** | .053 | .348 | .754 | .615 | .423 | .246 | .539 | .615 | never (clean null) |
| **high_γ** | .065 | *.032 | *.042 | .161 | .216 | .116 | *.019 | .097 | scattered (artifact) |

Reading: **the hierarchy is a selective filter.** low_γ is the clean demonstration — raw sees it, the cophenetic discards it at *all* 28 scales. δ's raw trace and its coarse-scale clearing are **disconnected** (dead in between), so the coarse δ is a surrogate-collapse artifact, not the same signal as raw δ. β survives every representation and every scale (the resolution-robust flagship); α survives from s≈1 up (a hierarchy-scale effect).

## CHECK 3 — the global methods, and why they differ (T_test trace, which bands clear)

| read-out | what it reads | bands that clear | placement |
|---|---|---|---|
| **dense raw FC** | pairwise edges (s→0 limit) | δ, α, β, low_γ | the fine limit — non-selective |
| **cophenetic (sweep)** | nested multiscale hierarchy | β (all s), α (s≥1) | selective; rejects low_γ/δ |
| **geodesic** | shortest paths on the backbone | δ, low_γ | path-length structure |
| **resistance** | all-paths (Laplacian pseudo-inverse) | none (dead) | global, multiscale-blind → averages everything away |
| **Grassmann** | leading global eigen-modes | β, low_γ | global-mode rotation |

The methods are **consistent, not contradictory** — each is sensitive to a different *carrier* of the reorganization, and the bands differ in which carrier holds their trace:
- **β** clears in raw, the cophenetic (all scales), AND Grassmann → its reorganization is carried redundantly (edges + nested hierarchy + global modes). That redundancy is *why β is the flagship*: every reasonable read sees it.
- **α** clears in raw and the cophenetic only (not Grassmann, not geodesic) → a **nested/hierarchical** effect, invisible to the leading global modes and to path length. This is the genuinely hierarchy-specific band.
- **low_γ** clears in raw, geodesic, AND Grassmann — but the cophenetic **kills it at every scale** → a global-mode/path/edge effect that does **not** form a resolution-robust hierarchy (consistent with §1's "low_γ carries the largest single tracers but is inconsistent"). This is the sharpest dissociation the multiscale read provides: it separates a global rotation from a held hierarchy.
- **δ** = edge/coarse-artifact; **θ** = null; **resistance** dead everywhere (a global multiscale-blind Laplacian read finds nothing — the strongest foil, as before).

So the answer to "make sense of the difference with the global methods": **β is universal, α is hierarchy-only, low_γ is everything-but-hierarchy, δ is edge/artifact.** The multiscale cophenetic is the only read-out that is both selective (α/β) and able to reject a global rotation (low_γ) that Grassmann/geodesic mistake for a trace.

## CHECK 2 — the proposed mechanism FAILED (report this, do not bury it)

I proposed that α/β survive because their task reorganization is *scale-coherent* while the discarded bands' is scale-local. **The data reject this metric as the discriminator.** Cross-scale coherence C of the per-pair task reorganization (mean cross-scale Spearman of D_task − D_pre), cohort medians:

| band | C (s≥1) | scales cleared |
|---|---|---|
| β | +0.691 | 28/28 |
| θ | +0.719 | 0/28 |
| low_γ | +0.684 | 0/28 |
| α | +0.684 | 17/28 |
| δ | +0.689 | 7/28 |
| high_γ | +0.755 | 7/28 |

Cross-band Spearman(C, scales_cleared) = **+0.06, p=0.91** — no relationship. θ (never clears) is *more* coherent than β (always clears). The metric is dominated by the trivial scale-smoothness of cophenetic distances, common to every band, so it neither confirms nor is the right probe. **Honest status: the selection is real and structured (CHECK 1) but its cause is not scale-coherence.** The likely cause, pointed to by CHECK 1 + §1's "largest single tracers" observation, is **edge-locality of the persistence**: low_γ's trace is carried by a few strong edges that the UPGMA clustering absorbs, while β's is distributed/nested and survives coarse-graining. That is a targeted follow-up probe (concentration of the per-pair trace, or fraction surviving projection onto the hierarchy), **not yet run**.

## What survives for the paper, and what does not

- **DEFENSIBLE (demonstrated):** raw FC is the non-selective fine limit; the multiscale cophenetic is a **selective filter** (α/β only) that additionally **dissociates bands by structural carrier** — uniquely separating a held hierarchy (α/β) from a global rotation (low_γ) that raw/Grassmann/geodesic cannot tell apart. β's cross-representation robustness is the flagship result.
- **NOT DEFENSIBLE (drop it):** "only the hierarchy sees the trace" (raw sees α/β; script 21: no new info). "Scale-coherence is why α/β are selected" (CHECK 2 failed). A smooth τ-morph from raw to α/β (the clustering selects already at s=0.05 — the filter is the *hierarchical clustering*, active at all scales, not a gradual τ effect).
- **OPEN (next probe):** the *cause* of the selection — the edge-locality hypothesis (low_γ = few strong edges washed by clustering; β = distributed/nested). Until this is shown, the selection is real but mechanistically unexplained; state it as "the hierarchy retains the reorganizations that are robust to coarse-graining and discards the edge-local ones," with the cause flagged as under test.

## Reproduce

```
PY=/home/giulio/Documents/miniconda3/envs/lapbrain/bin/python
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 SA_WORKERS=12 $PY \
  scripts/01_compute/sparsified_arc/28_raw_inclusion_scale_selection.py \
  --bands delta,theta,alpha,beta,low_gamma,high_gamma --R 200      # ~4 min
```
