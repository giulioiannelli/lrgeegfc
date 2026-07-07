---
name: task-persistence-investigation
type: guide
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-04-25
updated: 2026-07-04
pointers:
  - .agents/guides/task-persistence-investigation/2026-04-25_task-trace-canonical.md
  - .agents/preprint/headlines/README.md
  - .agents/guides/task-persistence-investigation/2026-06-25_per-node-trace-decomposition.md
---

# Task-persistence investigation

**Canonical home for every multiscale measure built to investigate task-induced
reorganization in the post-task resting state. New measures land here as
mathematically rigorous scope reports first; code follows. Mathematical
rigor is required so a careful reader can spot definition errors, edge cases,
and identifiability problems before any compute is launched.**

## Purpose

The standing scientific question (mirrored at the top of `CLAUDE.md`):

> Does `task_test` leave a band-specific, multiscale structural trace in
> `rest_post` LRG dendrograms that is cohort-wide (≥ 8/10 at n=10 under
> `imcoh_abs`)?
>
> Cohort returned to n=10 on 2026-04-25 with Pat_14 vendor-replacement.
> Canonical reformalization at the leaf-set level (regimes P/T/R/RA):
> [`2026-04-25_task-trace-canonical.md`](2026-04-25_task-trace-canonical.md).

This folder collects the descriptive multiscale measures we develop to surface
that exact question. Explicitly **not** scalar-gate hypothesis tests at
n=9 + FDR q < 0.05 — that loop was closed in the April scalar session (now in
`archive/scalar-vi-era/`). The measures that survived and now back the paper's
headlines are indexed under **Live measures** below; every dead-end is preserved
under **Archived measures (by era)**.

## File layout

```
task-persistence-investigation/
├── README.md                            # this file (index)
├── YYYY-MM-DD_<measure-name>.md         # LIVE scope reports (10 — see below)
└── archive/                             # retired measures, by research era
    ├── scalar-vi-era/                   #   April VI(k)/CBR/MRL/substrate/E1 dead-ends
    ├── kc-era/                          #   May KC / CTM / RF tree-distance family
    ├── epi-negatives/                   #   negative per-node epi-marker attempts
    ├── replay-n4/                       #   the pulled replay headline (former N4)
    └── opaque-null-era/                 #   rejected CRC / coordinated matrix nulls
```

Naming: `YYYY-MM-DD_<short-kebab-name>.md`. Date = day the scope was first
drafted; do not change on revisions, update the `updated` frontmatter field
instead. The date prefix is retained inside archive buckets so re-analyses of the
same idea remain distinguishable by date.

## Required scope-report structure

Every report under this folder MUST contain, in order:

1. **Frontmatter** — `name, type=scope, era, status ∈ {draft, current, superseded}, created, updated, pointers`.
2. **Renormalization head** — 1–2 sentences. The juice. A reader skimming the
   first paragraph must come away knowing what the measure is and what
   question it answers.
3. **Notation** — symbols, indices, sets, and the mathematical objects
   the measure operates on. Define once; never overload.
4. **Definitions** — each predicate / function / aggregator stated as a
   formula with its domain and range. Cite literature when borrowing
   (Kendall–Colijn 2016, Robinson–Foulds 1981, Maris–Oostenveld 2007, etc.).
5. **Properties** — range, monotonicity, invariance, identifiability,
   complexity. State *what the measure cannot detect* (negative properties
   are as important as positive ones).
6. **Caveats & failure modes** — probe bias at coarse `h_rel`, edge effects
   at the root, ties in heights, parameter sensitivity, finite-N artifacts,
   anatomy-dominated baselines. Each caveat with a mitigation.
7. **Pseudocode** — language-agnostic, line-oriented, every loop bound
   explicit. Not Python — algorithmic.
8. **Visualization spec** — figure axes, colormaps, annotation rules,
   what a reader is supposed to *see*. Include reading rules
   ("vertical bands of color = nested persistence"; "horizontal contrast = band heterogeneity").
9. **Connection to prior tools** — a paragraph or table mapping this
   measure to existing operationalizations.
   What does it complement; what does it subsume; what does it not replace.
10. **Implementation plan** — script path, output paths, library entry points
    (must reuse `lrg_eegfc.utils.metrics.tree*`; no new helpers without a
    promotion-to-library justification per `coding-rules.md`).
11. **Open questions** — parameter defaults, scope decisions deferred.
    Marked explicitly so a future reader knows what's still load-bearing.

## Status semantics

- `draft` — written, not yet implemented or sanity-checked.
- `current` — implemented and feeding into a live analysis.
- `superseded` — replaced by a later scope; keeps the file (don't delete);
  moved to the relevant `archive/<era>/` bucket with a top-of-body pointer.

## Live measures (back the current paper — headlines N1/N2/N3)

| Date | Measure | Backs | Tagline |
|------|---------|-------|---------|
| 2026-04-25 | [Task-trace canonical (P/T/R/RA)](2026-04-25_task-trace-canonical.md) | N1 foundation | Definitional umbrella — the four regimes (Persistence / Trace / Reset / Rearrange) at the leaf-set level; shared vocabulary for every downstream measure. |
| 2026-04-26 | [Continuous-trace matrix](2026-04-26_continuous-trace-matrix.md) | N1 measurement | Methodological ancestor of `ρ_split`: per-pair `Δ_task`/`Δ_rest` on the distance matrix, drift-floor null, same-probe control. (Operates on raw `D(τ_max)`; the headline uses the cophenetic `D_coph` refinement.) |
| 2026-05-30 | [Asymmetric pair-trace regression `s_TR`](2026-05-30_asymmetric-pair-trace-regression.md) | N1 method-lock | Comparative audit that keeps Spearman `ρ_split^coph` over slope / `R²` / Pearson alternatives under the shared matched-strength null. |
| 2026-06-05 | [Sampling-conditioned localization atlas](2026-06-05_localization-atlas.md) | N1.2 (β→OFC) | Sampling-conditioned positive-rate localization — the anatomical concentration claim (β→OFC system). |
| 2026-06-12 | [Cophenetic consolidation arc](2026-06-12_cophenetic-consolidation-arc.md) | **N2 backbone** | Encoding vs inference: `e = D_taskL − D_rsPre`, `f = D_taskT − D_taskL`; the flagship arc that separates inferred structure from encoded pairs. |
| 2026-06-12 | [Cross-phase cophenetic taxonomy](2026-06-12_cross-phase-cophenetic-taxonomy.md) | N1.6 bridge | Cut-free anchor / trace / reset / reorganize decomposition of `ρ^coph` (φ1 persistence ⊥ φ2 excursion). |
| 2026-06-18 | [Inference-mark localization](2026-06-18_inference-mark-localization.md) | N2.6 | Localizes the arc components (encoding→OFC; inference→cingulate, duration-downgraded to a directional hint). |
| 2026-06-22 | [τ-sensitivity of the cophenetic trace](2026-06-22_tau-sensitivity-cophenetic-trace.md) | N1.5 (τ-robust) | τ-sweep vindicating `τ=1/λ_max`: the genuine trace is fine-scale and τ-robust; coarse-τ "gains" are a collapse artifact. |
| 2026-06-22 | [τ-sweep of the consolidation arc](2026-06-22_tau-sweep-consolidation-arc.md) | N2.2 (mesoscale) | Scale-sweep of the inference-specific component — mildly mesoscale-favouring, robust not exclusive. |
| 2026-06-25 | [Per-node trace decomposition](2026-06-25_per-node-trace-decomposition.md) | N1 per-node | Bottom-up explanator of per-patient heterogeneity: `ρ_split` → per-node carrier / neutral / anti-trace scores, matched-strength calibrated. |

## Archived measures (by era)

Preserved for the record; retired, do not cite as live. Each file keeps its
original date prefix so the archaeology is legible.

- **`archive/scalar-vi-era/` (23)** — the April VI(k)/scalar-session lineage.
  The CBR/MRL single-module-matching family (`cbr-*`, `cohesion-cbr`,
  `containment-cbr`, `consensus-subtree`, `task-anchored-cbr`,
  `module-retention-landscape`, `mrl-vs-cbr-reconciliation`,
  `cophenetic-neighbourhood`, `trace-modules`, `multiscale-partition-coherence`);
  the substrate + rebuild-ladder diagnostics (`functional-tree-distance`,
  `residual-subspace-trace`, `edge-vs-hierarchy-discriminability`,
  `raw-fc-phase-distance`, `drift-triangle-null`, `cohort-coverage-matrix`,
  `decision-rules`, `measure-correctness-audit`, `e1-spectral-subspace-alignment`);
  and the eigenmode/Grassmann side-branches (`perspective-multidirectional-multiscale`,
  `eigenmode_embedding`, `epi-grassmann-embedding`). Superseded by the
  continuous cophenetic-trace + matched-strength approach.
- **`archive/kc-era/` (10)** — the May Kendall–Colijn / CTM / RF tree-distance
  family (`lrg-trace-investigation-plan`, `rf-k-multiscale-measure`,
  `kc-anchor-modules`, `kc-rearrangement-modules`, `kc-reset-modules`,
  `kc-leaf-topological-memory`, `taxonomy_mutual_exclusivity`,
  `kc-partition-merge-node-trace`, `epi-cross-phase-rigidity`,
  `band_agnostic_lrg`). KC retired 2026-05-18; the cophenetic ultrametric replaced it.
- **`archive/epi-negatives/` (4)** — negative per-node epi-marker attempts
  (`epi-eigenmode-localization`, `epi-node-trace-marker`, `occult-epi-node-marker`,
  `signed-magnetic-laplacian-epi`). Distinct from the *shipped* R3 marker, which
  is the relational diffusion-community propagator (see
  `.agents/preprint/headlines/03_epileptogenic_markers.md` + `epi-marker-analysis.md`).
- **`archive/replay-n4/` (4)** — the pulled replay headline (former N4)
  (`replay-states-multiscale-reinstatement`, `replay-states-v3-scale-target-pair`,
  `replay-sustained-and-timescale`, `intrinsic-scale-reorganization`). Removed
  from the paper 2026-06-23 (transient/bursty replay cohort-negative; sustained
  reinstatement folded into N1). Restore point only.
- **`archive/opaque-null-era/` (2)** — rejected matrix-level Haar-rotation nulls
  (`on-manifold-coherency-surrogate`, `coordinated-cross-phase-null`). Retired:
  opaque, construction-dependent, run before synthetic validation. Matched-strength
  is the default null. Postmortem: `scripts/archive/2026-06_opaque-matrix-nulls/`.

> **Note.** The R3 epileptogenic-marker family (the shipped relational-propagator
> SOZ detector) has **no scope report in this folder** — its process record lives
> in `.agents/preprint/headlines/epi-marker-analysis.md`. Only the *negative*
> per-node attempts (`archive/epi-negatives/`) were ever scoped here.

## Rule — every new task-trace tool lands here first

When developing a new descriptive multiscale measure for the task-trace question:

1. Write the scope report **first**. No code before scope.
2. Mathematical rigor is non-negotiable — every definition, every predicate,
   every aggregator stated as a formula. A reader must be able to spot
   errors from the math alone.
3. Cite primitives in `lrg_eegfc.utils.metrics.tree*` and
   `lrg_eegfc.utils.metrics.hypothesis`. No private re-implementations.
4. Cross-reference what the measure complements or supersedes.
5. Add a row to the **Live measures** table above (date, link, headline backed, tagline).
6. State explicitly what the measure does **not** measure — the negative
   property is part of the contract.

This rule is mirrored in the project-level enforcement list:
`.agents/guides/04_rules/never-always-list.md` and the seed list at the top
of `CLAUDE.md` / `AGENTS.md`.

## Cross-references

- Paper headline architecture (where these measures land): `.agents/preprint/headlines/README.md`
- Current results workflow: `.agents/reports/2026-07-04_results-readiness-map.md`
- Tree primitives: `src/lrg_eegfc/utils/metrics/tree.py`,
  `src/lrg_eegfc/utils/metrics/tree_distance.py`
- Hypothesis-test helpers (use, do not reinvent):
  `src/lrg_eegfc/utils/metrics/hypothesis.py`
- ImCoh method note: `.agents/guides/02_methods/imcoh-guide.md`
- Probe-bias controls: `.agents/guides/02_methods/probe-bias-guide.md`
