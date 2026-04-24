---
name: modular-trace-investigation
type: report
era: COHORT_N9
status: dead
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

# Modular task-trace investigation — handoff notes (updated 2026-04-24, post-Stage-1)

## The question

Is there a **multiscale, cohort-wide, band-specific** structural trace of
task (learn+test) in the rest_post LRG dendrogram, detectable across all
9 patients and distinguishable from directional drift alone?

Scope locked: 9 cross-phase patients, FC = `imcoh_abs`, LRG caches at
`data/cache/imcoh_lrg/` (2026-04-22 snapshot). Do not recompute FC / LRG.

## Already locked at n=9 (do not re-run)

- **H1** (task-learn ≈ task-test): 6/6 bands q<0.02.
- **H2c** (ρ>0, ultrametric drift co-varies with task drift, diffusion level): 6/6 bands q<0.005, **includes the 4 patients that fail module-level tests**.
- **H2d** (pairwise co-cluster persistence Δρ>0, discrete scales): 6/6 bands q<0.005.
- **H2d-θ** (θ has least block memory): θ<α p=0.010, θ<δ p=0.049.
- **H2e** (drift-floor split-half): 6/6 bands q=0.004.
- **H3** (rest vs task is a real axis): 6/6 bands q<0.05.

All unpaired partition-distance tests (H2a family) failed at n=9.

## Dead ends from the current investigation — report as negative evidence, do not re-open

| approach | what it tested | outcome |
|:---------|:---------------|:--------|
| CBR leaf-Jaccard + π + size-balanced weight, 4-pattern classifier | is the SAME CHANNEL SET grouped into a module in rpost-and-tt but not rpre? | band-selective (β, γ_l, α) but only **5/9 patients** contribute TRACE modules at any size. 4 patients have best-J_tt ≤ 0.38 — truly absent, not near-miss. Gallery of 85 PDFs at `data/outputs/figures/cbr_gallery/` documents the phenomenon but also its patient-sparseness. |
| Restricted cophenetic Spearman per module | does the INTERNAL BRANCHING of a matched module agree across phases? | previously considered as topology add-on; rejected by the user as insufficient structural evidence on its own (too local, still tied to leaf-set match). |
| Raw ultrametric matrix direct comparison (H2a-RAW) | does D^rpost align with D^tt more than with D^rpre? | failed 0/6 at n=9. |
| Frobenius ratio (H2a-FROB) | same, via ratio-of-norms | failed, θ in wrong direction. |
| Partition-level VI/H/NMI directed contrasts at each k | scale-resolved partition similarity | cluster-perm only finds δ k=23-31 and α k=2-4 at p<0.05. Visual unanimity pattern band-selective but sub-threshold on aggregate. |

All of these operate either on (a) raw pairwise distances or (b) leaf-set
matching at specific scales. Neither captures the **global topological
arrangement of the tree as a multiscale object**.

## 2026-04-24 — Stage 1 (14-metric rerun under `imcoh_abs`): FAILS pre-reg, β-only signal

Rerun of the 14 MSC-era tree-comparison scalars (CophPearson/Spearman,
BakersGamma, NormL1/L2 on ultrametric, TopK5/10/20 merges, WeightedARI
coarse/fine, MeanARI, WeightedVI coarse/fine, MeanVI) under post-reset
`imcoh_abs` LRG caches, 9 patients × 6 bands × 3 phase pairs. H2a contrast
`s = sim(tt, post) − sim(pre, post)` with Wilcoxon + FDR-BH(m=6) + rb +
10k-boot CI. See:
- Driver: `scripts/01_compute/stage1_14metrics_rerun.py`
- Inventory: `data/reports/imcoh_vi/stage0a_metric_catalog.md`
- Stats + report: `data/reports/imcoh_vi/stage1_14metrics_rerun.md`
- Figures: `data/reports/imcoh_vi/figures/stage1_14metrics_{heatmap,profiles}.pdf`

**Outcome — none of the 14 metrics meet the pre-registered criterion.**
Best metrics (NormL1_Ult, NormL2_Ult, MeanARI, WeightARI_Fine) reach 2 of 4
criteria. The failing criteria are c1 (β AND γ_l q<0.05) and c2 (β AND γ_l
n≥7/9). The passing criteria are c3 (θ AND α null, q>0.1) and c4
(rb(β) − rb(θ) ≥ 0.5).

**Key asymmetry discovered**: β behaves as hypothesized (n_positive 6-8/9,
rb_beta ≈ +0.5 to +0.7, p-values 0.03-0.16 pre-FDR), but **γ_l does not**
(n_positive 3-6/9, rb often negative). The band-specific hypothesis
"β AND γ_l" is not supported at single-scalar resolution — only β carries
directional trace signal in these metrics.

**dmax stability**: the LRG ultrametric is normalized so `dmax ≈ 0.99`
across every (patient, band, phase) — fractional cuts (`h_rel = h/dmax`)
are meaningfully comparable across phases without further normalization.
Diagnostic: `data/reports/imcoh_vi/stage0b_dmax.md`.

**Decision**: per simple-first strategy, advance to Stage 2 (literature
search for tree-similarity scalars that combine heights AND topology —
the category missing from the 14). If Stage 2/3 still find nothing on β
that survives FDR, the paper either (a) revises hypothesis to "β only"
with narrower claim + scale-resolved supplement, or (b) stops at Stage 4's
scale-resolved band × scale map for per-cell evidence.

## 2026-04-24 — Stage 2/3 (KC + MC + wRF scalar tree-distance tests): FAILS pre-reg, **key heights-vs-topology diagnostic**

Stage 2 literature menu (`.agents/reports/2026-04-24_stage2-literature-menu.md`)
identified Kendall-Colijn, Matching Cluster, and weighted Robinson-Foulds
as the canonical tree-distance scalars that blend branch-lengths AND
topology — the category the 14 lacked. Implementation in
`src/lrg_eegfc/utils/metrics/tree_distance.py` (KC, MC, wRF).

Stage 3 driver: `scripts/01_compute/stage3_tree_distance.py`. Eight metric
variants run on 9 patients × 6 bands: KC at λ ∈ {0.00, 0.25, 0.50, 0.75,
1.00}, MC unweighted and weighted, wRF. Outputs:
`data/reports/imcoh_vi/stage3_tree_distance{.md,_stats.csv,_verdict.csv,
_raw.csv,_contrast.csv}` + figures incl. `kc_lambda.pdf`.

**Outcome — none of the 8 variants meet the pre-registered criterion.**
Best are at 2/4 (c3 θ/α null + c4 rb separation), same ceiling as Stage 1.

**Key scientific finding — KC(λ) diagnostic**: sweeping λ from 0 (pure
topology) to 1 (pure heights) shows the β-vs-θ rank-biserial separation
grows from 0.09 at λ=0 to 1.46 at λ=1. At pure topology (λ=0), θ actually
goes in the WRONG direction (rb_θ=+0.38, i.e., θ's rpost is closer to tt
than to rpre in topology only — noise level). At pure heights (λ=1),
θ is strongly negative (rb_θ=−0.82) and β is positive (rb_β=+0.64).

Implication: **at scalar resolution, the task trace signal lives in
continuous merge heights, not in branching topology.** Topology-based
scalars tell us nothing about the β-trace; height-based scalars see it
but are redundant with H2c. The 14-metric Stage 1 and topology-blend KC
Stage 3 agree: n=9 with FDR m=6 is underpowered even where rb~+0.65.

**Secondary finding**: wRF shows an unexpected **δ-band** signal
(rb=+0.82, 7/9 positive, p=0.014 pre-FDR, q=0.085). Not in the hypothesis
(β/γ_l) but worth noting as a bipartition-scale phenomenon. wRF is
correspondingly null on β (rb=−0.07) — its subtree-relocation sensitivity
hurts in the regime where KC/MC see something.

**Consequence for Stage 4 design**: Stage 4 was originally planned as a
*topology-based* scale-resolved band × scale map (each cell = proportion
of patients with a qualifying task-retained module at that scale). Stage
3 shows the scalar topology signal is null in β. So Stage 4 cannot be
a pure topology-at-multiple-scales analysis — it has to include heights
or become a scale-resolved cophenetic / H2c-refinement. Decision point
deferred to user.

## 2026-04-24 — Stage 4a (KC λ-sweep characterization): per-band signatures, still NOT cohort-wide

Fine 21-point λ grid over KC distances, 9 patients × 6 bands × 3 phase pairs.
Script: `scripts/01_compute/stage4_kc_exploration.py`. Outputs:
`data/reports/imcoh_vi/stage4_kc_{exploration.md, stats.csv, contrasts.csv,
distances.csv}` + 3 figures incl. `stage4_kc_lambda_per_band.pdf`.

Per-band λ peak summary (H2a rank-biserial):

| band | peak λ | peak rb | n_pos@peak | p@peak |
|:---|:---|:---|:---|:---|
| δ | 0.00 | +0.20 | 5/9 | 0.30 |
| θ | 0.40 | −0.07 | 4/9 | 0.57 |
| α | 0.55 | +0.16 | 5/9 | 0.34 |
| β | **0.30** | **+0.69** | **6/9** | **0.033** |
| γ_l | 0.15 | +0.56 | 7/9 | 0.069 |
| γ_h | 0.55 | +0.69 | 6/9 | 0.033 |

Per-band λ-profiles differ: γ_l peaks in the topology half, γ_h in the
blend, β is robust across much of the axis, θ is anti-trace at both
ends, δ/α weak. Descriptively informative but **still capped at 6–7/9
patients** — same ceiling as Stage 1.

## 2026-04-24 — STOP: do not run more scalar tests, surface existing VI(k) signal instead

User feedback (verbatim, 2026-04-24): *"all the signals were there with the
imcoherence framework and the VI(k) profiles... i see that there are traces
of persistence. we just have to find the way to make them emerge and you
are making me losing time."*

At n=9 no scalar measure — the 14 from Stage 1, KC/MC/wRF from Stage 3, or
the KC λ-sweep from Stage 4a — produces a cohort-wide (9/9 or ≥7/9
unanimous) band-specific β/γ_l trace. Best is 6–7/9 with Pat_15 consistently
anti-trace. Hypothesis-gated pre-reg is the wrong frame for this cohort
size. The signal was already present in the pre-existing |ImCoh| VI(k)
partition-multiscale analyses; further scalar work should be halted and
the next session should instead extract and characterize what those
existing artifacts already show.

**Do NOT (next session)**: launch Stage 4b, Stage 5, or any new scalar
cohort test with strict pre-reg FDR at q<0.05 m=6. That bar is unreachable
at this cohort size for effect sizes of rb ≈ 0.5–0.7.

**Do instead (next session)**: read the existing ImCoh artifacts and
surface what's already present. Priority order:
1. `.agents/reports/2026-04-24_h1-h4-vi-results.md` — canonical VI(k)
   unanimity tables under `imcoh_abs`. Where does band selectivity appear?
2. `scripts/01_compute/h2_partition_multiscale.py` outputs +
   `h2_partition_cluster_perm.py` outputs — multiscale partition-distance
   landscape. δ k=23-31 (Δ_VI, p=0.014★) and α k=2-4 (Δ_H, p=0.050★)
   are the two significant clusters found there. Does the visual landscape
   differentiate bands more strongly than cluster-perm cohort test admits?
3. `scripts/01_compute/h2d_coactivation_persistence.py` + H2d-θ finding
   (θ has least block memory, θ<α p=0.010, θ<δ p=0.049) — the one
   genuine band-specific cohort-wide result we already have. How to
   present it as the paper's central structural band claim?
4. `scripts/01_compute/h2c_ultrametric_drift.py` — cohort-wide directional
   drift (9/9, 6/6 bands). Universal, not band-specific, but the
   structural-claim backbone.
5. `.agents/reports/2026-04-24_multiscale-task-trace.md` — writing-agent
   handoff. Read this to understand how the existing results are being
   framed for the paper; the scalar-test narrative we just ran did NOT
   plug into this. The fix is to rewrite from the existing cohort-wide
   results, not to generate more tests.

**Stage 1-4 artifacts from today (do not delete, but do not cite as new
results either)**: files under `data/reports/imcoh_vi/stage{0a,0b,1,3,4}*`
and `src/lrg_eegfc/utils/metrics/tree_distance.py`. The KC / MC / wRF
module is usable for future work but none of the Stage 1-4 results pass
the cohort criterion at n=9.

## The core problem, stated precisely

We have three orthogonal levels of description:

1. **Diffusion-ensemble level** (cophenetic distances): universal drift
   direction (H2c, H2d) — passes for all 9 patients in all 6 bands.
2. **Module-membership level** (Jaccard on leafsets): band-selective
   (β, γ_l, α) but patient-sparse — only 5/9 patients exhibit it.
3. **Topology level** (branching structure of the tree as a graph, not
   the distances or the memberships): **not yet measured**.

The paper cannot land a cohort-wide band-specific structural claim on
level 2 alone. We need a level-3 measure that:

- produces ONE scalar per (patient, band) — suitable for cohort-wide
  Wilcoxon + FDR at n=9;
- is inherently MULTISCALE — aggregates information from the full
  branching hierarchy, not just at a single k or h;
- is sensitive to TREE TOPOLOGY, not to leaf-set identity alone and
  not to the pairwise continuous distances (that's H2c's job);
- shows BAND-SPECIFICITY — β / γ_l pass, θ / α do not, with effect
  sizes separated;
- is PAIRED within patient (like H2d) — task-drift direction compared
  against each patient's own rest-drift baseline.

## Candidate directions worth exploring (none yet implemented)

The next agent should not commit to any of these before reviewing all.
Listed roughly in increasing ceremony.

### (a) Kendall-Colijn-style rank-depth metric

For each ordered pair (i,j) in the tree, compute `m(i,j)` = number of
internal nodes on the path from the root to their MRCA (or equivalently
the rank of the MRCA-merge in the height-ordered sequence of merges).
This gives a pair-vector of length N(N-1)/2 that encodes *only* the
branching pattern — invariant to absolute merge heights but sensitive to
the order of merges. Per phase, this is a scalar vector; comparison via
Spearman.

Scalar per (patient, band): `κ(pat, band) = Spearman(m^rpost, m^tt) −
Spearman(m^rpost, m^rpre)`. Paired within patient. No free parameter.

Strengths: purely topological, well-defined distance, uses full tree.
Open question: is it multiscale in the required sense? It depends on
where in the tree the differences sit — may need to be scale-stratified.

### (b) Soft bipartition overlap with prominence weighting

Every internal node of a tree defines a bipartition of the leaves. A
bipartition "matches" another if their symmetric Jaccard ≥ τ. Count
prominence-weighted matches:

`B(P, Q) = Σ_{β∈P} π(β) · 1[∃β'∈Q : J(β,β') ≥ τ] / Σ_{β∈P} π(β)`

Scalar per (patient, band): `B(rpost, tt) − B(rpost, rpre)`.
Multiscale by construction (π weights stable splits across the tree).
Free parameter τ (pre-register at 0.8 or 0.85).

Strengths: topology-level, prominence-weighted, fully tree-aware.
Risks: leaf-level still via Jaccard — a better version might use pair-
containment overlap instead of membership Jaccard.

### (c) Multiscale induced-subtree comparison

For each h-bin on a log grid, take the INDUCED tree on the channels that
are in a non-trivial cluster at that h (those above a min-size). Compare
the induced subtrees across phases via Robinson-Foulds or by descendant-
based Jaccard across internal nodes. Aggregate via π-weighted h-integration.

Strengths: genuinely multiscale, full topological comparison at each scale.
Risks: heavier to implement; null model less obvious.

### (d) Full-tree Wasserstein on merge-event point clouds

Each tree → a 2D point cloud {(h(η), size(η)) : η internal}. Compare
across phases via Wasserstein-1 distance. Task-trace scalar:
`W(rpost, tt) − W(rpost, rpre)`.

Strengths: distribution-level comparison of merge events, highly robust.
Risks: discards identity of leaves entirely — may overweight bulk
statistics over structurally meaningful alignments.

## Pass/fail for the "cohort-wide band phenomenon" claim

Before running ANY new measure, pre-register:

- For the candidate scalar `s(pat, band)`:
  - Paired Wilcoxon `s > 0` per band, FDR-BH m=6.
  - β AND γ_l at q<0.05, with at least 7/9 patients having s>0 per band.
  - θ AND α remain null (q>0.1, median s ≈ 0).
  - Effect size (rank-biserial) separation between "memory" and "no-
    memory" bands: rb(β) − rb(θ) ≥ 0.5.

If a candidate scalar meets all four, we have the cohort-wide band
phenomenon. If only 2–3 criteria met, the paper claim weakens
accordingly (narrower band set or "partial memory").

## Priorities for the next session

1. **Stop iterating on Jaccard-leafset classifiers.** They've given what
   they can — 5/9 patients in β/γ_l/α at small sizes. Keep the gallery
   as descriptive support, not as the primary claim.
2. **Implement candidate (a) — Kendall-Colijn rank-depth κ — first.**
   It is the cheapest, most principled topology-only measure, and
   directly testable at n=9 with the H2c/H2d statistical pipeline.
3. **If (a) fails, try (b) — soft bipartition overlap.** If that also
   fails, (c) and (d) sequentially.
4. **Do not commit the paper to a "module-level trace" claim.** Until
   a level-3 measure passes cohort-wide band testing, the paper's
   structural claim remains open.

## Files of record for the next session

- This file: `.agents/reports/2026-04-24_modular-trace-investigation.md`
- Compact pastable handoff: `.agents/reports/MODULAR_TRACE_COMPACT_MESSAGE.md`
- Writing context: `.agents/reports/2026-04-24_multiscale-task-trace.md`
- Method definitions: `.agents/guides/02_methods/h2-metrics.md`
- Era index: `.agents/reports/2026-04-24_pipeline-status.md`
- Gallery (descriptive, not statistical): `data/outputs/figures/cbr_gallery/`
- Tree utilities already in place: `src/lrg_eegfc/utils/metrics/tree.py`
  (`tree_internal_nodes`, `jaccard_leafsets`, `h_log_grid`, `dmax_from_Z`,
  `fcluster_at_h_rel`)
- Stats helpers: `scripts/01_compute/_shared.py` and
  `scripts/01_compute/h2d_coactivation_persistence.py` (`cluster_perm`).

## Guardrails that do not change

- Pat_14 excluded (corrupt task_test).
- Pat_03 flagged as 1024 Hz outlier.
- Pat_10 is 113-channel (drop task rows [53,54,55] at load).
- FC method: `imcoh_abs` only.
- No H2a-family revival. No pooling of VI/H/NMI into consensus scalars.
- All new metrics must be paired within-patient. Unpaired three-phase
  tests fail at n=9 — that is the main lesson of this investigation.
