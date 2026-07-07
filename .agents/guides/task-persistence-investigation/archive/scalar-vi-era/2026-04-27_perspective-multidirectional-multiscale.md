---
name: perspective-multidirectional-multiscale
type: perspective
era: COHORT_N10
status: open
created: 2026-04-27
updated: 2026-04-27
pointers:
  - .agents/reports/2026-04-27_task-persistence-reconciliation.md
  - .agents/guides/task-persistence-investigation/2026-04-26_continuous-trace-matrix.md
  - .agents/guides/task-persistence-investigation/2026-04-25_cohesion-cbr.md
  - .agents/guides/task-persistence-investigation/2026-04-25_module-retention-landscape.md
  - .agents/guides/task-persistence-investigation/2026-04-26_multiscale-partition-coherence.md
---

# Perspective — multidirectional multiscale, what we measured vs what we should be measuring

**The trace question is intrinsically multidirectional in three independent
realms — topology, temporal, spatial — and most of our results come from
collapsing two of them onto an axis that has no physical meaning. The
inconsistencies across measures are not noise; they are diagnostic of the
collapse choice. Before any further analyses, we owe ourselves a deep
re-anchoring of the *cut parameter* to a physically meaningful quantity
(τ, the LRG diffusion timescale) and a multi-axis triangulation discipline
that refuses to claim cohort-wide structure unless every relevant axis
agrees. This file is the conceptual statement of that perspective; it
proposes nothing operational and contains no algorithms — its sole purpose
is to constrain the way we evaluate every existing and future measure.**

---

## 1. The three axes, and why our analyses keep colliding with axis-collapse

The trace claim is by construction multiscale and multidirectional. Three
independent physical/structural axes parametrise the search space:

| Axis | What it scales | What anchors cross-patient comparison |
|---|---|---|
| **Topological / hierarchical** | the LRG dendrogram (partitions, subtrees, ultrametric distances) | the diffusion timescale `τ` set by the network Laplacian |
| **Temporal / spectral** | the frequency band index `b ∈ {δ, θ, α, β, γ_l, γ_h}` | shared band definitions in Hz |
| **Spatial / anatomical** | the 3D layout of electrodes in the brain | physical millimetres + Desikan-Killany annotation |

Every concrete operationalisation we ran chooses *one* topology cut
parameter (integer `k`, fractional `h_rel`, or physical-scale `mm`),
fixes the spectral axis at a single band, and reduces spatial information
to a per-cluster summary. Different choices of cut parameter induce
different *cross-patient alignments* of the same underlying data — and
the inconsistencies we have been chasing are mostly inconsistencies of
the alignment, not of the underlying signal.

The pattern across all measures is consistent: any test that lets every
patient self-normalise (integer `k`, scalar tree distances) declares
universality across many bands; any test that forces cross-patient
alignment on a physically meaningful axis (h_rel, mm, KC) returns at most
two bands cohort-wide and never the same two. **Universality is an
artifact of forgiving alignment.** The honest bands are α and β under
controlled tests, with no cohort-wide anatomical convergence.

## 2. What we actually have

For brutal accounting we restate, per band, what each test currently says:

| | δ | θ | α | β | γ_l | γ_h |
|---|---|---|---|---|---|---|
| integer-k Δ_VI ≥ 8/10 (uncontrolled axis) | ✓ wide ridge | ✗ | ✓ mid-k | ✓ low-k fragile | ✗ | ✓ two ridges |
| Δ_VI vs drift floor (paired Wilcoxon) | pass | — | pass | **fail** | — | pass |
| H2c continuous controls (split + cross-probe + drift) | **fail** | fail | pass | pass | pass | **fail** |
| H2d block-pair (no real null) | universal | universal | universal | universal | universal | universal |
| h_rel-fixed Δ_VI cohort-wide | dissolves | — | — | — | — | scattered |
| spatial-scale (mean cluster diameter) | ✓ ~20 mm | ✗ | ✗ | ✗ | ✗ | ✓ ~25 mm |
| Strict subtree identity (MRL J=0.9) | null | null | null | null | null | null |
| Cohesion-CBR per-patient | varies, not yet cross-checked vs VI claims |

Two bands pass *every* test that has a control: **α and β** under H2c
continuous controls (paired Wilcoxon p = 0.008 / 0.014 vs drift floor)
and at integer-k Δ_VI. The bands declared cohort-wide by the spatial
axis (δ, γ_h) are exactly the ones suspect for artefact contamination
(slow drift / EMG / saccades) and they fail the H2c continuous controls.
The integer-k headline that originally promoted δ as *the* central
result is now seen to be a consequence of self-normalising k-alignment.

## 3. Why every cut parameter we tried failed cross-patient comparison

| cut parameter | what it does | why it fails cohort-wide |
|---|---|---|
| **integer `k`** | cuts each tree at the threshold producing exactly k clusters | each patient's "k=28" is at a different fractional depth, different effective n_eff, different physical scale; cohort claim is "every patient agrees at *their own* nominal cluster count", not at a shared scale |
| **`h_rel = h / dmax`** | cuts at fixed fractional depth | dmax is constant 0.9901 by LRG normalization, so h_rel is just `h` rescaled; cross-patient `h_rel(k)` spans 0.45–0.77 at fixed `k`; γ_h `dmin/dmax` reaches 0.71 in some patients (rake topology) so fine-scale comparison is structurally degenerate |
| **mean cluster diameter (mm)** | aligns by physical module size | trees are non-monotone in diameter; patients have heterogeneous diameter envelopes; cohort denominator non-monotonic by bin; only extremal bands (δ, γ_h) survive — exactly the suspect ones |

Common feature: **none of these is anchored to the diffusion physics
that built the dendrogram**. They are post-hoc reparametrisations of
the linkage matrix.

## 4. τ as the natural axis we lost

The LRG dendrogram's merge heights are not free parameters of the
clustering. They are derived from the network's Laplacian-diffusion
entropy `S(τ)` and the corresponding stability functional `C(τ) =
-dS/d(log τ)`. The cophenetic distance d(i, j) between two contacts in
the ultrametric is a strictly monotone function of the diffusion time
`τ` at which their associated random-walk equilibration domains merge.
**Cutting at `τ` is the only cut parameter with intrinsic physical
meaning that is approximately cross-patient comparable.**

What the LRG pipeline currently does that breaks τ-comparability:

```
linkage_matrix = linkage(dists / tmax, method)
tmax = max_distance + 0.01 * max_distance
```

The division by `tmax` rescales every tree to dmax ≈ 0.9901, which is
exactly what destroys absolute τ. After this normalisation, "h_rel"
is the only thing we have left, and h_rel does not align across
patients (different trees have different merge-height *distributions*
inside [0, 1]).

Two ways to recover τ:

- **Drop the normalisation** at the linkage step and propagate absolute
  merge heights through the pipeline. Re-derive every analysis on
  un-normalised trees. Costly but methodologically honest.
- **Bypass the linkage entirely** and work with the continuous diffusion
  observables `S(τ)` and `C(τ)` directly, which are the LRG framework's
  *primary* outputs. Phase-comparisons at the level of τ-resolved
  observables (e.g., `S^post(τ) − S^test(τ)` averaged over τ-windows;
  `C^post(τ) − C^pre(τ)`) test the same trace claim on the original
  diffusion physics, before any clustering discretisation. This is
  conceptually cleaner.

We have not done either. Until we do, no cut-parameter axis we
construct will be cohort-comparable in a defensible way. *This is the
single largest open methodological problem.*

## 5. KC λ-sweep — the topology ↔ metric memory dichotomy

Kendall–Colijn 2016 introduces a parametric tree distance
`d_KC(λ; T1, T2) = ||(1−λ) m + λ M||` where `m` carries pure topology
information (LCA edge counts) and `M` carries pure metric information
(branch lengths between LCAs). The two limits:

- `λ = 0`: agreement iff topology matches → **memory of *which* leaves
  cluster together** — what CBR / Cohesion-CBR captures.
- `λ = 1`: agreement iff branch lengths match → **memory of pairwise
  cophenetic distances** — what H2c continuous tests capture.

Sweeping λ on directed phase-pair distances would directly probe which
kind of memory dominates the trace claim per band. This is conceptually
the bridge we have been missing between the discrete-CBR and the
continuous-H2c tests; KC was previously implemented as a single scalar
distance only and dismissed as stale-numeric. The λ-sweep was never
attempted. **Whether the trace lives in topology, metric, or a mixture
is an open question that KC alone can answer cleanly.**

## 6. Cross-realm triangulation discipline (going forward)

No measure passes alone. The methodological compass for any future
result:

1. **Per-patient visual verification.** Every cohort claim must be
   verifiable on actual dendrograms with the trace-relevant subtrees
   highlighted. If the eye cannot see it on at least 8/10 patient
   panels, the cohort number is statistical noise inside an unanchored
   coordinate system.

2. **Multi-axis robustness.** A claim is acceptable iff it survives:
   - the topology-side test (CBR, MRL, KC λ ≈ 0),
   - the metric-side test (H2c controlled, KC λ ≈ 1),
   - some form of cross-patient alignment (mm, τ, percentile-of-implant),
   - per-band bootstrap and LOO leave-out.
   No "majority of 4/6 axes" — *every* axis must agree.

3. **No claim from an uncontrolled axis.** The integer-k headline at
   ≥ 8/10 was an uncontrolled-axis claim. We do not repeat that
   mistake. Self-normalising axes give descriptive material only,
   never headline cohort claims.

4. **Inconsistency is data.** When two measures disagree, the
   disagreement is a *finding*, not a problem to paper over. The
   disagreement between integer-k Δ_VI and h_rel-fixed Δ_VI is the
   most important finding of the past two days: the trace is
   partition-resolution-locked, not scale-locked.

5. **Strict reviewer mode for memory and prior claims.** We do not
   trust auto-memory entries. Every reuse of a previously-stated
   number is verified against the current cache. Memories that say
   "X was true at era Y" are pointers, not facts.

## 7. The δ / γ_h trap

The spatial-axis result that promoted δ and γ_h as the only cohort-wide
bands is a methodological trap, not a discovery. δ and γ_h are exactly
the bands routinely contaminated by:

- δ: slow drift, vigilance, infraslow trends, session-order effects.
- γ_h: EMG, saccades, electrode noise, line harmonics.

When a method's cohort-wide pattern is exactly the shape of the
dominant artefact, the conservative reading is that the method is
sensitive to the artefact, not that the artefact is the trace.
Cognitive bands α and β should be where the task signal lives;
they pass the controlled continuous test (H2c) and the integer-k test
but they do not survive the spatial axis. **The interpretation that
fits every test simultaneously is "α and β reorganization is real but
patient-heterogeneous; δ and γ_h are partly artefact"**, not "δ is the
true cohort headline".

## 8. What we won't do until we resolve the conceptual issues

- We won't add new tests to the table.
- We won't run new operationalisations of "the trace".
- We won't promote the spatial-axis δ / γ_h result without a τ-anchored
  re-test.
- We won't write the paper section that depends on uncontrolled-axis
  cohort claims.

## 9. What we owe ourselves to think about, in order

1. **Recover τ.** Either undo the dmax normalisation and re-run a
   minimal subset of analyses on absolute heights, or pivot to the
   pre-linkage `S(τ)` / `C(τ)` continuous observables.
2. **KC λ-sweep on directed phase-pair distances** (`rest_pre→rest_post`
   vs `task_test→rest_post`). One scalar per (patient, band, λ); cohort
   inference per band; identify the λ where the trace is best resolved.
3. **Cross-check Cohesion-CBR vs Δ_VI per-patient.** For every (patient,
   band) cell where Δ_VI > 0, does the corresponding Cohesion-CBR
   per-leaf figure show a non-trivial trace fraction? If not, what does
   the discrepancy tell us about the difference between partition-
   divergence and subtree-identity memory?
4. **Stratify by implant geometry.** Are there sub-cohorts of patients
   with similar coverage that show *consistent* α / β trace at the
   same scale, while heterogeneous-coverage patients dilute the
   cohort? This is a covariate we have not used.
5. **Explicit statistical-power discipline.** At n=10 with effect sizes
   r_rb ≈ 0.3–0.5, what is the maximum we can defensibly claim? The
   "≥ 8/10" cohort gate at α = 0.05 with cell-level FDR has known
   limits. Claiming more than the data allows is the easiest way to
   produce inconsistent results.

## 10. What this perspective is *not*

It is not a plan. It does not commit us to running anything. It is a
methodological constraint on what we will accept as a result and how
we will read inconsistencies between measures. Future scope reports
under this folder must conform to it; results that do not meet the
multi-axis robustness bar are descriptive notes only, not cohort
claims.

The next session should pick up here, not produce more measures.
