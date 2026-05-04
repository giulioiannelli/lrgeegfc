---
name: cophenetic-neighbourhood
type: scope
era: COHORT_N10
status: draft
created: 2026-04-25
updated: 2026-04-25
pointers:
  - .agents/guides/task-persistence-investigation/README.md
  - .agents/guides/task-persistence-investigation/2026-04-25_cohesion-cbr.md
  - .agents/guides/task-persistence-investigation/2026-04-25_cbr-classifier-failure-analysis.md
  - .agents/reports/2026-04-24_multiscale-task-trace.md
  - src/lrg_eegfc/utils/metrics/tree.py
---

# Cophenetic-Neighbourhood Persistence (CNP)

**Per-leaf classifier on the cophenetic distance vectors `D^φ_ℓ ∈ ℝ^N`
each leaf ℓ has in each phase φ. Asks the *hierarchy-level* question
"does ℓ keep its neighbours from task into rsPost, and were those
neighbours different in rsPre?". Sanity check 2 confirmed the FC
Laplacian spectra have **no clean spectral gap** in *any* (patient,
band) cell — the LRG dendrograms encode a continuous hierarchy, not
discrete modules. CNP is the natural per-leaf measure on a continuous
hierarchy: no SIZE_MIN floor, no anchor scan, no subtree-membership
ambiguity. Every leaf gets four soft affinities (trace / persist /
reset / rearrange) from corner-distance on the
`(ρ_pre_post, ρ_task_post)` plane. Local version of H2c.**

## 1. Notation

Inherits §1 of `2026-04-25_cohesion-cbr.md`. Recap:

- `P, B, Φ` — `COHORT_N10` patients, 6 bands, 4 phases.
- `T^{p,b,φ}` — binary tree induced by the LRG linkage on `imcoh_abs`,
  with `N_p` leaves and `N_p − 1` internal nodes.
- `Z^{p,b,φ}` — scipy linkage matrix.

For every (patient, band, phase) and every ordered pair of leaves
`(ℓ, i)`, define the **cophenetic distance**

    d^{p,b,φ}(ℓ, i) := h(LCA_{T^{p,b,φ}}({ℓ, i}))

i.e., the merge height at which ℓ and i first co-cluster in
`T^{p,b,φ}`. By construction `d^{p,b,φ}(ℓ, ℓ) = 0` and `d^{p,b,φ}` is
the unique ultrametric induced by the linkage. Implementation:
`scipy.cluster.hierarchy.cophenet(Z)` returns the condensed form;
`scipy.spatial.distance.squareform(...)` converts to a full N × N
matrix.

The **cophenetic-neighbourhood vector** of leaf ℓ in phase φ is the row

    D^{p,b,φ}_ℓ := (d^{p,b,φ}(ℓ, 1), …, d^{p,b,φ}(ℓ, N_p)) ∈ ℝ^{N_p}.

The diagonal entry `D^{p,b,φ}_ℓ[ℓ] = 0` is excluded from any
correlation involving `D^{p,b,φ}_ℓ`.

## 2. Definitions

### 2.1 Per-leaf hierarchy-similarity scores

For phases `a, b ∈ Φ` and leaf `ℓ`,

    ρ^{p,b}_{a,b}(ℓ) := Spearman(D^{p,b,a}_ℓ \ {ℓ}, D^{p,b,b}_ℓ \ {ℓ})  ∈ [-1, 1].

Spearman (rank-based) rather than Pearson because the cophenetic
distances span orders of magnitude (early merges are tiny, root merges
are large) and we care about the *order* of ℓ's neighbours, not
absolute distances.

The three load-bearing scores:

    ρ_pre_post(ℓ)  := ρ^{p,b}_{rsPre, rsPost}(ℓ)         ℓ keeps rsPre neighbours into rsPost
    ρ_task_post(ℓ) := ρ^{p,b}_{task, rsPost}(ℓ)          ℓ keeps task neighbours into rsPost
    ρ_pre_task(ℓ)  := ρ^{p,b}_{rsPre, task}(ℓ)           reported as a sanity sidecar

with the **task collapse**

    D^{p,b,task}_ℓ := (D^{p,b,task_learn}_ℓ + D^{p,b,task_test}_ℓ) / 2

— element-wise mean of the two task-phase cophenetic vectors. If
`task_test` is missing for the patient (Pat_14 case until 2026-04-25
vendor replacement; n/a today since cohort returned to N10), use
`task_learn` alone.

Note that `ρ_pre_post`, `ρ_task_post` are *not* anchored on rsPost in
the same sense as the CBR family — they are **leaf-symmetric** by
construction (the cophenetic vector of leaf ℓ is its own perspective).
The four soft affinities below are still anchored on the (rsPost vs
rsPre, rsPost vs task) decomposition because that is the scientific
contrast we care about.

### 2.2 Soft pattern affinities

Each leaf has a point `(ρ_pre_post(ℓ), ρ_task_post(ℓ)) ∈ [-1, 1] ×
[-1, 1]`. The four corners encode the four CBR patterns:

    TRACE corner     := (-1, +1)        ℓ moved away from rsPre, kept task-buddy
    PERSIST corner   := (+1, +1)        ℓ kept everyone
    RESET corner     := (+1, -1)        ℓ kept rsPre, lost task structure
    REARRANGE corner := (-1, -1)        ℓ kept no one

Affinity =  1 − distance to corner / max possible distance. The
geometry uses Euclidean distance on `[-1, 1]²` and divides by
`2√2` (the diameter of the square):

    a_p(ℓ) := 1 − ‖(ρ_pre_post(ℓ), ρ_task_post(ℓ)) − corner_p‖₂ / (2√2)
              ∈ [0, 1].

`p ∈ {trace, persist, reset, rearrange}`. Affinities are not mutually
exclusive and do not sum to 1 — they are independent axes.

### 2.3 Per-leaf classification

    dominant(ℓ)   := argmax_p a_p(ℓ).
    confidence(ℓ) := a_dominant(ℓ) − a_secondbest(ℓ)   ∈ [0, 1].

Colour rule for the figure:
- Hue: `PATTERN_COLOR[dominant(ℓ)]` (the user-pinned palette: trace
  red, persist black, reset gray, rearrange blue).
- Alpha: solid `1.0` for every leaf — the corner-distance geometry on
  `[-1, 1]²` is balanced (every corner is equidistant from its three
  non-opposite neighbours), so persistent under-saturation of trace
  (the pathology that bit cohesion-CBR's first pass) does not arise.
  Confidence is reported in the leaf CSV for downstream filtering.

## 3. Properties

| property | value |
|:---|:---|
| Range of each `ρ_*` | `[-1, 1]` |
| Range of each `a_p` | `[0, 1]` |
| Anchor-free | **Yes** — no rsPost subtree scan, no SIZE_MIN floor. |
| Robust to absence of spectral gap | **Yes** — works on continuous hierarchies (which sanity check 2 confirmed is what every (patient, band) cell has). |
| Robust to rake topology | **Yes** — Pat_10 δ rake-like dendrogram still gives meaningful per-leaf cophenetic vectors; the issue becomes "do the rake-leaf orderings agree across phases?" rather than "do the modules survive?". |
| Symmetric in `(rsPre ↔ rsPost)` swap | **No** — the affinity geometry is anchored on rsPost (corners encode "ℓ ended up like task / rsPre / etc."). |
| Detects fragmentation | **Yes** by construction — fragmentation in task = ℓ's task neighbourhood differs from its rsPost neighbourhood = `ρ_task_post` low. |
| Identifiability | A leaf with `(0, 0)` on the affinity plane gets all four affinities ≈ 0.5 → low confidence → genuinely ambiguous, not silently classified. |
| Complexity | Cophenetic matrix is `O(N²)` per phase from `cophenet(Z)`. Three Spearman correlations per leaf, each `O(N log N)` → `O(N² log N)` per (patient, band). For N≈120: ~80k Spearman calls cohort-wide; trivial. |

**Negative property — what CNP cannot detect:**
- Anchored subtree identity. CNP doesn't say "*these specific leaves*
  formed a module that persisted"; it says "leaf ℓ's hierarchical
  neighbours stayed the same / changed". The two are related but not
  identical. CBR-family complements CNP for "name the module"
  questions.
- Per-scale resolution. CNP collapses across all merge heights via
  the Spearman rank — it cannot say "ℓ's local-scale neighbourhood
  persisted but its global-scale neighbourhood reorganised". A
  `h_rel`-windowed CNP variant (restrict the cophenetic vector to
  pairs with `d ≤ h_threshold`) is straightforward and out-of-scope.

## 4. Caveats & failure modes

| caveat | mitigation |
|:---|:---|
| Spearman is sensitive to ties. The cophenetic distance can have ties (multiple leaves merging at the same height in scipy's UPGMA). | scipy's `spearmanr` uses average-rank for ties — standard, documented. |
| Probe-bias inheritance. Same-probe contacts have trivially small cophenetic distance in every phase → boost ρ_* uniformly. | This biases CNP toward `persist` for probe-clustered leaves. The bias is real and *interpretable*: those leaves *are* always-near-each-other. Probe-debiased FC pipeline (out of scope) would address it. |
| Pat_10 δ stress test. Sanity 2 showed Pat_10 δ has the lowest spectral structure and dendrograms that diverge across phases. CNP should report low correlations everywhere → mostly rearrange. | Predicted; we explicitly test for it in §6. |
| Leaf-self correlation. The diagonal `D^φ_ℓ[ℓ] = 0` would inflate Spearman if included. | Explicitly excluded in §2.1. |
| Continuous correlations have no "limbo zone" but the geometry of the affinity simplex still has *near-corners* where confidence is low. | Solid alpha; confidence reported separately for downstream filtering. |

## 5. Pseudocode

```
INPUT  Z_rpre, Z_tl, Z_tt, Z_rpost  (linkage matrices, same N leaves)

D_rpre  ← cophenet_matrix(Z_rpre)        # N × N
D_tl    ← cophenet_matrix(Z_tl)
D_tt    ← cophenet_matrix(Z_tt)          (None if missing)
D_rpost ← cophenet_matrix(Z_rpost)
D_task  ← (D_tl + D_tt) / 2 IF D_tt ELSE D_tl

per_leaf ← []
FOR ℓ in {1..N}:
    mask     ← {1..N} \ {ℓ}              # exclude self-distance
    v_rpre   ← D_rpre[ℓ, mask]
    v_task   ← D_task[ℓ, mask]
    v_rpost  ← D_rpost[ℓ, mask]
    ρ_pp     ← spearmanr(v_rpre,  v_rpost)
    ρ_tp     ← spearmanr(v_task,  v_rpost)
    ρ_pt     ← spearmanr(v_rpre,  v_task)        # sanity sidecar
    a_trace     ← affinity(ρ_pp, ρ_tp, corner=(-1,  1))
    a_persist   ← affinity(ρ_pp, ρ_tp, corner=( 1,  1))
    a_reset     ← affinity(ρ_pp, ρ_tp, corner=( 1, -1))
    a_rearrange ← affinity(ρ_pp, ρ_tp, corner=(-1, -1))
    dominant    ← argmax_p (a_p)
    confidence  ← top1(a_p) − top2(a_p)
    APPEND (ℓ, ρ_pp, ρ_tp, ρ_pt, a_*, dominant, confidence)

OUTPUT per_leaf
```

`affinity(ρ_pp, ρ_tp, corner)` returns `1 - sqrt((ρ_pp - corner[0])² +
(ρ_tp - corner[1])²) / (2*sqrt(2))`.

## 6. Visualization spec

One PDF per (patient, band) at
`data/audit/per_patient_hierarchy_cnp/{Pat_XX}/{band}_hierarchy-crossphase.pdf`.

Layout — *identical* to audit_08 / cohesion-CBR (4 phase dendrograms +
codebar) so all measures are directly comparable. The only meaningful
difference is the per-leaf input: each leaf's colour comes from its
CNP `dominant`, alpha is solid `1.0`.

Branches: an internal-node link gets the colour of its dominant
category iff every leaf below it has the same dominant category, else
background gray.

Auxiliary CSVs:
- `leaf_assignment.csv` columns:
  `patient, band, leaf_id, dominant, confidence,
   rho_pre_post, rho_task_post, rho_pre_task,
   a_trace, a_persist, a_reset, a_rearrange`.
- No `anchor_classification.csv` — CNP has no anchor scan.

**Reading rule.** Saturated red leaves = ℓ's hierarchical neighbours
in rsPost match its task neighbours and not its rsPre ones. Saturated
black leaves = ℓ's hierarchical neighbours never moved. The story is
*per leaf*; spatial clustering of same-colour leaves in rsPost
indicates modules whose hierarchical neighbourhoods all behaved the
same way (correlated trace/persist).

## 7. Connection to prior tools

| existing tool | how CNP relates |
|:---|:---|
| H2c (cophenetic ρ on Δ-vectors) | H2c is the *cohort-pooled, all-pairs* Spearman ρ between rsPre→task and rsPre→post ultrametric *shifts*. CNP is the *per-leaf, neighbourhood-vector* Spearman between phases. They share the same primitive (cophenetic distance) and will agree on direction by construction; CNP localises *which leaves* are doing the agreement. |
| H2d (block coactivation persistence) | H2d is per-pair-per-scale. CNP collapses across scales but works per-leaf. Complementary. |
| Cohesion-CBR (variant E) | Cohesion-CBR asks "does this *subtree* persist as a coherent group". CNP asks "does this *leaf*'s neighbourhood persist". The difference: cohesion-CBR requires a discrete subtree (failed on Pat_10 δ rake topology); CNP doesn't. Same figure layout; different per-leaf input. |
| MRL (Module-Retention Landscape) | MRL is the scalar-field version of CBR-trace. Disjoint from CNP — MRL counts modules, CNP describes leaves. |

## 8. Implementation plan

- **Library reuse**: `scipy.cluster.hierarchy.cophenet`,
  `scipy.spatial.distance.squareform`, `scipy.stats.spearmanr`. No
  new helpers in `lrg_eegfc.utils.metrics.tree`.
- **Script**: `scripts/01_compute/audit/audit_13_cnp.py`. Reuses
  `plot_dendrogram_with_strip` from `audit_08`. Per-leaf input is
  generated directly from cophenetic matrices; no anchor loop.
- **Outputs**:
  - `data/audit/per_patient_hierarchy_cnp/{Pat_XX}/{band}_hierarchy-crossphase.pdf`
  - `data/audit/per_patient_hierarchy_cnp/leaf_assignment.csv`
- **Pat_14**: cohort returned to N10 on 2026-04-25 with vendor
  replacement (per README frontmatter). If a phase is still missing
  for any patient, fall back to `D_task = D_tl`.
- **No promotion to library** in this pass; if a third caller appears
  the cophenetic loader goes to `src/lrg_eegfc/utils/metrics/tree.py`.

## 9. Open questions

- **Confidence-graded alpha**. Currently solid `1.0`; alternative is
  `0.4 + 0.6 · confidence` like cohesion-CBR's earlier pass. Decision
  after first visual review.
- **Task collapse via vector-mean vs correlation-mean**. Default is
  vector-mean (`D_task = (D_tl + D_tt) / 2`). Alternative is to
  compute `ρ_task_post = mean(ρ_taskL_post, ρ_taskT_post)`. Vector
  mean is simpler but blends timescales if `D_tl` and `D_tt` have
  different absolute height scales; correlation mean is invariant.
  Defer.
- **`h_rel`-windowed CNP**. Restrict the cophenetic vector to pairs
  with `d^φ_ℓ ≤ h_threshold` to get scale-resolved CNP. Separate
  scope.
- **Probe-debiased CNP**. Strip the contribution of same-probe
  pairs from the cophenetic vector before computing Spearman.
  Requires the probe annotation pipeline (out of scope for the first
  cut).
