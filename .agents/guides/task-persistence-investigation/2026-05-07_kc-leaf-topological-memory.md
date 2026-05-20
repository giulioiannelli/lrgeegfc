---
name: kc-leaf-topological-memory
type: scope
era: COHORT_N10
status: current
created: 2026-05-07
updated: 2026-05-07
pointers:
  - .agents/reports/2026-05-07_kc-cross-phase-taxonomy.md
  - src/lrg_eegfc/utils/metrics/tree_distance.py
  - scripts/01_compute/audit/audit_53_kc_leaf_topological_memory.py
---

# KC leaf-topological memory

**Per-leaf, per-pair visualization of KC topology persistence (λ=0). Surfaces *which leaves remember their partners* across `task_test → rest_post` directly on each phase's dendrogram, and shows that the visual asymmetry between trace bands (β, α) and ergodic bands (θ, γ_h) is a per-pair phenomenon, not just a module-level one.** The figure is descriptive — no null, no test — designed to make the topological component of the KC scalar trace (`§5.2`) directly readable.

## Notation

- `n` — number of leaves (LRG nodes / sEEG contacts), aligned across phases.
- `Z_φ ∈ ℝ^{(n-1)×4}` — scipy linkage matrix for phase `φ ∈ {pre, tt, post}` at `τ_min = 1/λ_max`.
- `m_φ(i,j) ∈ ℕ` — depth from the root of the MRCA of leaves `i, j` in `Z_φ`'s tree (number of edges from root to MRCA; 0 at root, larger for shallower MRCAs / closer pairs). Computed via `lrg_eegfc.utils.metrics.tree_distance.kc_vectors(Z_φ)[0]`.
- `MRCA_φ(i,j) ∈ {n, …, 2n-2}` — internal cluster id of the MRCA in `Z_φ`.
- `δ ∈ ℕ_{≥0}` — integer tolerance on the topology match. Default `δ = 0` (exact match); relaxed only on inspection.
- `M ∈ ℕ_{>0}` — number of top partners drawn per focal leaf. Default `M = 10`.
- `s_\max ∈ (0, 1]` — MRCA-subtree-size cap as a fraction of `n`. Default `s_\max = 0.30` (admit only fine-grained matches; see Caveats — *Trivial-MRCA pathology*).
- `\text{size}_φ(i, j) := |\text{subtree under } MRCA_φ(i, j)|`, with `\text{size}_φ \in \{2, …, n\}`.
- `i^⋆` — focal leaf, default `argmax_i f_i` with smallest-id tiebreak; pluggable list for top-N focals in a future iteration.

## Definitions

### Topological-trace pair

For ordered phases `(pre, tt, post)`:

```
T(i,j) := ( |m_tt(i,j) − m_post(i,j)| ≤ δ )
       ∧ ( |m_pre(i,j) − m_tt(i,j)|  > δ )
       ∧ ( size_tt(i,j)  ≤ s_max · n )
       ∧ ( size_post(i,j) ≤ s_max · n )
```

Domain: unordered leaf pairs `i ≠ j`. Range: `{0, 1}`. The first conjunct says "the pair sits at the same topological depth in task and post"; the second says "rest-pre disagrees with task"; the size cap restricts admission to pairs whose MRCA defines a fine-grained subtree in BOTH tt and post (load-bearing — see Caveats §). With `δ = 0` we require strict integer equality; with `δ = 1` we allow off-by-one. Setting `s_max = 1` recovers the un-filtered predicate.

### Per-leaf persistence score

```
f_i := |{ j ∈ [n] \ {i} : T(i,j) = 1 }| / (n − 1)   ∈ [0, 1]
```

Range `[0, 1]`. `f_i = 0` ↔ no partner satisfies the predicate; `f_i = 1` ↔ every partner does.

### Top-`M` partners of the focal leaf

```
P_{i^⋆} := { j : T(i^⋆, j) = 1 }                         (predicate set)
P_{i^⋆}^{(M)} := top-M of P_{i^⋆} ranked by m_tt(i^⋆, j) descending
```

Larger `m_tt` ↔ deeper MRCA ↔ closer pair ↔ shorter highlighted path. Ranking by `m_tt` descending picks the *closest* partners first, which keeps the colored path lengths short and the visual readable.

## Properties

- **Range.** `f_i ∈ [0, 1]`; `m_tt(i,j) ∈ {0, …, n-1}` for binary trees.
- **Symmetry.** `T(i,j) = T(j,i)`; `f_i ≠ f_j` in general.
- **Invariance.** `T` and `f_i` are invariant under leaf-label permutations applied identically to the three trees, and under root-preserving rotations of internal nodes.
- **Identifiability.** Two pairs `(i,j)` and `(i,k)` with the same `m_tt` value are indistinguishable in the predicate; the partner ranking has ties, broken by partner leaf id.
- **What it cannot detect.** Height-only memory (λ=1 component): two pairs sharing MRCA depth but at different heights are indistinguishable here. The counterpart (height memory) is a separate measure deferred to a future scope.
- **Not a test.** `f_i` is a descriptive scalar per leaf; the figure does not compute a null, a p-value, or a cohort fraction. Band asymmetry is read by *inspection across small-multiples*, not by significance.

## Caveats & failure modes

- **Trivial-MRCA pathology (load-bearing — fixed by `s_max`).** Without the size cap, a leaf sitting in a tiny subtree absorbed at one near-root merge produces O(n) m-matched pairs all routed through that single anchor merge. First-render diagnostic on Pat_07 θ leaf 38: 97 of 115 candidate partners satisfied `m_tt = m_post = 3 ≠ m_pre`, all because the `(38, j)` MRCA in tt + post is the same near-root cluster (depth 3 from root, subtree size > 100). The visual saturated with whole-tree-spanning red paths and the figure read "ergodic θ shows stronger trace than β" — the opposite of intent. Mitigation: `s_max = 0.30` cap on the MRCA-subtree fraction in BOTH tt and post; collapses the pathology (Pat_07 θ → 2 partners survive) while preserving the genuine signal (Pat_07 β → 10 partners). The cap is now part of the canonical predicate; sensitivity to `s_max ∈ {0.15, 0.30, 1.0}` is recoverable via the `--max-mrca-frac` CLI flag.
- **Tolerance sensitivity.** Exact `δ = 0` on a 117-leaf tree is harsh. If `f_i ≡ 0` for the bulk of leaves, the figure is empty. Mitigation: report a sensitivity panel showing `mean f_i` at `δ ∈ {0, 1, 2}` for each (patient, band) cell; relax `δ` only after inspection.
- **Argmax ties.** `f_i` is rational with denominator `n-1`; ties at the maximum are common. Mitigation: smallest-id tiebreak; pluggable list lets us render top-N focals later.
- **Pat_03 (1024 Hz outlier).** `m_φ(i,j)` is dimensionless (integer count), so the patient's sampling-rate quirk does not enter; Pat_03 panels are interpretable on the same scale as the rest.
- **Path overlap.** With `M = 10` partners, intersecting paths still occur near a busy MRCA. Mitigation: alpha-blended LineCollection, ranking by `m_tt` favours short paths.
- **Trivial-MRCA pairs.** Pairs whose MRCA is the root (`m = 0`) cannot be "trace-persistent" in any visually-meaningful sense (their entire path traverses the whole tree). The predicate still admits them; they get demoted by the `m_tt` ranking.
- **Inter-phase leaf-order independence.** The dendrogram leaf-order in each phase is arbitrary (scipy default). The predicate operates on cluster ids, not draw-order; coloring is computed in cluster-id space and only mapped to draw-order when rendering. No alignment step is required.

## Pseudocode

```text
INPUT   Z_pre, Z_tt, Z_post   linkage matrices, n leaves
        delta ∈ ℕ_{≥0}        tolerance on topology match
        M ∈ ℕ_{>0}            top-M partners per focal

# 1. KC m-vectors per phase
m_pre  ← kc_vectors(Z_pre)[0]    array shape (n*(n-1)/2,)
m_tt   ← kc_vectors(Z_tt)[0]
m_post ← kc_vectors(Z_post)[0]

# 2. Pairwise predicate (with size cap)
size_tt   ← pair_mrca_subtree_sizes(Z_tt)              # condensed
size_post ← pair_mrca_subtree_sizes(Z_post)            # condensed
T ← (|m_tt − m_post| ≤ delta)
  ∧ (|m_pre − m_tt| > delta)
  ∧ (size_tt ≤ s_max · n)
  ∧ (size_post ≤ s_max · n)

# 3. Per-leaf score (squareform-aware iteration)
f ← zeros(n)
for each pair index k mapping to (i, j) with i < j:
    if T[k]: f[i] += 1; f[j] += 1
f /= (n − 1)

# 4. Focal leaf
i_star ← argmin { i : f[i] = max(f) }                # smallest-id tiebreak

# 5. Top-M partners of focal
P ← { j : T[index(i_star, j)] = 1 }
rank P by m_tt(i_star, j) descending, take first M
P_M ← those j

# 6. Path tracing on each phase's dendrogram
for φ in {pre, tt, post}:
    parent_φ ← _build_parent_height_depth(Z_φ)[0]
    for j in P_M:
        mrca ← MRCA(i_star, j) in Z_φ                # last shared ancestor
        path_i ← chain(i_star → mrca) in Z_φ          # cluster ids
        path_j ← chain(j → mrca) in Z_φ
        record (path_i ∪ path_j, mrca) for rendering

# 7. Render: dendrogram + per-leaf histogram strip + colored paths
#    one panel per phase; focal bar marked, partner bars colored.
```

## Visualization spec

- **Layout (one figure = one (patient, band) cell):** 2 rows × 3 columns. Top row = dendrograms at `pre / tt / post`. Bottom row = per-leaf histogram strip (one bar per leaf, x = leaf order from top dendrogram, y = `f_i ∈ [0, 1]`). Histogram width matches dendrogram width column-wise.
- **Per-phase coloring on the dendrogram:**
    - Base dendrogram: light gray (`#bababa` lines, `#c8c8c8` leaf bars).
    - Focal leaf bar: saturated accent (hex `#d62728`), with caret marker `▼` above the bar.
    - Top-`M` partner leaf bars: same accent, slightly desaturated (60% saturation).
    - Path edges from `i^⋆` and from each `j ∈ P_M` up to their respective MRCA: accent color, alpha-blended LineCollection, line-width 1.4 px.
- **Histogram strip:** y-axis `[0, 1]` linear; bars colored gray by default, focal bar accent-saturated, partner bars accent-desaturated. Light horizontal grid at `f_i ∈ {0.25, 0.5, 0.75}`.
- **Reading rules:**
    - On `tt` and `post`, the colored paths from `i^⋆` to its top-`M` partners traverse the same number of edges by predicate. Visually: the colored region's depth profile is similar.
    - On `pre`, the same partners trace different-depth paths (by predicate). Visually: the colored region's depth profile diverges.
    - In a trace band, the histogram has a tall peak at `i^⋆` and a non-trivial tail of partner bars; in an ergodic band the histogram is uniform-low and the focal accent has no anchored neighbors.

## Connection to prior tools

| Prior tool | Relation |
|---|---|
| §5.2 KC scalar trace (audit_36 / audit_46) | Decomposes the SCALAR `Δ_KC^{tt-post}` against `Δ_KC^{pre-post}` into per-leaf, per-pair contributions. The scalar pools all `(i,j)`; this tool localises. |
| audit_47 KC trace-module view | Module-level (subtree-Jaccard gate); this tool is leaf-level (no gate, continuous score). The two views are complementary: modules pick *coherent leaf sets that move together*; this tool picks *individual leaves with persistent partner sets*. |
| audit_51b §5.6 cohort composite | Cohort-level taxonomy at module granularity. This is the cohort-level **per-leaf** picture; both views can sit in §5 to argue "trace lives at multiple granularities". |
| MRL / CBR (parked) | Subtree-identity vs continuous distance. This tool is per-pair-discrete (KC `m`), not per-pair-continuous (CBR rank). |
| RF(k) (parked) | Bipartition-based; this is path-based. |

What it does **not** replace: the §5.2 KC scalar tests significance; this tool surfaces the signal but does not test it. What it does **not** capture: height memory (λ=1), phase-symmetric effects (reset / anchor), and inter-leaf reorganization where partner identity changes but module identity is preserved.

## Implementation plan

- **Script:** `scripts/01_compute/audit/audit_53_kc_leaf_topological_memory.py`.
- **Library reuse:** `lrg_eegfc.utils.metrics.tree_distance.kc_vectors` for `m`; `_build_parent_height_depth` for parent/height arrays; `_ancestor_chain` for path walk to MRCA. No new helpers in the library at first pass; if audit_53 reuses helpers a second time elsewhere, promote them per the second-use rule.
- **Outputs:**
    - `data/audit/kc_leaf_topological_memory/figures/{patient}__{band}__delta-{δ}__topM-{M}.pdf`
    - `data/audit/kc_leaf_topological_memory/per_leaf_scores.csv` (long: patient, band, leaf, f_i, is_focal, in_top_M)
- **Initial test cells:** Pat_07 β (multiscale trace flagship; 6 trace modules at the §5.6 showcase) vs Pat_07 θ (ergodic counterexample).
- **Pluggability:** focal-leaf list is a `list[int] | None` argument; `None` → argmax. Top-`M` is a positional argument. Tolerance `δ` is a positional argument. δ = 0, M = 10 default.

## Open questions

- **`δ` default.** Start exact (`δ = 0`), inspect histogram density, decide if `δ = 1` is needed. Likely needed for low-`f_i` cells (γ_h, θ?).
- **`s_max` default.** Locked at 0.30 from the Pat_07 β/θ first render. Revisit if the cohort sweep produces empty cells at well-supported (β, β-α) cells; never raise to 1.0 (recovers the trivial-MRCA pathology).
- **Top-N focal small-multiples.** Once N=1 figure is read, N=3 (top-3 focals) gives a sense of cohort-internal heterogeneity. Layout becomes 6 columns (3 phases × 2 focal-leaf rows) or 3×3 (one focal per row). Defer until N=1 visual is locked.
- **Focal-anchored vs MRCA-anchored ranking.** `m_tt` ranking favors close pairs. Alternative: rank by `m_tt − m_pre` magnitude (the predicate "delta") to favor pairs that *changed* the most across the test. Defer to second iteration.
- ~~**Off-tree leaves.** Pairs with MRCA at root (`m_tt = 0`) are admitted by the predicate but visually empty (path = whole tree). Should they be excluded from the predicate set? Likely yes; revisit after first render.~~ **Resolved 2026-05-07**: subsumed by the `s_max` cap (root MRCA has subtree size = `n` > `s_max · n` for any `s_max < 1`).
