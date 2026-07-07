---
name: kc-anchor-modules
type: scope
era: COHORT_N10
status: draft
created: 2026-05-07
updated: 2026-05-07
pointers:
  - .agents/reports/2026-05-07_kc-trace-module-visualization.md
  - .agents/guides/task-persistence-investigation/2026-05-07_kc-reset-modules.md
  - .agents/guides/task-persistence-investigation/2026-05-07_kc-rearrangement-modules.md
  - .agents/guides/01_project/terminology.md
---

# KC anchor modules — coherent across all three phases

**Head.** Fourth class in the cross-phase taxonomy. An anchor module
(user term: "persist") is a leaf set that exists as a coherent
subtree in ALL THREE phases — rest_pre, task, AND rest_post. This is
the **null** of the trace / reset / rearrangement story: anchors are
the modules that the task does NOT touch. Surfacing anchors is the
explicit baseline against which trace, reset, and rearrangement are
the deviations. A high anchor count means the brain mostly preserves
its module structure across rest → task → rest; a low anchor count
means task-induced reorganization dominates.

## Notation

(Inherits from `2026-05-07_kc-reset-modules.md`,
`2026-05-07_kc-rearrangement-modules.md`, and audit_47.)

Let `H_φ` be the LRG dendrogram of patient `p` in band `b` at phase
`φ ∈ {pre, tt, post}`. `S_φ` = set of internal subtrees with leaf
sizes in `[3, ⌊N/2⌋]`. Jaccard `J(A, B) = |A ∩ B| / |A ∪ B|`.

## Definitions

**Anchor module candidate.** A leaf set `L_pre` indexed by an
internal node `u_pre ∈ S_pre`. Anchored at rPre (any phase
would do — the predicate is symmetric across phases — but rPre is
chosen as the natural baseline anchor).

**Anchor predicate.** `Anchor(L_pre)` holds iff:

1. **Triple matching at the same Jaccard floor.**
   `∃ u_tt ∈ S_tt, u_post ∈ S_post :`
   `   J(L_pre, L_tt(u_tt)) ≥ J_ANCHOR ∧`
   `   J(L_pre, L_post(u_post)) ≥ J_ANCHOR ∧`
   `   J(L_tt(u_tt), L_post(u_post)) ≥ J_ANCHOR`
   where `J_ANCHOR = 0.6` (mirror of `J_PERSIST` in trace and reset).

2. **Size matching.** All three matched subtrees within
   `SIZE_RATIO_RANGE = 1.5` of each other in size — restricts the
   search to comparable-scale modules.

3. **Size floor.** `|L_pre| ≥ MIN_SIZE = 3`.

4. **Maximizer pick.** If multiple `(u_tt, u_post)` candidates
   satisfy gate 1, pick the pair that maximizes the geometric mean
   of the three Jaccards: `(J(L_pre, L_tt) · J(L_tt, L_post) · J(L_pre, L_post))^(1/3)`.

The anchor module is the triple `(L_pre, L_tt*, L_post*)`.

## Properties

- **Range.** Jaccards in `[0, 1]`; anchor count `n_anchor ∈
  {0, 1, 2, …}`. Per-(patient, band) output is the count plus the
  triple-leaf-set list.

- **Disjoint with trace.** Trace requires `J(L_tt, L_pre) < 0.5`
  (gate 2). Anchor requires `J(L_tt, L_pre) ≥ 0.6`. Disjoint by
  construction.

- **Disjoint with reset.** Reset requires `∀ tt: J(L_pre, tt) < 0.5`.
  Anchor requires `∃ tt: J(L_pre, tt) ≥ 0.6`. Disjoint by
  construction.

- **Disjoint with rearrangement.** Rearrangement requires `∀ pre:
  J(L_post, pre) < 0.5`. Anchor requires `∃ pre: J(L_post, pre) ≥
  0.6`. Disjoint by construction.

- **Phase-symmetry.** The anchor predicate is symmetric under any
  permutation of phases. Anchor at rPre, anchor at tt, and anchor at
  rPost should produce the same module set up to leaf-set ordering
  (since gate 1 enforces a triple match). Sanity-check at
  implementation: anchor lists should overlap heavily under different
  anchor phases.

- **What the measure DOES NOT detect.**
  - **Modules that drift gradually** — small Jaccard shifts < 0.4
    accumulating across phases would still pass at `J ≥ 0.6` per
    pair. A "graded anchor" measure would need a strict per-pair
    `J ≥ 0.85` — out of scope here, document as an open question.
  - **Height-only changes within a stable topology** — same blind
    spot as the other measures.

- **Complexity.** Triple-loop over rPre × tt × rPost subtrees,
  filtered by size. `O(N³)` worst case, but the size-ratio gate
  prunes aggressively. Tractable at N ≈ 115 (≈ 100 ³ at most).

- **Cohort gate.** No pre-registered threshold. Anchor counts
  serve as the comparison baseline. We expect anchor count ≫ trace
  count ≫ reset count ≈ rearrangement count under the trace
  hypothesis. Defer pre-registration until first sweep.

## Caveats & failure modes

| | Caveat | Mitigation |
|--|--|--|
| 1 | Probe bias / IMCOH_ABS — same as the other measures. | No additional mitigation; era-level FC choice handles it. |
| 2 | Edge effects at the dendrogram root. | Size cap `≤ ⌊N/2⌋` blocks root-near matches. |
| 3 | Pat_03 outlier (1024 Hz). | Log Pat_03 anchor count separately; flag in cohort plot. |
| 4 | Greedy de-duplication picks one anchor per overlap class. | Same `MAX_OVERLAP_FRAC = 0.30`. |
| 5 | Anchor counts may be very high — most modules are anchors at the IMCOH_ABS LRG layer. | Output is a count + size distribution. If the count exceeds `MAX_ANCHOR_TOTAL = 12` per (patient, band), report only the top 12 by triple-Jaccard score and log the total in the CSV. |

## Pseudocode

```
input:  H_pre, H_tt, H_post, N
output: anchor_modules

J_ANCHOR    = 0.60
MIN_SIZE    = 3
SIZE_CAP    = floor(N / 2)
SIZE_RATIO  = 1.5
MAX_OVERLAP = 0.30
MAX_TOTAL   = 12

S_pre  = subtrees(H_pre,  size_lo=MIN_SIZE, size_hi=SIZE_CAP)
S_tt   = subtrees(H_tt,   size_lo=MIN_SIZE, size_hi=SIZE_CAP)
S_post = subtrees(H_post, size_lo=MIN_SIZE, size_hi=SIZE_CAP)

candidates ← []
for u_pre in S_pre:
    L_pre ← leaves(u_pre)
    s     ← |L_pre|
    s_lo  ← max(MIN_SIZE, s / SIZE_RATIO)
    s_hi  ← min(SIZE_CAP, s * SIZE_RATIO)

    best_score ← 0
    best_tt    ← null
    best_post  ← null
    for u_tt in S_tt with size in [s_lo, s_hi]:
        j_pre_tt ← Jaccard(L_pre, leaves(u_tt))
        if j_pre_tt < J_ANCHOR: continue
        for u_post in S_post with size in [s_lo, s_hi]:
            j_pre_post ← Jaccard(L_pre, leaves(u_post))
            if j_pre_post < J_ANCHOR: continue
            j_tt_post ← Jaccard(leaves(u_tt), leaves(u_post))
            if j_tt_post < J_ANCHOR: continue
            score ← (j_pre_tt * j_pre_post * j_tt_post)^(1/3)
            if score > best_score:
                best_score ← score
                best_tt    ← u_tt
                best_post  ← u_post

    if best_tt is null: continue

    candidates ← append (L_pre, leaves(best_tt), leaves(best_post),
                         h_pre(u_pre), h_tt(best_tt), h_post(best_post),
                         best_score)

# Greedy de-duplication on union(L_pre, L_tt, L_post)
sort candidates descending by (size, score)
anchor_modules ← []
used ← ∅
for cand in candidates:
    union_leaves ← cand.L_pre ∪ cand.L_tt ∪ cand.L_post
    overlap ← |union_leaves ∩ used| / |union_leaves|
    if overlap < MAX_OVERLAP:
        anchor_modules ← append cand
        used ← used ∪ union_leaves
    if |anchor_modules| ≥ MAX_TOTAL: break

return anchor_modules
```

## Visualization spec

3 × 3 grid identical in shape to audit_47 / audit_48:

- **Rows:** λ ∈ {0, 0.5, 1} — KC blend regimes. λ=0 is all anchors
  (topology); λ=1 adds a triple-height-match constraint
  `max_pair |h_a − h_b| / max(h_a, h_b) ≤ HEIGHT_REL_THR = 0.30`;
  λ=0.5 = intersection.
- **Columns:** rest_pre, task, rest_post.
- **Each cell:** dendrogram + FC network.
- **Coloring rules.** Each phase's dendrogram colors its own
  matched subtree (`L_pre`, `L_tt`, `L_post` respectively). Internal
  link coloring uses the smallest-containing-module rule.
- **Network:** phase-INVARIANT module membership = `L_pre ∪ L_tt ∪
  L_post`. Edge weights per phase. Anchor visualization is the
  baseline: all three networks should look similar (bright
  intra-module edges across all three phases). Visual contrast is
  minimal — and that's the point.
- **Layout:** `lrg_sfdp` seeded with the rest_pre partition at
  k = LAYOUT_K = 20.
- **Merge-height bars** (λ ∈ {0.5, 1}): three dashed `axhline`s
  per module, one per phase at `h_pre`, `h_tt`, `h_post`. If the
  height-match constraint is met (λ=1), the three bars cluster
  together vertically; large spread → λ=1 fails.

**Reading rules.** An anchor module reads as: colored coherent
subtree on ALL THREE dendrograms; bright intra-module edges in ALL
THREE networks; little visual contrast across phases. Anchor =
"didn't move".

## Connection to prior tools

| Prior tool | Relation to KC anchor |
|--|--|
| KC trace modules (audit_47) | Anchors are the trace's complement at the cohort level. The trace count is what's left after removing anchors and other patterns. |
| KC reset modules (audit_48 scope) | Disjoint by gate construction. Anchor + reset cover the "rPre + rPost coherent" space partitioned by tt behavior (anchor: tt also coherent; reset: tt fragments). |
| KC rearrangement modules (audit_49 scope) | Disjoint by gate construction. |
| Terminology guide (`anchor`) | Direct match. The user-term "persist" is a synonym; the project's terminology guide standardizes "anchor". |
| Cohesion-CBR (2026-04-25) | Anchor modules should appear in Cohesion-CBR with high affinity in ALL THREE phases. Sanity-check overlap. |
| §5.2 KC scalar trace | Anchors do NOT contribute to `T_d^KC < 0` because all three pairwise distances are similar (high Jaccard everywhere). Anchors are the null. |

## Implementation plan

- **Script:** `scripts/01_compute/audit/audit_50_kc_anchor_module_view.py`
  (sister to audit_47–49).
- **Output figures:**
  `data/audit/kc_anchor_module_view/figures/{patient}__{band}__option_C.pdf`
- **Per-module summary CSV:**
  `data/audit/kc_anchor_module_view/anchor_subtrees_summary.csv`
- **Library reuse:** same as audit_48/49.
- **Refactor recommendation:** by audit_50, four scripts share the
  same scaffolding. Library promotion to
  `lrg_eegfc.visuals.kc_module_view` is overdue.
- **Empty-cell skip:** same convention.

## Open questions

1. **Higher Jaccard floor for anchors?** `J_ANCHOR = 0.6` matches the
   trace/reset persistence threshold but might over-count "weak
   anchors" that are really partial overlaps. A stricter
   `J_ANCHOR = 0.75` would isolate the strict-anchor population.
   Defer to first sweep; if anchor counts saturate at MAX_TOTAL for
   most cells, raise the floor.

2. **Triple-Jaccard scoring vs averaged.** Currently uses the
   geometric mean of the three pairwise Jaccards. An alternative is
   the minimum (`min(J_pre_tt, J_pre_post, J_tt_post)`). Geometric
   mean is gentler — flag for sensitivity check.

3. **Anchor at one phase vs symmetric enumeration.** Currently
   enumerates only rPre subtrees as candidates. Could enumerate over
   the union `S_pre ∪ S_tt ∪ S_post`; expected to give the same
   modules up to ordering (anchors are symmetric). Defer.

4. **Cohort gate / report format.** Anchor counts may need to be
   reported as ratios to total non-trivial subtrees, not absolute
   counts (since `|S_φ|` varies by patient). Will inform the
   manuscript figure.
