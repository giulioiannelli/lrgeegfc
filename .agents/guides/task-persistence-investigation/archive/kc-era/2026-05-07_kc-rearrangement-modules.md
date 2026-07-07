---
name: kc-rearrangement-modules
type: scope
era: COHORT_N10
status: draft
created: 2026-05-07
updated: 2026-05-07
pointers:
  - .agents/reports/2026-05-07_kc-trace-module-visualization.md
  - .agents/guides/task-persistence-investigation/2026-05-07_kc-reset-modules.md
  - .agents/guides/01_project/terminology.md
  - .agents/guides/task-persistence-investigation/2026-04-25_task-trace-canonical.md
---

# KC rearrangement modules — rPost-only emergent subtrees

**Head.** Counter-example to trace. A rearrangement module is a leaf
set that is coherent **only in rest_post** — fragmented in BOTH
rest_pre and task_test. Where trace asks "what did task induce that
persisted?" and reset asks "what did task disrupt that revertedstill?",
rearrangement asks "what reorganized in rest_post without precedent
in rest_pre OR task_test?". A nonzero rearrangement count is a
strict-emergent signal, distinct from the task-driven trace
mechanism. Surfacing rearrangement modules is the natural cohort
control for the trace claim: if the cohort produces large numbers
of rearrangement modules, the rPost coherence we attribute to task
might be partly rPost intrinsic dynamics, not task-locked.

## Notation

(Inherits from `2026-05-07_kc-reset-modules.md` and audit_47.)

Let `H_φ` denote the LRG dendrogram of patient `p` in band `b` at
phase `φ ∈ {pre, tt, post}`, computed at `τ = 1/λ_max(L̂_φ)`. Each
internal node `u` of `H_φ` defines a subtree `S_φ(u)` with leaf set
`L_φ(u) ⊆ {1, …, N}` and merge height `h_φ(u) > 0`. Write the set
of *non-trivial* internal subtrees as

    S_φ = { S_φ(u) : 3 ≤ |L_φ(u)| ≤ ⌊N/2⌋ }.

Jaccard similarity:  `J(A, B) = |A ∩ B| / |A ∪ B| ∈ [0, 1]`.

Smallest containing subtree in phase φ:
`C_φ(L) = arg min_{u : L ⊆ L_φ(u)} |L_φ(u)|`.

## Definitions

**Rearrangement module candidate.** A leaf set `L_post` indexed by
an internal node `u_post ∈ S_post`. Anchored at rPost (since rPost
is the coherent phase). The test asks: was this leaf set fragmented
in BOTH rPre AND tt before becoming a coherent rPost subtree?

**Rearrangement predicate.** `Rearrange(L_post)` holds iff:

1. **rPre fragmentation, leaf-set side.**
   `∀ u_pre ∈ S_pre : J(L_post, L_pre(u_pre)) < J_DISRUPT = 0.5`.

2. **rPre fragmentation, containment side.**
   `|C_pre(L_post)| ≥ FRAG_FACTOR · |L_post|` with `FRAG_FACTOR = 2`.
   The smallest rPre subtree containing `L_post` must be at least
   twice the size of `L_post`. Mirrors the trace and reset
   fragmentation guards.

3. **tt fragmentation, leaf-set side.**
   `∀ u_tt ∈ S_tt : J(L_post, L_tt(u_tt)) < J_DISRUPT = 0.5`.

4. **tt fragmentation, containment side.**
   `|C_tt(L_post)| ≥ FRAG_FACTOR · |L_post|`.

5. **Size floor.** `|L_post| ≥ MIN_SIZE = 3`.

The rearrangement module is just `L_post` (no anchor pair, since
neither rPre nor tt provides a matching subtree).

## Properties

- **Range.** `J ∈ [0, 1]`, `|C| ∈ [3, N]`. `Rearrange` is a Boolean
  predicate; per-(patient, band) output is a count
  `n_rearrange ∈ {0, 1, 2, …}` with the rearrangement modules listed.

- **Anchored at rPost.** Unique among the three measures (trace
  anchored at tt, reset anchored at rPre): rearrangement is
  rPost-anchored because rPost is the only phase where the module
  is coherent.

- **No matching pair.** Trace and reset both produce a `(L_anchor,
  L_match)` pair (tt↔rPost for trace, rPre↔rPost for reset);
  rearrangement has no match because both rPre and tt fragment
  the leaf set. Output is just `L_post`.

- **Asymmetry with trace.** A trace module's leaf set `L_tt`
  fragments in rPre but matches in rPost (J ≥ 0.6). A rearrangement
  module's leaf set `L_post` fragments in BOTH rPre and tt. By
  construction the two predicates are disjoint: `Reset ∩ Rearrange
  = ∅`, `Trace ∩ Rearrange = ∅`. (A rearrangement module is by
  definition not also a trace, because trace requires tt-coherence
  via Jaccard ≥ 0.6.)

- **Asymmetry with reset.** A reset module `L_pre` is coherent in
  rPre AND rPost, but fragmented in tt. A rearrangement module
  `L_post` fragments in both rPre AND tt. The predicates are also
  disjoint by gate 1 (rearrangement requires rPre fragmentation,
  reset requires rPre coherence).

- **What the measure DOES NOT detect.**
  - **Trace modules** — they have task-coherence; gate 3 rules them
    out.
  - **Reset modules** — they have rPre-coherence; gate 1 rules them
    out.
  - **Anchors** (modules unchanged across all phases) — they have
    both rPre and tt coherence; gates 1 and 3 both fail.
  - **Trace-by-merge-height-only** — same blind spot as audit_47:
    the visualization does not see height-only patterns where
    leaves stay grouped but merge heights shift.

- **Complexity.** Same as trace and reset: enumerating subtrees is
  `O(N)` per phase; per-candidate Jaccard scan is `O(|S_φ|)`.
  Total `O(N²)` per (patient, band). Trivial at N ≈ 115.

- **Cohort gate.** No pre-registered threshold. The expected value
  under the null hypothesis ("rPost coherence is task-driven, not
  intrinsic") is small or zero rearrangement counts; large
  rearrangement counts in any band would weaken the trace
  interpretation. Decide pass/fail empirically after first sweep.

## Caveats & failure modes

| | Caveat | Mitigation |
|--|--|--|
| 1 | A subtree very common in all 3 phases (anchor-like) might pass if both rPre and tt happen to fail Jaccard while rPost succeeds. | The FRAG_FACTOR = 2 containment guard rules these out — true anchors have a small containing subtree in rPre and tt, not 2× larger. |
| 2 | Probe bias / IMCOH_ABS era — same as trace. | No additional mitigation needed beyond the era-level FC choice. |
| 3 | Edge effects at the dendrogram root. | Size cap `|L_post| ≤ ⌊N/2⌋` prevents root match. |
| 4 | Pat_03 outlier (1024 Hz, 3× higher MSC). | Log Pat_03 rearrangement count separately; flag in cohort plot. |
| 5 | Random-baseline expectation may be nonzero — even noise dendrograms can produce tt fragmentation by chance. | Optional null calibration: run on within-baseline split halves of rPre vs full rPre/tt/post; the rPost-only-coherence pattern should be rare under null. Defer to first sweep. |
| 6 | Greedy de-duplication picks one module per overlap class. | Same `MAX_OVERLAP_FRAC = 0.30` as trace and reset; flag in open questions. |

## Pseudocode

```
input:  H_pre, H_tt, H_post
        N
output: rearrange_modules                  # list of L_post

J_DISRUPT   = 0.50
FRAG_FACTOR = 2.0
MIN_SIZE    = 3
SIZE_CAP    = floor(N / 2)
MAX_OVERLAP = 0.30
MAX_TOTAL   = 12

S_pre  = subtrees(H_pre,  size_lo=MIN_SIZE, size_hi=SIZE_CAP)
S_tt   = subtrees(H_tt,   size_lo=MIN_SIZE, size_hi=SIZE_CAP)
S_post = subtrees(H_post, size_lo=MIN_SIZE, size_hi=SIZE_CAP)

candidates ← []
for u_post in S_post:
    L_post ← leaves(u_post)
    s      ← |L_post|

    # Gate 1+2: rPre fragmentation
    fail ← false
    for u_pre in S_pre:
        if Jaccard(L_post, leaves(u_pre)) ≥ J_DISRUPT: fail ← true; break
    if fail: continue

    smallest_containing_pre ← min { |leaves(u_pre)| : L_post ⊆ leaves(u_pre) }
    if smallest_containing_pre < FRAG_FACTOR * s: continue

    # Gate 3+4: tt fragmentation
    fail2 ← false
    for u_tt in S_tt:
        if Jaccard(L_post, leaves(u_tt)) ≥ J_DISRUPT: fail2 ← true; break
    if fail2: continue

    smallest_containing_tt ← min { |leaves(u_tt)| : L_post ⊆ leaves(u_tt) }
    if smallest_containing_tt < FRAG_FACTOR * s: continue

    candidates ← append (L_post, h_post(u_post))

# Greedy de-duplication
sort candidates descending by s
rearrange_modules ← []
for cand in candidates:
    if any(Jaccard(cand.L_post, m.L_post) > MAX_OVERLAP for m in rearrange_modules):
        continue
    rearrange_modules ← append cand
    if |rearrange_modules| ≥ MAX_TOTAL: break

return rearrange_modules
```

## Visualization spec

3 × 1 grid (single row, three columns) — KC λ regimes are not
applicable here because rearrangement has no matched pair to compare
heights between. λ-blend distinguishes topology vs heights of a
matched pair `(L_a, L_b)`; a rearrangement module has only `L_post`.

Layout per (patient, band):

- Columns: rest_pre, task, rest_post.
- Each cell: dendrogram (left) + FC network (right).
- **Coloring rules** (mirror of trace, with phase semantics shifted):
  - rPost dendrogram: colors `L_post` of each rearrangement module
    as a coherent subtree.
  - tt dendrogram: colors `L_post` as **scattered leaf-bars** —
    fragmentation visualization.
  - rPre dendrogram: colors `L_post` as scattered leaf-bars —
    second fragmentation visualization.
  - Internal links: smallest-containing-module rule, same as trace.
  - Explicit colored vertical bars at each leaf x-position
    (scipy `link_color_func` workaround, same as trace).
- **Network:** phase-INVARIANT module membership (`L_post`), edge
  rank-percentile per phase. Strong intra-module edges in rPost
  saturate; weak intra-module edges in rPre AND tt fade. The visual
  contrast distinguishes rearrangement from trace: trace networks
  show task as the *strongest* phase, while rearrangement networks
  show rPost as the *only* strong phase.
- **Layout:** `lrg_sfdp` seeded with the rest_post partition at
  k = LAYOUT_K = 20. (rPost is the anchor.)
- **Merge-height bar.** Single dashed `axhline` at `h_post(u_post)`
  in the rPost dendrogram only, in the matching color. Trace and
  reset both have two heights to compare; rearrangement has one.

**Reading rules.** A rearrangement module reads as: colored coherent
subtree on the rPost dendrogram; colored scattered leaf-bars on
BOTH the rPre and tt dendrograms. In the network: bright intra-module
edges in rPost; faded edges in rPre AND tt. Visually, rearrangement
= "didn't exist before, didn't exist during task, exists only after."

**Reading rules in cohort context.** A small rearrangement count
across the cohort is consistent with the trace interpretation
("rPost coherence is task-locked"). A large rearrangement count in
any band weakens the claim — it shows rPost reorganization can
happen without task involvement. Treat the rearrangement gallery as
a control on the trace gallery, not as primary evidence.

## Connection to prior tools

| Prior tool | Relation to KC rearrangement |
|--|--|
| KC trace modules (audit_47) | Disjoint by construction (rearrangement requires tt fragmentation; trace requires tt coherence). Rearrangement is the symmetric phase-coherence inverse: rPre + tt fragment vs rPre fragment alone. |
| KC reset modules (2026-05-07 scope) | Disjoint by construction (rearrangement requires rPre fragmentation; reset requires rPre coherence). Reset and rearrangement are exhaustive partitions of the "rPost-coherent + something fragmented" space. |
| MRL (Module Retention Landscape, 2026-04-25) | MRL counted tt-subtrees that are NEW in rPost (`¬present_pre ∧ present_post`). Rearrangement is the rPost-anchored variant: stronger constraint that the leaf set also fragments in tt. |
| Cohesion-CBR (2026-04-25) | Continuous Jaccard + corner-distance; rearrangement modules should appear in Cohesion-CBR with high rPost-anchored affinity AND low rPre + low task affinity simultaneously. Use Cohesion-CBR scores as a sanity check. |
| §5.2 KC scalar trace | Pools trace + reset components (both produce `T_d^KC < 0`). Rearrangement does NOT contribute to `T_d^KC < 0` because the tt-rPost distance is high (fragmented). Rearrangement appears in `d(rPre, tt)` ≈ `d(rPre, rPost)` ≈ `d(tt, rPost)` triangle, all comparable. |
| Terminology guide (`emergent`) | "emergent" is defined as "module that did not exist in RPre". Rearrangement is `emergent ∧ ¬tt_coherent`, a strict subcase. Document the distinction in `terminology.md` if any rearrangement modules survive the first sweep. |

## Implementation plan

- **Script:** `scripts/01_compute/audit/audit_49_kc_rearrangement_module_view.py`
  (sister to audit_47 + audit_48).
- **Output figures:**
  `data/audit/kc_rearrangement_module_view/figures/{patient}__{band}__option_C.pdf`
- **Per-module summary CSV:**
  `data/audit/kc_rearrangement_module_view/rearrangement_subtrees_summary.csv`
- **Library reuse:**
  - `lrg_eegfc.utils.metrics.tree.tree_internal_nodes`
  - `lrg_eegfc.workflow.fc.load_fc_matrix`,
    `lrg_eegfc.workflow.lrg.load_lrg_result`
  - `lrg_eegfc.visuals.network_templates.compute_layout`,
    `_finalize_axes`
  - `scipy.cluster.hierarchy.dendrogram, linkage`
- **Refactor recommendation.** Same as for audit_48: promote shared
  scaffolding (subtree enumeration, palette, dendrogram + network
  rendering, multi-cell composition) to
  `lrg_eegfc.visuals.kc_module_view`. With audit_47 + audit_48 +
  audit_49 sharing ≥ 60 % of structure, this is the second-use
  promotion threshold.
- **Empty-cell skip.** Same convention as audit_47: if the module
  list is empty, do not render a figure; preserve the CSV row with
  `option_C = ""`.
- **Cohort:** same n=9 subset that audit_47 sweeps (Pat_06, Pat_03,
  Pat_05, Pat_13, Pat_02, Pat_10, Pat_14, Pat_07, Pat_15) plus
  Pat_08 if any rearrangement cells are non-empty.

## Open questions

1. **Should rearrangement and reset share a script?** They both
   produce a single anchored leaf set with a 3 × 1 visualization
   layout, distinct only by which phase is the anchor. A combined
   script `audit_48_kc_phase_anchored_view.py` with a `pattern`
   argument might be cleaner than two parallel scripts. Defer to
   implementation time after the library promotion lands.

2. **Symmetric J_DISRUPT threshold.** Currently 0.5 for both rPre
   and tt fragmentation. A more permissive 0.55 might surface more
   rearrangement modules without compromising the visual. Defer
   pending first sweep.

3. **Cohort gate / null model.** Should we pre-register a maximum
   acceptable rearrangement count per band as a cohort sanity
   check (e.g., "if any band has > 5 rearrangement modules per
   patient on average, the trace claim is weakened")? Defer until
   we have first-pass numbers.

4. **Visualization at λ regimes.** Single row vs three rows (with
   λ blend annotations on heights). Default: single row, since
   λ-blend has no meaningful target without a matched pair. If a
   reviewer asks for "consistency with trace + reset visualizations",
   add a 3 × 3 layout that simply repeats the same single-row data
   three times — purely cosmetic.

5. **Anchor at rPost vs combined.** Could equivalently look for
   leaf sets whose smallest containing subtree in rPre AND tt is
   ≥ 2 × the leaf set's size, by enumerating *all* subsets — but
   that's combinatorially infeasible. Anchoring at rPost subtrees
   is the natural restriction.
