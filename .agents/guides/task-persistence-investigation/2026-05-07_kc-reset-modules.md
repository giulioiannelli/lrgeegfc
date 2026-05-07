---
name: kc-reset-modules
type: scope
era: COHORT_N10
status: draft
created: 2026-05-07
updated: 2026-05-07
pointers:
  - .agents/reports/2026-05-07_kc-trace-module-visualization.md
  - .agents/guides/01_project/terminology.md
  - .agents/guides/task-persistence-investigation/2026-04-25_task-trace-canonical.md
---

# KC reset modules — task-disrupted-then-restored subtrees

**Head.** Symmetric inverse of the KC trace-module visualization
(`audit_47`). A reset module is a leaf set that exists as a coherent
subtree in BOTH rest phases (rPre and rPost) but fragments during
the task. The cohort scalar trace `T_d^KC = d(tt, rPost) − d(rPre, tt) < 0`
pools two mechanisms: (i) emergence + persistence of new modules
(=trace, scope `2026-04-25_trace-modules.md` family) and
(ii) preservation of pre-existing modules across the task break
(=reset). This scope report defines the reset measure and pins the
visualization that surfaces these subtrees concretely on the FC
network and dendrogram.

## Notation

Let `H_φ` denote the LRG dendrogram of patient `p` in band `b` at
phase `φ ∈ {pre, tt, post}`, computed at `τ = 1/λ_max(L̂_φ)`. Each
internal node `u` of `H_φ` defines a subtree `S_φ(u)` with leaf set
`L_φ(u) ⊆ {1, …, N}` and merge height `h_φ(u) > 0` (linkage value
in `D_φ = 1/ρ̂_φ`). Write the set of *non-trivial* internal subtrees
as

    S_φ = { S_φ(u) : 3 ≤ |L_φ(u)| ≤ ⌊N/2⌋ }.

Jaccard similarity between two leaf sets is

    J(A, B) = |A ∩ B| / |A ∪ B|     ∈ [0, 1].

For a candidate leaf set `L`, the *smallest containing subtree* in
phase φ is

    C_φ(L) = arg min_{u : L ⊆ L_φ(u)} |L_φ(u)|       (size).

If no internal subtree contains `L` then `C_φ(L) = {1, …, N}` (the
root) and `|C_φ(L)| = N`.

## Definitions

**Reset module candidate.** A leaf set `L_pre` indexed by an
internal node `u_pre ∈ S_pre`. The module is *anchored* in rPre
(so the candidate enumeration ranges over rPre subtrees). The
test then asks: does this rPre module reappear in rPost (after
the task break) but fragment in tt?

**Reset predicate.** `Reset(L_pre)` holds iff:

1. **rPost match (anchor preservation).**
   `∃ u_post ∈ S_post : J(L_pre, L_post(u_post)) ≥ J_PERSIST`
   where `J_PERSIST = 0.65` (matched-pair locking; same threshold
   as KC trace). Let `L_post* = L_post(u_post*)` be the maximizing
   match.

2. **Task fragmentation, leaf-set side.**
   `∀ u_tt ∈ S_tt : J(L_pre, L_tt(u_tt)) < J_DISRUPT`
   with `J_DISRUPT = 0.5`. (Symmetric to the KC trace's
   anti-anchor gate but with tt and rPre swapped — here tt is the
   fragmenting phase.)

3. **Task fragmentation, containment side.**
   `|C_tt(L_pre)| ≥ FRAG_FACTOR · |L_pre|` with `FRAG_FACTOR = 2`.
   The smallest tt-subtree that contains `L_pre` must be at least
   twice the size of `L_pre`. Symmetric to the trace fragmentation
   guard: blocks the case where `L_pre` is fully contained in a
   slightly larger tt-subtree (Jaccard fails the dispersion test
   silently otherwise).

4. **Size floor.** `|L_pre| ≥ MIN_SIZE = 3`.

5. **Symmetric rPost-vs-tt gate.** `J_DISRUPT` is also enforced on
   `(L_post*, tt)` — i.e., the rPost match itself must fragment in
   tt. Equivalently, after locking `L_post*`, require
   `∀ u_tt ∈ S_tt : J(L_post*, L_tt(u_tt)) < J_DISRUPT`. Without
   this, a triple where rPre and rPost are both intact but tt only
   partially overlaps with one of them would pass on `(L_pre, tt)`
   alone while concealing that the rPost copy is still present in tt.

The **reset module** is the pair `(L_pre, L_post*)`; the algorithm
keeps both leaf sets so the visualization can render them
phase-specifically (matching the trace visualization's per-phase
coloring rule).

## Properties

- **Range.** `J ∈ [0, 1]`, `|C| ∈ [3, N]`. `Reset` is a Boolean
  predicate, so the per-(patient, band) output is a count
  `n_reset ∈ {0, 1, 2, …}` with the reset modules listed.

- **Identifiability.** Anchors at the *rPre* leaf set (not tt as in
  the trace): the question is "what existing rPre modules survive
  the task and reappear in rPost?", so the rPre subtree is the
  natural candidate-enumeration root.

- **Asymmetry with the trace.** The trace and reset are NOT
  partition-disjoint. A leaf set could in principle satisfy both
  predicates if its tt subtree both matches a rPost subtree (trace)
  and is the smallest containing tt-subtree of an rPre module
  (reset). In practice the gates J_PERSIST = 0.65 and J_DISRUPT =
  0.5 prevent overlap because a tt-subtree that matches rPost at
  ≥ 0.65 cannot simultaneously fragment rPre at < 0.5 *and* fragment
  rPre's content in tt at < 0.5 in the symmetric direction.
  Sanity-check during implementation by intersecting the two
  module lists per (patient, band).

- **What the measure DOES NOT detect.**
  - **Modules that exist in all three phases** (anchors, in the
    `terminology.md` taxonomy). They satisfy gate 1 but fail gate 2.
  - **Height-only resets** — modules with stable topology across
    pre / tt / post but a transient merge-height shift in tt. Same
    blind spot as the KC trace visualization for height-only
    traces (Pat_08-style cases).
  - **Modules that emerge during task and disappear in rPost**
    (transient). They have no rPre anchor — the candidate
    enumeration cannot find them.

- **Complexity.** Same as the trace measure. Enumerating
  `S_pre ∪ S_post ∪ S_tt` is `O(N)` per phase. The Jaccard scan
  per candidate is `O(|S_φ|)` — total `O(N²)` per (patient, band).
  Trivial at N ≈ 115.

- **Cohort gate.** Because reset is a different structural pattern
  than trace, "≥ 8/10" is not directly portable; reset cohort
  thresholds should be set after observing the per-patient counts
  rather than pre-registered against the trace cohort gate.

## Caveats & failure modes

| | Caveat | Mitigation |
|--|--|--|
| 1 | Probe bias — same-probe contacts have artificially high MSC at coarse `h_rel`. | Era is `IMCOH_ABS`, which kills the zero-phase-lag artefact by construction (Nolte 2004). No same-probe zeroing needed; remove only if a per-patient anatomy review shows reset modules dominated by single-probe contacts. |
| 2 | Edge effects at the dendrogram root. | The size cap `|L_φ(u)| ≤ ⌊N/2⌋` prevents the root and near-root subtrees from being matched. |
| 3 | Ties in heights across nearby internal nodes. | Use scipy linkage with `method='average'` consistently; ties are broken deterministically. |
| 4 | Pat_03 outlier (1024 Hz, 3× higher MSC). | IMCOH_ABS-era so MSC is not an issue, but log Pat_03 reset count as a secondary row in the output table; flag it distinctly in the cohort plot per the never-always-list. |
| 5 | Greedy de-duplication picks one module per overlap class — could miss a near-duplicate that is structurally distinct. | Same `MAX_OVERLAP_FRAC = 0.30` as the KC trace measure; flag in the open-questions section. |
| 6 | High `J_PERSIST` thresholds yield zero modules in some patient-band cells (sparsity). | Mirror the trace's `J_PERSIST = 0.65` (validated to give 25/27 visible cells in the trace gallery); if the reset gallery is too sparse, document and lower to 0.6 in a follow-up scope report. |

## Pseudocode

```
input:  H_pre, H_tt, H_post                # SciPy linkage matrices
        N                                   # number of leaves
output: reset_modules                      # list of (L_pre, L_post)

J_PERSIST   = 0.65
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
for u_pre in S_pre:
    L_pre  ← leaves(u_pre)
    s      ← |L_pre|

    # Gate 1: rPost match (anchor preservation)
    best_post  ← null
    best_jpost ← 0
    for u_post in S_post:
        if size_ratio(L_pre, leaves(u_post)) outside [1/1.5, 1.5]: continue
        j ← Jaccard(L_pre, leaves(u_post))
        if j > best_jpost: best_jpost ← j; best_post ← u_post
    if best_jpost < J_PERSIST: continue

    L_post ← leaves(best_post)

    # Gate 2: tt fragmentation, leaf-set side (no size restriction)
    fail ← false
    for u_tt in S_tt:
        if Jaccard(L_pre, leaves(u_tt)) ≥ J_DISRUPT: fail ← true; break
    if fail: continue

    # Gate 3: tt fragmentation, containment side
    smallest_containing_tt ← min { |leaves(u_tt)| : L_pre ⊆ leaves(u_tt) }
    if smallest_containing_tt < FRAG_FACTOR * s: continue

    # Gate 5: symmetric rPost-vs-tt gate
    fail2 ← false
    for u_tt in S_tt:
        if Jaccard(L_post, leaves(u_tt)) ≥ J_DISRUPT: fail2 ← true; break
    if fail2: continue

    candidates ← append (L_pre, L_post, h_pre(u_pre), h_post(best_post))

# Greedy de-duplication on L_pre (modules sharing > 30 % leaves collapse)
sort candidates descending by (s, best_jpost)
reset_modules ← []
for cand in candidates:
    if any(Jaccard(cand.L_pre, m.L_pre) > MAX_OVERLAP for m in reset_modules):
        continue
    reset_modules ← append cand
    if |reset_modules| ≥ MAX_TOTAL: break

return reset_modules
```

## Visualization spec

3 × 3 grid identical in shape to the KC trace visualization (audit_47):

- **Rows:** λ ∈ {0, 0.5, 1} — KC blend regimes. λ=0 is all reset
  modules (topology); λ=1 adds the height-match constraint
  `|h_pre − h_post| / max(h_pre, h_post) ≤ HEIGHT_REL_THR = 0.30`;
  λ=0.5 = intersection.
- **Columns:** rest_pre, task, rest_post.
- **Each cell:** dendrogram (left) + FC network (right).
- **Coloring rules** (mirror of trace, with phase semantics
  swapped):
  - rPre dendrogram: colors `L_pre` of every reset module as a
    coherent subtree (since rPre is the anchor).
  - rPost dendrogram: colors `L_post` (the rPost match leaves) as
    a coherent subtree.
  - tt dendrogram: colors `L_pre ∪ L_post` as **scattered
    leaf-bars**. The fragmentation visualization is the centerpiece
    of the figure: leaves that are coherent in pre and post burst
    apart in tt.
  - Internal links color by smallest module M with
    `link.leaves ⊆ M.{L_pre|L_post|union}` under that phase's
    leaf-set rule.
- **Network:** phase-INVARIANT module membership (`L_pre ∪ L_post`),
  edge rank-percentile per phase. Strong intra-module edges in
  pre/post saturate; weak tt intra-module edges fade.
- **Layout:** `lrg_sfdp` seeded with the rest_pre partition at
  k = LAYOUT_K = 20 (the rPre-anchored module enumeration motivates
  using the rPre partition rather than the tt partition that the
  trace visualization uses).
- **Merge-height bars** (λ ∈ {0.5, 1} only): `axhline` at `h_pre`
  in the rPre dendrogram and at `h_post` in the rPost dendrogram,
  in the matching color. tt has no anchor module so no axhline.

**Reading rules.** A reset module reads as: colored coherent
subtree on the rPre and rPost dendrograms; colored *scattered*
leaf bars on the tt dendrogram. In the network: bright intra-module
edges in rPre and rPost, faded edges in tt. Visually, reset =
"opens up during task, closes again after". Trace = "opens up
during task and stays open" — and the comparison should make
clear that resets and traces are not the same phenomenon despite
the cohort scalar `T_d^KC < 0` being driven by both.

## Connection to prior tools

| Prior tool | Relation to KC reset |
|--|--|
| KC scalar trace (Section 5.2 / 5.3) | Pools trace + reset + height components into one number. The reset measure decomposes one of those components into per-patient identifiable modules. |
| KC trace modules (audit_47) | Symmetric — both detect cross-phase Jaccard pattern but at different anchors (tt vs rPre) and gates (anti-rPre vs anti-tt). Reset is *not* a special case of trace, even at λ=0; the predicates are different. |
| MRL (Module Retention Landscape, 2026-04-25) | MRL counts tt-subtrees that are NEW in rPost (`¬present_pre ∧ present_post`). Reset is the orthogonal regime: `present_pre ∧ ¬present_tt ∧ present_post`. Together they cover the trace and reset corners of the `(p, t, r)` regime cube from `2026-04-25_task-trace-canonical.md`. |
| Cohesion-CBR (2026-04-25) | Continuous Jaccard + corner-distance affinity per module. Reset modules should be a subset of high-affinity Cohesion-CBR modules with the additional constraint of low task-affinity. Use Cohesion-CBR scores as a sanity check rather than a replacement. |
| §5.4 Grassmann principal angles (audit_46) | Subspace-level — sees rotations of eigenspaces. Reset is module-level — sees specific leaf subsets. Complementary; the cohort signal in §5.4 is consistent with reset (eigenspace rotation back to baseline). |

The reset measure does **not** replace the cohort scalar test — it
adds a per-patient mechanistic decomposition.

## Implementation plan

- **Script:** `scripts/01_compute/audit/audit_48_kc_reset_module_view.py`
  (sister to audit_47).
- **Output figures:**
  `data/audit/kc_reset_module_view/figures/{patient}__{band}__option_C.pdf`
- **Per-module summary CSV:**
  `data/audit/kc_reset_module_view/reset_subtrees_summary.csv`
- **Library reuse (no private re-implementations):**
  - `lrg_eegfc.utils.metrics.tree.tree_internal_nodes` — subtree
    enumeration
  - `lrg_eegfc.workflow.fc.load_fc_matrix`,
    `lrg_eegfc.workflow.lrg.load_lrg_result` — data IO
  - `lrg_eegfc.visuals.network_templates.compute_layout`,
    `_finalize_axes` — layout + network rendering (already used
    by audit_47)
  - `scipy.cluster.hierarchy.dendrogram, linkage` — dendrogram
- **Refactor opportunity:** if audit_47 and audit_48 share more
  than 60 % of their structure (likely), promote the shared
  scaffolding to `lrg_eegfc.visuals.kc_module_view` in the same
  PR as audit_48, per the project's "second-use library promotion"
  rule (`coding-rules.md`).
- **Cohort:** same n=9 subset that audit_47 uses, plus Pat_08 if
  any reset cells are non-empty (Pat_08 had 0 trace modules; reset
  may or may not produce visible modules — log either way).

## Open questions

1. **`J_PERSIST` value.** Mirror the trace at 0.65, or relax to 0.6
   for the reset gallery if 0.65 produces too few modules?
   Decide after the first sweep.
2. **Symmetric rPost-vs-tt gate (Gate 5).** Is this a strict
   requirement or should the algorithm count an asymmetric reset
   (rPre fragments in tt but rPost partly overlaps tt) as a
   half-credit module? Default: strict (Gate 5 enforced); revisit
   if the cohort is too sparse.
3. **Anchor at rPre vs rPost.** Could equivalently anchor at
   rPost subtrees and ask "does this rPost module exist in rPre
   but fragment in tt?" Should produce the same modules up to
   ordering. Default: anchor at rPre for prose consistency
   ("modules existed before the task and survived it").
4. **Cohort gate.** Cohort threshold for "trace exists" was ≥ 8/10.
   Reset is a different phenomenon and the literature offers no
   default. Defer pre-registration until after the first sweep.
5. **Height-only reset.** Out of scope for this measure (same blind
   spot as the trace visualization). A separate scope report
   covering merge-height-only patterns would be needed if the
   gallery is sparse.
