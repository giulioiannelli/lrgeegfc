---
name: multiscale-partition-coherence
type: scope
era: COHORT_N10
status: draft
created: 2026-04-26
updated: 2026-04-27
pointers:
  - .agents/guides/task-persistence-investigation/README.md
  - .agents/guides/task-persistence-investigation/2026-04-25_cophenetic-neighbourhood.md
  - .agents/guides/task-persistence-investigation/2026-04-25_cohesion-cbr.md
  - .agents/reports/2026-04-24_multiscale-task-trace.md
  - src/lrg_eegfc/utils/metrics/tree.py
---

# Multi-Scale Partition Coherence (MSPC)

**Per-leaf, per-scale Jaccard of cluster-mate sets across phases. At every
dendrogram cut `k ∈ [2, K_max]`, each leaf has a set of cluster-mates in
each phase; MSPC compares those sets across phases at every scale and
classifies each (leaf, scale) cell into trace / persist / reset /
rearrange via corner-distance affinities. Directly exploits the
multiscale character of the LRG ultrametric — does not require spectral
gaps or discrete modules. Two visualisations: (1) the standard 4-phase
dendrogram layout with per-leaf k-averaged dominant pattern;
(2) a multiscale strip per (patient, band) with leaves × scales coloured
by per-cell dominant pattern, so a reader can see how each leaf's
classification changes with scale.**

## 1. Notation

Inherits the cohort/band/phase notation of the parent scope reports.
Recap and additions:

- `P, B, Φ` — `COHORT_N10` patients, 6 bands, 4 phases.
- `Z^{p,b,φ}` — scipy linkage on `imcoh_abs`; `T^{p,b,φ}` the induced
  tree with `N_p` leaves.
- For each scale `k ∈ {2, …, K_max(p)}` and phase `φ`,

      c^{p,b,φ}_k : L_p → {1, …, k}

  is the flat partition produced by `scipy.cluster.hierarchy.fcluster(
  Z^{p,b,φ}, k, criterion='maxclust')`. Each leaf gets a cluster id.

- The **cluster-mate set** of leaf `ℓ` at scale `k` in phase `φ` is

      coclust^{p,b,φ}_k(ℓ) := { j ∈ L_p \ {ℓ} : c^{p,b,φ}_k(j) = c^{p,b,φ}_k(ℓ) }

  i.e. the leaves in the same flat cluster as `ℓ`, excluding `ℓ` itself.
  `|coclust^{p,b,φ}_k(ℓ)| = (size of ℓ's cluster) − 1`. Singleton
  clusters give the empty set.

- **Scale window**: default `K_max(p) := ⌊N_p / 2⌋`. Caps at
  half-tree to avoid trivial fine-scale partitions where most leaves
  are singletons.

- **Task collapse** at each (k, ℓ):

      J_task_k(ℓ) := mean( J(coclust^{task_learn}_k(ℓ), coclust^{rest_post}_k(ℓ)),
                              J(coclust^{task_test}_k(ℓ),  coclust^{rest_post}_k(ℓ)) )

  If `task_test` is missing for this patient (legacy Pat_14 case
  pre-2026-04-25 vendor replacement; n/a today), use task_learn alone.

## 2. Definitions

### 2.1 Per-(leaf, scale) Jaccard of cluster-mates

For phases a, b ∈ Φ and (leaf, scale) (ℓ, k):

    J^{a,b}_k(ℓ) := |coclust^a_k(ℓ) ∩ coclust^b_k(ℓ)| / |coclust^a_k(ℓ) ∪ coclust^b_k(ℓ)|

with `J = 0` by convention when both sets are empty.

`J^{a,b}_k(ℓ) ∈ [0, 1]`. `J = 1` iff `ℓ` has the exact same cluster-mates
in phases a and b at scale k. `J = 0` iff no overlap.

The two load-bearing scores per (ℓ, k):

    J_pp(k, ℓ) := J^{rest_pre, rest_post}_k(ℓ)
    J_tp(k, ℓ) := mean(J^{task_learn, rest_post}_k(ℓ), J^{task_test, rest_post}_k(ℓ))

### 2.2 Per-(leaf, scale) classification

For each (ℓ, k) compute corner-distance affinities on the unit square
`[0, 1]²` (Jaccard ≥ 0, unlike Spearman):

    TRACE corner     := (0, 1)        low J_pp, high J_tp
    PERSIST corner   := (1, 1)        both high
    RESET corner     := (1, 0)        high J_pp, low J_tp
    REARRANGE corner := (0, 0)        both low

    a_p(ℓ, k) := 1 − ‖(J_pp(k, ℓ), J_tp(k, ℓ)) − corner_p‖₂ / √2          ∈ [0, 1]

(`√2` is the diameter of the unit square, so a_p reaches 0 only at the
opposite corner.)

    dominant(ℓ, k)   := argmax_p a_p(ℓ, k)
    confidence(ℓ, k) := a_dominant(ℓ, k) − a_secondbest(ℓ, k)            ∈ [0, 1]

### 2.3 K-averaged per-leaf classification

For the standard 4-phase layout we collapse across scales:

    Ā_p(ℓ) := mean_k a_p(ℓ, k)                                            ∈ [0, 1]
    dominant(ℓ)   := argmax_p Ā_p(ℓ)
    confidence(ℓ) := Ā_dominant(ℓ) − Ā_secondbest(ℓ)

This gives the "what is this leaf's dominant story across scales"
classification used for the dendrogram colour.

### 2.4 Singletons, empty sets

If at scale k both `|coclust^a_k(ℓ)| = 0` and `|coclust^b_k(ℓ)| = 0`
(leaf ℓ is a singleton cluster in both phases at this k), define
`J = 0`. The leaf is "alone everywhere" at this scale — corresponds to
the rearrange / unclassified end of the spectrum. Filtered into
rearrange via the corner-distance affinity at `(0, 0)`.

## 3. Properties

| property | value |
|:---|:---|
| Range of J | `[0, 1]` |
| Range of a_p | `[0, 1]` |
| Multiscale | **Yes by construction** — every k from 2 to ⌊N/2⌋ contributes. |
| Anchor-free | Yes — no rsPost subtree scan, no SIZE_MIN floor. |
| Robust to absence of spectral gap | Yes — works on continuous hierarchies; we just cut at every scale and look at cluster-mate sets. |
| Robust to rake topology | Yes — Pat_10 δ rake just gives rapidly shrinking cluster sizes as k grows; the leaf still has a well-defined cluster-mate set at each k. |
| Sensitive to fragmentation | Yes — fragmentation in task = ℓ's task cluster-mates at scale k differ from rsPost cluster-mates = `J_tp` low at that k. |
| Detects band-specific scale signatures | Yes — different bands can show different patterns of dominant-vs-k via the multiscale strip. |
| Complexity | `O(K · N²)` per (patient, band) for the all-pair cluster-mate matrices. For N=120, K=60: ~860k ops × 4 phases = 3.4M per cell. ~200M cohort-wide; trivial. |

**Negative property — what MSPC cannot do:**
- Doesn't name a *specific* subtree that persists. MSPC says "this leaf
  ℓ has trace-like cluster-mates at scales k=10–25 and persist-like at
  k=2–8"; it does not extract the leaf-set of a single coherent
  subtree the way CBR does.
- Doesn't separate functional from anatomical co-clustering. If two
  leaves are always-near-each-other for non-task reasons (assuming the
  ImCoh-vs-MSC reduction of probe-bias is real, which we have not
  verified in this session), MSPC will still mark them persist.
  Probe-debiased FC pipeline would address it.
- The dominant-vs-k pattern can be noisy at fine scales where most
  clusters are singletons. We expect saturation around mid-k.

## 4. Caveats & failure modes

| caveat | mitigation |
|:---|:---|
| K-averaged dominant blurs scale-specific signal. | Always emit the multiscale strip alongside; never report Ā_p as the only result. |
| At very small k (k=2,3), the partition is dominated by anatomy / global structure; J's tend to be high everywhere → persist is over-represented at coarse scales. | Acknowledged. The K-averaged collapse weights all k equally; if the user wants to emphasise mid-scale they can apply a k-window. Open question §9. |
| Task collapse via mean blends timescales. | Inherits the same issue from CNP / cohesion-CBR. Defer mean-vs-min decision to first visual review. |
| Singleton-singleton case (J = 0 at fine k for sparse clusters) inflates rearrange. | Solid rearrange at fine k is honest if ℓ is alone; filter at projection time if needed. |

## 5. Pseudocode

```
INPUT  Z_rpre, Z_tl, Z_tt, Z_rpost  (linkage matrices, same N leaves)
       K_max = floor(N / 2)

dominant_at_k[ℓ, k]  ← 0          # for the strip
A_table[ℓ, p]        ← 0          # for the K-averaged dominant
n_k_seen[ℓ]          ← 0

FOR k in {2, …, K_max}:
    labels_rpre  ← fcluster(Z_rpre,  k, criterion='maxclust')
    labels_tl    ← fcluster(Z_tl,    k, criterion='maxclust')
    labels_tt    ← fcluster(Z_tt,    k, criterion='maxclust')   IF Z_tt
    labels_rpost ← fcluster(Z_rpost, k, criterion='maxclust')

    same_pre   ← labels_rpre[:, None] == labels_rpre[None, :]
    same_tl    ← labels_tl[:, None]   == labels_tl[None, :]
    same_tt    ← labels_tt[:, None]   == labels_tt[None, :]      IF Z_tt
    same_rpost ← labels_rpost[:, None] == labels_rpost[None, :]

    # Jaccard per leaf ℓ between coclust sets
    J_pp_k[ℓ] ← |same_pre[ℓ] & same_rpost[ℓ]| / |same_pre[ℓ] | same_rpost[ℓ]|
              evaluated row-wise (subtract 1 for self-match before dividing)
    J_tl_k[ℓ] ← same way with same_tl, same_rpost
    IF Z_tt: J_tt_k[ℓ] ← same way with same_tt, same_rpost
    J_tp_k[ℓ] ← (J_tl_k[ℓ] + J_tt_k[ℓ]) / 2  IF Z_tt ELSE J_tl_k[ℓ]

    FOR each leaf ℓ:
        a_trace     ← affinity(J_pp_k[ℓ], J_tp_k[ℓ], corner=(0, 1))
        a_persist   ← affinity(...,                          corner=(1, 1))
        a_reset     ← affinity(...,                          corner=(1, 0))
        a_rearrange ← affinity(...,                          corner=(0, 0))
        dominant_at_k[ℓ, k] ← argmax_p a_p
        A_table[ℓ, p]       += a_p
        n_k_seen[ℓ]         += 1

A_table[ℓ, p] /= n_k_seen[ℓ]
dominant[ℓ]   ← argmax_p A_table[ℓ, p]
confidence[ℓ] ← top1(A_table[ℓ]) − top2(A_table[ℓ])

OUTPUT dominant_at_k, A_table, dominant, confidence
```

`affinity(x, y, corner)` returns `1 - sqrt((x-corner[0])² + (y-corner[1])²) / sqrt(2)`.

## 6. Visualization spec

Two figures per (patient, band):

### 6.1 Standard layout (audit_08 compatible)

`data/audit/per_patient_hierarchy_mspc/{Pat_XX}/{band}_hierarchy-crossphase.pdf`

Identical to audit_08 / cohesion-CBR / CNP — 4 phase dendrograms with
codebar. Per-leaf colour from §2.3 dominant; solid alpha (`1.0`).
Branches: dominant colour iff every leaf below shares the same
dominant; else background gray.

### 6.2 Multiscale strip

`data/audit/per_patient_hierarchy_mspc/{Pat_XX}/{band}_multiscale-strip.pdf`

A single `N × K_max` image:
- y-axis: leaves, ordered by rsPost dendrogram display order (so
  leaves clustering nicely in rsPost are adjacent).
- x-axis: scale k from 2 to K_max.
- Colour at cell (ℓ, k): `PATTERN_COLOR[dominant(ℓ, k)]`.
- Solid alpha (no confidence-fading at this resolution; the reader
  should pick out coloured "vertical bands" for scale-coherent
  patterns and "horizontal stripes" for leaves with consistent
  cross-scale stories).

A right-margin codebar duplicates the K-averaged dominant per leaf
so the two figures' classifications are visually keyed.

**Reading rule (multiscale strip):**

- Saturated red horizontal stripe for a leaf = trace at most scales.
- Saturated black for a leaf = persist at most scales (anatomy-like).
- Vertical red column = many leaves are trace at the same specific k
  (a localised modular task imprint).
- Diagonal red transition (trace at coarse k, rearrange at fine k) =
  scale-specific reorganisation.
- Pat_10 δ would predictively show mostly persist + rearrange;
  Pat_06 γ_l would show clear vertical or horizontal red bands.

Auxiliary CSVs:
- `leaf_assignment.csv`: `patient, band, leaf_id, dominant, confidence,
  A_trace, A_persist, A_reset, A_rearrange`.
- `multiscale_assignment.csv`: `patient, band, leaf_id, k, dominant_at_k,
  J_pp, J_tp` (one row per (leaf, scale)).

## 7. Connection to prior tools

| existing tool | how MSPC relates |
|:---|:---|
| H2d (block coactivation persistence) | H2d aggregates over leaf pairs and over patients. MSPC is the **per-leaf-per-scale** version of H2d — same primitive (cluster co-membership at scale k), different aggregation. They will agree on direction; MSPC localises *where* and *at what scale* H2d's signal lives. |
| CBR family | CBR picks one anchor subtree per (patient, band, anchor_row). MSPC sweeps every scale of every leaf. CBR names "this 11-leaf module"; MSPC says "this leaf's cluster-mates at scale k=10 are task-like, at k=20 are persist-like". |
| CNP | CNP collapses across scales via Spearman on full cophenetic vectors. MSPC keeps the scales explicit. They share the per-leaf level but differ on scale handling. CNP is faster; MSPC carries more signal. |
| H2c | H2c is global per pair, k-averaged, cohort-pooled. MSPC is per-leaf, scale-resolved, per-cell. Same statistical family, opposite direction of localisation. |
| MRL | MRL operates on subtree leaf-sets at native heights. MSPC operates on flat-partition cluster-mate sets at integer k. Disjoint. |

## 8. Implementation plan

- **Library reuse**: `scipy.cluster.hierarchy.fcluster`. No new helpers
  needed.
- **Script**: `scripts/01_compute/audit/audit_14_mspc.py`. Reuses
  `plot_dendrogram_with_strip` from audit_08 for figure-1 and writes a
  custom imshow for figure-2.
- **Outputs**: as listed in §6.
- **No promotion to library** in this pass; if a third caller appears
  the cluster-mate Jaccard goes to `lrg_eegfc.utils.metrics.tree`.

## 9. Open questions

- **K-window weighting.** Currently uniform mean over all k from 2 to
  ⌊N/2⌋. A mid-scale weighting (e.g., gaussian peak at k = log₂N)
  would emphasise the scales where the trace lives in the gallery
  exemplars. Defer until first visual review.
- **Task collapse rule** (mean vs min vs separate). Default mean. Defer.
- **Confidence-graded alpha on the standard figure.** Currently solid
  to match user feedback on previous variants. The strip is solid at
  cell level; confidence could fade individual cells if needed.
- **K_max choice.** ⌊N/2⌋ vs N−1. The half-tree cap excludes very
  fine partitions where most clusters are singletons; if those scales
  carry interesting signal we'd want to extend. Defer.

## 10. MSPC vs CBR — granularity dichotomy

CBR (`2026-04-25_cbr-investigation.md`) and MSPC ask the same
scientific question — *did some structure appear during task and
persist into rsPost?* — but operate on different geometric objects.
CBR thinks in **subtrees**; MSPC thinks in **flat partitions**. The
distinction is sharp enough that the two measures can disagree on
the same data without contradicting each other.

### 10.1 Atomic objects

| | **CBR family** | **MSPC** |
|---|---|---|
| atomic object | a subtree (= leafset under one internal node of `T_φ`) | a single leaf `ℓ` |
| atomic question | "is leafset `L` coherent in phase φ?" | "what are leaf `ℓ`'s cluster-mates in phase φ at scale `k`?" |
| meaning of "coherent in φ" | `L` IS the leafset under some internal node of `T_φ`, possibly with slack (Jaccard ≥ τ, containment ≥ τ, or exact equality) | trivially coherent — `coclust_k(ℓ, φ)` is *defined* by the partition at every `k`; there is no coherence test |
| anchor | yes — pick rsPost or taskT internal nodes; iterate | none — every leaf participates |
| scale | implicit — each candidate subtree carries its native height | explicit — sweep `k = 2 … ⌊N/2⌋`; classify every `(ℓ, k)` cell |
| similarity | leafset-to-leafset Jaccard / containment / exact equality | row-wise Jaccard of cluster-mate sets at fixed `k` |
| threshold | hard cutoffs on similarity (legacy `J ≥ 0.5`) | soft 4-corner affinity, argmax (no cutoff) |
| output unit | a list of *modules* (leafsets) labelled trace/persist/reset/rearrange | a label per *leaf* (K-averaged) and per *(leaf, scale)* (strip) |
| detects | task-induced *modules* — explicit "this 11-leaf group appeared in task and persisted" | task-induced *neighbourhood-identity flips* — "this leaf has the same neighbours in task and post, but not in pre" |
| misses | distributed reorganisation (no single subtree captures it); partial overlaps below `τ`; modules that exist at integer `k` but never coalesce into a native subtree at any height | the *identity* of which leaves form a coherent group — MSPC says "ℓ₁ traces, ℓ₂ traces" without saying whether they are co-clustered |

### 10.2 Worked example

Suppose at `k = 10` the partitions are (only the involved clusters shown):

| phase | clusters |
|-------|----------|
| pre   | `{a,b,c}`, `{d,e,f,g}`, `{h,i,j,k,l,m,n,o,p,q}`, … |
| task  | `{a,b,c}`, `{d,h,i,j}`, `{e,f,g,k,l,m,n,o,p,q}`, … |
| post  | `{a,b,c}`, `{d,h,i,j}`, `{e,f,g,k,l,m,n,o,p,q}`, … |

Two events: cluster `{d,e,f,g}` split — `{d}` joined `{h,i,j}`, while
`{e,f,g}` merged into the larger cluster.

- **CBR (legacy on rsPost anchor):** scans rsPost subtrees. The subtree
  `{d,h,i,j}` exists in both task and post (`J = 1.0` against itself)
  and *does not* exist in pre (best Jaccard against any pre subtree is
  `0.25`, against `{d,e,f,g}`). `{d,h,i,j}` is labelled **trace**. CBR
  names the module: "this 4-leaf group is the task-induced module."

- **MSPC at `k = 10`:** evaluates every leaf —
  - `d`: pre-mates `{e,f,g}`, task-mates `{h,i,j}`, post-mates `{h,i,j}`.
    `J_pp = 0`, `J_tp = 1` → **trace** (corner T).
  - `e`: pre-mates `{d,f,g}`, task-mates `{f,g,k,l,m,n,o,p,q}`,
    post-mates `{f,g,k,l,m,n,o,p,q}`. `J_pp = 2/10 = 0.2`, `J_tp = 1.0`
    → **trace**.
  - `f, g`: same as `e` by symmetry → **trace**.
  - `h, i, j`: pre-mates ⊂ `{h,…,q}`, task and post mates `{d, others}`.
    `J_pp ≈ 0.3`, `J_tp = 1.0` → **trace** (weaker).
  - `a, b, c`: identical neighbours every phase → **persist**.

Same event, different reports: **CBR returns one named 4-leaf module**;
**MSPC returns 7 trace leaves but does not say which 7 belong
together**.

### 10.3 Why both exist in the investigation folder

They fail in different ways:

- **CBR is too strict at hard thresholds.** Legacy CBR returned zero
  trace exemplars in 4/9 patients (Pat_08, Pat_10, Pat_13, Pat_15)
  at `J_rpre ≤ 0.5 ∧ J_tl, J_tt ≥ 0.75`. Cohesion-CBR
  (`2026-04-25_cohesion-cbr.md`) softens this with the same
  corner-distance affinity trick MSPC uses, but keeps the rsPost
  subtree anchor.

- **MSPC is too granular.** It never names a module, so a reader
  cannot write "this DLPFC subtree traces" — only "these 7 leaves
  trace at scale `k = 10`". The trace-coherence diagnostic
  (`audit_14b`) confirms that MSPC trace leaves at the right `k`
  *do* co-cluster (purity ≈ 1.0 in 35/60 cohort cells), so the
  module is recoverable post-hoc — but not by MSPC itself.

### 10.4 Operational consequence

Use CBR (or its Cohesion variant) when the writeup needs to point at
*a specific subtree* — gallery exemplars, named cluster, "this
4-leaf group is the trace module." Use MSPC when the writeup needs
to point at *a leaf-level scale-resolved map* — multiscale strip,
"these leaves trace at this `k`, those at that `k`." For the
γ_l-vs-α/β reconciliation in §11 (TBD), the leaf-level map is
load-bearing because the dilution test requires per-leaf masks,
not pre-named modules.

### 10.5 Open theoretical question

A clean unified measure would (1) classify per-leaf like MSPC,
(2) auto-discover the module containing the trace leaves like CBR,
(3) avoid hard thresholds like Cohesion-CBR. We do not have such a
measure. CBR + MSPC together cover the question; neither alone does.
