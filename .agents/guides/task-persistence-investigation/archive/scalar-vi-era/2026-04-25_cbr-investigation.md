---
name: cbr-investigation
type: investigation
era: COHORT_N9
status: active
created: 2026-04-25
updated: 2026-04-25
pointers:
  - .agents/guides/task-persistence-investigation/README.md
  - .agents/guides/task-persistence-investigation/2026-04-25_module-retention-landscape.md
  - .agents/guides/task-persistence-investigation/2026-04-25_task-anchored-cbr.md
  - .agents/guides/task-persistence-investigation/2026-04-25_containment-cbr.md
  - .agents/guides/task-persistence-investigation/2026-04-25_consensus-subtree.md
  - data/audit/per_patient_hierarchy/
  - data/outputs/figures/cbr_gallery/
---

# CBR investigation — single-module-matching family

**Umbrella entry for the cluster-birth-retention (CBR) family of
descriptive measures. Each measure asks the same shape of question — *is
this subtree a coherent module in some phase and not in another?* — and
differs only in (i) which dendrogram is the anchor, (ii) what similarity
score replaces leafset Jaccard, (iii) how strict the threshold is. This
file indexes the four CBR variants drafted on 2026-04-25 and tracks the
session-level decisions; the individual scope reports carry the math.**

## 1. The shape of the question

For every (patient, band), the LRG produces four dendrograms — one per
phase. A "task-induced module that persisted" is a leafset `L ⊆ L_p`
that:

- **was not** a coherent subtree of `T^{p,b,rest_pre}`, AND
- **is** a coherent subtree of `T^{p,b,task_test}` (and ideally
  `T^{p,b,task_learn}` too), AND
- **is still** a coherent subtree of `T^{p,b,rest_post}`.

Every CBR variant operationalises the three predicates differently:

- **What counts as "coherent in this phase"?** Set-membership in the
  leafset family (CON), best-Jaccard ≥ τ (legacy / TA-CBR), or
  best-containment ≥ τ (B).
- **Where do we look for candidates?** Iterate over `task_test`'s
  internal nodes (TA), `rest_post`'s internal nodes (legacy / B-rpost),
  or both (B-task and CON).

## 2. Why we are revisiting CBR right now

The legacy CBR (`scripts/archive/2026-04_failed-scalar-session/gallery_cbr_modules.py`)
plus its leaf-projection (`scripts/01_compute/audit/audit_08_per_patient_hierarchy.py`)
returned **zero CBR-trace exemplars** in 4 of 9 patients (Pat_08, Pat_10,
Pat_13, Pat_15) at the strict thresholds `J_rpre ≤ 0.5 ∧ J_tl, J_tt ≥
0.75`. The visual census (`data/audit/per_patient_hierarchy/{Pat_XX}/`)
confirms — those four patients show only PERSIST and RESET modules in
the colored dendrograms. No red.

User concern (paraphrased): the *algorithm* might be missing real trace
subtrees because of the symmetric Jaccard, the rest_post anchor, and
the strict thresholds. Five concrete failure modes were enumerated
(see chat 2026-04-25): anchor asymmetry, Jaccard's containment
penalty, single-best-match, hard-threshold cliffs, size window. The
CBR family this scope catalogues is the systematic exploration of
those failure modes, one variant per mode.

## 3. Variants drafted in this session

| variant | scope file | anchor | similarity | thresholds | failure modes addressed |
|---:|:---|:---|:---|:---|:---|
| Legacy | (no scope; `audit_08`) | `T_rsPost` | symmetric Jaccard | `J ∈ {0.5, 0.75, 0.8, 0.4}` | none — baseline |
| **A. TA-CBR** | [`2026-04-25_task-anchored-cbr.md`](2026-04-25_task-anchored-cbr.md) | `T_taskT` | symmetric Jaccard | identical to legacy | failure 1 (anchor asymmetry) |
| **B. Containment-CBR** | [`2026-04-25_containment-cbr.md`](2026-04-25_containment-cbr.md) | both `T_taskT` and `T_rsPost` | directional `C(S, T) = max\|S∩leaves(v)\|/\|S\|` | identical to legacy (re-interpreted on `C`) | failures 1 (via dual-anchor) + 2 (containment) |
| **D. Consensus-subtree** | [`2026-04-25_consensus-subtree.md`](2026-04-25_consensus-subtree.md) | both | exact set equality (no score) | none — strict membership | strictest possible; tests whether the absence in Pat_08/10/13/15 is real |
| MRL | [`2026-04-25_module-retention-landscape.md`](2026-04-25_module-retention-landscape.md) | `T_taskT` | symmetric Jaccard via `match_τ` | swept τ | scalar-field aggregate of TA-CBR; sibling, not nested |

Variants C, E, F, G, H from the parent menu (chat 2026-04-25) are *not*
drafted in this session — they leave the CBR/single-module-matching
shape of the question and require their own scope reports if/when
needed.

## 4. Run order (this session)

User instruction (2026-04-25): *"do A and B first and then while I judge
the results also do D"*.

1. ✅ Scope A drafted → `2026-04-25_task-anchored-cbr.md`.
2. ✅ Scope B drafted → `2026-04-25_containment-cbr.md`.
3. ✅ Scope D drafted → `2026-04-25_consensus-subtree.md`.
4. ⏳ **Implement A** → `scripts/01_compute/audit/audit_09_task_anchored_cbr.py`,
   outputs under `data/audit/per_patient_hierarchy_taskanchored/`.
5. ⏳ **Implement B** → `scripts/01_compute/audit/audit_10_containment_cbr.py`,
   outputs under `data/audit/per_patient_hierarchy_containment/{rpost_anchor,task_anchor}/`.
6. ⏳ **Implement D** (after user judges A/B) →
   `scripts/01_compute/audit/audit_11_consensus_subtree.py`, outputs
   under `data/audit/per_patient_hierarchy_consensus/`.

All four CBR variants reuse the audit_08 figure layout (4 phase
dendrograms, branches/leaves coloured, codebar) so they can be
inspected side by side. Colour palette is fixed:
`trace=#d62728, persist=#000000, reset=#7f7f7f, rearrange=#1f77b4`
(user-pinned).

## 5. What we expect to learn

| measure | informative if … |
|:---|:---|
| TA-CBR vs. legacy | If task-anchoring surfaces trace exemplars in Pat_08/10/13/15, the legacy "no trace" verdict is an anchor-asymmetry artifact. |
| Containment vs. Jaccard | If containment surfaces exemplars in the same patients but Jaccard does not, the legacy verdict is a containment-penalty artifact. |
| CON | If CON returns zero across the entire cohort, exact subtree-equality is too strict and we know the load-bearing predicate must be a similarity score, not equality. If CON returns even one exemplar, we have an undeniable proof-of-existence. |

If A, B, AND D all confirm the legacy verdict ("no trace in
Pat_08/10/13/15"), the cohort heterogeneity is a *real* property of
the data, not an algorithmic artifact, and the next move is no longer
about the metric — it is about the patients (probe coverage, signal
quality, task engagement).

## 6. Where this lives

- **Scope reports**: this folder (`task-persistence-investigation/`).
- **Implementation**: `scripts/01_compute/audit/audit_{09,10,11}_*.py`.
- **Outputs**: `data/audit/per_patient_hierarchy_{taskanchored,containment,consensus}/`.
- **Legacy comparison**: `data/audit/per_patient_hierarchy/`
  (audit_08, already on disk).
- **Cohort context**: `data/audit/REPORT.md` (audit_07).

The README's Active-scope-reports table carries one line for the entire
CBR family pointing here; per-variant scope files are linked from §3
above, not from the README, to keep the index uncrowded.
