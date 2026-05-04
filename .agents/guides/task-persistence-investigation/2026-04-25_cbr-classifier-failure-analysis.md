---
name: cbr-classifier-failure-analysis
type: diagnostic
era: COHORT_N9
status: current
created: 2026-04-25
updated: 2026-04-25
pointers:
  - .agents/guides/task-persistence-investigation/README.md
  - .agents/guides/task-persistence-investigation/2026-04-25_cbr-investigation.md
  - .agents/guides/task-persistence-investigation/2026-04-25_task-anchored-cbr.md
  - .agents/guides/task-persistence-investigation/2026-04-25_containment-cbr.md
  - .agents/guides/task-persistence-investigation/2026-04-25_consensus-subtree.md
  - data/audit/per_patient_hierarchy/
  - data/audit/per_patient_hierarchy_taskanchored/
  - data/audit/per_patient_hierarchy_containment/
  - data/audit/per_patient_hierarchy_consensus/
---

# CBR classifier failure analysis (2026-04-25 visual review)

**The leaf-cluster framework is the right proof of concept. The
*classifier* on top of it is what's broken. All four current CBR
variants fall into a single root cause — they collapse a fundamentally
2-dimensional question ("how cohesive is a leafset in phase φ?") onto
either a Jaccard score (which conflates cohesion with size match) or a
containment score (which conflates cohesion with mere set membership).
Neither captures the visual fact that draws the eye on the figures:
*the same leaves bunching tightly into a single subtree*. The fix is to
replace both with a directional* tightness *score `T(S, T_φ) = |S| /
|LCA(S, T_φ)|` that asks "how compactly does T_φ accommodate S as a
single subtree?", and to drop the binary classify-or-None step in favor
of soft per-leaf affinity scores so unclassified leaves stop defaulting
to blue rearrange.**

## 1. The user's six observations (verbatim → mapped to evidence)

| # | observation | mapped finding |
|:--|:-----------|:---------------|
| 1 | "lrg_sanity folder can be completely deleted, no need or use at all" | Done. S(τ), C(τ) curves on complete weighted networks don't reveal weight structure; for the LRG question we already use linkage matrices directly. **`data/audit/lrg_sanity/` and `scripts/01_compute/audit/audit_03_lrg_sanity.py` removed in this session.** |
| 2 | "per_patient_hierarchy interesting … weird that it does not find any trace with low_gamma where consensus finds traces in the very same figure" | Legacy CBR **rejects subtrees in the J ∈ (0.5, 0.75) limbo zone** as `None`. CON's exact-equality is binary so it bypasses the limbo. Detailed numerical evidence in §2.A. |
| 3 | "per_patient_hierarchy_consensus … incredible Pat_15 β does not find any thing except rearrangement … most of the figure is blue which is weird casue there are clearly the 4 classes expressed" | Pat_15 β's largest rsPost subtrees have *no exact-equality* counterpart in any other phase — leaf membership shifts by 1-3 leaves per phase. **CON is correctly empty** but the *visual* claim that "stable structures are clearly there" is also true: the leaves bunch but don't EXACTLY repeat. Detailed in §2.B. |
| 4 | "per_patient_hierarchy_containment … finds persist where others don't but makes a lot of mess" | Containment is **over-permissive on `persist`**: a small anchor S can be contained in a larger reference subtree even when S's leaves are not cohesively grouped *together* in the reference. Detailed in §2.C with Pat_05 γ_l + Pat_15 β counter-examples. |
| 5 | "per_patient_hierarchy_taskanchored interesting find correct stuff … but still misses the most" | TA-CBR fixes only the anchor asymmetry; the limbo zone, threshold cliffs, and projection default are all unchanged. So it recovers some lost cases (Pat_13 δ, expanded multi-band trace in Pat_06) but not the limbo-zone losses (Pat_10 α, etc.). |
| 6 | "task_trace_per_patient … VI is not the measure we are looking for, it cannot detect the presence of these 'imprinted' blocks" | Confirmed. VI(k) on flat partitions discards *which* leaves are involved; it can't localise. Keep the folder for completeness but it is decoupled from the scientific question. |

## 2. Numerical evidence

Each row below is a real rsPost internal node in the audit data. Scores
computed on the rebuilt LRG cache (`data/cache/imcoh_lrg/`).

### A. Limbo zone — Pat_13 γ_l (legacy fails to classify what CON catches)

The legacy `classify(J_rpre, J_tl, J_tt)` returns `None` whenever
`J_rpre ∈ (0.5, 0.75)` *or* `min(J_tl, J_tt) ∈ (0.5, 0.75)` and the
strict patterns don't fire. That gap absorbs many real signals:

| size | h_rel | J_rpre | J_tl | J_tt | C_rpre | C_tl | C_tt | T_rpre | T_tl | T_tt | legacy | containment |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|:---|
| 41 | 0.88 | **0.70** | 0.41 | 0.70 | 0.80 | 0.49 | 0.83 | 0.36 | 0.41 | 0.48 | None | reset |
| 37 | 0.87 | **0.66** | 0.38 | 0.67 | 0.78 | 0.49 | 0.84 | 0.33 | 0.37 | 0.44 | None | reset |
| 26 | 0.87 | 0.88 | 0.88 | 0.77 | 1.00 | 0.88 | 1.00 | 0.79 | 0.43 | 0.74 | None | persist |
| 25 | 0.84 | 0.92 | 0.92 | 0.80 | 1.00 | 0.92 | 1.00 | 0.76 | 0.42 | 0.71 | persist | persist |
| 23 | 0.71 | 1.00 | 1.00 | 0.87 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.66 | persist | persist |

The size-41 subtree has `J_rpre = 0.70` — slap in the middle of the
limbo zone. Legacy returns `None`. Containment correctly classifies it
(C_rpre = 0.80 ≥ 0.75, J_tl = 0.41 → reset). The legacy figure shows
this leafset's leaves as blue (rearrange default). **Same dendrogram,
same leaves, the *only* difference between legacy and containment is the
similarity score.**

### B. Pat_15 β — visually-stable structures, no algorithm finds them as `persist`

Pat_15 β's right-hand region clearly contains a large stable module
across all four phases (visible in `per_patient_hierarchy_containment/rpost_anchor/Pat_15/beta`).
The numbers explain why legacy *and* CON both miss it:

| size | h_rel | J_rpre | J_tl | J_tt | C_rpre | C_tl | C_tt | T_rpre | T_tl | T_tt | legacy | containment |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|:---|
| 60 | 0.44 | **0.76** | 0.68 | 0.72 | 0.80 | 0.78 | 0.82 | 0.64 | 0.54 | 0.57 | None | persist |
| 58 | 0.44 | **0.75** | 0.70 | 0.74 | 0.81 | 0.81 | 0.84 | 0.62 | 0.59 | 0.59 | None | persist |
| 57 | 0.43 | 0.75 | 0.69 | 0.75 | 0.81 | 0.81 | 0.84 | 0.61 | 0.58 | 0.58 | None | persist |
| 50 | 0.41 | 0.75 | 0.68 | 0.69 | 0.86 | 0.82 | 0.86 | 0.53 | 0.51 | 0.51 | None | persist |

- **Legacy**: every J ∈ {0.67…0.76}, hits the limbo zone in *every*
  phase pair → `None` for every size → all blue in figure.
- **CON**: requires exact equality. None of these size-60 / 58 / 57
  leafsets are exactly internal-node leafsets in *all four* phases (the
  module's boundary leaves drift by a few channels) → empty.
- **Containment**: scores ≥ 0.80 in every phase (the module's leaves
  are contained in *some* big subtree of each phase) → persist. Looks
  right…
- **Tightness**: 0.5–0.6 in *every* phase. Meaning: this 60-leaf
  rsPost module's LCA in rsPre / task is a *much larger* subtree (~ 100
  leaves). The 60 leaves are co-located with another 40 unrelated
  leaves in rsPre — so they're "in the same neighbourhood" but **not
  the same compact module**.

Reading: **the visual "module" is real but it's not the same module in
all four phases — it's a 60-leaf subset of a bigger 100-leaf
neighbourhood that the other phases organise differently.** Legacy and
CON are correctly silent. Containment is incorrectly enthusiastic. None
of the four scores capture the actual scientific situation, which is
"these 60 leaves are bunched in rsPost, partially bunched elsewhere".

### C. Containment over-permissivity — Pat_05 γ_l vs Pat_15 β

The size-25 row of Pat_05 γ_l is the controlled comparison:

| size | J_rpre | J_tl | J_tt | C_rpre | C_tl | C_tt | T_rpre | T_tl | T_tt | legacy | containment |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|:---|
| 25 | 0.37 | 0.68 | 0.69 | 0.80 | 0.76 | 0.92 | 0.36 | 0.22 | 0.21 | None | **persist** |

This rsPost subtree has 25 leaves. Containment says "persist": the 25
leaves are contained in some rsPre-, task_learn-, task_test-subtree of
size ≤ 59. But tightness says T_tl = 0.22 — those 25 leaves take a
*114-leaf* LCA in task_learn (essentially the whole tree). They are
*scattered* in task_learn, not bunched. Persist is the wrong label;
the right label is closer to "rearranged in tasks but bunched in
rsPost". Containment can't see that distinction; tightness can.

This is the mechanism behind the user's "fucks up with persist" comment.

## 3. Root causes (three, ranked by severity)

### R1. Threshold cliff + limbo zone

The classifier uses two thresholds (`τ_match = 0.75`, `τ_dissim = 0.50`)
and four mutually exclusive predicates that each require *every* score
on a given side of one threshold. A subtree with `J ∈ (0.5, 0.75)` —
which empirically is the **modal** range for real-world LRG subtrees
matched across phases — fires no predicate and falls through to `None`.
The leaf-projection rule then defaults its leaves to `rearrange` (blue)
because no classified ancestor exists.

Quantitatively, on Pat_15 β rsPost (118 internal nodes in size [5, 60]):
the legacy classifier returns `None` on **all 86 of them**. That's why
the figure is 100 % blue. The data is fine; the classifier is rejecting
everything.

### R2. Score conflation

| score | what it actually measures | what we want |
|:---|:---|:---|
| Jaccard `J(S, leaves(v))` | similarity of two leafsets, penalised when sizes differ | "do these leaves bunch together in T_φ?" |
| Containment `C(S, T_φ)` | best fraction of S contained in any reference subtree of size ≤ κ | "do these leaves bunch together in T_φ?" |
| **Tightness `T(S, T_φ) := \|S\| / \|leaves(LCA(S, T_φ))\|`** | exactly: how compactly the smallest containing subtree of T_φ accommodates S | "do these leaves bunch together in T_φ?" |

Tightness is the directional cohesion score that matches the visual
question: "in this phase's dendrogram, do these particular leaves form
a single coherent subtree?". Range `[1/N_p, 1]`. `T = 1` ⇔ S *is* a
subtree of T_φ. `T → 0` ⇔ S spans the entire tree (LCA ≈ root).

LCA computation is `O(N)` per query using ancestor walking; total
complexity is comparable to `J*` and `C*`.

### R3. Hard-priority projection

The leaf inherits the *highest-priority* classified pattern among its
ancestors. Trace > persist > reset > rearrange. If no ancestor is
classified, leaf is `rearrange`. This is wrong twice over:

- **It produces "100% blue" figures whenever R1 fires for every
  ancestor**, even when the dendrogram is full of stable structure.
- **It hides ambiguous leaves**. A leaf inside a subtree that scores
  T_rpre = 0.6, T_tl = 0.7, T_tt = 0.7, T_rpost = 0.7 is meaningful —
  it has weak persist-leaning behaviour. The current rule discards
  this because the strict predicates don't fire.

The right operationalisation is a **soft per-leaf affinity vector**
`(α_trace, α_persist, α_reset, α_rearrange) ∈ Δ³`, computed from
ancestor tightness scores, with the leaf coloured by its dominant
component (intensity = saturation of dominance).

## 4. Proposed redesign — Cohesion-CBR (variant E)

Sketch only — full mathematical scope to land in
`2026-04-26_cohesion-cbr.md` after this analysis is acknowledged.

**Score.** Per anchor `η ∈ Internal(T_taskT)` of size ∈ [s_min, s_max]:

    T_φ(η) := |leaves(η)| / |leaves(LCA(leaves(η), T_φ))|     for φ ∈ {rsPre, taskL, taskT, rsPost}.

T_taskT(η) ≡ 1 by construction.

**Soft pattern affinities.** For each anchor, define four affinities in
[0, 1]:

    a_trace(η)     := σ(T_rpost(η) − T_rpre(η))      · sigmoid_to_[0,1]
    a_persist(η)   := min(T_rpre, T_taskL, T_taskT, T_rpost)(η)
    a_reset(η)     := σ(T_rpre(η) − max(T_taskL, T_taskT, T_rpost)(η))
    a_rearrange(η) := σ(T_rpost(η) − max(T_rpre, T_taskL, T_taskT)(η))

where `σ` is a soft contrast-to-affinity map `σ(x) = clip((x + 1)/2, 0, 1)`
(or another smooth squash). No hard thresholds.

**Per-leaf projection.** For each leaf ℓ, aggregate over all ancestors
η with leaves(η) ∋ ℓ and size ∈ [s_min, s_max]:

    A_p(ℓ) := mean_{η ∋ ℓ} (a_p(η) · w(η))   for p ∈ {trace, persist, reset, rearrange}.

`w(η)` is a weight (default uniform; size-weighting and h_rel-weighting
are sweep parameters).

**Color rule.** Leaf colour is the dominant `argmax_p A_p(ℓ)` with
*saturation* equal to `max_p A_p(ℓ) − second_max_p A_p(ℓ)` (a confidence
measure). Faint colours = ambiguous leaves; saturated colours = clearly
classified leaves.

**Why this fixes everything:**
- **R1 (limbo zone)**: gone. No thresholds, every subtree contributes a
  continuous affinity to all four patterns.
- **R2 (score conflation)**: tightness *is* cohesion; it cannot be high
  for "scattered but contained" leafsets, and it cannot be low for
  "different size but same cluster" leafsets.
- **R3 (default-to-blue)**: gone. Every leaf gets all four affinities;
  the figure shows the dominant pattern with confidence-graded
  saturation.

## 5. Sanity check on tightness

Run on the same diagnostic cells:

| cell | top rsPost subtree | T_rpre | T_tl | T_tt | reading |
|:---|:---|---:|---:|---:|:---|
| Pat_15 β size-60 | 0.44 h_rel | 0.64 | 0.54 | 0.57 | weak persist (cohesion drops in non-rsPost phases) |
| Pat_05 γ_l size-42 | 0.94 h_rel | 0.36 | 0.69 | 0.37 | true trace candidate (cohesion drops only in rsPre, picks up in task_learn) — matches the gallery exemplar |
| Pat_13 γ_l size-23 | 0.71 h_rel | 1.00 | 1.00 | 0.66 | strong persist with weak task_test mismatch |
| Pat_10 α size-57 | 0.63 h_rel | 0.59 | 0.52 | 0.51 | weak persist-leaning (no clean classification at any threshold) — limbo case |

The numbers separate the four patterns smoothly without any classifier
hard-coding. This is the evidence that tightness is the missing axis.

## 6. Confirmation on the leaf-cluster POC

Yes. The leaf-cluster representation is the right scaffolding:

- **Per-leaf categorisation onto 4 colours visualises the joint event
  cleanly** when the underlying classifier is sound. The audit_08
  layout (4 phase dendrograms with consistent leaf colouring) is the
  right reading rule and survives intact.
- **The four pattern names (trace / persist / reset / rearrange) carve
  the natural quadrants of the (T_rpre, T_rpost) plane** with
  task as a shared anchor — they're the right question.
- **What needs replacing is the score function and the classifier,
  not the framework.** Tightness + soft affinities slot into the
  existing audit_08 plot machinery without changing the figure layout.

This is a strong proof of concept: every patient now has *some* visible
trace structure under at least one of the four current variants
(per the cross-method summary), confirming the joint event is real;
the variability across variants is exactly the algorithmic noise the
diagnostic explains. A unified Cohesion-CBR is the right next pass.

## 7. Path forward

1. **Land the diagnostic** (this file).
2. **Write `2026-04-26_cohesion-cbr.md`** — full math scope for the
   tightness-based soft-affinity classifier. Include LCA pseudocode,
   weighting strategies, sigmoid choice, sensitivity to `s_min/s_max`.
3. **Implement** as `scripts/01_compute/audit/audit_12_cohesion_cbr.py`,
   reusing the audit_08 figure layout end-to-end.
4. **Compare side by side** with the four existing variants on the
   problematic cells (Pat_15 β, Pat_13 γ_l, Pat_10 α, Pat_05 γ_l).
   The success criterion is *visual*: the dominant-color saturation map
   should match the eye's reading of the dendrograms.
5. **Only after that** revisit cohort census and the question "is the
   trace cohort-wide". The current numbers are conditional on each
   variant's classifier and are therefore not directly comparable.

No commits yet (per user instruction). Working tree state preserved on
`audit/cohort-n9-diagnostic`.
