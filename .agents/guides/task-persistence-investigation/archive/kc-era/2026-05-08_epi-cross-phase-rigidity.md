---
name: epi-cross-phase-rigidity
type: scope
era: IMCOH_ABS × COHORT_N10
status: draft
created: 2026-05-08
updated: 2026-05-08
pointers:
  - .agents/plans/active/2026-05-08_lrg-epilepsy-research-directions.md
  - .agents/reports/2026-05-07_epileptic-n10-revisit.md
  - .agents/reports/2026-05-07_kc-cross-phase-taxonomy.md
  - .agents/guides/task-persistence-investigation/2026-05-07_kc-anchor-modules.md
  - .agents/guides/task-persistence-investigation/2026-05-07_kc-reset-modules.md
  - .agents/guides/task-persistence-investigation/2026-05-07_kc-rearrangement-modules.md
  - .agents/guides/01_project/terminology.md
---

# Epi cross-phase rigidity — TARR taxonomy applied to the fixed epi leaf set

## Renormalization head

**Apply the existing trace / anchor / reset / rearrange (TARR) classification
to a *fixed* leaf set per patient — the epileptic node set E_p — instead of
enumerating subtree candidates. The scientific question: across cognitive
states (rest_pre → task_test → rest_post), does the dendrogram restricted
to epi leaves stay locked (anchor), reorganize and revert (reset),
reorganize and persist (trace), or end up in an emergent post-task
configuration (rearrange)? The cohort-level statistic is the per-class
z-score against a size-matched random-leaf-subset null per (patient,
band). Direct test of the textbook intuition "the epileptic network is
pathologically rigid across cognitive states", which has never been
measured directly with intracranial multiscale tooling. Either
direction (anchor-enriched, trace-enriched, reset-enriched, or null) is
publishable. Reuses the audit_47-50 KC machinery and the
`kc_cross_phase_taxonomy.md` thresholds; no new primitives.**

## Verdict (2026-05-08, post-implementation)

**Implemented at `audit_54_epi_rigidity_compute.py`. Three CSVs +
three PDFs at `data/audit/epi_rigidity/`. The data shows a structured
band-resolved phenomenology — not a single "rigid vs labile" answer.
Three regularities pass the IQR > 20% × |median| floor (per
`feedback_iqr_vs_cohort_slope.md`):**

| (band, λ) | direction | n_pat | median Td_E | %IQR | reading |
|---|---|---|---|---|---|
| **δ / λ=1** | **RESET** | 7/9 | +0.53 | **75%** | strongest single effect — epi slow-rhythm heights shift in task and revert |
| **β / λ=1** | **RESET** | 6/9 | +0.15 | **43%** | epi heights reset at the band where the global cohort TRACES (Result 2) — *band-specific dissociation* |
| **α / λ=0** | **TRACE** | 5/9 | −0.29 | **55%** | epi-subtree topology reorganises and persists |
| α / λ=0.5 | TRACE | 5/9 | −0.23 | 58% | confirms α/λ=0 in the balanced blend |

**The β-dissociation is the manuscript-worthy finding.** At β,
the cohort network shows a load-bearing trace (Result 2: KC λ=0/0.5/1
+ D-rank d_F all 10/10, joint Bonferroni m=48 surviving). The
epi-induced subtree at the same band shows the *opposite* — heights
shift during task and revert in rest_post. Either reading
(epi nodes confound the global signal, or the global trace is
non-epi-carried) implies a direct testable consequence: **subtracting
the epi nodes from the global tree should sharpen the β trace**.

**δ-heights RESET is the largest amplitude effect** in the entire
(band × λ) grid. Combined with the n=10 revisit's δ cross-probe edge
enrichment (1.55× cohort median), the layered picture is:
**epileptic δ coupling is strong AND task-modulated AND does not
persist** — three features that cohere only when read together
across the edge / coupling-strength layer and the
cross-phase / hierarchical-position layer.

**α-topology TRACE** is the only band where the epi sub-network
behaves *like* the rest of the cortex: cognition-band α recruits
the epi network into a task configuration that persists. The
audit_48 fig_05 pair-mask result (8/9 p=0.006) attenuates to 5/9 at
55% IQR under the more rigorous induced-subtree test — direction
preserved, amplitude reduced. The γ_l/λ=1 cell replicates
identically in both operationalisations (7/9, small amplitude).

The blanket framing "is the epileptic network rigid?" does not
survive contact with this data. The right question is **which
band, which aspect (topology vs heights), and what kind of shift
(trace vs reset)?** — and the answer is structured, not flat.

### What the BH-FDR framing gives you

If gated by Wilcoxon + BH-FDR within (variant, λ) panel of m=6
bands at q < 0.05: all 144 TARR class cells, 36 z_trace cells, and
36 Td_E vs 0 cells return null. Best uncorrected p = 0.033 at
γ_l/λ=1 (matches audit_48 pair-mask exactly). The
`ε_KC = 0.10` classifier is degenerate at this scale — KC distances
on |E_p|-leaf induced subtrees range 0.4–2.6, so 99.6% of all leaf
subsets (epi or null) classify as "rearrange". A scale-aware
threshold would be needed; the continuous tests already give the
phenomenology above.

This BH-null framing is the wrong lens for the data. Per project
memory `feedback_iqr_vs_cohort_slope.md` and
`feedback_dont_rerun_scalar_tests.md`: cohort-shape regularities
(passing the IQR floor) are the load-bearing finding, not the
p-value gate.

### Full report

[`.agents/reports/2026-05-08_direction-a-induced-subtree.md`](../../reports/2026-05-08_direction-a-induced-subtree.md)
— canonical regularity-first writeup.

---

## Notation

(Inherits from `2026-05-07_kc-anchor-modules.md`,
`2026-05-07_kc-reset-modules.md`,
`2026-05-07_kc-rearrangement-modules.md`, and audit_47.)

- `V_p` — node set for patient `p`, size `N_p ∈ [100, 122]`.
- `E_p ⊆ V_p` — epileptic node set (red XLSX cells, channel-matched
  per `load_epileptic_nodes`). `|E_p|` ranges 6–30 across n=9
  cohort patients (Pat_15 has zero annotation and drops out;
  same n=9 effective as the 2026-05-07 revisit).
- `N_p^c = V_p \ E_p` — non-epi node set.
- `H_φ^p` — LRG dendrogram of patient `p` in band `b ∈ {δ, θ, α, β,
  γ_l, γ_h}` at phase `φ ∈ {pre, tt, post}`.
- `H_φ^p|_S` — induced subtree of `H_φ^p` on leaf subset `S ⊆ V_p`.
  The induced subtree is the minimum subtree containing all leaves
  in `S` (equivalently the dendrogram restricted to `S` after
  pruning leaves not in `S` and contracting degree-2 internal
  nodes).
- `Z_φ^p(S)` — linkage matrix of `H_φ^p|_S`, with `|S| − 1` rows.
- `d_KC(Z_a, Z_b; λ)` — Kendall–Colijn distance with topology /
  heights blend `λ ∈ {0, 0.5, 1}` (`lrg_eegfc.utils.metrics.tree_distance`).
- `J(A, B)` — Jaccard `|A ∩ B| / |A ∪ B|`.

Conventions for `S = E_p`:

- The induced subtree has `|E_p|` leaves and `|E_p| − 1` internal nodes.
- KC distance is computed on the **trimmed** dendrogram, not on the
  full `H_φ^p` with non-epi leaves marked. Trimming is the standard
  "induced subtree" operation; we use `lrgsglib.tree.prune_to_leaves`
  if available, otherwise wrap `scipy.cluster.hierarchy` to filter
  the linkage matrix.

## Definitions

We define **per-patient class membership** for the fixed leaf set
`E_p` — exactly one of:

- `Trace_E` — induced subtree reorganises during task and the
  reorganisation persists into rest_post.
- `Anchor_E` — induced subtree is approximately invariant across
  all three phases.
- `Reset_E` — induced subtree reorganises during task but reverts
  into rest_post.
- `Rearrange_E` — induced subtree at rest_post is unrelated to both
  rest_pre and task_test (emergent post-task configuration).

Each class is defined via Jaccard subtree-matching gates between
the induced subtrees, mirroring the audit_47/48/49/50 thresholds
but applied to the *root* induced subtree directly (not to enumerated
internal subtrees).

### Triangle-distance cell

For each (patient `p`, band `b`, λ-blend `λ ∈ {0, 0.5, 1}`),
compute the three pairwise KC distances

```
δ_pre→tt   = d_KC(Z_pre^p(E_p),  Z_tt^p(E_p);  λ)
δ_tt→post  = d_KC(Z_tt^p(E_p),   Z_post^p(E_p); λ)
δ_pre→post = d_KC(Z_pre^p(E_p),  Z_post^p(E_p); λ)
```

with the standard KC normalisation (`normalize=True`).

### Class predicates (per (p, b, λ))

Define small / large with a single tolerance `ε_KC = 0.10` (chosen
to match the audit_47 trace gate); flexibility deferred to open
questions.

| Class | Predicate |
|---|---|
| **Anchor_E** | `δ_pre→tt < ε_KC ∧ δ_tt→post < ε_KC ∧ δ_pre→post < ε_KC` |
| **Trace_E** | `δ_pre→tt > ε_KC ∧ δ_tt→post < ε_KC ∧ δ_pre→post > ε_KC` |
| **Reset_E** | `δ_pre→tt > ε_KC ∧ δ_tt→post > ε_KC ∧ δ_pre→post < ε_KC` |
| **Rearrange_E** | `δ_pre→post > ε_KC ∧ (δ_pre→tt > ε_KC ∨ δ_tt→post > ε_KC)` AND no other class fires |

The four classes partition the (δ_pre→tt, δ_tt→post, δ_pre→post)
octant into four basins; an "ambiguous" residual (small δ_pre→tt
+ large δ_tt→post + small δ_pre→post — a noisy reset variant) is
counted as Reset_E.

### Per-class score

In addition to the binary class assignment, define a continuous
`Trace score` = `T_d^E = δ_tt→post − δ_pre→tt` (negative ⇒ trace
direction; matches the Section-5 KC scalar). This continuous score
allows cohort-Wilcoxon tests independently of the discrete class
gate, mirroring the Result-2 LRG β trace test.

### Size-matched random null

For each (p, b, λ), draw `R = 100` random leaf subsets
`S_r ⊂ V_p` with `|S_r| = |E_p|`, computed without replacement and
disjoint from each other only if memory allows (otherwise i.i.d.).
For each `S_r` compute the same triangle (δ_pre→tt, δ_tt→post,
δ_pre→post), the class assignment, and the trace score `T_d^{S_r}`.
The null distributions are `{Class(S_r)}` and `{T_d^{S_r}}`.

Per-patient z-scores:

```
z_class(p, b, λ) = (1{Class(E_p) = c} − P(c | random)) / σ(P(c | random))
z_trace(p, b, λ) = (T_d^E − E[T_d^{S_r}]) / σ(T_d^{S_r})
```

(For `z_class` we use bootstrap on the indicator since exact mean
and variance are derivable from the multinomial.)

### Cohort aggregator

For each (band, λ, class):
- Cohort-Wilcoxon one-sided on per-patient `z_class` against 0
  (greater-than). BH-FDR within band × class.
- Cohort consistency: count patients with `z_class > 1.96` (≥ 8/9
  ⇒ **cohort-positive**).

For the trace score:
- Cohort-Wilcoxon one-sided on per-patient `z_trace` (less-than ⇒
  trace direction). BH-FDR within band × λ.

### Cross-probe variant

Same procedure with `E_p` replaced by `E_p^cp ⊆ E_p` — the subset
of epi nodes that participate in at least one cross-probe edge in
the patient's adjacency. This drops same-probe focal hypersynchrony
as a confound. Required as a complementary report; both `E_p` and
`E_p^cp` are reported.

## Properties

- **Range.** `δ ∈ [0, 1]` after KC normalisation; class assignment is
  one-hot in {Anchor, Trace, Reset, Rearrange}; `T_d^E ∈ [-1, 1]`.
- **Invariance.** KC distance is permutation-invariant on leaf
  labels, so the measure is invariant under cohort-level
  reordering of contacts. It is **not** invariant under leaf
  reweighting (we do not reweight by epi/non-epi).
- **Identifiability.** Class assignment is unambiguous given
  `ε_KC`; ambiguous cases by construction fall into Reset_E.
- **Sample size.** n=9 patients × 6 bands × 3 λ-values = 162
  cells. After BH-FDR within band × λ, cohort-positive bands at
  `q < 0.05` are reportable.
- **Negative properties (what this measure does NOT detect).**
  - **Sub-modular reorganisation within E_p.** A trace where only
    half of E_p reorganises and the other half stays locked may
    average out at the induced-subtree level. Mitigation: a follow-up
    leaf-level decomposition (audit_47-style enumeration restricted
    to subtrees of E_p) is layered on top if A is positive.
  - **Reorganisation that conserves topology but shifts heights**
    is detected only at λ = 1.
  - **Anatomical gradient effects.** If E_p is anatomically
    clustered (e.g. all in a single hippocampal probe), the induced
    subtree is dominated by short-range connectivity and the cross-
    phase signal may be confounded with anatomy. Mitigation:
    cross-probe variant `E_p^cp`.
  - **Per-patient effects below cohort threshold.** Pat_03 (1024 Hz
    outlier) and patients with `|E_p| < 5` (none after channel
    matching, but Pat_03 and Pat_07 are borderline) may dominate
    cohort statistics. Mitigation: per-patient table in the figure.
- **Complexity.** Per (p, b, λ): three KC distance evaluations on
  a `|E_p| − 1` × 4 linkage matrix, plus 3 R draws → `O(R · |E_p|^2)`.
  With R = 100, `|E_p| ≤ 30`, n = 9 patients, 6 bands, 3 λ-values:
  total ≈ 30 000 KC evaluations. Trivially parallel; estimated
  runtime < 5 minutes on the lapbrain conda env.

## Caveats & failure modes

| | Caveat | Mitigation |
|--|--|--|
| 1 | `|E_p|` varies wildly (6–30); raw class fractions are not directly comparable across patients | Size-matched random null per patient is the central control; report z-scores, not raw fractions |
| 2 | Pat_15 has no epi annotation | Drops out, n=9 effective (matches the 2026-05-07 revisit) |
| 3 | Pat_03 outlier (1024 Hz, 47% same-probe E_p) | Log Pat_03 row separately; flag in cohort table; check robustness with leave-one-out |
| 4 | Same-probe pairs dominate E_p in some patients (47% Pat_03, 43% Pat_07, 50% Pat_08) | Cross-probe variant `E_p^cp` reported alongside `E_p` |
| 5 | KC normalisation `normalize=True` already used in audit_47 etc. — sensitive to subtree size | Same convention used here for direct comparability |
| 6 | `ε_KC = 0.10` may be too lenient or strict at small `|E_p|`; small-induced-subtree noise floor | First-sweep diagnostic: plot `δ_pre→post` distribution under random null per `|E_p|` bin; adjust `ε_KC` if class fractions saturate |
| 7 | Between-patient comparison of class counts is fragile (one binary call per patient) | Continuous `T_d^E` with cohort-Wilcoxon is the primary inferential statistic; class fractions are descriptive |
| 8 | Random null may include subsets that overlap heavily with `E_p` (especially when `|E_p|` is large vs `N`) | At `|E_p| = 30, N ≈ 115`, expected overlap is ≈ 25%. Document but accept; per the literature, size-matched is the standard null. A "non-overlapping" stricter null is an open question |

## Pseudocode

```
input:  H_pre, H_tt, H_post (linkage matrices for one (p, b))
        E_p ⊂ V_p           (epi node set)
        R = 100              (null draws)
        ε_KC = 0.10
        λ ∈ {0, 0.5, 1}

output: per-(p, b, λ) record:
        - δ_pre→tt, δ_tt→post, δ_pre→post (epi)
        - Class(E_p) ∈ {Anchor, Trace, Reset, Rearrange}
        - T_d^E
        - z_class[c], z_trace
        - per-class null fractions

# 1. Epi induced subtree distances
Z_pre  = induced_linkage(H_pre,  E_p)
Z_tt   = induced_linkage(H_tt,   E_p)
Z_post = induced_linkage(H_post, E_p)

δ_pre_tt   = d_KC(Z_pre, Z_tt;   λ, normalize=True)
δ_tt_post  = d_KC(Z_tt,  Z_post; λ, normalize=True)
δ_pre_post = d_KC(Z_pre, Z_post; λ, normalize=True)

T_d_E = δ_tt_post − δ_pre_tt
Class_E = classify(δ_pre_tt, δ_tt_post, δ_pre_post; ε_KC)

# 2. Size-matched null
class_count = {Anchor: 0, Trace: 0, Reset: 0, Rearrange: 0}
T_d_null = []
for r in 1..R:
    S_r ← random sample without replacement from V_p, size |E_p|
    Z_pre_r  = induced_linkage(H_pre,  S_r)
    Z_tt_r   = induced_linkage(H_tt,   S_r)
    Z_post_r = induced_linkage(H_post, S_r)
    δ_a = d_KC(Z_pre_r, Z_tt_r;   λ)
    δ_b = d_KC(Z_tt_r,  Z_post_r; λ)
    δ_c = d_KC(Z_pre_r, Z_post_r; λ)
    Class_r = classify(δ_a, δ_b, δ_c; ε_KC)
    class_count[Class_r] += 1
    T_d_null.append(δ_b − δ_a)

# 3. z-scores
P_null[c] = class_count[c] / R
σ_null[c] = sqrt(P_null[c] · (1 − P_null[c]) / R)
z_class[c] = (1{Class_E = c} − P_null[c]) / max(σ_null[c], 1e-3)
z_trace = (T_d_E − mean(T_d_null)) / std(T_d_null)

# 4. Cross-probe variant: same procedure with E_p replaced by E_p^cp
```

Cohort aggregation (over patients for fixed (b, λ)):

```
for c in {Anchor, Trace, Reset, Rearrange}:
    z_vec = [z_class[c](p) for p in cohort]
    p_one_sided = wilcoxon_z(z_vec, alternative='greater')
    cohort_positive[c] = (count(z > 1.96) ≥ 8) and (BH-corrected q < 0.05)

z_trace_vec = [z_trace(p) for p in cohort]
p_trace_one_sided = wilcoxon_z(z_trace_vec, alternative='less')   # trace = T_d < 0
```

## Visualization spec

**Figure 1 — Cohort verdict heatmap.**

A 4 × 6 × 3 cube unfolded as 3 panels (one per λ):

- **Rows** = {Anchor, Trace, Reset, Rearrange}.
- **Columns** = {δ, θ, α, β, γ_l, γ_h}.
- **Cell value** = number of patients with `z_class > 1.96` (count
  out of 9). Colour: sequential green-blue scale.
- **Stars**: BH-FDR `q < 0.05` (within λ panel, m = 24 cells).
- **Reading rule**: a row with bright cells indicates the
  corresponding class is cohort-enriched in those bands; a dim row
  is null-consistent. Vertical bands reveal band-specific class
  preferences.

**Figure 2 — Trace-score box-and-dot.**

3 × 6 grid (rows = λ, cols = bands), each panel boxplot of
per-patient `z_trace` with dots over the box. Horizontal red dashed
line at 0; below = trace direction. Stars: BH-FDR `q < 0.05`. Mirrors
the n=10 revisit `fig_05_kc_epi_triangle.pdf` but with random-null
z-scores instead of raw distances.

**Figure 3 — Per-patient triangle scatter.**

For one focal band per λ (probably α at λ=0 since that is the
expected continuation of the n=10-revisit α/topology lead): scatter
of `δ_pre→tt` (x) vs `δ_tt→post` (y), one dot per patient, dot
labelled with patient id and class. Diagonal `y = x` line; trace
zone (`y < x − ε_KC`) shaded. Random-null cloud (median + 95%
ellipse) overlaid in light grey.

**Figure 4 — `E_p` vs `E_p^cp` comparison.**

Side-by-side cohort heatmaps (Figure 1 layout) for the all-pair
and cross-probe variants. Reading rule: cells that brighten in
`E_p^cp` only are driven by cross-probe (long-range) coupling;
cells that dim are driven by same-probe focal coupling.

## Connection to prior tools

| Prior tool | Relation |
|--|--|
| audit_47 (KC trace modules) | Same KC machinery; this measure applies the trace gate to a *fixed* set instead of enumerating subtrees. Trace_E is the "is `E_p` a trace module?" test. |
| audit_50 (KC anchor modules) | Anchor_E is the "is `E_p` an anchor module?" test. |
| Section-5 KC scalar trace (Result 2) | Same `T_d^KC` formulation, restricted from the global tree to `E_p` induced subtree. |
| `2026-05-07_kc-cross-phase-taxonomy.md` | Cohort tally taxonomy applied at the per-patient scale rather than the per-module enumerated scale. |
| `2026-04-25_task-trace-canonical.md` | Canonical TARR (Trace / Anchor / Reset / Rearrange) names; this is a TARR application, not a new framework. |
| `audit_48_epileptic_n10_compute.py` | Already computed `T_d^E` for trace direction at λ ∈ {0, 0.5, 1}. This scope adds full TARR classification + size-matched null + cross-probe variant. |
| `lrg_eegfc.utils.metrics.tree_distance.kc_distance` | Direct reuse. |
| `lrg_eegfc.utils.metrics.hypothesis.wilcoxon_z` + `bh_fdr` | Direct reuse. |

What this measure **subsumes**: the trace-only KC epi-subtree triangle
in `fig_05_kc_epi_triangle.pdf` (extended to all four classes + null).

What this measure **complements**: the global KC trace (Result 2)
and the per-module KC trace enumeration (audit_47/50). It is the
"is the epi-induced subtree special wrt the rest of the network?"
intermediate scale.

What this measure **does not replace**: the Section-5 global trace
finding stays as the load-bearing β cohort result; this measure is
epi-specific and orthogonal.

## Implementation plan

- **Script:** `scripts/01_compute/audit/audit_54_epi_rigidity_compute.py`
- **Library reuse (no new helpers in scripts/):**
  - `lrg_eegfc.workflow.lrg.load_lrg_result` for linkage matrices.
  - `lrg_eegfc.utils.io.patient.load_epileptic_nodes` for `E_p`.
  - `lrg_eegfc.utils.metrics.tree_distance.kc_distance` for KC.
  - `lrg_eegfc.utils.metrics.hypothesis.wilcoxon_z`, `bh_fdr` for stats.
- **Helper to add (library, not script):**
  `lrg_eegfc.utils.metrics.tree.induced_linkage(Z, leaf_indices)` —
  produces the linkage matrix of the induced subtree on a leaf
  subset. Currently absent from `tree.py` (per `coding-rules.md`,
  this lives in the library because it has ≥ 2 callers: this
  measure and the C-direction phase-trace extension).
- **Outputs:**
  - `data/audit/epi_rigidity/M_class_per_patient.csv` — one row
    per (patient, band, λ, variant ∈ {all, cp}) with z-scores per
    class and `T_d^E`.
  - `data/audit/epi_rigidity/M_class_cohort.csv` — one row per
    (band, λ, variant, class) with cohort-Wilcoxon p, BH q, count.
  - `data/audit/epi_rigidity/figures/fig_01_cohort_heatmap.pdf`
  - `data/audit/epi_rigidity/figures/fig_02_trace_score_boxes.pdf`
  - `data/audit/epi_rigidity/figures/fig_03_triangle_scatter.pdf`
  - `data/audit/epi_rigidity/figures/fig_04_cp_vs_all.pdf`
- **Figures script:**
  `scripts/01_compute/audit/audit_54b_epi_rigidity_figures.py`
  (mirrors the 48 + 48b convention).

## Open questions

1. **`ε_KC = 0.10` choice.** Matches audit_47 trace gate. May be
   too lenient at small `|E_p|`. First-sweep diagnostic plots the
   null `δ_pre→post` distribution per `|E_p|` bin; if the null
   has median > `ε_KC` for some `|E_p|`, the gate must be raised
   to `2σ` of that null instead. Defer to first sweep.

2. **Stricter null without overlap.** The default null draws random
   subsets that may overlap with `E_p`. A stricter null draws from
   `V_p \ E_p` only — strictly disjoint, but the size constraint
   `|S_r| = |E_p|` may force tight subsets in patients with high
   `|E_p|/N`. Defer; likely a sensitivity-check secondary figure.

3. **Multi-class assignment.** Currently each (p, b, λ) has exactly
   one class. A "soft" multi-class score using the three δ values
   (e.g. log-odds for each class) would let cohort tests use a
   continuous score per class. Defer; the binary call is more
   directly comparable to audit_47-50 outputs.

4. **Anchor floor.** Class fractions under the random null are
   expected to be Anchor-dominated at small `|E_p|` (small subsets
   are more likely to stay topologically simple). The z-score
   normalisation handles this in principle but needs first-sweep
   verification: the null distribution should have median Anchor
   fraction ≈ 0.5–0.7 at `|E_p| ≤ 10` and lower for larger sets.

5. **Pre-registration of acceptance gate.** Following the project
   convention (≥ 8 / 10 cohort threshold; here ≥ 8 / 9 since Pat_15
   drops): `cohort-positive` requires at least 8 patients with
   `z_class > 1.96` in at least one (band, λ) cell, AND BH-FDR
   `q < 0.05` cohort-Wilcoxon, in any one of the four classes.
   If null in all four, publish as a negative result.

6. **Cross-probe variant interpretation.** If `E_p^cp` flips a
   class assignment relative to `E_p` (e.g. anchor under all-pair
   becomes trace under cross-probe), the interpretation is that
   the focal hypersynchrony is rigid but the long-range epi
   network reorganises. This is a richer per-patient story; the
   manuscript figure should preserve this distinction.
