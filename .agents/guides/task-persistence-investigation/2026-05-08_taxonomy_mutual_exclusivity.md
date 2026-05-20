---
name: taxonomy-mutual-exclusivity
type: report
era: IMCOH_ABS / COHORT_N10
status: current
date: 2026-05-08
created: 2026-05-08
updated: 2026-05-08
scope: anchor-anatomy-baseline + mutual-exclusivity tie-breaker
pointers:
  - .agents/guides/task-persistence-investigation/2026-05-07_kc-anchor-modules.md
  - .agents/guides/task-persistence-investigation/2026-05-07_kc-trace-modules.md (audit_47 scope)
  - .agents/reports/2026-05-08_section5_6_verify.md
  - scripts/01_compute/audit/audit_60_anchor_anatomy_baseline.py
  - scripts/01_compute/audit/audit_60b_anchor_anatomy_baseline_figure.py
  - data/audit/anchor_anatomy_baseline/per_module_same_probe_fraction.csv
  - data/audit/anchor_anatomy_baseline/cohort_aggregate.csv
  - data/reports/notes_verification_2026-05-08/figures/anchor_anatomy_baseline.pdf
---

# §5.6 anchor anatomy baseline + mutual-exclusivity tie-breaker

**Head.** Anchor-class same-probe pair fraction sits **6.2x cohort baseline at
beta** (median 0.524 vs baseline 0.084, n=43 anchor modules across 10 patients)
and 6.9x / 11.9x at alpha / low-gamma — anchors are **not** anatomy-neutral
nulls; they are partially anatomy-anchored under |ImCoh|, and the §5.6 prose
should add a probe-geometry caveat to any "modules invariant across phases"
language. Trace at the same beta cell sits at 4.0x cohort baseline (median
0.333, n=25), so the trace class also carries a non-trivial probe-geometry
component but anchor enrichment dominates trace enrichment by ~1.6x at beta.
The mutual-exclusivity tie-breaker rule is formalized below; the canonical
class-ordering is **anchor > trace > reset > rearrange > diffuse** (tightest
gate first), with leaves classified at level X removed from the candidate
pool of all looser classes before the next pass runs.

---

## Notation

Inherits from `2026-05-07_kc-anchor-modules.md` and `audit_47`.

- `p ∈ COHORT_N10`, `b ∈ {δ, θ, α, β, γ_l, γ_h}`, `φ ∈ {pre, tt, post}`.
- `H_φ` — LRG dendrogram of phase `φ` (linkage matrix `Z_φ`).
- `S_φ` — set of internal subtrees at phase `φ` with leaf-set size in
  `[MIN_SIZE = 3, ⌊N_p / 2⌋]`.
- `J(A, B) = |A ∩ B| / |A ∪ B|` — Jaccard.
- `pid: {0, …, N_p − 1} → Σ` — the per-contact probe label, derived
  from the canonical channel labels CSV via
  `lrg_eegfc.utils.io.patient.parse_seeg_label`. Σ is the alphabet
  of probe identifiers; for sEEG these are the leading letters
  (e.g. `A`, `H'`, `OF`).
- `same_probe(i, j) = 1[pid(i) = pid(j)]`.
- For a leaf set `L ⊆ {0, …, N_p − 1}`:
  - `n_pairs(L) = binom(|L|, 2)`.
  - `n_same(L) = |{(i, j) : i < j ∈ L, pid(i) = pid(j)}|`.
  - `f_sp(L) = n_same(L) / n_pairs(L)` if `|L| ≥ 2`, undefined otherwise.
- Cohort baseline (per patient, anchored only by the implant geometry):

  ```
  f_sp^baseline(p) = ( Σ_q binom(n_q, 2) ) / binom(N_p, 2)
  ```

  where `n_q` is the number of contacts on probe `q ∈ Σ`. This is the
  expected `f_sp` of a random `|L|`-leaf draw at the limit `|L| → N_p`,
  and the natural anchor-geometry null at any leaf-set size.

## Definitions

### Per-module same-probe statistics

For a class `c ∈ {anchor, trace, reset, rearrange}` and a module
`m ∈ M_c(p, b)`:

| symbol | formula | semantics |
|--------|---------|-----------|
| `L^strict_c(m)` | per §5.6: `L_pre ∩ L_tt ∩ L_post` for anchor; `L_tt ∩ L_post` for trace; `L_pre ∩ L_post` for reset; `L_post` (raw, pre-residual) for rearrange | strict-intersection leaf set |
| `f_sp(m)` | `n_same(L^strict_c(m)) / n_pairs(L^strict_c(m))` | within-module same-probe pair fraction |
| `Δ(m, p)` | `f_sp(m) − f_sp^baseline(p)` | anatomy-corrected enrichment |
| `R(m, p)` | `f_sp(m) / f_sp^baseline(p)` | enrichment ratio (use only when `f_sp^baseline(p) > 0`) |

Cohort aggregate per `(b, c)`:
- `med_f_sp(b, c) = median over modules of f_sp(m)` — central tendency.
- `med_R(b, c) = med_f_sp(b, c) / med_baseline(b)` where
  `med_baseline(b)` is the median per-patient cohort baseline at band `b`
  (the baselines vary slightly across patients — Pat_03 has the smallest
  at 0.081; Pat_15, Pat_02 the largest at 0.091).

### Mutual-exclusivity tie-breaker

The four detectors (`audit_47`–`audit_50`) emit candidate modules with
overlapping strict leaves: a single contact may appear in (e.g.) both
an anchor module and a trace module. The §5.6 manuscript treats the
five classes as a **partition of the leaf set**, requiring an explicit
tie-breaker.

**Class-ordering (canonical).** From tightest gate to loosest:

1. **anchor** — `J(L_pre, L_tt) ≥ 0.6 ∧ J(L_pre, L_post) ≥ 0.6 ∧
   J(L_tt, L_post) ≥ 0.6` (three-way matching; gate strength ≈ 0.6³).
2. **trace** — `J(L_tt, L_post) ≥ 0.6 ∧ J(L_pre, ·) < 0.5 ∧ FRAG_FACTOR=2`
   containment guard against rPre.
3. **reset** — symmetric to trace, `J(L_pre, L_post) ≥ 0.6 ∧ J(L_tt, ·) < 0.5`.
4. **rearrange** — rPost-anchored only, both rPre and tt fragmented
   (`J(L_post, L_pre) < 0.5`, `J(L_post, L_tt) < 0.5`, `FRAG_FACTOR=2`
   on both sides).
5. **diffuse** — residual leaves not classified by any of the above.

**Predicate.** For each leaf `i ∈ {0, …, N_p − 1}`, define
`cls(i) ∈ {anchor, trace, reset, rearrange, diffuse}` recursively:

- If there exists `m ∈ M_anchor` with `i ∈ L^strict_anchor(m)` and the
  residual `L^strict_anchor(m) \ used` still has size ≥ 3 after
  subtracting all earlier-classified leaves: `cls(i) ← anchor`.
- Else if there exists `m ∈ M_trace` analogously: `cls(i) ← trace`.
- Else if there exists `m ∈ M_reset` analogously: `cls(i) ← reset`.
- Else if there exists `m ∈ M_rearrange` analogously: `cls(i) ← rearrange`.
- Else: `cls(i) ← diffuse`.

The "size ≥ 3" floor mirrors the existing `assemble_class_modules`
(`audit_52:188–198`) rearrange residual rule. Modules whose residual
falls below the floor are dropped, and their leaves fall through to
the next class — they do not migrate elsewhere directly.

## Properties

### What the measure detects

- `f_sp(m)` is `[0, 1]`-valued, monotone in same-probe pair count for
  fixed `|L|`. `f_sp = 1` ⇔ all leaves on a single probe; `f_sp = 0` ⇔
  every pair crosses probes.
- `R(m, p) = 1.0` is the anatomy-neutral null: a module whose leaves
  are drawn uniformly from the implant. `R(m, p) > 1` flags
  anatomy-enrichment; `R(m, p) < 1` flags anti-anatomy structure
  (modules deliberately spanning probes).
- The mutual-exclusivity tie-breaker partitions `{0, …, N_p − 1}`
  cleanly: every leaf is assigned exactly one class, and the strict
  class-leaf sets are pairwise disjoint by construction.

### What the measure does NOT detect

- `f_sp(m)` is silent on the **spatial extent** of probes. Two probes
  3 mm apart count as cross-probe; two contacts on the same probe
  10 mm apart count as same-probe. It is a coarse anatomy proxy.
- It is silent on **anatomy other than implant geometry**. Two contacts
  in different cytoarchitectonic areas but on the same straight insertion
  trajectory will count as same-probe; two contacts in the same gyrus on
  different trajectories will count as cross-probe.
- The cohort baseline `f_sp^baseline(p)` is the null expected for a
  random leaf draw at the limit `|L| → N_p`. For finite `|L|` (typical
  module sizes 3–10), the random-draw null fluctuates above and below
  `f_sp^baseline(p)`. Modules at the lower-quartile boundary should be
  read against a `|L|`-matched bootstrap null, not the cohort scalar.
  This is documented as an open question, not implemented in the
  first pass.

### Range, identifiability, complexity

- `f_sp ∈ [0, 1]`; `R ∈ [0, 1/f_sp^baseline(p)]` (typically `[0, 12]`
  given baselines `≈ 0.08`).
- The tie-breaker is unique up to ordering of modules within a class
  (it is not order-independent across modules of the same class — first
  encountered wins). With deterministic input (sorted candidate list
  per class detector), the output is reproducible.
- Per-module: `O(|L|²)` for `n_same(L)`. Per-cell with all four classes:
  `O((MAX_TOTAL_per_class)² × Σ_classes)` — bounded by 4 × 12² = 576
  pair lookups per cell, fully tractable.
- Tie-breaker complexity: `O(N_p × 5)` — a single linear pass per leaf
  per class.

## Caveats & failure modes

| | Caveat | Mitigation |
|--|--|--|
| 1 | Probes vary in number of contacts (4–18 typical). A single 18-contact probe inflates `f_sp^baseline(p)` substantially. | Report per-patient baseline (variation 0.072–0.091 across cohort). |
| 2 | Pat_03 outlier (1024 Hz) has the smallest baseline (0.081) — the smallest implants. | Pat_03 included unflagged; the per-patient baseline correction handles the geometry, the band-resolved figure shows Pat_03 in the same boxplot. |
| 3 | Small modules (`|L| = 3`) produce `f_sp ∈ {0, 1/3, 2/3, 1}` — only four values, dominated by single-pair coincidence. | Cohort aggregate is over modules; the median across n=43 β anchors is robust. Per-module noise is high. |
| 4 | Anchor modules may overlap on contacts in their **leaf-set union** but NOT in their strict three-way intersection. The cohort aggregate is computed on the strict set. | Document choice; both options reasonable. Strict-set is consistent with §5.6 prose. |
| 5 | The tie-breaker is order-dependent: a leaf classified as anchor cannot later be reclassified to trace. The intent (anchor > trace) is correct given gate strength, but the implementation must commit early. | Class-ordering is registered; do not change without a sensitivity audit. |
| 6 | `MIN_RESIDUAL = 3` may drop borderline modules whose strict leaves are all consumed by a higher-priority class. | Document; report dropped count per cell as a sensitivity column in the CSV (added to `n_residual_leaves_after_tiebreak`). |

## Pseudocode — first-pass anchor anatomy baseline + tie-breaker

```
input:  Z_pre, Z_tt, Z_post, N_p, probe_ids[0..N_p-1]
output: per-module records, aggregate same-probe fractions

# 1. Cohort baseline geometry
counts_per_probe ← histogram(probe_ids)
n_same_baseline  ← Σ_q binom(counts_per_probe[q], 2)
f_sp_baseline    ← n_same_baseline / binom(N_p, 2)

# 2. Run the four canonical class detectors
M_anchor    ← find_anchor_triples(Z_pre, Z_tt, Z_post, N_p)[0.0]
M_trace     ← find_trace_pairs   (Z_tt, Z_post, Z_pre, N_p)[0.0]
M_reset     ← find_reset_pairs   (Z_pre, Z_post, Z_tt, N_p)[0.0]
M_rearrange ← find_rearrangement_modules(Z_pre, Z_tt, Z_post, N_p)

# 3. Mutual-exclusivity tie-breaker (canonical class-ordering)
used   ← ∅
T      ← {anchor: [], trace: [], reset: [], rearrange: [], diffuse: []}
for cls in [anchor, trace, reset, rearrange]:
    for m ∈ M_cls:
        L_strict ← strict_leaves_for_class(cls, m)
        L_resid  ← L_strict \ used
        if |L_resid| ≥ MIN_RESIDUAL = 3:
            T[cls] ← append (m, L_resid)
            used   ← used ∪ L_resid
diffuse_leaves ← {0..N_p-1} \ used
T[diffuse] ← [(none, diffuse_leaves)] if |diffuse_leaves| ≥ 3 else []

# 4. Per-module same-probe statistics (using strict leaves, NOT residuals)
records ← []
for cls in [anchor, trace, reset, rearrange]:
    for idx, m ∈ enumerate(M_cls):
        L      ← strict_leaves_for_class(cls, m)
        n_pair ← binom(|L|, 2)
        n_sp   ← |{(i,j): i<j ∈ L, pid(i) = pid(j)}|
        records ← append (cls, idx, |L|, n_pair, n_sp,
                          n_sp / n_pair if n_pair > 0 else NaN,
                          f_sp_baseline)

# 5. Aggregate per (band, class)
aggregate by (band, class) over records:
    n_modules, median(f_sp), mean(f_sp), q25(f_sp), q75(f_sp),
    median(|L|), median(f_sp_baseline)
```

## First-pass cohort numbers (locked 2026-05-08)

`data/audit/anchor_anatomy_baseline/cohort_aggregate.csv`. Median
same-probe-pair fraction at strict-intersection leaves vs cohort
baseline (median across patients):

| band | class | n_modules | median f_sp | baseline | enrichment R |
|------|-------|-----------|-------------|----------|--------------|
| **β** | **anchor** | **43** | **0.524** | 0.084 | **6.22x** |
| β | trace     | 25 | 0.333 | 0.084 | 3.96x |
| β | reset     | 18 | 0.348 | 0.081 | 4.28x |
| β | rearrange | 87 | 0.184 | 0.081 | 2.25x |
| α | anchor    | 16 | 0.548 | 0.080 | 6.87x |
| α | trace     | 23 | 0.288 | 0.081 | 3.54x |
| α | reset     | 13 | 0.333 | 0.087 | 3.85x |
| α | rearrange | 100 | 0.167 | 0.083 | 2.01x |
| δ | anchor    | 11 | 0.333 | 0.078 | 4.27x |
| δ | trace     | 17 | 0.308 | 0.084 | 3.65x |
| δ | reset     | 10 | 0.333 | 0.082 | 4.08x |
| δ | rearrange | 100 | 0.167 | 0.083 | 2.01x |
| θ | anchor    | 12 | 0.317 | 0.083 | 3.82x |
| θ | trace     | 16 | 0.322 | 0.077 | 4.19x |
| θ | reset     | 15 | 0.333 | 0.081 | 4.09x |
| θ | rearrange | 103 | 0.167 | 0.081 | 2.05x |
| γ_l | anchor   | 45 | 1.000 | 0.084 | 11.87x |
| γ_l | trace    | 23 | 0.400 | 0.078 | 5.12x |
| γ_l | reset    | 26 | 0.348 | 0.085 | 4.11x |
| γ_l | rearrange | 81 | 0.214 | 0.081 | 2.63x |
| γ_h | anchor   | 17 | 1.000 | 0.084 | 11.87x |
| γ_h | trace    | 11 | 0.333 | 0.081 | 4.09x |
| γ_h | reset    | 1  | 0.333 | 0.091 | 3.68x |
| γ_h | rearrange | 89 | 0.000 | 0.081 | 0.00x |

**Headline reading.**

1. **Anchors are NOT anatomy-neutral.** At every band, anchor enrichment
   is `≥ 3.8x` and at low-/high-gamma the median anchor module sits at
   `f_sp = 1.0` — entirely on a single probe. The §5.6 prose framing
   "anchors are the cross-phase null" needs a probe-geometry caveat:
   anchors are the *functional-task null*, but they are partly the
   *anatomy attractor* under |ImCoh|.
2. **Trace enrichment is real but smaller than anchor.** β trace sits
   at 3.96x baseline, beneath β anchor's 6.22x. Trace is anatomy-tinted
   but more cross-probe than anchor — consistent with the trace class
   capturing a genuine functional reorganization signal that crosses
   probe boundaries more often than anchors do.
3. **Rearrange is the anatomy-neutral class.** β rearrange enrichment
   is 2.25x; γ_h rearrange median sits at 0 (more pairs cross probes
   than within). Rearrange is the cleanest "spans the implant" class.
4. **Gamma anchors are extreme** — `f_sp = 1.0` median means that more
   than half of low-gamma anchor modules are entirely on one probe,
   driven by residual same-probe broadband coupling that |ImCoh|
   suppresses but does not eliminate. The gamma-band anchor result
   should be reported with the strongest anatomy caveat.

## Visualization spec

`data/reports/notes_verification_2026-05-08/figures/anchor_anatomy_baseline.pdf`,
2-row layout, vector PDF, no PNG sibling.

- **Top row:** boxplot per (band × class), 6 bands × 4 classes = 24
  boxes. Class colors: anchor purple `#807dba`, trace blue `#2c7fb8`,
  reset green `#41ab5d`, rearrange orange `#fc9272`. Cohort baseline
  drawn per-band as a thin gray dotted line at the cohort-median
  `f_sp^baseline(p)` plus a translucent band spanning the inter-patient
  min-max range.
- **Bottom row:** β-only per-patient scatter. x = size-weighted-mean
  anchor `f_sp`; y = size-weighted-mean trace `f_sp`. Diagonal `y = x`
  drawn as a thin reference. Cohort baseline cross drawn at the median
  baseline. Markers by `tab10`; per-patient label annotations offset
  to avoid stacking.

**Reading rules.**

- Top panel: anchor boxes well above the dotted line at every band ⇒
  anchors are anatomy-tinted. Trace and reset boxes only modestly
  above the dotted line. Rearrange boxes touch the line at γ_h.
- Bottom panel: most patients **below** y = x ⇒ at β, anchors are
  more probe-localized than traces in the same patient. Patients
  above the diagonal (Pat_10, Pat_13) have very few β trace modules
  (n=1, n=2 respectively) so the comparison is noisy.

## Connection to prior tools

| Prior tool | Relation |
|--|--|
| `audit_50_kc_anchor_module_view.py` | Defines the anchor module set this audit consumes; we add the per-module same-probe statistic and the per-patient baseline. |
| `audit_47_kc_trace_network_view.py` | Defines the trace module set; same role, different class. |
| `audit_48_kc_reset_module_view.py` | Defines the reset module set. |
| `audit_49_kc_rearrangement_module_view.py` | Defines the rearrange candidate set; we use the raw `leaves_post`, not the audit_52 residual. |
| `audit_52_section5_class_showcases.py` | Implements the §5.6 strict-intersection definitions and the rearrange-residual subtraction; **this** audit adopts the same strict definitions but **does not** apply the rearrange-residual rule (we want the raw rearrange anatomy signal). |
| `audit_15_anatomy_mspc.py` | Anatomy-MSPC at the partition level; not directly comparable (partitions ≠ modules) but the same Desikan-Killany / `(x, y, z)` axis matters for follow-up beyond probe geometry. |
| Probe-bias guide (`.agents/guides/02_methods/probe-bias-guide.md`) | The MSC-era probe-bias mitigation; under |ImCoh| the issue is *residual* anatomy, not the dominant zero-phase-lag bias. This audit quantifies how much residual remains. |
| §5.6 prose | Adds the "anchor anatomy caveat" requested by the question — should appear as a sentence in the §5.6 paragraph that introduces the anchor class. |

## Open questions

1. **Size-matched bootstrap null instead of cohort scalar baseline.** A
   3-leaf module's expected `f_sp` under random draw isn't exactly
   `f_sp^baseline(p)` — it's an integer ratio out of `binom(3, 2) = 3`
   trials, so the null expectation is `f_sp^baseline(p)` but the null
   variance is large. A `|L|`-matched bootstrap null (R=1000 random
   leaf draws of the right size, per patient) would give a per-module
   z-score. Out of scope for first pass; document as the natural V2.

2. **Probe-aware FC weighting at LRG construction.** A more principled
   anatomy mitigation than reporting the residual: zero (or down-weight)
   same-probe edges in `A` before the LRG eigendecomposition, then
   re-run `audit_47`–`50`. This was the MSC-era plan and was retired
   under |ImCoh| (memory `probe_bias_critical.md`). The β-anchor 6.2x
   number suggests the retirement was premature — same-probe `|ImCoh|`
   isn't trivially zero and may need attention. Track as a follow-up
   if §5.6 prose requires it.

3. **Anatomy beyond probe geometry.** Pair-level Desikan-Killany region
   labels (`implant_pat_NN.csv`) would let us replace `pid(i) = pid(j)`
   with `dk(i) = dk(j)`, and `(x, y, z)` Euclidean distance would let us
   replace the binary same-probe with a continuous proximity score.
   These are the canonical anatomy axes (memory
   `feedback_implant_anatomy_not_letters.md`). Out of scope for the
   first-pass anchor caveat but the right next step.

4. **Tie-breaker order sensitivity.** The canonical order is
   anchor > trace > reset > rearrange. Alternative orders
   (trace > anchor > … or strict-Jaccard-rank ordering) would change
   the residual leaf counts. A sensitivity sweep is straightforward:
   run the tie-breaker under the 4! = 24 class permutations; report
   the leaf-classification confusion matrix. Defer to §5.6 polish.

5. **Cohort-statistics vs per-cell statistics.** The aggregate above
   pools modules across patients and bands. A per-cell `R(p, b, c)`
   matrix would let us check whether anchor enrichment is band-specific
   or patient-specific. The per-module CSV carries the data; the
   matrix view is a follow-up figure.
