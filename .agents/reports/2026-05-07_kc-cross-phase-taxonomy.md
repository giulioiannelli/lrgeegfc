---
name: KC cross-phase module taxonomy
description: Four-class cohort tally of cross-phase module behavior (trace / reset / rearrange / anchor) on FC dendrograms — ties the §5 KC scalar trace to a per-class structural inventory
type: result
era: IMCOH_ABS × COHORT_N10
status: locked_2026-05-07
created: 2026-05-07
updated: 2026-05-07
scripts:
  - scripts/01_compute/audit/audit_47_kc_trace_network_view.py
  - scripts/01_compute/audit/audit_48_kc_reset_module_view.py
  - scripts/01_compute/audit/audit_49_kc_rearrangement_module_view.py
  - scripts/01_compute/audit/audit_50_kc_anchor_module_view.py
related:
  - .agents/reports/2026-05-07_kc-trace-module-visualization.md
  - .agents/guides/01_project/terminology.md
  - .agents/guides/task-persistence-investigation/2026-05-07_kc-reset-modules.md
  - .agents/guides/task-persistence-investigation/2026-05-07_kc-rearrangement-modules.md
  - .agents/guides/task-persistence-investigation/2026-05-07_kc-anchor-modules.md
  - data/audit/kc_trace_network_view/trace_subtrees_summary.csv
  - data/audit/kc_reset_module_view/reset_subtrees_summary.csv
  - data/audit/kc_rearrangement_module_view/rearrangement_subtrees_summary.csv
  - data/audit/kc_anchor_module_view/anchor_subtrees_summary.csv
---

# Four-class cross-phase module taxonomy on FC dendrograms

**Head.** Across the COHORT_N10 era (n=10 patients × 6 bands) we now
have a complete inventory of how cross-phase modules behave:
**trace** (107 modules at n=9; Pat_08 excluded as topology-blind by
construction), **reset** (83 modules at n=10), **rearrange** (560 at
n=10), **anchor** (144 at n=10). Three points stand out. First,
**rearrangement is the abundant default**: 560 modules in 60/60
cells, ~9.3 per cell — most rPost coherence is emergent, not
task-seeded. Second, **β co-locates trace and anchor** (23 trace +
43 anchor at n=10, the latter universal at 10/10 cohort cells) — the
load-bearing trace band rides on top of a *larger* phase-invariant
backbone, so the §5.2 signal is a small reorganization on a
structurally conserved hierarchy, not a wholesale rewrite. Third,
**Pat_10 inverts the trace↔reset ratio** (3 trace vs 17 reset; T:R =
0.18) and **Pat_08 has the most extreme inversion** (T:R = 0/7,
undefined; he is the height-only-trace patient by construction —
zero topology trace, but 7 reset + 50 rearrange + 10 anchor).

## Cohort scope: n=10 sibling, n=9 trace

The four audits run on different patient subsets, by design:

- **trace (audit_47)**: n=9 = COHORT_N10 \ Pat_08. Pat_08 is the
  canonical *height-only-trace* counterexample — 5/5 cohort scalar
  trace at β yet zero topology-trace modules at any band under the
  J(tt, rPre) < 0.5 fragmentation gate. Including him in the trace
  audit would only render 6 empty cells. The trace report
  (`.agents/reports/2026-05-07_kc-trace-module-visualization.md`)
  documents this exclusion as load-bearing for the trace ↔ height
  decomposition of §5.2.
- **reset (48), rearrange (49), anchor (50)**: n=10 = full
  COHORT_N10. There is no a priori reason Pat_08 lacks reset,
  rearrangement, or anchor modules; his trace blind-spot is
  topology-driven (J(tt, rPre) cannot fall below 0.5 in his β
  hierarchy), and the sibling audits use *different* anchor phases
  (rPre / rPost / triple-match) where the gate doesn't apply.

The n=10 sibling counts confirm the prediction: Pat_08 has 67
non-trace modules across the three sibling audits (7 reset + 50
rearrange + 10 anchor), with low_γ anchor showing 6 modules
(including a 38-leaf one) — substantial stable topology consistent
with his "trace via height shifts on a stable backbone" mechanism.

## The four classes

Per the project's [terminology
guide](../guides/01_project/terminology.md), modules across phases
{rest_pre, task_test, rest_post} are partitioned into four classes
by a Jaccard-fragmentation contract on subtree leaf sets. All four
audits use the same gates: Jaccard ≥ 0.6 = "coherent", < 0.5 =
"fragmented", with `FRAG_FACTOR = 2` containment guard against false
matches by inclusion.

| class | rPre | task_test | rPost | anchor / script | n_modules | cells / max |
|:--|:--:|:--:|:--:|:--|:--:|:--:|
| **trace** | fragmented | coherent | coherent | tt → rPost match (audit_47) | **107** | 45/54 (n=9) |
| **reset** | coherent | fragmented | coherent | rPre → rPost match (audit_48) | **83** | 39/60 (n=10) |
| **rearrange** | fragmented | fragmented | coherent | rPost-only emergent (audit_49) | **560** | 60/60 (n=10) |
| **anchor** | coherent | coherent | coherent | triple-match (audit_50) | **144** | 50/60 (n=10) |

cells / max = #(patient, band) cells with at least one module at λ=0
(topology) over total cells. Rearrangement uses `n_modules` (size ≥
3 stratum); the others use `n_modules_lambda0`.

## Per-band module counts at λ = 0

| | δ | θ | α | β | low_γ | high_γ | total |
|:--|--:|--:|--:|--:|--:|--:|--:|
| trace (n=9) | 15 | 15 | 23 | **23** | 21 | 10 | **107** |
| reset (n=10) | 10 | 15 | 13 | 18 | **26** | 1 | **83** |
| rearrange (n=10) | **100** | **103** | **100** | 87 | 81 | 89 | **560** |
| anchor (n=10) | 11 | 12 | 16 | 43 | **45** | 17 | **144** |

Bold = per-class band-leader. Cells per band:

| | δ | θ | α | β | low_γ | high_γ |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|
| trace (/ 9) | 7 | 7 | 8 | **9** | 8 | 6 |
| reset (/ 10) | 6 | 8 | 6 | 8 | **10** | 1 |
| rearrange (/ 10) | **10** | **10** | **10** | **10** | **10** | **10** |
| anchor (/ 10) | 7 | 7 | 9 | **10** | **10** | 7 |

**β is the only band where three classes (trace, rearrange, anchor)
are universal across the cohort** — trace 9/9, rearrange 10/10,
anchor 10/10 (reset 8/10 — high but not universal). This is the
structural fingerprint of the §5.2 load-bearing β trace: the
topology has a stable conserved backbone (anchor 10/10), a wide
reservoir of emergent rPost coherence (rearrange 10/10), AND a
discrete population of task-seeded modules (trace 9/9 modulo Pat_08
height-only). β is special not because it lacks competitors but
because it carries all of them.

**high_γ has 1 reset module across the entire cohort** — Pat_15
alone, even with Pat_08 added. γ-band rPre topology essentially
never rebuilds itself after a tt disruption. Consistent with the
§5.2 result that γ_h is ergodic on KC scalars: at fast scales the
rest_pre arrangement does not return.

## Per-patient cohort sums

Across all 6 bands, by class:

| patient | trace | reset | rearrange | anchor | T : R |
|:--|--:|--:|--:|--:|:--:|
| Pat_02 | 15 | 7 | 57 | 24 | 2.14 |
| Pat_03 | 12 | 8 | 65 | 9 | 1.50 |
| Pat_05 | 11 | 2 | 54 | 14 | **5.50** |
| Pat_06 | 16 | 5 | 70 | 5 | 3.20 |
| Pat_07 | 21 | 5 | 44 | 19 | 4.20 |
| **Pat_08** | **0** | **7** | **50** | **10** | **0.00** |
| Pat_10 | 3 | 17 | 47 | 16 | **0.18** |
| Pat_13 | 14 | 16 | 52 | 24 | 0.88 |
| Pat_14 | 11 | 8 | 56 | 16 | 1.38 |
| Pat_15 | 4 | 8 | 65 | 7 | 0.50 |
| **n=10** | 107 | 83 | 560 | 144 | **1.29** |

(trace column for Pat_08 is structurally 0, not measured — see scope
caveat above.)

**Three extremes:**

- **Pat_05** — most trace-heavy by ratio (T:R = 5.5; 11 trace, 2
  reset). His β trace is a single 47-leaf cluster covering ~40 % of
  the network (discovered when J_POST_MIN was lowered to 0.6).
- **Pat_07** — richest absolute trace count (21 across all 6 bands,
  including 6 at β). Manuscript flagship.
- **Pat_08 / Pat_10** — opposite-tail reset specialists. Pat_08
  carries 0 trace by construction (height-only) but has 7 reset.
  Pat_10 carries the cohort's largest reset count (17 modules; 22 %
  of cohort total at n=9, 20 % at n=10).

## Three observations

### 1. Rearrangement is the abundant default, not a special pattern

560 rearrangement modules vs 107 trace + 83 reset + 144 anchor =
334 matched-across-phase. **Rearrangement is ~1.7 × larger than all
three matched classes combined.** Cohort-wide every cell (60 / 60)
is rearrangement-positive, with mean ~9.3 modules per cell.

This reframes the trace question. The structural baseline is *not*
"a stable resting hierarchy that task occasionally perturbs". It is
"each phase's hierarchy is mostly its own, with a small population
of modules that match across phases". The trace pattern is the
narrow sub-class where the rPost match is **task-seeded** (the same
leaves were already coherent in tt). The reset pattern is the
narrow sub-class where the rPost match is **rPre-seeded** (the
leaves were coherent before task, broken by task, restored after).
The anchor is the narrowest sub-class where the match is
**phase-invariant**.

A useful sanity check is that trace + reset + anchor = 334 ≈ 60 % of
rearrange. Most rPost subtrees do not have a sibling in tt or rPre.

### 2. β co-locates trace and anchor (the load-bearing band carries both)

β is the only band where trace, rearrange, and anchor are each
maximally cohort-wide (9/9, 10/10, 10/10 — reset is 8/10 at β, high
but not universal). At β the cohort totals are:

- 23 trace modules (the §5.2 signal — task seeds new rPost coherence),
- 43 anchor modules (the largest non-low_γ anchor band-count;
  phase-invariant backbone),
- 18 reset modules (rPre topology that survives task disruption),
- 87 rearrange modules (emergent rPost coherence).

**The §5.2 trace lives on top of a substantially larger anchor
backbone** — the cohort has 43 phase-invariant β subtrees vs 23
task-seeded ones. This is a more honest statement of the §5.2
mechanism than "task reorganizes β topology": task **adds** a
specific population of new modules to a backbone that mostly does
not change. The KC scalar trace at λ = 0 is the difference of two
small numbers (modules that appear in tt+rPost but not rPre) on top
of a much larger constant (modules that appear in all three phases).

Adding Pat_08 (n=9 → n=10 sibling) bumped β anchor from 42 to 43
(only 1 anchor module: a 10-leaf cluster) while contributing 12
rearrangement and 2 reset modules at β — consistent with the
height-only-trace reading: his topology is largely stable in the
small region that matters for his scalar trace, with most cross-
phase variation expressed as rearrangement.

### 3. Pat_10 inverts the trace↔reset ratio; Pat_08 maxes the inversion

Pat_10 T:R = 0.18, an order of magnitude below the cohort median
(1.29). Pat_08 T:R is undefined (0 trace by construction, 7 reset)
— the most extreme reset specialist by ratio, but for a different
reason than Pat_10: Pat_10 has measurable trace and just very large
reset; Pat_08 has *no* topology trace at all.

Both are corroborating evidence for cohort heterogeneity in the
trace mechanism:

- Pat_10's task does disrupt resting topology (high reset — 17
  modules) but does not seed new persistent modules (low trace —
  3). The four-class inventory makes the structural origin of
  Pat_10's H2a / H4 ergodic-θ outlier behavior visible.
- Pat_08's task either does not reorganize topology or reorganizes
  it in a way that does not survive J(tt, rPre) < 0.5 — but the §5.2
  scalar trace fires at β anyway, because his β trace is carried by
  *merge-height shifts within stable topology* (the λ=1 KC component).

The cohort range of T:R from undefined / 0.18 (Pat_08, Pat_10) to
5.50 (Pat_05) covers nearly two orders of magnitude. The §5.2
scalar trace tests for the cohort-aggregate signal but the per-
patient mechanism splits across at least these two regimes
(trace-rich vs reset-rich, with a height-only third regime for
Pat_08).

### 4. high_γ resets do not exist (still)

high_γ has 1 reset module across the entire n=10 cohort (Pat_15
alone). Pair this with high_γ trace at 10 modules (the band's
lowest count) and rearrange at 89 (the third-highest count). high_γ
topology behaves as if each phase rebuilds a fresh hierarchy that
has almost no leaf overlap with the others, except for whatever is
in the anchor (17 modules, also low). **At γ_h, "the previous
resting hierarchy" effectively does not propagate** — rPre never
resurfaces in rPost (reset = 1) AND task seldom seeds it (trace =
10). This is the structural correlate of γ_h's ergodicity on the KC
scalar trace.

Adding Pat_08 to the n=10 cohort did NOT change this: Pat_08 has 0
high_γ reset and 4 high_γ rearrange modules.

## Why the four-class table matters for the manuscript

The §5.2 KC scalar trace at λ = 0 is a single number per
(patient, band) that is significant cohort-wide for β. The
four-class taxonomy unpacks that number into a structural inventory:

- Trace modules are **what the scalar tests for** — discrete
  identifiable subtrees that fragment in rPre and crystallize in
  tt+rPost. The cohort has 107 of them at λ = 0 across n=9; 23 are
  at β; the cohort-level signal at β is statistically detectable.
- Anchor modules are **what the scalar does NOT test for** — but
  they are 1.9 × more abundant at β (43 vs 23 at n=10). A β-trace
  manuscript that omits anchor modules implicitly suggests the
  hierarchy is mobile under task. The data say it mostly is not;
  task adds a small task-seeded layer on top of a larger phase-
  invariant backbone.
- Reset modules are **the inverse of trace** (rPre coherent, tt
  fragments, rPost recovers). They exist (83 modules at n=10), are
  patient-specific (Pat_10 carries 20 % of cohort total), and
  cluster at low_γ rather than β. They give a vocabulary for
  Pat_10's outlier behavior without resorting to "noise".
- Rearrangement modules are **the abundant default** that calibrates
  the others. The fact that rearrangement is ~3.8 × more common
  than trace at β (87 vs 23) tells us most rPost coherence at β is
  NOT task-seeded. The trace signal is the small task-seeded slice.

A future writeup of §5 should at minimum cite the trace + anchor +
rearrangement table for β to defend against a "the brain mostly
reorganizes" reading and to localize what the trace actually is.

## Pointers

- Trace report (canonical): `.agents/reports/2026-05-07_kc-trace-module-visualization.md`
- Trace scope (Module Retention Landscape): `.agents/guides/task-persistence-investigation/2026-04-25_module-retention-landscape.md`
- Reset scope: `.agents/guides/task-persistence-investigation/2026-05-07_kc-reset-modules.md`
- Rearrangement scope: `.agents/guides/task-persistence-investigation/2026-05-07_kc-rearrangement-modules.md`
- Anchor scope: `.agents/guides/task-persistence-investigation/2026-05-07_kc-anchor-modules.md`
- Terminology guide: `.agents/guides/01_project/terminology.md`
- §5.2 KC scalar trace: `.agents/reports/2026-05-05_result-2-lrg-beta-trace.md`

### Audit scripts and outputs

| class | script | figure dir | summary CSV | rows |
|:--|:--|:--|:--|:--:|
| trace | `audit_47_kc_trace_network_view.py` | `data/audit/kc_trace_network_view/figures/` | `trace_subtrees_summary.csv` | 9×6=54 |
| reset | `audit_48_kc_reset_module_view.py` | `data/audit/kc_reset_module_view/figures/` | `reset_subtrees_summary.csv` | 10×6=60 |
| rearrange | `audit_49_kc_rearrangement_module_view.py` | `data/audit/kc_rearrangement_module_view/figures/` | `rearrangement_subtrees_summary.csv` | 10×6=60 |
| anchor | `audit_50_kc_anchor_module_view.py` | `data/audit/kc_anchor_module_view/figures/` | `anchor_subtrees_summary.csv` | 10×6=60 |

Each script renders 3 × 3 grids (rows = λ regime ∈ {0, 0.5, 1} for
trace/reset/anchor, or size strata ∈ {≥3, ≥5, ≥8} for rearrangement
where λ-blend does not apply; cols = phase). Leaf-coloring rule
(2026-05-07): colored vertical bars at each leaf x-position extend
from y = 0 up to the leaf's parent merge height, overriding scipy's
default leaf-link gray. The earlier short-bar rendering (4 % of
y-range) is retired across all four scripts. Empty cells skipped.

## Open

1. **Library promotion.** audit_47/48/49/50 share ~60 % of helper
   scaffolding (linkage compute, module module-coloring, leaf-bar
   rendering, network drawing). Per the second-use library promotion
   rule these helpers should land in `lrg_eegfc.visuals.kc_module_view`
   in a follow-up commit.
2. **Anatomical mapping.** Cross-class module leaves should be mapped
   to (x, y, z) + Desikan-Killany regions via `implant_pat_NN.csv`
   to test whether trace modules at β concentrate on
   Hippocampus + left fusiform (the §5.2 anatomical claim) and
   whether anchor modules localize to a different anatomy.
3. **Cross-class summary figure.** A per-band stacked bar of T / R /
   Re / A counts — or a per-(patient × band) grouped bar chart —
   would surface the cohort distribution at a glance. Defer until
   library promotion lands.
