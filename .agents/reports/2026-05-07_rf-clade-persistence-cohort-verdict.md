---
name: 2026-05-07_rf-clade-persistence-cohort-verdict
type: report
era: IMCOH_ABS x COHORT_N10
status: parked
created: 2026-05-07
updated: 2026-05-07
pointers:
  - .agents/guides/task-persistence-investigation/2026-05-06_rf-k-multiscale-measure.md
  - .agents/guides/task-persistence-investigation/2026-05-04_lrg-trace-investigation-plan.md
  - .agents/guides/task-persistence-investigation/2026-05-07_kc-anchor-modules.md
  - .agents/guides/task-persistence-investigation/2026-05-07_kc-rearrangement-modules.md
  - .agents/guides/task-persistence-investigation/2026-05-07_kc-reset-modules.md
  - data/reports/section_5_lrg_trace/14_rf_clade_persistence/
  - data/audit/kc_trace_network_view/
  - data/audit/kc_reset_module_view/
  - data/audit/kc_rearrangement_module_view/
  - data/audit/kc_anchor_module_view/
  - scripts/01_compute/audit/audit_48_rf_k_clade_persistence.py
  - scripts/01_compute/audit/audit_48b_rf_k_soft_jaccard.py
  - scripts/01_compute/audit/audit_48c_rf_real_vs_null_scatter.py
  - scripts/01_compute/audit/audit_48d_trace_clades_on_dendrograms.py
---

# RF(k) clade-persistence — cohort numbers (PARKED, not in manuscript)

**The RF clade-persistence probe (hard threshold + soft mean-Jaccard +
per-clade strict 4-mode taxonomy) is set aside as a Section-5 measure.
Numbers are computed and documented below for the record; the probe
is NOT included in the manuscript. KC λ=0 (per-pair MRCA, audit_36 +
audit_46) and the audit_47-50 KC module-view family (trace / reset /
rearrangement / anchor) are the load-bearing tools for β's graded
clade reorganization. RF's Jaccard-based formulation has structural
problems with size-sensitivity and threshold-binarity that prevent it
from cleanly capturing graded phenomena. We keep the RF output on disk
as a strict-identity sensitivity check (16/60 cells across cohort
satisfy the strict emergent-then-persists predicate) in case it becomes
useful later, but it is parked, not load-bearing.**

## Why parked

RF reads clade leafsets via Jaccard. Jaccard is symmetric and size-
sensitive: a 10-leaf task clade matched against a 15-leaf rs_pre
containing subtree gives J = 10/15 = 0.67 (looks like presence) while
the same 10-leaf clade matched against a 5-leaf rs_pre subset gives
J = 5/10 = 0.50 (looks like partial presence). Neither value cleanly
distinguishes "leaves are tightly nested in a larger module" (which
is graded reorganization, not absence) from "leaves are partially
clustered" (which is genuine partial absence).

Soft RF (mean best-Jaccard, no threshold) inherits the same shape
problem — averaging a metric whose values are misleading does not
fix the underlying issue. Hard RF (Jaccard ≥ 0.70) compounds the
problem by binarizing it.

KC operates one geometric layer below: pair-level MRCA depths. Each
leaf-pair's MRCA depth is an integer count of internal nodes, not a
size-sensitive set comparison. KC cumulates graded pair-MRCA shifts
robustly across thousands of pairs. β at KC λ=0 is 10/10 vs null with
p = 0.000977 (audit_46); soft RF at β is 9/10 with p = 0.0068. KC was
already cleaner; soft RF was a redundant confirmation, not a new
finding.

The audit_47-50 KC module-view family (Section-5 4-mode taxonomy) is
ALSO Jaccard-based but uses size-matched best-J + containment
fragmentation factor — gates that are calibrated to the size-sensitivity
issue and find clade-with-drift modules that show the trace narrative
visually. The RF strict per-clade catalog (audit_48d) and audit_47-50
catalogs answer related but different questions; audit_47-50 is what
the manuscript will cite.

## Numbers documented for the record

## Headline numbers

Within-baseline-null gate (paired Wilcoxon real > null; 10 patients;
within-probe BH-FDR at m=6):

| band | hard RF n_above_null | hard p | hard q | **soft RF n_above_null** | **soft p** | **soft q** | Pat_03 drop p_null |
|:--|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| δ | 6/10 | 0.246 | 0.25 | 6/10 | 0.246 | 0.25 | 0.41 |
| θ | 7/10 (anti) | 0.032 | 0.10 | 8/10 (anti) | 0.019 | 0.037 | 0.037 |
| α | **9/10** | **0.005** | **0.029** | **9/10** | **0.0020** | **0.012** | **0.0039** |
| **β** | 7/10 | 0.053 ✗ | 0.105 ✗ | **9/10** | **0.0068** | **0.021** | **0.014** |
| low_γ | 7/10 | 0.080 | 0.12 | 7/10 | 0.042 | 0.063 | 0.064 |
| high_γ | 7/10 | 0.138 | 0.17 | 7/10 | 0.116 | 0.14 | 0.21 |

Controlled trace under soft RF: **α + β**. Soft RF brings β from BH-fail
to BH-pass; α is controlled under both variants.

## Why the soft variant works

The hard threshold counts a clade as "persisted" only if its best
Jaccard match crosses θ ≥ 0.70. β's task-induced reorganization
shifts ~20–50% of clade leaves between phases — clades reshuffle
without dissolving (best-J stays in [0.5, 0.85]) and without
crystallizing (best-J doesn't cross 0.85). Hard RF assigns these
graded reorganizations a discrete 0; soft RF assigns them their
actual best-J value and sums.

Per-pair Spearman ρ across 10 patients between T_KC λ=0 and
T_RF_hard at β: −0.47 (sign convention flipped; perfect agreement
would be ≈ −1). The two probes share some per-patient signal but
disagree at the band-cohort level because of the threshold.
Soft RF and KC λ=0 should track each other much more closely
(KC's m vector and soft RF's mean best-J both cumulate graded
shifts).

## Per-clade dendrogram visualization (audit_48d) — 4-way taxonomy

The strict gate is bestJ ≥ 0.65 (presence) and bestJ ≤ 0.30 (absence),
where bestJ = max Jaccard against ALL internal-node leafsets of the
target tree (NOT same-cut-level only — using same-cut-level mis-classifies
PERSIST as RESET when a small clade stays grouped INSIDE a larger
task-test cluster). Top-3 disjoint clades per (patient, band, mode)
are rendered on all 4 phase dendrograms (rs_pre, task_learn, task_test,
rs_post) with distinct colors.

Four mutually exclusive modes:

| mode | anchor | predicate | visual narrative |
|:--|:--|:--|:--|
| **trace** | task_test | bestJ(C, post) ≥ 0.65 AND bestJ(C, pre) ≤ 0.30 | scatter → cluster → persist |
| **reset** | rs_pre | bestJ(C, post) ≥ 0.65 AND bestJ(C, task) ≤ 0.30 | cluster → scatter → cluster |
| **persist** (= anchor) | rs_pre | bestJ(C, task) ≥ 0.65 AND bestJ(C, post) ≥ 0.65 | cluster across all phases |
| **rearrange** | rs_post | bestJ(C, pre) ≤ 0.30 AND bestJ(C, task) ≤ 0.30 | scatter → scatter → cluster |

Cohort coverage (cells with ≥ 1 matching clade / 10):

| band | trace | reset | persist | rearrange |
|:--|:--:|:--:|:--:|:--:|
| δ | 3 | 1 | 5 | **10** |
| θ | 3 | 3 | 4 | **9** |
| α | 3 | 3 | 5 | **9** |
| **β** | 3 | 4 | **9** | 7 |
| low_γ | 2 | 2 | 7 | 8 |
| high_γ | 2 | 0 | 2 | 5 |
| total / 60 | **16** | **13** | **32** | **48** |

Three striking observations:

1. **β is dominantly PERSIST (9/10 cells), NOT TRACE (3/10).**
   At the strict per-clade level, β's modules STAY THE SAME across
   phases — they don't dissolve in task and re-emerge. This is a
   meaningful reframe of the soft-RF cohort result: the 9/10
   above-null finding is asymmetric drift WITHIN persistent modules,
   not task-induced new modules persisting into rest.

2. **REARRANGE is near-universal (48/60 cells)**: rs_post carries
   spontaneous emergent structure absent from both rs_pre and
   task_test. δ, θ, α all show 9-10/10 patients with at least one
   rearrange clade. rs_post is NOT a "frozen task state"; it has
   its own band-specific organization.

3. **Strict TRACE (emergent-in-task AND persists-in-post) is rare**:
   16/60 cells; no band exceeds 3/10. Strict trace is the ideal
   mechanistic narrative but it is the LEAST common phenomenon at
   the per-clade level. This is consistent with the canonical
   2026-04-25 trace-modules J=0.9 cohort-null and the MRL cohort-null,
   reinterpreted: strict emergent-then-persist clades exist sparsely
   across the cohort; the cohort-level signal is in graded asymmetry
   (soft RF, KC), not in strict per-clade emergence.

Pat_06 β trace example: one 35-leaf clade scattered across rs_pre,
contiguous subtree in task_test, mostly contiguous in rs_post
(J_post=0.62, J_pre=0.36).
Pat_02 β persist example: two clades (9 + 8 leaves) clustered in ALL
4 phases (J_tt=0.78/0.70, J_post=0.82/0.89).
Pat_06 θ reset example: 2 clades (6 + 6 leaves) clustered in rs_pre
+ rs_post but scattered through task_learn + task_test.
Pat_07 θ rearrange example: 3 clades (9, 8, 7 leaves) emergent only
in rs_post (J_pre ≈ J_tt ≈ 0.10–0.16).

## Convergence with substrate + other LRG probes

| probe | α | β | low_γ |
|:--|:--|:--|:--|
| substrate d_S (raw FC, audit_25) | 8/10 trace, drift caveat | 7/10 trace, controls pass | 7/10 trace, controls pass |
| KC λ=0 (per-pair MRCA, audit_36) | passes BH | **passes BH (10/10 vs null)** | passes BH |
| soft RF (per-clade graded, audit_48b) | **9/10 above null, q=0.012** | **9/10 above null, q=0.021** | 7/10 above null, q=0.063 |
| per-clade dendrogram (visual, audit_48d) | 5/10 cells | **8/10 cells** | 6/10 cells |

α + β are convergent across substrate + per-pair + per-clade
+ per-clade-visual. β has the cleanest convergence at the per-clade
visual level (8/10 cells). Low_γ is borderline at every probe.

## What's no longer load-bearing

The 2026-04-25 MRL cohort-null and the 2026-04-25 Trace-Modules
strict-J=0.9 cohort-null are reinterpreted. They did NOT measure
"no clade-level reorganization exists." They measured "no clade-level
reorganization crosses our hard threshold (J=0.9 strict identity)."
The continuous soft variant + the per-clade dendrogram catalog at
J_post ≥ 0.65 / J_pre ≤ 0.30 show that clade-level reorganization
exists at α + β + low_γ — it's just GRADED rather than discrete.

The MRL ↔ CBR reconciliation report
(`.agents/guides/task-persistence-investigation/2026-04-25_mrl-vs-cbr-reconciliation.md`)
already anticipated this in 2026-04-25; the soft RF result
empirically confirms its diagnosis.

## Status: parked

RF is NOT scheduled for further development as a Section-5 measure.
The audit scripts and outputs remain on disk for reproducibility:

- audit_48 (hard RF threshold)
- audit_48b (soft RF mean best-Jaccard)
- audit_48c (real-vs-null scatter)
- audit_48d (per-clade strict 4-mode taxonomy on 4-phase dendrograms)

If the strict per-clade catalog (16/60 trace, 13/60 reset, 32/60
persist, 48/60 rearrange) becomes useful later as a strict-identity
sensitivity check against the audit_47-50 lenient catalogs, the data
is here. Otherwise this entire probe family is set aside.

Section 5.4's load-bearing measures are KC λ=0 (audit_36 cohort scalar +
audit_46 controls) and the audit_47-50 KC module-view 4-mode catalogs
(`2026-05-07_kc-{trace-network,reset,rearrangement,anchor}-modules.md`).
RF does not appear in the manuscript.

## Pointers

Tables: `data/reports/section_5_lrg_trace/14_rf_clade_persistence/tables/`
Figures: `data/reports/section_5_lrg_trace/14_rf_clade_persistence/figures/`
Per-clade catalog: `figures/trace_clades_on_dendrograms/{Pat_XX}_{band}.pdf` (29 files)
Real-vs-null scatter: `figures/rf_real_vs_null_hard_vs_soft.pdf`
Cohort heatmaps: `figures/rf_k_cohort_heatmap.pdf` (hard) +
                  `figures/rf_k_soft_cohort_heatmap.pdf` (soft)
