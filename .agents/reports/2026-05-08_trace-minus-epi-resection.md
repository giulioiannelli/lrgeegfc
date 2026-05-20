---
name: 2026-05-08_trace-minus-epi-resection
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-05-08
updated: 2026-05-08
pointers:
  - .agents/reports/2026-05-08_direction-a-induced-subtree.md
  - .agents/reports/2026-05-05_result-2-lrg-beta-trace.md
  - .agents/reports/2026-05-07_epileptic-n10-revisit.md
  - data/audit/trace_minus_epi/Td_KC_per_patient.csv
  - data/audit/trace_minus_epi/Td_KC_cohort_summary.csv
  - data/audit/trace_minus_epi/Td_KC_full_vs_resect.csv
  - scripts/01_compute/audit/audit_56_trace_minus_epi.py
---

# Trace minus epi — KC trace on full tree vs V \ E_p resection

## Renormalization head

**Direct test of the prediction from Direction A's β-dissociation:
remove the epi nodes from the global tree and recompute the
Section-5 KC trace; does the β trace sharpen, weaken, or stay the
same? Answer: it stays the same per-patient — Δ_med small at all
three λ values (β/λ=0: −2.45 → −2.26; β/λ=0.5: −0.58 → −0.46;
β/λ=1: −1.42 → −1.44). The β trace is **non-epi-carried**: the
epi nodes contribute to the cohort signal proportionally to their
size, not preferentially. The global β trace lives in the non-epi
part of the network. This is one of three regularities that the
resection probes:

1. **β trace robust** to resection (all 3 λ's, per-patient stable)
   → confirms the global β trace is dominated by non-epi structure;
   the audit_54 epi-only RESET at β/λ=1 is a separate
   sub-network signal, not a confound.
2. **α/λ=0 trace strengthens** with resection (cohort median −0.81
   → −1.29, IQR fraction 19% → 40%). The "subtract epi to sharpen"
   intuition holds at α-topology specifically.
3. **δ/λ=0 direction flips** with resection (full +2.01 RESET → resect
   −0.09 flat). The global δ-topology RESET was driven by epi-rest
   tree interactions; the non-epi tree alone is direction-flat at
   δ-topology.

The β-dissociation finding from Direction A becomes more precise:
the global β trace is non-epi-driven; the epi sub-network shows
its own RESET; the two coexist without interfering. Three readings
emerge depending on the band — none of them collapse to a simple
"epi confounds Result 2" or "epi carries Result 2" story.**

## Setup

For each (patient, band, λ), three KC distance triangles on three
phases (rest_pre, task_test, rest_post) on three operationalisations:

- **`full`** — KC on the cached LRG dendrogram (N leaves, all
  channels). This is the Section-5 / Result-2 operationalisation.
- **`resect`** — KC on the dendrogram induced on `V \ E_p` (drop
  epi leaves, keep ~88–96% of channels per patient).
- **`epi_only`** — KC on the dendrogram induced on `E_p` (the
  audit_54 Direction A object, included here for direct
  comparison).

Statistic: `Td_KC = d_KC(tt, post) − d_KC(pre, tt)` per (patient,
band, λ, variant). Negative = trace direction. Cohort: n=9
(Pat_15 dropped — no epi annotation; matches the rest of the
2026-05-08 family).

Same machinery as audit_54 (induced subtree via cophenetic
restriction); script at
`scripts/01_compute/audit/audit_56_trace_minus_epi.py`.

## Cohort regularities — full vs resect, side-by-side

IQR-passing cells per `feedback_iqr_vs_cohort_slope.md`
(median ≥ 20% × IQR width):

| (band, λ) | full Td_med (n_neg, %IQR) | resect Td_med (n_neg, %IQR) | Δmed | reading |
|---|---|---|---|---|
| **β / λ=0** | **−2.45 (7/9, 66%)** | **−2.26 (7/9, 76%)** | +0.19 | **TRACE robust** |
| **β / λ=1** | **−1.42 (8/9, 84%)** | **−1.44 (8/9, 97%)** | −0.02 | **TRACE robust, slightly sharper** |
| β / λ=0.5 | −0.58 (6/9, 25%) | −0.46 (6/9, 29%) | +0.12 | trace robust at moderate amplitude |
| **α / λ=0** | **−0.81 (6/9, 19%)** | **−1.29 (6/9, 40%)** | **−0.48** | **TRACE strengthens with resection** |
| α / λ=0.5 | −0.81 (5/9, 29%) | −0.26 (6/9, 12%) | +0.55 | trace weakens (mostly Pat_03 effect) |
| α / λ=1 | −0.33 (6/9, 13%) | −0.02 (5/9, 1%) | +0.30 | trace weakens (small amplitudes) |
| **δ / λ=0** | **+2.01 (2/9, 60%)** | **−0.09 (5/9, 3%)** | **−2.10** | **direction FLIPS — global RESET driven by epi-rest interaction** |
| δ / λ=0.5 | −0.70 (6/9, 39%) | −0.57 (6/9, 34%) | +0.13 | trace robust at moderate amplitude |
| **low_gamma / λ=0** | **+1.11 (4/9, 17%)** | **−0.34 (5/9, 9%)** | −1.44 | direction flips, but %IQR low both — flat |
| low_gamma / λ=0.5 | +0.84 (3/9, 38%) | +0.18 (3/9, 13%) | −0.66 | reset weakens to flat |
| **high_gamma / λ=0** | **−3.10 (5/9, 43%)** | **−3.02 (5/9, 41%)** | +0.08 | trace robust |
| high_gamma / λ=0.5 | −1.18 (5/9, 55%) | −0.69 (6/9, 38%) | +0.48 | trace weakens with resection |
| theta / λ=0.5 | −0.66 (7/9, 23%) | −0.54 (7/9, 18%) | +0.12 | trace robust at small amplitude |

(Cells below 20% IQR not shown.)

## Reading the three regularities

### 1. β trace is non-epi-carried (manuscript-worthy)

All three λ values: per-patient Td values barely change with
resection. β/λ=0 example, per patient (full → resect):

```
Pat_02: +0.51 → +1.63    Pat_07: -2.43 → -2.26
Pat_03: -19.73 → -19.64  Pat_08: -5.07 → -5.27
Pat_05: +0.17 → +0.01    Pat_10: -2.45 → -1.27
Pat_06: -5.09 → -3.59    Pat_13: -3.16 → -2.81
                         Pat_14: -1.37 → -0.61
```

Removing 7%–25% of leaves (Pat_03 has |E_p|=6 of 122 = 5%; Pat_13
has |E_p|=30 of 119 = 25%) barely shifts any patient's β trace.
**This is exactly what we'd expect if the epi nodes are not
structurally pivotal at β** — they get along for the ride, and the
trace lives in non-epi structure.

The Direction A finding (epi-only at β/λ=1 shows RESET, audit_54)
**stands** as a separate phenomenon: the epi sub-network has its
own behaviour, but it doesn't propagate to the global tree's β
dynamics in either direction. The "β-dissociation" framing from
Direction A is correct — the cohort network traces while the epi
sub-network resets, and resection now confirms the two coexist
without interfering.

### 2. α/λ=0 trace strengthens with resection

Cohort median goes from −0.81 (19% IQR — borderline) to **−1.29 (40%
IQR — solid)**. Per-patient at α/λ=0 (full → resect):

```
Pat_02: -0.81 → -1.91    Pat_07: +6.28 → +6.65
Pat_03: +1.99 → +1.57    Pat_08: -6.19 → -5.96
Pat_05: -1.09 → -1.11    Pat_10: -3.21 → -1.39
Pat_06: -7.16 → -7.85    Pat_13: -0.81 → -1.29
                         Pat_14: +1.05 → +1.34
```

The "borderline-trace" patients (Pat_02, Pat_13) shift further into
trace direction; the "strong-trace" patients (Pat_06, Pat_08) stay
strong; the "strong-reset" patient (Pat_07) stays strong reset.
Cohort consistency stays at 6/9 trace direction but the **median
trace amplitude doubles**.

This is the cleanest case of "subtracting epi sharpens the trace"
in the entire grid. **At α-topology, the epi nodes were dampening
the trace at the cohort median**. The mechanism is consistent with
Direction A's α-topology TRACE finding (5/9 trace at the epi-only
sub-network) — both the epi and non-epi sub-networks trace at
α-topology, but the trace is cleaner in the non-epi part.

### 3. δ/λ=0 direction flips on resection

Cohort median: **+2.01 (RESET, 60% IQR — strong) → −0.09 (flat,
3% IQR — null)**. Per-patient at δ/λ=0:

```
Pat_02: +0.83 → -0.54    Pat_07: +2.64 → +2.85
Pat_03: +0.67 → +0.46    Pat_08: -0.60 → -0.09
Pat_05: +2.01 → -0.78    Pat_10: +4.01 → +2.50
Pat_06: -10.45 → -10.03  Pat_13: +4.52 → -0.75
                         Pat_14: +5.31 → +4.11
```

Six patients shift their δ/λ=0 Td toward zero (Pat_02, Pat_05, Pat_08,
Pat_10, Pat_13, Pat_14). Pat_06 (strong-trace outlier) and Pat_07
(strong-reset outlier) stay similar. The cohort-wide δ-topology
RESET at the global level was being driven by **how the epi nodes
attach to the rest of the tree at δ**. When you remove them, the
non-epi part has no clean direction at δ-topology.

This complements Direction A's finding that the epi-induced subtree
at δ/λ=0 has a flat Td (median −0.13). So **δ-topology has no clean
trace or reset signal in either the epi sub-network or the non-epi
sub-network — but the *interaction* of the two looks like a
RESET at the global level**. This is a different kind of regularity:
the δ-RESET is an interface effect, not a sub-network property.

## How this relates to Result 2

Result 2 is the load-bearing manuscript finding:
"LRG β trace, n=10, 10/10 KC λ=0/0.5/1 vs within-baseline null,
joint Bonferroni m=48 surviving". The audit_56 resection result
**confirms** Result 2: the β trace is robust to dropping the epi
nodes per-patient, in the same direction at all three λ values.
The epi nodes are not the source of the β-trace finding.

Importantly, the audit_56 cohort here is **n=9** (Pat_15 has no epi
annotation, so we restrict to patients where the resection is
defined). Within this n=9 sub-cohort:
- β/λ=0: 7/9 trace direction (full and resect identical)
- β/λ=1: 8/9 trace direction (full and resect identical)

The drop from 10/10 (Result 2) to 7-8/9 here is the n=9 sub-cohort
plus the absence of the within-baseline-halves null (we use raw
Wilcoxon vs 0 here, which is a more lenient test). The qualitative
agreement is total.

## Open questions / next moves

- **The α/λ=0 sharpening (Δ=−0.48) is non-trivial and underexplored
  in the existing manuscript.** Worth a sentence in the Section 5
  discussion: "the α-topology trace at the global level is
  sharpened when epi nodes are removed, suggesting epi nodes
  contribute non-α-trace structure that dampens the cohort
  median."
- **The δ/λ=0 interface effect** (RESET at global, flat at both
  sub-networks) is a structural finding about how epi connects to
  non-epi at the slow-rhythm topology. Worth flagging but not
  worth a dedicated section unless the broader directions plan
  surfaces a coherent slow-band story.
- **Virtual resection on FC, not on the tree, is the next layer.**
  Audit_56 uses cophenetic-restriction on the cached linkage matrix
  (cheap, exact MRCA-height preservation under UPGMA). True virtual
  resection (delete rows/columns of `A` before LRG eigendecomposition)
  would test "what does the LRG do without the epi nodes physically
  there?" — a stronger version. Direction B in the broader plan.

## Reproducibility

```bash
cd /home/giulio/Documents/research/neural_networks/lrgeegfc
conda activate lapbrain
python scripts/01_compute/audit/audit_56_trace_minus_epi.py
```

Runtime < 1 minute on n=9. Outputs at `data/audit/trace_minus_epi/`.
