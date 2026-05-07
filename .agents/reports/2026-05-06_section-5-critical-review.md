---
name: 2026-05-06_section-5-critical-review
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-05-06
updated: 2026-05-06
pointers:
  - .agents/reports/2026-05-06_section-5-writing-brief.md
  - .agents/reports/2026-05-05_result-2-lrg-beta-trace.md
  - .agents/reports/2026-05-05_section-5-manuscript-draft.md
  - data/audit/section5_v2_round2/README.md
  - data/audit/section5_v2_round3_redo/README.md
---

# Section 5 — Critical review of the writing-agent brief

**Head.** Of the 8 items the brief had to defend, **2 stand as written, 5
must be reframed, 1 must be retracted entirely**. The headline claims that
collapse on contact with the data are: (i) the **"67× dendrogram amplification"**
is a unit-conversion artifact between two distance functions on different
scales — both probes have *identical* scale-invariant effect strength
(rank-biserial r_rb = −1, n_real_below_null = 10/10); (ii) the brief's
**"5/6 bands have no uncorrected global-probe trace"** silently excludes
CTM (measure 05) where α/β/low_γ are all BH-surviving controlled traces
at p ≤ 0.014 — α has been load-bearing for §5.3 of the manuscript; (iii)
the **β-band anatomy "Hippocampus 3.5×, p = 0.011" claim** does not
survive Bonferroni across the 48 cohort regions tested (threshold 0.001)
and **collapses entirely without Pat_03** (Hip 2.4× p = 0.124, fusiform
1.1× p = 0.534) — Pat_03, the 1024 Hz outlier, contributes 2/5 Hip and
3/5 β fusiform leaves, making the β anatomy headline a Pat_03-driven
finding that contradicts the brief's own "Pat_03 not load-bearing"
assertion. The **p = 0.000977 numbers** are the test-floor of an n = 10
one-sided paired Wilcoxon (1/2¹⁰), not measurements of effect strength:
5 of 42 cells max out the same floor and are not differentiable. The
**"partition-structure not edge-weight"** dichotomy is wrong: β trace
appears at every layer tested (edge-rank, edge-Pearson, edge-Frobenius,
LRG d_F, KC λ=0/0.5/1) — the LRG layer adds tree views, it does not
amplify a separate phenomenon. The **anti-aligned patients Pat_15 / Pat_07
each carry 9 β trace-leaves** under the per-leaf calibrated null —
the brief's implicit assumption that they are excluded from the anatomy
catalog is wrong.

**One claim survives all gates intact: low_γ ctx-lh-fusiform** (4.0×,
p = 0.00001 full cohort; 4.03×, p = 0.00002 without Pat_03; Pat_02
carries 9/13 leaves; survives Bonferroni m=48). This is the load-bearing
anatomy headline. The β anatomy headline must shift to either
ctx-lh-superiortemporal (uncorrected, Pat_03-robust) or be dropped.

The locked framing decisions §1a–1e of the brief need rewriting before any
manuscript pass. This document specifies the find/replace edits.

---

## Summary table

| # | item | verdict | required edit |
|:-:|:--|:--|:--|
| 1 | "67× dendrogram amplification" | **RETRACT** | Replace with: rank-biserial r_rb = −1 on both KC β λ=0 and d_F β; identical scale-invariant strength. Drop the volcano callout. |
| 2 | "p = 0.000977" / "joint Bonferroni m = 48, 6 cells surviving" | **REFRAME** | Clarify that p = 0.000977 = 1/2¹⁰ is the n=10 paired Wilcoxon floor; 5 cells max out the same floor. The headline is "all 5 cells reached the maximum-strength outcome of the n=10 one-sided test against the within-baseline null" — not "p < 0.001 effect size differentiates them." |
| 3 | "5/6 bands have no global-probe trace, only β + low_γ" | **RETRACT** | CTM has α/β/low_γ at BH-q = 0.027 within m=6 (manuscript §5.3 commitment). Measure 11 alone has α with 3 uncorrected p < 0.05 cells (Grassmann k=13 p=0.003 BH-surviving, d_P p=0.007, d_S p=0.019). Restore α to the controlled-trace headline. |
| 4 | "β enriched at Hippocampus 3.5× / fusiform 2.5×" | **REFRAME** | The hypergeometric controls per-region only, not over-region. Across 48 regions tested, Bonferroni p ≤ 0.001 — no β region survives. **Without Pat_03**: Hip → p = 0.124, fusiform → p = 0.534, both COLLAPSE. ctx-lh-superiortemporal at p = 0.030 (without Pat_03) survives uncorrected. **low_γ ctx-lh-fusiform survives both gates** (full p = 0.00001, no-Pat_03 p = 0.00002; Pat_02 carries 9/13 leaves) — the only anatomy claim that does. β anatomy is fragile; low_γ fusiform is robust. |
| 5 | "β = partition-structure shift, not edge-weight shift" | **REFRAME** | Substrate β shows on edge-rank d_S 7/10, edge-Pearson d_P 7/10, edge-Frobenius d_F 8/10. LRG β shows on KC topology, KC heights, d_F LRG matrix distance. β is multifacet-supported across every geometric layer, not unique to the LRG-partition view. Drop the dichotomy. |
| 6 | Anti-aligned patients (Pat_15, Pat_07) handled? | **REFRAME** | Pat_15 carries **9 β trace-leaves** in the per-leaf null catalog; Pat_07 carries **9 β trace-leaves**. The brief's framing that they are excluded is wrong. The trace-leaf catalog is built per-(patient, band) against each patient's own null — pro-cohort and anti-cohort patients both produce trace-leaves. Anatomy claim must acknowledge this contributor split or restrict to pro-cohort patients. |
| 7 | Non-independence in §5.6 limitations vs §5 head | **REFRAME** | Move "this is not independent replication; it is a methodologically distinct geometric layer" to the *first* sentence of §5 results (or the §5 head paragraph). Limitations is the wrong place. |
| 8 | Provenance: scripts cited per claim | **EXTEND** | Brief cites CSVs but not the script that produced each. Add `scripts/01_compute/audit/audit_NN_*.py` column to the defensible-claims table. |

---

## Item 1 — "67× dendrogram amplification" (RETRACT)

### Issue

§1a of the brief leads with "the LRG dendrogram amplifies the same
underlying signal ~67× over matrix distance at β (effect 5.3 vs 0.08 for
the same Wilcoxon p = 0.001)." The two effect sizes are computed as
`median(T_d_null) − median(T_d_real)` on:

- **KC β λ = 0**: a tree-distance metric. KC distance sums squared
  differences in pair-depth across N(N−1)/2 leaf pairs; for N ≈ 120 leaves
  with integer-valued common-ancestor depths up to N, the metric scales
  with N². Effect = 2.85 − (−2.44) = **5.29 in tree-distance units.**
- **d_F β**: a self-normalized Frobenius distance on triu(D̂(τ)),
  bounded in approximately [0, 1]. Effect = 0.053 − (−0.026) = **0.079
  (dimensionless, normalized).**

Dividing one by the other (5.29 / 0.079 ≈ 67) is a unit-conversion
artifact. The two metrics live on incomparable scales by construction.

### Verification

Required: scale-invariant comparison on the same statistical test, same
patients, both probes. The brief's own test was a one-sided paired
Wilcoxon `T_d_real < T_d_null` with n = 10 patients per cell. The natural
scale-invariant strength is **rank-biserial r_rb** (the normalized
Wilcoxon W statistic, bounded [−1, +1]).

```
KC β λ=0   (real vs within-null):  r_rb = −1.000,  n_real<null = 10/10
d_F β       (real vs within-null):  r_rb = −1.000,  n_real<null = 10/10
```

Both probes have **rank-biserial = −1** (every patient is in trace
direction against their own null) and Wilcoxon p = 0.000977 (the n = 10
one-sided floor). Scale-invariantly, **the two probes have identical
cohort effect strength.** The "67×" is a unit-conversion artifact that
disappears under any properly normalized comparison.

Source CSV: `data/audit/lrg_global_probe_controls/cohort_controls_summary.csv`
(rows: `probe=drank, band=beta, variant=d_F` and `probe=kc, band=beta,
variant=lambda=0.0`). Verification command in
`scripts/01_compute/audit/audit_41_lrg_global_probe_controls.py:240–296`.

### Required edit

**Brief §1a head paragraph: replace**

> "The most novel methodological contribution of Section 5 is **the LRG
> dendrogram amplifies the same underlying signal ~67× over matrix
> distance at β** (effect 5.3 vs 0.08 for the same Wilcoxon p=0.001)."

**with**

> "The β trace at the LRG layer is detected at maximum strength of the
> n=10 paired Wilcoxon (real vs within-baseline null) by both the
> KC tree-distance probe (λ=0/0.5/1) and the d_F LRG-matrix-distance
> probe — every patient is in trace direction against their own null
> in both views. The LRG layer does not amplify the substrate effect; it
> renders the same underlying β reorganization as a tree-topology +
> tree-heights + matrix-distance signature in addition to the edge-level
> signature visible in §4."

**Brief §3 "Do NOT write" table**: the row "Effect size 0.08 (D-rank) or
Effect size 5 (KC) without context" must be replaced with: "**Never cite
a numerical 'amplification ratio' between probes on different scales.**
The two probes' effect-size ratios are unit-conversion artifacts; cite
rank-biserial or n_real_below_null instead."

**Volcano figure (5A)**: the inline callout "the dendrogram (KC) amplifies
effect ~67× over matrix distance (D-rank) at β" must be removed. The
volcano remains useful to show *which* cells survive the Bonferroni line,
but the amplification narrative is wrong and should not appear in the
caption.

---

## Item 2 — "Joint Bonferroni m=48, 6 cells surviving at p=0.000977" (REFRAME)

### Issue

Round-2 reported KC β λ=0 absolute Wilcoxon p = 0.019 (n_trace = 7/10),
KC β λ=1 p = 0.042 (n=8/10), d_F β p = 0.042 — none of which survive
Bonferroni m=48 (threshold 0.00104). The brief claims four β cells
"survive joint Bonferroni m=48 at p=0.000977." The two p-values come
from different statistical tests and are not interchangeable.

### Verification

Two distinct statistical frameworks produced the two sets of p-values:

| framework | test | H0 | source | β KC λ=0 p |
|:--|:--|:--|:--|:--|
| **absolute Wilcoxon** (round-2) | one-sample paired Wilcoxon on `T_d` | `T_d ≥ 0` | `tables/Td_per_patient_per_band_lambda.csv` (cohort statistic) | **0.019** (n_trace = 7/10) |
| **within-baseline null** (measure 11) | paired Wilcoxon on `T_d_real − T_d_null` | `T_d_real ≥ T_d_null` | `data/audit/lrg_global_probe_controls/cohort_controls_summary.csv` | **0.000977** (n_real<null = 10/10) |

The within-baseline null version compares the real triangle scalar
against a per-patient null triangle constructed from rest halves
(RPreA → RPreB → RPostA, no task data). Under H0 "real ≈ null", the
paired Wilcoxon with n = 10 patients in the same direction gives the
floor p-value:

```
n=10 one-sided paired Wilcoxon, all 10 patients with real < null
→ W+ = 0,  W− = sum(1..10) = 55
→ p = 1 / 2^10 = 0.0009766
```

Five cells (KC β λ=0, KC β λ=0.5, KC β λ=1, d_F β, Grassmann low_γ k=13)
**all max out at this floor**. Reporting "p = 0.000977 across 5 cells"
means "5 cells exceed every permutation of the null with n=10 patients" —
not "5 cells have differentiated p-values." Beyond p = 0.000977 the test
cannot resolve further; you'd need n > 10 for a stricter floor.

### Required edit

The number m = 48 is structurally correct (42 controls + 6 per-leaf
bands) but the manuscript narrative around p = 0.000977 must distinguish:

**Manuscript §5 methods (existing draft, the "Methods paragraph")**: the
sentence "Wilcoxon paired test over n=10 patients tested whether the real
T_d was below the patient's own T_d^null" is correct but must be followed
by:

> "With n = 10 paired observations the one-sided paired Wilcoxon has a
> floor p-value of 1/2¹⁰ = 0.000977 (achieved when all 10 patients are
> in trace direction against their own null); five of the 42 cells
> tested reach this floor. Joint Bonferroni at m = 48 (42 controls
> cells + 6 per-leaf bands; threshold 0.05/48 = 0.00104) is satisfied
> for all 5 controls cells (β at KC λ ∈ {0, 0.5, 1}, β at d_F, low_γ at
> Grassmann k=13) and for low_γ per-leaf — 6 cells. The p-values do not
> differentiate strength among the 5 controls cells; the Wilcoxon test
> at n = 10 is at its resolution floor."

**Brief §1b "β is load-bearing, low_γ is secondary"**: keep the two-band
picture but add the framing:

> "Across the within-baseline null framework (measure 11), β has 4
> cells at the n=10 Wilcoxon floor (p = 0.000977): KC λ=0/0.5/1 and d_F.
> low_γ has 1 cell at the floor: Grassmann k=13. Within-probe BH-FDR
> q ≤ 0.05 (m=42) adds 7 more cells across α/β/low_γ. **β at the absolute
> Wilcoxon (T_d < 0) is at p = 0.019 (KC λ=0) / 0.042 (KC λ=1) / 0.042
> (d_F)** — does not survive Bonferroni m=42 in this stricter framing.
> The within-baseline null is the more sensitive test because it uses
> each patient's own drift floor as the comparator; the manuscript cites
> it as the controlled headline."

**Brief §2 defensible-claims table**: the row "β all 4 facets pass joint
Bonferroni m=48 — p=0.000977, 10/10 patients each" is correct *iff*
qualified as "(within-baseline null framework; n=10 paired Wilcoxon
floor)". Otherwise it is misleading by suggesting the four cells have
independently small p-values.

---

## Item 3 — "5/6 bands have NO uncorrected global-probe trace" (RETRACT)

### Issue

§1b of the brief asserts: "α / δ / θ / high_γ show no global-probe trace
under controls." §2 row says "5/6 bands have NO uncorrected global-probe
trace; only β + low_γ have any cell with p<0.05." This contradicts the
manuscript Section 5.3 commitment to CTM α/β/low_γ controlled trace.

### Verification

Two control frameworks, both valid, give different per-band verdicts:

**Framework A — CTM split-vs-drift (measure 05):** the LRG ultrametric
distance changes Δ_task and Δ_rest are correlated per-pair across N(N−1)/2
contact pairs; the Wilcoxon tests whether `ρ_split > ρ_drift` cohort-wide.
This is the load-bearing controlled headline cited in the manuscript
draft (§5.3) and in round-2 README §5.3.

| band | ρ_split med | n_split>0 | p split>drift | BH-q (m=6) | verdict |
|:-:|:-:|:-:|:-:|:-:|:--|
| **α** | **+0.115** | **8/10** | **0.0068** | **0.027** | **controlled trace** |
| **β** | **+0.222** | **8/10** | **0.0137** | **0.027** | **controlled trace** |
| **low_γ** | **+0.140** | **7/10** | **0.0098** | **0.027** | **controlled trace** |
| δ | +0.031 | 6/10 | 0.246 | 0.334 | marginal |
| θ | −0.049 | 3/10 | 0.278 | 0.334 | null |
| high_γ | −0.014 | 4/10 | 0.423 | 0.423 | null |

α at p = 0.007 is the strongest controlled cohort signal in any LRG
measure of this section. **The brief drops it without justification.**

**Framework B — within-baseline null (measure 11):** real vs null
triangle, paired one-sided Wilcoxon. α has uncorrected p < 0.05 cells:

| probe | band | variant | n_real<null | Wilcoxon p | BH-q |
|:--|:-:|:--|:-:|:-:|:-:|
| Grassmann | α | k=13 | 9/10 | **0.0029** | **0.010** ✓ |
| D-rank | α | d_P | 9/10 | **0.0068** | **0.021** ✓ |
| D-rank | α | d_S | 8/10 | **0.0186** | **0.035** ✓ |

α has **3 uncorrected p < 0.05 cells** and 3 BH-FDR surviving cells
within m = 42. The brief's "α has no global-probe trace" is false even
under its own (measure 11) framework.

Sources:
- `data/reports/section_5_lrg_trace/05_ctm_sigma_aggregate/tables/cohort_summary.csv`
- `data/audit/lrg_global_probe_controls/cohort_controls_summary.csv` (full 42-row table reproduced in measure 11 report)
- `scripts/01_compute/audit/audit_33_ctm_triangle.py` (CTM re-mining)
- `scripts/01_compute/audit/audit_41_lrg_global_probe_controls.py` (within-baseline null)

### Required edit

**Brief §1b "β is load-bearing, low_γ is secondary via different route"**:
must become a **three-band** picture, not two:

> "Three bands carry controlled cohort trace at the LRG global-probe
> layer:
> - **β** is load-bearing across multiple geometric facets: 4 cells at
>   the n=10 Wilcoxon floor (KC λ=0/0.5/1, d_F) under within-baseline
>   null (measure 11); BH-q = 0.027 in the CTM split-vs-drift framework
>   (measure 05); 7/10 trace at the absolute T_d < 0 Wilcoxon (KC λ=0
>   p = 0.019, KC λ=1 p = 0.042). Bonferroni-surviving and
>   Pat_03-dropout-robust at the topology view (KC λ=0, p=0.019 → 0.037
>   without Pat_03).
> - **α** is **CTM-load-bearing**: ρ_split = +0.115, p split-vs-drift =
>   0.0068, BH-q = 0.027 within m=6 (measure 05) — the strongest single
>   p-value in any LRG controlled framework. Within-baseline null
>   (measure 11) BH-surviving at Grassmann k=13 (p = 0.003) and at
>   D-rank d_P (p = 0.007) and d_S (p = 0.019). KC silent — α trace
>   does not crystallize into tree distance.
> - **low_γ** is multi-mode + heights: 1 cell at the n=10 Wilcoxon
>   floor (Grassmann k=13, measure 11); BH-surviving at KC λ=1
>   (p = 0.014, heights only) and d_F (p = 0.019); per-leaf calibrated
>   localization 10/10 patients ≥3 trace-leaves (measure 12, p = 0.001).
> δ / θ / high_γ are null at every controlled framework."

**Brief §2 defensible-claims table**: row "5/6 bands have NO uncorrected
global-probe trace" is wrong. Replace with:

> "Three bands (α, β, low_γ) show controlled cohort trace at the LRG
> layer; geometric mechanism differs by band (β = tree, α = multi-mode +
> rank, low_γ = multi-mode + heights). δ / θ / high_γ are null at every
> framework tested."

**Brief §6 Result-1 → Result-2 narrative arc**: the second paragraph
must be revised to include α, *not just β + anatomy*:

> "Section 5 (LRG): at the LRG dendrogram layer, three bands carry
> controlled cohort trace. β is the most multifacet-supported (every
> geometric view, n=10 Wilcoxon floor on 4 cells) and Pat_03-robust at
> the topology view. α has the strongest single CTM controlled p
> (0.007). low_γ shows multi-mode and heights signatures plus per-leaf
> calibrated localization. The bands separate by geometric mechanism:
> β = tree, α = multi-mode + rank, low_γ = multi-mode + heights." ...

The "left fusiform / Hippocampus" anatomy claim is reframed in Item 4
below.

---

## Item 4 — Anatomy enrichment per-leaf null + region multiple-comparison + Pat_03 sensitivity (REFRAME)

### Issue

§1c of the brief: "After hypergeometric correction for cohort-wide
contact density: β: Hippocampus (3.5×, p=0.011, 3/9 patients) + left
fusiform cortex (2.5×, p=0.045, 2/9)." The hypergeometric is a per-region
test against random implant coverage given the trace-leaves themselves.
It does **not** correct for (a) multiple regions tested simultaneously,
nor (b) the per-leaf null trace classification, nor (c) which patients
contribute the leaves.

### Verification

#### 4a — Multiple-comparison across regions

`scripts/01_compute/audit/audit_45_section5_headline_figures.py` filters
to regions with `enrichment > 1, n_contacts ≥ 5, n_pat_trace ≥ 2` for
the figure. The hypergeometric is computed per filtered region, not
across the full set of regions tested.

Re-running on the **full 48 cohort regions with ≥ 5 contacts** (β eligible
set: 48 regions; total cohort contacts N = 866; trace-leaves K = 46;
baseline rate K/N = 5.31 %):

| region | n_contacts | n_trace | n_pat | enrich | p_hyper |
|:--|:-:|:-:|:-:|:-:|:-:|
| Hip | 27 | 5 | 3 | 3.49× | 0.011 |
| ctx-lh-fusiform | 38 | 5 | 2 | 2.48× | 0.045 |
| ctx-lh-superiortemporal | 53 | 6 | 3 | 2.13× | 0.055 |
| ctx-lh-entorhinal | 10 | 2 | 1 | 3.77× | 0.095 |
| ctx-rh-medialorbitofrontal | 22 | 3 | 1 | 2.57× | 0.106 |
| ctx-lh-middletemporal | 85 | 5 | 3 | 1.11× | 0.478 |
| ... 42 more regions | | | | | |

Bonferroni m = 48 threshold = **p ≤ 0.001**. **Zero β regions survive
Bonferroni-across-regions.** 2 regions survive uncorrected (Hip, fusiform).

#### 4b — Pat_03 sensitivity for anatomy

The brief asserts Pat_03 (1024 Hz outlier) is not load-bearing for the β
finding (defensibly, for KC λ=0: p=0.019 → 0.037 without Pat_03). But
Pat_03 contributes:

```
Hippocampus β leaves by patient:
  Pat_02: 2,  Pat_03: 2,  Pat_13: 1,  others: 0
  → Pat_03 = 2/5 of cohort Hip leaves

ctx-lh-fusiform β leaves by patient:
  Pat_03: 3,  Pat_13: 2,  others: 0
  → Pat_03 = 3/5 of cohort fusiform leaves
```

Re-running the hypergeometric **without Pat_03** (β: K = 40,
N = 866, baseline 4.62 %):

| region | n_trace | enrichment | p_hyper |
|:--|:-:|:-:|:-:|
| Hip | 3 | 2.41× | **0.124** (collapses) |
| ctx-lh-fusiform | 2 | 1.14× | **0.534** (collapses entirely) |
| ctx-lh-superiortemporal | 6 | 2.45× | **0.030** (survives uncorrected) |

**The β anatomy headline (Hippocampus + fusiform) is Pat_03-driven.**
Without Pat_03 the only surviving claim is left superior temporal
(uncorrected, n = 3 patients).

This contradicts the brief's §3 row "Pat_03 drives the β finding —
NOT defensible: KC λ=0 survives Pat_03 dropout." That row is correct
*for the global probe* but **wrong for the anatomy claim.**

#### 4c — Per-leaf null

The 611 trace-leaves are classified by a per-(patient, band) τ_95(p, b)
null floor (measure 12, audit_42). The classification rule is:

```
ρ_ℓ_demeaned ≥ τ_95(p, b)   where τ_95 = 95th percentile of ρ_ℓ_null over leaves
```

The within-leaf null is computed from `Spearman((D^pre_B − D^pre_A)[ℓ,·],
(D^post_B − D^post_A)[ℓ,·])`. Source:
`scripts/01_compute/audit/audit_42_per_leaf_rho_null.py`,
`data/reports/section_5_lrg_trace/12_per_leaf_rho_null/tables/cohort_calibration_summary.csv`.

This null **is** in place and is the right calibration. The independent
question (4a + 4b above) is whether the *region-level* enrichment built
on top of the trace-leaves survives multiple-region correction and
Pat_03 sensitivity.

#### 4d — Patient-overlap of contributors

Hippocampus β contributors: Pat_02, Pat_03, Pat_13.
ctx-lh-fusiform β contributors: Pat_03, Pat_13.
ctx-lh-superiortemporal β contributors: 3 patients (full names not
re-extracted; the catalog shows `n_pat_trace = 3` from the figure
script's filter).

**Pat_13 contributes to both Hip and fusiform; Pat_03 contributes to
both.** The "Hip + fusiform are independent regions co-enriched at β"
reading is partially false: 1 of 3 Hip patients and 1 of 2 fusiform
patients overlap (Pat_13), and Pat_03 is in both.

### Required edit

**Brief §1c "Anatomy = Hippocampus + left fusiform"**: must become

> "Anatomy at β is fragile under (a) multiple-region Bonferroni and (b)
> Pat_03 dropout. Across 48 cohort regions tested, no β region survives
> Bonferroni m=48 (threshold p ≤ 0.001). Two regions survive
> uncorrected: Hippocampus (3.5×, p = 0.011, 3 patients) and
> ctx-lh-fusiform (2.5×, p = 0.045, 2 patients). **Both collapse without
> Pat_03**: Hip → 2.4×, p = 0.124; fusiform → 1.1×, p = 0.534. The only
> region whose enrichment **strengthens** without Pat_03 is
> ctx-lh-superiortemporal (3 patients, 6 leaves; full cohort p = 0.055
> borderline; without Pat_03 p = 0.030, uncorrected only). The honest
> β anatomy claim is: 'left superior temporal cortex carries
> sampling-corrected trace-leaves at uncorrected p = 0.030 (without the
> 1024 Hz outlier patient); Hippocampus and fusiform claims are
> Pat_03-driven and do not survive his removal.' low_γ enrichment at
> ctx-lh-fusiform survives Bonferroni m=48 — verify with
> Pat_03-dropout sensitivity before citing as load-bearing."

**Brief §3 "Do NOT write" table**: add new row

> | "Hippocampus drives the β trace localization" | The Hip enrichment
> at p=0.011 fails Bonferroni m=48 across 48 cohort regions (threshold
> 0.001) and collapses to p=0.124 without Pat_03 (the 1024 Hz outlier).
> Pat_03 contributes 2 of the 5 cohort Hip leaves. The defensible β
> anatomy claim, after Pat_03 sensitivity, is left superior temporal
> cortex only. |

**Brief §2 defensible-claims table**: rows "β trace enriched at
Hippocampus 3.5× p=0.011" and "β trace enriched at left fusiform 2.5×
p=0.045" must be marked as **NOT defensible** under multiple-region or
Pat_03 controls. Replace with one row:

> | β trace anatomy: ctx-lh-superiortemporal (uncorrected p=0.030 without Pat_03) | Sampling-corrected hypergeometric, 3 patients, 6 leaves; only β region whose enrichment strengthens under Pat_03 dropout. | `data/audit/lrg_localization_anatomy/region_enrichment.csv` (recompute without Pat_03 row) |

**Manuscript anatomy paragraph (§5.5)**: must be rewritten to lead with
the multiple-comparison + Pat_03 sensitivity caveats, not bury them.

**low_γ ctx-lh-fusiform Pat_03 sensitivity (verified 2026-05-06):**
Per-patient breakdown of the 13 cohort low_γ fusiform trace-leaves:
**Pat_02: 9 leaves**, Pat_03: 2 leaves, Pat_13: 2 leaves. Pat_02 carries
the bulk; Pat_03 contributes ≈ 15 %.

Hypergeometric without Pat_03 (K = 70, N = 770, fusiform cohort
contacts = 30):
```
low_γ fusiform without Pat_03: n_trace=11, enrichment=4.03×, p=0.00002
```

vs full cohort (3.95×, p = 0.00001). **The low_γ ctx-lh-fusiform
enrichment is robust to Pat_03 dropout** — it is the ONE anatomy claim
in either band that survives BOTH multiple-region Bonferroni m=48 (full
cohort p = 0.00001 < threshold 0.001) AND Pat_03 sensitivity. It can
stand as the load-bearing anatomy headline for low_γ. The β anatomy
claims (Hip, fusiform) do not have a comparable robustness profile.

**Figure 5D anatomy_2d.pdf**: caption must declare the Pat_03 sensitivity.
Either re-render the bar chart twice (full cohort and without Pat_03
side-by-side) or remove Hip + fusiform as headlines and lead with
ctx-lh-superiortemporal.

---

## Item 5 — "β = partition-structure, NOT edge-weight" (REFRAME)

### Issue

§1d, §6 narrative arc, and the manuscript draft assert: "the β trace at
LRG is a *partition-structure* shift, not just an edge-weight shift."
This dichotomy implies that LRG reveals something that edge-level analysis
misses. It does not. β shows on every layer.

### Verification

| layer | probe | β cohort | source |
|:--|:--|:-:|:--|
| **Substrate edge-rank** | d_S (1 − Spearman(triu(A))) | 7/10, p = 0.042 | `data/audit/raw_fc_phase_distance/final_verdict_table_with_controls.csv` (§4 Result 1) |
| **Substrate edge-Pearson** | d_P | 7/10 | same |
| **Substrate edge-Frobenius** | d_F | 8/10 | same |
| **LRG matrix-rank** | d_S on D̂(τ) | 8/10, p = 0.042 (round-2) | `02_d_rank_triangle/tables/Td_per_patient_per_band.csv` |
| **LRG matrix-Pearson** | d_P on D̂(τ) | 8/10, p = 0.042 (round-2) | same |
| **LRG matrix-Frobenius** | d_F on D̂(τ) | 8/10, p = 0.042 (round-2); 10/10 vs within-baseline null at floor (m11) | same + `lrg_global_probe_controls/cohort_controls_summary.csv` |
| **LRG tree topology** | KC λ = 0 | 7/10, p = 0.019 (absolute); 10/10 vs within-baseline null at floor | `03_kc_lambda_triangle/tables/Td_per_patient_per_band_lambda.csv` |
| **LRG tree heights** | KC λ = 1 | 8/10, p = 0.042 (absolute); 10/10 vs within-baseline null at floor | same |
| **LRG tree balanced** | KC λ = 0.5 | n/r at absolute; 10/10 vs within-baseline null at floor | same |
| **LRG eigenstructure** | Grassmann k = 13 | 5/10 absolute (NULL); 8/10 vs null (BH-surviving) | `04_grassmann_triangle/tables/Td_per_patient_per_band_k.csv` |
| **LRG ultrametric ρ_split** | CTM | 8/10, p split-vs-drift = 0.014, BH-q = 0.027 | `05_ctm_sigma_aggregate/tables/cohort_summary.csv` |

**β is multifacet-supported**, not partition-unique. Every probe except
Grassmann at low k is positive for β at the cohort level, and Grassmann
at higher k gets β to 7-8/10 (round-2 k-sweep).

### Required edit

**Brief §1d "Result 2 = structural enrichment of Result 1"**: keep this
framing (it is correct that LRG is the same data, different layer), but
**delete** the sentence "The novel content of Section 5 over Section 4:
partition-structure shifts (LRG-layer) vs edge-weight shifts (raw-FC).
Different question, related answer." Replace with:

> "The novel content of Section 5 over Section 4 is multifacet support:
> β trace at LRG passes through tree topology (KC λ=0), tree heights
> (KC λ=1), balanced topology+heights (KC λ=0.5), LRG matrix distance
> (d_F, d_P, d_S), and CTM ultrametric per-pair correlation.
> Section 5 does not amplify the substrate effect; it shows that the
> substrate β edge-shift produces a coherent restructuring detectable at
> every coarser-grained geometric layer. The LRG-as-test is therefore an
> internal-consistency check that the β reorganization is not localized
> to a specific noise mode of one distance metric. The
> *partition-structure-only* framing is wrong — β reorganizes both
> partition (tree topology + heights) AND magnitude (d_F at LRG and
> substrate)."

**Brief §3 "Do NOT write" table**: add row

> | "β shows up at LRG that was invisible at edge level" or "β = LRG-layer-only finding" | β shows on edge-rank d_S 7/10, edge-Pearson d_P 7/10, edge-Frobenius d_F 8/10 (substrate, §4); the LRG layer adds tree-topology + tree-heights views of the same effect. The right framing is "multifacet-supported across every layer", not "partition-only". |

**Manuscript §5.2 results paragraph (existing draft)**: the sentence
"β trace lives in EVERY geometric feature" (in measure 11 report,
copied verbatim into the Result 2 writeup) is correct and should
anchor the §5.2 narrative. The Section 5 head paragraph (currently:
"To test whether this trace also imprints on the multiscale partition
structure of the network") needs to broaden:

> "To test whether the β edge-rank trace from §4 carries into the
> multiscale geometry of the network, we applied the LRG framework
> ... [as before] ... and confirmed that β reorganization is detectable
> at every coarser-grained layer of the resulting dendrogram: tree
> topology (KC λ = 0), tree heights (KC λ = 1), balanced (λ = 0.5),
> matrix-distance (d_F), and per-pair ultrametric correlation (CTM)."

---

## Item 6 — Anti-aligned patients in the trace-leaf catalog (REFRAME)

### Issue

The brief mentions cohort heterogeneity (6 pro / 2 mixed / 2 anti) but
does not explicitly state how the 2 anti patients (Pat_15, Pat_07) are
handled in the per-leaf trace-leaf catalog. The user's review prompt asks
for verification that the anatomy excludes them from the catalog (the
brief implicitly assumes this) and asks what their anti-trace anatomy
looks like.

### Verification

The per-leaf calibration uses each patient's **own** within-baseline null
(τ_95(p, b)). A leaf is flagged trace if `ρ_ℓ_demeaned ≥ τ_95(p, b)`, regardless
of whether the patient is pro or anti at the global-probe level.
Pat_15 and Pat_07 each have **per-patient τ_95(p, b) thresholds** and produce
trace-leaves above them.

Counts in the per-leaf catalog at all bands (`per_trace_leaf.csv`):

| patient | β | δ | θ | α | low_γ | high_γ | total |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| Pat_15 (anti-aligned) | **9** | 3 | 9 | 3 | 4 | 11 | 39 |
| Pat_07 (anti-aligned) | **9** | 9 | 8 | 14 | 10 | 12 | 62 |

**Both anti-aligned patients carry 9 β trace-leaves each.** They are not
excluded from the catalog. The catalog measures localized leaves above
each patient's own drift floor — it can flag leaves in patients whose
*global* (cohort-aggregated triangle T_d) signal is anti-trace.

This means the cohort anatomy bar charts at β (87 leaves total over 9
contributing patients; Pat_05 has 0 leaves) **include** 18 leaves from
the two anti-aligned patients (≈ 21 % of the cohort β catalog).

The brief's §3 row "Pat_03 drives the β finding — NOT defensible: KC λ=0
survives Pat_03 dropout (p=0.037)" is correct for the global probe, but
the analogous question for the anti-aligned patients **has not been
asked**: does the per-leaf catalog at β survive removing Pat_15 and
Pat_07? And what regions do their 18 β trace-leaves cluster in?

The anti-aligned patients' anatomy is not characterized.

### Required edit

**Brief §1c anatomy framing**: must add

> "The per-leaf trace catalog includes leaves from the two anti-aligned
> patients (Pat_15: 9 β leaves; Pat_07: 9 β leaves; combined ~21 % of
> the cohort β catalog). Their inclusion is by design — the per-leaf
> calibration uses each patient's own within-baseline null, and it can
> flag localized leaves in patients whose global triangle T_d is
> anti-trace. Two consequences:
> 1. The cohort anatomy bar charts at β count 18 leaves (Pat_15 + Pat_07)
>    that come from globally anti-aligned patients. Their regional
>    pattern has not been characterized — anti-aligned anatomy is an
>    open question.
> 2. Sampling-corrected enrichment at, e.g., ctx-lh-superiortemporal
>    (3 patients, 6 leaves cohort) does not isolate which patients are
>    pro vs anti-aligned; restrict-to-pro-cohort sensitivity test
>    needed before manuscript citation."

**Manuscript §5.5 anatomy paragraph (existing draft)**: must declare

> "The per-leaf calibrated catalog (measure 12) flags leaves above
> each patient's within-baseline drift floor; it does not require the
> patient to be pro-cohort at the global triangle level. The catalog
> therefore includes contributions from Pat_15 and Pat_07 (the two
> globally anti-aligned patients): 9 β trace-leaves each, 18 of the 87
> cohort β trace-leaves. We do not characterize their anatomy
> separately; restrict-to-pro-cohort sensitivity tests of the regional
> enrichment are deferred to revision."

**Optional follow-up compute**: re-run audit_45's regional enrichment
with the 8-patient pro-cohort subset (n=10 minus Pat_15, Pat_07) and
report side-by-side. This is a quick re-run if the user wants the brief
locked at full sensitivity.

---

## Item 7 — Non-independence framing position (REFRAME)

### Issue

The brief proposes putting the non-independence statement in the
limitations paragraph (§5.6). The user's read is that this should be the
*first* sentence of §5 results, since the manuscript already commits in
the §5 head to "we applied the LRG framework to the same imaginary-coherency
matrices ... compared dendrogram triangles."

### Verification

The manuscript draft Section 5 head paragraph
(`2026-05-05_section-5-manuscript-draft.md` lines 24–40) begins:

> "The raw functional-connectivity analysis (Section 4) revealed a
> β-band task-shaped trace in resting-state networks at the channel-pair
> level. To test whether this trace also imprints on the multiscale
> partition structure of the network, we applied the Laplacian
> Renormalization Group (LRG) framework to each (patient, phase, band)
> connectivity matrix..."

This implicitly says "LRG is applied to the same matrices". It does
**not** explicitly say "Section 5 is a methodologically distinct view of
the same data, not an independent replication." The current limitations
paragraph (§5.6, lines 121–144) does say it explicitly — but only at the
end.

### Required edit

**Manuscript §5 head paragraph**: insert after the first sentence (after
"...at the channel-pair level."):

> "Because the LRG analysis is applied to the same imaginary-coherency
> matrices analyzed in §4, Section 5 is a methodologically distinct
> geometric view of the same data — not an independent replication.
> Cross-cohort replication remains required for a generalizability claim."

The current §5.6 sentence on this is no longer needed in the limitations
paragraph (or can be condensed).

**Brief §1d**: change "The reviewer-defense move is in the *limitations
paragraph* of the manuscript draft" to "The non-independence is stated
up front in the §5 head paragraph; the limitations paragraph adds the
specific statistical caveats but no longer carries the framing
disclosure."

---

## Item 8 — Provenance: scripts cited per claim (EXTEND)

### Issue

§2 of the brief gives a defensible-claims table with CSV paths but no
script paths. A reviewer should be able to trace each number to (a) the
CSV cell, (b) the analysis script, and (c) the pipeline that produced
the inputs to that script.

### Verification

For each row of the brief's §2 table, the script path and the input it
consumes:

| brief row | CSV | producing script | consumes |
|:--|:--|:--|:--|
| β all 4 facets at p=0.000977 | `data/audit/lrg_global_probe_controls/cohort_controls_summary.csv` | `scripts/01_compute/audit/audit_41_lrg_global_probe_controls.py` | `data/cache/imcoh_lrg_halves/{Pat_NN}/{band}_{phase}_lrg_imcoh-abs.npz` (half-LRG cache) + measures 02/03/04 real T_d tables |
| β KC λ=0 Pat_03-dropout p=0.037 | `data/audit/lrg_global_probe_controls/pat03_dropout.csv` | same audit_41 (lines 304–331) | `data/reports/section_5_lrg_trace/03_kc_lambda_triangle/tables/Td_per_patient_per_band_lambda.csv` |
| 67× amplification | (computed inside figure) | `scripts/01_compute/audit/audit_45_section5_headline_figures.py` | same controls CSV — **but the claim is wrong, see Item 1** |
| 5/6 bands no uncorrected — only β + low_γ | (filtered subset of controls) | filter operation in audit_45 — **but the claim is wrong, see Item 3** | controls CSV |
| low_γ Grassmann k=13 p=0.000977 | controls CSV | audit_41 | half-LRG cache |
| 9/10 patients ≥3 calibrated trace-leaves at β | `data/reports/section_5_lrg_trace/12_per_leaf_rho_null/tables/cohort_band_summary.csv` | `scripts/01_compute/audit/audit_42_per_leaf_rho_null.py` | half-LRG cache |
| Bidirectional cross-validation at β (enrichment p=0.004, concentration p=0.004) | `data/reports/section_5_lrg_trace/09_global_local_cross_validation/tables/cohort_band_summary.csv` (or `data/audit/lrg_cross_validation/`) | `scripts/01_compute/audit/audit_43_global_local_cross_validation.py` | per-leaf and CTM σ outputs |
| β Hippocampus 3.5× p=0.011 | `data/audit/lrg_localization_anatomy/region_enrichment.csv` | `scripts/01_compute/audit/audit_45_section5_headline_figures.py` (function `make_anatomy_2d`) | `data/reports/section_5_lrg_trace/13_anatomical_localization/tables/per_trace_leaf.csv` (output of `audit_44_anatomical_mapping.py`) — **but multiple-region Bonferroni and Pat_03 sensitivity not applied, see Item 4** |
| β fusiform 2.5× p=0.045 | same | same | same — **same caveats** |
| low_γ fusiform 4.0× p<0.001 | same | same | same — **Pat_03-sensitivity unverified, see Item 4** |

### Required edit

**Brief §2 defensible-claims table**: add a fourth column "producing
script" populated as above. **Replace** the wrong rows (67×, "5/6 bands",
"Hippocampus", "fusiform") per Items 1, 3, 4.

**Manuscript reproducibility appendix (§7 of brief)**: the
verification snippet provided at brief §7 prints β Bonferroni m=42
surviving cells, Pat_03 dropout, and β enriched regions. It does not
print **without-Pat_03 anatomy** — must be added:

```python
# After the brief's existing snippet:
import pandas as pd
from scipy.stats import hypergeom

leaves = pd.read_csv(
    "data/reports/section_5_lrg_trace/13_anatomical_localization/"
    "tables/per_trace_leaf.csv"
)
leaves = leaves[~leaves.region.isin(["Wm", "Unk"])]
beta_no03 = leaves[(leaves.band == "beta") & (leaves.patient != "Pat_03")]
# (re-compute hypergeometric without Pat_03; expect Hip p~0.124, fusiform p~0.534)
```

---

## Lock-in summary (after all 8 items)

The §1 framing decisions of the brief are now revised as follows. **These
are the locked numbers and framings the writing pass must build on. No
further negotiation.**

### 1a (revised) — Lead with multifacet support, not amplification

β trace is detected at the n=10 Wilcoxon floor (within-baseline null) by
4 cells: KC λ=0 (topology), KC λ=0.5 (balanced), KC λ=1 (heights), d_F
(LRG matrix distance). Rank-biserial r_rb = −1 on each. The LRG layer
adds tree-distance and ultrametric matrix-distance views of the same β
reorganization detectable at the substrate edge level (§4); it does not
amplify the underlying signal. **No "67×" claim anywhere in the
manuscript.**

### 1b (revised) — Three-band picture, not two

- **β** load-bearing across multiple geometric facets (4 cells at
  Wilcoxon floor; CTM BH-q = 0.027; KC λ=0 Pat_03-dropout-robust at
  p = 0.037).
- **α** CTM-load-bearing (p = 0.0068, the strongest single controlled
  cohort p in any LRG framework); within-baseline-null BH-surviving at
  Grassmann k=13 + D-rank d_S/d_P.
- **low_γ** multi-mode + heights (Grassmann k=13 at Wilcoxon floor; KC
  λ=1 BH-surviving; per-leaf calibrated 10/10 patients at floor).

### 1c (revised) — Anatomy split: β fragile, low_γ robust

**β** has zero regions surviving Bonferroni m=48 (threshold p = 0.001).
Without Pat_03, the only surviving uncorrected β claim is **left
superior temporal cortex (3 patients, p = 0.030)**. Hippocampus and
fusiform β headlines collapse without Pat_03 (Hip p = 0.124; fusiform
p = 0.534). **No β "Hippocampus" or "fusiform" claims anywhere in the
manuscript without explicit Pat_03 caveat.**

**low_γ ctx-lh-fusiform** is the load-bearing anatomy headline: 4.0×
enrichment full cohort (p = 0.00001, n=3 patients with Pat_02 carrying
9/13 leaves), survives Bonferroni m=48 (p ≪ 0.001 threshold), and
**survives Pat_03 dropout** (4.03×, p = 0.00002 without Pat_03). This
is the only anatomy claim in either band that meets both robustness
gates. low_γ inferior temporal (uncorrected p = 0.026) is a secondary
borderline claim; without dropout test should not be cited as
manuscript-load-bearing.

### 1d (revised) — Section 5 is a multifacet view of §4, not an independent replication

Non-independence stated up front in §5 head paragraph. β multifacet
support across every layer (substrate edge-rank/Pearson/Frobenius +
LRG matrix-rank/Pearson/Frobenius + KC topology/heights/balanced + CTM
+ Grassmann at higher k) is the §5 novel content. **Drop the
"partition-not-edge" dichotomy.**

### 1e — Voice and tone (unchanged)

Third person plural, past tense; renormalization-style heads;
TRACE / ANCHOR / RESET / EMERGENT taxonomy; never "persistence".

### Joint MTC (revised numbers)

| scope | m | threshold | surviving cells |
|:--|:-:|:-:|:--|
| controls only (measure 11, 6 bands × 7 probes) | 42 | 0.00119 | **5** (β KC λ=0/0.5/1 + d_F; low_γ Grassmann k=13) — all at Wilcoxon floor |
| controls + per-leaf (m11 + m12) | 48 | 0.00104 | **6** (above + low_γ per-leaf) |
| add cross-validation (+m09) | 60 | 0.000833 | **5** (loses low_γ per-leaf) |
| **β at the absolute T_d < 0 Wilcoxon (round-2)** | 42 | 0.00119 | **0** (best: KC β λ=0 at p = 0.019) |
| **CTM split-vs-drift (measure 05) within m=6 BH-FDR** | 6 | q ≤ 0.05 | **3** (α p=0.0068, β p=0.014, low_γ p=0.010) |
| **Joint BH-FDR across 5 LRG probe families (round-2 absolute, m=84)** | 84 | q ≤ 0.05 | **0** (min q = 0.287) |

The honest joint-MTC narrative for the manuscript: the within-baseline
null framework at n=10 is a *more powerful* test than the absolute
Wilcoxon (paired against each patient's own drift), and it surfaces 6
joint-Bonferroni-surviving cells. **The 5 controls cells + 1 per-leaf
cell at p = 0.000977 are at the n=10 paired Wilcoxon floor; they do not
differentiate strength**. The absolute Wilcoxon (T_d < 0) does
differentiate but does not survive m=84 joint BH-FDR.

The manuscript should cite both frameworks: within-baseline null as the
primary controlled framework (m=42 controls + m=6 per-leaf, 6 cells
joint-Bonferroni surviving); CTM as the second controlled framework
(m=6 within-probe, 3 cells BH-surviving — α/β/low_γ); the absolute
Wilcoxon and joint BH-FDR across 5 probes as the conservative
honesty statement (no individual cell joint-significant).

### Anti-aligned patients

Pat_15 and Pat_07 carry 9 β trace-leaves each in the per-leaf catalog.
The 87 β cohort trace-leaves include 18 from anti-aligned patients
(≈ 21%). Their regional pattern is uncharacterized; restrict-to-pro
cohort sensitivity test needed before submission of any anatomy
claim.

---

## Files this review references

- Brief: `.agents/reports/2026-05-06_section-5-writing-brief.md`
- Manuscript draft: `.agents/reports/2026-05-05_section-5-manuscript-draft.md`
- Result 2 writeup: `.agents/reports/2026-05-05_result-2-lrg-beta-trace.md`
- Round-2 README: `data/audit/section5_v2_round2/README.md`
- Round-3 redo README: `data/audit/section5_v2_round3_redo/README.md`
- Measure 05 (CTM): `data/reports/section_5_lrg_trace/05_ctm_sigma_aggregate/report.md`
- Measure 11 (controls): `data/reports/section_5_lrg_trace/11_lrg_global_probe_controls/report.md`
- Measure 12 (per-leaf null): `data/reports/section_5_lrg_trace/12_per_leaf_rho_null/report.md`
- Measure 13 (anatomy): `data/reports/section_5_lrg_trace/13_anatomical_localization/report.md`
- Per-leaf catalog: `data/reports/section_5_lrg_trace/13_anatomical_localization/tables/per_trace_leaf.csv`
- Region enrichment: `data/audit/lrg_localization_anatomy/region_enrichment.csv`
- Controls CSV: `data/audit/lrg_global_probe_controls/cohort_controls_summary.csv`
- Audit 41 script: `scripts/01_compute/audit/audit_41_lrg_global_probe_controls.py`
- Audit 42 script: `scripts/01_compute/audit/audit_42_per_leaf_rho_null.py`
- Audit 43 script: `scripts/01_compute/audit/audit_43_global_local_cross_validation.py`
- Audit 44 script: `scripts/01_compute/audit/audit_44_anatomical_mapping.py`
- Audit 45 script: `scripts/01_compute/audit/audit_45_section5_headline_figures.py`

## Verification commands (re-runnable in `lapbrain`)

```bash
# Item 1: rank-biserial verification (KC β λ=0 vs d_F β under within-baseline null)
/home/giulio/Documents/miniconda3/envs/lapbrain/bin/python -c "
import pandas as pd
import numpy as np
kc_real = pd.read_csv('data/reports/section_5_lrg_trace/03_kc_lambda_triangle/tables/Td_per_patient_per_band_lambda.csv')
kc_null = pd.read_csv('data/audit/lrg_global_probe_controls/kc_null.csv')
sub_real = kc_real[(kc_real.band=='beta') & (kc_real.lam==0.0)].set_index('patient')['T_KC']
sub_null = kc_null[kc_null.band=='beta'].set_index('patient')['null_T_lam0.0']
real_kc, null_kc = sub_real.values, sub_null.loc[sub_real.index].values
print(f'KC β λ=0 n_real_below_null = {(real_kc<null_kc).sum()}/{len(real_kc)}')
"

# Item 4: anatomy without Pat_03
/home/giulio/Documents/miniconda3/envs/lapbrain/bin/python -c "
import pandas as pd
from scipy.stats import hypergeom
leaves = pd.read_csv('data/reports/section_5_lrg_trace/13_anatomical_localization/tables/per_trace_leaf.csv')
leaves = leaves[~leaves.region.isin(['Wm','Unk'])]
PATIENTS = ['Pat_02','Pat_03','Pat_05','Pat_06','Pat_07','Pat_08','Pat_10','Pat_13','Pat_14','Pat_15']
N = 866  # cohort cortical+subcortical contacts
beta_no03 = leaves[(leaves.band=='beta') & (leaves.patient!='Pat_03')]
K = len(beta_no03)
for reg in ['Hip','ctx-lh-fusiform','ctx-lh-superiortemporal']:
    n_trace = (beta_no03.region==reg).sum()
    n_contacts = {'Hip':27,'ctx-lh-fusiform':38,'ctx-lh-superiortemporal':53}[reg]
    p = float(hypergeom.sf(n_trace-1, N, K, n_contacts))
    print(f'{reg}: n_trace={n_trace}, p_no_pat03={p:.3f}')
"
```

Expected output for the second block:
```
Hip: n_trace=3, p_no_pat03=0.124
ctx-lh-fusiform: n_trace=2, p_no_pat03=0.534
ctx-lh-superiortemporal: n_trace=6, p_no_pat03=0.030
```

---

**End of critical review. The brief at
`.agents/reports/2026-05-06_section-5-writing-brief.md` should be edited
to reflect the lock-in summary above before any writing pass begins.**
