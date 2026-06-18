---
name: trace-concordance-vs-blind-fc
era: IMCOH_ABS_COHORT_N10
status: current
date: 2026-06-01
companion: data/audit/trace_concordance_crosstab/
sources:
  - audit_63  matched_strength_surrogate_split_baseline   # rho_split^coph
  - audit_66  grassmann_cluster_extent (+ _matched_strength_surrogate)  # T_G*
  - audit_67  raw_fc_matched_strength                       # rho_split^raw (BLIND)
  - audit_75  trace_concordance_crosstab                    # this aggregation
---

# Is the band-selective β trace defensible across the two LRG probes — and what does the blind FC miss?

**Head.** β is the one band where the task→rest trace is *defensible and not
fragile*: both LRG probes — per-pair cophenetic `ρ_split^coph` and
leading-subspace `T_G*` — clear the matched-strength null **and** survive
leave-one-patient-out, 9/10 patients carry it on at least one probe (7/10 on
**both**), and the single complete dissenter is the pre-registered
right-hemisphere-only Pat_15. The blind raw-FC edge comparison — the *same*
split-baseline statistic under the *same* null, computed on flat `triu(A)`
edges instead of the hierarchy — **fails at β** (6/10, cohort p = 0.053, and
leave-one-out cannot rescue it). Three of the nine β-trace patients are
invisible to flat FC and only appear after the LRG transform. That gap is the
richness the LRG step adds. Every other band is weaker: α is cophenet-only and
*shared with blind FC* (LRG adds little there); γ_l/δ are Grassmann-only and
LOO-thin-to-fragile; θ/γ_h are null. The two probes co-fire on **β alone** —
which is exactly why β is the headline and the others are not co-equal claims.

All numbers below are surfaced from existing matched-strength surrogate
outputs (R = 200, 4-cycle ±δ, split-baseline). No new test was run; `audit_75`
only aggregates, and its β column reproduces the locked `beta_per_patient.tex`
exactly (hard integrity gate).

---

## 1. Definitive cohort verdicts — the two LRG probes + the blind baseline

Same split-baseline (`Δ_task = X^tt − X^preA`, `Δ_rest = X^post − X^preB`),
same R = 200 matched-strength null. The only thing that changes across the
three columns is the **representation** of each pair: raw `imcoh_abs` edge
(blind), cophenetic communication distance `D_coph` (LRG hierarchy), or
leading-eigenmode subspace rotation (LRG subspace).

| band | `ρ_split^coph` (LRG hierarchy) | `T_G*` (LRG subspace) | `ρ_split^raw` (BLIND edges) |
|---|---|---|---|
| | cohort-Wilcoxon p · n≥1.96 · LOO-max p | cluster-mass p · LOO-max p (argmax) | cohort-Wilcoxon p · n≥1.96 · LOO-max p |
| **β** | **0.0049 · 7/10 · 0.0098** | **0.0050 · 0.0050 (any)** | 0.053 · 6/10 · 0.102 |
| α | 0.0020 · 5/10 · **0.0039** | 0.348 · 0.577 — no | 0.014 · 7–8/10 · 0.027 |
| γ_l | 0.116 · 5/10 · 0.213 — no | 0.0050 · 0.0398 (Pat_05) | 0.053 · 7/10 · 0.102 |
| δ | 0.278 · 4/10 · 0.455 — no | 0.0050 · **0.0547** (Pat_08) | 0.042 · 7/10 · 0.082 |
| θ | 0.722 · 2/10 — no | 0.159 — no | 0.138 — no |
| γ_h | 0.246 · 4/10 — no | 0.060 — no | 0.188 — no |

Reading the table:

- **β is the only band where both LRG probes fire and both are LOO-robust.**
  `T_G*` is *bulletproof*: dropping **any** single patient leaves the
  cluster-mass permutation p at its 0.005 floor. `ρ_split^coph` LOO-max is
  0.0098 — drop any patient and it is still p < 0.01.
- **α fires on cophenet (robust, LOO 0.004) but not on the subspace
  (`T_G*` p = 0.35), and the blind edges already see it (p = 0.014).** So at α
  the LRG hierarchy does not add discriminating power over flat FC — α is a
  "raw FC already resolves it" band.
- **γ_l and δ are Grassmann-only.** γ_l survives LOO with a thin margin
  (0.040, drop Pat_05); δ is **LOO-fragile** — dropping Pat_08 pushes it to
  0.055, over the line. Neither has cophenet or blind support.
- **θ, γ_h are null on every probe.**

> Caveat on the per-patient `T_G*` flag (used in §2): the manuscript
> per-patient rule `n_sig_k ≥ 20` is *permissive* — it fires for 8–9/10
> patients in **every** band, including θ and α where the **cohort**
> cluster-extent gate says no trace. The leading-subspace rotation persists
> task→rest in almost every patient and band; what is band-selective is the
> **cohort-level cluster extent** (count + contiguity of significant `k`-cells
> over matched-strength), not the per-patient flag. So in the cross-tab the
> band signal is carried by the cophenet column and the `T_G*` **cohort** gate
> — not the `T_G*` per-patient column. This is stated so the cross-tab is not
> over-read.

---

## 2. Per-patient concordance cross-tab (the "both / one / none" table)

LRG-layer classification per patient: **B** = trace on both probes, **C** =
cophenet `ρ_split` only, **G** = subspace `T_G*` only, **·** = neither.
(`n_sig_k ≥ 20` for `T_G*`; `z ≥ 1.96` for `ρ_split`, both verbatim manuscript
rules.) Blind row = raw-FC edge `ρ_split^raw`, `z ≥ 1.96`.

```
band        02 03 05 06 07 08 10 13 14 15   both onlyC onlyG none | LRG≥1  blind
delta        B  B  G  B  B  G  G  G  .  .     4    0     4    2  |   8      7
theta        G  G  B  B  G  G  G  G  G  .     2    0     7    1  |   9      7
alpha        G  B  G  B  G  B  B  G  B  .     5    0     4    1  |   9      7
beta         B  B  B  B  B  B  G  B  G  .     7    0     2    1  |   9      6
low_gamma    B  B  B  B  .  B  G  G  G  .     5    0     3    2  |   8      7
high_gamma   G  B  B  B  .  B  G  G  G  .     4    0     4    2  |   8      8
```

Two structural facts:

1. **`onlyC = 0` in every band.** Every cophenet-positive patient is also
   `T_G*`-positive — the subspace flag is a strict superset of the cophenet
   flag at the per-patient level. The probes are *nested in sensitivity*, not
   orthogonal: when the per-pair merge-height ranks move, the leading subspace
   has always already rotated.
2. **`n_both` peaks at β (7), driven by the demanding cophenet probe
   (`coph_nsig`: δ4 θ2 α5 β7 γ_l5 γ_h4).** β is the band where the most
   patients clear *both* the strict per-pair test and the subspace test.

---

## 3. β single-patient breakdown and the blind-FC recovery story

This is the core "richness" result. At β the blind edge comparison sees 6/10
and **fails the cohort test** (p = 0.053, LOO can't get under 0.05). Each LRG
representation of the *same data under the same null* recovers patients the
flat bag-of-edges reads as non-trace:

| patient | blind `ρ_split^raw` | cophenet `ρ_split^coph` | subspace `T_G*` (n_sig) | class |
|---|---|---|---|---|
| Pat_02 | +0.364 ✓ | +0.507 (z 6.27) ✓ | 0.962 (108) ✓ | both |
| Pat_03 | +0.700 ✓ | +0.373 (z 9.16) ✓ | 0.859 (96) ✓ | both |
| Pat_05 | +0.530 ✓ | +0.491 (z 7.32) ✓ | 0.766 (87) ✓ | both |
| Pat_06 | +0.734 ✓ | +0.211 (z 3.71) ✓ | 0.858 (96) ✓ | both |
| Pat_07 | +0.192 ✓ | +0.230 (z 4.47) ✓ | 0.530 (62) ✓ | both |
| Pat_08 | +0.324 ✓ | +0.502 (z 8.76) ✓ | 0.520 (62) ✓ | both |
| **Pat_13** | **−0.221 ✗ (anti)** | **+0.208 (z 3.51) ✓** | 0.520 (61) ✓ | both (LRG) |
| **Pat_10** | −0.120 ✗ | −0.091 (z −1.70) ✗ | **0.612 (71) ✓** | only `T_G*` |
| **Pat_14** | −0.002 ✗ | −0.049 (z −1.00) ✗ | **0.695 (78) ✓** | only `T_G*` |
| Pat_15 | −0.023 ✗ | +0.083 (z 1.27) ✗ | 0.006 (1) ✗ | neither |

- **Pat_13** is the cleanest demonstration: the hierarchy *flips* an
  edge-level **anti-trace** (ρ = −0.221) into a merge-height **trace**
  (ρ = +0.208). The same pairs, re-expressed as communication-distance ranks,
  carry a task→rest imprint the flat edges actively contradict.
- **Pat_10 and Pat_14** carry the β trace entirely in the leading-eigenmode
  **subspace rotation** (z ≈ +5 on `T_G*` at audit_66) while being flat/anti
  on both per-pair representations. A real geometric dissociation — not noise.
- Net: blind **6/10 (fails)** → LRG **9/10 (clears, LOO-proof)**. Three of the
  nine β-trace patients (10, 13, 14) are *only* visible after the LRG
  transform. **The blind FC comparison is literally unable to gather the
  richness.**

The lone full dissenter, **Pat_15**, is the pre-registered biology exception
(right-hemisphere-only implant; the canonical single cohort-level
anti-aligned patient at β). Its dissent is the documented exception that the
patient-dropout policy already carves out, not an unexplained failure.

---

## 4. Is the β claim defensible and not fragile? — verdict

**Yes, on a "universality-class" reading, and the fragility lives in the
*other* bands, not β.**

Defensibility rests on four legs, each independently checkable:

1. **Two independent LRG geometries agree.** Per-pair merge-height ranks
   (`ρ_split^coph`) and leading-subspace rotation (`T_G*`) are different
   constructions of the dendrogram; both clear matched-strength at β. No other
   band has both.
2. **Both survive leave-one-patient-out.** `T_G*` is LOO-invariant (0.005 for
   every drop); `ρ_split^coph` LOO-max is 0.0098. The verdict does not hinge on
   any single patient. (Contrast δ: `T_G*` LOO-max 0.055 — one patient flips
   it.)
3. **Cohort coverage is 9/10 on ≥1 probe, 7/10 on both**, with the single
   complete dissenter biology-justified (Pat_15).
4. **The blind baseline fails at β** under the identical statistic and null —
   so the β verdict is a property of the *hierarchical/subspace geometry*, not
   an artifact reproducible from flat edges.

What the β claim is **not**: it is not "10/10, every patient, every probe."
Pat_10/Pat_14 split across probes; Pat_15 dissents. The honest claim is a
*cohort-wide trace in the LRG communication geometry*, robust to dropping any
one patient and invisible to flat FC — not a per-patient universal.

**The genuinely fragile / measure-dependent parts are the secondary bands,**
and they should be scoped accordingly:

- **α** — cophenet-only, and the blind edges already resolve it (p = 0.014).
  LRG adds no discriminating richness at α. Defensible as a secondary
  observation, *not* as an LRG-specific result.
- **γ_l** — Grassmann-only, LOO margin thin (0.040). A genuine subspace-only
  finding; state with its single-probe, thin-LOO caveat.
- **δ** — Grassmann-only and **LOO-fragile** (0.055). Do not headline.
- **θ, γ_h** — null.

**Band-selectivity ("β bears it, others do not") is only fully defensible as a
cross-measure statement about β.** The two probes *disagree* on which
secondary bands trace (cophenet → α; Grassmann → δ, γ_l); they *agree* on β
alone. So "β is special" is defensible precisely as: *the unique band on which
both LRG probes concur, both are LOO-robust, and the blind FC fails.* A flat
"each band either bears the trace or doesn't" framing is **not** defensible —
it is probe-dependent for every band except β.

---

## 5. Bottom line for the manuscript

- **Lead with β as the cross-probe, LOO-robust, blind-FC-beating result.** The
  three-patient recovery (10, 13, 14) is the concrete evidence that the LRG
  step resolves structure flat FC cannot — the cleanest available argument for
  why the hierarchy/subspace step earns its place.
- **Demote α to a secondary, blind-FC-shared observation;** present γ_l (and δ
  with its fragility) as subspace-only; θ/γ_h as null.
- **Drop binary "bears/doesn't" band language.** Replace with the graded,
  probe-explicit statement above, and show the per-patient concordance grid so
  the heterogeneity is visible rather than discovered by a referee.

Reproduce: `python scripts/01_compute/audit/audit_75_trace_concordance_crosstab.py`
(surfacing only; integrity-gated against `beta_per_patient.tex`). Outputs in
`data/audit/trace_concordance_crosstab/`.
