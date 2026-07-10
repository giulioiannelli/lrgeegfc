---
name: cophenetic-gate-presence-vs-consistency
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-06-25
updated: 2026-06-25
pointers:
  - scripts/02_preprint/preprint_23_bands_perpair_rho_null_forest.py
  - data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv
  - scripts/01_compute/audit/audit_63_split_baseline_surrogate.py
  - .agents/preprint/locked/CONTROLS.md
  - .agents/preprint/locked/VERDICT_LEDGER.md
  - .agents/preprint/bands/03_gammalow.md
  - .agents/preprint/bands/00_cohort.md
  - data/preprint/figures/all_bands/patients_forest_plot/fig_bands_perpair_rho_null_forest_coph_sort-patient.pdf
---

> **Head.** The per-band cophenetic verdict hangs on **one number** — the C3
> one-sided paired Wilcoxon p < 0.05 — and that number answers the wrong
> question. It is a cohort-**consistency** test (does the *typical* patient
> shift in the trace direction?), but its failure is reported as **absence**
> ("no trace"). Those are not the same thing. **low-γ is the proof:** it carries
> the **same presence signature as α** (5/10 patients clear their own
> matched-strength null; binomial p = 6×10⁻⁵; obs/surr ratio 30× — *larger than
> β's*) and the **single largest per-patient trace among the candidate trace
> bands** (Pat_05 ρ=0.85 in γ_l — bigger than anything in β, whose max is 0.51,
> and tied with α's top; Pat_02 ρ=0.75 close behind). It fails the gate only
> because it is **directionally split**
> (4 patients lean anti, Pat_14 strongly), which a one-sided Wilcoxon punishes.
> The honest verdict is **"present but heterogeneous,"** not "no trace." The fix
> is to report the cophenetic trace on **two axes — presence × consistency —**
> instead of collapsing both into a single pass/fail. This is not an
> over-claim: it makes the cophenet probe finally **agree** with the already-locked
> Grassmann verdict (γ_l = `strong trace, only Grassmann`), and it leaves the one
> genuinely empty band — **θ** — exactly where it is.
>
> **The consistency key (Pat_06) is real.** Pat_06 clears its null in **all six
> bands** and is the single strongest tracer overall (mean z = 10.6, next 6.5).
> It is **not** a size or variance artifact (115 nodes, *below* cohort average;
> surrogate widths at the 70–80th percentile, i.e. *wider* than typical; ρ_split
> is a baseline-referenced double-difference, so global stability cannot inflate
> it). And the **same core five patients {02,03,05,06,08} carry both β and γ_l**
> — γ_l's clearers are an *exact subset* of β's. The trace is a real,
> patient-anchored phenomenon; bands differ in **how much of the cohort shares
> it**, not in whether the structure can carry it.

---

## 1. The trigger

`fig_bands_perpair_rho_null_forest_coph_sort-patient.pdf` (preprint_23) plots,
per band and per patient, the observed cophenetic trace ρ_split against that
patient's own matched-strength surrogate null (grey p5–p95 span; filled dot =
clears its own null; red = anti). The user's objection, looking at it:

> "γ_low seems as traced as the beta band while our headline tells that
> gamma_low has no trace … there are more 'anti patients' but surely we cannot
> assert 'no multiscale trace in gamma_low.'"

This is correct, and the figure's **own docstring already says so** (lines
21–24): it calls γ_l a **"split cohort"** (Pat_05/02/06 strongly pro,
Pat_14/15/07/13 anti → fails the gate), explicitly distinct from a **"non-trace"**
(δ, θ, γ_h, stems straddling the null both ways). The split/absent distinction
is known at the figure layer and **lost at the verdict layer.**

`ρ_split = Spearman(Δ_task, Δ_rest)`, with `Δ_task = D^tt − D^pre_A` and
`Δ_rest = D^post − D^pre_B` on **disjoint** baseline halves (audit_63). It is a
baseline-referenced double-difference — it measures whether the way task
displaced the dendrogram from baseline is preserved into rest_post. (This
matters for §6: a globally self-similar patient does **not** trivially score high.)

## 2. The numbers (ground truth: `per_patient_per_band.csv`, n=10)

| band | pos-clear (obs>p95) | neg-clear (obs<p5) | n_pro / n_anti | median ρ_split | **gate: Wilcoxon p (1-sided)** | **presence: binomial p** | obs/surr ratio |
|------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **β**   | **7** | 1 | 8 / 2 | 0.221 | **0.0049 ✓** | **4×10⁻⁷ ✓** | 23.7× |
| **α**   | **5** | 0 | 9 / 1 | 0.105 | **0.0020 ✓** | **6×10⁻⁵ ✓** | 8.3× |
| **γ_l** | **5** | 1 | 6 / 4 | 0.083 | **0.116 ✗** | **6×10⁻⁵ ✓** | **30.2×** |
| δ       | 4 | 2 | 5 / 5 | 0.008 | 0.278 ✗ | 0.0010 ✓ | 1.7× |
| γ_h     | 4 | **4** | 5 / 5 | 0.000 | 0.246 ✗ | 0.0010 ✓ | ≈0 |
| θ       | 2 | 2 | 5 / 5 | −0.040 | 0.722 ✗ | **0.086 ✗** | anti |

**Presence test** = `binomtest(pos_clear, 10, 0.05, alternative="greater")`. Under
the global null (no patient traces), each patient clears its own one-sided 95th
percentile with probability exactly 0.05, independently, so #clearing ~
Binomial(10, 0.05). It asks *"do more patients individually beat their own
connectivity-matched null than chance allows?"* — a calibrated **presence**
test, not an arbitrary count threshold.

Read the two right-hand columns against each other. **They disagree only for γ_l
and δ**, and that disagreement is the whole story: the gate (Wilcoxon) says "no,"
the presence test says "yes, overwhelmingly." γ_l's presence p (6×10⁻⁵) is
**identical to α's**, a locked trace band. Even under BH correction across all 6
bands (which the project does *not* apply), **only θ fails the presence test.**

## 3. Diagnosis — a one-axis gate doing a two-axis job

The locked decision rule (`CONTROLS.md:228`):

> **no trace** — C3 Wilcoxon p ≥ 0.05. All other controls become moot.

The C3 statistic is `wilcoxon(obs_rho − surr_median, alternative="greater")`.
A one-sided signed-rank test asks **"is the cohort *centre* shifted positive?"**
It is a **consistency** test. It necessarily fails in *two structurally opposite*
situations that it cannot tell apart:

1. **Absence** — every patient sits near its null (θ: median −0.04, nothing clears).
2. **Heterogeneity** — large effects exist but split in sign, so the rank-sum
   cancels (γ_l: median +0.08, but Pat_14 at z = −10.4 drags the signed ranks
   down; δ; γ_h, which is *symmetric* bimodal — 4 clear up, 4 clear down).

Reporting (2) as "no trace" is a **type-of-failure error**: the absence of a
*cohort-consistent* trace is not the absence of a *trace*.

**Root cause is datable.** On **2026-05-19** the project *retired* the
patient-count input to the verdict (`CONTROLS.md:237`):

> Cohort-agreement language retired … patient-count thresholds are arbitrary;
> the relevant question — *"does the per-patient distribution sit on the trace
> side?"* — is answered by the C3 paired Wilcoxon already.

The motivation (5/10 vs 6/10 thresholds *are* arbitrary) was sound, but the
remedy threw away the **presence axis** entirely. The count `n_above` *was* the
presence signal; folding it into the Wilcoxon assumed the Wilcoxon subsumes it.
For a unimodal cohort it nearly does. For a **bimodal** cohort it does the
opposite — it hides exactly the patients that carry the effect. The binomial
presence test recovers that axis **without** reintroducing an arbitrary cutoff:
its threshold is the calibrated 5% per-patient false-positive rate, and its
null is exact.

## 4. The fix — report presence × consistency, never collapse them

Replace the single pass/fail with a 2×2 read of the cophenetic substrate:

| | **consistency PASS** (Wilcoxon p<0.05) | **consistency FAIL** (Wilcoxon p≥0.05) |
|---|---|---|
| **presence PASS** (binomial p<0.05) | **cohort-wide trace** — β, α | **subgroup / heterogeneous trace** — γ_l, δ, γ_h |
| **presence FAIL** | (impossible) | **absent** — θ |

- A band is a **cohort-wide trace** only if it passes both (β strong, α modest).
- A band that passes presence but fails consistency is **present but
  heterogeneous** — the trace is real and often large, but carried by a
  subgroup rather than shared cohort-wide. This is γ_l, δ, and (symmetrically
  bimodal) γ_h.
- Only a band that fails **presence** is **absent** — θ, and θ alone.

This is the brutally-honest middle position. It does **not** promote γ_l to β's
status (γ_l is *not* cohort-consistent, full stop). It stops the false claim
that γ_l is empty.

**It also dissolves an internal contradiction.** γ_l already holds the locked
verdict `strong trace, only Grassmann ↑`. The Grassmann cluster-mass gate is, in
effect, a **presence/magnitude** test (it fires on the few strong reorganizers —
Pat_05 is logged as a "massive reorganizer," m=18, in the displacement
taxonomy). So the Grassmann probe was *already telling us γ_l is present*; the
cophenet probe's binary Wilcoxon was the only voice saying "no trace," and it
was saying so for a reason (split cohort) that the word "no trace" misrepresents.
The presence axis makes the cophenet substrate agree with Grassmann: **both
probes credit γ_l with a real, subgroup-carried trace.** The two probes then
*dissociate cleanly on geometry, not on existence*: α = distributed per-pair
trace (cophenet-present, Grassmann-absent), γ_l = low-rank subspace trace
(Grassmann-present, cophenet-present-but-split). That dissociation is a result,
not an embarrassment.

## 5. The unified band ladder (the message, per band)

Ordering by (presence, then consistency) — the honest hierarchy:

1. **β — cohort-wide trace, both probes (flagship).** 7/10 clear, 8/10 pro,
   only 1 anti. Most consistent band; survives at both LRG probes.
2. **α — cohort-wide trace, per-pair only (modest).** 5/10 clear, 9/10 pro,
   *zero* anti. Cleanest direction; distributed (no low-rank subspace signature).
3. **γ_l — present but split; the largest per-patient trace of any candidate
   band.** 5/10 clear (= α), obs/surr 30×, Pat_05 ρ=0.85 (> β's max 0.51), but
   4 patients lean anti → no cohort direction. Low-rank subspace trace
   (Grassmann-strong). **Not "no trace."**
4. **δ — present but split (weak).** 4/10 clear, balanced by 2 anti; median ≈ 0.
   Weak Grassmann-only.
5. **γ_h — bidirectional (two populations).** 4 clear up *and* 4 clear down —
   genuinely opposite-signed subgroups, net zero. Present, but with no shared
   direction at all.
6. **θ — absent.** The only band that fails the presence test (2/10, p=0.086).
   The genuine null.

The cophenetic trace is therefore **present in 5 of 6 bands** and
**cohort-consistent in 2**. The bands differ in *cohort penetrance and
directional agreement*, not in whether the LRG hierarchy can carry a task trace.

## 6. The consistency key — Pat_06 and the core subgroup

**Pat_06 is a genuine broadband super-tracer, not an artifact.**

- Clears its own null in **all six bands**; mean z = 10.6 (next: Pat_03 = 6.5).
- **Not a size artifact:** 115 nodes / 6555 pairs, *below* cohort average;
  Spearman(mean z, N_nodes) = −0.02.
- **Not a variance artifact:** its surrogate widths sit at the **70–80th
  percentile** in δ/θ/α (i.e. *wider* than typical), yet obs_rho ≈ 0.78–0.83.
  The large z comes from a large numerator, not a small denominator.
- **Not a global-stability artifact:** ρ_split is a baseline-referenced
  double-difference on disjoint baseline halves, so a merely reproducible
  recording does not score high.

**The positive mechanism — implant coverage of the OFC hotspot (not an
artifact, a localization).** Ruling out artifacts is necessary but not an
explanation. The explanation is anatomical: **Pat_06 has the heaviest OFC
implant in the cohort — 19 OFC contacts (16.5 % of its montage) vs a cohort
median of 0**, plus 60/115 (52 %) frontal. Its trace *concentrates* in OFC
(frac_trace_raw = 0.64) and cingulate (0.59) — **the known β-trace hotspot**
(audit_83/92, OFC q≈0.01). Pat_06 lights up across every band because it
densely samples the tissue that carries the trace **and** is a strong
consolidator. This ties the broadband super-tracer directly to the localization
result rather than leaving it as an unexplained outlier.

**But coverage gates *visibility*, not *sign* — the honest twist.** The cohort
correlation between OFC/frontal coverage and trace strength is only weak-positive
(Spearman ≈ +0.19), and one patient breaks it cleanly: **Pat_14 has the
*second*-heaviest OFC/frontal implant (15 OFC, 45 frontal) yet is the cohort's
*worst* tracer** (mean z = −3.2, strongly anti in γ_l/γ_h). So OFC coverage is
roughly *necessary to see* an OFC trace (the temporal-only patients Pat_03/13
have 0 OFC and trace weakly) but **not sufficient** — among the patients who
sample OFC, whether they trace or anti-trace is a separate consolidator-vs-not
biology. Pat_06 = heavy OFC coverage **×** strong-consolidator. Pat_14 = heavy
OFC coverage **×** non-consolidator.

**What it means.** A single patient *can* carry a strong, reproducible,
multiscale reorganization that persists across every band, and **where** it sits
(OFC) is the same hotspot the cohort localization finds. Combined with the
**moderate per-patient coupling** across bands (mean band-band Spearman of obs_z
≈ 0.29), this says the trace is a **patient-anchored, anatomically-localized
phenomenon** — patients who trace in one band tend to trace in others, in the
same tissue — and that band-level cohort verdicts are about *how many patients
share a band's trace*, not about which bands are capable of carrying one.

**The core-subgroup nesting is the cleanest evidence γ_l is real:**

- β clearers: {02, 03, 05, 06, 07, 08, 13}
- γ_l clearers: {02, 03, 05, 06, 08} — an **exact subset** of β's.

The *same five core consolidators* carry both β and γ_l; β merely recruits two
extra marginal patients (07, 13). γ_l's trace is not a different, noisy
phenomenon — it is **the same trace in the same patients, just not extended to
the cohort margin.** (α's clearers {03,06,08,10,14} overlap only partially —
consistent with α tapping a somewhat different population and the α/γ_l geometry
dissociation in §4.)

**Honest caveat — Pat_06 does NOT rescue θ.** Pat_06 traces strongly even in θ
(z=11.3), where the *next* patient is z=1.99 and does not clear. So Pat_06's
trace is **broadband / frequency-nonspecific** — a per-patient property, not
evidence of a θ cohort trace. θ stays the genuine null. This is the same brutal
honesty the gate fix demands in the other direction: a single super-responder is
not a cohort effect, just as a split cohort is not an absence.

## 7. What this does and does not license

**Does:**
- Retire the phrase "no trace" for **γ_l, δ, γ_h** at the cophenet probe;
  replace with "**present but heterogeneous / no cohort-consistent direction**"
  (γ_l, δ) and "**bidirectional**" (γ_h).
- Add a **presence axis** (binomial on own-null clearing) to the cophenet
  verdict, reported beside the consistency Wilcoxon.
- State the **core-consolidator subgroup** {02,03,05,06,08} and Pat_06 as a
  positive consistency result.

**Does not:**
- Promote γ_l (or δ, γ_h) to a **cohort-wide** trace. They are *not*
  directionally consistent. β remains the only band strong on both axes and both
  probes; α the only other cohort-consistent band.
- Touch β's or α's verdicts (both already pass both axes).
- Rescue **θ** — it fails presence and stays the null.
- Introduce any new p-value as a *gate*. The binomial presence test is reported
  as **a count with its chance baseline**, not a new pass/fail hurdle. The
  verdict is the 2×2 cell, read off two transparent numbers (how many cleared;
  did they share a direction).

## 8. The figure fix

The forest figure currently colours each band letter **binary** — green if the
Wilcoxon gate passes, grey otherwise. That binary *is* the one-axis collapse: it
paints γ_l grey, which reads as "nothing," contradicting its own dots. Update
(implemented in `preprint_23`, `--verdict two-axis`):

- **Band letter, 3-way colour** by the 2×2 cell: **green** = cohort-wide trace
  (β, α); **amber** = present-but-heterogeneous (γ_l, δ, γ_h); **grey** = absent
  (θ).
- **A presence/consistency glyph** beside the letter encoding both axes at a
  glance (filled wedge = presence ✓; full vs split fill = consistency).
- Counts (pos-clear / neg-clear), binomial p and Wilcoxon p stay in the
  **caption / stdout** per the no-in-axes-text rule — but now both numbers are
  reported, not just the gate p.

The per-patient forest (null span + filled/open/red dots) is already honest and
unchanged; the only thing that was lying was the binary band-letter colour.

## 9. Cascade across the preprint

Load-bearing (must change to remove the contradiction):

- **`locked/CONTROLS.md:228-238`** — the decision rule. Add the presence axis;
  redefine "no trace" as **presence-fail**, and add a
  "**present, not cohort-consistent**" tier for presence-pass / consistency-fail.
  This is the keystone; everything else follows.
- **`locked/VERDICT_LEDGER.md`** — γ_l (`:215-221`), δ (`:324`), γ_h (`:~345`)
  cophenet rows: "no trace" → "present, not cohort-consistent (k/10 clear,
  binom p, Wilcoxon p)". β/α unchanged. θ unchanged ("absent").
- **`bands/03_gammalow.md:58-64,113-124`** — rewrite the C3 reading from "no
  trace at the cophenet probe" to "present-but-split"; note γ_l's cophenet
  substrate now *agrees* with its Grassmann verdict. (This file already contains
  the correct mechanism at `:64` — "the cophenet step disperses γ_l per-pair
  structure" — it just mislabels the consequence.)
- **`bands/00_cohort.md:29,41,72-76`** — the "cophenet does band resolution"
  thesis must be restated: the cophenet step resolves **consistency**, not
  **presence**. It demotes δ/θ/γ_h *from cohort-consistency*, and θ from
  presence; it does not erase the γ_l/δ/γ_h subgroup traces. The three-layer
  table γ_l row (5/10, p=.12) gains a presence column.
- **`HANDOFF_INDEX.md:116`** — γ_l cheatsheet line.

Secondary (propagate the wording): `README.md:62`,
`directives/writing_directive_2026-06-04_grassmann-subsection-results.md:70-82`
("do NOT promote γ_l to a trace band" → "γ_l is a subgroup/heterogeneous trace,
not a cohort-wide one"), `directives/methods_displacement_taxonomy_2026-06-03.md:91`,
`responses/2026-05-26_preprint_folder_audit.md`, `EVALUATION_PROTOCOL.md:119`
(presence axis is per-band, still no cross-band BH).

Unaffected: the localization six-rung gate (β→OFC), the encoding/inference arc
(N2), the low-γ *focal cingulate* encoding trace (a region-level result,
orthogonal to the whole-brain band verdict), all anatomy ledgers.

## 10. Open decisions for the PI

1. **Adopt the two-axis verdict?** (presence × consistency, replacing the binary
   "C3 p<0.05 else no trace"). Recommended — it is strictly more honest and
   removes the γ_l cophenet/Grassmann contradiction.
2. **Naming.** "present but heterogeneous" vs "subgroup trace" vs "split trace"
   for the presence-pass / consistency-fail tier. (Report uses these
   interchangeably; pick one for the ledger.)
3. **How hard to feature the core-consolidator subgroup** {02,03,05,06,08} +
   Pat_06 in the manuscript — a genuine new consistency narrative, or a
   supplementary observation.
4. **Execute the cascade now or stage it.** The locked-ledger edits are real
   scientific changes; this report stages them but does not apply them.
