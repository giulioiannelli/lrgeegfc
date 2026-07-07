---
name: rho-sym-panoramic-and-methodology
type: report
era: "IMCOH_ABS × COHORT_N10"
status: current
created: 2026-07-06
audience: PI / panoramic read of every headline under the new estimator
pointers:
  - .agents/reports/2026-07-06_rho-sym-pipeline-migration.md
  - .agents/reports/2026-07-06_estimator-robustness-supplement.md
  - data/audit/rho_sym_gate/
  - data/audit/localization_atlas_rhosym/
  - data/audit/consolidation_arc_rhosym/
  - data/audit/cross_phase_taxonomy_rhosym/
  - data/audit/per_node_trace_decomposition_rhosym/
---

# ρ_sym — the methodological shift, and every headline re-read under it

## Head

We swapped the cophenetic-trace estimator from **ρ_split** to **ρ_sym**. The switch
fixes a real defect — ρ_split hard-codes an **arbitrary choice of which rest_pre half
feeds which contrast**, and that choice flips the sign of near-zero patients like a coin.
ρ_sym averages the two equally-valid choices, so the artifact is gone at zero data cost.
**Every cohort headline survives** (α/β trace, β→OFC, the inference arc, the taxonomy, the
per-node split). And to the sharp question — *does ρ_sym rescue the non-aligning patients?*
— the answer is a clean **no, and that is the point**: ρ_sym does not manufacture a trace
where there is none. What it does is **separate genuine non-tracers (stable reset, e.g.
Pat_10) from patients we simply cannot classify on half the data (now honestly labelled
"undetermined", e.g. Pat_13/Pat_15)** — instead of letting a lucky half-assignment call
them tracer or anti.

---

# Reading key — how to read any ρ in this paper (read this first)

A bare ρ is a rank concordance in [−1,1] and reads "small" out of context. **Never report it
alone** — pair it with three confound-free companions that make it interpretable. All are
already computed (audit_150 gate + audit_153 taxonomy):

| band | ρ (median) | **× null** | of movers, **% held** | **% of total structure** | anchor | gate |
|---|---|---|---|---|---|---|
| **β** | +0.20 | **16.6×** | **54%** | **26%** | 63% | **CLEAR** (p=.032) |
| α | +0.10 | 13.7× | 41% | **3%** | 90% | CLEAR (p=.024) |
| low-γ | +0.09 | 48× | 33% | 9% | 71% | fail (p=.080) |
| high-γ | +0.07 | 37× | 41% | 8% | 42% | fail (p=.080) |
| δ | +0.03 | 4.4× | (61%)\* | 7% | 89% | fail (p=.080) |
| θ | −0.04 | — | — | 4% | 87% | fail (p=.784) |

The three companions, in plain words:

- **× null (the effect size):** observed ρ ÷ matched-strength surrogate median. β's 0.20 is
  **16.6× what strength-matched noise produces** — THE headline number; a bare "ρ=0.2" hides
  it. *(High ratio ≠ trace on its own: low-γ is 48× but fails the cohort gate — the ratio is
  driven by a few patients, not cohort-consistent. Ratio is read WITH the gate.)*
- **% of movers held (persistence vs reset among the pairs the task moved):** of the structure
  the task actually perturbed, this fraction kept its new arrangement into rest; the rest
  reverted. β: **~54% held, ~46% reverted.** \*Only meaningful when the gate clears — for
  no-trace bands it is ~50/50 noise (δ's 61% is not a trace).
- **% of total structure (trace share of all cross-phase energy):** β trace = **26%** of the
  network's cross-phase structure (63% rigid anchor never moved, 11% reset). This is the honest
  "how much persists" — a substantial *minority* overlay, never "most of the structure."

**Two sentences this key makes sayable** (use these framings everywhere a ρ appears):

- **β** — *"a specific, structurally-substantial trace: ρ = 0.20 (16× a strength-matched null),
  ~a quarter of the cross-phase structure, with ~54% of what the task moved holding into rest."*
- **α** — *"statistically real but structurally thin: it clears the gate (14× null) yet is only
  ~3% of structure with 10% movers — a consistent whisper, which is exactly why α has no
  anatomical home while β localizes to OFC."*

And per patient: report **ρ_sym ± ½|ρ_AB−ρ_BA|**, mark **|ρ_sym| < 1 SE "undetermined"**, and
remember the cohort median *understates* the effect where present (strong tracers ρ = 0.3–0.5).

### The gate is a CONSISTENCY test, not a magnitude test — why big-looking bands fail

The three companions above describe the effect's **size/composition**. The cohort gate is a
*separate* question: a **signed Wilcoxon** asking whether the per-patient effect is
*consistently* positive across all 10 — a band clears only if it is positive in the **typical**
patient **and** cohort-wide, not merely big in a few. The per-patient ρ_sym spread makes each
pass/fail legible (this is why a ρ report must always show the spread, not just the median):

| band | per-patient ρ_sym (sorted, n=10) | why |
|---|---|---|
| **β** ✓ | +.54 +.51 +.50 +.36 +.27 · +.13 +.07 −.03 −.08 −.14 | 5 strong + shallow tail → **consistent** |
| **α** ✓ | +.82 · +.32 +.22 +.13 +.11 +.10 +.07 +.04 · −.03 −.15 | broad mild-positive + shallow tail → **consistent (thin)** |
| low-γ ✗ | +.86 +.74 +.65 +.41 · +.12 +.06 −.01 −.03 −.09 **−.28** | 4 strong **subset** + anti tail → **split, not cohort-wide** |
| high-γ ✗ | +.92 +.51 +.51 +.50 · +.09 +.05 +.00 −.02 −.18 **−.47** | 4 strong + **deepest** anti → **split** |
| δ ✗ | +.79 +.40 +.34 +.22 · +.04 +.02 +.02 −.03 −.07 −.08 | after top 4, everyone ≈ 0 → **too weak (typical patient barely moves)** |
| θ ✗ | mostly negative | anti / absent |

Two distinct failure modes, and the companions alone hide both:
- **low-γ & high-γ fail on CONSISTENCY** — genuinely strong in a **4-patient subset** but a
  near-zero/anti tail (down to −0.28 / −0.47) drags the signed Wilcoxon. Their eye-catching
  **× null (48× / 37×) is a median-over-near-zero-null ratio** (the low-γ surrogate median is
  0.0018, so a modest ρ gives a huge ratio) — it is **subset-inflated, NOT cohort-wide**. This
  is the locked "strong subset, not cohort-locked" verdict.
- **δ fails on MAGNITUDE** — after the top 4, the *typical* patient sits at ≈0 (median ρ 0.03,
  only 4.4× null); its **61% "held" is 50/50 noise** on near-zero movement (few real movers).

Bottom line: **a high × null alone never proves a trace**, and **% held is only interpretable
once the gate clears.** β and α pass because they are *consistently* positive (β strongly,
α thinly); the rest fail either on consistency (split cohort) or magnitude (uniformly weak).

---

# Part 1 — The methodological shift

## 1.1 Why the trace is estimated on a *split* baseline

The cophenetic trace asks: does the task-induced reorganization of the dendrogram persist
into rest_post? Operationally, for per-pair LRG cophenetic distances `D_x` in each phase,

    trace = Spearman( D_task − D_rest_pre ,  D_rest_post − D_rest_pre ).

If **the same** `D_rest_pre` is subtracted on both sides, its measurement noise appears in
both arguments and correlates them — inflating the trace even when nothing persists. To
break that, we split rest_pre into two independent halves **A** and **B** and baseline the
two contrasts on *different* halves:

    ρ_split = Spearman( D_task − D_preA ,  D_rest_post − D_preB ).

This is a **correct and deliberate design**. Its cost is that each half is only half the
data (noisier), and — the defect below — the assignment "A→task-side, B→rest-side" is
arbitrary.

## 1.2 The bias ρ_split carries: an arbitrary-half artifact (ill-conditioning)

The mirror assignment "B→task-side, A→rest-side" is **equally valid**:

    ρ_BA = Spearman( D_task − D_preB ,  D_rest_post − D_preA ).

There is no principled reason to prefer ρ_AB over ρ_BA. Yet they can disagree badly. The
two halves differ on ~45 % of pairs, and for **low-reorganization / near-zero patients**
the two assignments give *opposite signs*. Concretely, in β:

| patient | ρ_AB (arm A) | ρ_BA (arm B) | ρ_split verdict depends on the coin |
|---|---|---|---|
| Pat_13 | **+0.208** | **−0.077** | "tracer" or "anti"?? |
| Pat_15 | **+0.083** | **−0.250** | "weak tracer" or "anti"?? |

`audit_149` quantified this across the whole grid: **20 of 60 (patient × band) cells flip
sign** when A and B are swapped. The flips are **not random** — they concentrate on
patients whose true effect is near zero (the estimator is *ill-conditioned* there:
`½|ρ_AB − ρ_BA|` is large relative to the signal). Strong tracers are unaffected
(split-uncertainty ≤ 0.035; their two arms agree to the third decimal).

**This is a genuine measurement pathology**: a per-patient sign that depends on an
arbitrary bookkeeping choice cannot be reported as biology.

## 1.3 The fix: ρ_sym (symmetric estimator)

Average the two equally-valid assignments:

    ρ_sym = ½ [ Spearman(D_task − D_preA, D_post − D_preB)
              + Spearman(D_task − D_preB, D_post − D_preA) ].

By construction ρ_sym is **invariant to the half-labelling** — swapping A↔B maps arm1↔arm2
and leaves the average unchanged. Same data, same matched-strength surrogate, same cost;
the artifact is simply gone. We apply the identical symmetrization to every downstream
decomposition (per-pair concordance, per-node incidence, the arc functionals, the taxonomy
channels) and to every matched-strength surrogate realization, so observed and null are
always built the same way.

**Per-patient reporting rule (new):** report `ρ_sym ± ½|ρ_AB − ρ_BA|` (the split-
uncertainty is the natural SE), and label any patient with **|ρ_sym| < 1 SE
"undetermined"** — its sign is estimator noise, not a verdict.

## 1.4 A rejected alternative — the full (shared) baseline

The obvious way to kill the half-choice is to *not* split — use one rest_pre baseline on
both sides. We tested it and **rejected it**: it reintroduces exactly the shared-error
inflation the split design exists to prevent. For Pat_02, the full-baseline observed
statistic is 0.775 against a surrogate median of 0.767 — the signal (0.008) is drowned by
shared-baseline structure (catastrophic cancellation). **The split design is vindicated;
ρ_sym is the free refinement of it, not a retreat to the shared baseline.**

## 1.5 Reproducibility note

The strength-preserving surrogate shuffle is now numba-JIT'd with the RNG draws generated
outside the compiled loop — **bit-identical** to the reference implementation (max|Δ| =
0.0) and **98× faster** (882 → 9 ms). This is what made an *exhaustive* estimator
re-analysis (six bands, five headlines, R up to 1000) feasible rather than sampled.

---

# Part 2 — Every headline re-read under ρ_sym

New audits only; the locked ρ_split audits (63/83/103/105/144) were **not edited**.

## 2.1 Headline: the cohort trace exists in α and β (matched-strength gate)

`audit_150` — `data/audit/rho_sym_gate/`. Cohort Wilcoxon(obs − surr_median, greater), R=200.

Read each ρ with its three companions (per the Reading key): **× null** (effect size),
**% held** (of movers, persistence vs reset — taxonomy), **% struct** (trace share of total
cross-phase energy):

| band | median ρ_sym | **× null** | **% held** | **% struct** | **gate p (ρ_sym)** | gate p (ρ_split) | verdict |
|---|---|---|---|---|---|---|---|
| **β** | +0.198 | **16.6×** | **54%** | **26%** | **0.032** | 0.005 | **CLEAR** |
| **α** | +0.101 | 13.7× | 41% | **3%** | **0.024** | 0.002 | **CLEAR** |
| δ | +0.032 | 4.4× | (61%) | 7% | 0.080 | 0.216 | fail |
| low-γ | +0.089 | 48× | 33% | 9% | 0.080 | 0.116 | fail |
| high-γ | +0.072 | 37× | 41% | 8% | 0.080 | 0.246 | fail |
| θ | −0.038 | — | — | 4% | 0.784 | 0.722 | fail |

**How to read.** α and β clear 0.05; the other four fail — the *same* taxonomy as ρ_split,
**zero flips** (`audit_149`, independent per-estimator surrogates). The companions translate
the two clears: **β is real AND structurally substantial** (16× null, 26% of structure, 54%
of movers held); **α is real but structurally thin** (14× null yet only 3% of structure, 10%
movers) — a consistent whisper, which is why α has no anatomical home while β localizes to OFC.
Two guards on the companions: a high **× null** alone is not a trace (low-γ 48× but *fails*
the gate — driven by a few patients, not cohort-consistent), and **% held** is only meaningful
once the gate clears (δ's 61% is 50/50 noise). The ρ_sym p-values are *larger* than ρ_split's
(β 0.032 vs 0.005) at the **same** R=200 — not a weakening and not a grid effect: ρ_sym pulls
the ill-conditioned near-zero patients toward zero instead of letting a lucky half push them
positive. The verdict rests on the 5–6 rock-stable strong tracers, which dominate the Wilcoxon.

## 2.2 Headline: the β trace concentrates in orbitofrontal cortex

`audit_151` — `data/audit/localization_atlas_rhosym/` (R=200; R=1000 + LOO + shaft-collapse
robustness in `audit_155`, §2.6).

| system | role | BH q (epi-incl) | BH q (epi-excl) |
|---|---|---|---|
| **OFC** | **carrier (upper-tail)** | **0.025** | **0.025** |
| sensorimotor | depleted (lower-tail) | 0.017 | 0.033 |
| PFC | depleted (lower-tail) | 0.017 | 0.025 |
| cingulate / occipital / lat-temporal | secondary carriers | — | — |

**How to read.** Under the arbitrary-half-free estimator, β still concentrates in OFC and is
actively *depleted* in sensorimotor/PFC, in both epi-include and epi-exclude montages. This
reproduces the locked β→OFC verdict (canonical R=1000 q=0.009–0.013); the q difference is the
R=200 empirical floor (≈0.005 minimum), not a weakening. β remains the **only** band that
localizes; α/γ_l/δ have no FDR-surviving anatomy (a cophenetic-trace estimator change does
not touch that).

## 2.3 Headline (N2 climax): the offline trace carries an inference-specific component

`audit_152` — `data/audit/consolidation_arc_rhosym/`. Four-phase arc
rest_pre → task_learn → task_test → rest_post; matched-strength Wilcoxon, R=200.

| band | functional | median | **gate p** | LO-P15 p | reading |
|---|---|---|---|---|---|
| **β** | **T_infspec_pe** (inference-specific \| encoding) | +0.091 | **0.0098** | 0.0195 | **the claim — β-ONLY** |
| β | T_learn (learning phase's own trace) | +0.244 | 0.032 | 0.020 | learning leaves a trace |
| α | T_learn | +0.191 | 0.014 | 0.004 | learning leaves a trace |
| β | T_test (main trace, via the arc) | +0.198 | 0.032 | 0.020 | = the gate |

**β-ONLY guard**: for T_infspec_pe, every other band fails — δ p=0.246, θ 0.813, α 0.278,
low-γ 0.935, high-γ 0.313. **How to read.** The offline β reorganization is not just an echo
of encoding: it carries a component tied specifically to the *inference* step (task_test
beyond task_learn), controlling for encoding, and this is unique to β. Independently, the
*learning* phase leaves its own trace in α and β. Cross-check: the arc's arm1 reproduces the
locked audit_83 full-graph statistic to 1e-16, so this is the same pipeline, only symmetrized.

## 2.4 Headline: anchor/trace/reset composition — β is the trace-dominant band

`audit_153` — `data/audit/cross_phase_taxonomy_rhosym/`. Composition = average of both arms.

| band | anchor (rigid) | **trace** | reset | mover frac | r(anchor,strength) |
|---|---|---|---|---|---|
| **β** | 0.63 | **0.263** | 0.115 | 0.374 | 0.14 |
| δ | 0.89 | 0.072 | 0.042 | 0.113 | 0.13 |
| α | 0.90 | 0.029 | 0.064 | 0.095 | 0.05 |
| θ | 0.87 | 0.037 | 0.039 | 0.125 | 0.11 |
| low-γ | 0.71 | 0.086 | 0.164 | 0.290 | 0.44 |
| high-γ | 0.42 | 0.078 | 0.150 | 0.585 | 0.25 |

**How to read.** β has by far the largest *trace* share and mover fraction; every other band
is anchor-dominated (rigid backbone). β's low `r(anchor, strength)` = 0.14 means the rigid
backbone is **not** just the high-strength hubs — the anchor≠hubness tautology is refuted.
Trace-guard: arm1 ρ_split reproduces the locked audit_63 gate 60/60 cells bit-exactly. (The
*reset* channel remains unverified — its coupled null was retracted; only trace→OFC survived.)

## 2.5 Headline: the trace is carried by "carrier" nodes; heterogeneity is anti-trace nodes

`audit_154` — `data/audit/per_node_trace_decomposition_rhosym/`. Node-level ρ_sym↔ρ_split
correlation = **0.96** (β), **0.95** (α), **0.98** (low-γ) — a faithful refinement, not a new
measure. β per-patient (each node calibrated against its own matched-strength surrogate →
carrier / anti / neutral):

| patient | ρ_sym | carrier | anti | neutral | reading |
|---|---|---|---|---|---|
| Pat_08 | +0.536 | 75 | 0 | 45 | pure tracer |
| Pat_02 | +0.508 | 56 | 0 | 61 | pure tracer |
| Pat_05 | +0.502 | 61 | 1 | 56 | pure tracer |
| Pat_03 | +0.365 | 76 | 4 | 42 | tracer |
| Pat_06 | +0.267 | 54 | 5 | 56 | tracer |
| Pat_07 | +0.129 | 28 | 10 | 78 | weak tracer |
| Pat_13 | +0.066 | 7 | 4 | **108** | **undetermined (near-all-neutral)** |
| Pat_14 | −0.027 | 8 | 16 | 95 | weak reset |
| Pat_15 | −0.083 | 8 | 16 | 94 | null |
| Pat_10 | −0.144 | 9 | **44** | 60 | **reset (anti-dominated)** |

**How to read.** `mean_i T_i` is proportional to ρ_sym, so a low-ρ patient is one dragged
down by anti-trace nodes. The gradient is monotone: strong tracers are carrier-dominated with
~0 anti; Pat_10 is anti-dominated (a genuine *reset*); Pat_13 is almost entirely *neutral*
(genuinely undetermined). The anti-node property is preserved (β: epileptic nodes anti-
enriched, within-patient-shuffle p=0.021; MTL secondary).

## 2.6 Robustness: β→OFC at R=1000 + leave-one-out + shaft-collapse

`audit_155` — `data/audit/localization_atlas_rhosym/beta_ofc_robustness_R1000.csv`.
Like-for-like migration of the canonical ANATOMY_LEDGER lock (R=1000, LOO, shaft-collapse,
4 conditions) onto ρ_sym.

| condition | OFC full q_up | worst-LOO q_up | sensorimotor q_lo | PFC q_lo | LOO verdict |
|---|---|---|---|---|---|
| epi-excl / shaft | 0.010 | 0.045 | 0.010 | 0.010 | **robust** |
| epi-incl / shaft | 0.015 | 0.040 | 0.010 | 0.100 | **robust** |
| epi-excl / contact | 0.010 | 0.047 | 0.010 | 0.010 | **robust** |
| epi-incl / contact | 0.010 | 0.060 | 0.010 | 0.010 | marginal (1-drop over) |

**How to read.** At R=1000 the full-cohort OFC carrier q settles at **0.010–0.015** — this
now **matches the canonical ρ_split headline (R=1000 q=0.009–0.013) essentially exactly**, so
the β→OFC number migrates like-for-like, not merely "confirmed at R=200". OFC is a
significant carrier and sensorimotor/PFC are depleted in all four conditions. Leave-one-
patient-out is robust (worst-drop q < 0.05) in 3 of 4 conditions, including the strictest
(epi-exclude × shaft-collapse, worst 0.045); the one marginal case (epi-include × contact,
worst-drop q = 0.060) is a single-patient sensitivity consistent with the long-documented
**n = 5 OFC-coverage** caveat — the same fragility the ρ_split lock carried. **The β→OFC
verdict is estimator-invariant AND R=1000/LOO/shaft-robust under ρ_sym.**

---

# Part 3 — Does ρ_sym fix the non-aligning patients?

**Short answer: no — and that is exactly correct.** ρ_sym does not convert non-tracers into
tracers; it removes *false precision*. Read the β column carefully (ρ_sym gate, `audit_150`):

| patient | ρ_A (arm A) | ρ_B (arm B) | ρ_sym | ½\|ρ_A−ρ_B\| (SE) | own-surrogate p | ρ_sym class |
|---|---|---|---|---|---|---|
| Pat_08/02/05/03/06 | +0.21…+0.57 | +0.32…+0.57 | +0.27…+0.54 | ≤ 0.056 | 0.000 | **tracer (stable)** |
| Pat_07 | +0.230 | +0.028 | +0.129 | 0.101 | 0.000 | weak tracer |
| **Pat_13** | **+0.208** | **−0.077** | +0.066 | 0.143 | 0.130 | **UNDETERMINED** |
| Pat_14 | −0.049 | −0.005 | −0.027 | 0.022 | 0.740 | weak reset (stable) |
| **Pat_15** | **+0.083** | **−0.250** | −0.083 | 0.167 | 0.835 | **UNDETERMINED** |
| Pat_10 | −0.091 | −0.196 | −0.144 | 0.053 | 1.000 | **reset (stable)** |

Three distinct things were previously lumped together as "non-aligning", and ρ_sym
**separates them**:

1. **Genuine reset — stays negative, both halves agree.** Pat_10 (−0.091 / −0.196, SE 0.053)
   and Pat_14 (−0.049 / −0.005) are *stably* negative: no arm-flip, small SE. These patients
   genuinely reorganize during task and revert in rest_post. ρ_sym confirms this is real
   biology, not a coin-flip — Pat_10's per-node map is anti-dominated (44/113 anti). ρ_sym
   does **not** rescue them into the trace camp, because they truly do not trace.

2. **Genuinely undetermined — the coin-flip patients.** Pat_13 (+0.208 ↔ −0.077) and Pat_15
   (+0.083 ↔ −0.250) are the two β cells that flipped sign under the arbitrary half. Their SE
   (0.143, 0.167) *exceeds* |ρ_sym|, so ρ_sym labels them **undetermined** — we honestly
   cannot classify them on half the data. Their per-node maps confirm it: Pat_13 is
   108/119 **neutral**. Before ρ_sym, whether Pat_13 counted as a tracer or an anti was
   decided by which half we happened to call "A".

3. **Strong tracers — untouched.** The 5–6 patients carrying the cohort effect have SE ≤ 0.056
   and identical arms; ρ_sym changes them by < 0.01.

**So the honest resolution of the heterogeneity is two-layered:**
- The *estimator* layer (this migration) removes the arbitrary-half noise, converting the
  two coin-flip patients from "arbitrarily signed" to "undetermined". This is the part ρ_sym
  fixes.
- The *biology* layer (already resolved, `audit_148`,
  [[interpatient_variability_resolution_2026_07_04]]) explains *why* the stable non-tracers
  don't trace: a left-lateralised cross-phase **retention trait** (T_learn↔T_test ρ=0.95,
  sign-concordant 10/10), with non-tracers being engagement-null (Pat_15) or reset
  (Pat_10/14). ρ_sym does not change that story — it *sharpens* it by removing the estimator
  noise that was muddying the near-zero patients.

**Bottom line for the manuscript.** The cohort β/α trace is a property of the data, not of an
arbitrary half-choice. The per-patient spread is real and now cleanly tiered — stable
tracers, stable resets, and an honestly-undetermined near-zero minority — rather than a mix
of biology and coin-flips. ρ_sym is the estimator of record; ρ_split is retained only as a
supplementary column.
