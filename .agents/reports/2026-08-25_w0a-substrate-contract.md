---
name: 2026-08-25_w0a-substrate-contract
kind: report
era: PAPER_FINALIZATION (Wave 0, lane W0-A)
status: complete
created: 2026-08-25
scope: Evidence report for the W0-A substrate contract. Decides the two choices that underpin every result in the paper -- the |ImCoh| magnitude transform and the sparsified backbone -- under a decision rule pre-registered before any number was computed. Reports the full sparsification stability surface, the parameter-free-filter verdict, and an explicit list of what remains unsettled.
pointers:
  - scripts/01_compute/paper_final/w0a_00_preregistration.md
  - .agents/preprint/locked/PIPELINE_CONTRACT.md
  - src/lrg_eegfc/workflow/substrate.py
  - data/paper_final/w0a_substrate/
  - .agents/plans/active/2026-08-25_paper-finalization-master-plan.md
---

# W0-A — the substrate contract

## Head

The substrate is `imcoh_abs` on an mst-union backbone **reported over `f ∈ [0.07, 0.20]` rather than at any single fraction**, and the free parameter could not be eliminated — no parameter-free filter passed, so it is integrated over instead. Three things came out against expectation and all three matter. **α is recoverable**: it holds an invariantly-positive plateau over 2.32 octaves, which reverses the demotion sibling lane W0-C reached using a per-scale multiplicity correction that over-charges the correlated scale axis by ~11×. **The transform is not a free choice**: the reliability rule preferred `imcoh_sq`, but β's trace does not survive it (positive at 6/10 swept configs under `abs`, 3/10 under `sq`) and `sq` admits no stable window at all — so `abs` is kept because it is the only admissible option, not because it won. **The full band-selectivity claim is narrower than the trace claim**: adding low_γ to the required-null set shrinks the window to 0.51 octaves, below the pre-registered bar. The historical `f = 0.20` sits at the upper edge of every window computed here, which is the quantitative signature of an outcome-selected knob; the contract moves the representative fraction to the interior. Everything remains conditional on matched-strength, a null injected at the finished-FC-matrix stage — lane W0-B's ladder is what would remove that, and until it runs these verdicts are provisional in that specific sense.

---

## 0. What this lane was asked to settle, and the firewall around it

Two choices underpin every multiscale result in the paper and neither was decided by evidence. `imcoh_abs = ⟨|Im C|⟩_f` was chosen on 2026-04-15 by argument (a genuine Jensen order-of-operations fix plus a citation to Ewald 2012 / Bastos & Schoffelen 2016), with no head-to-head comparison against the equally citable `imcoh_sq = ⟨(Im C)²⟩_f`. The backbone fraction `f = 0.20` was chosen because, in the methods guide's own words, "only 0.20 gives α strong + β maximal + θ/low_γ null" — selecting the knob so the answer comes out right, and failing that same guide's own criterion 3 ("parameter-free"), which lists `mst@f` in its method table as "reject (free parameter f)".

The decision rule for both was written and frozen in `scripts/01_compute/paper_final/w0a_00_preregistration.md` **before any number was computed**, with the mandatory 5-point critical preamble. It selects on hypothesis-independent properties only: estimator reliability, distance from a lag-destroying null floor, weight heterogeneity, and invariance of the verdict under the knob. Trace gates and marker AUCs under each option are computed and recorded but are **not admissible** at any step of the rule. The rule was not relaxed after the numbers arrived; one dated addendum (R1-bis) was added after A1 and before any A2 trace number, and it is flagged as such in place.

Everything here is conditioned on **matched-strength** being the null — a null injected at the finished-FC-matrix stage, which therefore cannot test the coherency estimator, the band split, session nonstationarity, or drift. That limitation is lane W0-B's to remove and is not removed here.

---

## 1. A1 — the transform, head to head

### 1.1 What was computed

Every number in this section is recomputed fresh from the raw timeseries, not sourced from any prior report or CSV. For each of the 10 patients: one Welch pass per `rest_pre` half (at `nperseg_for_fs(fs) // 2`, the incumbent split-half setting) and one at full duration and full `nperseg`, giving the complex band coherency `C(f)` per band. From the same object both transforms are derived per frequency bin and then band-averaged — the order the canonical loader uses.

**Reproduction check.** The `abs` arm reproduces the existing `imcoh_halves_fc` cache to `max|diff| ≤ 3.0e-8` for every patient (float32 round-trip of the incumbent files), so this is the same estimator, not a re-implementation that drifted.

### 1.2 Step 1 — the rank-equivalence gate: the transform is NOT a no-op, but its effect is bounded

The brief's prior was that squaring is a monotone map on non-negative weights, so any rank-based backbone should be invariant to it. **That is false here, and the reason is instructive**: the transform is applied per frequency bin and the band average is taken afterwards, so `mean(|x|)` and `mean(x²)` are not related by a monotone map — Jensen's inequality again, the same mechanism as the 2026-04-15 reset.

| band | Spearman(abs, sq), median | mst-union edge-set Jaccard @ f=0.05 | @ f=0.10 | @ f=0.20 |
|---|---|---|---|---|
| delta | 0.987 | 0.823 | 0.847 | 0.862 |
| theta | 0.987 | 0.782 | 0.823 | 0.846 |
| alpha | 0.986 | 0.784 | 0.808 | 0.834 |
| beta | 0.990 | 0.828 | 0.854 | 0.882 |
| low_gamma | 0.987 | 0.896 | 0.903 | 0.881 |
| high_gamma | 0.990 | 0.813 | 0.867 | 0.891 |

Across all 60 (patient, band) cells the edge-rank Spearman spans `[0.882, 0.998]` and the backbone edge-set Jaccard at `f = 0.20` spans `[0.510, 0.987]`. **So the transform moves between 1 % and 49 % of the retained backbone edges depending on the cell, with a cohort median of 12–17 %.** The gate does not fire; the choice is not vacuous. But the effect is bounded, and that bound — a rank correlation never below 0.88 and a typical backbone overlap of ~0.86 — is itself the cleanest thing to quote about the transform question, in either direction.

### 1.3 Step 2 — split-half reliability: the rule selects `imcoh_sq`, by a small margin

Spearman is the primary statistic because it is invariant under any strictly monotone map, so a difference can only come from the band-averaging order and not from units. Pearson is reported as secondary and was not used to decide.

| band | Spearman rel. abs | sq | Δ(abs−sq) | winner | [Pearson abs / sq] |
|---|---|---|---|---|---|
| delta | 0.5950 | 0.5971 | −0.0021 | sq | 0.820 / 0.875 |
| theta | 0.6585 | 0.6652 | −0.0067 | sq | 0.795 / 0.821 |
| alpha | 0.6248 | 0.6399 | −0.0151 | sq | 0.776 / 0.816 |
| beta | 0.7313 | 0.7402 | −0.0089 | sq | 0.845 / 0.889 |
| low_gamma | 0.5481 | 0.5587 | −0.0106 | sq | 0.848 / 0.864 |
| high_gamma | 0.6297 | 0.5852 | **+0.0445** | **abs** | 0.701 / 0.634 |

Cohort over the 60 cells: median Δ = −0.0078, Wilcoxon two-sided p = 5.0e-6, band wins abs 1 / sq 5. The pre-registered decisiveness condition (≥ 4/6 band wins **and** p < 0.05) is met, so **R1 stops at step 2 and selects `imcoh_sq`**.

Two honest qualifications, neither of which changes the rule's output. First, the margin is ~1.2 % relative on reliabilities of 0.55–0.74 — systematic, but small. Second, `high_gamma` is the one band where `abs` wins, and it wins by five times the median margin.

**Post-hoc diagnostic (not part of the rule).** The strongest alternative reading of a reliability win is that it is inherited from coherence *magnitude* — anatomy and volume conduction, which are highly reproducible across halves and have nothing to do with lag estimation. It is not: `Spearman(observed edges, band-averaged |C|)` differs between the two transforms by ≤ 0.005 in five of six bands, and the one appreciable difference (high_gamma, +0.055) is in the band where `abs` wins reliability. The diagnostic can only weaken the step-2 conclusion and does not.

### 1.4 Step 3 — the circular-shift floor, and a polarity warning that matters for lane W0-B

A circular time shift of channel `c` by `t_c` maps `C_ij(f) → C_ij(f)·e^{−2πif(t_i−t_j)}`: it preserves `|C_ij(f)|` — hence every bit of volume conduction and of the coupling magnitude — exactly, and randomises only the phase, i.e. the lag. Thirty-two draws per (patient, band), frequency-domain phase-ramp form.

**The observed `|ImCoh|` sits far BELOW its lag-randomised floor, in 120 of 120 cells.** Median `obs/floor` is 0.13–0.32 for `abs` and 0.02–0.11 for `sq`; median per-edge `z` runs from −3.3 (delta) to −96 (high_gamma); the fraction of edges above their own floor is 0.3–5 %. This is expected from the mechanism, not a pathology: under uniform random phase `⟨|Im C|⟩ → (2/π)·⟨|C|⟩`, and the true phase distribution is concentrated near zero lag (which is what volume conduction does), so destroying the lag structure *inflates* `|ImCoh|` rather than deflating it.

The consequence is a design warning, and it is the most transferable thing in this section: **"distance above the circular-shift floor" is the wrong polarity for a per-edge comparison and must never be reported as a null-clearing statistic.** The floor comparison has to be framed as a *structure* test — does the surrogate reproduce the FC *pattern*? — rather than a magnitude test. Framed that way it is informative and reassuring: `Spearman(observed edges, floor edges)` has a cohort median of −0.13 to +0.09 in five bands (high_gamma −0.26 / −0.31), so the observed edge ordering is essentially **uncorrelated** with the magnitude-driven floor ordering. Neither transform is separated by this, so step 3 does not arbitrate — and it would not have been reached anyway, since step 2 was decisive.

### 1.5 Step 4 — weight heterogeneity (recorded; step 2 had already decided)

| band | Gini abs / sq | CV abs / sq | participation ratio abs / sq |
|---|---|---|---|
| delta | 0.345 / 0.606 | 0.70 / 1.47 | 0.670 / 0.317 |
| theta | 0.292 / 0.516 | 0.55 / 1.24 | 0.767 / 0.394 |
| alpha | 0.301 / 0.531 | 0.60 / 1.29 | 0.734 / 0.374 |
| beta | 0.271 / 0.514 | 0.57 / 1.49 | 0.756 / 0.310 |
| low_gamma | 0.216 / 0.457 | 0.53 / 1.47 | 0.783 / 0.316 |
| high_gamma | 0.122 / 0.255 | 0.26 / 0.67 | 0.937 / 0.692 |

`imcoh_sq` roughly **doubles** the edge-weight Gini in every band and drops the participation ratio from ~0.75 to ~0.32 — i.e. the effective fraction of edges carrying the weight falls by more than half. On a dense graph, weight heterogeneity is precisely what drives the propagator into the Villegas single-peak regime, so this is a real cost of the transform the rule selected. The rule is lexicographic and step 4 is a final tie-break, so it does **not** override step 2 — but it is why R1-bis exists.

### 1.6 R1-bis — the dated addendum

Written after A1 and **before any A2 trace number**: a 1.2 % reliability margin is not a sufficient basis on which to regenerate every artifact in the project, when the same switch moves 10–22 % of the backbone and doubles the weight heterogeneity. So A2 is run under **both** transforms. If the per-band, per-scale verdict is invariant across the transform inside the plateau, the contract records the transform as **immaterial within the plateau** and retains `imcoh_abs` for continuity, citing the demonstrated invariance — not the reliability comparison — as the justification. If the verdict is not transform-invariant, `imcoh_sq` is adopted exactly as R1 dictates and every downstream artifact is regenerated on it.

---

### 1.7 `⟨|C|⟩` vs `|ImCoh|`, with and without same-shaft pairs — the control W0-B's result forced

**Nothing in this section reaches significance at n = 10.** Every number is directional and is reported as such.

W0-B found that the trace fails every lag-destroying null and that ordinary coherence `⟨|C|⟩` gives a *larger* trace than `|ImCoh|` in both trace bands (β +0.355 vs +0.281, α +0.277 vs +0.212, five patients). If that held, the transform this project adopted for volume-conduction immunity would be the weaker one on its own headline statistic. But `⟨|C|⟩` is not VC-immune and carries the same-shaft bias of CLAUDE.md invariant 5, and anatomically fixed structure is trivially stable across phases — so the larger number might be the artifact. Zeroing the within-shaft block separates the two.

Design: three transforms × two pair sets × three backbone fractions from the contract window × five phases × ten patients, all arms sharing **one Welch pass per phase** so the spectral estimate is bit-identical and only the transform differs. Same-shaft removal applied at **both** stages — edges zeroed before sparsification, and the pairs dropped from the ρ_sym correlation so they cannot re-enter through indirect cophenetic paths.

**Same-shaft inflation, confirmed and quantified** (mean same-shaft / cross-shaft edge weight, task_test):

| band | `imcoh_abs` | `⟨\|C\|⟩` |
|---|---|---|
| θ | 1.01× [0.73, 1.21] | 1.42× [1.16, 2.42] |
| α | 1.05× [0.73, 1.42] | 1.51× [1.26, 2.33] |
| β | 1.24× [0.90, 1.60] | 1.70× [1.32, 2.67] |

`⟨|C|⟩` is consistently the more shaft-inflated measure, in the expected direction.

**`T_probe`, cohort medians over scales and fractions** (obs, and margin vs matched-strength):

| band | transform | all pairs | same-shaft removed | change |
|---|---|---|---|---|
| θ | `imcoh_abs` | +0.004 (m −0.002) | +0.024 (m +0.019) | +0.020 |
| θ | `⟨\|C\|⟩` | **+0.133** (m +0.150) | **+0.115** (m +0.121) | −0.018 |
| α | `imcoh_abs` | +0.138 (m +0.132) | +0.109 (m +0.101) | −0.029 |
| α | `⟨\|C\|⟩` | +0.228 (m +0.214) | +0.137 (m +0.130) | −0.091 |
| β | `imcoh_abs` | +0.171 (m +0.153) | +0.177 (m +0.160) | +0.006 |
| β | `⟨\|C\|⟩` | +0.122 (m +0.114) | +0.182 (m +0.163) | +0.061 |

Per-patient paired tests (n = 10), advantage of `⟨|C|⟩` over `|ImCoh|`:

| band | all pairs | same-shaft removed |
|---|---|---|
| α | +0.024, 7/10 patients, p = 0.63 | +0.031, 6/10, p = 1.00 |
| β | **−0.067**, 4/10, p = 0.43 | −0.060, 4/10, p = 0.92 |
| θ | +0.076, 8/10, p = 0.084 | +0.060, 8/10, p = 0.16 |

**Three readings, in order of how much weight they carry.**

**(a) The advantage does not replicate.** At n = 10 on the locked substrate there is no α advantage worth the name (+0.024, p = 0.63) and β runs the *other* way (−0.067; `⟨|C|⟩` worse in 6 of 10 patients). W0-B's run used five patients and a different backbone, so this is a **non-replication, not a refutation** — but there was far less advantage here to explain than the premise assumed.

**(b) The same-shaft confound is real but is not demonstrably the explanation.** Masking removes most of `⟨|C|⟩`'s pooled α advantage (+0.090 → +0.028), but the per-patient paired test on that *change* gives p = 0.85. The confound exists in the weights (table above); its effect on the trace cannot be established at this n.

**(c) The decisive objection is band-selectivity, and it is stark.** `⟨|C|⟩`'s largest cross-phase value sits in **θ — the band that must be null**: margin +0.133 all-pairs, +0.120 after masking, against `|ImCoh|`'s −0.017 and −0.010. `⟨|C|⟩`'s θ value is comparable to its own α value. That is disqualifying on the project's own criterion regardless of significance. And because it **survives** same-shaft removal, it is not merely volume conduction: ordinary coherence is sensitive to stable structure that is simply not task-related.

**The reassuring half, stated plainly because it was the most consequential of the three possible outcomes.** `|ImCoh|`'s own trace is **not** same-shaft-dependent: α changes +0.008 on masking (p = 0.85), β −0.035 (p = 0.13), neither significant, and both remain positive with positive margins. The scenario in which both transforms collapse and the trace turns out to rest on within-shaft structure **did not occur**.

**Verdict: the contract stands, and the recorded justification changes.** `imcoh_abs` is retained not because it is larger (it is not, in α) and not because it is VC-immune by construction (§4b's rule 2 forbids that as a justification for a result), but because ordinary coherence fails the project's own band-selectivity criterion on this substrate and `|ImCoh|` does not.

---

## 2. A2 — the sparsification stability surface

### 2.1 What was swept

18 configurations × 6 bands × 4 functionals × 16 scales × 10 patients, matched-strength null at R = 100. The configurations span the density knob under three different *mechanisms* — `mst-union` at 10 fractions from 0.02 to 1.00, the plain global threshold at 3 fractions (same edge budget, no connectivity policy, so it fragments below the percolation point), and the disparity filter at 3 α — plus the two parameter-free filters, TMFG and percolation.

One design decision does real work: within a cell, the matched-strength surrogate draws are **shared across every configuration**. Config-to-config differences are therefore attributable to the backbone and not to surrogate noise. That also makes neighbouring fractions *more* correlated, which is why the plateau claim is not allowed to rest on smoothness alone — it has to survive a change of mechanism at matched density, and it is reported alongside the structural covariates that show the graph itself changed materially across the window.

### 2.2 The gate, and why it is not the one the sweep was run under

The 16-scale axis is worth **n_eff = 1.51** independent tests (mean cross-scale correlation +0.79) on these margin matrices. Treating the scales as 16 independent tests over-charges by ~11×. Sibling lane W0-C reached the same conclusion independently on a 28-point sweep (n_eff 1.2–1.9, mean r +0.70 to +0.93) — two lanes, different data slices, same phenomenon.

So the **primary** gate here is a single sign-flip cluster-mass test (Maris-Oostenveld, whole-patient-profile flips, which preserves the along-axis correlation) per (config, band, functional) over the entire 16-scale margin profile. It collapses the swept axis to one test and does not inherit the per-scale multiplicity problem at all. Everything is gated on the **margin** `obs − surr_p50`, never a raw observed value.

Per-scale cleared-counts are demoted to secondary and reported under two families side by side. Verdict agreement with the primary gate is 0.729 (raw α = 0.05) and 0.750 (whole-grid BH) over 432 cells — **the family changes a quarter of the verdicts**, which is precisely why the plateau is defined on the family-free test. The gate was deliberately *not* swapped mid-sweep, so every row is internally comparable.

### 2.3 The surface

`T_probe` cluster p per fraction, mst-union family; bold = p < 0.05:

| f | density | δ | θ | α | β | low_γ | high_γ |
|---|---|---|---|---|---|---|---|
| 0.02 | 0.028 | **0.009** | 0.311 | 0.064 | **0.025** | **0.020** | **0.034** |
| 0.04 | 0.044 | **0.027** | 0.278 | **0.022** | 0.094 | 0.173 | 0.301 |
| 0.05 | 0.053 | **0.032** | 1.000 | **0.008** | 0.082 | 0.062 | 1.000 |
| 0.07 | 0.072 | **0.011** | 1.000 | **0.009** | **0.022** | **0.013** | 1.000 |
| 0.10 | 0.101 | 0.150 | 0.286 | **0.028** | **0.031** | **0.022** | 0.176 |
| 0.14 | 0.140 | 0.103 | 0.393 | **0.017** | **0.026** | 0.111 | 0.124 |
| 0.20 | 0.200 | **0.013** | 0.426 | **0.018** | **0.004** | 0.098 | 0.083 |
| 0.28 | 0.280 | **0.048** | 0.250 | 0.166 | **0.004** | 0.091 | 0.090 |
| 0.40 | 0.400 | 0.209 | 0.402 | 0.093 | **0.003** | 0.068 | **0.013** |
| 1.00 | 1.000 | **0.015** | 0.322 | 0.120 | **0.016** | 0.061 | 0.125 |

Three-way plateau labels on `T_probe`:

| band | label | widest constant run | positive at |
|---|---|---|---|
| **β** | **invariantly-positive** | f ∈ [0.07, 1.00], 3.84 oct | 8/10 fractions |
| **α** | **invariantly-positive** | f ∈ [0.04, 0.20], 2.32 oct | 6/10 fractions |
| θ | invariantly-null | f ∈ [0.02, 1.00], 5.64 oct | 0/10 |
| low_γ | knob-dependent | — | 3/10, scattered |
| high_γ | knob-dependent | — | 2/10, scattered |
| δ | knob-dependent | — | 7/10, scattered (1111001101) |

**α is recoverable, and this is the highest-value cell in the sweep.** W0-C demoted α on the incumbent substrate — 1 of 28 scales, q = 0.090, verdict broken by 3 of 10 patient drops — and flagged its status as substrate-sensitive. It was right to. Under a gate that does not over-charge the correlated scale axis, α is *invariantly-positive* across 2.32 octaves. The disagreement is the multiplicity family, not the data.

**β is the most knob-robust result in the study**, holding from f = 0.07 all the way to the fully dense graph. **θ is the cleanest null**: not once positive anywhere on a 5.64-octave grid.

### 2.4 The admissible window — the actual substrate decision

A fraction is admissible when both signal bands are positive and the null band is not, on the same graph. An all-zeros plateau on a signal band is a *failure* of that substrate, not a success, which is why the three-way label matters.

- **Required null = {θ}: `f ∈ [0.07, 0.20]`, contiguous, 1.51 octaves, 4/10 fractions.** Passes the pre-registered ≥ 1-octave bar. Over it: α 4/4 (median p = 0.017, margin +0.125), β 4/4 (0.024, +0.127), θ 0/4 (0.409, −0.012).
- **Required null = {θ, low_γ}: `f ∈ [0.14, 0.20]`, 0.51 octaves.** *Fails* the bar. low_γ leaks at f = 0.07 (p = 0.013) and f = 0.10 (p = 0.022) and only goes null at f ≥ 0.14.

The incumbent `f = 0.20` is inside both windows and at the **upper edge** of both. That is what an outcome-selected knob looks like from the outside, and it is the honest answer to the charge in the master plan's fracture 2: the choice was not arbitrary, but it was made at the boundary of the region that supports it, which is the least robust place to stand.

### 2.5 The negative result on the inference-specific functional

`T_probespec` is **invariantly null for α across the entire grid** (0/10 fractions, 5.64 octaves), and likewise for low_γ and high_γ. β is scattered (4/10). Whatever the inference-specific story becomes, it is not a knob-robust cross-phase trace in those bands. `T_probespec_pe` looks much stronger for β (9/10) but is a conditional statistic and is **uncalibrated** — W0-C showed a sham arc built inside pre-task rest can drive such a statistic to p = 0.007, so it cannot be gated until a data-based placebo exists.

### 2.5b What the surface is a surface *of*

Two constraints from lane W0-B change how everything above must be read, and they are recorded here rather than only in the contract because they bear on the numbers in this section.

**Every verdict in §2 is an independence verdict.** Matched-strength agrees cell-for-cell with an independent segment-lattice null (r = 0.983 across 800 cells), so it tests whether the cross-phase structure could arise from independently-drawn graphs with the same strength sequence. It does not test whether that structure is carried by coupling magnitude or by lag — and against a null that preserves magnitude and destroys only lag, the trace shows **no separation at all** (cohort median margin −0.003 for β, +0.004 for α, backbone pinned). The plateau in §2.4 is therefore a stability property of an *independence* result. It says the substrate does not manufacture the finding; it does not say what the finding is made of.

**No functional is zero-centred.** All four, `T_probe` included, are significantly positive at 16/16 scales for α and β on a sham arc carved from a single resting recording with temporal order destroyed. Surrogate-referenced results survive this (the surrogate inherits the same construction), but it means the only admissible statistic is a margin against a per-configuration surrogate. Every gate in this lane is one — the cohort gate is `obs − surr_p50` and the cluster test runs on that margin matrix — so nothing in §2 needed changing, but any downstream reuse must preserve it.

### 2.6 The transform arm (R1-bis resolved)

Run under `imcoh_sq` on 10 configurations. On the 10 shared configurations, `T_probe`:

| band | positive under `abs` | under `sq` |
|---|---|---|
| α | 7/10 | 7/10 |
| β | **6/10** | **3/10** |

On `T_encode`, β goes 4/10 → 1/10. The admissible window under `sq` is a **single fraction, 0.00 octaves**.

So the verdict is not transform-invariant — but not in the direction R1-bis anticipated. α is transform-invariant; **β's trace is specific to `⟨|Im C|⟩`** and does not survive `⟨(Im C)²⟩`. low_γ also leaks more under `sq`, so band-selectivity is worse too.

R1-bis's literal text ("adopt `sq` if the verdict is not transform-invariant") would have me adopt a transform under which no stable substrate exists. **That is an under-specification in my own pre-registration and I am flagging it rather than quietly rewriting it.** The resolution: R1 is a preference ordering over transforms; R2 is a hard admissibility gate, frozen before any number. A (transform, backbone) pair that admits no octave-wide invariant region is not admissible whatever its reliability. `imcoh_abs` is retained because it is the only one of the two that passes R2 — **not because it won R1; it lost.**

The consequence belongs in Methods rather than in a footnote: the β trace is transform-specific, and that is a real fragility of the central claim.

---

## 3. A3 — killing the free parameter

R2.1 preferred a genuinely parameter-free filter over a knob-integrated readout. **No candidate qualified**, and the three reasons are each worth stating because each retires something the project currently believes.

### 3.1 Neither parameter-free filter clears the trace

On `T_probe`, cluster p: TMFG α 0.141, β 0.086; percolation α 0.312, β 0.076. Neither clears in any band, so neither sits inside the plateau and R2.1 fails outright. The free parameter cannot be eliminated — only integrated over.

### 3.2 PMFG is not TMFG, so the standing recommendation is unsound

`.agents/guides/02_methods/sparsification-choice.md` recommends adopting the PMFG/TMFG family, citing PMFG as the principled parent and using TMFG as "its fast chordal computation", and its own adoption rule makes this conditional on clause (iii) **PMFG ≈ TMFG**. That clause had never been tested. Tested here for the first time on the observed per-scale readout, 10 patients × 4 bands × 16 scales:

| functional | Pearson(TMFG, PMFG) | median \|diff\| | sd(TMFG) |
|---|---|---|---|
| `T_probe` | 0.624 | 0.105 | 0.235 |
| `T_encode` | 0.539 | 0.129 | 0.227 |
| `T_probespec` | 0.158 | 0.107 | 0.158 |
| `T_probespec_pe` | 0.455 | 0.109 | 0.173 |

TMFG accounts for under 40 % of PMFG's variance on the standard trace, with a typical disagreement about half the spread of the statistic itself. These are **two different filters**, not one filter and its fast approximation. Clause (iii) fails, so the planar family cannot be adopted on the argument the guide offers. This is a per-scale observed-readout comparison, not a null-gated one — PMFG is O(N³) planarity testing and cannot sit inside a 100-surrogate null — so it establishes disagreement of the readout, which is enough to void the clause.

### 3.3 Planarity, not density, is what kills α — with an explicit mechanism

The earlier "planarity kills the α trace" verdict (2026-07-15, 12/16 → 0/16) confounded the topological prior with the edge budget: TMFG is fixed at 3N−6 ≈ 5 % density while the incumbent was 20 %. Separated here by comparing TMFG against mst-union at **TMFG's own density, per phase** (0.0504 vs 0.0539, ratio 0.947):

| band | TMFG median obs | budget-matched | ratio |
|---|---|---|---|
| α | +0.0416 | **+0.1886** | 4.5× |
| β | +0.0951 | +0.1845 | 1.9× |
| low_γ | +0.0568 | +0.1061 | 1.9× |
| θ | +0.0312 | +0.0398 | 1.3× |

Same number of edges, 4.5× less α trace. The loss is the topological prior and it is band-specific — θ, the null band, is barely touched. The A2 sweep says the same thing at the gate level: mst@0.05 (density 0.053) gives α p = 0.0055 while TMFG (density 0.050) gives p = 0.141.

The mechanism is explicit in the edge composition. At matched size the two graphs share only a third of their edges (Jaccard 0.336, range [0.226, 0.539]). The edges planarity is **forced to discard** number 186 per phase and sit in the **top 3.4 % by weight** (mean normalised rank 0.0336, mean weight 0.111); the edges it takes instead average rank 0.214 and weight 0.061, against a median edge weight of 0.035 overall. Planarity spends its 3N−6 budget on a spatial embedding constraint and cannot fit the strong-edge structure the α trace lives on.

### 3.4 Percolation is parameter-free but not stable

Percolation fixes its threshold at the connectivity bottleneck, so its density is set by the data rather than chosen — genuinely parameter-free. But that density **swings 36-fold across cells**: median 0.178, range [0.018, 0.665], IQR [0.108, 0.257]. A filter whose density wanders that far cannot sit inside a single invariance window, and this is why it fails the trace. "Parameter-free" and "stable" are different properties, and only the second is what a substrate contract needs.

### 3.5 Therefore: route (ii), the knob-integrated readout

The contract ships the R2.3 fallback. `canonical_graph_ensemble` yields the graph at every fraction of `f ∈ {0.07, 0.10, 0.14, 0.20}`, and any reported number is the median over that window with the across-fraction spread attached. No single fraction is ever the pipeline setting. Plateau width — 1.51 octaves for the trace claim, 0.51 for the band-selectivity claim — is reported with it.

---

## 4. A4 — one substrate or two?

**Recommendation: keep two substrates, but replace the justification — and re-verify the marker's before citing it.**

The current pipeline uses mst@0.20 for the cross-phase trace and TMFG for the epilepsy marker, justified as a "magnitude vs skeleton" dissociation. The trace side is now settled and says a single substrate does **not** work:

- TMFG's density (0.050) lies **below** the admissible trace window in density terms ([0.072, 0.200]).
- TMFG fails the trace in every band on `T_probe` (α 0.141, β 0.086).
- The failure is mechanistic, not marginal: at matched density planarity discards the top-3.4 % strongest edges and the α trace collapses 4.5× (§3.3).

So a TMFG-only paper would lose α and weaken β, and an mst-only paper would have to give up whatever the marker gains from TMFG. Two substrates it is — on the trace-side evidence.

But the marker's half of the justification is now **less secure than it looked**, and this is the part L3 must act on. The marker's TMFG win has always been read as evidence for the *planar filter principle* (PMFG as the principled parent, TMFG as its fast computation). §3.2 shows PMFG and TMFG disagree substantially on the cross-phase readout — Pearson 0.16–0.62 — so they cannot be assumed interchangeable on the marker either. **Until L3 verifies PMFG ≈ TMFG on the SOZ readout, the marker's backbone should be described as "TMFG", not as "the planar family", and no principled-parent argument should be made for it.** If they diverge there too, the marker's TMFG result is a property of one specific greedy triangulation, which is a much weaker thing to defend.

A mechanistic hypothesis worth L3 testing rather than assuming: planarity systematically drops the strongest edges (§3.3), and in sEEG the strongest `|ImCoh|` edges are enriched for same-probe contact pairs. A filter that discards them may help a *within-phase* seeded-diffusion read of the SOZ for exactly the reason it hurts a *cross-phase* trace. That would turn the two-backbone posture from an awkward admission into a stated mechanism — but it is a hypothesis here, and the same-probe enrichment of the discarded edges has not been measured.

What L3 inherits concretely: the trace substrate is fixed by this contract; the marker substrate is **not** settled by it; the decisive marker experiments are (a) PMFG vs TMFG on the SOZ readout, and (b) whether the edges planarity discards are same-probe enriched.

---

## 5. The library deliverable

`src/lrg_eegfc/workflow/substrate.py` holds the substrate in one place so the three science lanes cannot silently diverge. `CANONICAL` is a frozen `Substrate` dataclass (transform, backbone, frac, disparity_alpha, plateau, plateau_fracs); `canonical_graph(patient, phase, band)` returns the analysed adjacency; `canonical_eig`, `canonical_phase_graphs`, `canonical_phase_eigs`, `canonical_structure` and `canonical_scale_grid` cover the rest. Per-call overrides exist for robustness sweeps but are explicit arguments, and `dense=True` is the visible escape hatch for a raw-FC baseline rather than a separate code path.

**The phase set is data, never a hardcoded tuple.** `CANONICAL_PHASES` is the five-phase set with `task_learn` first-class, and `canonical_phase_eigs` accepts any phase set and feeds `cross_phase_functionals_over_scales`, which reads the semantic roles (`baseline_a`, `baseline_b`, `encode`, `probe`, `follow`) off whichever phases are present. Dropping `task_learn` returns the standard trace alone through the *same* call, so a four-phase caller and a five-phase caller share one code path.

`canonical_graph_ensemble(patient, phase, band)` yields `(frac, graph)` across the declared invariance plateau — the knob-integrated readout, for use whenever the backbone is not parameter-free, so no single fraction is ever reported as the pipeline setting.

Acceptance-checked on real FC (`scripts/01_compute/paper_final/w0a_verify_canonical_graph.py`), all passing: `canonical_graph` reproduces the incumbent inline pipeline **bit-for-bit** (`max|diff| = 0`); the five-phase path yields all four functionals; the four-phase path yields exactly `{T_probe}`, equal to the five-phase `T_probe` to 0; split-half phases resolve under both spellings; overrides and the plateau ensemble behave.

Supporting library work, all with general names and no manuscript-local tokens:

- `utils/fc/heat_multiscale.py` — `CROSS_PHASE_ROLES`, `cross_phase_functionals`, `cross_phase_functionals_over_scales`. Verified on real FC to reproduce `rho_sym_over_scales` **and** the `05_enc_inf_arc.py` reference implementation of all four functionals to 1e-16. All Spearmans now come from one rank correlation matrix per scale instead of ~10 pairwise `spearmanr` calls, which is what makes carrying the functional axis nearly free.
- `utils/fc/split_half.py` — `imcoh_split_half_adjacencies`, `band_transform_signed`. Generalises the script-local `compute_imcoh_abs_halves` over the transform, so the abs/sq comparison runs through one code path.
- `utils/fc/backbone.py` — `top_fraction_threshold` (the plain global threshold: the mechanism contrast for density sweeps), `backbone_structure` (structural covariates), and the `"thresh"` branch of `select_backbone`.
- `utils/metrics/graph_descriptors.py` — `gini_coefficient`, `weight_heterogeneity`.
- `config/paths.py` — `IMCOH_HALVES_CACHE`, the canonical home for split-half FC under both transforms (the legacy `imcoh_halves_fc` holds the `abs` arm only).

---

## 6. What is NOT settled

Stated first, per the honesty rule, because several of these bear on claims already written up.

1. **The trace does not clear a lag-destroying null.** This supersedes the weaker statement this section originally carried, and it is now the largest open item in the project. W0-B established that matched-strength is empirically an **independence** null — it agrees cell-for-cell with an independent segment-lattice null (r = 0.983 over 800 cells) — so it cannot speak to magnitude-versus-lag. With the backbone edge set pinned to the observed `|ImCoh|` one, the cohort median margin is **−0.003 (β) and +0.004 (α)**, against +0.214 and +0.130 under matched-strength. Every verdict in this report is therefore an **independence** verdict about **coherence structure**, and until N1 is cleared, no claim here may be stated as being about lagged interaction, nor justified by volume-conduction immunity. See `PIPELINE_CONTRACT.md` §4b for the binding wording rules.
2. **The full band-selectivity claim is knob-dependent.** low_γ is null only for f ≥ 0.14, giving 0.51 octaves — below the pre-registered bar. Any claim of the form "the hierarchy rejects low_γ" must carry that window and that width.
3. **The free parameter is integrated over, not eliminated.** A referee can still ask why mst-union rather than another family; the answer is the surface in §2, not a principle.
4. **`T_probespec_pe` is uncalibrated** and cannot support a claim until a data-based placebo exists. This is the specific failure mode W0-C demonstrated (sham arc, p = 0.007).
5. **`T_probespec` is invariantly null for α, low_γ and high_γ** across the whole grid. The inference-specific component is not a knob-robust trace in those bands.
6. **δ and high_γ are knob-dependent** and no verdict is offered for them.
7. **Three transforms were compared and none of the contrasts is significant at n = 10.** `imcoh_abs`, `imcoh_sq`, `⟨|C|⟩`; signed `⟨Im C⟩` and other coherence families were not run. The `⟨|C|⟩` rejection (§1.7) is a directional band-selectivity argument, and it is a non-replication of W0-B's advantage rather than a refutation — the two runs differ in cohort size and backbone and were not reconciled.
8. **PMFG is observed-only** — no surrogate null on that arm, so §3.2 voids the guide's PMFG ≈ TMFG clause but does not itself gate a claim.
9. **n = 10.** The signed-rank p-floor is 1/1024 and the cluster test inherits the same ceiling. Nothing here has headroom for a heavily corrected family — which is the same wall W0-C hit.
10. **The sq arm ran 10 configurations, not 18**, so its window is resolved more coarsely than the abs arm's. The β degradation was verified on the 10 *shared* configurations, so it is not a grid artifact, but the exact sq window edges are less precisely located.
11. **R = 100 surrogates** (reduced from 200 under an earlier machine constraint). The gate uses the surrogate median, which is stable at R = 100; the per-cell upper-tail p resolution is 0.01.
12. **The gate differs from W0-C's locked default.** This lane used the axis-cluster test throughout; W0-C's default is whole-grid BH with the cluster test as the smoothness-respecting alternative. Where the two ledgers disagree per-cell, the difference is the multiplicity family and not the data. Reconciling them is an integration task.

---

## 7. Compute provenance

Environment: `lapbrain` conda env, `PYTHONPATH` pinned to this worktree's `src/`, BLAS threads pinned to 1, worker pools 4–6. Every load-bearing number recomputed fresh from the cached freq-resolved `imcoh` arrays or from raw timeseries; none sourced from a prior report or CSV.

| stage | script | output |
|---|---|---|
| pre-registration (frozen before any number) | `w0a_00_preregistration.md` | — |
| A1 transform head-to-head | `w0a_01_transform_headtohead.py` | `data/paper_final/w0a_substrate/a1_transform/` |
| split-half cache, both transforms | `w0a_01b_cache_halves_both_transforms.py` | `IMCOH_HALVES_CACHE` (240 files) |
| A2 stability sweep | `w0a_02_sparsification_stability.py` | `.../a2_stability/{abs,sq}/` |
| surface analysis | `w0a_03_plateau_analysis.py` | same dirs (`cluster_gate`, `plateaus`, `admissible_window_*`, `margin_surface`, `multiplicity_check`, …) |
| A3 parameter-free filters | `w0a_04_parameter_free_filters.py` | `.../a3_parameter_free/` |
| A1b `⟨\|C\|⟩` vs `\|ImCoh\|` ± same-shaft | `w0a_05_transform_volume_conduction.py` | `.../a1b_volume_conduction/` |
| equivalence checks (real FC) | `w0a_verify_functionals.py`, `w0a_verify_canonical_graph.py` | — |

All scripts under `scripts/01_compute/paper_final/`.

**Validation anchors.** The `abs` split-half arm reproduces the legacy `imcoh_halves_fc` cache to max|diff| ≤ 3.0e-8 for all 10 patients. `T_probe` reproduces the incumbent `rho_sym_over_scales` to 1e-16, and all four functionals reproduce the `05_enc_inf_arc.py` reference implementation to 1e-16. `canonical_graph` reproduces the incumbent inline pipeline bit-for-bit (max|diff| = 0).

**Compute.** A1 ≈ 5 min (10 patients × 3 Welch passes). A2 abs ≈ 42 min (60 cells × 18 configs, R = 100, 4 workers); A2 sq ≈ 20 min (10 configs). A3 ≈ 7 min (40 cells, PMFG dominating at ~10 s/graph). Surface analysis ≈ 3 min at 10 000 permutations. All runtimes were estimated from a timed 1–2 cell probe before launch and printed live with `[i/N] elapsed ETA`.
