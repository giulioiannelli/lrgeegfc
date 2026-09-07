# Lane E — pre-registration addendum: drift control and the paired contrast

**Frozen before any number from the new estimators was computed.** It supersedes nothing in `lane_e_00_preregistration.md`; it adds three estimators that the original did not anticipate, because the original did not know what the sham would show.

## Head

The ordered sham — five pseudo-phases carved in true temporal order out of a single resting recording, containing no task — returns positive cross-phase functionals. That is a task-free explanation for part of what this lane has been measuring, and it lands hardest on exactly the functional the lane was built around. This addendum freezes, in advance, (a) an estimator that conditions the drift away, (b) a contrast that differences it away, and (c) the criterion by which each is judged to have worked or failed.

## What the sham established, and what it did not

Cohort-median margin (observed minus that patient's own matched-strength surrogate median, knob-averaged over the four plateau fractions, median over the sixteen scales), reproduced from disk before writing this document.

| arm | what it destroys | worst cell over the four bands measured |
|---|---|---|
| ordered | nothing (contiguous windows, true temporal order) | `T_infspec` β **+0.061**, θ **+0.044**; `T_infspec_pe` β +0.041 |
| shuffled | temporal order, via 30 s block reassignment | ≤ \|0.027\| everywhere |
| eqdur | temporal order **and** the task-block duration asymmetry | ≤ \|0.044\| everywhere, mostly ≤ 0.02 |

The mechanism is temporal adjacency: phases that sit next to each other in a session resemble each other more than distant ones, and a five-phase arc laid out in time inherits that similarity whether or not a task happened. The matched-strength null cannot see it, because matched-strength shuffles the *finished* connectivity matrix and therefore reproduces the strength distribution while destroying the temporal relationship that produced the effect. This is not a defect of the matched-strength null; it is a null aimed at a different alternative.

What the sham does **not** establish is that everything is contaminated. `T_learn` runs *negative* on the ordered β sham (−0.041) while the real β value is +0.112, and is ≈0 on the ordered α sham (+0.006) against a real +0.149. A no-signal arc does not reproduce `T_learn` at all. Establishing that cleanly is a deliverable of this addendum in its own right, and given how much of the project rests on the persistence headline it is the deliverable I rank first.

## Notation

Five phases: `A`, `B` (contiguous first and second halves of `rest_pre`), `task_learn`, `task_test`, `rest_post`. `D_x` is the condensed cophenetic distance vector of the UPGMA tree of the diffusion communication distance at dimensionless scale `s = τ·λ_max`, on the locked `imcoh_abs` mst-union backbone.

Reorganisation vectors: `e = D_learn − D_A`, `e2 = D_learn − D_B`, `f = D_test − D_learn`, `g = D_test − D_A`, `g2 = D_test − D_B`, `p = D_post − D_B`, `p2 = D_post − D_A`, and newly `d = D_B − D_A`.

`d` is the **task-free drift vector**: same subject, same session, same estimator, a comparable elapsed timespan, and no task inside it. It is already implied by every arc the pipeline has ever computed and costs nothing to extract. `split_half.py`'s docstring fixes the geometry — "a recording is cut in two contiguous halves" — so `d` points from earlier to later in time, which is the direction the ordered sham showed is contaminating.

## Estimator 1 — drift-controlled functionals

`T_x^{|d}` is `T_x` with `d` partialled out of both arguments of every constituent correlation: `T_test_d = ½[pr(g,p|d) + pr(g2,p2|d)]`, `T_learn_d = ½[pr(e,p|d) + pr(e2,p2|d)]`, `T_infspec_d = ½[pr(f,p|d) + pr(f,p2|d)]`, `T_infspec_ped = ½[pr(f,p|e,d) + pr(f,p2|e2,d)]`.

Why the existing symmetric estimator cannot escape drift on its own: the A/B arm average was built to cancel *shared-baseline* bias, and it does. Under a monotone session drift, though, both arms pick up a drift component with the **same sign** — `task_test` is after `A`, `rest_post` is after `B` — so the symmetrisation averages two contaminated arms instead of one clean and one dirty. Conditioning on `d` is the smallest change that reaches the actual alternative.

Two caveats are recorded here so they are not discovered later. `A` and `B` are half-length, so `d` is noisier than the full-phase vectors, and partialling out a noisy regressor **under-corrects** (regression dilution) — a surviving drift-controlled value is therefore an **upper bound** on the drift-free effect, not an unbiased estimate of it. And `d` spans one rest recording while the arc spans much longer, so under non-linear drift `d` has the right direction but the wrong magnitude; the equal-duration sham arm is the comparison that bears on this.

## Estimator 2 — the paired contrast

`C = T_learn − T_test`, per patient, per scale, per plateau fraction, knob-averaged over fractions.

`T_learn = ρ(e,p)` asks whether the encoding reorganisation persists; `T_test = ρ(g,p)` asks whether the *total* task reorganisation persists; and `g = e + f`. So `C > 0` means encoding predicts what persists better than the whole task does, and `C < 0` means the probe-specific increment `f` is carrying persistence. Either sign is a dissociation; a **sign difference between bands** is the strongest form of one.

The two arguments share `e` and share `p` entirely, so `C` differences out a large common component rather than summing two independent errors. Whether it is quieter or noisier than its inputs is empirical and will be measured, not assumed, and Lane S is measuring the same quantity's variance independently.

## Estimator 3 — the paired real-minus-sham difference

For any per-patient statistic `S` (a functional, a drift-controlled functional, or `C`), the reference is that patient's **own** ordered sham:

`ΔS_i(s) = margin_real,i(s) − margin_shamOrdered,i(s)`

where each `margin` is the locked `patient_margin` (observed minus that patient's matched-strength surrogate median), knob-averaged over the four plateau fractions, and the sham side is averaged over the two independent no-signal sources (`rest_pre` and `rest_post`).

This is **not** a test against zero. The comparison object is a surrogate arc that inherits the construction exactly — same patient, same electrodes, same estimator, same phase-duration profile, same backbone, same scale grid, same matched-strength null applied on top — and differs only in that no task happened inside it. It is a stronger reference than matched-strength, not a substitute for it: matched-strength is still applied first, on both sides.

Because the subtraction is within patient, patient-level nuisances (electrode count, graph size, SNR, null height) cancel to first order rather than being averaged over.

## Critical preamble (five points)

**1. The claim.** Encoding and inference leave *separable* traces in post-task rest: there is at least one measurable axis — band, scale, magnitude, spatial loading, or which contact pairs carry it — on which `task_learn`-referenced and `task_test`-referenced persistence differ, and the difference is not reproduced by an arc containing no task.

**2. The null.** For every statistic, that patient's own ordered sham arc, itself referenced to matched-strength. Formally `ΔS = 0`.

**3. The strongest plausible alternative the null must control for.** That the arc's five phases are ordered in time and adjacent phases resemble each other for reasons having nothing to do with cognition — electrode drift, impedance change, arousal decline, sleepiness, medication. Under this alternative every functional is positive, `T_infspec` most of all, because `f = D_test − D_learn` compares the two *temporally adjacent* task blocks.

**4. Does the null control for it, by mechanism?** Yes, and this is the point of the construction: the ordered sham *is* five contiguous windows of a real recording in true temporal order, with window sizes proportional to the real phase durations, so it carries the identical adjacency structure, the identical duration profile and the identical spectral degrees of freedom, and it contains no task. What it **cannot** reject: (i) a task effect that is itself a drift effect of larger magnitude — the sham fixes the *shape* of contamination, not a ceiling on it; (ii) anything about lagged interaction, because `imcoh_abs` and the phase-randomisation null leave `|C_ij(f)|` intact, so every claim here is about **coherence-magnitude structure**, not about lagged coupling; (iii) drift that is non-linear over the hours the real arc spans, since both `d` and the sham are measured inside single rest recordings.

**5. What would falsify the claim, and what remains.** `ΔC` failing its cluster gate in every band after correction falsifies separability *on this axis*. `ΔT_learn` and `ΔT_test` failing would falsify something larger — the persistence headline itself — and would outrank everything else in this lane. Remaining regardless of outcome: the lag null is uncleared for the headline trace as well as for enc/inf, so neither is a claim about lagged interaction; `task_test` is longer than `task_learn` in all 10 patients (ratio 0.40–0.74), and while the equal-duration sham arm addresses the mechanism, the project forbids truncation so no length-matched replication exists; and `n = 10` bounds everything.

## Declared gate and multiplicity family

**Gate.** `axis_cluster_gate` (Maris–Oostenveld sign-flip cluster mass over whole per-patient profiles, 10 000 flips, cluster-forming `z ≥ 1`) applied to the per-patient × per-scale matrix of `Δ`. This collapses the 16-scale axis to **one** test per band, which is the honest charge: W0-A, W0-C and Lane S independently measured that axis as worth ~1.1–1.5 independent tests, so per-scale BH over the full grid over-charges by roughly tenfold.

**Family — primary.** Six bands within one statistic; BH across the six. Chosen because the scientific question is per-band ("in which bands, if any"), and the six bands are the coordinated family that question interrogates.

**Family — secondary, reported alongside.** The full 6 bands × 4 functionals grid (24 cells) with whole-grid BH, as the conservative bound. Any cell that clears the primary but not the secondary is reported as clearing the primary only, explicitly.

**Directionality.** `ΔT_test` and `ΔT_learn` are tested **one-sided** (the claim is directional: real persistence exceeds no-task drift). `ΔC` is tested **two-sided**, implemented as the gate run on `+ΔC` and on `−ΔC` with `p = min(2·min(p₊, p₋), 1)`, because either sign of the contrast is a dissociation and no direction is predicted per band.

**Leave-one-out.** Every reported cell is recomputed under each single-patient drop, including its BH step, and the number of verdict flips is reported. Per-patient counts are descriptive and never the gate.

## Validation criteria, frozen in advance

A drift-controlled estimator is accepted only if it passes **both**. These are tests, not magnitude filters; the magnitudes below them are reported as descriptive companions, never as gates.

- **V1 — calibration.** On the **ordered sham**, `T_x^{|d}` must fail its own axis-cluster gate (p > 0.05, uncorrected, in every band). A statistic that still fires on an arc with no task in it has not been fixed. Descriptive companion reported alongside: the cohort-median sham margin, which for an accepted estimator should be small relative to the uncorrected value it replaces.
- **V2 — retention.** On the **real arc**, `T_x^{|d}` must clear its axis-cluster gate under the declared primary family in at least one band. Descriptive companion: retention, `margin(T_x^{|d}) / margin(T_x)`, reported per band per functional.

Passing V1 only means the estimator has been broken rather than fixed. Passing V2 only means it has not been corrected. **If drift control kills `T_test` and `T_learn` as well, that is the finding, it outranks everything else in this lane, and it is reported first.**

## Analyses this addendum authorises

| id | question | statistic | family |
|---|---|---|---|
| E5 | Is the persistence headline a drift artifact? | `ΔT_test`, `ΔT_learn`, one-sided | 6 bands × 2, BH |
| E6 | Do encoding and inference dissociate? | `ΔC`, two-sided | 6 bands, BH |
| E7 | Does drift control preserve anything? | `T_x^{|d}` under V1 and V2 | 6 bands × 4, BH |
| E8 | Is `T_infspec` recoverable at all? | `ΔT_infspec`, `ΔT_infspec_pe` | reported, not rescued |
| E9 | Is the contrast scale-dependent? | Lane S's null-referenced scale readout applied to `C` | cited, not rebuilt |

## Deviations register

Anything I change after this file is committed is appended here with a reason, including changes that make a result weaker. Every variant tried is reported, including the ones that failed; no construction is searched over until one produces significance.

- **Sham arms.** The re-run uses 1 ordered, 2 shuffled and 1 equal-duration realization per patient per source (previously 1/3/2), to buy the two missing gamma bands within the same wall-clock. The ordered arm — the load-bearing one — is unchanged, and the previous 4-band cells are retained on disk as an independent cross-check of the arms whose realization count dropped.
- **Bands.** Extended from the 4 bands the sham originally covered to all 6, because without low γ and high γ no claim that a dissociation is β-specific or δ-specific is admissible, and the user has explicitly said another band is an acceptable outcome.
