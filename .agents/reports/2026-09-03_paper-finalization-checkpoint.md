---
name: 2026-09-03_paper-finalization-checkpoint
kind: checkpoint
era: IMCOH_ABS × COHORT_N10 → PAPER_FINALIZATION (Wave 0 complete + integrated)
status: current
created: 2026-09-03
scope: State of play for the paper-finalization effort. Wave 0 (substrate / null ladder / cohort gate) is complete and merged onto `integration/wave0`. Records what is settled, what died, what is running, and the infrastructure changes. Written as a compaction-survivable handoff.
pointers:
  - .agents/plans/active/2026-08-25_paper-finalization-master-plan.md
  - .agents/preprint/locked/PIPELINE_CONTRACT.md
  - .agents/reports/2026-08-25_w0a-substrate-contract.md
  - .agents/reports/2026-08-25_w0b-null-ladder.md
  - .agents/reports/2026-08-25_w0c-cohort-gate-and-tau.md
---

# Paper finalization — checkpoint, 2026-09-03

## Head

Wave 0 is complete and integrated. The **trace is real and better controlled than it has ever been** — it clears two independence nulls built on unrelated mechanisms that agree at r = 0.983 across 800 cells, it survives session repartitioning, and it now sits on a substrate reported over a knob-*integrated* window rather than a tuned threshold. Almost every interpretive layer above it failed its controls: the effect is **not lag-specific**, the hierarchy detects **nothing** raw connectivity cannot, and the scale axis looks **flat**. Two lanes are running to decide whether the two things the user requires — a real encoding-versus-inference dissociation, and genuine variability across scales — are recoverable or genuinely absent.

## 1. Settled by Wave 0

| question | verdict | source |
|---|---|---|
| Is the trace real? | **Yes.** Clears matched-strength (β 16/16 scales) and an independent segment-lattice null (β 16/16); the two agree at r = 0.983 over 800 cells. | W0-B |
| Is it tied to the task or to the session? | **Task.** Order-preserving rotation: β 10/16 scales, best p = 0.0049. Not BH-significant (q = 0.186) — suggestive, not decisive. | W0-B |
| Is it about time-lagged interaction? | **No.** Every lag-destroying, magnitude-preserving null fails at every scale; the fixed-backbone control gives cohort median margin **−0.003 (β)**, i.e. no separation. The effect is carried by coupling **strength**. | W0-B |
| Which FC transform? | `imcoh_abs`, kept because it is the **only admissible** option — `imcoh_sq` has a 0.00-octave stability window and β does not survive it (6/10 → 3/10 configs). Reliability alone marginally preferred `sq` (~1.2%). | W0-A |
| Which backbone? | mst-union, reported over **f ∈ [0.07, 0.20]** (1.51 octaves, contiguous). The free parameter could not be eliminated — it is **integrated over**. No parameter-free filter qualifies (TMFG α p = 0.141, percolation α p = 0.312). | W0-A |
| Is α real? | **Yes** — invariantly-positive plateau over 2.32 octaves. This **reverses** W0-C's demotion; the lanes disagreed on multiplicity family, not on data. | W0-A |
| Does the hierarchy see what raw FC cannot? | **No.** 0/168 cells under two different residualizers; not a power artifact (75-90% of variance remains). Closed. | W0-C |
| β scale-invariant vs α scale-tuned? | **Does not reproduce.** Classifier 0.200 vs 0.167 chance (p = 0.31); α not more scale-tuned than β (p = 0.90). The old claim rested on a Friedman *failing to reject*. | W0-C |
| Band-selectivity ("hierarchy rejects low_γ")? | **Knob-dependent** — window shrinks to 0.51 octaves, below the pre-registered bar. The trace claim is robust; this one is not. | W0-A |
| PMFG ≈ TMFG? | **Refuted** (r = 0.16-0.62). Kills the sparsification guide's core clause; guide marked SUPERSEDED. | W0-A |

## 2. Two rules that now bind everything

**No cross-phase functional may ever be tested against zero.** All four — `T_test` included, not only the conditional `T_infspec_pe` — are significantly positive at 16/16 scales for α and β on a sham arc carved from one resting recording *with temporal order destroyed*. White-noise input returns −0.031..+0.030, so it is not the algebra; the construction manufactures positives from shared structure. This does **not** void surrogate-referenced results (the gate has always been observed-vs-surrogate, and the surrogate inherits the construction), but "T > 0 hence a trace" is dead permanently.

**Matched-strength alone is insufficient.** It agrees cell-for-cell with an independent segment-lattice null, which shows it is empirically an *independence* null — it tests whether any genuine coupling is required and nothing more. Every cohort claim now needs three rungs: matched-strength/N1b, **N1 with the backbone fixed**, and **N3 order-preserving**. Where N1 is not cleared (currently α, β, and all four functionals) the claim must be stated as being about **coherence structure**, and volume-conduction immunity must not be its justification.

Also locked: the scale sweep is worth only **1.2-1.9 independent tests** (both W0-A and W0-C measured this independently), so per-scale BH over the full grid overcharges roughly tenfold and is the wrong family. Family choice is the difference between "β holds robustly" and "nothing holds".

## 3. Integration

Branch **`integration/wave0`** merges all three lanes. Only `heat_multiscale.py` collided: W0-A and W0-B had each independently implemented the five-phase cross-phase functionals under different names. **The two implementations were verified to agree to full float precision on identical input before unification** (`T_test` −0.0658031755908399 under both), which cross-validates both lanes' results. W0-A's roles-based version is the single source; it now emits both naming schemes and accepts W0-B's `phases=` tuple, so scripts from either lane run unchanged. All lane entry points import clean: `workflow.substrate` (`canonical_graph`, `canonical_graph_ensemble`, `canonical_phase_eigs`), `utils.surrogate.timeseries_nulls`, `utils.metrics.cohort_gate`, `workflow.phase_graphs`.

**Every reported number must come from `canonical_graph_ensemble` over the plateau fractions** — a single-fraction number is a tuned number.

## 4. Running now

- **Lane S — is the scale axis really flat, or is the readout hiding it?** The motivating contradiction: the tree resolves **118 clusters at the finest scale and 4 at the coarsest**, yet the statistic barely moves. Strongest lead: the cophenetic measure takes only **~20 distinct values** across 6 903 pairs, which is structurally forced (UPGMA on N leaves gives at most N−1 merge heights) and plausibly explains *both* the flatness and W0-C's independent finding that the hierarchy expresses only **1-6%** of the raw reorganization. Closed sub-question: the grid-sampling explanation is **falsified** — W0-C's 28-point sweep already had 14 points below s = 3.
- **Lane E — make encoding vs inference real.** Core task is the calibration problem: `T_infspec_pe` reads positive on no-task data, so no p-value from it is quotable until a formulation sits at zero on a no-signal input. Current positive evidence to re-derive, not assume: encoding persists in β (q = 0.037); inference-specific clears in δ/α/β conditioned on encoding (q = 0.043). Enc/inf is the **one place scale structure has appeared** (2.3-3.2 independent scales vs 1.2-1.9 for the whole-task functional).
- **W0-A residual** — `⟨|C|⟩` vs `imcoh_abs`, **with and without same-probe pairs**. Ordinary coherence gives a *larger* trace (β +0.355 vs +0.281) but carries the 2-8× same-probe bias, so the advantage may be volume conduction. Decides whether the transform choice must be revisited.

## 5. Infrastructure changes made this session

- **Compute cage.** All analysis runs in a cgroup v2 slice (`user@1000.service/compute.slice`) with `memory.high` 12G, `memory.max` 16G, `cpu.max` 10 of 16 cores, plus a `compute-cgroup-herder.service` (user unit, enabled) that captures new `lapbrain` processes every 10 s. Record since installation: **`oom_kill 0`**, 657k soft-throttle events, zero system-wide OOM kills. It cannot kill a process — only throttle.
- **Agent compute must be launched detached** (`setsid ... < /dev/null > log 2>&1 &`) and waited on with a bounded foreground polling loop. Lanes lost compute repeatedly because child processes died with the agent session; this was misdiagnosed as the cage killing them.
- **TeamViewer.** Root cause of the drops was **`TeamViewer_Desktop` crashing ~15 s into each session while the GNOME screen was locked** (`LockedHint=yes`), which broke IPC to the daemon and terminated the session with `IPC_Error`. Screen now unlocked; crashes fell from 4-in-27-minutes to 3 in four days. GUI watchdog hardened (2 s poll, core dumps off). Memory pressure was an aggravator, never the cause — the crash loop predates this work by weeks.

## 6. Open decisions for the user

1. **The narrative.** The multiscale framing is not currently supportable. Lanes S and E decide whether it is recoverable. If both come back negative, the honest paper is: a robust, well-controlled cross-phase trace carried by coupling strength, plus the band dissociation, plus the epilepsy marker — with the null ladder itself as a methods contribution.
2. **L3 (epilepsy) has not been re-tested at all** and is structurally independent of everything that failed — it is a *within-phase* readout that never touches the cross-phase statistic, lag-specificity, or the scale axis. It may now be the strongest result in the paper. Not yet launched.
3. **`T_infspec_pe` cannot be quoted** until Lane E calibrates it.

## Update 2026-09-03 — W0-A and Lane S closed and merged into `integration/wave0`

**W0-A (substrate, closed).** `imcoh_abs` stands, but the recorded justification changed. Ordinary coherence `⟨|C|⟩` does NOT beat it at n=10 (α +0.004 paired, p=0.77; β −0.036, worse in 6/10) — W0-B's 5-patient advantage is a non-replication. What disqualifies `⟨|C|⟩` is band-selectivity: its θ margin is +0.15, larger than its own β (+0.11), and it survives same-shaft masking, so it is not merely volume conduction. `|ImCoh|`'s own trace is NOT same-shaft-dependent (α +0.001 p=0.92, β −0.025 p=0.23) — the worst available outcome did not occur. Contract §4b now binds every lane: matched-strength is an INDEPENDENCE null; three rungs on every cohort claim; where N1 is uncleared (α, β, all four functionals) the claim is about coherence STRUCTURE and VC-immunity may not be its justification. Numbers above are my independent recompute from `data/paper_final/w0a_substrate/a1b_volume_conduction/per_patient_scale.csv`, not the lane's summary; they reproduce the lane's directionally and in magnitude.

**Lane S (scale, closed, NEGATIVE).** The scale axis is flat and the readout is not the reason. Verified from `data/paper_final/lane_s_scale/gate_grid.csv`: δ/α/β clear with the supra-threshold cluster spanning the ENTIRE axis s = 0.02→180 (q = 0.034 each), θ null at 0/32 scales. The 32-point sweep is worth **n_eff 1.06–1.26** independent tests against a **noise floor of 2.9–3.4** — the axis is more redundant than noise. Five scale-local readouts all failed, every one below its own noise floor on the statistic meant to vindicate it. Mechanism measured: four decades of scale perturb the cophenetic vector about as much as splitting one rest recording in half (ρ 0.46–0.59 vs 0.44–0.59). Bonus: the UPGMA tree compresses ~6900 pairs into 117 cells, a hard ceiling of ~11–14% on expressible raw reorganisation — **this IS W0-C's unexplained "hierarchy expresses only 1–6%"**, no orthogonality story needed.

**Two things I found on top of the lane reports, neither yet acted on:**

1. **β is the one live scale-dependence lead.** Its margin rises monotonically +0.052 (fine) → +0.180 (coarse). The lane's null-referenced slope test gives q = 0.078, not significant, and that is the correct test — a naive Spearman over the 32 points gives p = 0.0004 but the axis is worth ~1 independent test, so that p is inflated ~30× and must not be quoted. Directional hypothesis for a pre-registered replication only.

2. **Band-selectivity is weaker than the project has been claiming, and two independent lanes now say so.** Lane S: low_γ margin median **+0.096** against β's **+0.105**, cluster spanning the whole axis, q = 0.063 — i.e. low_γ's effect size is on par with β and it misses the gate only narrowly. W0-A independently: the full band-selectivity window (requiring low_γ null too) is 0.51 octaves, below the pre-registered one-octave bar. Only θ (−0.031) and high_γ (+0.028) are genuinely low. **This is a live threat to paper result #2 (per-band dissociation) and nobody has been assigned to it.**

**Consequence for the paper's three results.** Result 1's "scale-dependent" half is dead — the phenomenon is scale-INVARIANT, which is a different (defensible, narrower) claim and matches `lrg_relational_lens_not_scale_superset`. Result 2 is weakened as above. Result 3 (epilepsy/SOZ) remains untested in this wave and is structurally independent of everything that failed.

**Still running:** Lane E (encoding vs inference).

## Update 2026-09-03 (later) — Lane E sham verdict, computed by me from its cells after it died on a 529

Lane E's agent was killed by a server-side API error mid-run, leaving `data/paper_final/lane_e_encinf/sham/cells/*.npz` complete (10 patients × 4 bands × 2 source recordings). I analysed them directly rather than waiting. Aggregation: margin = obs − median over R=100 matched-strength surrogates, averaged over realizations and the 4 plateau fractions, median over the 16 scales.

**The contamination is temporal ORDER, and only temporal order.** Ordered sham (fake 5-phase arc carved from one rest recording, blocks in natural time order, no task anywhere) returns non-zero margins; the shuffled arm (order destroyed) returns ≤ |0.026| and the equal-duration arm ≤ |0.046|. Matched-strength cannot see this because it shuffles the finished FC.

Ordered-sham margins:

| band | T_test | T_learn | T_infspec | T_infspec_pe |
|---|---|---|---|---|
| alpha | 0.039 | −0.009 | 0.014 | 0.035 |
| beta | 0.020 | **−0.063** | **0.062** | 0.051 |
| delta | 0.025 | 0.029 | 0.010 | 0.009 |
| theta | 0.036 | 0.004 | 0.046 | 0.025 |

Percentage of the REAL margin (`verdict/per_scale_grid.csv`, same aggregation) that the ordered sham alone reproduces:

| band | T_test | T_learn | T_infspec | T_infspec_pe |
|---|---|---|---|---|
| alpha | 30% | −6% | **175%** | 64% |
| beta | **14%** | **−52%** | **89%** | 61% |
| delta | 25% | 36% | 34% | 18% |
| theta | −103% | −10% | 383% | 179% |

**T_test and T_learn are SAFE — this is the most important line in this file.** β T_test keeps +0.123 of +0.143; β T_learn keeps +0.185 because the sham runs *negative* there; α T_learn keeps +0.150. The persistence headline is not a within-session drift artifact.

**T_infspec as constructed is dead.** `f = D_test − D_learn` correlates two adjacent-in-time phases and is measuring block adjacency. T_infspec_pe retains ~40% but was already un-calibrated; 61% sham-reproduced is not a foundation.

### Why the symmetric estimator cannot escape this, and the fix

The swap `½[ρ(D_test−D_A, D_post−D_B) + ρ(D_test−D_B, D_post−D_A)]` cancels *shared-baseline* bias, not drift: under a monotone session drift both arms acquire a positive drift component in the same direction, so the swap averages two contaminated arms. `split_half.py` states the recording is "cut in two contiguous halves", so **A = first half of rest_pre, B = second half**, and therefore `d := D_B − D_A` is a task-free drift vector already computed in every run at zero cost. Proposal handed to Lane E: recompute every functional as a partial correlation controlling for `d`, with a pre-registered two-sided acceptance criterion — **≈ 0 margin on the ordered sham AND retained margin on the real arc**. Caveats to carry: A/B are half-length so `d` is noisy and partialling under-corrects (any residual is an upper bound); and `d` spans one rest recording while the arc spans far longer, so if drift is non-linear `d` has the right direction but the wrong magnitude.

### The dissociation lead, band-agnostic per the user's instruction

Real minus ordered-sham, cohort-median (a LEAD, not a result — not the paired per-patient test, no error bars): α T_learn +0.150 vs T_test +0.091 (encoding ahead +0.059); β +0.185 vs +0.123 (encoding ahead +0.062); δ flips, +0.051 vs +0.075 (test ahead +0.024). Lane E instructed to gate the paired per-patient contrast `C = T_learn − T_test` against each patient's own ordered sham, and to extend the sham to low_γ and high_γ (currently unsham'd, real T_test +0.097).

### The τ question reopened — contrasts, not functionals

Lane S tested scale-locality of one functional at a time and never tested **contrasts between functionals**; two flat quantities can have a non-flat difference. On the real grid `C(s) = margin[T_learn] − margin[T_test]` changes sign along τ: δ −0.053→+0.004 (ρ=+0.64), θ −0.069→+0.028 (ρ=+0.87), β −0.055→+0.013 (ρ=+0.62), α +0.002→+0.038 (ρ=+0.46); and T_infspec alone has δ ρ=−0.83, α ρ=−0.64. **These ρ are hypotheses only** — cohort medians over an axis worth ~1.1 independent tests (p inflated ~15×), and a difference of two noisy quantities is noisier than either, which inflates scale-structure statistics in exactly the direction that trapped four of Lane S's five candidate readouts. Lane S reopened to run the contrast through its own noise-floor and cluster-gate machinery, referenced to the sham's contrast profile, and to state a verdict on the low_γ band-selectivity near-miss it tabulated but never discussed.

## Update 2026-09-07 — Lane S part B closed (merged). Both of the coordinator's leads are dead, and I was the one who was wrong.

Both lanes died on a session limit on 2026-09-03 having finished their compute and most of their commits. Lane S part B is merged; Lane E's verdict script had never been run and I ran it (below).

**The τ contrast lead is an artefact, and my caveat about it was correct in direction but far too generous.** I flagged that the `ρ(log s)` values I handed the lane were inflated "by roughly an order of magnitude" because Part A measured the axis as worth ~1.1 independent tests. Lane S built the actual null — promote a held-out matched-strength realization into the observed slot, rebuild the cohort-median profile, retake the Spearman, 200×. **Under that null the cohort-median |ρ| has a 95th percentile of 0.69–0.82.** So a |ρ| of 0.8 on this axis is a ~1-in-20 event on data with no task in it. Inflation of the naive p: **1200× for `T_learn` β (4.4e-05 → 0.055), 550× for `T_infspec` δ (6.3e-05 → 0.035)**. 0 of 30 (object × band) cells clear at q < 0.05; smallest q anywhere is 0.060. The inflation *grows with the apparent strength of the correlation*, so the most persuasive-looking numbers were the most inflated — which is exactly why the lead looked good to me.

**Contrasts are flat too.** Verified by me from `contrast_scale_structure.csv`: `q_neff = 1.000` for every one of the 36 (object × band) cells — the contrast carries no more independent cross-scale information than its own noise floor produces (n_eff 1.50–2.84 against a floor of 2.97–3.68). No band shows the two clusters a reversal requires; per-patient zero-crossing scales are no tighter than chance (p ≥ 0.44). And `C_learn_minus_test` has `cluster_q = 1.000` in **all six bands** — under matched-strength the contrast does not clear the cohort gate anywhere, so there is no dissociation whose scale-dependence could be asked about. **Answer to "is there anything multiscale only τ can detect": no, not in this family, and the contrast route is now closed too.**

**One nuance not to flatten.** Several cells do show a significant *absolute* slope after BH (`q_slope_abs` ≈ 0.02–0.04 for T_test δ/low_γ/high_γ and T_learn δ/β). Part A already explains it: the margin grows with scale while the surrogate spread grows with it in step, so the effect size in margin units rises (β ~6× along the axis) while standardized detectability does not (per-scale p stays 0.019–0.12 throughout). Report both halves; a reader shown only the margin curve would conclude the opposite of the verdict.

**Band-selectivity: not supported, in either direction.** From `band_separation_verdict.csv`: low_γ margin **+0.095** against β's **+0.092** — low_γ is if anything the larger of the two. low_γ's own clearing is a **gate artefact** (LOO p rises 0.045 → 0.121, `loo_collapses = True`, clears only at f ≤ 0.10, knob-dependent), *but* the direct paired test `β − low_γ` is null in **both** directions (p = 1.000 each way). So at n = 10 we can neither establish low_γ's trace nor demonstrate that β differs from it. **The paper's second result cannot currently be stated as a band dissociation.** This is the third independent line pointing the same way (W0-A's 0.51-octave window, Lane S part A's tabulated near-miss, now the direct paired test).

## Update 2026-09-07 — Lane E closed (merged). I ran its verdict script; it had never been run.

Lane E died on a session limit having committed everything but never executed `lane_e_10_dissociation_verdict.py`. Compute was complete on disk (60 real cells, 180 sham cells, 6 bands, knob-integrated over the four plateau fractions). I ran it (~35 min, 10 000 sign-flip permutations × 132 gate calls) and cross-checked every headline with my own paired Wilcoxons on the per-patient table before trusting it.

### 1. The encoding-vs-inference dissociation does not exist — in any band

`C = T_learn − T_test`, each patient referenced to **their own ordered sham**, two-sided, BH over 6 bands:

| band | Δ C | p | q |
|---|---|---|---|
| delta | −0.011 | 0.220 | 0.661 |
| theta | +0.041 | 0.172 | 0.661 |
| alpha | −0.009 | 0.775 | 0.930 |
| beta | +0.022 | 0.521 | 0.782 |
| low_gamma | +0.029 | 0.455 | 0.782 |
| high_gamma | +0.030 | 1.000 | 1.000 |

The user's instruction was that another band was acceptable. All six were tested. **Nothing.** Lane S independently reached the same place from the scale side: `C_learn_minus_test` has `cluster_q = 1.000` in all six bands under matched-strength. Two lanes, two references, one answer.

### 2. `T_infspec` is demonstrably contaminated and must not be used

On the ordered sham — data with no task in it — `T_infspec` is **significantly positive**: β median +0.053, **p = 0.001, q = 0.006**; `T_infspec_pe` β +0.056, q = 0.046. Descriptive share of the real margin that the sham alone reproduces: **β 84%, α 69%, low_γ 152%, θ 484%**. `f = D_test − D_learn` correlates two phases adjacent in time and is reading adjacency. Retire the construct.

### 3. The persistence result survives in substance but is NOT formally established against this reference

Descriptive, knob-integrated, share of the real margin reproduced by the ordered sham: **α 12%, β 25%, δ 30%, low_γ 27%**, and high_γ 116% / θ negative-real — i.e. the two bands the project treats as null behave as null against this reference too. That is the reassuring half and it is real.

But the paired test does not clear BH. Best cells: α `Δ T_test` p = 0.043 **q = 0.146**; β `Δ T_learn` p = 0.044 q = 0.157; β `Δ T_test` p = 0.073 q = 0.146. Nothing at q < 0.05, and `loo_n_lost` is 7–10 of 10, so it is fragile to dropping a single patient. **This is a power limit, not a refutation** — the sham supplies only two realizations per patient (one per source recording), so the per-patient Δ is noisy. Against matched-strength alone the same functionals are fine (real `T_test` α q = 0.032, β q = 0.032). **Honest status: the persistence headline clears independence nulls and is directionally clean against the drift sham, but is not yet significant against it. Deepening the sham ensemble is the cheapest way to settle it and should be the next compute.**

### 4. My drift-regressor proposal failed, and the reason is structural

Partialling out `d = D_B − D_A` fails **both** frozen criteria. V1 (must be ≈0 on the sham): it **overcorrects** to −0.123 (α `T_test_d`), −0.153 (β `T_learn_d`). V2 (must be retained on the real arc): β `T_test` +0.148 → +0.016 (q = 0.458); β `T_learn` +0.112 → −0.094.

The diagnosis is not the noise caveat I flagged when proposing it. It is a **shared-term artefact**: `p = D_post − D_B` and `d = D_B − D_A` both contain `D_B`, so they are negatively correlated by construction and partialling `d` out removes signal, not just drift. The confirming detail: `T_infspec` — the one functional built from neither A nor B — is the only one the correction leaves alone (β +0.063 → +0.060 on the sham, +0.082 → +0.066 on the real). Anyone revisiting this must build the drift regressor from phases disjoint with the functional it corrects. Recorded so it is not retried in the same form.

### Where the three paper results now stand

| | status |
|---|---|
| 1. Persistence, scale-dependent | persistence clears independence nulls and is directionally clean against drift, **not yet significant against drift**; "scale-dependent" is **dead** (Lane S parts A and B) |
| 2. Per-band dissociation | **unsupported at n = 10** — β vs low_γ null in both directions; and no encoding/inference dissociation in any band |
| 3. Epilepsy / SOZ | **still untested in this wave**, and structurally independent of every failure above |

Result 3 has never been re-derived and is now, by elimination, the most likely surviving contribution. It is a within-phase readout: it never touches ρ_sym, the cross-phase construction, the drift sham, the lag null, or the scale axis. Nothing in Waves 0–1 bears on it either way.
