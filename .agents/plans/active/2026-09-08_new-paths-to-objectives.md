---
name: new-paths-to-objectives
kind: plan
era: IMCOH_ABS × COHORT_N10 → PAPER_FINALIZATION (Wave 2 design)
status: active
created: 2026-09-08
scope: After Waves 0–1 closed every route we had tried, this is the map of the routes that remain to the four objectives (scale-specific on τ; band / spatial specificity; encoding vs inference; SOZ marker), with today's feasibility numbers, the structural reasons the old routes failed, and the lanes to run next.
pointers:
  - .agents/reports/2026-09-03_paper-finalization-checkpoint.md
  - .agents/preprint/locked/PIPELINE_CONTRACT.md
  - .agents/guides/task-persistence-investigation/2026-09-08_hierarchy-level-decomposition.md
  - data/paper_final/feasibility/
---

# New paths to the four objectives (2026-09-08)

**Unverified throughout.** Every number in §2–§5 is a same-day feasibility probe with no null, no knob integration and no leave-one-out. Nothing here is a result. The point of this document is to say which constructions can still reach the objectives, why the old ones could not, and what each new one needs before it may be called a result.

## Head

The old route to "multiscale" was structurally incapable of delivering it, and the two feasibility probes run today say why: on a ~100-node sparse backbone the heat kernel's pair ranking barely changes between s = 1 and s = 200 (effective number of independent scales 1.1–1.5 out of 12), and the one characteristic scale the specific heat shows scatters 3–7× across patients on our working backbone. Whole-vector readouts at different τ were the same statistic wearing different clothes. The construction that escapes this is to take the scale resolution from the diffusion hierarchy itself: bin node pairs by the level of the rest_pre tree at which they merge, so that per-level statistics use disjoint pair sets and are independent by construction. On that construction the β trace is not flat: near zero for pairs joining at the top 2–4 modules, ρ_sym ≈ 0.25–0.30 for pairs joining when 16 or more modules remain, with a within-patient trend in 9/10 patients (p = 0.010) that survives dividing by each level's own split-half reliability and holds across backbone fractions. Encoding vs inference does not separate by level (contrast within ±0.05 at every level), so that objective has exactly one principled construction left, a within-block early/late design that is immune to monotone drift. The SOZ marker remains the least-tested and most promising piece, and the level construction gives it a scale axis too.

## 0. Where the four objectives stand, and whether each failure was structural

| objective | what was tried (Waves 0–1) | outcome | structural or empirical |
|---|---|---|---|
| 1. depends on specific τ scales, across patients | ρ_sym of the whole cophenetic vector at 16–32 scales; contrasts between functionals along s; five scale-local readouts | flat: n_eff 1.06–1.26 vs noise floor 2.9–3.4; q_neff = 1.000 in 36/36 cells; my τ-contrast lead had a 1200× inflated p | **structural** — see §2: the kernel has ~1 independent scale and the tree compresses ~6900 pairs into ~117 values |
| 2. band or spatial specificity | paired β vs low_γ; localization audits (July) | β vs low_γ null both ways at n = 10; β delocalised, α placeless, low_γ focal to SOZ | **empirical at the whole-vector level**; untested at the level-resolved one |
| 3. encoding vs inference | T_learn vs T_test, T_infspec, partial-out-encoding, drift regressor, six bands, ordered sham | no dissociation in any band (q ≥ 0.66); T_infspec reads temporal adjacency (positive on no-task data, β q = 0.006); drift regressor fails structurally (shared D_B) | **structural for block contrasts**: learn always precedes test; today's probe adds that the level profile of learn and test is identical |
| 4. SOZ marker | nothing in this wave; July numbers AUC δ 0.83 / low_γ 0.82 / β 0.75 unverified on the locked substrate | untested | independent of every failure above |

The persistence result itself (task state retained into rest_post at α, β) clears every independence null and is directionally clean against the temporal-drift sham (sham reproduces α 12 %, β 25 % of the margin) but does not reach significance against it (best q = 0.146 with two sham realizations per patient). That is a power limit and is handled by the sham-deepening step in §7.

## 1. The four fundamental issues, and the decision on each

**FC construction.** Keep `imcoh_abs` on the mst-union backbone integrated over f ∈ [0.07, 0.20] (contract). W0-B showed the trace is coherence structure, not lag-specific interaction, so volume-conduction immunity is not its justification; the same-shaft mask (`pair_mask`) is, and every level-resolved statistic below is reported with same-shaft pairs excluded. Bipolar re-referencing is not an option here (it makes the same-shaft bias worse; probe-bias guide §3). One genuinely different FC that is worth one lane as a contrast, not a replacement: orthogonalised amplitude-envelope correlation (Hipp 2012), because it is the coupling mode that carries β resting networks in the literature and a trace present in phase coupling but absent in amplitude coupling (or the reverse) would be a mechanistic statement. It is listed in §7 as optional.

**LRG with τ.** Today's audit (`spectral_range_audit.csv`, 840 cells) fixes what the τ axis can and cannot do on these graphs. Sparse backbones give 2.1–2.4 decades of spectral range but a single dominant specific-heat peak in most cells; the heat kernel's pair ranking has 1.1–1.5 effective independent scales over s ∈ [1, 200]; the peak position is reproducible across patients only for dense, mst0.07 and tmfg (q75/q25 ≈ 1.5–1.8) and scatters 3.4–5.6× for mst0.14 and mst0.20. Backbone density explains 70 % of the variance in these quantities, band 1 %. Consequence: **τ is not a usable independent variable for a cohort claim on this substrate**. The scale coordinate that is usable is the hierarchy level of the diffusion tree, which is comparable across implants when expressed as the number of modules remaining (k) rather than as τ. The LRG still supplies the coordinate system (the tree at s = 1 is the diffusion hierarchy); it stops being asked to supply a τ-resolved readout it cannot give.

**Nulls.** Three rungs stay mandatory for every cohort claim (matched-strength, backbone-fixed N1, order-preserving N3). Two additions: (i) every level-resolved statistic is referenced to the same statistic on the ordered sham, per level, because the finest levels are the strongest pairs and could carry a drift signature of their own; (ii) the ordered-sham ensemble must be deepened beyond two realizations per patient before any margin against drift is called significant. Whole-vector n_eff testing (Lane S) is retired for the level construction because the levels are disjoint by design; the cross-level question becomes a within-patient trend or a coarse-vs-fine paired contrast, both of which cancel the positive bias of ρ_sym.

**Cross-patient aggregation.** Align on k (modules remaining), never on τ or s. Report per-patient profiles and test the profile shape within patients (trend, coarse-vs-fine), then aggregate with Wilcoxon and LOO. For band specificity, replace the single paired test by a repeated-measures design across the six bands (Friedman, then post-hoc), which is what n = 10 can actually power. Pooling pairs across patients into a mixed model is admissible only as a sensitivity analysis, never as the headline, because pair-level dependence within a patient is not modelled by the tools we have.

## 2. Objective 1 — scale-specific on τ, across patients

**Why the old readout could not do it.** Feasibility probe 1 (`feasibility_wavelet_scale.csv`) compared three readouts over s ∈ [0.02, 180]: the cophenetic vector, the heat kernel itself, and the kernel's log-derivative (its graph wavelet, the band-pass object of spectral graph wavelet theory). Below s ≈ 2 all three are identical to the raw FC (split-half reliability constant at 0.72 α / 0.80 β from s = 0.02 to 2.6); above s ≈ 10 the wavelet isolates the lowest Laplacian modes and its reliability collapses (α 0.72 → 0.37 by s = 53). Effective independent scales over the 16-point grid: 2.0–2.4 for every readout. There is about one decade of non-trivial, reliable diffusion content on these graphs, and inside it the readouts are strongly correlated with the raw FC. A whole-vector statistic cannot be scale-specific on this substrate.

**The construction that can.** Build the diffusion tree of the full rest_pre backbone at s = 1 (UPGMA of 1/ρ̂, the existing library tree). For each node pair, record k = the number of clusters remaining at the merge that joins them. Bin pairs by k into log₂ levels: k ∈ {2}, {3–4}, {5–8}, {9–16}, {17–32}, {33–64}, {> 64}. Levels are disjoint pair sets of 170–1100 pairs each, with same-shaft fraction 6–23 % (removed by the mask anyway). Within each level compute ρ_sym on a non-ultrametric readout, either the dense |ImCoh| or the heat kernel at a fixed s. The cophenetic vector itself cannot be used within levels because it is constant within a merge node (this is the ~117-value compression W0-C found). Scale specificity is then a statement about which levels carry the retained change, tested within patients.

**Today's numbers** (`level_decomposition.csv`, cross-shaft pairs, cohort medians):

| β, heat kernel s = 1 | k = 2 | 3–4 | 5–8 | 9–16 | 17–32 | 33–64 | > 64 |
|---|---|---|---|---|---|---|---|
| ρ_sym | 0.02 | −0.07 | 0.21 | 0.17 | 0.25 | 0.28 | 0.30 |
| ρ_sym / split-half reliability | 0.07 | −0.12 | 0.33 | 0.30 | 0.36 | 0.41 | 0.48 |

Within-patient trend (Spearman of level index vs ρ_sym, Wilcoxon one-sided across patients): β heat s = 1 p = 0.010 (9/10 positive), fine-minus-coarse contrast +0.30, p = 0.005; β raw p = 0.024; θ heat s = 4.7 p = 0.005; α flat (p 0.05–0.70); δ and high_γ flat. Robustness (`level_robustness.csv`): β raw p = 0.014–0.024 at every backbone fraction with the reference tree at s = 1, and 0.04–0.12 with the reference tree at s = 4.7; β heat p = 0.003–0.010 at f = 0.14 and 0.20, 0.138 at f = 0.07. The gradient does not change sign anywhere in the sweep.

**What this would mean if it survives.** The retained β change lives in the fine structure of the diffusion hierarchy (pairs that already communicate at short diffusion times) and is absent between the top-level modules. That is a scale-specific, cross-patient, hierarchy-referenced statement, in the units the LRG defines (modules remaining), and it does not depend on the reader believing that τ is comparable across implants.

**What kills it.** (a) The ordered sham shows the same level gradient (then it is a property of drift or of pair strength, not of the task). (b) Matched-strength surrogates rebuilt per level show it (then it is degree structure). (c) The gradient is a reliability artefact not removed by the ratio above; the cleanest check is to reference each level to its own sham, not to its reliability. (d) It disappears under knob integration over f ∈ [0.07, 0.20] and reference-tree scale.

**Lane M1 (scale specificity by hierarchy level).** Scope report first (`task-persistence-investigation/2026-09-08_hierarchy-level-decomposition.md`, written today), then: real arc and ordered sham per level, matched-strength per level, knob integration over f and reference scale, LOO, six bands. Compute: the real-arc part runs in seconds; the sham part reuses `sham_dense` from Lane E and is the same cost as Lane E's sham (hours, not days). Falsifier as above.

## 3. Objective 2 — band or spatial specificity

**Band.** The single paired β vs low_γ test was the wrong instrument at n = 10. The level construction gives each band a profile over seven levels rather than one number, and the profile shape is a within-patient quantity. The question becomes a level × band interaction: does β's gradient differ from α's flat profile? Today's probe says β and θ rise with level, α does not, δ and high_γ carry nothing. If that holds against the sham, the band statement is "β and θ traces are fine-level; α's is level-independent", which is a specificity claim that a repeated-measures test across six bands can carry. Report it with Friedman across bands on the fine-minus-coarse contrast, then post-hoc pairs, then LOO.

**Spatial.** The level construction localises for free: the fine-level pairs that carry the β trace belong to identifiable modules of the rest_pre tree, and each module has a Desikan-Killiany footprint. The July verdict that β has no fine home was reached on whole-vector node scores; a module-level statement ("the trace is carried by pairs inside these k ≈ 16–32 modules, whose anatomy is X") is a different object and has not been tested. Pool across patients by region using the existing `cohort_localization.py` machinery, resolution-preserving, whole-grid BH, LOO. What kills it: modules carrying the trace have no consistent anatomy across patients, which would return us to "delocalised" honestly.

## 4. Objective 3 — encoding vs inference

**What is closed, and why.** Block-level contrasts (learn vs test) at the whole-vector level are reproduced by a no-task arc in the correct temporal order (Lane E), and today's probe shows they do not separate by hierarchy level either: the level profiles of T_learn and T_test are the same, and their contrast is within ±0.05 at every level with sign counts of 2–8/10. Any statistic of the form f(learn) − f(test) on whole-block FC is confounded with temporal order by construction, and no regressor built from phases that share a term with the functional can fix it (the shared-D_B artefact).

**The one construction left.** Learning is front-loaded inside task_learn, whereas task_test is stationary in what it demands. A 2 × 2 within-block design uses time inside each block as its own control: encoding signature = (late-learn − early-learn) − (late-test − early-test). Both blocks have the same within-block time direction, so monotone drift cancels in the difference. Inference then has two readings, and both are testable: (i) the retained change in rest_post resembles the late-learn state (consolidation of what was encoded) versus the test state; (ii) the test block adds a component absent from late-learn, measured as the test-minus-late-learn change projected on the levels that carry the retained trace. Neither needs trial markers, which we do not have and will not get. The null is the same 2 × 2 on the ordered sham (windows carved from one resting recording in true order). What kills it: no within-learn change in excess of the within-test change at any band or level, in which case the honest statement is that the learn and test blocks leave the same offline signature, a task-general trace, and the paper says so.

**Cost.** Windowed |ImCoh| on task blocks (Welch per window) is the heavy step; the existing sliding-window infrastructure (time-window guide, dynamics gallery) covers it. Estimate before launch, as the rule requires.

## 5. Objective 4 — the SOZ marker

**Status.** Never re-derived on the locked substrate. July's node-level AUC δ 0.83 / low_γ 0.82 / β 0.75 with a leave-one-patient-out calibrated detector at median 0.85 in 9/10 patients is the strongest number this project has ever produced, and it touches none of the machinery that failed (no ρ_sym, no cross-phase construction, no drift sham, no τ contrast).

**Why it belongs to objective 1 as well.** July also recorded that the marker is τ-dependent (peaks near s ≈ 10 rather than at the raw-FC limit). On the level construction, the node-level version is direct: a node's SOZ score at level k is built from its pairs inside that level. If SOZ discrimination peaks at a specific k across patients while raw strength does not, that is a scale-specific result for the marker, and one where the cross-patient alignment problem is already solved by k.

**Physiology meets pathology.** July's "β trace spares SOZ" was unverified. On the level construction it becomes a testable, interesting statement: the fine-level pairs that carry the trace exclude SOZ nodes more than chance. If true, one operator reads both the consolidation signature and the pathological tissue, and shows them to be disjoint. That is the sentence a high-tier paper would be built on.

**Lane L3.** Re-derive from scratch on the locked substrate; node-level AUC per band and per level; baselines that a reviewer will demand (node strength, degree, band power, same-shaft removed, distance control); leave-one-patient-out calibration; matched-strength null for the node scores. The benchmark table and the reviewer-demand list are being compiled from the literature and will be appended in §6.

## 6. Literature grounding

(to be appended when the three literature briefs return: LRG scale definitions and GSP scale-resolved statistics; transitive-inference electrophysiology and offline-reinstatement designs immune to time; SOZ-marker benchmarks and reviewer demands)

## 7. Order of work

1. **Lane M1** — scale specificity by hierarchy level, with sham and matched-strength per level, knob-integrated, six bands. Cheap, and it decides whether objective 1 exists.
2. **Deepen the ordered sham** — more realizations per patient (sliding carve positions). Needed by M1's sham rung and by the persistence result's q = 0.146.
3. **Lane L3** — SOZ marker re-derivation with level axis and the reviewer baselines; the trace-spares-SOZ test rides on M1's output.
4. **Lane B1** — band × level interaction and module-level localisation, on M1's output.
5. **Lane E2** — within-block 2 × 2 for encoding vs inference. Heaviest compute; launch after M1 tells us which levels to project on.
6. Optional: amplitude-envelope FC as a mechanistic contrast, only if M1 and L3 stand.

Only after 1–3 is there anything to write. The results and methods sections are not touched until then.

## 8. What this plan does not promise

The level gradient is a same-day probe on one readout family with no null. The ordered sham could reproduce it. The within-block design could return nothing, in which case the encoding/inference claim is dropped and replaced by "task-general". The SOZ marker could fall against strength and power baselines. If all three fail, the honest paper is the persistence trace at α/β against independence nulls plus whatever the marker survives, and it is not a Nature Neuroscience paper. The plan is built so that each of those outcomes is known within days, not weeks.
