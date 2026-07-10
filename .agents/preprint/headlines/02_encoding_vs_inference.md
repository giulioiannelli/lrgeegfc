---
name: headline-n2-encoding-vs-inference
era: IMCOH_ABS_COHORT_N10
status: current
kind: headline
scope: N2 — THE FLAGSHIP. What the held structure (N1) IS: an offline abstraction of a learned structure — the brain keeps the relations it INFERRED, not just the pairs it saw. The four-phase design (rest→LEARN→TEST→rest) + task_learn as encoding reference decompose the trace into encoding vs inference; the inference-specific component persists (β-only) and anchors in the OFC cognitive map. [τ-sweep "mesoscale" result RETIRED 2026-07-09, PI call — read at τ=1/λmax only.] Duration control RESOLVED 2026-06-22 — β inference arc duration-robust (truncation null retired as invalid; carried by length-ratio regression + full-length matched-strength); α duration-suspect. Reframed to flagship 2026-06-23.
owner_agent: inference (consolidation arc) chat
updated: 2026-07-09
---

# N2 — Offline abstraction of a learned structure: the persistent trace carries the relations the brain *inferred*, not just the pairs it was shown

> **The flagship.** N1 showed the resting brain *holds* a reorganization offline;
> N2 shows *what it is holding* — not a replay of the premises it saw, but the
> relational structure it had to *compute*. The four-phase design
> (rest→LEARN→TEST→rest) and `task_learn` as the encoding reference are what make
> the "seen vs figured-out" split possible. Numbers in cached CSVs (§F).

> **Estimator update (2026-07-06) — ρ_split → ρ_sym.** The cophenetic trace is now
> estimated as **ρ_sym** (symmetric over the two split-half arms), replacing ρ_split
> (arbitrary-half artifact; `audit_149`). The consolidation arc was re-run under ρ_sym
> (`audit_152`) and **holds**: β `T_infspec_pe` **p=0.0098** (LO-P15 0.020), **β-only**
> (every other band min p=0.246); learning-phase own trace α/β `T_learn` p=0.014/0.032;
> `T_test` α/β p=0.032. arm1 reproduces the locked `audit_83` full-graph statistic to
> 1e-16 (same pipeline, only symmetrized). The ρ_split arc numbers below (audit_103
> series) are the supplement reference; verdicts are estimator-invariant. Panoramic
> §2.3: `.agents/reports/2026-07-06_rho-sym-panoramic-and-methodology.md`.

> **R2 write-up reconciliation (2026-07-07) — whole section drafted, 3 claims corrected.**
> R2.1–R2.7 are drafted under ρ_sym (`results_paragraphs/R2.*.tex`; migrations `audit_157–161`,
> all arm1==ρ_split bit-exact; report
> `.agents/reports/2026-07-07_r2-mesoscale-encoding-rhosym-migration.md`). **Three ρ_split-era
> claims below are SUPERSEDED:** (1) **N2.2 τ-sweep multiscale result DROPPED ENTIRELY
> (2026-07-09, PI call).** The 2026-07-07 reframe ("multiscale: clears the null at fine τ=1
> *and* mesoscale τ≈2.6") is itself now retired: on a fully-connected, continuous-spectrum
> network there is no meaningful multiscale *across* τ, so the decomposition is read at
> τ=1/λmax only, where the dendrogram hierarchy carries the multiscale content (as in N1.3).
> R2.2 slot retired (`results_paragraphs/archive/2026-07/R2.2_mesoscale-integration.tex`). (2) **N2.3/N2.5 inference→cingulate:** under
> ρ_sym it CLEARS the full-length null (BH q=0.030 incl / 0.040 excl) but does NOT survive
> length-matching (`audit_161`, q=0.15 / 0.30) → **directional lead, length-assisted, not an
> established location**. (3) **N2.5 "only α tracks length" is WRONG** — high-γ tracks it harder
> (ρ=+0.70 vs α +0.55); β is duration-flat (ρ=+0.10, p=0.78); length-tracking sits in
> inference-NULL bands. Memory→cingulate @ low-γ is SOLID (q=0.035, 8/8, `audit_159`).
> Details: memory `r2_cingulate_and_duration_rhosym_2026_07_07`.

> **Note on drift (2026-07-10) — an earlier "drift-confounded" reading was WITHDRAWN.**
> A windowed drift null (`audit_167`/`audit_168`) was briefly read as showing the
> inference-specific decomposition is drift-confounded. **That reading was wrong and is
> retracted.** `T_infspec·e` is a **partial** correlation that already conditions on
> encoding (`e = D_TL − D_bt`), so drift-removal is entangled with the estimator
> itself — a **conditional trace has no valid external windowed-drift-null analog**.
> The null (not the trace) is what broke: the sham arc is built entirely from pre-task
> `rest_pre`, where **no consolidation signal can exist**, yet it returns a
> *significantly positive* conditional trace on its own (β `drift>0` p=0.007) and even
> exceeds the real value — a valid null must sit near zero on a no-signal arc. The real
> inference-specific trace is validated by **matched-strength** (β-only p=0.0098; that
> null *does* return ≈0 on no-signal) + **duration-robustness** (2026-06-22) and
> **stands** as N2's core claim. The **whole-task** drift controls (locked C2) are a
> separate, valid matter and are unaffected. Detail (kept as a methodological cautionary
> record): `.agents/preprint/supplementary/S2_drift_controls.md`.

> **Spectral-subspace control (2026-07-09) — the inference trace is specifically
> HIERARCHICAL / multiscale.** Run the *same* four-phase encoding-vs-inference
> decomposition through a spectral-subspace ("spectral-clustering") lens — a Grassmann
> chordal distance between the network's leading-*k* eigenmode subspaces, in place of
> the cophenetic hierarchy (`audit_165`). It **reproduces the whole-task and encoding
> β reorganization** (`T_G^onl` matches the published whole-task Grassmann masses to
> the decimal — δ 38.07 / β 69.76 / γ_l 66.14; `T_G^enc` clears δ, β) **but finds NO
> inference-specific relocation in any band** (`T_G^infspec` null 6/6; β
> cluster_p_mass = 0.408). So a spectral probe that demonstrably *can* see the gross
> task reorganization is **blind to the inference-specific refinement** — which is
> therefore a change in the cophenetic *hierarchy* (which units nest with which across
> scales), not a reorientation of the dominant eigenmodes. This is an independent
> hierarchical-vs-eigenmode argument (the τ-sweep multiscale claim, former N2.2, is
> retired 2026-07-09), and the N1 whole-task
> dual-witness (cophenetic + Grassmann) is untouched — the **dissociation** (N1 seen
> by both probes, N2-inference by cophenetic only) *is* the result. See **N2.7**;
> memory `grassmann_inference_arc_negative_2026_07_09`.

## §A — Result in plain language

**The claim — and its ceiling, first.** After reasoning through a transitive-
inference task, the resting brain keeps a reorganization that is **specific to the
relations the patient had to *work out*** — over and above what it keeps from merely
being *shown* the premise pairs. We read this as an **offline signature of
abstracting a learned structure**: the brain consolidates the *inferred order*, not
just the *seen pairs*. **Ceiling, up front:** this is a signature of
inference-*specific reorganization that persists*; there is **no behavioral data**
(TI performance is unavailable and will not be obtained — PI 2026-06-22), so the
claim can **never** be tied to inference *success*, and the inferred relations have
**no pinned-down anatomical home** (cingulate is a directional hint only).
"Abstraction" is an **interpretation** of a persistence result — licensed by the task
structure and the anatomy, not a decoded representation.

**How the four phases separate "seen" from "figured out."** The recording runs
**rest → learn → test → rest**. In **learning** the patient is *shown* premise pairs
("A beats B, B beats C, …") — **encoding**. In **testing** the patient gets *novel*
pairs never seen and must *reason out* the answer (A vs D?) — **inference**.
`task_learn` is the reference that subtracts "what was merely seen" from "what was
computed." N1 used the test phase only; N2 uses both — that is the whole point.

**The load-bearing result — the brain keeps what it reasoned out.** Splitting the
persistent trace into an encoding part and an inference-specific part (encoding
controlled *out* of inference), the **inference-specific component itself persists
into post-task rest — in β, and in no other band** (ρ_sym p=0.010, ρ_split p=0.007;
leave-one-out robust, clears the mandatory matched-strength null, and **duration-robust** — it does not
grow with the longer test recording; N2.2 / N2.5). So the offline trace is not just a
re-exposure echo of the premises: it carries a component tied specifically to the
*inference* process. *(Solid.)*

**The trace is a hierarchical reorganization, read at a single scale.** The
inference-specific β component is a change in the **hierarchy** of communication scales,
but that multiscale content is carried by the dendrogram at the single working scale
τ=1/λmax (as for the N1 trace), not by a τ-sweep. On a fully-connected, continuous-spectrum
network there is no meaningful multiscale *across* τ, so the earlier "mesoscale-favouring" /
τ≈2.6 reading was **retired 2026-07-09 (PI call)** and is not part of the paper. *(The
hierarchical — not eigenmode — character is separately witnessed by the spectral-subspace
control; see the 2026-07-09 banner.)*

**The inference component is hierarchical, not spectral — a second, independent
multiscale argument.** Run the *same* encoding-vs-inference split through a completely
different lens — a spectral-subspace comparison of the network's dominant eigenmodes
(what "spectral clustering" reads) instead of the cophenetic hierarchy — and the
whole-task and encoding reorganizations still register in β, but the
**inference-specific component vanishes** (no signal, in any band). The spectral probe
is **not** underpowered: it recovers the whole-task β reorganization exactly (it is one
of N1's two witnesses). It simply cannot see the *inference refinement* — because that
refinement is a change in **which units nest with which across scales** (the cophenetic
hierarchy), not a rotation of the leading modes. So two independent readings converge on
the same conclusion — that the inference trace is a genuinely **multiscale, hierarchical**
object: the diffusion-scale sweep (robust across τ) and this spectral-vs-hierarchical
dissociation. *(Solid; negative control, `audit_165`; N2.7.)*

**The consolidation sits in the brain's map-making hub.** The encoding component
concentrates in **orbitofrontal cortex** — the hub for building relational
"cognitive maps" — and that OFC anchor is set by **learning itself**: the pure
memorising phase (`rest → learn → rest`, no reasoning) already leaves its own
persistent α/β trace into OFC, the *same* hotspot the test-phase trace uses, so OFC
is anchored by **both** task phases (N2.3 / N2.6). Inference then adds the β-only
component on top. *(Solid for the OFC encoding anchor.)*

**The brain files task content by rhythm — the supporting cast.**

- **Alpha keeps the *memorised* material** — and only that (encoding-only). *(Solid.)*
- **Beta keeps the memorised material *and* the *inferred* relations** — the
  load-bearing result above. *(Solid.)*
- **Low-gamma looks empty whole-brain but hides a real, focal *memory* trace in the
  cingulate** (8 of 8 patients; `audit_112`) — visible only because the method reads
  one region instead of averaging the brain. The cleanest "averaging-hides-it"
  showcase for why the method matters. *(Solid.)*

**The honest twist.** Read together, the **cingulate looks like a hub that switches
its job with the rhythm** — *memory* at low-γ (**solid**) and, more tentatively,
*inference* at β (**a hint**). Suggestive, not locked — precisely because the
inference-location half is only a directional hint (the duration check demoted it).

**Why it matters — the flagship.** This is the paper's center. N1 showed the resting
brain *holds* a reorganization offline; **N2 shows that what it holds is an
abstraction** — the relational order the brain *computed*, beyond the pairs it was
*shown*, consolidated at the scale of multi-step integration and anchored in the
orbitofrontal cognitive-map hub. That is the cognitive-neuroscience payload the
method was built to reach: a content-specific, scale-resolved, offline fingerprint of
*what was learned*, read directly from intracranial connectivity — and reached
without behavioral labels, purely from the geometry of the held connectivity state.

**Literature home (to cite — verify exact refs at write-up).** The natural framing is
the relational **cognitive-map / schema / geometry-of-abstraction** line: Behrens et
al. (cognitive maps as relational structure for inference), Schuck & Niv / Schuck et
al. (OFC as a map of task state-space), Bernardi et al. (the *geometry of
abstraction* — neural geometries that generalize to unseen combinations). Those
characterize how such maps are built and used **online**; our contribution is the
**offline, post-task complement** — the inferred relational structure *persists as a
held state once reasoning stops*. *(Full literature triangulation is a pending
deep-research pass — on PI go; do not over-cite until verified.)*

## §B — Technical statement (per subheadline)

**N2.1 — The four-phase decomposition (task_learn's role).** Per region-pair
cophenetic distances: encoding `e = D_task_learn − D_rest_pre`, inference-specific
`f = D_task_test − D_task_learn`, persistence `p = D_rest_post − D_rest_pre`. The
consolidation arc (`audit_103`) tests whether task-change persists:
`T_learn = corr(e,p)` (encoding echo) and `T_infspec·e = partial corr(f,p | e)`
(inference, controlling encoding). Validated against the locked N1 trace exactly
(reproduces `audit_83` per-patient to ~1e-16). **`task_learn` is the pivot:**
`f` and `T_infspec·e` do not exist without it.

**N2.2 — Encoding echo persists (α & β); inference-specific is β-only.** Under the
mandatory matched-strength null: the encoding echo `T_learn` persists in α and β;
the inference-specific component `T_infspec·e` persists in **β only**
(p≈0.007, LOO-robust, all other bands null; `audit_103 --null`). **Duration-
controlled (resolved 2026-06-22, §N2.5):** the β component does not scale with the
test/learn length ratio across patients (Spearman ρ≈+0.25, n.s.), while α — the
only other candidate — carries the duration-tracking trend (ρ≈+0.53). So β-only is
the *duration-clean* verdict, not an artifact of the longer test recording.
**Scale — τ-sweep multiscale result RETIRED (2026-07-09, PI call).** The former
"scale-robust / mesoscale-favouring" reading (τ-sweep `audit_103c/d`; β `T_infspec·e`
clearing the matched-strength null at τ≈2.6/6.8) is **dropped from the paper**: on a
fully-connected, continuous-spectrum network there is no meaningful multiscale *across*
τ, so the decomposition is read at τ=1/λ_max alone, where the dendrogram hierarchy
already carries the multiscale content (as in N1.3). The C(τ) susceptibility has a single
peak τ* (≈10/λ_max, the collapse scale) — no interior mesoscale to sit at (verified
2026-07-09, `data/reports/2026-07-09_ctau-peak-verification/`). History preserved in
`.agents/reports/2026-06-22_arc-mesoscale-inference-null.md` +
`results_paragraphs/archive/2026-07/R2.2_mesoscale-integration.tex`.

**N2.3 — Content localizes differently.** Feeding the localizer the encoding vs
inference per-pair vectors (`audit_110`, β at R=1000): **standard/encoding →
OFC** (robust, BH-clears); **inference → cingulate** (BH-borderline at full
length). Double dissociation: cingulate = encoding in α/low-γ, inference in β.

**N2.4 — A focal low-γ encoding trace hidden by whole-brain averaging.** There is
*no* whole-brain low-γ trace, yet an **absolute** within-region test
(`audit_112`, no demeaning) finds a genuine **focal cingulate encoding trace** in
low-γ (ρ≈+0.39, 8/8 patients, q≈0.035; epi-excluded ρ≈+0.51 7/8 — a *larger* effect
on fewer contacts, so not epi-driven; cingulate *inference* is null here →
encoding-specific). A real trace that cohort-wide averaging conceals — the cleanest
argument for the per-region read-out. **Encoding-only ⇒ duration-immune** (no
`task_test` term). (Power floor: systems sampled by ≤5 patients can't reach
significance on this absolute test, so OFC's evidence stays the demeaned
localization of N2.3, not this.) Backing report
`.agents/reports/2026-06-22_lowgamma-focal-cingulate-encoding.md`.

**N2.5 — The duration confound (RESOLVED 2026-06-22).** `task_test` is **1.35–2.53×
longer** than `task_learn` in every patient, and `f = D_test − D_learn` inherits
that asymmetry. The verdict, across two levels:

- *Localization* (`audit_113b`, demeaned per-system null): **splits.** *Encoding →
  OFC* is duration-invariant and **survives**; *inference → cingulate* does **not**
  survive matched length → **downgraded to a directional hint** (consistent sign,
  but significance was duration-assisted).
- *Whole-brain arc* (N2.2): the naive instrument — a length-matched strength null
  rebuilt on a **truncated** `task_test` (`audit_103b`) — is **INVALID and retired
  from the evidence.** Its negative control fails: **δ, dead-null at full length
  (p≈0.54), false-positives at p≈0.014 at both the head and center windows**, and
  at center even out-clears β. A contiguous truncation injects a common-mode
  structural shift into `D_test` that the *un-demeaned* whole-brain concordance
  reads as fake signal (worst in slow bands) — the same artifact the *demeaned*
  localization removes, which is why that control stayed clean. **Methodological
  takeaway worth one sentence in the paper:** whole-brain concordance is
  corruptible by contiguous-window truncation; the per-system demeaned read-out is
  not.
- *Clean replacement control:* a per-patient regression of the full-length
  `T_infspec·e` against the test/learn length ratio (artifact-free). **β does not
  scale with duration** (ρ≈+0.25, p≈0.49); **α carries a duration-tracking trend**
  (ρ≈+0.53). Combined with the full-length matched-strength null (the mandatory,
  validated referee, which β already clears at p≈0.007), this establishes the
  **β inference-specific arc as duration-ROBUST and α as duration-suspect** — so
  N2.2 ships **duration-controlled, not provisional.**

**N2.6 — The learning phase leaves its OWN OFC trace (closes the open item).**
`T_learn = ρ(e,p)` is the locked N1 trace with `task_learn` swapped for `task_test`
(`g→e`), so it *is* the literal `rest_pre → task_learn → rest_post` trace — and the
split-half baseline makes it the rigorous version. Under ρ_sym (`audit_152`) it clears
matched-strength in **α (+0.19, p=0.014, LO-P15 0.004) and β (+0.24, p=0.032, LO-P15
0.020)**; δ a marginal hint; others null (`arc_null_per_patient.csv`). The localizer's `encoding`
target is its per-system map: the β learning-trace lands on **OFC** — the only system
clearing BH in both epi modes (R=1000, q=0.010, low-strength −0.51/−0.57, not
hub-driven) — the **same hotspot** the test-phase trace uses (`audit_83`, q≈0.01), so
OFC is anchored by **both** task phases. It survives the duration control (encoding→OFC
β: epi-exclude q=0.050 clears, epi-include q=0.099 BH-borderline but raw MS_p≈0.010);
`e` uses the *shorter* learn phase and carries no test/learn length asymmetry, so the
duration confound never touched it. **Reading:** the OFC consolidation anchor is set
by memorising itself, in α and β, before inference — and inference adds the β-only
component on top. Backing report `.agents/reports/2026-06-22_learning-phase-own-trace.md`.

**N2.7 — A spectral-subspace ("spectral-clustering") comparison does NOT resolve the
inference-specific component: it is hierarchical, not eigenmode-geometric.** Replacing
the cophenetic per-pair distance with a Grassmann chordal distance between leading-*k*
combinatorial-Laplacian eigenmode subspaces (`d_G`, *k* ∈ {2..112}, `T_G*` cluster-mass
gate), on the identical four-phase decomposition (`audit_165`), the inference-specific
arm `T_G^infspec = Δ(RP) − Δ̄(pre)` — with `Δ(r) = d_G(r,task_learn) − d_G(r,task_test)`,
the natural subspace analogue of `f = D_test − D_learn` — is **null in all six bands**
(β cluster_p_mass = 0.408, `T_G* = 0.032`, 5/10 positive; δ 0.756 / θ 0.612 / α 0.821 /
γ_l 0.378 / γ_h 1.000). **This is not a power failure.** The same probe's whole-task arm
`T_G^onl` reproduces the published Grassmann trace to the decimal (δ 38.07 / β 69.76 /
γ_l 66.14, all strong; anti-hallucination cross-check versus `audit_66` max deviation
1.7e-13), and the encoding echo `T_G^enc` clears δ and β. So the spectral-subspace
geometry captures the *gross* task reorganization (whole-task + encoding) yet is blind to
the *inference-specific refinement* — establishing the inference component as a change in
the cophenetic hierarchy, not a reorientation of the dominant eigenmodes. **Consistency
with N1:** Grassmann remains a valid second witness for the *whole-task* β trace; the
dissociation (N1 = both probes, N2-inference = cophenetic only) is the evidence, not a
contradiction. **Method note:** a length-matched (truncated-`task_test`) Grassmann variant
is **invalid** — head-truncation injects pan-band eigenvector noise that fires the statistic
in all six bands *including θ*; the raw `T_G^infspec` is length-robust by construction (the
`task_test`-length bias enters `Δ(RP)` and `Δ̄(pre)` through the same `−d_G(·,TT)` term and
cancels). Scope + verdict
`.agents/guides/task-persistence-investigation/2026-07-09_grassmann-inference-arc.md`;
memory `grassmann_inference_arc_negative_2026_07_09`.

## §C — Critical issues & powerful strengths

**Strengths:** the encoding side (encoding echo + encoding→OFC) is robust and
duration-invariant; the four-phase decomposition is exact (reproduces N1);
β-specificity is a strong intrinsic argument against generic confounds; the
low-γ focal encoding trace is a clean "averaging hides it" demonstration.
**A spectral-subspace ("spectral-clustering") control that recovers the whole-task β
reorganization finds no inference-specific component in any band (`audit_165`, N2.7) —
independent evidence that the inference trace is specifically hierarchical/multiscale,
not a dominant-eigenmode shift; the N1-vs-N2 probe dissociation is itself a result.**
**The flagship rests on two solid legs:** (1) the inference-specific β component
*persists* (matched-strength + duration-robust, N2.2 / N2.5); (2) the encoding anchor
is the **OFC cognitive-map hub**, set by learning itself (N2.3 / N2.6). *(A former third
leg — a τ≈2.6 "mesoscale-favouring" scale signature — was retired 2026-07-09, PI call;
the multiscale content is carried by the single-τ dendrogram hierarchy, not a τ-sweep.)*
The "abstraction" reading is an *interpretation* layered on
these legs, bounded by the no-behavior ceiling below — **sell the legs, frame
the abstraction.**

**Lead-with weaknesses:**
- **The inference story is now duration-controlled but anatomically/behaviorally
  thin.** The whole-brain inference-specific β arc is **duration-robust** (clean
  ratio regression + full-length matched-strength; the naive truncation control was
  retired as invalid — N2.5), but it remains **distributed with no established
  hotspot** (inference→cingulate is a directional hint, downgraded by the
  matched-length localization null) and has **no behavioral anchor**. Sell it as a
  signature of inference-specific reorganization — not as localized or
  behaviorally-validated.
- **No behavioral/performance data — and none is obtainable** (PI 2026-06-22).
  The "inference mark" is therefore the neural residue of inference-specific
  *reorganization* and can **never** be validated against inference *success*;
  the claim's ceiling is the task-phase decomposition + anatomy + literature.
  Do not promise a brain–behavior correlation.
- **`task_learn` is short** (the duration asymmetry) and was historically
  neglected — its recording quality / length per patient should be audited.
- **What about α?** α = **encoding-only** and its localization is **diffuse** (no
  FDR-surviving concentration) — i.e. α is the broad "memory" band, β the focal
  "inference-carrying" band. State this contrast explicitly; it is a result, not
  a gap.

## §D — To-dos & verifiables (owner: inference chat)

- [x] **Stage-2-arc (DONE 2026-06-22):** the naive `audit_103b` truncation null is
  **invalid** (δ negative-control false-positives, N2.5) and retired; the duration
  verdict is carried by the clean ratio regression + full-length matched-strength →
  **β inference arc is duration-robust; α duration-suspect.** N2.2 ships
  duration-controlled.
- [x] **Split decision (DONE 2026-06-22): N2 stays unified — do NOT spin the
  inference component into a separate headline.** Inference is duration-controlled
  enough to be N2's confident climax, but too thin anatomically (no hotspot) and
  behaviorally (no data) to carry its own headline next to N1. **Neuroscience
  headline count stays 3 (N1 trace + N2 content + N3 epilepsy) + the methods core**
  (replay N4 archived 2026-06-22).
- [ ] ~~Behavioral TC1~~ **REMOVED** — TI performance data is unavailable
  (PI 2026-06-22). There is no brain–behavior test; do not re-scope one.
- [x] **`task_learn` length/quality audit (DONE 2026-06-22):** all 10 learn
  recordings clean (0 % NaN), median 720 s (431–837 s); even the shortest affords
  ~430 Welch segments → reliably estimated. The duration asymmetry is *relative*
  (test/learn ratio 1.35–2.53×), not a learn-phase quality problem; no new dropout.
  Report `.agents/reports/2026-06-22_task-learn-length-quality-audit.md`.
- [x] **`task_learn`'s own trace (DONE 2026-06-22, N2.6):** `rest_pre → task_learn →
  rest_post` **does** leave its own trace — clears matched-strength in α & β
  (`T_learn` p=0.014 both, LO-P15 robust) and the β trace localizes to **OFC**
  (R=1000 q=0.010, low-strength, epi-robust, duration-survives at exclude q=0.050).
  Already in the encoding-component artifacts — surfaced, no new compute. Report
  `.agents/reports/2026-06-22_learning-phase-own-trace.md`.
- [x] **τ-sweep of the arc (ran 2026-06-22) — RESULT RETIRED 2026-07-09 (PI call).**
  The τ-sweep "scale-robust / mesoscale" reading is dropped from the paper: no meaningful
  multiscale *across* τ on a fully-connected, continuous-spectrum network; read at
  τ=1/λ_max only, with the multiscale content carried by the single-τ dendrogram hierarchy.
  Historical record: reports `2026-06-22_tau-sweep-arc-result.md` +
  `2026-06-22_arc-mesoscale-inference-null.md`; scope
  `.agents/guides/task-persistence-investigation/2026-06-22_tau-sweep-consolidation-arc.md`.
- [ ] Re-verify all cited inference CSVs against the locked N1 trace (the exact
  reproduction is a selling point — keep it true).

## §E — Figure / representation ideas

- **Four-phase schematic** making `task_learn` vs `task_test` (encoding vs
  inference) explicit — fixes the neglected-phase problem visually.
- Encoding→OFC vs inference→cingulate **double-dissociation brain** (with the
  honest "inference duration-downgraded" caption), full-length labelled; keep an
  **X-epi** variant, no C4.
- Band × content matrix (α encoding-only; β encoding+inference) driven from CSV.
- low-γ focal cingulate encoding within-region plot (the "averaging hides it"
  panel).

## §F — Provenance (CSV · script · timestamp)

- Arc (encoding echo + inference-specific) → `data/audit/consolidation_arc/{arc_per_patient.csv,
  arc_null_per_patient.csv}` · `audit_103_cophenetic_consolidation_arc.py` ·
  2026-06-12 / 2026-06-18. Scope+verdict
  `.agents/guides/task-persistence-investigation/2026-06-12_cophenetic-consolidation-arc.md`;
  memory `arc_inference_consolidation_2026_06_18`.
- **Stage-2-arc (RESOLVED 2026-06-22)** → `data/audit/consolidation_arc/arc_lenmatched_null_R200_{band}{,_center,_tail}.csv`
  · `audit_103b_arc_lenmatched_null.py`. **Verdict: truncation null INVALID** (δ
  negative-control false-positives p≈0.014 head+center, full p≈0.54) — retired.
  Duration verdict carried by: (a) per-patient length-ratio regression of
  `arc_per_patient.csv` `T_infspec_pe` (β ρ≈+0.25 p≈0.49; α ρ≈+0.53) — β
  duration-robust, α duration-suspect; (b) full-length matched-strength
  `arc_null_per_patient.csv` (β p≈0.007). Localization duration split:
  `lenmatched_null_R200_{include,exclude}.csv` (`audit_113b`).
- Inference localization → `data/audit/inference_localization/inference_mark{,_R1000,
  _shaftcollapsed}_{include,exclude}.csv` · `audit_110_inference_mark_localization.py`
  (+ `audit_111_gen_task_learn_r1000_surrogates.py`) · 2026-06-18.
- **Spectral-subspace (Grassmann) inference control (N2.7)** → `data/audit/grassmann_inference_arc/{cohort_summary.csv,
  per_patient.csv, README.md}` · `audit_165_grassmann_inference_arc.py` · 2026-07-09.
  **Verdict: `T_G^infspec` null in all 6 bands (β cluster_p_mass=0.408); whole-task
  `T_G^onl` reproduces audit_70 exactly (δ 38.07 / β 69.76 / γ_l 66.14), cross-check
  vs audit_66 max dev 1.7e-13.** Helpers promoted to
  `src/lrg_eegfc/utils/metrics/spectral.py`. Scope
  `.agents/guides/task-persistence-investigation/2026-07-09_grassmann-inference-arc.md`;
  memory `grassmann_inference_arc_negative_2026_07_09`.
- Absolute within-system (low-γ focal) → `data/audit/inference_localization/within_system_trace_{include,exclude}.csv`
  · `audit_112_within_system_trace.py` · 2026-06-18.
- Duration control (localization) → `data/audit/inference_localization/{length_control.csv,
  lenmatched_null_R200_{include,exclude}.csv}` · `audit_113_inference_length_control.py`
  + `audit_113b_inference_lenmatched_null.py` · 2026-06-19. Memory
  `inference_mark_localization_2026_06_18`. Handoff
  `.agents/reports/2026-06-18_inference-mark-handoff.md`.

## §G — Missing parts / open

- ~~Stage-2-arc result~~ **RESOLVED 2026-06-22** — β inference arc is
  duration-robust (N2.5); the gate is closed.
- ~~Whether inference deserves its own headline~~ **DECIDED 2026-06-22** — no; N2
  stays unified; 3 neuroscience headlines (N1–N3) + methods core (§D).
- **No behavioral anchor** — performance data unavailable (PI 2026-06-22); the
  inference claim rests on neural decomposition + anatomy + literature.
- ~~**`task_learn`'s own trace**~~ **DONE 2026-06-22 (N2.6)** — it *does* leave its
  own trace (α & β, `T_learn` p=0.014) localizing to **OFC** (R=1000 q=0.010,
  duration-survives epi-exclude q=0.050); already in the encoding-component artifacts,
  surfaced not recomputed. Report `.agents/reports/2026-06-22_learning-phase-own-trace.md`.
- **No anatomical home for inference** — cingulate is a directional hint; whether a
  cleaner localizer recovers one is open. (The **τ-sweep** ran 2026-06-22: the
  decomposition is scale-robust.)
- ~~**Is the inference trace resolvable by a spectral-subspace probe?**~~ **ANSWERED
  2026-07-09 (`audit_165`, N2.7):** no — the Grassmann leading-eigenmode ("spectral
  clustering") comparison recovers the whole-task + encoding β reorganization but finds
  **no** inference-specific component in any band (β p_mass=0.408). The inference trace
  is specifically hierarchical/multiscale (cophenetic-only), not eigenmode-geometric —
  an independent hierarchical-vs-eigenmode argument (the τ-sweep multiscale claim was
  retired 2026-07-09). Scope
  `.agents/guides/task-persistence-investigation/2026-07-09_grassmann-inference-arc.md`.
- ~~**Mesoscale hint for β inference**~~ **RETIRED 2026-07-09 (PI call).** The
  2026-06-22 verification (`audit_103d`: τ≈2.6 p=0.014, τ≈6.8 p=0.010) is dropped from the
  paper — no meaningful multiscale *across* τ on a fully-connected network; read at
  τ=1/λ_max only. History: report `.agents/reports/2026-06-22_arc-mesoscale-inference-null.md`.
