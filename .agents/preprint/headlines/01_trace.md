---
name: headline-n1-trace
era: IMCOH_ABS_COHORT_N10
status: current
kind: headline
scope: N1 — the multiscale persistence trace left by a transitive-inference task in post-task rest; band-specific by consistent localization (β→OFC; γ_low→PFC DEMOTED to an aside, PI 2026-07-04 — γ_low's cohort trace fails the gate and PFC is the β-depleted cortex) + per-band probe dissociation (α cophenet-only, δ Grassmann-only, β/γ_low both); cohort-level over a variable cohort (detectability control); carried by healthy gray-matter coupling in β; emergent (not raw-strength).
owner_agent: localization + white-matter + preprint-general-questioning(null)
updated: 2026-07-07
---

# N1 — A transitive-inference task leaves a held, multiscale connectivity trace in post-task rest that consolidates to a consistent cortical home in β → orbitofrontal cortex

> Emerges from CORE. Numbers live in cached CSVs (§F). Raw FC is the baseline,
> never the result.

> **Estimator update (2026-07-06) — ρ_split → ρ_sym.** The cophenetic-trace
> estimator of record is now **ρ_sym**, the symmetric average of the two
> equally-valid split-half arm assignments. Bare ρ_split hard-coded *which*
> rest_pre half baselined *which* contrast — an arbitrary choice that flipped the
> sign of near-zero patients (20/60 patient×band cells flip under A↔B swap;
> `audit_149`). ρ_sym is invariant to that choice at zero data cost. **Every
> verdict below is estimator-invariant (0/6 bands flip).** The ρ_sym numbers are
> deliberately *more conservative* (β gate p=0.032, α p=0.024 at R=200) because the
> estimator pulls ill-conditioned near-zero patients toward zero instead of letting
> a lucky half push them positive — that is the estimator doing its job, not a
> weakening. Per-patient reporting is now **ρ_sym ± ½|ρ_AB−ρ_BA|**, with
> `|ρ_sym| < 1 SE` labelled **"undetermined."** ρ_split is retained only as a
> supplement. Full read: `.agents/reports/2026-07-06_rho-sym-panoramic-and-methodology.md`.

## §A — Result in plain language

Patients learned an order from pairs ("A beats B, B beats C, …") and were then
tested on pairs they had never seen, which they had to reason out (transitive
inference). We compared the brain's resting network organization **before** and
**after** the task. The finding: **the task reorganizes the network and the
reorganization does not fully wash out — it persists into the post-task rest.**

Three things make this a real result rather than a measurement artifact:

1. **It is band-specific — in where it lands, not in whether it exists.** The trace
   is a cohort-level effect, clearest in **β** (cohort-median ρ_sym +0.20; 6/10 patients
   clear their own strength-matched null; gate p=0.032). The median only *ranks* the
   rhythms; the **control certifies** them — γ_low carries the largest single-patient
   traces yet fails the cohort null, so peak size is not the verdict. The **size of +0.20**
   is set by the null, not by 1: ρ_sym ignores coupling strength, and strength-matched
   surrogates reproduce almost none of it (median ρ ≈ 0.01), so +0.20 sits an order of
   magnitude above the strength-only baseline — and it *understates* the effect, since half
   the cohort at the measurement floor drags the median down while the patients that do trace
   reach ρ = 0.3–0.5. The between-patient spread is
   real, not a measurement-reliability artifact (detectability control, Methods), so
   we say *cohort-level over a variable cohort* — never "every patient shows it."
   Several rhythms carry a trace, but it consolidates to a **consistent anatomical
   home** most clearly in **β→OFC**; in α it is real but spatially scattered, **γ_low**
   shows a PFC concentration on a band whose cohort trace fails the gate (a **demoted
   aside**, PI 2026-07-04), and **θ shows nothing**. An **independent global-mode readout confirms
   β**; the two readouts are otherwise complementary and the full readout × band map
   is a **Methods** detail — the biology (OFC consolidation), not the probe
   comparison, is the result.

2. **It is an emergent, multiscale property, not a strength reshuffle.** At the
   level of raw connectivity, β is actually the *weakest* band; the clean
   persistence only appears once we read the multiscale hierarchy, and it
   **survives the strength-matched surrogate** — so it is not just "the strong
   edges moved around." Plain edge-by-edge comparison would miss it — **and so
   would an off-the-shelf spectral-clustering / PCA pipeline.** The α companion of
   the trace is the cleanest proof: it is recovered by the full multiscale read-out
   *only*, invisible to the leading-mode spectral view however you read it (N1.2b).

3. **It has an address and a tissue type.** The persistent β reorganization is
   brain-wide but **over-expresses in orbitofrontal cortex (OFC)** — a hub for
   building relational "cognitive maps," exactly what a transitive-inference task
   should leave behind — and it is carried by **healthy gray-matter cortical
   coupling**, *not* by the epileptic seizure core. So it is neither everywhere-
   and-meaningless nor a pathology artifact. **And the tissue rule is
   band-dependent** — β *spares* the diseased core while α *recruits* it (N1.7) —
   the bridge from this cognitive trace to the epileptogenic read-out (N3).

4. **It is a state the brain *holds*, not a fading echo.** Watching rest unfold in
   time (20 s windows), the post-task rest does not merely drift near the task
   pattern on average — it **dwells** there, window after window (**10/10 patients,
   p=0.001**), at *every* scale of the hierarchy. We tested hard for the opposite —
   that the trace flashes back in brief, replay-like *bursts* — across every
   representation, scale, band, and timescale down to the ripple range; there are
   **no bursts**. So the trace is a **sustained reinstatement** the brain settles
   into and holds offline, not a transient replay process (N1.6).

**Why it matters:** an offline, multiscale, cortical trace of a reasoning task,
sitting in the β band and the orbitofrontal map system, that the resting brain
**actively holds as a state** — the resting network does **not** return to where it
started (a non-ergodic, held reorganization), read directly from intracranial
connectivity. **This sets up the flagship (N2):** N1 establishes that the brain
*holds a structure offline* and that our method can *see* it; **N2 asks what that
held structure *is* — and finds it is an abstraction: the relations the brain
reasoned out, not just the pairs it saw.**

## §B — Technical statement (per subheadline; reference CSVs, no tables)

**N1.1 — The trace exists (cophenetic, the primary probe), cohort-level.** On the
per-pair multiscale cophenetic probe `ρ^coph` (estimated as **ρ_sym**), post-task
rest sits closer to the task state than pre-task rest does, at the cohort level.
**There is no one-axis descriptor of this result: the median ρ_sym only *orders* the
bands; the binding verdict is the control — per-patient consistency against each
patient's own strength-matched null (C3).** By that gate the trace holds for
**β** (median +0.20, 6/10, gate p=0.032; passes C1–C4) and **α** (median +0.10,
p=0.024); **γ_low** carries the *largest* individual traces and localizes to PFC (N1.3)
yet **fails** the cohort gate (ρ_sym p=0.080, as under ρ_split) — the living proof that
peak/median is not the verdict; **θ** shows nothing. **Reading the magnitude of +0.20:**
ρ_sym ignores coupling strength, so its scale is set by the null, not by 1 — strength-matched
surrogates reproduce almost none of it (median ρ ≈ 0.01), placing +0.20 an order of magnitude
above the strength-only baseline, and the cohort median is itself a floor (half the cohort at
the measurement floor; the patients that trace reach ρ = 0.3–0.5, §C). Anchor the magnitude on
the null **value** (≈0.01), **not** the ×null **ratio** — the ratio does not rank bands
(γ_low's 48× is larger and fails the gate). The between-patient spread is a **detectability axis** but
**not** purely measurement reliability — a residual who-traces axis survives full
per-phase reliability correction (Methods) — so the trace is **cohort-level over a
variable cohort**, not universal. **Confirmation:** an independent **global-mode
Grassmann** readout of the same hierarchy **also clears at β** (cluster-extent
permutation, LOO-robust) — the two-probe convergence at the flagship band. The
probes are otherwise complementary (α cophenetic-led, δ Grassmann-led); the full
readout × band map is a **Methods** table, not a headline (and the Grassmann arm is
not yet baseline-matched to the cophenetic split-baseline — §D). *(C3 matched-strength
is the binding control; see CORE §B.)*

**N1.2 — Emergent multiscale property, not raw strength.** The three-layer
contrast (raw FC → raw `D(τ_max)` → cophenet) shows β is the *weakest* band at
raw FC yet clears matched-strength at the cophenet layer — the multiscale wrap
resolves it. The matched-strength pass rules out a degree-reorganization account.

**N1.2b — Invisible to spectral clustering / PCA, not just edge-wise (the α proof).**
A head-to-head (`audit_143`) puts a literal **leading-k spectral embedding** (the
representation k-means clusters) and the subspace Grassmann through the **same
`ρ_split` statistic** as cophenetic, swapping only the per-pair distance (the
cophenetic arm reproduces the locked α/β `ρ_split` bit-exact). **α is recovered by
the full multiscale cophenetic distance alone** — *both* spectral read-outs (per-pair
embedding *and* subspace Grassmann) miss it — so a standard spectral-clustering / PCA
analysis would not see the α trace at all. Cophenetic is also uniquely **selective**:
it fires on exactly the two signal bands (α, β) while the textbook embedding
over-detects into the null band (θ). *(Honest limit: **β is a tie** — the leading-k
embedding recovers β as well, so the spectral-superiority claim rests on α +
selectivity, never on β.)* Numbers: `data/audit/spectral_distance_swap/cohort_summary.csv`.

**N1.3 — Consistent localization: β→OFC (γ_low→PFC demoted to an aside, PI 2026-07-04).** The trace
consolidates to a reproducible anatomical home most clearly in β. **β → OFC:**
the brain-wide β cophenet trace concentrates, **above each patient's own demeaned
baseline**, in the OFC system (matched-strength R=1000, rank concordance, canonical
observed trace; BH-clears across the 9 a-priori systems within β, in all four
contact/shaft × epi-incl/excl conditions; **leave-one-out immovable** — OFC
top-ranked in 10/10 folds at q≤0.05; shaft-collapse robust; low-strength/non-hub;
bilateral; peaks at the **system** scale, washes out at lobe/hemisphere), with
**sensorimotor and PFC depleted** (lower-tail matched-strength) as the anatomical
complement. It is a **hotspot on a distributed trace, not a container** (within-OFC
lean ≈ 0.59). **Estimator-invariant:** the ρ_sym migration reproduces this
like-for-like (`audit_151`/`audit_155`) — OFC carrier at R=1000 q=0.010–0.015
(matching the canonical ρ_split q=0.009–0.013), sensorimotor + PFC depleted
q_lo ≤ 0.033, LOO-robust in 3/4 conditions; the one marginal case (epi-incl × contact,
worst-drop q=0.060) is the long-documented n=5 OFC-coverage fragility, not an
estimator effect. **γ_low → PFC (demoted aside, PI 2026-07-04):** within γ_low a PFC concentration does clear
the localization null (all four conditions; carrier p<0.05 in 10/10, top-ranked 8/10, occipital
depleted) — **but γ_low's cohort trace itself fails the gate (N1.1, ρ_sym p=0.080)**, and the
concentration sits in the cortex where **β is depleted**, so it does **not** corroborate the OFC
cognitive-map story. **Prefrontal cortex is a reproducible consolidation home in no band.** Do
not present γ_low→PFC as a second localization. In **α** the cohort trace
is real (ρ_sym p=0.024) but has **no surviving carrier** — spatially scattered,
patient-specific. θ and γ_high localize (MTL, parietal) on a **net-null** band —
reproducible spatial structure, **not a trace**; do not report as one. *(Per-band
consistency taxonomy: consistent {β} / γ_low PFC demoted to an aside (trace fails gate) /
patient-specific {α, δ, γ_high} / absent {θ}. The OFC verdict flipped twice before locking — state that.)*

**N1.4 — The β trace is carried by healthy gray-matter coupling, not the epi core.**
Stated in the LRG framework, **for β**: the genuine pair-count hotspot is **gray↔gray
cortical coupling** (WM↔WM depleted); pairs internal to the clinically-labelled
seizure core do **not** carry the β trace (fail their matched-strength null; "core
spared" is directional only). *(Band-dependent — α instead recruits the core; see N1.7.)* Excluding SOZ contacts leaves β/α intact. **Two prior
sub-claims withdrawn (brutal honesty):** "epi-exclusion *strengthens* the trace"
(generic node-count, decimation control) and "the epi↔healthy interface is the
strongest carrier" (not a pair-count hotspot). The WM result is parallel: the
trace **survives** a gray-only montage but the apparent "sharpening" was
node-count, not WM-specific.

**N1.5 — τ-robust, fine-scale.** A τ-sweep (`audit_121`) confirms the trace is
fine-scale and stable around τ = 1/λ_max; coarse-τ "gains" are a collapse
artifact. The result is not a scale cherry-pick.

**N1.6 — The trace is a SUSTAINED STATE the brain holds (dynamical reinstatement),
not a transient replay.** Resolving rest into 20 s windows, the per-window task-
likeness contrast `g = ρ^coph(window, task) − ρ^coph(window, rest_pre)` is higher in
post-task rest than in pre-task rest across the cohort (**10/10, p=0.001, LOO 0.002**;
`audit_125`/`audit_129`) — i.e. N1 expressed *dynamically*: rest_post **dwells** in
the task configuration. The hold is **scale-invariant and representation-invariant**:
the per-window level is positive at every diffusion scale τ and in **five distinct
Laplacian-propagator read-outs** — cophenetic, magnetic/directional, raw propagator,
subspace/Grassmann, normalized (`audit_133`–`139`, shift +8–10/10, p≤0.007 throughout).
The complementary **transient/replay** reading is a *rigorous cohort negative*: no
burstiness, no separable states, no isolated flashes, no sequence — across every
collapsed axis (scale, target, per-pair/OFC), every band, and the whole resolvable
timescale from 20 s down to the **high-γ coherence floor at ~0.2 s** (the ripple
regime; `audit_124`–`141`), each gated by a time-shuffle + node-permuted-placebo null.
So the reinstatement is a *held state*, not a flashing process. **Inherits N1's
matched-strength** (it is the same `ρ^coph` effect at finer temporal resolution, not a
new FC quantity); a dedicated window-level matched-strength is the one optional control
to make N1.6 fully standalone (§D). The directional-flow variant was tested and
**retracted** (magnitude-inherited, `audit_140`).

**N1.7 — BRIDGE to the clinic: β spares, α recruits the epileptic core.** The
persistence mechanism treats diseased tissue *oppositely across bands*, and this
is the natural bridge from the cognitive trace (N1/N2) to the epileptogenic
read-out (N3). In **β**, the trace is carried by healthy gray↔gray cortical
coupling and the seizure core is **spared** — pairs internal to the clinically
labelled SOZ do not carry it (N1.4). In **α**, by contrast, the same per-pair
persistence **recruits** the epileptic core — epi↔epi coupling is concentrated and
clears matched-strength (C5). So the focal, inference-carrying band (β) avoids
diseased tissue while the diffuse memory band (α) pulls it in: the cognitive trace
and the epileptogenic network are **not independent**. Present this as the bridge,
not an aside. *(Verification brief: `verification/verify_beta_spares_alpha_recruits.md`.)*

## §C — Critical issues & powerful strengths

**Airtight:** band-specificity; matched-strength pass on both probes (β);
emergent-not-strength; OFC localization survives the full six-rung gate incl.
shaft-collapse + R=1000 + LOO; volume-conduction-immune substrate; τ-robust.

**Lead-with weaknesses / provisional:**
- **OFC coverage is n=5** (irreducible) and the verdict flipped 3× historically —
  present as *concentration above baseline*, never "the trace lives in OFC."
- **Per-patient spread is now cleanly tiered (ρ_sym, `audit_154`).** β splits into
  **6 stable tracers** (Pat_08 +0.54, Pat_02 +0.51, Pat_05 +0.50, Pat_03 +0.37,
  Pat_06 +0.27, Pat_07 +0.13 weak; SE ≤ 0.10, own-surrogate p<0.05), **2 undetermined**
  (Pat_13 +0.07, Pat_15 −0.08 — `|ρ_sym| < 1 SE`, sign is estimator noise; Pat_13 is
  108/119 neutral nodes), and **2 stable resets** (Pat_10 −0.14, 44/113 anti nodes;
  Pat_14 −0.03). ρ_sym does **not** manufacture a trace where there is none — it
  separates genuine resets (Pat_10/14) from genuinely-undetermined near-zero patients
  that ρ_split's arbitrary half had been calling tracer-or-anti by coin-flip. **Note:**
  Pat_15 is no longer a clean "β-anti sanctioned dropout" — under ρ_sym it is
  **undetermined** (−0.08 ± 0.17), not anti; revisit that dropout framing. Pat_02 is
  still the β LOO argmax driver (partly epi-coupled). Report LOO-max; no single-patient
  "strong" tag.
- **`task_learn` was not used** for the core trace — "the task" = `task_test`
  only. This reference-phase choice must be justified (see §G); it is the bridge
  to N2.
- **θ is not simply "anti-trace"** — characterize, don't assert (README §5).
- **N1.6 sustained-state framing**: the *level* (rest_post dwells in the task config)
  is essentially N1 dynamically — its strength is robustness (5 representations, all
  scales) + the clean transient-negative boundary, **not** a new independent effect.
  The genuinely-beyond-static-N1 piece ("rest_post is a *tighter* attractor than
  rest_pre") is only a weak trend (A2b p=0.08) — do **not** claim "tighter," claim
  "held." Window-level matched-strength not yet run (inherits the static C3); see §D.

## §D — To-dos & verifiables

- [ ] **(cross-probe rigor, PARKED)** Grassmann uses the single-rest_pre triangle
  (`audit_66`/`audit_70`), not the split-baseline `ρ_split` (`audit_63`) — the
  cross-check is not fully baseline-matched. A split-baseline Grassmann would firm up
  (or retire) the δ Grassmann-only cell; not on R1's critical path (Methods caveat).
- [ ] (localization) Re-confirm OFC CSVs still match the locked verdict; produce
  the system-scale concentration figure with X-epi panels (not C4).
- [ ] (white-matter) Re-state the gray-matter-dominant pair-count result as the
  *significant* finding; keep the decimation-controlled withdrawals explicit.
- [ ] (null model) θ characterization; outlier-leverage panel (LOO-max per band).
- [ ] (null model) τ-robustness figure/supplement from `audit_121`.
- [ ] **(`task_learn` gap)** does `rest_pre → task_learn → rest_post` leave its
  own β trace, and does it localize to OFC too? Decide whether N1's "task" should
  be test-only, learn-only, or both — coordinate with N2.
- [ ] **(N1.6 optional bulletproof)** window-level matched-strength surrogate on the
  per-window `g` (make the sustained-state claim standalone, not inheriting C3).
  Cheap; only needed if a referee treats N1.6 as independent of static N1.

## §E — Figure / representation ideas

- Cohort β trace forest/gate plot (per-band verdict), **X-epi variant** alongside.
- OFC system-scale concentration brain (bilateral; system not region), with the
  honest "hotspot not container" annotation in caption (not in-figure).
- Gray↔gray vs WM↔WM vs epi-core pair-class contrast (drives the tissue claim).
- Three-layer band-resolution strip (shared with CORE).
- **No C4 cross-probe figure.**

## §F — Provenance (CSV · script · timestamp · locked-in)

- ρ^coph C1/C2/C4 → `data/audit/ctm_triangle/cohort_summary.csv` ·
  `audit_33_ctm_triangle.py` · 2026-05-26 (+ `c4_wilcoxon_cohort.csv` ·
  `audit_71_c4_wilcoxon_cohort.py` · 2026-05-19). Locked: `../locked/VERDICT_LEDGER.md`.
- ρ^coph C3 matched-strength (ρ_split, retained as supplement) →
  `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` ·
  `audit_63_split_baseline_surrogate.py` · 2026-05-15.
- **ρ_sym gate (estimator of record)** → `data/audit/rho_sym_gate/` ·
  `audit_150_rho_sym_gate.py` · 2026-07-06 (β p=0.032, α p=0.024; rest fail).
  Estimator-robustness (split=sym, 0/6 flips; full-baseline rejected) →
  `data/audit/estimator_regate/` · `audit_149_estimator_robustness_regate.py`.
- **β→OFC under ρ_sym** → `data/audit/localization_atlas_rhosym/` ·
  `audit_151_localization_rhosym.py` (R=200) + `audit_155_beta_ofc_robustness_rhosym.py`
  (R=1000/LOO/shaft; OFC q=0.010–0.015) · 2026-07-06.
- **ρ_sym downstream** → cross-phase taxonomy `data/audit/cross_phase_taxonomy_rhosym/`
  (`audit_153`; β trace-dominant, comp_trace 0.263), per-node carrier/anti
  `data/audit/per_node_trace_decomposition_rhosym/` (`audit_154`; node ρ_sym↔ρ_split
  0.96), consolidation arc `data/audit/consolidation_arc_rhosym/` (`audit_152`) ·
  2026-07-06. Panoramic: `.agents/reports/2026-07-06_rho-sym-panoramic-and-methodology.md`.
- Grassmann → `data/audit/grassmann_cluster_extent/cohort_summary.csv` ·
  `audit_70_grassmann_cluster_extent.py` · 2026-05-26.
- OFC localization → `data/audit/localization_atlas/matched_strength_R1000_{include,exclude}.csv`
  · `audit_83_localization_matched_strength.py` + `audit_92_localization_R1000_surrogates.py`
  · 2026-06-09; README `data/audit/localization_atlas/README.md`. Locked:
  `../locked/ANATOMY_LEDGER.md` (2026-06-10 reinstatement).
- Gray-matter / epi pair-class → `data/audit/epi_stratified/{cophenetic_cohort.csv,
  pairclass_decimation_cohort.csv}` · `audit_77_epi_stratified_cophenetic.py`
  (2026-06-16, n=10) + `audit_85_wm_decimation_control.py` (decimation control,
  2026-06-12). Cascaded in `../bands/01_beta.md` §3.2.5.
- WM survival → `data/audit/wm_stratified/` · `audit_83_wm_stratified_cophenetic.py`
  / `audit_85_wm_decimation_control.py` / `audit_86_wm_grassmann_cluster_extent.py`
  · 2026-06-08…12.
- τ-robustness → `audit_121_tau_sweep_cophenetic_trace.py` · 2026-06-22; memory
  `tau_sensitivity_trace_2026_06_22`.
- Spectral head-to-head (N1.2b: α cophenetic-only; selectivity; β tie) →
  `data/audit/spectral_distance_swap/cohort_summary.csv` ·
  `audit_143_spectral_distance_swap_headtohead.py` · 2026-06-22.
- N1.6 sustained reinstatement (dynamical trace; held-not-flashed) →
  `data/audit/replay_states/{sustained_reinstatement_cohort,cohort_verdict,
  tau_resolved_cohort,directional_cohort,raw_propagator_cohort,subspace_cohort,
  normlap_cohort,subsecond_gamma_cohort}.csv` · `audit_124`–`141` · 2026-06-22/23.
  Full verdict ledger: `.agents/reports/2026-06-22_replay-states-verification.md`.
  Transient-replay investigation CLOSED (negative); sustained = the headline finding.

## §G — Missing parts / open

- **`task_learn`** as a trace probe in its own right (above) — the single
  biggest under-explored dimension of N1.
- **No behavioral anchor exists** — TI performance data is unavailable and will
  not be obtained (PI 2026-06-22); the "consolidation" framing rests on the
  task-phase decomposition (N2) + anatomy + literature, never a performance
  correlation.
- **Coordinated cross-phase null** (CORE §G) applies here too.
