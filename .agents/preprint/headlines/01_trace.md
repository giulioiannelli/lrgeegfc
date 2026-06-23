---
name: headline-n1-trace
era: IMCOH_ABS_COHORT_N10
status: current
kind: headline
scope: N1 — the band-specific, multiscale persistence trace left by a transitive-inference task in post-task rest; β on both probes; over-expressed in OFC; carried by healthy gray-matter cortical coupling; emergent (not raw-strength); invisible to edge-wise comparison.
owner_agent: localization + white-matter + preprint-general-questioning(null)
updated: 2026-06-23
---

# N1 — A transitive-inference task leaves a multiscale connectivity trace that persists into rest as a sustained state the brain holds, peaks in β, and over-expresses in orbitofrontal cortex

> Emerges from CORE. Numbers live in cached CSVs (§F). Raw FC is the baseline,
> never the result.

## §A — Result in plain language

Patients learned an order from pairs ("A beats B, B beats C, …") and were then
tested on pairs they had never seen, which they had to reason out (transitive
inference). We compared the brain's resting network organization **before** and
**after** the task. The finding: **the task reorganizes the network and the
reorganization does not fully wash out — it persists into the post-task rest.**

Three things make this a real result rather than a measurement artifact:

1. **It is band-specific.** The persistence is strongest and most complete in the
   **β rhythm** — and it is the *only* band our two complementary read-outs both
   confirm. θ shows nothing. So this is a selective signature, not a generic
   property of all connectivity. (The other bands tell a probe-dissociation story
   that *is* the methodology's punchline — see CORE / N2.)

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
   band-dependent** — β *spares* the diseased core while α *recruits* it (N1.6) —
   the bridge from this cognitive trace to the epileptogenic read-out (N3).

4. **It is a state the brain *holds*, not a fading echo.** Watching rest unfold in
   time (20 s windows), the post-task rest does not merely drift near the task
   pattern on average — it **dwells** there, window after window (**10/10 patients,
   p=0.001**), at *every* scale of the hierarchy. We tested hard for the opposite —
   that the trace flashes back in brief, replay-like *bursts* — across every
   representation, scale, band, and timescale down to the ripple range; there are
   **no bursts**. So the trace is a **sustained reinstatement** the brain settles
   into and holds offline, not a transient replay process (N1.7).

**Why it matters:** an offline, multiscale, cortical trace of a reasoning task,
sitting in the β band and the orbitofrontal map system, that the resting brain
**actively holds as a state** — the resting network does **not** return to where it
started (a non-ergodic, held reorganization), read directly from intracranial
connectivity. **This sets up the flagship (N2):** N1 establishes that the brain
*holds a structure offline* and that our method can *see* it; **N2 asks what that
held structure *is* — and finds it is an abstraction: the relations the brain
reasoned out, not just the pairs it saw.**

## §B — Technical statement (per subheadline; reference CSVs, no tables)

**N1.1 — The trace exists and is band-specific.** On the per-pair multiscale
cophenetic probe `ρ^coph`, the post-task rest sits closer to the task state than
the pre-task rest does, at the cohort level, for **β** (passes C1, C2, C3, C4)
and **α** (cophenet-only). On the global-mode Grassmann probe, **β** also passes
(cluster-extent permutation, Decision-8 mass gate, LOO-robust). β is the unique
"both-probes" band; θ/γ_high pass neither. *(C3 is the binding control; see CORE
§B.)*

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

**N1.3 — Concentrates in orbitofrontal cortex.** The brain-wide β cophenet trace
concentrates, **above each patient's own demeaned baseline**, in the OFC system
(matched-strength R=1000, rank concordance, canonical observed trace; BH-clears
across the 9 a-priori systems within β, in all four contact/shaft × epi-incl/excl
conditions; 4/5 implanted patients positive; leave-one-out + shaft-collapse
robust; low-strength/non-hub tissue; bilateral; peaks at the **system** scale,
washes out at lobe/hemisphere). It is a **hotspot on a distributed trace, not a
container** (within-OFC lean ≈ 0.59). Hippocampus/MTL = real but sub-threshold
hint; occipital = low-coverage (K=3) secondary; the Grassmann subspace probe does
**not** localize. This verdict **flipped twice before locking** — state that.

**N1.4 — Carried by healthy gray-matter coupling, not the epi core.** Stated in
the LRG framework: the genuine pair-count hotspot is **gray↔gray cortical
coupling** (WM↔WM depleted); pairs internal to the clinically-labelled seizure
core do **not** carry the trace (fail their matched-strength null; "core spared"
is directional only). Excluding SOZ contacts leaves β/α intact. **Two prior
sub-claims withdrawn (brutal honesty):** "epi-exclusion *strengthens* the trace"
(generic node-count, decimation control) and "the epi↔healthy interface is the
strongest carrier" (not a pair-count hotspot). The WM result is parallel: the
trace **survives** a gray-only montage but the apparent "sharpening" was
node-count, not WM-specific.

**N1.5 — τ-robust, fine-scale.** A τ-sweep (`audit_121`) confirms the trace is
fine-scale and stable around τ = 1/λ_max; coarse-τ "gains" are a collapse
artifact. The result is not a scale cherry-pick.

**N1.6 — BRIDGE to the clinic: β spares, α recruits the epileptic core.** The
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

**N1.7 — The trace is a SUSTAINED STATE the brain holds (dynamical reinstatement),
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
to make N1.7 fully standalone (§D). The directional-flow variant was tested and
**retracted** (magnitude-inherited, `audit_140`).

## §C — Critical issues & powerful strengths

**Airtight:** band-specificity; matched-strength pass on both probes (β);
emergent-not-strength; OFC localization survives the full six-rung gate incl.
shaft-collapse + R=1000 + LOO; volume-conduction-immune substrate; τ-robust.

**Lead-with weaknesses / provisional:**
- **OFC coverage is n=5** (irreducible) and the verdict flipped 3× historically —
  present as *concentration above baseline*, never "the trace lives in OFC."
- **Outliers / leverage:** Pat_15 (right-hemi-only, β-anti — sanctioned dropout),
  Pat_10 (β-anti), Pat_02 (β LOO argmax driver, partly epi-coupled). Report
  LOO-max; no single-patient "strong" tag.
- **`task_learn` was not used** for the core trace — "the task" = `task_test`
  only. This reference-phase choice must be justified (see §G); it is the bridge
  to N2.
- **θ is not simply "anti-trace"** — characterize, don't assert (README §5).
- **N1.7 sustained-state framing**: the *level* (rest_post dwells in the task config)
  is essentially N1 dynamically — its strength is robustness (5 representations, all
  scales) + the clean transient-negative boundary, **not** a new independent effect.
  The genuinely-beyond-static-N1 piece ("rest_post is a *tighter* attractor than
  rest_pre") is only a weak trend (A2b p=0.08) — do **not** claim "tighter," claim
  "held." Window-level matched-strength not yet run (inherits the static C3); see §D.

## §D — To-dos & verifiables

- [ ] (localization) Re-confirm OFC CSVs still match the locked verdict; produce
  the system-scale concentration figure with X-epi panels (not C4).
- [ ] (white-matter) Re-state the gray-matter-dominant pair-count result as the
  *significant* finding; keep the decimation-controlled withdrawals explicit.
- [ ] (null model) θ characterization; outlier-leverage panel (LOO-max per band).
- [ ] (null model) τ-robustness figure/supplement from `audit_121`.
- [ ] **(`task_learn` gap)** does `rest_pre → task_learn → rest_post` leave its
  own β trace, and does it localize to OFC too? Decide whether N1's "task" should
  be test-only, learn-only, or both — coordinate with N2.
- [ ] **(N1.7 optional bulletproof)** window-level matched-strength surrogate on the
  per-window `g` (make the sustained-state claim standalone, not inheriting C3).
  Cheap; only needed if a referee treats N1.7 as independent of static N1.

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
- ρ^coph C3 matched-strength →
  `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` ·
  `audit_63_split_baseline_surrogate.py` · 2026-05-15.
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
- N1.7 sustained reinstatement (dynamical trace; held-not-flashed) →
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
