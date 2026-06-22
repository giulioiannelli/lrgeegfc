---
name: headline-n1-trace
era: IMCOH_ABS_COHORT_N10
status: current
kind: headline
scope: N1 — the band-specific, multiscale persistence trace left by a transitive-inference task in post-task rest; β on both probes; over-expressed in OFC; carried by healthy gray-matter cortical coupling; emergent (not raw-strength); invisible to edge-wise comparison.
owner_agent: localization + white-matter + preprint-general-questioning(null)
updated: 2026-06-22
---

# N1 — A transitive-inference task leaves a multiscale connectivity trace that persists into rest, peaks in β, and over-expresses in orbitofrontal cortex

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
   edges moved around." Plain edge-by-edge comparison would miss it.

3. **It has an address and a tissue type.** The persistent β reorganization is
   brain-wide but **over-expresses in orbitofrontal cortex (OFC)** — a hub for
   building relational "cognitive maps," exactly what a transitive-inference task
   should leave behind — and it is carried by **healthy gray-matter cortical
   coupling**, *not* by the epileptic seizure core. So it is neither everywhere-
   and-meaningless nor a pathology artifact. **And the tissue rule is
   band-dependent** — β *spares* the diseased core while α *recruits* it (N1.6) —
   the bridge from this cognitive trace to the epileptogenic read-out (N3).

**Why it matters:** an offline, multiscale, cortical trace of a reasoning task,
sitting in the β band and the orbitofrontal map system — a candidate signature of
**memory/inference consolidation** read directly from intracranial connectivity.

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

## §G — Missing parts / open

- **`task_learn`** as a trace probe in its own right (above) — the single
  biggest under-explored dimension of N1.
- The **dynamic complement** — does this static signature recur as transient
  *states* in time-resolved rest? — is now headline **N4 (replay states)**; N1
  establishes the multiscale-multiband signature N4 searches for.
- **No behavioral anchor exists** — TI performance data is unavailable and will
  not be obtained (PI 2026-06-22); the "consolidation" framing rests on the
  task-phase decomposition (N2) + anatomy + literature, never a performance
  correlation.
- **Coordinated cross-phase null** (CORE §G) applies here too.
