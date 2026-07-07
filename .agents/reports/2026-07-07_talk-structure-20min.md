---
name: talk-structure-20min
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-07
pointers:
  - .agents/reports/2026-07-06_rho-sym-pipeline-migration.md
  - .agents/reports/2026-07-07_r2-mesoscale-encoding-rhosym-migration.md
  - .agents/preprint/headlines/02_encoding_vs_inference.md
  - data/reports/results_section1/
  - data/preprint/figures/all_bands/
---

# 20-minute talk — structure & figure plan

## Head

A 20-slide arc (7 intro + 13 methods/results/conclusion) that goes from
"the brain is a network" to our flagship result: a transitive-inference
task leaves a **held, multiscale connectivity trace** in post-task rest
that (N1) concentrates in **orbitofrontal cortex** in the **β band**, and
(N2, the climax) carries specifically **the relations the brain *inferred*,
not just the pairs it *saw*** — an offline abstraction of a learned order.
The method thread is the spine: imaginary coherence (kills volume
conduction) → per-band → Laplacian Renormalization Group diffusion
hierarchy (reads the mesoscale) → cophenetic cross-phase similarity vs a
matched-strength null. Every headline result already has a rendered PDF
on disk **except the N2 encoding-vs-inference figure**, which needs one
regeneration pass (§Gaps).

Pace: ~20 slides / 20 min ≈ 1 min/slide. Intro fast (~40 s each),
results slower (~75 s each). Budget ~2 min of slack for one demo/anecdote.

---

## Accuracy flags — lock these before building slides

1. **Use ρ_sym numbers, not ρ_split.** Current estimator is ρ_sym
   (symmetric split-half). If a draft slide shows **β p=0.005 / 23.7×**
   or **α p=0.002**, those are the older ρ_split values. Current:
   **β p=0.032, 16.6× null; α p=0.024, 13.7× null.** Verdicts are
   estimator-invariant (0/6 bands flip); ρ_sym is just more conservative.
2. **The N2 mesoscale is NOT a "peak."** Say **multiscale robustness**
   (β inference clears the null at fine *and* meso scales); the 13-point
   τ-sweep has no interior maximum. The "scale of multi-step integration"
   peak claim was dropped 2026-07-07. Keep relational-chaining as
   *motivation*, not as a located scale.
3. **Never rank bands by the ×null ratio.** γ_low's ratio is *larger*
   than β's yet it fails the cohort gate. ρ_coph is magnitude-blind;
   its scale is set by the null (~0.01), and the **matched-strength
   control** — not the ratio — certifies a band. Anchor "+0.20 is large"
   on the null (~0.01) and on per-patient tracers reaching 0.3–0.5.
4. **No behavioral data — state it once, plainly.** TI performance is
   unavailable and will not be obtained. "Abstraction/consolidation" is
   an interpretation licensed by task structure + anatomy + literature,
   never a brain–behavior correlation. Don't imply one exists.
5. **OFC is a hotspot, not a container.** Implanted in 5/10 (4 positive);
   bilateral; within-OFC ~59% pro. Present as *concentration above
   baseline*, not "the trace lives only in OFC."
6. **Don't feature specific heat C(τ)/entropy curves as "the method."**
   Our FC graphs are fully connected with a near-continuous spectrum →
   C(τ) has no informative peak here. The multiscale content comes from
   the **dendrogram hierarchy at τ=1/λ_max**, not from a τ-sweep of C(τ).

---

## PART I — INTRODUCTION (7 slides, ~5 min)

### I-1 · Title
*Message:* set the frame — "Multiscale network analysis of brain
connectivity: what a diffusion view reveals about learning and disease."
- Title, authors, affiliation, one-line hook.
- *Figure:* optional — the glass-brain cohort overlay as a backdrop.
  `data/outputs/figures/implant_in_brain/cohort_implants_overlay.pdf`

### I-2 · The brain is a network
*Message:* functional connectivity turns a recording into a graph;
cognition and epilepsy are *network* phenomena, so the graph is the object.
- Nodes = recording sites (sEEG contacts, ~110–122/patient); edges =
  functional coupling between two contacts' signals.
- We don't study channels, we study the organization *between* them.
- *Figure:* one FC adjacency matrix + its network.
  `data/outputs/figures/fc_templates/single_adjacency/` (pick a clean
  β/rest example) or `data/outputs/figures/network_templates/matrix_plus_network/`.

### I-3 · The problem of scale — why "mesoscale"
*Message:* the biology lives between the single edge and the global
average, in **nested communities** — and there is no single "correct"
number of communities.
- micro (one edge: noisy, misses organization) — macro (one number:
  throws structure away) — **meso (modules → systems → hierarchy).**
- Cognition/disease *reorganize* this mesoscale; a flat edge list can't see it.
- *Figure:* same network partitioned at two scales (5 vs 20 communities).
  `data/outputs/figures/section3/fig_I/fig_I2_brain_communities_ImCoh_n5_n20.pdf`

### I-4 · The idea: read all scales at once (the thesis)
*Message:* the Renormalization Group is physics' tool for "zooming out";
the **Laplacian RG** uses **diffusion as the zoom knob** — diffusion time
τ *is* the scale. One operator, all scales, a principled hierarchy.
- Short τ → local/fine structure; long τ → global/coarse. Sweep τ = sweep scale.
- This is the talk's promise: what does a *multiscale* view reveal that
  edge-by-edge FC cannot? (Answer arrives in Results.)
- *Figure:* schematic (build in slides) — a graph with a heat blob
  spreading at increasing τ, resolving into a dendrogram. Optionally pair
  with a real dendrogram: `data/outputs/figures/network_templates/circular_dendrogram_network/`.

### I-5 · The paradigm: transitive inference & cognitive maps
*Message:* patients learn an ordered chain (A>B>C>D…) then judge **novel**
pairs (B>D) they never saw — answerable only by *inferring* the latent
order. The brain builds relational **cognitive maps** in **OFC**, so a
relational trace should land there.
- Encoding (premises shown) vs inference (novel pairs reasoned out) — the
  distinction the flagship result exploits.
- Literature home: OFC as the substrate for task/state-space maps
  (Wilson–Niv–Schoenbaum; Behrens/Schuck geometry-of-abstraction).
- *Figure:* schematic (build in slides) — the chain A>B>C>D and the
  inferred B>D.

### I-6 · The experiment
*Message:* 10 stereo-EEG epilepsy patients, four phases bracketing the
task with rest.
- Cohort n=10; depth electrodes (11–14 probes); 2048 Hz (Pat_03 1024 Hz,
  handled at config, a full member).
- **rest_pre → task_learn → task_test → rest_post** — the two rests bracket
  learning + inference.
- *Figure:* the 4-phase timeline (build in slides) + cohort implants.
  `data/outputs/figures/implant_in_brain/cohort_implants_overlay.pdf`

### I-7 · The question
*Message:* does the task leave a **band-specific, held, multiscale** trace
in rest_post? A *joint* event — reorganization occurs during task **and**
is still visible after it ends.
- Vocabulary (state once): **trace** (changed & stuck) · **anchor**
  (never changed) · **reset** (changed & reverted) · **emergent** (new).
  Our claim is a *trace*.
- Elegance: the task is *learning a hierarchy*; the measure is *a persisting
  change in a connectivity hierarchy*. Form mirrors form.
- **Caveat up front:** no behavioral data — the read is offline
  consolidation by structure+anatomy, never performance.

---

## PART II — METHODS (4 slides, ~4 min)

### M-1 · Why imaginary coherence (the killer methods slide)
*Message:* intracranial FC is dominated by **volume conduction** —
instantaneous, zero-lag leakage. Imaginary coherence is **blind to it by
construction**, because a zero-lag component is purely real and its
imaginary part is exactly zero.
- Standard coherence (MSC) uses the full cross-spectrum → captures the
  artifact → same-probe coupling **2–8× inflated**, driving the hierarchy
  by electrode geometry, not brain function.
- ImCoh = Im(S_ij)/√(S_ii·S_jj) ∈ [−1,1] (Nolte 2004). We feed
  **imcoh_abs = ⟨|ImCoh(f)|⟩** over in-band bins (non-negative for the
  Laplacian; |·| per-bin *then* average, Jensen).
- Number to quote: same/cross-probe ratio **5.5× (MSC) → 1.13× (imcoh_abs)**.
- *Figure:* MSC-vs-ImCoh communities or dendrograms on the brain.
  `data/outputs/figures/section3/fig_I/fig_I1_brain_communities_MSC_n5_n20.pdf`
  (pair with I2 from I-3) **or**
  `data/outputs/figures/section3/fig_H/fig_H1_dendrograms_MSC_vs_ImCoh_3patients.pdf`.

### M-2 · Why frequency bands
*Message:* different rhythms index different processes; **band-specificity
is both the signal and a control** — a real effect is frequency-selective,
an artifact is band-blind.
- Table (flag γ_high = 80–300 Hz, intracranial, **not** textbook):
  δ 0.53–4 · θ 4–8 · α 8–13 · β 13–30 · γ_low 30–80 · γ_high 80–300 Hz.
- Broadband would average the signal band (β) with the null band (θ) and
  wash it out.
- Foreshadow: at the *raw edge* level the drift is positive in *every*
  band incl. θ (band-blind) — specificity only *emerges* through the
  multiscale read (→ Result R1).
- *Figure:* a 6-band row of FC matrices.
  `data/outputs/figures/fc_templates/row_per_band/` (one patient/phase).

### M-3 · The LRG pipeline — from matrix to hierarchy
*Message:* build the diffusion operator, let heat flow for time τ, and read
the resulting **communication hierarchy**. The heat kernel sees the *whole
web* (all indirect paths), not one edge.
- Chain (keep it visual, minimal algebra):
  Â (imcoh_abs) → L̂ = D̂ − Â → K̂(τ)=e^{−τL̂} → ρ̂(τ)=K̂/Z →
  communication distance D_ij(τ)=(1−δ_ij)/K_ij (ultrametric, Villegas 2025)
  → UPGMA **dendrogram** at τ=1/λ_max → **cophenetic distance** = merge
  height of the lowest common ancestor.
- "Multiscale" here = the dendrogram's own nested cuts (we commit to one τ;
  see caveat — don't show C(τ)).
- Grounding: Villegas et al., Nat. Phys. 2023 (L̂,K̂,ρ̂); PRR 2025
  (communication distance, ultrametric, UPGMA).
- *Figure:* matrix → network → dendrogram in one row.
  `data/outputs/figures/network_templates/matrix_plus_network/` +
  `data/outputs/figures/section3/fig_H/fig_H1_dendrograms_MSC_vs_ImCoh_3patients.pdf`.

### M-4 · How we test for a trace — and the null that matters
*Message:* similarity of two phases = **cophenetic correlation ρ_coph**
(Spearman on cophenetic distances). A trace must clear a **matched-strength
surrogate null** — the mandatory control that killed our earlier headline.
- ρ_sym = symmetric split-half estimator (fixes sign flips on near-zero
  patients); scale set by the null (~0.01), not by 1.
- The gate is the **Wilcoxon vs matched-strength**, not patient counts,
  not the ×null ratio. Cross-phase drift alone ≠ trace.
- *Figure:* signal-vs-null triangle per band.
  `data/preprint/figures/all_bands/fig_bands_null_triangle_coph.pdf`
  (optionally beside `..._null_triangle_rawfc.pdf` to preview the
  raw-vs-multiscale contrast).

---

## PART III — RESULTS (8 slides, ~10 min)

### R-1 · The reveal: band-specificity is *emergent*
*Message:* raw edge-level drift is positive in **every** band (θ +0.12,
indistinguishable) and raw β doesn't even clear the null (p=0.053). Read
through the diffusion hierarchy, the non-tracing bands **collapse** (θ →
−0.04) and **β is retained and now clears** (p=0.032). The multiscale read
*earns its keep* here.
- Guardrail: the hierarchy **compresses** magnitude — it separates, it
  doesn't amplify. So this is filtering, not inflation.
- *Figure:* `data/reports/results_section1/fig_trace_c_raw_vs_multiscale.pdf`.

### R-2 · The gate: which bands hold the trace
*Message:* cohort-level, **β and α** hold it; γ_low/δ/γ_high fail; **θ is
the clean negative.** Peak ≠ verdict.
- ρ_sym: **β +0.20, 6/10, p=0.032, 16.6×** · **α +0.10, 6/10, p=0.024,
  13.7×** · γ_low +0.08 p=0.080 (fails) · θ −0.04 p=0.784 (anti).
- γ_low carries the *largest single-patient* traces (Pat_05 +0.86 vs β max
  +0.54) yet fails the cohort gate → why we certify by the null, not the peak.
- *Figure:* `data/reports/results_section1/fig_trace_a_band_forest.pdf`
  (per-band per-patient forest) **or**
  `data/reports/rho_sym_band_map/fig_rho_sym_band_map.pdf`.

### R-3 · β → orbitofrontal cortex
*Message:* the β trace **concentrates in OFC** — the cognitive-map hub —
and OFC is the *only* enrichment survivor; sensorimotor & PFC are depleted.
- BH **q = 0.010–0.015** (R=1000) across 9 a-priori systems, in **all four**
  conditions (contact/shaft × epi-in/out); LOO-robust.
- Hotspot not container: bilateral; within-OFC ~59% pro; OFC in 5/10 (4+).
  The one marginal LOO (q=0.060) is the known n=5 coverage fragility.
- *Figure:* `data/reports/results_section1/fig_trace_b_ofc_localization.pdf`
  (glass brain) or polished
  `data/preprint/figures/beta/anatomy/fig_beta_anatomy_brain.pdf`.

### R-4 · It's real: two independent probes, strength-immune
*Message:* β is the **only** band cleared by **both** an ultrametric-geometry
probe (cophenetic) and a global-mode probe (**Grassmann**) — two different
strength-independent lenses agreeing.
- Grassmann β cluster-mass **p=0.005, LOO 0.005**; α fails Grassmann
  (0.348); γ_low passes Grassmann only; θ,γ_high none.
- Band × probe dissociation is the core "why you should believe it" panel.
- *Figure:* `data/preprint/figures/all_bands/fig_bands_coph_grassmann_dissociation_map.pdf`
  (optionally `data/preprint/figures/beta/grassmann/fig_beta_grassmann_heatmap.pdf`).

### R-5 · Held, not replayed — and who carries it
*Message:* post-task rest **dwells** in the task configuration (not a
transient replay), and the cohort splits into a stable structure.
- Reinstatement: post>pre task-likeness **10/10, p=0.001, LOO 0.002**, at
  every scale; transient/replay is a clean negative down to the ~0.2 s
  ripple floor. Say "held," never "tighter" (that's only a trend).
- 3 tiers: **6 stable tracers** (Pat_08 +0.54 … Pat_07 +0.13) · **2
  undetermined** (Pat_13, Pat_15) · **2 resets** (Pat_10 −0.14, Pat_14).
  (Pat_15 is now *undetermined*, no longer a clean β-anti dropout.)
- *Figure:* `data/reports/results_section1/fig_trace_e_reinstatement.pdf`
  (+ `data/outputs/figures/section_5_lrg_trace/headline/per_patient_slopes.pdf`
  for the tiers).

### R-6 · The flagship (N2): encoding vs inference
*Message:* the four-phase design lets us split the trace into what was
**shown** (encoding, task_learn) and what was **reasoned out** (inference-
specific, task_test with encoding controlled out) — and ask which one persists.
- Decomposition (one clean line): e = D_learn−D_pre (encoding),
  f = D_test−D_learn (inference-specific), p = D_post−D_pre (persistence);
  T_infspec = partialcorr(f, p | e).
- Motivation (not a scale claim): relational chaining A>B,B>C ⇒ A>D
  integrates multi-step paths — why we look across scales.
- *Figure:* schematic (build in slides) — the 4-phase design mapped to
  e / f / p. **Rendered N2 figure needs regeneration — see §Gaps.**

### R-7 · The flagship (N2), result: the *inference* persists — β only
*Message:* the inference-specific component **itself persists into rest** —
in **β and no other band** — and the consolidation **anchors in OFC**,
set by learning. This is the offline abstraction of a learned order.
- β inference-specific persistence: **+0.091 vs null +0.005, p=0.0098,
  LO-P15 0.020, 6/10**; **multiscale-robust** (also τ≈2.6, p=0.0049 — not
  a peak); **duration-immune** (vs length ρ=+0.10, p=0.78); every other
  band null (next p≈0.25).
- Encoding→OFC both phases (q=0.010 excl / 0.040 incl, PFC depleted);
  cingulate multiplexes — **low-γ → cingulate memory** (q=0.035, 8/8, solid);
  inference→cingulate is a **directional lead only** (fails length-matching).
- *Figure:* `data/audit/inference_localization/figures/systems_inference_dissociation.pdf`
  or `.../brain_inference_dissociation_beta.pdf` (legacy — regenerate for
  polish); summarize the taxonomy with
  `data/reports/per_band_phenomenology/fig_per_band_phenomenology.pdf`.

### R-8 · Coda (N3): the same operator reads epileptogenic tissue
*Message:* the identical diffusion view doubles as a **clinical marker** —
seizure contacts form a strength-independent **diffusion community**, and
the cognitive bands split over it: **β spares** the core, **α recruits** it.
- Co-diffusion AUC (matched-strength gate): **δ 0.80, γ_low 0.74, β 0.69**;
  off-shaft δ LOSO **0.72, 8/10**. Seed-based 6-band detector **mean AUC
  0.81, precision@5 60% (~7× chance), 9/10** above chance.
- Scope (state it): clinical SOZ labels, **not** surgical outcome.
- *Figure:* `data/reports/results_section1/fig_trace_f_soz_divergence.pdf`
  (β spares / α recruits — the cognitive→clinical bridge) +
  `data/preprint/figures/all_bands/fig_epi_compound_probability.pdf` or
  `fig_epi_soz_marker.pdf`.
- **Optional/cuttable** if time is tight — the talk stands on N1+N2 alone.

---

## PART IV — CONCLUSION (1 slide, ~1 min)

### C-1 · Take-homes
*Message:* a multiscale diffusion view of intracranial FC reveals what
edges and averages miss: a **held, band-specific hierarchical trace** of
inference, in **OFC**, plus a clinical read-out for free.
- **Method:** ImCoh (no volume conduction) × per-band × LRG diffusion
  hierarchy × matched-strength null — band-specificity *emerges* from the
  multiscale read.
- **N1:** β cophenetic trace, cohort-level, → OFC (two probes, strength-immune).
- **N2 (flagship):** the *inferred* relations persist offline in β — an
  abstraction of a learned order (interpretation, no behavior).
- **N3:** same operator localizes epileptogenic tissue.
- Outlook: behavioral link, larger OFC coverage, prospective clinical test.
- *Figure:* reuse the dissociation map or per-band phenomenology as a
  one-glance summary.

---

## Figure gaps — what to produce before Sunday

**Must-fix (1 regeneration):**
- **N2 encoding-vs-inference rendered figure.** Scripts exist
  (`scripts/01_compute/audit/audit_158_encoding_localization_rhosym.py`,
  `audit_160_inference_localization_rhosym.py`, `audit_157/159/161_*`) but
  their ρ_sym output dirs are empty/missing
  (`data/audit/inference_localization_rhosym/` has 0 figures). The only N2
  figures on disk are **legacy (pre-ρ_sym)**:
  `data/audit/inference_localization/figures/brain_inference_dissociation_beta.pdf`,
  `systems_inference_dissociation.pdf`. For the flagship slide, regenerate a
  clean ρ_sym version (or, if time-boxed, use the legacy figure and label
  numbers from R-7 by hand).

**Build in the slide editor (trivial, no code):**
- 4-phase timeline (rest_pre → learn → test → rest_post) — for I-6/R-6.
- TI task schematic (chain A>B>C>D; inferred B>D) — for I-5.
- LRG "diffusion zoom" schematic (heat blob → dendrogram) — for I-4/M-3.
- e / f / p decomposition mapped onto the 4 phases — for R-6.

**Deliberately NOT shown:** specific-heat C(τ)/entropy curves — flat for our
fully-connected spectrum, not load-bearing (would invite a wrong Villegas
analogy). Only exploratory PNGs exist anyway.

**Nice-to-have if time:** encoding→OFC vs inference→cingulate double-
dissociation brain (not on disk); `fig_alpha_anatomy_brain.pdf` (α has no
anatomical home — could be shown as a *negative*).

---

## Verified figure appendix (repo-relative paths)

**Headline N1 set** — `data/reports/results_section1/`
- `fig_trace_a_band_forest.pdf` — band gate (which bands hold)
- `fig_trace_b_ofc_localization.pdf` — β trace → OFC on the brain
- `fig_trace_c_raw_vs_multiscale.pdf` — the emergent-not-strength reveal
- `fig_trace_d_tissue_class.pdf` — gray/WM/epi carrier
- `fig_trace_e_reinstatement.pdf` — held-state (post dwells in task config)
- `fig_trace_f_soz_divergence.pdf` — β spares / α recruits the SOZ

**ρ_sym band map** — `data/reports/rho_sym_band_map/`
- `fig_rho_sym_band_map.pdf`, `fig_rho_sym_detectability_decomposition.pdf`,
  `fig_raw_vs_cophenet_resolving_power.pdf`

**Per-band phenomenology** — `data/reports/per_band_phenomenology/fig_per_band_phenomenology.pdf`

**Manuscript tree** — `data/preprint/figures/`
- `all_bands/fig_bands_coph_grassmann_dissociation_map.pdf` — core band×probe panel
- `all_bands/fig_localization_matched_strength.pdf`, `fig_localization_beta_two_null.pdf`
- `all_bands/fig_bands_null_triangle_{coph,rawfc}.pdf` — signal vs null
- `all_bands/fig_epi_*.pdf` (18 files: `fig_epi_soz_marker.pdf`,
  `fig_epi_compound_probability.pdf`, `fig_epi_distant_soz_marker.pdf`,
  `fig_epi_occult_discovery.pdf`, `fig_epi_stratified_*` …)
- `beta/anatomy/fig_beta_anatomy_brain.pdf`, `beta/grassmann/fig_beta_grassmann_heatmap.pdf`,
  `beta/dendrogram_persistence/fig_beta_dendrogram_persistence_Pat_05_*.pdf`

**Method-context figures** — `data/outputs/figures/`
- `implant_in_brain/cohort_implants_overlay.pdf` — cohort implants (intro)
- `section3/fig_I/fig_I{1,2}_brain_communities_{MSC,ImCoh}_n5_n20.pdf` — MSC vs ImCoh, two scales
- `section3/fig_H/fig_H1_dendrograms_MSC_vs_ImCoh_3patients.pdf` — dendrograms differ
- `section2/fig_G/fig_G1_eigenvalue_spectrum_Pat_0{2,5,8}_MSC_vs_ImCoh.pdf` — spectra
- `fc_templates/single_adjacency/*.pdf`, `fc_templates/row_per_band/*.pdf`, `fc_templates/row_per_phase/*.pdf`
- `network_templates/{matrix_plus_network,circular_dendrogram_network}/` — Pat_05 β exemplars
- `section6/task_trace_band_k_n10_imcoh_abs.pdf` — band×k cohort map

**Legacy N2 (pre-ρ_sym, regenerate)** — `data/audit/inference_localization/figures/`
- `brain_inference_dissociation_beta.pdf`, `systems_inference_dissociation.pdf`
- `data/audit/consolidation_arc/figures/arc_tau_sweep.pdf` — τ sweep
