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

**Design: a "wonder-talk" arc — open on the puzzle, not the field. NO
generic "brain as a network." Lead with the inference you were never told,
then turn our method into its answer. Every beat hangs on a precise
citation (full refs + DOIs in "Intro references" at the end — DOI pass in
progress). Two threads braid on I-6: the cognition is hierarchy-shaped, the
method is hierarchy-shaped — form mirrors form. (Alternative ordering for a
network-physics audience: method-first — I-4/I-5 scale+diffusion, then TI as
the perfect testbed. Default below is cognitive-first.)**

### I-1 · Cold open — "the inference you were never told"
*Message / spoken hook:* "I show you Anna > Ben, Ben > Carla, Carla > Dan.
I never say a word about Anna and Dan. Yet you already know. Where did that
knowledge come from — and where does it go when you stop?"
- No data, no field intro — just the puzzle and an evocative title.
- *Title options:* **"The shape of an inference"** · "Offline abstraction" ·
  "What the brain keeps after it stops thinking."
- *Figure:* none (or one faint chain graphic). Let the words land.

### I-2 · That was abstraction, not memory
*Message:* you didn't store pairs — you built a *structure* (an order) and
read a new fact off it. That is a **cognitive map**: knowledge organized for
inferences you were never trained on.
- *Attach:* Behrens et al. 2018, *Neuron* ("What is a cognitive map?"); the
  abstraction-geometry line — Bernardi et al. 2020, *Cell*; Whittington et
  al. 2020, *Cell* (Tolman–Eichenbaum Machine). TI as the classic probe:
  Dusek & Eichenbaum 1997, *PNAS*.
- *Figure:* schematic (build) — pairs shown → a line assembled → B>D read off.

### I-3 · Two precise predictions — where, and when
*Message:* state our whole experimental bet as literature-derived,
falsifiable predictions (not vibes).
- **Where:** a relational map should live in **orbitofrontal cortex** →
  Wilson et al. 2014, *Neuron*; Schuck et al. 2016, *Neuron* (OFC = map of
  task/state space).
- **When:** the map is cemented **offline, after** the task — post-task rest
  connectivity carries the just-learned structure → Tambini, Ketz & Davachi
  2010, *Neuron*.
- So we look in **rest_post**, in **frontal cortex**, for a persisting hierarchy.
- *Figure:* schematic (build) — an OFC pin + a "look after the task" arrow
  onto the phase timeline.

### I-4 · The catch — a map is multiscale (the scale slide)
*Message:* a learned order is **nested** (items → groups → the whole chain)
and there is **no privileged scale** to read it at. Standard FC forces a
choice: single edges (finest), or a partition at a *k* you had to pick.
- *Attach:* **Betzel & Bassett 2017, *NeuroImage* ("Multi-scale brain
  networks")** — brain organization spans scales with no correct resolution;
  Bassett & Sporns 2017, *Nat. Neurosci.* (network neuroscience).
- *Figure:* same network partitioned at two scales (5 vs 20 communities).
  `data/outputs/figures/section3/fig_I/fig_I2_brain_communities_ImCoh_n5_n20.pdf`

### I-5 · The turn — borrow renormalization from physics
*Message:* let **heat diffuse** across the connectivity graph; diffusion
time **τ becomes a continuous zoom** from local to global — coarse-graining
*by diffusion*. This **Laplacian Renormalization Group** turns a connectivity
matrix into a **tree** — an ultrametric hierarchy — read at all scales at
once. Spoken line: "a hierarchy is exactly what diffusion sees."
- *Attach:* Villegas et al. 2023, *Nature Physics* (LRG); lineage —
  De Domenico & Biamonte 2016, *PRX* (spectral-entropy / density-matrix);
  De Domenico 2017, *PRL* (diffusion geometry → clusters).
- *Figure:* schematic (build) — heat blob spreading at increasing τ →
  dendrogram; optionally
  `data/outputs/figures/network_templates/circular_dendrogram_network/`.

### I-6 · Form mirrors form (thesis / emotional peak)
*Message:* braid the threads — the cognition is a learned **hierarchy**; our
method **reads** hierarchy. We measure a **tree-shaped change** in
connectivity that the task induces and rest keeps — and can ask whether it
carries the pairs **shown** or the relations **inferred**. That
correspondence is *why this works*.
- *Figure:* two-panel concept — a learned-order tree beside a connectivity
  dendrogram, same shape.

### I-7 · The stage — rare data, a design built to catch it
*Message:* now, and only now, the concrete setup.
- **sEEG** — direct intracranial human recordings, the rare mm/ms window;
  **n=10**; depth electrodes (11–14 probes; 2048 Hz, Pat_03 1024 Hz at
  config, a full member).
- Four phases **rest_pre → task_learn → task_test → rest_post** — the bracket
  that lets us subtract encoding from inference and catch what *persists*.
- Restate the precise question: a **band-specific, held, multiscale** trace
  in rest_post — a *joint* event (reorganizes during task AND survives it).
- Vocabulary (once): **trace** (changed & stuck) · **anchor** (never
  changed) · **reset** (changed & reverted) · **emergent** (new). Ours is a
  *trace*.
- **One caveat, once:** no behavioral data → "abstraction/consolidation" is
  read from structure + anatomy + this literature, never from performance.
- *Hands off into M-1:* "we have the data — but can we even measure
  connectivity without the intracranial artifact?"
- *Figure:* 4-phase timeline (build) +
  `data/outputs/figures/implant_in_brain/cohort_implants_overlay.pdf`.

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

---

## Intro references (attaches by slide)

Full citations, verified against Crossref (authoritative DOI metadata).
Links are canonical `doi.org` resolvers.

**I-2 — abstraction / cognitive map (the wonder):**
- Behrens et al. (2018). What is a cognitive map? Organizing knowledge for flexible behavior. *Neuron* 100:490–509. https://doi.org/10.1016/j.neuron.2018.10.002
- Bernardi et al. (2020). The geometry of abstraction in the hippocampus and prefrontal cortex. *Cell* 183:954–967. https://doi.org/10.1016/j.cell.2020.09.031
- Whittington et al. (2020). The Tolman–Eichenbaum Machine: unifying space and relational memory through generalization in the hippocampal formation. *Cell* 183:1249–1263. https://doi.org/10.1016/j.cell.2020.10.024
- Dusek & Eichenbaum (1997). The hippocampus and memory for orderly stimulus relations. *PNAS* 94:7109–7114. https://doi.org/10.1073/pnas.94.13.7109 *(classic TI, rodent)*
- Heckers et al. (2004). Hippocampal activation during transitive inference in humans. *Hippocampus* 14:153–162. https://doi.org/10.1002/hipo.10189 *(human TI counterpart)*

**I-3 — where = OFC:**
- Wilson, Takahashi, Schoenbaum & Niv (2014). Orbitofrontal cortex as a cognitive map of task space. *Neuron* 81:267–279. https://doi.org/10.1016/j.neuron.2013.11.005
- Schuck, Cai, Wilson & Niv (2016). Human orbitofrontal cortex represents a cognitive map of state space. *Neuron* 91:1402–1412. https://doi.org/10.1016/j.neuron.2016.08.019

**I-3 — when = offline consolidation in rest:**
- Tambini, Ketz & Davachi (2010). Enhanced brain correlations during rest are related to memory for recent experiences. *Neuron* 65:280–290. https://doi.org/10.1016/j.neuron.2010.01.001
- Tambini & Davachi (2019). Awake reactivation of prior experiences consolidates memories and biases cognition. *Trends Cogn. Sci.* 23:876–890. https://doi.org/10.1016/j.tics.2019.07.008
- Tse et al. (2007). Schemas and memory consolidation. *Science* 316:76–82. https://doi.org/10.1126/science.1135935 *(prior schema governs offline consolidation)*

**I-4 — multiscale / the scale problem:**
- Betzel & Bassett (2017). Multi-scale brain networks. *NeuroImage* 160:73–83. https://doi.org/10.1016/j.neuroimage.2016.11.006
- Bassett & Sporns (2017). Network neuroscience. *Nat. Neurosci.* 20:353–364. https://doi.org/10.1038/nn.4502

**I-5 — diffusion / Laplacian renormalization group (the method):**
- Villegas, Gili, Caldarelli & Gabrielli (2023). Laplacian renormalization group for heterogeneous networks. *Nature Physics* 19:445–450. https://doi.org/10.1038/s41567-022-01866-8
- Villegas, Gabrielli, Poggialini & Gili (2025). Multi-scale Laplacian community detection in heterogeneous networks. *Phys. Rev. Research* 7:013065. https://doi.org/10.1103/PhysRevResearch.7.013065
- De Domenico & Biamonte (2016). Spectral entropies as information-theoretic tools for complex network comparison. *Phys. Rev. X* 6:041062. https://doi.org/10.1103/PhysRevX.6.041062
- De Domenico (2017). Diffusion geometry unravels the emergence of functional clusters in collective phenomena. *Phys. Rev. Lett.* 118:168301. https://doi.org/10.1103/PhysRevLett.118.168301

**M-1 — imaginary coherence (FC measure):**
- Nolte et al. (2004). Identifying true brain interaction from EEG data using the imaginary part of coherency. *Clin. Neurophysiol.* 115:2292–2307. https://doi.org/10.1016/j.clinph.2004.04.029
- Bastos & Schoffelen (2016). A tutorial review of functional connectivity analysis methods and their interpretational pitfalls. *Front. Syst. Neurosci.* 9:175. https://doi.org/10.3389/fnsys.2015.00175

**Optional add-ons:**
- Park, Miller & Boorman (2021). Inferences on a multidimensional social hierarchy use a grid-like code. *Nat. Neurosci.* 24:1292–1301. https://doi.org/10.1038/s41593-021-00916-3 *(neural code for an inferred hierarchy)*
- Ellenbogen et al. (2007). Human relational memory requires time and sleep. *PNAS* 104:7723–7728. https://doi.org/10.1073/pnas.0700094104 *(human relational inference + offline consolidation in one cite)*

---

## Handoff — parallel-worktree workflow & the interactive HTML deck

**Status (2026-07-09):** this file is the single source of truth for the
20-min talk. Structure is locked — intro rebuilt as the "wonder" arc
(I-1…I-7); methods (M-1…4), results (R-1…8), conclusion (C-1) complete.
Remaining polish: (a) ✅ verified DOIs folded into "Intro references"
(Crossref-checked); (b) regenerate the N2 encoding-vs-inference figure
(see Gaps); (c) build the interactive HTML deck (below).

**Workflow (travel-ready):**
- Presentation work happens on branch `presentation/talk-html-deck` in a git
  **worktree**, in parallel with the main `audit/cohort-n10-diagnostic`
  checkout on this laptop.
- The branch is **pushed to the remote**; on the travel laptop, clone the
  repo and check out that branch to continue seamlessly.

**Deliverable to build in the worktree chat — interactive HTML deck:**
- Self-contained and portable (must open offline on the travel laptop, no
  server / no network). Suggested: a single-file **reveal.js** deck with
  JS/CSS inlined, or plain HTML/CSS/JS. Everything embedded.
- 20 slides following the arc in this file (7 intro + 13). Speaker notes =
  the *Message* lines here.
- Figures: current assets are **PDF (vector)**. For the web, convert the
  chosen headline PDFs to **SVG** (keeps vector) or high-res PNG and embed.
  Headline set: `fig_trace_c` (raw-vs-multiscale), `fig_trace_a` (band gate),
  `fig_trace_b` (OFC), `fig_bands_coph_grassmann_dissociation_map`,
  `fig_trace_e` (reinstatement), the N2 inference figure (once regenerated),
  `fig_epi_*` + `fig_trace_f` (coda), `fig_I2` (multiscale communities),
  `fig_H1` (MSC vs ImCoh). Paths in the figure appendix above.
- "Cool/interactive" ideas (pick a few — don't overload):
  - I-5: a **τ-slider** morphing a network → dendrogram (the diffusion zoom).
  - I-1: the A>B>C>D chain assembling, then **B>D lighting up**.
  - R-2: hover a band → its per-patient forest / gate p.
  - R-3: rotatable/zoomable OFC glass-brain (or an animated reveal).
  - I-6: two trees (learned order ‖ connectivity dendrogram) snapping to the
    **same shape** — the "form mirrors form" beat.
- Keep the **accuracy flags** (top of this file) authoritative: ρ_sym numbers
  (β p=.032/16.6×, α p=.024/13.7×), N2 = robustness *not* a peak, no-behavior
  caveat, OFC hotspot-not-container.

**Open decisions to resolve while building:**
- Opening order: cognitive-first (default) vs method-first (network-physics
  audience). See Part I design note.
- Cold-open for I-1: transitive-inference story (recommended) vs method-
  wonder ("watch a thought become geometry").
- N3 epilepsy coda (R-8): keep as the "same tool, second payoff" bridge, or
  cut for a deeper-cognition talk. Cuttable without breaking N1+N2.
