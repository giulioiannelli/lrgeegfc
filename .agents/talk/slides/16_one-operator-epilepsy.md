---
name: talk-slide-16-one-operator-epilepsy
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
slide: 16
status: draft
updated: 2026-07-13
canva: NEW — merges old 19 (band dissociation) + old 20 (epilepsy marker). α-recruits-SOZ SOFTENED
  (β spares the SOZ is the confirmed result; α does not clearly recruit on the backbone).
---

# Slide 16 — One operator, two readouts: cognition and the seizure zone

1. TITLE
One operator, two readouts — cognition in α/β, epilepsy in δ/low-γ/β, β the bridge

2. MAIN CONCEPT
- **Head (the juice):** the *same* diffusion operator that carries the cognitive trace — nothing
  retuned — localizes epileptogenic tissue, read inside a single recording as diffusion affinity from a
  seed. And across bands the two jobs dissociate cleanly: **cognition lives in α/β, epilepsy in
  δ/low-γ/β, and β is the bridge**. A cognitive lens that pays a clinical dividend for free.
- **The band dissociation (one glance, don't go band-by-band):**
  - **cognition = α · β** — α the single-scale carrier, β the scale-invariant one (held, coarse-left,
    SOZ-independent — tissue pair-class + node-exclusion controls).
  - **epilepsy = δ · low-γ · β** — strongest in δ; δ marks the seizure zone but carries **no** cognitive
    trace, the cleanest face of the split.
  - **β = the bridge** — the one band in both families.
  - **θ silent** in both — the selectivity control: a global artefact would light every band.
- **The seizure-zone marker.** The onset contacts form a **strength-independent co-diffusing community** —
  relational, not "where connectivity is strong": it beats a matched-strength fake-SOZ null. Seed-affinity
  AUC, strongest in δ: **δ 0.83 · low-γ 0.82 · β 0.745 · α 0.61**, each above its own fake-SOZ null in most
  patients (δ 7/10, low-γ 8/10, β 8/10).
- **Over the seizure zone the cognitive bands part ways.** **β spares the SOZ** (its trace lives on
  healthy cortex — three ways: pair-class, node-exclusion, size-matched decimation). α does **not** clearly recruit the SOZ on the backbone (the old
  "α recruits +0.41" was a whole-graph effect — softened). So the cognitive trace and the disease marker
  sit on *different tissue*.
- **What the sweep buys, honestly.** τ lifts **ranking, not precision**: multiscale beats single-scale AUC
  in 7–8/10 patients, but top-of-list precision comes from **fusing the bands**, not from τ — single-band
  prec@5 = 0.40, **fused prec@5 = 60% (~7× base rate)**.
- **Two populations + a calibrated detector.** 8 patients are the co-diffusing community (community-level
  AUC up to **0.99**); 2 are right-hemisphere hub cases where the community read fails (0.49 / 0.57) and node
  strength rescues them. A **calibrated leave-one-patient-out detector** lands at **median AUC 0.87, 9/10
  above chance**. The marker is **robust to sparsification** (dense → mst@0.20: δ .80→.83, low-γ .74→.82,
  β .69→.745) — a property of the operator, not the thinning.
- **Scope:** this reads clinically-labelled **SOZ, not surgical outcome** — a **triage marker, not a
  diagnostic**.

3. ON-SLIDE TEXT
same operator, read within one recording → the seizure zone

band dissociation:  cognition = α · β   ·   epilepsy = δ · low-γ · β   ·   β = the BRIDGE
  θ silent (both) = the control · δ marks the SOZ but carries NO cognitive trace

SOZ = a strength-INDEPENDENT co-diffusing community (relational, beats a fake-SOZ null)
seed-affinity AUC:  δ 0.83 · low-γ 0.82 · β 0.745 · α 0.61   (above own null: δ7 · γ_low8 · β8 /10)
β SPARES the SOZ (cognitive ≠ epileptic tissue) · α does not clearly recruit it (backbone)
τ = ranking, not precision → band FUSION gives prec@5 = 60% (~7× base)
two populations: 8 community (AUC↑.99) + 2 right-hemi hubs → LOPO detector median AUC 0.87 (9/10 > chance)
robust dense→mst@0.20 · scope: SOZ labels, not surgical outcome (triage, not diagnostic)

4. SPEECH
And the same tool, with nothing retuned, pays a second dividend — a clinical one. But first the picture
that ties the whole talk together: read across every band, the one operator sorts the frequencies into two
families. Cognition lives in alpha and beta; epilepsy in delta, low-gamma, and beta; and beta is the bridge
between them. Delta is the sharpest case — it marks the seizure zone but carries no cognitive trace at all,
and theta is silent in both, which is our control that none of this is a global artefact. Now the clinical
read. Inside a single recording, the seizure-onset contacts light up as their own co-diffusing community —
and crucially it's strength-independent, it beats a null that scrambles the seizure zone but keeps every
node's connectivity. Strongest in delta, area under the curve point-eight-three, then low-gamma, then beta.
And notice the tissue split: beta spares the seizure zone — its cognitive trace lives on healthy cortex —
while low-gamma is the seizure zone. Sweeping tau reorders the candidate list but doesn't sharpen the top;
precision comes from fusing the bands — sixty percent in the top five, about seven times chance. There are
really two populations: eight patients are this diffusion community, up to point-nine-nine at the community
level, and two are right-hemisphere hubs where plain strength does the work. Calibrated, leave-one-patient
-out, the detector sits at a median area-under-the-curve of point-eight-seven, nine of ten above chance.
One operator, one null — a cognitive lens that doubles as a clinical one.

Careful: scope is clinically-labelled SOZ, NOT surgical outcome — triage marker, not diagnostic. The SOZ is
a strength-INDEPENDENT co-diffusing community (relational), never a per-contact classifier. δ is an epilepsy
band with NO cognitive trace — do NOT call it cognitive (that IS the dissociation, not a contradiction). β
SPARES the SOZ (confirmed 3 ways: pair-class, node-exclusion, size-matched decimation); α does NOT clearly recruit the SOZ on the backbone — do NOT
lead with the old "α recruits +0.41" (whole-graph, softened). τ = ranking, NOT precision; precision is band
FUSION — don't credit τ with the hit-rate. Null = matched-strength ONLY (no drift, no second null); never
rank bands by the ×-null ratio; read AUCs as reported (per-band), not best-scale. low-γ→SOZ may be an
anchor not a trace (matched-strength controls strength, not tissue stability) — fine to present as the disease marker regardless. β is the bridge
because it is in BOTH families — coherent, not a contradiction. Cuttable if you run long (keep the band
dissociation + the AUC line; drop the two-populations/LOPO detail).

5. FIGURES
- **MAIN — SOZ = strength-independent diffusion community.** `data/outputs/figures/talk/fig_soz_diffusion_community.pdf`
  (relational-community layout — still valid).
- **Per-band AUC vs fake-SOZ null.** `data/reports/results_section3/fig_epi_a_relational_marker.pdf`
  ⚠ REBUILD to the mst@0.20 AUCs (δ .83 / low-γ .82 / β .745) — the script/panel still carry the dense
  δ .80 / low-γ .74 / β .69.
- **Two populations (community vs right-hemi hub).** `data/reports/results_section3/fig_epi_c_two_populations.pdf`
  (verify community AUC↑.99 + the 2 hub cases .49/.57).
- **Calibrated LOPO detector.** `data/reports/results_section3/fig_epi_b_calibrated_detector.pdf` (verify median .87, 9/10 > chance).
- **Scope / ceiling (τ=ranking; fused prec@5 = 60%).** `data/reports/results_section3/fig_epi_d_scope_ceiling.pdf`.
- **Band-dissociation table (build/relabel):** two columns (cognitive trace · epileptogenic AUC) across
  θ/α/β/low-γ/δ/high-γ, β flagged in BOTH columns as the bridge; NO OFC cell for β. PNG fine for the talk.
- ⚠ RETIRE the "α recruits the SOZ" panel framing; if reused, relabel to "β spares the SOZ."

6. REFERENCES
- Ours (audit_132 all-contacts propagator SOZ marker; audit_89 diffusion community; mst@0.20 recovery arc 2026-07-12; low-γ→SOZ verdict 2026-07-13).

7. CANVA STATUS
NEW slide (merges old 19 band-dissociation + old 20 epilepsy). Band dissociation (cognition α/β · epilepsy
δ/low-γ/β · β bridge) + the strength-independent SOZ community marker (AUC δ.83/low-γ.82/β.745; fake-SOZ null
δ7/low-γ8/β8), β spares the SOZ, α-recruits SOFTENED, τ=ranking not precision (fused prec@5 60%), two
populations (community↑.99 + 2 right-hemi hubs), LOPO detector median .87, dense→mst@0.20 robustness. Missing
on deck: rebuild fig_epi_a to the mst@0.20 AUCs, place the section-3 figures + band-dissociation table,
compressed on-slide text, presenter notes. Cuttable if running long.
