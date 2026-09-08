---
name: talk-slide-17-one-operator-epilepsy
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 ONLY — one backbone for trace AND marker)
slide: 17
status: draft
updated: 2026-07-16
canva: REBUILT on mst@0.20 ONLY (TMFG dropped 2026-07-16 — one backbone, nothing to justify).
  Marker AUCs → mst@0.20 (δ.83/low-γ.82/β.75). Cohort-recovery / rescue / two-backbone all REMOVED
  (they were TMFG-specific). Honest: detector median 0.85, 9/10 above chance, ONE right-hemi hub
  (Pat_15) at chance. Rescue panel deleted; slide = 3D hero + 3 curves + bullet list.
---

# Slide 17 — From Cognition to Pathology

1. TITLE
From Cognition to Pathology
SUBTITLE (the thesis, on-slide small or spoken): one operator, two readouts — the same diffusion propagator that builds the cognitive hierarchy also localizes the seizure-onset zone.

2. MAIN CONCEPT
- **Head (the juice):** the *same* diffusion propagator — the heat kernel \(e^{-\tau L}\) on the *same*
  mst@0.20 connectivity backbone that carries the cognitive trace — has two readouts. Read **across
  phases**, how much its task reorganization persists is the cognitive trace. Read **within a single
  recording** as diffusion affinity from a seed, it localizes epileptogenic tissue. Nothing retuned —
  same operator, same backbone, read a different way.
- **The seizure-zone marker is relational, not a per-contact label.** A contact-by-contact SOZ
  marker does not exist here: within a patient any contact-level signal collapses into how strongly
  the contact is coupled, and no contact-level rule transfers across patients. Read *relationally*,
  the seizure-onset contacts **diffuse together** — a co-diffusing community. With node strength
  removed, a contact's seed-affinity to the rest of the onset set separates seizure from healthy
  tissue **above a matched-strength fake-SOZ null**: strongest in **δ (AUC 0.83)** and
  **low-γ (0.82)**, then **β (0.75)** — clearing its own null in **low-γ 8/10 · β 8/10 · δ 7/10**
  patients (α and θ weak, 0.61 / 0.74, 5/10 — not disease-band markers).
- **The protocol (minimal).** Put heat on the known onset seeds and read the propagator: contact
  \(i\)'s score is its seed-affinity \(a_i = \langle e^{-\tau L}\rangle_{\mathrm{seed}\to i}\) with node
  strength regressed out (so it is *relational*, not "where coupling is strong"); rank contacts → the
  SOZ shortlist. Deploy either **one champion band (δ, AUC 0.81)** or **fuse all six** into the
  leave-one-patient-out logistic (median 0.85).
- **The role of τ.** τ is the diffusion *time* — how far the heat spreads from the seeds. The marker
  reads a **coarse** scale (\(s=\tau\lambda_{\max}\approx 10\)): heat travels several steps, so it
  reaches the onset contacts on *other* electrodes (the distant discovery / the arcs), not just the
  seeds' immediate neighbours. Fine τ keeps it local; coarse τ is what makes the seizure zone read as
  one distributed co-diffusion group. (A single fixed operating scale, same for every patient — not
  tuned per case.)
- **Distant discovery — proximity removed by construction.** Seed the known onset contacts on *some*
  electrodes, and the same affinity ranks the onset contacts on *different* electrodes above healthy
  ones (**δ AUC 0.72, 8/10**, a label-shuffle collapsing to chance). Because the \(|\mathrm{ImCoh}|\)
  substrate is immune to zero-lag volume conduction, this grouping cannot be spatial proximity.
- **A calibrated detector.** Fuse the six bands into a leave-one-patient-out logistic probability of
  seizure onset: **median AUC 0.85**, **9 of 10 patients above chance**, probability **0.69 at true
  onset vs 0.33 at healthy** contacts, label-shuffle null **0.47**. The single strongest deployable
  band is **δ alone (AUC 0.81)**; fusion lifts top-five precision to **60%** (≈ 7× the base rate). A
  nonlinear model does not beat it (AUC 0.82) — we keep the interpretable logistic.
- **The honest exception.** One right-hemisphere implant (**Pat_15**) sits at chance (**AUC 0.50**):
  its seizure zone is a network *hub* rather than a co-diffusing community, so plain node strength —
  not the propagator — marks its onset. It is the *same* patient that is the cognitive β-trace outlier
  (right-hemisphere-only implant): the one anti-aligned patient is anti-aligned in both reads —
  coherent, not two separate failures.
- **The cognitive and disease reads sit on different tissue.** **β spares the SOZ** — its cognitive
  trace lives on healthy cortex (three ways: pair-class, node-exclusion, size-matched decimation).
  And **δ is the cleanest split**: it marks the seizure zone strongly but carries **no** cognitive
  trace at all.
- **Scope:** this reads clinically-labelled **SOZ, not surgical outcome** — a **triage marker, not a
  diagnostic** — and it spans tissue classes (**41% of onset contacts lie in white matter**, ranked
  without a gray-matter prior).

3. ON-SLIDE TEXT — FINAL (user-edited 2026-07-16; this is the ONLY text on the Canva slide)

Seed the known onset →  aᵢ = ⟨ e^{−τL} ⟩_{seed→i}  (strength removed)  → rank
τ = diffusion time: coarse (s≈10) → heat reaches distant onset, not just neighbours

one champion band (δ) — or fuse 6 → detector
Seizure zone = a co-diffusing community  (onset set = a high-affinity block of K)
Calibrated detector: median AUC 0.85  — 9/10 above chance  [ROC]
Triage marker, not a diagnostic  — top-5 ≈ 7× base rate  [precision@k]

⚠ NOTE: the user's draft parenthetical "(sum over rows of Propagator K)" was CORRECTED — the full
row-sum Σⱼ K[i,j] is the heat-centrality ≈ the node-strength baseline that is REGRESSED OUT
(contradicts "strength removed"). The marker is the PARTIAL sum over the ONSET columns only,
aᵢ = Σ_{j∈onset} K[i,j] (line 1). Use "(onset set = a high-affinity block of K)" — the structural
co-diffusing-community meaning — or "(sum the ONSET rows of K)", never the bare full row-sum.

4. SPEECH
And the same tool pays a second dividend — a clinical one. Same propagator, same backbone that built the cognitive hierarchy — only now, instead of asking how the hierarchy persists across phases, we drop heat on a few known seizure-onset contacts in one recording and watch where it spreads. Score each contact by the seed-heat it receives, node strength removed so it's relational, read at a coarse scale so the heat reaches the onset contacts on *other* electrodes, not just neighbours — the distant discovery you see on the brain. The seizure zone lights up as its own co-diffusing community, above a null that scrambles the zone but keeps each node's connectivity — strongest in delta and low-gamma, AUC around point-eight.

Deploy delta alone, or fuse all six bands into a leave-one-patient-out detector, and it lands at median AUC point-eight-five, nine of ten patients above chance — a triage marker, not a diagnostic, its top-five shortlist about seven times base rate. The one honest exception, a single right-hemisphere implant at chance, is the *same* patient that's our cognitive outlier — a hub, not a community. And the two reads sit on different tissue: beta spares the seizure zone, delta marks it and carries no cognitive trace at all. One operator, one backbone — a cognitive lens that doubles as a clinical one.

Careful: scope is clinically-labelled SOZ, NOT surgical outcome — triage marker, not diagnostic. The
SOZ is a strength-INDEPENDENT co-diffusing community (relational), never a per-contact classifier. δ
is an epilepsy band with NO cognitive trace — do NOT call it cognitive (that IS the dissociation, not
a contradiction). β SPARES the SOZ (confirmed 3 ways). ONE backbone now: trace AND marker both on
mst@0.20 (TMFG dropped 2026-07-16) — if asked "same graph as the trace?" the answer is YES, literally
the same backbone and operator, read within-phase instead of cross-phase (nothing to justify). Lead
the marker as "strongest in δ / low-γ" — β is the strongest *cognition* band but only a MODEST marker
(0.745); α/θ are weak markers (0.61 / 0.74, 5/10). Be honest that it is 9/10, NOT cohort-wide: Pat_15
sits at chance (state it — the ROC shows it hugging the diagonal). Null = matched-strength ONLY.
low-γ→SOZ may be an anchor not a trace (matched-strength controls strength, not tissue stability) —
fine as the disease marker regardless. Cuttable if you run long (keep the marker AUC line + the
detector line; drop the distant-discovery / δ-only detail). NOTE: cohort-recovery / rescue / "all 10
detect" were TMFG-specific and are RETIRED — do NOT say them on mst@0.20.
prec@5 multiplier: keep ONE on the deck. On-slide "~7×" is median-based (prec@5 ≈ 0.60 / base ≈ 0.086 → ~7×); the §5 figures note's "~5×" is mean-based (0.54 vs a random ~0.10). Speech says "about seven times" to match the on-slide; reconcile §5 if you standardize.

5. FIGURES  (dark/transparent talk PNGs, rebuilt 2026-07-16 on mst@0.20; scripts in scripts/07_figures/)
Slide set = ONE visual + THREE curve panels + the on-slide BULLET LIST (§3). All words live on the
slide, NOT in the figures (only axis labels). Every curve panel = per-patient fluctuations (thin)
under a bold cohort average; all on the mst@0.20 backbone against the same matched-strength null.
- **VISUAL — SOZ = a co-diffusing community (3D brain, zoomed to the contacts).**
  `data/outputs/figures/talk/fig_epi_seed_spread_3d.png` — exemplar Pat_14 δ rest_post: heat placed
  on the known onset contacts of one electrode (gold diamonds) diffuses across the
  volume-conduction-immune \(|\mathrm{ImCoh}|\) graph and re-reaches the *unseeded* onset contacts on
  *other* electrodes — **solid bright-green discs = correctly-identified onset** (the seizure zone
  rediscovered; disc size ∝ affinity), warm arcs threading the co-diffusing onset web; healthy cortex
  stays cool. A few **steel ✕ = false positives**. Solid green not open rings — plotly's 3D renderer
  draws thin strokes almost invisibly, so the onset marker must be a filled glyph. Zoomed to the
  electrode bbox; out-of-pial contacts snapped just inside the schematic shell (cosmetic — coords are
  display-only). Script `talk_fig_epi_seed_spread_3d.py` (glyph legend GREEN=hit / ✕=false-alarm /
  gold=seed; conceptual illustration on the dense graph — distant-discovery δ 0.72, 8/10 is robust).
- **CURVE 1 — seizure-zone AUC across bands.** `data/outputs/figures/talk/fig_epi_auc_by_band.png`
  — per-patient AUC profile (thin) + bold cohort-median + IQR; band-coloured nodes; chance at 0.5.
  Reads: peaks δ .83 / low-γ .82, β .75, dips α (.61). (source: `epi_arc_mst020`.)
- **CURVE 2 — detector ROC (leave-one-patient-out).** `data/outputs/figures/talk/fig_epi_detector_roc.png`
  — thin per-patient node ROC + bold mean ROC (AUC ≈ 0.85); 9 patients bow above the diagonal, one
  right-hemi hub (Pat_15) hugs it — the honest exception, shown not stated. (source: `marker_detector_mst020`.)
- **CURVE 3 — precision @ k (vs a random-ranking band).** `data/outputs/figures/talk/fig_epi_precision_at_k.png`
  — thin per-patient + bold cohort-mean precision vs shortlist size, over a shaded **random band**
  (5–95% of cohort-mean precision under random ranking, Monte-Carlo; narrows toward the base rate).
  Reads: the detector's top-list sits clearly ABOVE random — prec@5 = 0.54 mean / 0.60 median vs a
  random 5–95% of [0.04, 0.18] (≈ 5× the base rate). Curves from `talk_fig_epi_curves.py`.

Retired / do NOT place: fig_epi_rescue (cohort-recovery — TMFG-specific, deleted), fig_epi_band_dissociation,
fig_epi_two_populations, fig_epi_relational_marker (in-figure text / dense dot-matrices).

6. REFERENCES
- Ours (audit_132 all-contacts propagator SOZ marker; audit_89 diffusion community; detector on
  mst@0.20 `27_marker_detector_backbone.py` SA_BACKBONE=mst020; low-γ→SOZ verdict 2026-07-13). NOTE:
  the 2026-07-15 sparsifier verdict recommended TMFG for the marker — SUPERSEDED 2026-07-16 by the
  one-backbone (mst@0.20-only) decision.

7. CANVA STATUS
REBUILT on mst@0.20 ONLY (TMFG dropped 2026-07-16 — one backbone). Marker AUCs → δ.83/low-γ.82/β.75
(beats-null low-γ8/β8/δ7; α/θ weak 5); calibrated detector median 0.85, 9/10 above chance, P(SOZ)
.69 vs .33, shuffle .47, δ-only .81, prec@5 60%. Honest exception: Pat_15 (right-hemi hub) at chance
0.50 — the same cognitive outlier. Cohort-recovery / rescue / two-backbone ALL REMOVED (TMFG-specific).
β spares the SOZ + δ = the clean split kept. Slide = 3D hero + 3 curves + 5-bullet list. Missing on
deck: place CURVE 1/2/3 (mst@0.20) + the bullets, presenter notes. Cuttable if running long.

⚠ DOWNSTREAM (not this slide): §3 (`results_sec_3.tex`) is currently written on TMFG — cohort
recovery, two-backbone, fig_epi2 rescue, β 0.87. If "mst@0.20 only" is project-wide, §3 needs the
same revert (bigger job). FLAGGED, not yet done.
