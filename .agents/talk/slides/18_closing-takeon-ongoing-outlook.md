---
name: talk-slide-18-closing-takeon-ongoing-outlook
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery, ρ_sym)
slide: 18
status: draft
updated: 2026-07-16
canva: NEW compound CLOSER — absorbs old 17 (take-homes) + old 18 (outlook-thanks) into ONE
  three-part slide: TAKE-HOMES · ONGOING · OUTLOOK (+ thanks footer). Each take-home closes a
  slide-06 "why it matters" anchor ("oh — they solved that"); the AI/outlook point closes the arc
  back to the slide-06 opening. β localization is OPEN (ongoing), never a result (delocalized 2026-07-14).
---

# Slide 18 — Take-homes · ongoing · outlook  (the compound closer)

1. TITLE
Take-homes, what's ongoing, and where it goes

(Alt punchier title if wanted: "From a change we can detect to structure we can read — and what's next".)

2. MAIN CONCEPT
- **Head (the juice):** one sparsified diffusion operator, one null (matched-strength), reads the brain
  across its three axes at once — and every "why it matters" problem we opened with, pinned in the
  literature only by *behaviour and theory*, we can now **read quantitatively in the network**. Cognitive
  maps → a multiscale lens; offline consolidation → encoding vs inference held at rest; replay → a held,
  sustained trace; data-efficiency → the reason the human solution is worth reading out. Detecting a
  change became **characterizing** it, and the same lens doubles as a clinical one. Three parts close the
  talk: what we take home, what we're working on now, and where it goes — ending back at the opening
  (few-shot human cognition → inspiration for data-efficient machines).
- **The compound closer is three columns:** TAKE-HOMES (the messages, ~6), ONGOING (open threads on our
  bench, ~6), OUTLOOK (where it goes, ~6) — plus a thin thanks footer. Deliberately messages, not the
  detailed result tables (those were delivered on 13 / 15 / 16); the numbers live in presenter notes.

### TAKE-HOMES — each one closes a slide-06 "why it matters" anchor
(The arrow ↩ names the slide-06 problem it answers; keep the arrows in presenter notes, not on-slide.)
- **T1 · One lens, three axes.** A single multiscale diffusion framework reads the brain across time
  (frequency bands), space (high-resolution sEEG), and topology (the diffusion hierarchy) at once.
  ↩ *cognitive maps / multiscale* (Behrens 2018, Son 2021).
- **T2 · From phenomenology to structure — detect became characterize.** Simple pairwise connectivity
  *detects* that something changed; the multiscale hierarchy *characterizes* what kind — **which band**
  (α/β), **at which scale** (β scale-invariant across all 16, α mesoscale), and a **higher-order coarse
  grouping the edges miss**. Many relevant features hide in this higher-order embedding.
  ↩ *the gap: read the change quantitatively, not just by behaviour.*
- **T3 · A held, offline trace of reasoning.** The task's reorganization outlives the task — a **sustained
  resting state** (reinstatement 10/10), an offline *consolidator* of the learned structure, not a
  transient replay burst. ↩ *replay / offline reactivation* (Liu–Behrens 2019).
- **T4 · Rest keeps the order it inferred, not just the pairs it saw.** Encoding and inference both
  persist; the inferred order is **mesoscale-emergent** — visible only as the network coarse-grains, the
  shape of an *abstraction*, not a memorized detail. ↩ *offline consolidation / encoding-vs-inference*
  (Ellenbogen 2007).
- **T5 · One operator, two readouts.** The *same* propagator that builds the cognitive hierarchy localizes
  the **seizure-onset zone** from one recording — a **triage-grade** marker (median AUC 0.85, 9/10 above
  chance). The bands dissociate cleanly: **cognition α/β · epilepsy δ/low-γ/β · β the bridge**, and the two
  reads sit on different tissue — **β spares the seizure zone, δ marks it**. ↩ *epilepsy is a multiscale
  network disorder* (Kramer & Cash 2012).
- **T6 · Only if it beats the boring explanation.** Every claim is discounted against a **matched-strength**
  null — the reorganization is not merely "where connectivity is strong". ↩ *the guardrail, repeated all talk.*

### ONGOING — open threads currently on the bench
- **O1 · Per-patient trace phenotyping in the heterogeneous bands.** δ / θ / γ carry **strong traces per
  patient that don't yet cohere** across the cohort (γ_low reaches Δ≈0.6 individually); mapping what drives
  the spread — laterality, patient-specific structure. (Fluctuations are signal, not noise.)
- **O2 · Higher-order, beyond-pairwise markers.** Pushing from second-order |ImCoh| edges into
  hypergraph / higher-order structure — both to sharpen the cognitive trace and for a stronger epilepsy
  marker.
- **O3 · Is low-γ→SOZ a task-trace or a stable-tissue anchor?** Matched-strength controls connection
  strength, not tissue *stability* — settling whether the disease-band localization is task-induced or a
  fixed property of that tissue.
- **O4 · What exactly do simpler metrics miss?** Pinning the multiscale value-add over raw connectivity and
  simple descriptors, band by band (the raw-vs-cophenetic comparison) — so "characterize" is quantified,
  never over-claimed as "only multiscale detects it".
- **O5 · Reconciling the resting task-dynamics with the held trace.** How the windowed state-space
  trajectory through the encoding/inference space gives rise to the sustained cross-phase state.
- **O6 · Resolving β's spatial signature.** Genuinely **delocalized**, or a real *coarse* lateralization
  (a marginal left tilt) — on a powered, resolution-preserving cohort statistic. (Open — no anatomical
  home is claimed.)

### OUTLOOK — where it goes (nicest, arc-closing)
- **F1 · Close the loop to behaviour.** Acquire transitive-inference performance to link the *strength and
  shape* of the characterization to cognitive ability — the brain–behaviour bridge we deliberately don't
  claim yet.
- **F2 · Fluctuations as a biomarker.** Even without behaviour, a patient's large **departure from the
  cohort** — especially in bands where leave-one-out flips a band's encoding- or inference-specificity —
  may itself flag pathology or a cognitive deficit. The spread could be the diagnostic. *(Hypothesis, not
  a result.)*
- **F3 · From triage to clinic.** A prospective clinical test of the diffusion seizure-onset marker (and
  settling its trace-vs-anchor status) toward a deployable triage tool.
- **F4 · Richer, nonlinear decompositions.** Beyond the second-order edge and the linear diffusion read —
  nonlinear / higher-order decompositions resolving *specific timescales*, and deeper task-dynamics.
- **F5 · Bigger, more diverse data — and other tasks.** Larger, more varied cohorts, other cognitive
  paradigms, and a substrate-independent version of the multiscale trace: the lens is reusable wherever
  cross-scale signatures are expected.
- **F6 · Back to intelligence — and to machines.** Reading how the human brain builds rich relational
  structure from *little* data feeds straight back to where we started: mechanisms of few-shot human
  cognition as inspiration for **data-efficient AI architectures**. ↩ *closes the slide-06 opening
  (data-efficiency / the "AI" kick).*

3. ON-SLIDE TEXT   (condensed — single phrases; three labelled blocks + a thanks footer)

TAKE-HOMES
• one lens, three axes — time (bands) · space (sEEG) · topology (diffusion hierarchy)
• detect → characterize — which band (α/β) · which scale (β invariant, α mesoscale) · a higher-order grouping the edges miss
• a held, offline trace of reasoning — a sustained resting state, an offline consolidator (not a transient burst)
• rest keeps the order it inferred, not just the pairs it saw — inference mesoscale-emergent = the shape of an abstraction
• one operator, two readouts — same propagator localizes the seizure zone (triage-grade, AUC 0.85, 9/10); cognition α/β · epilepsy δ/low-γ/β · β the bridge; β spares the SOZ
• real only if it beats the boring explanation — every claim vs a matched-strength null

ONGOING
• per-patient traces in the other bands (δ/θ/γ) — strong but not yet cohort-coherent
• higher-order / beyond-pairwise markers — cognition and epilepsy
• is low-γ→SOZ a task-trace or a stable-tissue anchor?
• what exactly do simpler metrics miss — the multiscale value-add, band by band
• reconciling the resting task-dynamics with the held trace
• resolving β's spatial signature — delocalized, or a coarse lateralization?

OUTLOOK
• close the loop to behaviour — link characterization to inference performance
• fluctuations as a biomarker — a patient far from the cohort may be flagging pathology
• from triage to clinic — a prospective test of the seizure-onset marker
• richer, nonlinear decompositions — specific timescales, deeper task-dynamics
• bigger / more diverse data · other tasks · a substrate-independent trace
• back to intelligence — few-shot human cognition → inspiration for data-efficient AI

— thanks —  Sapienza collaborators (data + task) · the LRG authors (Villegas, Gabrielli, Poggialini, Gili)

4. SPEECH  (~100 s closing)
Let me pull it together. One diffusion operator, one null, reads the brain across its three axes at once — time in the frequency bands, space in the electrodes, topology in the hierarchy. And remember how we opened: cognitive maps, offline consolidation, replay, learning from little data — all of it established by behaviour and by theory, none of it read directly in the network. That's the move. We turned detecting a change into characterizing it: which band carries it — alpha and beta; at which scale — beta at every scale, alpha at the mesoscale; and a higher-order, coarse structure the edges simply miss. The trace is held — a sustained resting state, an offline consolidator, not a fleeting burst. And what rest keeps isn't just the pairs the brain saw, it's the order it inferred — a re-grouping that only appears when you zoom out, the shape of an abstraction. Then the same operator pays a clinical dividend: read inside one recording it localizes the seizure zone, triage-grade, and the bands split cleanly — cognition in alpha and beta, epilepsy in delta, low-gamma and beta, on different tissue. And none of it counts unless it beats the boring explanation — connection strength.

What are we working on now? Per-patient traces in the other bands, which are strong but don't yet line up across people; higher-order markers beyond pairwise edges; whether the low-gamma seizure signal is a trace or just stable tissue; exactly what the simpler metrics miss; and resolving whether beta has any spatial signature at all.

And where it goes — this is the part I'm most excited about. We want behaviour, to tie the characterization to how well someone actually reasons. But here's the twist: even without behaviour, a patient who sits far from the cohort — especially where dropping them flips a band's signature — might be telling us about their own pathology. The fluctuation could be the diagnostic. Then richer, nonlinear decompositions; bigger and more varied data; a version of the trace that doesn't depend on this substrate. And finally, full circle: understanding how the human brain builds so much structure from so little data is exactly the kind of mechanism that could inspire more data-efficient machines. Thank you — to our Sapienza collaborators for the data and the task, and to the LRG authors whose framework made all of this possible.

Careful (guardrails — this is the CLOSER; every claim must match slides 14/16/17 exactly):
- **Matched-strength is the SOLE null.** Never "drift", "second null", "survives both nulls". Every take-home
  is discounted against connection strength — T6 states it explicitly.
- **detect → CHARACTERIZE, never "only multiscale sees it".** Raw pairwise FC DOES detect the β change
  (robust, complementary); the edge is discrimination (α/β) + scale (β invariant / α mesoscale) + a
  higher-order grouping. O4 is exactly the honest "quantify the value-add" thread — don't pre-empt it.
- **NO β anatomical home. β is DELOCALIZED** (localization slide dropped 2026-07-14). β localization is
  ONGOING (O6), framed as an OPEN question, NEVER a result. No β→OFC, no "coarse left home" as a finding.
  The only focal address anywhere is the DISEASE one (low-γ→SOZ), and even that may be an anchor (O3).
- **Held / sustained, NOT transient replay.** T3 = sustained reinstatement (10/10); replay is the slide-06
  *motivation concept* (Liu–Behrens), our version is the held state — "offline consolidator", not "we
  replicate replay".
- **Inference is α AND β** (7/16 each, δ cut), mesoscale-emergent. Never "β-only inference".
- **Trace = cohort CONSISTENCY, not presence/absence.** The other bands are NOT "trace-free" — they trace
  per patient but don't cohere (that's O1). Never call any band silent.
- **Epilepsy: triage marker, NOT diagnostic; SOZ labels, not surgical outcome.** β SPARES the SOZ, δ MARKS
  it — a dissociation, not a contradiction. Be honest it's 9/10 (Pat_15, the right-hemi hub, at chance —
  the same cognitive outlier). Lead the marker as "strongest in δ / low-γ" (β is the strongest *cognition*
  band, only a modest marker at 0.75).
- **NO brain–behaviour correlation is claimed** — we have no per-patient behaviour. F1 is the OUTLOOK to
  get it; F2 (fluctuation-as-biomarker) is an explicit HYPOTHESIS, flag it as speculative, never a result.
- **The AI point (F6) is soft** — few-shot human cognition as *inspiration* for data-efficient AI. Never
  "we beat AI", "we build AI", or "AI has solved this". It closes the slide-06 data-efficiency anchor.
- Do NOT re-teach the three axes or the Betzel schematic (slide 5 owns them); T1 names them, doesn't
  re-derive.

5. FIGURES
Text-forward, three columns. At most ONE compact synthesis glyph — keep the words dominant.
- OPTIONAL synthesis glyph — the band-dissociation one-liner as a tiny table/strip: **cognition α/β ·
  epilepsy δ/low-γ/β · β the bridge** (a compressed copy of the slide-16 dissociation). PNG fine. This is
  the single cleanest visual for T5 if a visual is wanted.
- OPTIONAL "one lens · three axes" 3-icon glyph (bands / sEEG / diffusion) built in Canva for T1.
- NO new result plots here — the results were delivered on 13 (trace), 15 (encoding/inference), 16
  (epilepsy). The closer is messages.
- ⚠ RETIRE any "which method wins" / whole-graph Grassmann dissociation panel / β→OFC anatomy panel
  (superseded framings).

6. REFERENCES
On-slide names carried from earlier slides (no new citations needed on the closer):
- Cognitive map / task: Behrens 2018 · Son 2021 (slide 6). Replay / consolidation: Liu–Behrens 2019 ·
  Ellenbogen 2007 (slide 6). Epilepsy: Kramer & Cash 2012 (slide 6). Method: Villegas 2023/2025 (LRG) ·
  Nolte 2004 (imaginary coherence).
- All results are ours (mst@0.20 backbone; matched-strength null; ρ_sym cophenetic trace; diffusion SOZ
  marker). No new external citation.

7. CANVA STATUS
NEW compound CLOSER — absorbs old slide 17 (take-homes) + old 18 (outlook-thanks) into ONE three-part
slide. Old 17 / 18 are SUPERSEDED (kept on disk for their detailed numbers + speech, which feed the
presenter notes here). Deck now ends on this single slide.

Structure: THREE labelled blocks (TAKE-HOMES ~6 · ONGOING ~6 · OUTLOOK ~6) + a thin thanks footer.
Each take-home closes a slide-06 "why it matters" anchor (one lens↩cognitive maps · detect→characterize↩the
gap · held trace↩replay · encoding/inference↩offline consolidation · one operator↩epilepsy · matched-strength↩guardrail);
the final outlook point (F6, AI inspiration) closes the arc back to the slide-06 opening (data-efficiency).

Settled numbers behind the message labels (presenter notes only, not on-slide): β trace 16/16 scale-invariant,
α 12/16 mesoscale; reinstatement 10/10 sustained (~88% co-move); inference α AND β 7/16 each, mesoscale-emergent;
epilepsy marker δ .83 / low-γ .82 / β .75, detector median AUC 0.85, 9/10 above chance, Pat_15 at chance.
Matched-strength is the sole null. β delocalized (localization OPEN, ongoing — never a result).

Missing on deck: build the Canva page (three text columns + optional dissociation glyph + thanks footer),
paste presenter notes. Thanks can be spun to its own final slide if preferred (keep the acknowledgements
either way).
