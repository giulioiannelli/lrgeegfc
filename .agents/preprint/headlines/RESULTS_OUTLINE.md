---
name: results-outline
era: IMCOH_ABS_COHORT_N10
status: working_draft
kind: results-outline
scope: table of contents for the paper's RESULTS section — 3 result headlines + subheadlines, phrased as findings (methods are not DESCRIBED here — the framework, the two read-outs, and the surrogate controls are explained in Methods; method names in titles/text are fine). Titles set 2026-06-22 (noun-phrase headers; one alt each). R2 reframed to the flagship — offline abstraction of a learned structure — 2026-06-23.
updated: 2026-06-23
---

# Results — outline / table of contents

**How to use.** This is the **skeleton** of the Results section: **three result
headlines (R1–R3)**, each with sub-results. **Titles are set** — short
noun-phrase headers (detail lives in the sub-results), one alternative under each.
Rules held here: **don't describe the methodology in these sections** (results
are phrased as findings — the framework, the read-outs, and the controls are
explained in Methods, not here; method *names* in titles/text are fine); robustness facts appear as *properties of
the finding* ("not explained by connection strength", "stable across scales");
status flags in brackets where a claim is a *hint* or *in progress*. NN-fit
judgement at the foot.

---

## R1 — the persistence trace
**Title:** *A multiscale β-band afterimage of learning and reasoning.*  
*(alt: A β-band afterimage of the task experience.)*

- **R1.1 — The reorganization persists into rest, and only in some rhythms.** After the task, the resting network sits measurably closer to the task configuration than it did beforehand — the reorganization does **not** wash out. The persistence is **band-specific**: strongest and fully confirmed in **β**, while **θ shows nothing** — a selective signature, not generic drift.
- **R1.2 — A multiscale reorganization, not a reshuffle of the strongest links.** The change is invisible to a connection-by-connection comparison — **and to an off-the-shelf spectral-clustering / PCA analysis** (the α component is recovered only by the full multiscale read-out, missed by both the spectral embedding and the leading-mode subspace) — and is not explained by which links are strongest; it emerges only when the network is read across many scales at once, and it survives a control that holds each node's total connection strength fixed. *(β is recovered by the spectral view too — the method-superiority point rests on α + band-selectivity.)*
- **R1.3 — It concentrates in orbitofrontal cortex.** Though brain-wide, the trace is over-expressed — above each patient's own baseline — in **orbitofrontal cortex**, bilaterally, the hub the brain uses to build relational maps. A hotspot on a distributed trace, not a container *(implanted in 5 patients)*.
- **R1.4 — It rides healthy cortex, not the disease.** The persistence is carried by coupling between **healthy gray-matter cortex**; the epileptic seizure core does **not** carry it — so it is neither a pathology artifact nor an everywhere-and-nowhere effect.
- **R1.5 — *(bridge to R3)* Diseased tissue is treated oppositely by band.** The focal inference band (**β**) **spares** the epileptic core, while the diffuse memory band (**α**) **recruits** it — the first sign that the cognitive trace and the epileptogenic network are linked.

## R2 — *(flagship)* offline abstraction of a learned structure
**Title:** *Offline abstraction of a learned structure: the trace carries the relations the brain inferred, not just the pairs it saw.*  
*(alt: The brain keeps the relations it reasoned out.)*

- **R2.1 — The brain keeps the relations it reasoned out, not just the pairs it saw.** Splitting the persistent trace into what was **encoded** (premises shown) and what was **inferred** (novel pairs reasoned out, with encoding controlled out), the **inference-specific component itself persists** into post-task rest — in **β alone**. The offline trace carries the *computed* structure, not merely a replay of experience: an offline signature of abstracting a learned order. *(Solid, duration-controlled.)*
- **R2.2 — The consolidation lives at the scale of multi-step integration.** The inference-specific β component is **mesoscale-favouring** — it clears the strength-matched control most strongly at intermediate scales, the natural scale of the multi-step relational paths that transitive inference is built from, while every other band stays null there. *(Verified — mesoscale-robust, not exclusive.)*
- **R2.3 — The consolidation anchors in the orbitofrontal cognitive map, and learning sets the anchor.** The encoding component concentrates in **orbitofrontal cortex** — the relational-map hub — and the pure memorising phase, with no reasoning demanded, already leaves its own α/β trace there, so OFC is anchored by **both** task phases. *(Solid.)*
- **R2.4 — The two components live in different rhythms.** The **memory** component persists in **both α and β**; the **inference-specific** component persists in **β alone**. Not an artifact of the longer test session — the effect does not scale with recording length, and the only length-tracking trend rides α, not β. *(Solid, duration-controlled.)*
- **R2.5 — A hidden memory trace surfaces in the cingulate at low-γ.** A focal, encoding-related trace appears in the **cingulate in low-γ** that is completely invisible at the whole-brain level — the clearest case of the method reading a single region. *(Solid.)*
- **R2.6 — Inference points toward the cingulate.** The inference-specific component leans toward **cingulate cortex** rather than orbitofrontal cortex. *(Directional hint only — not an established location.)*
- **R2.7 — The cingulate seems to switch its content by rhythm.** Read together, the cingulate carries **memory at low-γ and inference at β** — one region multiplexing different cognitive content across frequencies. *(Suggestive — the inference half is a hint.)*

## R3 — the same organization identifies epileptogenic tissue
**Title:** *Propagator-inspired markers of epileptogenic tissue.*  
*(alt: Epileptogenic tissue read from the network propagator.)*

- **R3.1 — Seizure contacts form a co-organized network group.** Given only a few known seizure contacts as seeds, the same connectivity organization points to the others that belong with them — including **distant** contacts on entirely separate electrodes.
- **R3.2 — Combined across rhythms, it yields a seizure-zone probability.** Pooling the read-out across frequency bands produces a **calibrated, cross-patient probability** that any given contact belongs to the seizure-onset zone.
- **R3.3 — Doing it without seeds is the open frontier.** Whether the seizure group can be recovered with **no known seizure contacts at all** is **[under investigation]**.
- **R3.4 — Patients fall into two regimes.** Detection behaves one way when the seizure zone forms a tight community and differently when it sits as a **network hub** — two populations the method must treat separately.
- **R3.5 — Honest scope.** Everything here is validated against **clinical seizure-onset labels, not surgical outcome** — a triage/shortlist tool to focus clinical attention, not a standalone localizer.

---

## Title workshop (next pass)

*(R1–R3 titles set 2026-06-22 — short and expressive; one alt kept under each.)* Rationale (what each title leads on):
- R1 — "afterimage" of **both** learning (encoding) and reasoning (inference); multiscale + β; brain-wide, so OFC stays in R1.3 — out of the title.
- R2 — *(flagship)* leads on **offline abstraction of a learned structure** (the brain keeps the inferred relations); the per-band dissociation (encoding vs inference) + the mesoscale scale-signature + the OFC anchor are the *evidence*; the anatomy of inference itself stays out of the title (cingulate is a hint, OFC is a concentration).
- R3 — *propagator-inspired* markers; spans all bands (δ marker strongest); label-validated triage (R3.5).

## Notes held

- **Methods excluded by design** — the framework / two read-outs / surrogate
  controls / volume-conduction-immune substrate are a Methods section; they appear
  in Results only as result-properties (R1.2).
- **Replay narrative removed (2026-06-22).** The dynamic replay thread — the bursty
  test (former R4) and the "sustained reinstatement" framing once carried on R1 — is
  **pulled**: the bursty/transient-replay test came back **cohort-negative** and the
  "sustained replay" relabel of R1 was just the persistence trace under a riskier
  word. To be restored only if something new surfaces.
- **No behavioral data** — performance is unavailable (PI 2026-06-22); nothing in
  this outline may be tied to task performance.
- NN-fit (2026-06-23): NN core = **R1 (the method reveals a held, non-ergodic
  structure) → R2 (the flagship — what the held structure *is*: an offline
  abstraction of a learned order)**; R3 = clinical coda. The
  "abstraction/consolidation" framing rests on the phase decomposition + the
  mesoscale scale-signature + the OFC anchor + literature, **not** on replay or
  behavior. The no-behaviour ceiling is the binding constraint — the abstraction
  claim is an *interpretation* of an inference-specific persistence result, never a
  behavioral or mechanistic proof. All remaining claims sit on proven results.
