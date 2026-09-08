---
name: ellenbogen-relation-and-extension-2026-07-15
type: reference
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
status: current
created: 2026-07-15
scope: How Ellenbogen et al. 2007 (PNAS, "Human relational memory requires time and sleep")
  relates to our study, and precisely what we add and — honestly — what we give up. Written to
  back the slide-06 Ellenbogen anchor (`.agents/talk/slides/06_why-it-matters.md`) and reusable
  for the preprint related-work / discussion. Brutal-honesty posture: we do NOT "beat" their
  result; we open the box at a different level of analysis, with real trade-offs.
pointers:
  - .agents/talk/slides/06_why-it-matters.md            # the slide this backs
  - .agents/talk/2026-07-13_multiscale-narrative-audit.md   # settled results / do-not list
  - .agents/preprint/established_results/2026-07-13_settled-three-results.md
references:
  - Ellenbogen, Hu, Payne, Titone & Walker (2007) PNAS 104(18):7723–7728, DOI 10.1073/pnas.0700094104
---

# Ellenbogen 2007 ↔ our study: link, extension, and honest limits

## Head (the juice)

Ellenbogen measured the **behavioural output** — people can do the transitive inference
*later*, not right after learning — and from that inferred, indirectly, that the relational
order is bound "offline." We measure the **neural substrate** they could only hypothesize,
and we **characterize its form** (which rhythm, at what scale, distributed or focal). So we do
not "beat" their result: we open the black box behind it, at a different **level of analysis**
(neural vs behavioural), and with real trade-offs — no behaviour in our cohort, wake not sleep,
patients not healthy adults.

---

## What Ellenbogen 2007 actually showed (behaviour only)

- **Same task class.** 56 healthy young adults learned five premise pairs (A>B, B>C, C>D, D>E,
  E>F → hidden hierarchy A>B>C>D>E>F) to criterion with feedback, then were tested after an
  offline delay (20 min, 12 h wake, 12 h sleep, or 24 h) on the premises *and* on novel
  **inference pairs** — B>D and C>E (one degree of separation) and **B>E** (two degrees, the
  hardest).
- **The dissociation.** Premise-pair retention was near-identical across all delays (>85%) —
  the *building blocks* are equally held. But **inference was at chance right after learning**
  (20-min group 52%, not different from chance) and became **highly significant only after an
  offline delay** (12/24 h groups >75%, P<0.001).
- **Sleep boosts the hardest inference specifically.** For the two-degree pair B>E: **93% after
  sleep vs 69% awake** (P=0.03) — a ~35% relative advantage for the most distant judgment.
- **Below awareness.** The benefit was *not* accompanied by increased subjective confidence —
  it "operates below the level of conscious awareness."
- **Their interpretation.** Premise learning is "necessary but not sufficient"; the relational
  hierarchy is *bound offline* into an abstract **"metamemory representation,"** especially
  during sleep.

**Level of analysis = purely behavioural** (accuracy + confidence). The offline binding is a
**black box**: they observe the *output* (inference ability appears later) and posit an internal
reorganization they never measure. They *speculate* a hippocampal-replay mechanism; no neural
data.

---

## The link — what maps to what

| Ellenbogen (behaviour) | Ours (neural network, sEEG) |
|---|---|
| The inferred order is built **offline, post-task** | We probe exactly that window — **rest_post** (rest_pre → task_learn → task_test → rest_post) |
| Premises retained immediately; **inference develops later** | Neural **encoding vs inference-specific** decomposition — the inference component is the residual after partialling out the shown pairs (T_infspec_pe) |
| "Metamemory" = an **abstraction** of the hierarchy | **Mesoscale-emergent** trace: invisible at the finest diffusion scale, appears only after coarse-graining — an abstraction is *coarser* than the surface pairs |
| Benefit is **implicit** (no confidence gain) | A resting-state trace — an offline *state*, not deliberate recall |

Ellenbogen is our **external license**: it is *why* probing the resting brain after this task is
principled, and *why* we expect an offline, abstract object rather than a copy of the shown pairs.
It is motivation, not our data.

---

## What we add (the honest extension)

1. **Behavioural inference → a direct neural observable.** They *inferred* offline binding from
   later behaviour; we *see* the reorganization in the rest_post functional network.
2. **We separate encoding from inference at the neural level** and show the *inference-specific*
   component persists into rest (β 7/16 scales, p_meso 0.005; α 7/16, 0.005) — a neural counterpart
   to "the inferred order, not just the shown pairs." (Encoding fires already at the fine scale,
   β s₁ 0.032; inference is silent at s₁ (β 0.097) and emerges at the mesoscale.)
3. **We characterize the FORM, not just presence** — the core value-add (the slide's thesis).
   Ellenbogen has one number (accuracy). We can state **which band** (α/β), **at what scale**
   (β *scale-invariant*, 16/16, Friedman p 0.17 ns; α *single-scale/mesoscale*, 12/16, Friedman
   0.004), and that it is **distributed**, not focal. A behavioural score cannot express any of
   this. "Mesoscale-emergent" is, concretely, the neural signature of the abstraction they posit.
4. **Held as a state**, not merely inferable at a later test — sustained reinstatement (rest_post
   dwells in the task configuration more tightly than rest_pre, 10/10, p 0.001). ⚠ whole-graph
   τ_min companion, not a mst@0.20 backbone finding — present as corroboration.
5. **Reach they lacked** — intracranial access (they could only speculate about mechanism), and
   the *same* diffusion operator also reading the **epileptogenic network** (δ/low-γ/β), entirely
   outside their scope.

**In one line:** Ellenbogen turned "inference develops offline" into a behavioural fact; we turn
the *offline binding* itself into a measured, multiscale neural object.

---

## The onset-vs-maturation bridge (the talk's central hook)

> **Correction (2026-07-15).** An earlier draft of this section pinned "our rest_post window"
> onto Ellenbogen's **20-min test bar** and captioned it "behaviour still at chance → but our
> trace has already begun." That framing is **retracted**. It welded our recording *window*
> (a span, no % coordinate) onto their delayed *test result* (a value on a % axis) and — the
> fatal error — **hid that our `task_test` precedes `rest_post`**. In our design the patient
> has *already attempted the inference* (task_test) before the window we analyse; Ellenbogen's
> delay contains *no test*. So "our neural onset precedes behaviour" is **not** a claim we can
> make about our own cohort — they produced behaviour in task_test first. The corrected bridge
> below keeps the two studies on **separate lanes with no shared axis**.

The honest link is a **level-of-analysis** bridge over one **shared offline period**, not a
shared coordinate:

- Both paradigms concern the **offline period after the same relational task**, during which the
  inferred order is bound. Ellenbogen reads that period by **testing behaviour AFTER it** (a single
  delayed test; inference is at chance at 20 min and matures only over hours + sleep). We read the
  *same kind of* period by **recording the network INSIDE it** — `rest_post`, wake, minutes.
- The contrast is **onset vs maturation**: our window sits at the **onset** end (minutes, wake),
  where the reorganization is beginning; their **matured** inference lives at the far end (hours,
  sleep). Our lane makes **no claim past `rest_post`** — no 12 h / 24 h / sleep on our side.

Why this is the right hook — and its honesty rails:
- It **links to their result** (same task, same offline logic) and **foreshadows ours** (a network
  trace forming in rest_post) in one move — via a shared *period*, never a shared axis or value.
- It **does not touch behaviour** — no brain–behaviour correlation; "the trace itself is the
  evidence." The behavioural maturation is *theirs*; the onset trace is *ours*.
- We never say "neural precedes behaviour" **for our subjects** — they already produced behaviour in
  task_test. "Onset before behaviour" stays Ellenbogen's *between-subjects behavioural* fact.
- Do NOT claim we captured their *matured* inference or the sleep-dependent B>E effect: their endpoint
  needs hours + sleep; our window is early wake. We claim **onset**, not completion — "the
  reorganization has *begun* at rest," never "we see the inference they measured."
- The figure (`slide06_ellenbogen_offline_inference`, redrawn 2026-07-15 as *"Onset vs Maturation"*)
  bakes this in: **two separate rows** on a broken time-after-learning axis. TOP = their behaviour
  (% correct, blue/coral); BOTTOM = our protocol drawn as `rest_pre → task_learn → task_test →
  rest_post` with **task_test explicitly BEFORE** the teal read-window ("inference attempted here,
  before our window") and **no y-quantity** (network state, not a score). Their offline delay is
  drawn *empty* ("no test given"); our bottom-right is left *empty* ("beyond our recorded window —
  maturation not measured"). The only cross-lane element is the labelled relationship "they read the
  OUTPUT by testing AFTER it · we record the PROCESS beginning INSIDE it — no shared subjects, no
  shared axis, no brain–behaviour correlation." The design + red-team rationale is in the workflow
  synthesis that produced it (three independent lenses converged on *replace, do not salvage*).

---

## Where we do NOT improve on them (state this too)

- **No behaviour in our cohort.** Their entire design rests on a behavioural anchor; we have the
  substrate but not per-patient performance tied to it. Our "inference" reading is licensed by
  **task design + mesoscale-emergence + literature (them)**, *never* a brain–behaviour correlation
  (the project forbids claiming one). The "performance benchmark" is an **opportunity**, not
  something we have done. This is the largest asymmetry — and precisely why Ellenbogen is such
  valuable external license.
- **Wake rest, not sleep.** Their strongest effect *requires sleep* (12–24 h). Our rest_post is
  minutes of **wake** rest immediately post-task. We are likely catching the **early initiation**
  of offline consolidation, not the sleep-matured endpoint. Do NOT claim we captured "the same
  thing" — this is a genuine difference (and an interesting extension: the reorganization may begin
  immediately, in wake rest).
- **Patients, not healthy adults.** Ours are drug-resistant-epilepsy patients — diseased brain,
  no healthy baseline; the pathology confounds the cognitive claim (the Beat-3 point cuts both
  ways).
- **Correlational, offline.** Neither study manipulates the reorganization causally.

Net: the "improvement" is that we **open the black box** — from a behavioural inference of offline
binding to a direct neural measurement and a multiscale *characterization* of its form — at the
cost of their clean behavioural anchor and healthy-subject baseline, and in an early wake window
rather than across sleep. It is complementary, not a better version of the same measurement.

---

## Ready-to-use framings

- **Talk (slide 06 anchor):** "Ellenbogen showed *behaviourally* that this inferred order is built
  offline, at rest — but the reorganization itself was a black box. We record that window directly
  and **characterize** the trace: which rhythm, at what scale, distributed or focal — turning their
  inferred offline binding into a measured, multiscale neural object." (Q&A caveat kept ready: no
  behavioural link in our cohort; wake, not sleep.)
- **Preprint (related work / discussion):** lead with the shared paradigm and the encoding/inference
  dissociation as the behavioural precedent; position our contribution as the **neural substrate +
  scale-characterization** of the offline abstraction; declare the wake-vs-sleep and no-behaviour
  limits up front (per the results-not-story + brutal-honesty rules).

---

## Guardrails carried from the settled narrative (do not violate when citing this)

- Ellenbogen is **motivation, not our data**; no brain–behaviour correlation is claimed anywhere.
- Detect → **CHARACTERIZE**, never "only multiscale detects it" (raw FC detects the β trace).
- Inference is **α/β**, not β-only (δ inference-specific was cut as a partial-corr artifact).
- Single null: **matched-strength** (drift retired).
- No cognition **localization** ("distributed," network-scale — not a hotspot).
- "Held, not replayed" is a **whole-graph companion**, not a backbone result.
