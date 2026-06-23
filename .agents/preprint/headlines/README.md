---
name: headlines-readme
era: IMCOH_ABS_COHORT_N10
status: current
kind: headline
scope: master blueprint for the paper's headline architecture — ONE methodology core + 3-4 neurophysiological headlines that emerge from it; routing, file spec, global rules. This chat crystallizes headlines; per-headline agents do the investigation. N2 reframed to the flagship (offline abstraction of a learned structure), kept in slot 2, 2026-06-23.
updated: 2026-06-23
---

# `headlines/` — the paper's headline architecture (blueprint)

**Head.** The paper is **one methodology core** — a multiscale Laplacian
diffusion read-out of sEEG functional connectivity, ρ̂(τ) = e^{−τL̂}/Z, with two
non-redundant probes (per-pair cophenetic `ρ^coph` + global-mode Grassmann `d_G`)
and a locked C1–C5 null battery — **from which 3 neurophysiological headlines
emerge**: (1) a band-specific, multiscale **persistence trace** left by a
transitive-inference task — the method reveals a *held, non-ergodic structure* in
post-task rest (β, orbitofrontal, gray-matter cortical coupling); (2) **the flagship
— what that held structure *is*: an offline abstraction of a learned order** (the
brain keeps the relations it *inferred*, not just the pairs it saw; evidenced by an
inference-specific component that persists in β, favours the mesoscale of multi-step
integration, and anchors in the orbitofrontal cognitive map — this is where
`task_learn` finally earns its place); (3) the **epileptogenic read-out** of the same
propagator (seed-based, now seed-free, SOZ markers).
"One operator, multiple read-outs." This folder holds one crystallized file per
headline, plus the core; each is a scaffold a dedicated agent fleshes out.

> **This chat's job (and the only job here):** crystallize headlines into
> paper-ready scaffolds and write the routing/instructions below. **Investigation,
> validation, and number-crunching happen in the per-headline agent chats**, not
> here. When a headline needs work, write the task into that headline file's
> *To-dos & verifiables* section and hand the file to the owning agent.

---

## 1. Architecture — core + headlines (the amalgam)

```
CORE  ── multiscale Laplacian framework  [00_CORE_methodology.md]
        · ρ̂(τ)=e^{−τL̂}/Z, one operator, scale-continuous renormalization
        · two probes: ρ^coph (per-pair, multiscale, interpretable, localizable)
                       d_G(k)  (global leading-mode subspace ≈ PCA/spectral-clustering)
        · the multiscale-beats-global-spectral argument (probe dissociation)
        · band-resolution-not-amplification (3-layer: raw FC → raw D(τ) → cophenet)
        · C1–C5 null battery; matched-strength is the referee
        · |ImCoh| substrate (volume-conduction immune); τ-robustness
            │
            ├── N1  THE PERSISTENCE TRACE                    [01_trace.md]
            │     the task leaves a band-specific multiscale trace in rest_post;
            │     β both probes; over-expressed in OFC; gray-matter cortical;
            │     emergent (not raw-strength); invisible to edge-wise comparison
            │       N1.1 trace exists & is band-specific (β both / α cophenet-only)
            │       N1.2 emergent multiscale property, survives matched-strength
            │       N1.3 concentrates in orbitofrontal cortex (localization)
            │       N1.4 carried by healthy gray-matter coupling, not the epi core
            │       N1.5 τ-robust, fine-scale (audit_121)
            │       N1.6 BRIDGE→clinic: β spares / α recruits the epileptic core
            │
            ├── N2  (FLAGSHIP) OFFLINE ABSTRACTION OF A LEARNED STRUCTURE   [02_encoding_vs_inference.md]
            │     what the held structure IS: the brain keeps the relations it
            │     *inferred*, not just the pairs it saw. Evidence — four phases
            │     (rest→LEARN→TEST→rest), task_learn = encoding ref; α keeps encoding
            │     only, β keeps encoding + inference; the inference-specific component
            │     persists (β), favours the mesoscale of multi-step integration, and
            │     anchors in the OFC cognitive map
            │       N2.1 the four-phase decomposition (task_learn's role)
            │       N2.2 inference-specific component persists (β-only, mesoscale-favouring); encoding echo in α & β
            │       N2.3 content localizes differently (encoding→OFC; inference→cingulate*)
            │       N2.4 low-γ focal cingulate encoding trace (hidden by whole-brain avg)
            │       N2.5 duration confound RESOLVED — β arc duration-robust (clean
            │            ratio regression + matched-strength); naive truncation null
            │            retired as invalid (δ neg-control fails); α duration-suspect
            │       N2.6 the learning phase leaves its own OFC trace (both phases anchor OFC)
            │     * inference→cingulate DOWNGRADED by the duration control (directional hint)
            │
            └── N3  EPILEPTOGENIC READ-OUT                    [03_epileptogenic_markers.md]
                  the same propagator carries SOZ information
                    N3.1 relational δ diffusion community marks distant SOZ
                    N3.2 seed-based 6-band compound detector (calibrated P(SOZ))
                    N3.3 seed-FREE markers (live: audit_118/119/120)
                    N3.4 the two-population split (community vs hub patients)
                    N3.5 honest scope (clinical-label not outcome; precision ceiling)
```

**Headline count: 3 (replay narrative removed 2026-06-22).** The dynamic replay
thread — both the bursty test (former N4) and the "sustained reinstatement" framing
once carried on N1 — has been **pulled from the active structure**: the
bursty/transient-replay test came back **cohort-negative** (audit_124–131) and the
"sustained replay" relabeling of N1 was not earning its keep (it was the persistence
trace under a riskier word). The former N4 headline is **archived**
(`archive/2026-06/`), to be restored only if something new surfaces.
**N2 stays unified and is the paper's flagship** (the stage-2-arc duration control
resolved — β inference arc is duration-robust; the inference *component* is
anatomically thin and has no behavioral anchor, so it is N2's climax not a
standalone headline — but the **abstraction** framing built on it is the
headline-of-headlines, kept in slot 2 after N1, before N3). Core = **N1 (the method
reveals a held structure) + N2 (the flagship — offline abstraction of a learned
structure) + N3 (epilepsy)**. (The abstraction/consolidation framing rests on the
phase decomposition + the mesoscale scale-signature + the OFC anchor + literature,
not on replay or behavior.)

**Where did the "spectral critique" go?** It is **not** a standalone
neurophysiological headline — it is the **methodological core** (`00_CORE`).
The point that global eigenmode / PCA / spectral-clustering descriptors miss
multiscale per-pair reorganization *is* the justification for the whole method,
so it threads through N1–N3 and is argued in full in the core. (Earlier drafts
floated it as "Headline IV"; folded into the core 2026-06-22 per PI.)

---

## 2. The narrative spine (how to sell it)

1. **Method** (core): connectivity is usually read with a single threshold/scale
   or with global spectral embeddings that are uninterpretable and scale-blind.
   We read the *whole* diffusion hierarchy. Two probes; matched-strength referee.
2. **It reveals a held, non-ergodic structure** (N1): the resting network does not
   return to where it started — a multiscale persistence trace, band-specific (β),
   spatially concentrated (OFC), in healthy cortex, that the brain *holds as a
   state*. An emergent hierarchical property edge-wise/global-spectral tools miss,
   not a strength reshuffle. *(What the method can see.)*
3. **And what it holds is an abstraction** (N2 — the flagship): the held structure
   is the relational order the brain *inferred*, not just the pairs it saw. An
   inference-specific component persists (β, duration-robust), favours the mesoscale
   of multi-step integration, and anchors in the orbitofrontal cognitive map — **an
   offline signature of abstracting a learned structure.** *(The center of the paper;
   bounded by the no-behavior ceiling — an interpretation of a persistence result,
   not a behavioral proof.)*
4. **And the same lens is clinically useful** (N3): the propagator flags
   epileptogenic tissue. One framework, cognition → clinic.

---

## 3. The neglected dimension — `task_learn` (fix this everywhere)

The four phases are `rest_pre → task_learn → task_test → rest_post`. **Most of
the core trace work (N1) used only `task_test`** as "the task"; `task_learn` was
ignored. That is a real gap, and N2 is where it gets repaired:

- `task_learn` = **encoding** (patient is *shown* premise pairs).
- `task_test`  = **inference/retrieval** (patient *works out* novel pairs).
- The consolidation arc (`audit_103`) is the only place both phases are used:
  encoding `e = D_task_learn − D_rest_pre`, inference-specific
  `f = D_task_test − D_task_learn`. `task_learn` is the **reference that
  separates inference from encoding** — without it there is no N2.

**Open questions to seed in N1/N2 to-dos:** does `task_learn` leave its *own*
multiscale trace measured `rest_pre → task_learn → rest_post`? Is the OFC
localization the same for learn vs test? Is the β "both-probes" verdict
test-specific or does learn show it too? Treat the test-only core trace as a
**reference-phase choice that must be justified**, not a default.

---

## 4. Per-headline file spec (every headline `.md` must have these)

Each headline file is a **paper-section scaffold**, brutally honest, complete on
its own, in this order:

1. **§A Result in plain language** — result-oriented, interpretation-first.
   What we found and *why it matters*. Minimal numbers, **no p-value dumps**.
2. **§B Technical statement** — the measures, the C1–C5 controls applied, one
   block per subheadline. **Reference cached results; never hardcode tables**
   (see §5).
3. **§C Critical issues & powerful strengths** — outliers/leverage, confounds,
   retractions, what is airtight vs provisional. Lead with the weakness.
4. **§D To-dos & verifiables (for the owning agent)** — concrete, checkable,
   each tagged with the agent/chat that owns it (§6).
5. **§E Figure / representation ideas** — honor the figure rules in §5.
6. **§F Provenance** — bullet list (not a table): result → CSV path · generating
   script · timestamp · where the verdict is locked.
7. **§G Missing parts / open** — what we don't yet know (incl. `task_learn`).

---

## 5. Global rules (apply to every headline file)

- **No hardcoded result tables in `.md`.** State the result qualitatively in
  prose; for the numbers, **point to the cached CSV, name the script that
  generated it, and give the timestamp**. One source of truth = the CSV. This
  also fixes the "wall of unreadable p-values" problem. (If a number must appear
  inline for readability, e.g. a headline AUC, keep it to one and still cite the
  CSV.) Avoid markdown tables entirely — use bullets.
- **Frame the C1–C5 battery completely** (defined once in `00_CORE`; reference
  it, don't re-derive):
  - **C1** within-baseline split-half null (ρ^coph): is the cross-phase shift
    bigger than within-`rest_pre` noise?
  - **C2** drift-floor null (ρ^coph): bigger than slow within-session drift?
  - **C3** **matched-strength surrogate — MANDATORY, both probes**: bigger than
    what degree-preservation alone produces? *This is the referee; it has killed
    several headline claims.*
  - **C4** cross-probe restriction (ρ^coph): does restricting to cross-shaft
    pairs degrade the trace? (same-shaft control).
  - **C5** epi-zone exclusion (both probes, **secondary/mechanistic**): does the
    trace survive/strengthen when SOZ contacts are dropped? + the mandatory
    **node-count decimation control** before any "exclude subset → strengthens".
  - **Localization battery** (N1.3/N2.3): A1 hypergeometric + A3 matched-strength,
    R=1000, shaft-collapse, leave-one-out, within-band multiplicity (a-priori
    systems).
- **Figure rules:**
  - **Do NOT build figures around C4 cross-probe controls** — not informative to
    look at.
  - **DO feature the C5 epi-exclusion (X-epi)** comparison — it reads cleanly and
    sells the healthy-tissue point. Keep X-epi panels in the figure set.
  - Follow the project plotting guide (PDF-only, vector, no suptitle, mplstyle).
- **Outliers / per-patient leverage are a cross-cutting critical issue** —
  surface them, never hide:
  - **Pat_15** — right-hemisphere-only implant, β-LRG anti-aligned (the one
    sanctioned dropout exception); 2nd epi hub-patient.
  - **Pat_10** — β-anti at full β; detector hub-patient (marginal AUC).
  - **Pat_08** — δ Grassmann LOO leverage.
  - **Pat_02** — β LOO argmax-p driver (large, partly epi-coupled effect).
  - Policy: report LOO-max p; never let a single patient carry a "strong" tag
    (`feedback_no_single_patient_p_driven`).
- **"θ is anti-trace" is too simple — do not assert it.** At the cophenet probe
  θ's cohort median sits in the anti direction, but this is a single-layer,
  fragile, possibly leverage-driven sign; at raw FC θ behaves like every other
  band. Characterize θ (outlier sensitivity, layer-dependence) before saying
  anything stronger than "no trace on either probe." Flagged as a to-do.
- **Brutal honesty / results-in-the-Laplacian-framework**: raw FC is the
  comparison baseline, **never** a result; a claim significant-on-raw but only a
  trend on cophenet is reported as the cophenetic trend
  (`feedback_results_only_in_laplacian_framework`). Mark unverified measures
  "unverified"; lead each section with its binding limitation.

---

## 6. Agent routing (who owns what)

The PI runs one investigation chat per headline. When a headline file's §D has
open items, hand that file to its owner:

- **N1 trace** → the **localization** chat (OFC, systems, shaft-collapse) +
  the **white-matter** chat (gray/WM stratification, decimation controls) +
  **preprint-general-questioning** chat (now on the **null model** — C3 design,
  θ nuance, outlier leverage).
- **N2 encoding vs inference** → the **inference** investigation chat (arc,
  inference localization, low-γ focal encoding; duration control **RESOLVED**
  2026-06-22 — β arc duration-robust, α duration-suspect; τ-sweep of the arc next).
  No behavioral test (performance data unavailable).
- **N3 epileptogenic markers** → the **epi-marker-trace-analysis** chat
  (6-band detector, seed-free markers, cross-phase routing rigidity).
- **CORE methodology** → **preprint-general-questioning** chat (probe
  non-redundancy, band-resolution argument, C1–C5 framing, reviewer-proofing).

When you hand off, the agent should: finish the section, validate the cited CSVs
still match, spot critical issues, and propose figures — then return an updated
file. Do not run those investigations in *this* (crystallization) chat.

---

## 7. File index

- [`00_CORE_methodology.md`](00_CORE_methodology.md) — the methodology spine.
- [`01_trace.md`](01_trace.md) — N1, the multiscale persistence trace.
- [`02_encoding_vs_inference.md`](02_encoding_vs_inference.md) — N2, cognitive content.
- [`03_epileptogenic_markers.md`](03_epileptogenic_markers.md) — N3, SOZ read-out.
- *(replay narrative removed 2026-06-22 — former `04_replay_states.md` archived under `archive/2026-06/`; restore only if something surfaces.)*
- [`verification/`](verification/README.md) — verification briefs to hand to agents (one per investigation).
- [`epi-marker-analysis.md`](epi-marker-analysis.md) — **detailed backing doc** for
  N3 (full audit-by-audit epi investigation history). N3 is the crystallized
  headline; this is its evidence appendix.

Operational sources of truth remain the per-band briefs (`../bands/`) and locked
ledgers (`../locked/`). Headlines synthesize; they do not override the ledgers.

---

## 8. Verification briefs (hand these to agents)

The [`verification/`](verification/README.md) folder holds one self-contained brief
per investigation flagged as most interesting for neurophysiology / epilepsy. Each
follows the **mandatory 5-point critical preamble** (claim · null · strongest
alternative · does-the-null-control-it · falsification) + exact steps + data/scripts
+ pass/fail + caveats + what to return. Current set:

- **`verify_beta_spares_alpha_recruits`** (N1.6 bridge).
- **`verify_virtual_resection`** (epilepsy; mechanistic validation of N3).
- ~~`verify_stage2_arc_duration`~~ **RESOLVED 2026-06-22** — β inference arc is
  duration-robust (clean ratio regression + matched-strength); naive truncation
  null retired as invalid. N2.5 closed.
- **`verify_seedfree_epi_and_rigidity`** (N3.3 + cross-phase routing rigidity).

**No brain–behavior verification exists** — TI performance data is unavailable
(PI 2026-06-22). Do not scope one.
