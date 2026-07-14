---
name: talk-review-tau-sweep-propagation
type: review
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery + τ-sweep)
status: current
created: 2026-07-13
updated: 2026-07-13
pointers:
  - .agents/reports/2026-07-12_mst020-recovery-arc.md
  - .agents/talk/PLAN.md
  - .agents/talk/slides/
scope: >
  Overall review of whether the mst@0.20 + τ-sweep recovery propagated onto the
  22-slide talk deck. NOT a fix pass — this is the dispatch doc for parallel
  per-slide fixing. Verdicts + ranked weak spots + self-contained per-slide punch-list.
---

# Talk review — did the mst@0.20 / τ-sweep recovery propagate?

**Head.** The recovery landed on the **spine** (slides 11 pipeline · 15 controls
ladder · 16 OFC · 21 take-homes are current; 14 nulls is drift-clean) but **not
on the whole deck**. Propagation ≈ **60%**. There is **one hard self-contradiction**
(slides 13/14 run "three co-equal reads / who's right" while 15 crowns cophenetic
the *sole* band-carver), **one now-false headline** (18 "inference persists, β-only"),
**one internal contradiction** (19's low-γ→cingulate *encoding* claim vs the ledger,
where cingulate = β *inference*), and a **drift leitmotif** still live in 06/08/PLAN.
The arc itself is sound; it is over-weight on adjudication (13–15) and on post-trace
results (16–19), and the `Slides.md` index (19 stubs) is out of sync with the 22-slide
deck. All fixes below bake in the three locked fork decisions.

## Locked fork decisions (2026-07-13)
1. **Grassmann — KEEP, made honest.** 13/14 keep the three-read "who's right" hook, but
   Grassmann is a **complementary at-a-scale (subspace) read** whose only defensible
   band under matched-strength is **β** (drop the "γ_low, δ" carve). It is *not* a rival
   the null must break — 15's spectral foil is a *different* read (effective resistance),
   so no contradiction if 13/14 stop calling Grassmann a band-carver that "disagrees".
2. **Encoding vs inference — inference stays flagship, DE-STALE.** Keep inference as the
   climax (18); fix the falsehoods: inference clears **δ/α/β** (not β-only), is raw-visible,
   and cophenetic's value-add is **localization (cingulate), not detection**. Encoding also
   traces (β) — acknowledge in one clause, don't cast it as pure nuisance.
3. **Sparsification — REVEAL AT 11.** Slides 5/7 stay clean conceptual motivation; the
   "dense FC is the degenerate LRG case → sparsify to mst@0.20 backbone" step appears at
   slide 11 (already there). 5/7 only need to *not contradict* it.

## Canonical numbers (source of truth for every fixer)
- **Trace (cophenetic, matched-strength, τ-swept 16 scales s≡τλmax∈[1,180], READ PER-SCALE):**
  - β = **scale-invariant** trace, clears **16/16** scales, cohort **p=0.001**.
  - α = **mesoscale** trace, clears **12/16** scales (s≈1.5), **p=0.007** (= the old τ_min 0.024, sharpened).
  - θ & low_γ = clean null (**0/16**). δ & high_γ = scale-collapse artifacts — **NOT claimed** as trace.
- **Encoding (T_learn, bug-fixed):** β **16/16** scales, cophenetic-UNIQUE (coph .003/.024;
  raw .053, clustering/geodesic/spectral FAIL). α readable from a BASIC metric (raw .032, clustering .010).
- **Inference (T_infspec_pe):** clears **δ/α/β** @meso (β 7/16, α 7/16, δ 10/16). β inference raw-visible;
  coph value-add = **cingulate localization** (audit_171). Raw T_infspec clears β only (fine scale).
- **Epilepsy marker (mst@0.20, AUC):** δ **.83** / γ_low **.82** / β **.745** (best-scale .87/.85/.90).
  Trace = α/β · epi = δ/γ_low/β · **β the only overlap**. Marker, not an improved detector (prec@5 flat).
- **Null:** matched-strength is the **SOLE** null. **Drift RETIRED** (directional task → trace monotonic
  → drift degenerate with the alternative). No "second null", no "double bar", no "strength · time".

---

## Global verdicts

### 1 · Propagation ≈ 60%
- **Done (spine):** 11 (prose; figure pending), 14 (drift-clean), 15, 16, 21.
- **Not propagated:** the *result numbers* (14/18/19/20 carry old ρ_sym-era p-values / AUCs),
  the *enc/inf flip* (17/18/19 + README + PLAN still "β-only inference"), the *Grassmann framing*
  (13/14 co-equal), and the *sparsification touchpoint* in the method slides (5/7 imply raw dense FC → LRG).

### 2 · Message cleanliness — 3 breaks + 1 noise source
- **HARD CONTRADICTION:** 13/14 "three co-equal reads, who's right" ⟂ 15 "cophenetic is the sole
  band-carver; raw/strength/clustering/geodesic/spectral all fail". Fork 1 resolves it (Grassmann =
  complementary, not rival) but 13/14 must be reworded so the contradiction actually disappears.
- **FALSEHOOD:** 18 "inference persists, β-only" — false under the bug-fixed pipeline (δ/α/β).
- **INTERNAL CONTRADICTION:** 19 "low-γ → focal cingulate ENCODING trace (8/8, q=0.035)" contradicts
  the ledger (cingulate = β INFERENCE; low-γ is a 0/16 null). Re-derive or cut.
- **NOISE:** drift leitmotif in 06/08/PLAN muddies the otherwise-clean "one null" message.
- Otherwise the core message ("one operator → multiscale hierarchy → band-specific held trace →
  same operator flags epilepsy") is consistent and lands.

### 3 · Structure — sound arc, two weight problems + a stale index
- Arc (question→lens→bet→machine→trial→localization→clinical→close) escalates well; thesis
  "it must be multiscale" is well-placed at 15; title/through-line coherent. **Keep.**
- **Over-weight A:** 13→14→15 = ~3 slides on "trust the measure", and "who's right?" is the hook on
  BOTH 13 (closer) and 14 (title) — de-duplicate.
- **Over-weight B:** 16→17→18→19 = four post-trace result slides; 19 partly re-treads 15/16.
- **Index drift:** `Slides.md` = 19 stubs for a 22-slide deck; 07/09/10 H1 titles are off-by-one from
  their filenames. Fix the index FIRST so parallel per-slide fixers don't collide on numbers.

---

## Weak spots (ranked)

1. **[HIGH] 13/14/15 adjudication bloat + duplicated "who's right?" hook.** Three slides, one idea.
   Minimum: de-dupe the hook (13 sets up the tension, 14 owns "who's right", or merge the framing).
   Consider whether 13's "three reads" + 14's null can tighten to two slides.
2. **[HIGH] 18 "β-only" is a live falsehood** on the current deck — top priority content fix.
3. **[MED-HIGH] 19 low-γ→cingulate encoding contradicts the ledger** AND 19 re-treads 15/16. Decide if
   19 earns its slot or folds into 16/18; either way the low-γ encoding claim must go/cite-fresh.
4. **[MED] Sparsification is load-bearing but single-line at 11** (Fork 3). Ensure that line carries the
   weight (dense = degenerate specific-heat single-peak; Villegas outlier) and that 05's "no modelling
   beyond the pairwise graph" doesn't read as "raw dense FC feeds LRG".
5. **[MED] `Slides.md` index ≠ deck (19 vs 22); 07/09/10 H1 off-by-one.** Blocks clean parallelism — fix first.
6. **[LOW-MED] Drift leitmotif in 06/08** — mechanical de-drift for a clean one-null message.
7. **[LOW] 08 and 11 both claim the 3D brain-tree hero** — assign ownership (bet=teaser still, pipeline=build).

---

## Per-slide punch-list (dispatchable; each row self-contained)

Legend: **✓** no change · **○** n/a · **⚠** fix. Priority P1 (content-critical) → P3 (housekeeping).

| # | Slide | Directive | Pri |
|---|---|---|---|
| 01 | title | ✓ none | — |
| 02 | ti-puzzle | ✓ none | — |
| 03 | seeg-dataset | ✓ none | — |
| 04 | functional-connectivity | ✓ none | — |
| 05 | higher-order-multiscale | ⚠ Soften "no modelling beyond the pairwise graph" / "just the graph you already have" so it means *no higher-order modelling*, NOT "raw dense FC feeds LRG". Keep the hypergraph-vs-multiscale contrast. Do NOT introduce sparsification here (Fork 3). One-line touch. | P3 |
| 06 | why-it-matters | ⚠ De-drift: "beats the nulls: strength · time" → "beats the null: matched-strength". Drop the "passage of time" half throughout (concept/on-slide/speech). | P2 |
| 07 | diffusion-lrg | ⚠ Numbering: H1 "Slide 6"→"Slide 7", frontmatter name→talk-slide-07. Keep conceptual (Fork 3 — no sparsify step here). Ensure the τ-sweep video doesn't *claim* a multiscale ladder on dense FC (lens-level illustration only). | P3 |
| 08 | the-bet | ⚠ "(strength · time)" → "(matched-strength)" ×2 (main concept + on-slide). Assign the 3D hero: bet = teaser, pipeline(11) = build. | P2 |
| 09 | imaginary-coherence | ⚠ Numbering only: H1 "Slide 8"→"Slide 9", name→talk-slide-09. | P3 |
| 10 | frequency-bands | ⚠ Numbering: H1 "Slide 9"→"Slide 10", name→talk-slide-10; reconcile 35%/45% Canva figure. | P3 |
| 11 | pipeline | ✓ prose current. ⚠ FIGURE rebuild (task #13, the RESUME point): insert |ImCoh|→mst@0.20 backbone panel between panels 2 and 3; convert single-τ D(τ) panel to a 2–3-snapshot τ-sweep morph. Confirm the one-line "dense degenerate → sparsify" motivation carries the weight (Fork 3). | P1 (fig) |
| 12 | tanglegram-categories | ⚠ Swap the TRACE exemplar from Pat_06 **δ** (artifact band) to an α or β patient (e.g. a β left-lateralized tracer); regenerate that tanglegram; keep the disclosed-selected caveat. | P2 |
| 13 | lasting-trace | ⚠⚠ (Fork 1) Keep the 3-read hook but: (a) Grassmann band claim "β, γ_low, δ" → **β only**, framed as a COMPLEMENTARY subspace read (not a disagreeing rival); (b) add per-scale τ-sweep to the coph verdict (β 16/16 p=0.001 scale-invariant; α 12/16 p=0.007 mesoscale); name mst@0.20; (c) FIG `fig_before_nulls_three_measure`: trim Grassmann over-claim + fix footer "strength / drift nulls" → "matched-strength null", regen PNG; (d) soften descriptive "blunt drift"; (e) ensure no wording survives that 15 will contradict. | P1 |
| 14 | the-nulls | ⚠ (Fork 1) matched-strength-sole-null is already correct — keep. Fix: (a) Grassmann band list "δ, β, γ_low" → defensible **β**, complementary framing; (b) replace old cophenetic gate p's (α .024 / β .032) with the per-scale trace numbers (β 16/16 .001, α 12/16 .007); (c) de-duplicate the "who's right?" hook shared with 13; (d) `fig_whos_right_three_methods` stays a 3-method matrix (Fork 1) but with corrected Grassmann bands — keep DISTINCT from 15's controls-ladder figure. | P1 |
| 15 | survivors-table | ✓ current (controls ladder, drift-clean, mst@0.20). | — |
| 16 | beta-ofc-consolidation | ✓ current. | — |
| 17 | encoding-vs-inference | ⚠ (Fork 2) Keep inference-flagship setup, but stop casting encoding as pure nuisance-to-partial-out: add one clause that encoding also traces (β, cophenetic). Keep the "seen pairs vs inferred order?" hook. No number changes (setup slide). | P2 |
| 18 | inference-persists | ⚠⚠ (Fork 2) "β-only" → inference clears **δ/α/β** (raw-visible for β). Replace single scalar (+0.091/+0.005, p=0.0098) with per-scale counts across the 16-scale sweep (β 7/16, α 7/16, δ 10/16 @meso). KEEP the localization dissociation (encoding→OFC vs inference→cingulate) + "value-add = localization NOT detection". Re-derive or CUT the low-γ→cingulate "8/8 q=0.035" (old ρ_sym; low-γ is 0/16 null). Scrub the residual "same reason drift is" analogy. | P1 |
| 19 | per-band-taxonomy | ⚠ (Fork 2 + numbers) α "gate p=0.024" → mesoscale 12/16 **p=0.007** (state as sharpened); add scale axis (β scale-invariant 16/16 .001; α mesoscale; θ/low-γ clean 0/16; δ/γ_high scale-collapse, not trace); "δ: clinical" → epi marker = **δ/γ_low/β**; DROP or re-derive the low-γ→cingulate ENCODING claim (contradicts ledger: cingulate = β INFERENCE); clarify "only β clears the full bar" = holds AND localizes. Consider folding into 16/18 (weak-spot #3). | P1 |
| 20 | epilepsy-marker | ⚠ AUCs δ .80 / γ_low .74 / β .69 → **δ .83 / γ_low .82 / β .745** (points 18, 25, 30). Keep the τ-buys-ranking-not-precision + mst@0.20-robustness framing (current). | P2 |
| 21 | take-homes | ✓ current (drift retired, controls ladder, band-specificity emergent). | — |
| 22 | outlook-thanks | ✓ none | — |

### Hub docs
| File | Directive | Pri |
|---|---|---|
| `Slides.md` | ⚠ Stale 19-stub outline for a 22-slide deck. Rebuild to the 22-slide index; fix the slide-6/10/11 stubs (sparsify at 11, Grassmann complementary, controls ladder). **Do this FIRST** — it's the scaffold that keeps parallel per-slide fixers from colliding. | P1 |
| `README.md` | ⚠ Slide-18 index line "18 · Inference persists (β only)" → drop "(β only)". | P3 |
| `PLAN.md` | ⚠ Board re-sync: delete drift refs (Verify L273, CANVA L230, Blockers L300); swap old gate p's for per-scale numbers; reframe 14/15 figs (Grassmann complementary, not co-equal); insert FC→sparsify(mst@0.20)→τ-scan into the slide-11 ladder/sequence; correct slide-18 to inference-not-β-only. Large — its own cleanup task. | P2 |

## Suggested parallelization order
1. **First, solo:** `Slides.md` index rebuild (unblocks clean numbering for everyone).
2. **P1 content wave (parallel):** 13, 14, 18, 19 (+ 11 figure rebuild). These carry the contradiction,
   the falsehood, and the ledger conflict — highest risk if a viewer sees them stale.
3. **P2 wave (parallel):** 06, 08, 12, 17, 20, PLAN.
4. **P3 housekeeping (parallel/batched):** 05, 07, 09, 10 numbering + README.
