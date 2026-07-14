---
name: talk-results-restructure-10slot-2026-07-13
kind: talk-structure-plan
era: IMCOH_ABS_COHORT_N10 (mst@0.20 recovery)
status: current
created: 2026-07-13
updated: 2026-07-13 PM   # grew 9→10 slots (higher-order β = own slide); folded the PM checkpoint pull
scope: THE current results-act restructure for the talk (pipeline-onward, 10 slots), built on the settled three results + the PM checkpoint's 4-pillar articulation + the narrative audit + this session's methodological Q&A (pedestal-vs-peak selectivity, detection-vs-form, higher-order β complementarity). The raw-dense/backbone graph confound is RESOLVED (same-graph apples control). RESUME ANCHOR after compaction — iterate on this. (Filename keeps "9slot" for pointer stability; it is now 10 slots.)
supersedes: the results-act structure of .agents/talk/REFOCUS_2026-07-13_three-results.md (the 8-slot in-place refocus)
pointers:
  - .agents/reports/2026-07-13_CHECKPOINT-PM-narrative-and-cascade.md            # THE governing latest state (4 pillars + through-line)
  - .agents/preprint/established_results/2026-07-13_settled-three-results.md     # the three results, locked (R1/R2 tables mid-cascade)
  - .agents/talk/2026-07-13_multiscale-narrative-audit.md                        # the 8 talk points + climax, graded
  - data/sparsified_arc/controls_ladder/cohort_gate.csv                          # pedestal-vs-peak, raw DENSE (the headline choice)
  - data/sparsified_arc/controls_ladder_apples/cohort_gate.csv                   # SAME-GRAPH control (everything on backbone) — the S4 rigor backup
  - data/sparsified_arc/coph_beyond_raw_s*/                                      # S6 higher-order β (coph⊥raw), script 21_coph_beyond_raw.py
  - data/sparsified_arc/figures/fig_scale_signatures.pdf                         # S5 + S7 curves (built)
---

# Talk results restructure — 10 slots (pipeline-onward)

**Head.** The results act is a **10-slot machine → number → meaning → [discriminate · scale ·
higher-order] → payoff → one-operator → closer** arc. The spine: **simple FC *detects* the
reorganization but responds *broadly* (a pedestal across all bands); the multiscale pipeline
earns its place three ways — it *discriminates* the bands (sharp α/β peak), *characterizes*
their scale (β scale-invariant vs α single-peak), and adds a *higher-order* coarse-scale view
the pairwise edges miss — and the same operator frames cognition and epilepsy.** Framing is
**information-DEPTH, NOT exclusive detection.** Through-line (PM checkpoint, say it): **raw FC
is a robust, *complementary* detector — never "fragile / blind / insufficient." The hierarchy
adds discrimination + characterization + a higher-order lens; it does not replace the edges.**

The three method-value pillars land as three consecutive beats (S4/S5/S6); then the scientific
payoff they enable (S7 inference held offline); then external validation (S8 epilepsy).

## The 10-slot sequence (S1 = pipeline, already built; S2–S10 to build/rewrite)

| Slot | Slide | Distinct job / message | Pillar / verb | Figure | Absorbs (old) |
|---|---|---|---|---|---|
| S1 | **The pipeline** *(DONE)* | the **machine** — FC → sparsify (mst@0.20) → diffusion swept over τ → cophenetic hierarchy. Dense degenerate; sparsify so τ resolves scales; propagator sums all paths (higher-order graph feature). | — | 7-panel seq *(built)* | 11 |
| S2 | **The measure: ρ^coph_sym** | the **number** — Spearman on cophenetic distance vectors, symmetrized over rest_pre split-halves; plants the **two orthogonal axes** (scale τ × cross-phase). *(Absorbs the ρ formula off the crowded pipeline slide.)* | — | formula + schematic | (ρ formula now on 11) |
| S3 | **The 4-fold taxonomy** *(KEPT standalone — user values it)* | the **meaning of the cross-phase axis** — a population's fate pre→task→post = anchor / trace / reset / reorganized, read off the tanglegram (parallel = kept place; crossing = scrambled). Conceptual menu, illustrative single-patient tanglegrams (disclosed-selected), NO anatomy/cohort claim. | — | tanglegram | 12 |
| S4 | **Detect ≠ discriminate: pedestal vs peak** | trace holds (10/10). **Raw FC detects it (β .024) — NOT "only multiscale sees it."** Raw = **broad pedestal** (positive all 6 bands, θ a top-3 effect +0.165, clears 4/6); cophenetic = **sharp α/β peak, θ negative**. θ = internal control. **Same-graph control lands (apples): even reading pairwise edges on the identical backbone they stay non-selective → the selectivity is the HIERARCHY, not the thinning.** | **1 · discriminate** | **pedestal-vs-peak profile** *(TO BUILD)* | 13 + 14 |
| S5 | **The scale-signature: only β is genuinely multiscale** | the **meaning of the scale axis** — τ-sweep grades the two survivors: β **scale-invariant** (16/16, .014→.001; Friedman p=.17 ns) vs α **single mesoscale peak** (12/16, .024→.007; Friedman p=.004). A single-scale probe can't tell them apart. | **4 · characterize** | `fig_scale_signatures` *(built)* | 15 |
| S6 | **The higher-order view: complementary, not redundant** *(NEW slide)* | is the hierarchy just a fancy repackaging of the edges? **No.** At the **coarse scale** the cophenetic read carries β grouping structure the pairwise edges lack (coph⊥raw β meso/coarse p=.032, effect +0.15→+0.24, BH-marginal q=.064; δ@coarse clean q=.039). And symmetrically **raw carries persistent structure the hierarchy lacks in *all* bands** → two **complementary** lenses. This is "higher-order graph feature" made concrete: the propagator sums all paths, so coarse τ reads multi-step community regrouping edges can't express. **The positive form of the anti-"raw is fragile" guardrail.** | **2 · higher-order** | coph⊥raw incremental panel *(TO BUILD, from `coph_beyond_raw_s*`)* | (was implicit in 11) |
| S7 | **The payoff: premises vs the order inferred, held offline** | encoding + inference **both persist**. **Inference-specific = α/β (NOT β-only, NOT δ — δ-inference CUT as a partial-corr artifact).** It is **mesoscale-emergent** (ns at s₁ .097 → sig at meso .005). **Detection-vs-form:** raw detects inference (.010); in the hierarchy it is a **mesoscale coarse re-grouping** = the shape of an abstraction. Held not replayed (10/10). **Distributed** — no hotspot (~88% co-move). | payoff | `fig_scale_signatures` (inference curve) | 16 + 17 + 18 |
| S8 | **One operator, two readouts** | band dissociation: cognition α/β, epilepsy δ/low-γ/β, **β the bridge**. Same operator reads SOZ (AUC δ.83, strength-independent, robust to sparsification). Contrast: **cognition distributed/placeless, but epilepsy keeps the one anatomical address in the whole talk — low-γ→SOZ (q.010)**. β **spares** SOZ; **α does NOT recruit SOZ** on the backbone (softens old "α recruits +0.41"). | 3 · external | epi figs | 19 + 20 |
| S9 | **Take-homes** | multiscale = **discrimination + scale-characterization + a higher-order lens, NOT exclusive detection**. Simple FC detects broadly (pedestal), robust + complementary; the pipeline discriminates bands (peak), grades scale, adds the coarse higher-order view. Inferred order = distributed offline abstraction. One operator frames cognition + epilepsy. | — | — | 21 |
| S10 | **Thanks** | — | — | — | 22 |

Setup reads **machine (S1) → number (S2) → meaning (S3)**; the two-axes idea S2 plants is
cashed by **S3 (cross-phase axis) + S5 (scale axis)**. The three method-value pillars are
consecutive (S4 discriminate · S5 characterize · S6 higher-order), then the payoff (S7), then
external validation (S8).

## The pedestal-vs-peak (S4) — numbers locked; graph confound RESOLVED

**Headline = raw DENSE** (`controls_ladder/cohort_gate.csv`) — the honest "simple FC as anyone
would actually compute it." T_test, med_obs · gate_p:

| band | raw_fc (DENSE, headline) | coph_meso (backbone) |
|---|---|---|
| δ | +0.072 · .024 ✓ | +0.122 · .116 |
| θ | **+0.165 · .053** | **−0.017 · .50** |
| α | +0.190 · .019 ✓ | +0.197 · .007 ✓ |
| β | +0.244 · .024 ✓ | +0.202 · .001 ✓ |
| low-γ | +0.116 · .032 ✓ | +0.039 · .216 |
| high-γ | +0.139 · .065 | +0.046 · .161 |

Raw = **broad positive pedestal** (all 6 bands positive; θ a top-3 effect that *misses only on
variance*; lights δ/low-γ non-carriers) → clears 4/6 (δ/α/β/low-γ). Cophenetic = **sharp α/β
peak, θ negative** → 2/6. **θ is the internal control** that makes "raw non-selective /
cophenetic band-sharp" non-circular.

**SAME-GRAPH control (apples), the rigor backup** (`controls_ladder_apples/cohort_gate.csv` —
every descriptor read on the identical mst@0.20 backbone B): raw-on-backbone STILL fires
broadly (α .042 ✓ / β .024 ✓ / high-γ .042 ✓; δ/θ borderline .053; and its **largest** cohort
effect is θ +0.214, which still misses), while coph-on-backbone is unchanged sharp α/β. ⇒ with
the graph held fixed, the selectivity is attributable to the **hierarchy read-out, not the
sparsification.** (This retires the earlier "pipeline-vs-pipeline, can't attribute" caveat.)
`fig_method_ladder` stays demoted to backup (its one unique-detection claim — β encoding .053
vs .003 — is LOO-fragile; see correction 3).

## The methodological corrections that MUST survive compaction

1. **GRAPH CONFOUND — RESOLVED (was correction #1).** The same-graph apples run
   (`controls_ladder_apples`, `14_controls_ladder_mst020.py` fixed 2026-07-13 to read every
   descriptor on backbone B) confirms the attribution: on the identical graph, pairwise edges
   stay non-selective, cophenetic stays sharp α/β. S4 now says "same graph, different read-out"
   — an upgrade, not a caveat. Keep the dense-raw pedestal as the intuitive headline and the
   apples run as the rigor backup.
2. **DETECTION-vs-FORM (S7).** Do NOT say inference is "invisible at the fine scale / you must
   climb the hierarchy to see it" — raw FC DETECTS β inference (.010). Correct: raw detects
   *that* it persists; in the cophenetic HIERARCHY it is a **mesoscale coarse re-grouping** (the
   *form*), the shape of an abstraction. Within-cophenetic mesoscale-emergence (s₁ ns .097 →
   meso sig .005) is clean (both scales on the same backbone).
3. **"ENCODING COPHENETIC-ONLY" IS DEAD → "robust reader," never "exclusive."** Raw FC encoding
   .053 (R=1000-confirmed), LOO 4/10 fragile; cophenetic 10/10 LOO robust. The hierarchy is the
   *robust* reader, not the *only* one. This is the S6/S9 guardrail: raw FC is robust +
   complementary, the overclaim "raw FC is insufficient/fragile" is BANNED.
4. **LOCALIZATION IS DEAD → DISTRIBUTED (incl. enc/inf).** β trace delocalized every scale
   (nothing clears BH; ~88% co-move); **tissue-DISTRIBUTED** (gray_gray/cross_wm/wm_wm all
   q≈.015; β spares SOZ). The enc→OFC / inf→cingulate double dissociation ALSO dissolves
   (enc→OFC β p=.156, 0/16; inf→cingulate β p=.23; only a sub-threshold α lean). **DROP
   OFC/cingulate "different cortex" entirely** — the honest "where" is DISTRIBUTED. The R2
   dissociation is by **SCALE**, not anatomy: encoding fine-legible + coarse-strengthening;
   inference mesoscale-emergent. (`17_localization_arc_mst020`, `19_tissue_pairclass_mst020`.)
5. **EVERYTHING IS SCALE-TESTED.** Laterality (β left, sig 6/16, peak ρ.74 @s16 — scale-robust ✅),
   duration (β inference ns everywhere ✅), tissue, SOZ all read per-scale. Present laterality as
   corroboration.
6. **DECOMPRESS the pipeline:** moving the ρ_sym formula onto S2 frees the (built, crowded)
   pipeline slide.

## The 4 pillars (PM checkpoint) → slot mapping (honest grades)

1. **Band-selectivity** (non-circular) → **S4**: coph 2/6 (α/β) vs raw 4/6; cross-validates R3 band dissociation.
2. **Higher-order β** (non-circular) → **S6**: coph⊥raw β persists at coarse scale, complementary (`21`).
3. **Epilepsy** (non-circular, external ground truth) → **S8**: community beats strength on SOZ AUC; low-γ→SOZ survives.
4. **Scale-invariance** (partly self-referential — acknowledge) → **S5**: β Friedman .17 vs α .004.

## Figures
- **Built (this session, PNG in `data/outputs/figures/talk/`; gen in `scripts/07_figures/talk_fig_*.py`):**
  - **`fig_measure_rhocoph.png` (S2)** — the measure as COLOURED radial dendrograms (reference /
    same-rest-half ρ=0.53 green / strength-matched-shuffle ρ=0.23 pale), colour = per-contact
    ρ^coph preservation. This is the "put the colour on the tree" visual the user asked for
    (like `fig_reinstatement_learn_mst020`); REPLACES the abstract two-axes plane as S2's main.
    Exemplar Pat_08 β. ⚠ optional: the low anchor is matched-strength (a hard null → pale-green
    0.23, not pure grey); swap to a full scramble for a punchier grey if wanted.
  - **`fig_rhocoph_tau_bands.png` (lasting-trace / S5 band-comparator)** — cohort ρ^coph_sym(s),
    all 6 bands over the 16-scale sweep + matched-strength floor; β held every scale (16/16),
    α mesoscale (12/16), θ/low-γ on the floor (0/16), δ/high-γ fine-only. Source
    `ms_mst020/cohort_gate.csv`. (User-requested band-comparison tool.)
  - **`fig_pedestal_vs_peak.png` (S4)** — raw FC pedestal (fires 4/6, θ a top effect) vs
    cophenetic sharp α/β peak (2/6, θ negative); θ shaded as internal control; y = obs − matched-
    strength null. Source `controls_ladder/cohort_gate.csv` (dense headline; apples printed for caption).
  - **`fig_coph_beyond_raw.png` (S6)** — L: β coph|raw grows coarse (higher-order, +0.15→+0.24);
    R: coarse-scale bars, each read-out sees what the other misses (complementary). Source `coph_beyond_raw_s*`.
- **Also built earlier:** `fig_scale_signatures` (S5, S7), 7-panel pipeline (S1), epi figs (S8), tanglegram (S3).
- **Demoted to backup:** `fig_method_ladder`.
- **S3 fix:** swap the *trace* tanglegram exemplar off **Pat_06 δ** (artifact band) → an α or β patient.

## NOT for the talk (newest, not yet load-bearing)
- `ms_marker_exploration/` (post-checkpoint): tests whether clever multiscale use beats the
  fixed-τ SOZ marker. FusShape_L2 prec@5 .74 (vs .58 fixed) is intriguing but **no
  label-shuffle-LOPO null verdict yet**; the PM checkpoint doesn't cite it; the scope's honest
  prior is "won't break the precision ceiling." Keep OUT of the deck until the null lands.

## Do-not list (guardrails)
No "only multiscale sees it" / "invisible to simple methods" / "multiscale-exclusive"; no "raw
FC fragile/blind/insufficient" (raw is robust + complementary); no "β trace → OFC" (delocalized;
OFC = encoding only, and even that dissolves on the backbone); no "inference β-only" (α/β) and no
δ inference-specific (cut); no drift / "second null" / "survives both nulls"; never rank bands by
×-null ratio; never call a band "false"; read per-scale, never best-scale. Null = matched-strength
ONLY.

## Flex points
- → 8 slots: merge S4+S5 or S5+S6 ("what the pipeline buys") if the act runs long.
- Upstream ripple (slide 5): diffusion-time ≡ topology (not 3 independent axes); "space" = distributed/network-scale.

## NEXT (resume here)
Structure LOCKED at 10 slots (user: higher-order β = own slide, 2026-07-13 PM). Craft order from
"the measure": **S2 (measure) → S3 (taxonomy, exemplar fix) → S4 (pedestal-vs-peak, build fig) →
S5 (scale-signature) → S6 (higher-order, NEW, build coph⊥raw fig) → S7 (inference payoff) →
S8 (one operator) → S9 (take-homes)**. Graph confound is resolved — no external dependency left.
