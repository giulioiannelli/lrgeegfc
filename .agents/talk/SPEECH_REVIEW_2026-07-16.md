---
name: talk-speech-review-2026-07-16
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
status: active
created: 2026-07-16
updated: 2026-07-16
pointers:
  - .agents/talk/PLAN.md
  - .agents/talk/slides/
---

# Whole-deck speech review — timing, message coverage, critical spots (2026-07-16)

## Head

The 18 live slides carry **~4,640 spoken words ≈ 35.7 min at a measured 130 wpm / 30.9 min at a brisk 150 wpm — 6 to 11 minutes over a 25-min cap.** Nearly all the overage (≈ 640 s of the 641 s to cut) sits in **six slides that were written at roughly 2× the time budget their own SPEECH headers declare** — 13, 14, 15, 07, 06, 17. The fix is not invention: cut each of those back to the budget its header already states, plus light trims on 03/08/11b/16, and the deck lands at **~25.0 min @130 / ~21.7 @150 (~3,350 words)**. Every big message survives — the cuts remove restatement, cross-slide duplication, and nuance-overload, not content. Separately, there are **8 critical spots** (one is a hard numbering bug; two are pending user decisions).

## 1. Timing — now vs each slide's own stated budget vs target

Delivery-rate assumption: a dense technical talk with figure-pointing runs **130–150 wpm**; plan at ~135 and rehearse. "budget" = the `~N s` the slide's own SPEECH header declares (blank = none declared).

| slide | now | budget | target | cut |
|---|---:|---:|---:|---:|
| 01 title | 0s | — | 7s | add a 1-line open |
| 02 ti-puzzle | 93s | — | 90s | keep |
| 03 seeg-dataset | 91s | — | 75s | −16s |
| 04 functional-conn | 57s | — | 55s | keep |
| 05 higher-order | 70s | 65s | 65s | on target |
| 06 why-it-matters | 178s | 140s | 130s | **−48s** |
| 07 diffusion-lrg | 180s | 100s | 105s | **−75s** |
| 08 form-mirrors | 110s | 95s | 95s | −15s |
| 09 imaginary-coh | 78s | — | 75s | keep |
| 10 frequency-bands | 54s | — | 54s | keep |
| 11 pipeline | 73s | — | 70s | keep |
| 11b the-measure | 109s | — | 85s | −24s |
| 12 tanglegram | 67s | — | 65s | keep |
| 13 lasting-trace | 234s | 110s | 110s | **−124s** |
| 14 detect-v-disc | 206s | 110s | 105s | **−101s** |
| 15 encoding-v-inf | 264s | 130s | 125s | **−139s** |
| 16 epilepsy | 100s | — | 90s | −10s |
| 17 closer | 178s | 100s | 100s | **−78s** |
| **TOTAL** | **35.7 min** | | **25.0 min** | **−10.7 min** |

Release valve if a live run-through still spills: slide 16 (epilepsy) is marked cuttable in its own notes — dropping it buys ~1.5 min and takes the deck to ~22–23 min.

## 2. Big-message coverage — complete; the 06→17 arc closes

All twelve core messages are present and land in the right place, and the "why it matters" arc (slide 06) is paid off one-for-one in the closer (slide 17): cognitive-maps→multiscale lens, offline-consolidation→encoding-vs-inference, replay→held resting trace, data-efficiency→AI-inspiration outlook, the detect→characterize gap→the whole results act, matched-strength→the running guardrail. No message is missing. The only coverage issues are three **redundancies** (reclaim time, below) and one **micro-gap** (no scripted opening line on the title slide).

Redundancies to reclaim: (a) the **symmetrised-baseline** explanation is given in full on both 11b and 15 — keep it on 11b, cut to a half-clause on 15; (b) **held-not-replayed** is stated at length on 15 and again on 17 — full on one, a clause on the other; (c) **"each contact ≈ a ~mm population"** appears on 03, 12, and 16 — keep on 12 (where it does work), trim elsewhere; (d) the **matched-strength refrain** appears on 6+ slides — that is a deliberate refrain, fine, but only slide 06 needs the full plain-language sentence.

## 3. Critical spots (ranked)

1. **[Blocks the 25-min goal] Length: 35.7 min @130 / 30.9 @150 — 6–11 min over.** ~640 of the 641 s to cut are in six slides (13/14/15/07/06/17) that overran their own declared budgets by ~2×. Fix = §4 cut plan.
2. **[Must-fix — numbering bug] Two slides both title themselves "Slide 9."** `09_imaginary-coherence` reads "Slide 9" (correct) and `10_frequency-bands` also reads "Slide 9" (should be 10); there is no in-file "Slide 10", and `11b` reads "Slide 11b (S2)". The in-file `# Slide N` headers are stale/duplicated and will mis-map to Canva pages and confuse the presenter. Renumber the headers 1–17 to match filename order (09→9, 10→10, 11→11, 11b→11b or fold into 11, …).
3. **[Structural — biggest single lever] Slide 15 does three jobs.** Setup + encoding result + inference result + a held-not-replayed coda → it is the longest slide (264 s, ~2× its 130 s budget) and it duplicates 11b's symmetrisation and 17's held-not-replayed. Cutting the duplicated beats alone reclaims ~90 s.
4. **[Density] Three dense method slides run back-to-back (07, 11, 11b ≈ 6 min) heavy with vocabulary** — specific-heat C(τ) mesoscale ladder, lowest-common-ancestor, symmetrised partial Spearman, propagator series. Keep the intuition, shed the jargon; the audience needs the *what*, not the derivation. 07 alone is 80 s over budget.
5. **[Pending decisions — 2, block a clean lock]** (a) **Slide 06 Kramer citation** — the on-slide footer cites a different Kramer 2012 (PNAS "seizures self-terminate…") than the "Kramer & Cash" the speech names; reconcile to one. (b) **Slide 07 reach figure** — the on-slide `talk_fig_spatial_reach.py` curve floors at the electrode pitch (3.5 mm = a single contact), so the delocalisation point must be anchored at the s = 1 dashed line, not at a "~15 mm finest" claim (that number is a *different*, max-diameter figure); decide whether to keep the pitch-flooring figure with s = 1 framing or swap to the max-diameter version.
6. **[Watch — weakest spoken claim] Slide 14's low-γ "mis-read" is partly confounded** (dense Grassmann vs the mst@0.20 backbone, per the slide's own Careful). It is the one place a sharp questioner can pull a thread. Lead the clean, same-backbone α spine; demote the γ-low mis-read to a single sentence.
7. **[Minor] Title slide has no scripted first sentence.** A talk's opening line should be chosen, not improvised — add one line (≤ 10 s).
8. **[Minor — internal consistency] prec@5 multiplier "~7×" (slide 16 speech/on-slide) vs "~5×" (§5 figures note)** — the talk is internally consistent at ~7×; just don't let the §5 note leak onto a slide.

## 4. Cut plan — per slide (specific)

Heavy cuts (to the header budget):

- **15 (264→125 s, −139).** Cut the full symmetrisation paragraph to a half-clause (it is on 11b). Compress the held-not-replayed coda to one clause (it is on 17). Say the amber/teal encoding-vs-inference contrast once, not twice. Keep: the disc-golf callback, the three figure beats (space / tinted tree / scale curves), the "on the null at fine scale → lifts at mesoscale, β AND α 7/16 each → it was never β alone", the abstraction close, and the detection-tied→value-add-is-FORM point. Worked rewrite drafted — see chat.
- **13 (234→110 s, −124).** Three figures keep their three jobs, but compress each beat. The Fig 2 violin "reads-both-ways" passage (the subtlest, longest) collapses to: "grey ≠ silent — individual patients trace strongly, low-γ out to 0.6, but not the same way, so the cohort median doesn't hold." The closing "so — the result, and what it means" paragraph restates the figures; halve it.
- **14 (206→105 s, −101).** Compress the two-paragraph low-γ / spectral-clustering beat to ~2 sentences and demote it below the α spine. Tighten the opener and the "relational lens, not a superset" close (keep the point, lose the length).
- **07 (180→105 s, −75).** Keep the community-size / effective-spatial-scale beat (user-required) but cut the "left of the line / cross it / push further" walk to two sentences. Tighten the micro-meso-macro illustration and the higher-order-for-free paragraph.
- **06 (178→130 s, −48).** Compress each of the three literature anchors to two sentences; shorten the data-efficiency opener; trim the detect→characterise enumeration (flat/peaked/emerging) to a phrase.
- **17 (178→100 s, −78).** The take-home sentence-chain is a good rapid recap — trim the clinical restatement (already slide 16) and the detect→characterise re-explanation. Trim the outlook enumeration; keep the AI-inspiration close (the payoff).

Light trims: **03** (−16, tighten the phase/kHz list), **08** (−15), **11b** (−24, shed LCA/symmetrisation vocabulary to intuition), **16** (−10). **01** add a one-line open.

## 5. Method

Word counts and budgets extracted mechanically from each slide's `4. SPEECH … 5. FIGURES` block (bracketed stage cues and the `Careful`/audit rails excluded — those are notes-to-self, not spoken). Rate band 130–150 wpm. This file is the cutting checklist; apply slide-by-slide and re-measure.
