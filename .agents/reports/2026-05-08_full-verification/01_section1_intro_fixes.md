---
type: report
status: current
date: 2026-05-08
era: IMCOH_ABS / COHORT_N10
section: 1
---

# Section 1 — fix list

**Head.** §1 is mostly clean: ImCoh swap framing, n=10 reads as the live cohort
through the second-page summary, and the band-accurate result preview at p.3
matches the audit. The remaining issues are minor: the `\textsc{Notes\textsubscript{MSC}}`
reference is consistent throughout (no fix), but the closing paragraph's claim
"controlled cohort claim at α, β, and low-γ via the per-pair correlation on D(τ)"
should be checked against the §5.3 anatomy line that already disclaims the
α qualification (mixed CTM alignment after the §5.5 Pat_07/Pat_14 disclosure).

## Confirmed (no action)

- "n=10" framing implicit in the reference to NOTES_MSC pipeline (the patient cohort and statistical protocol are stated as `unchanged` from the predecessor); cohort identity is reconfirmed in §2 via per-figure context.
- ImCoh swap framing (Nolte 2004, Ewald 2012, Bastos & Schoffelen 2016) is consistent with `imcoh-guide.md` and the reference list.
- The four-phase rest-task-rest naming (rsPre / taskLearn / taskTest / rsPost) is consistent with `terminology.md`.
- The three-axis multiscale framing (Topology / Frequency / Space) is the canonical §1 framing carried into §3 and §5.
- p.3 closing-paragraph result preview: "directionally consistent task-shaped pattern in four of six bands that survives drift controls" → matches Tab 2 / Tab 3 (δ soft trace, α / β / low-γ trace bands, θ drift-only, γ_h borderline). ✓
- p.3 closing-paragraph LRG headline: "controlled cohort claim at α, β, and low-γ via the per-pair correlation on D(τ) (within-probe BH-q = 0.027 at m = 6)" → matches §5.3 result table. ✓
- p.3 anatomy headline: "γ_l ctx-lh-fusiform Bonferroni-survived (m = 48)" → matches §5.5 result. ✓

## Numerical corrections (action: writing agent)

- None required at the §1 numerical level. All preview numbers (BH-q = 0.027 at m=6, m=48 anatomy correction) reproduce in `_section5_joint_bh.md` (m=48 family) and `_section5_5_verify.md` (anatomy m=48 Bonferroni).

## Framing rewrites (action: writing agent)

- **§1 closing paragraph, p.3, "with a per-band geometric decomposition (β multifacet, low-γ heights-only, α per-pair only)"** → currently consistent with §6.1 reading. Keep as-is. No rewrite.
- **§1 closing paragraph, p.3, "an anatomical headline at low-γ ctx-lh-fusiform that survives Bonferroni correction across the (band, region) cell space"** → §5.5 final prose discloses the Pat_02-double-dependence (47% fusiform implant + 50% per-patient trace rate). §1 preview is one level higher and does not need that disclosure, but a one-clause hedge would prevent over-reading: "Bonferroni correction across the (band, region) cell space (with a disclosed Pat_02 anatomical-density dependence in §5.5)". Optional.
- **NOTES_MSC framing** (p.2): "scientific question, experimental design, patient cohort, LRG framework, inter-phase comparison methodology, and statistical testing protocol are all `unchanged`. The single revision concerns the FC estimation step." → accurate; the entire MSC vs ImCoh comparison machinery in §2 / §3.2 is a `comparison`, not a re-derivation. Keep.

## Caveats to add (action: writing agent)

- None required at §1 level. The preview-paragraph caveats (drift controls, cross-baseline, joint multiple-comparison) are all developed in §4.6 and §5.3, which is the right place.

## Figure actions (coding agent)

- None. §1 contains no figures.

## Deferred / questions

- The §1 / §3.1 phrase "n=10" is **consistent**, but a careful reader will notice that the §2/§3.2 figures show "all six patients" (Pat_02, 03, 05, 06, 07, 08). The contradiction is at the §2 figure-caption layer, not at §1. See `02_section2_fc_fixes.md` for the cohort-relabel action items.
- The "scientific question, ... cohort ... unchanged" framing in §1 implicitly assumes the n=10 lock is already established in NOTES_MSC. Pat_14 was vendor-replaced 2026-04-25 (per `CLAUDE.md` + memory `audit_recovery_2026_04_25.md`); if the predecessor NOTES_MSC pre-dates the replacement, a sentence acknowledging "current cohort: Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15 (n=10, after Pat_14 vendor-replacement on 2026-04-25)" would close the loop. Optional and depends on whether NOTES_MSC has been updated.
