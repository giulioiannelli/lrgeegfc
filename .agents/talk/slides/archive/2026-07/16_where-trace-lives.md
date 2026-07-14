---
name: talk-slide-16-where-trace-lives
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
slide: 16
status: draft
updated: 2026-07-13
canva: NEW — the REVIVED localization slide (user 2026-07-13: "add the localization slide back, the
  audit found something nice"). Source: .agents/reports/2026-07-13_localization-cohort-statistic-verdict.md
  ([[localization_cohort_statistic_audit_2026_07_13]]). Supersedes old 16 (β held/placeless/lateralized).
---

# Slide 16 — Where the trace lives: a coarse home, and an honest null

1. TITLE
Where does the trace live? A coarse, left home — the fine address belongs to the disease

2. MAIN CONCEPT
- **Head (the juice):** we stress-tested our own result. The multiscale framework had seemed to *wash
  out* the anatomy — nothing localized. Was that real, or an artifact of the test? We ran the same
  localization four different ways on the identical graph and null. The answer is honest and, in the
  end, positive: our first "nothing localizes" was **partly a broken test**, but fixing it does **not**
  bring β→OFC back — β has **no fine anatomical home**. What it *does* have, once you match the
  parcellation to the scale, is a **coarse** one: the **left hemisphere**. And it is a **cognitive**
  trace, not an epileptic one — it *spares* the seizure zone. The only band that localizes to a
  specific place is the **disease** band.
- **The self-check (why this slide exists).** A single cohort statistic can hide a real effect *or*
  invent a fake one. So we held the backbone and the matched-strength null fixed and varied **only the
  cohort statistic** — four of them (pooled median · per-patient Wilcoxon · Stouffer · a
  consistency-favouring sign test), each read under whole-grid BH + leave-one-out. Disagreement between
  them is diagnostic.
- **Finding 1 — the wash-out was *part* artifact.** Our first localizer used a per-patient Wilcoxon
  across the 5–8 patients that cover each region; at 5 patients its p-floor is 0.031, so a region seen
  in five patients **cannot** clear correction no matter how strong the effect. A genuinely
  under-powered test.
- **Finding 2 — …but β→OFC does not return.** Under the powered statistics β sits at the **~2% noise
  floor** in every region and every target (trace, encoding, inference), and it is **density-invariant**
  — zero survivors at 10%, 20%, and 50% backbone density. The pre-registered β→OFC, encoding→OFC,
  inference→cingulate all fail. **β has no fine home; that verdict now stands for the right reason.**
- **Finding 3 — β is NOT placeless: a coarse, left home.** The fixed fine atlas was the wrong ruler —
  the diffusion coarse-grains space, so a coarse-scale trace must be tested at a coarse parcellation.
  Do that, and β **enriches in the left hemisphere** (z = +1.38, pooled p = 0.005, sign-consistency
  p = 0.005, **leave-one-out worst = 0.005**, 5/7 patients, at the mesoscale s ≈ 11; the right
  hemisphere is *depleted*). This is the between-patient laterality (ρ ≈ 0.69) surfacing as anatomy. β
  has a **coarse — hemisphere-level — home (left), and no fine one**: you have to zoom the *anatomy*
  out to a whole hemisphere to see it, which is why the fine atlas missed it. (Coarse in *space*; the
  signal itself peaks at the mesoscale s ≈ 11, so we do **not** claim a clean scale-to-granularity
  diagonal — the coarse spatial map is a good *lens*, not a proven one-to-one.) *(Honest: marginal
  under the strictest whole-grid correction — only two hemisphere units,
  q = 0.119 — so we lean on the a-priori left-lateralization hypothesis + the LOO robustness, not on a
  strict-BH win.)*
- **Finding 4 — β is a *cognitive* trace, not an epileptic one (three ways).** It **spares** the seizure
  zone: on non-SOZ pairs the β trace clears the null at 16/16 scales, in the SOZ core only 1/16; physically
  removing the SOZ contacts and rebuilding leaves it essentially intact (14/16, p = 0.002); and a
  size-matched random-removal control shows dropping the SOZ is no more disruptive than dropping the same
  number of random contacts. **β lives on healthy cortex.**
- **Finding 5 — α is genuinely placeless**, at every resolution. And the **one** specific anatomical
  address in the whole talk belongs to the **disease** band: **low-γ → SOZ** (q ≈ 0.04; it matches the
  seizure-marker AUC of 0.82, and it *vanishes* when you remove the SOZ). The broken Wilcoxon had hidden
  it. That is the hand-off to the next slide. *(Caveat: low-γ→SOZ may be a stable-tissue **anchor** rather
  than a task-trace — matched-strength controls connection strength, not tissue stability.)*

3. ON-SLIDE TEXT
did the scale framework wash out the anatomy? — we ran the localization 4 ways (same graph + null)

the old test was under-powered (Wilcoxon floor at K=5 patients) — so "nothing localizes" was part-artifact
BUT β → OFC does NOT return: β at the ~2% noise floor everywhere, at 10/20/50% density (β has NO fine home)

β is NOT placeless — a COARSE (hemisphere-level) home:  LEFT hemisphere  (z=+1.38, p=.005, LOO=.005, 5/7 pat, meso s≈11)
  = a coarse SPATIAL home (fine atlas missed it) + the laterality (ρ≈.69) as place  [no clean scale↔granularity diagonal]
β is COGNITIVE, not epileptic — SPARES the SOZ (non-SOZ 16/16 · SOZ-core 1/16 · survives SOZ removal 14/16)
α — placeless at every resolution
the ONE focal address is the DISEASE band → low-γ → SOZ (q≈.04, = the seizure marker)  → next slide

4. SPEECH
Now, where does the trace live — and here I want to show you how we check ourselves. The multiscale
framework had seemed to wash the anatomy out: nothing localized. That could be real, or it could be the
test. So we ran the same localization four different ways, on the identical graph and the identical null,
changing only the cohort statistic. First thing we found: our original test was under-powered — a
per-patient Wilcoxon across five-to-eight patients literally cannot reach significance at five patients,
regardless of the effect. So "nothing localizes" was partly a broken test. But — and this matters —
fixing it does not bring beta-to-orbitofrontal back. On every valid statistic beta sits at the noise
floor, in every region, at every backbone density. Beta has no fine anatomical home, and now that stands
for the right reason. Here's the nice part. The fixed fine atlas was the wrong ruler: the diffusion
coarse-grains space, so a coarse trace needs a coarse map. Test at the hemisphere level, and beta
lights up — in the left hemisphere, robustly, leave-one-out and all, while the right is depleted. That's
the between-patient laterality showing up as actual anatomy. Beta isn't placeless; it has a coarse home —
you have to zoom the anatomy out to a whole hemisphere to see it — a left one. And it's a cognitive trace, not an
epileptic one: it spares the seizure zone three different ways — it lives on the non-seizure pairs, and it
survives even when we physically cut the seizure contacts out. Alpha, by contrast, is genuinely placeless
at every resolution. And the one sharp anatomical address in this entire talk belongs to the disease band
— low-gamma marks the seizure zone — which is exactly where we're going next.

Careful: NEVER say "localization entirely dead / everything delocalized" — state PER BAND (β coarse-left,
α placeless, low-γ→SOZ). β→OFC / encoding→OFC / inference→cingulate stay RETIRED — do NOT reinstate a
FINE β localization; β's home is COARSE (hemisphere) only. β-left is MARGINAL under strict whole-grid BH
(q=.119, 2 units) — present it as robust-under-pooled+LOO + a-priori laterality, NOT as a strict-BH win;
say "we lean on the a-priori left hypothesis." Do NOT claim a clean scale↔granularity DIAGONAL (β's
hemisphere signal is at the mesoscale s≈11, not the coarsest). low-γ→SOZ MAY be an ANCHOR not a trace
(matched-strength controls strength, not tissue stability) — flag it. Null = matched-strength ONLY. The
SOZ-independence (spares the SOZ) is the β result; low-γ is the mirror (its trace IS the SOZ, vanishes
without it). Do NOT front-run the epilepsy AUCs / detector — that is slide 17.

5. FIGURES
- **MAIN (BUILD) — the scale-matched laterality panel.** β enrichment vs diffusion scale for LEFT vs RIGHT
  hemisphere: L rises to z ≈ +1.4 at the mesoscale (p=.005, LOO=.005), R depleted (negative every scale).
  Source `data/sparsified_arc/localization_bakeoff_mst020_glad/hemisphere_bakeoff.csv` (band=beta,
  grouping=hemisphere, target=trace). NOT YET BUILT — needs a small talk figure (band-color left vs grey
  right, sig markers). PNG fine.
- **SUPPORTING — the 4-statistic bake-off (delocalization is rigorous).**
  `data/outputs/figures/localization_bakeoff/fig_localization_statistic_bakeoff.pdf` (gen
  `scripts/07_figures/fig_localization_statistic_bakeoff.py`): system × statistic, showing β at the noise
  floor while low-γ survives. Use as the "we checked four ways" backup / the finding-2 anchor.
- **OPTIONAL — β spares the SOZ.** A small non-SOZ (16/16) vs SOZ-core (1/16) contrast, or the existing
  `data/outputs/figures/talk/fig_trace_f_soz_divergence.pdf` relabelled to "β spares the SOZ" (drop the
  old "α recruits" framing — softened on slide 17).
- **OPTIONAL — a left-hemisphere brain glyph** (β trace weight, left>right) to carry the "coarse home"
  beat visually. Build in Canva or reuse `data/outputs/figures/talk/fig_trace_laterality_brain_beta.png`.

6. REFERENCES
- The four-statistic cohort-localization bake-off + scale-matched parcellation ladder + SOZ-independence
  controls are ours (`23_localization_statistic_bakeoff_mst020`, `24_localization_rawcoph_bakeoff`,
  `19_tissue_pairclass`, `25_node_exclusion`; lib `utils/metrics/cohort_localization.py`;
  verdict 2026-07-13). Laterality driver = audit_148 (ρ ≈ 0.69).

7. CANVA STATUS
NEW slide, REVIVING localization (user directive 2026-07-13 PM). This is the methodological-honesty +
positive-reframe beat: we stress-tested the "wash-out" (4 statistics, same graph/null), found it was
part-artifact, confirmed β has no FINE home (β→OFC stays retired, density-invariant), but recovered a
COARSE left-hemisphere home (scale-matched) that is SOZ-independent (β = cognitive, spares the seizure
zone), with α placeless and the one focal address being low-γ→SOZ (the disease band). Missing on deck:
BUILD the scale-matched laterality panel (hemisphere_bakeoff.csv), place the 4-stat bake-off backup,
compress on-slide text, presenter notes. Hands off to slide 17 (the disease band). ⚠ All localization-audit
outputs are UNCOMMITTED (new work off audit/cohort-n10-diagnostic).
