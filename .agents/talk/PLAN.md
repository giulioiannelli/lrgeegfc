---
name: talk-sprint-plan
type: plan
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
status: active
created: 2026-07-11
updated: 2026-07-11
pointers:
  - .agents/reports/2026-07-07_talk-structure-20min.md
  - .agents/talk/slides/README.md
---

# Talk sprint — execution plan & coordination board (2026-07-11)

## Head

**One-day sprint to a mostly-ready 20-min Canva deck.** Five lanes work in
parallel off *this one file*. The **conductor** (main slide-editing chat) and
the user review the deck slide-by-slide, locking each slide's content + figure
choices; as each slide locks, the conductor writes concrete tasks into the
**Figure**, **Reference**, and **Verify** lanes, which run in their own chats.
The deliverable is **not the Canva build** (the user does that) — it is: every
slide's content finalized, **every figure gathered into one folder with its
generating code identified**, and every number verified — all under `data/` +
`scripts/` so it rides the travel-bundle / git sync to the portable machine.

**How to read this file:** find your lane in **§Lanes**, read **§Protocol**,
then work your rows in the **§Slide board** + **§Figure manifest**. Mark your
own boxes done. Don't touch other lanes' rows.

---

## Goal (end of today)

- [ ] All 22 slides content-locked (structure fixed + text/speech finalized).
- [ ] **Every figure gathered in `data/outputs/figures/talk/`** (one folder).
- [ ] **Each figure's generating script identified + path recorded** here.
- [ ] External (from-paper) figures fetched into `data/outputs/figures/talk/_external/`.
- [ ] Every load-bearing number verified against the **latest ρ_sym-era** artifacts.
- [ ] Nothing outside `data/` + `scripts/` (so it all syncs).

Deck arc + numbers SSOT: `.agents/reports/2026-07-07_talk-structure-20min.md`.
Per-slide specs: `.agents/talk/slides/*.md`. This file is the **coordination
board**, not the content — content lives in the slide files.

---

## Lanes

| Lane | Runs in | Mandate | Primary output |
|------|---------|---------|----------------|
| **M — Conductor** | **this chat** | Review structure, then lock each slide's content + figure list *with the user*. **Score each slide against the Canva deck** (page + % + what's missing). Write F/R/V tasks into the board as slides lock. Keep the board current. | Updated `slides/*.md`; filled board |
| **W — Writing** | lateral chat | Finalize each **locked** slide's on-slide text + spoken script so it's clear and paste-ready. Polish prose only — do not change the locked scientific content. | Final text in `slides/*.md` |
| **F — Figures (ours)** | lateral chat | For each **locked** figure: finalize/regenerate it, **copy the final PDF into `data/outputs/figures/talk/`**, and record its generating script path here. | Gathered figs + code map |
| **R — References (external)** | **Conductor, inline** (user 07-11: "do the work of lane R, don't delegate") | Verify/pin external citations (WebSearch→DOI) + extract published figures from local PDFs (pdftoppm+crop, verify visually) into `…/talk/_external/`; record source + license/attribution in CITATIONS.md. | External figs + citations |
| **V — Verify** | lateral chat | Re-confirm every load-bearing number **fresh from the latest artifacts** (never trust prior prose/CSV — recompute/cross-check). Flag any mismatch. | Verified ✓ / flags |

The conductor (M) is the only lane that *locks* a slide. W/F/R/V act **only on
locked slides.**

---

## Protocol (read before editing this file)

1. **Work your lane only.** Edit only your lane's cells/rows. The conductor owns
   the slide board's `Locked` column and adds new figure rows.
2. **The gate is `Locked`.** A slide's W/F/R/V work is actionable **only when its
   `Locked` box in the §Slide board is `[x]`.** `[ ]` = conductor still working it → wait.
3. **Mark done in place.** Flip `[ ]`→`[x]` and append `— <lane> <date> · <one-line note>`.
   Example: `[x] fig in talk/ — F 07-11 · regenerated from fig_arc_a.py`.
4. **Small, localized edits.** Multiple chats share this file; touch one row at a
   time, don't reflow the whole doc.
5. **Blockers don't stall you.** Log them in **§Blockers** and move to your next row.
6. **Figures ride `data/`; code rides `scripts/`/`src/`.** Never put binaries in
   `.agents/`. Record code as repo-relative paths.
7. **Always score against Canva.** Every slide's part-7 (CANVA STATUS) and the board's
   `Canva pg · %` column track its state on the deck **The Multiscale Shape of Neural
   Inference** (design `DAHO6jbUdrY`). Re-read the deck when scoring — don't guess.
8. **Slide files are PLAIN 7-part** (title · main concept · on-slide text · speech ·
   figures · references · Canva status). See `slides/README.md`; exemplar = R-7.
9. **Methods slides carry their formulas as compile-ready LaTeX** in the ON-SLIDE
   TEXT block (the user compiles them with a Canva LaTeX tool → image). Define
   symbols before using them (build the ladder: e.g. power → coherency → the
   quantity we use). Match the LaTeX to what the code actually computes — verify
   against the implementation, don't paraphrase a textbook. On-slide prose stays
   in **bullet-list** form alongside the formulas. (Convention set 2026-07-11 on
   slide 09; applies to slides 07/09/10/11 and any future methods slide.)

---

## Slide board

`Locked` = conductor froze the slide's content + figure list (the gate for all
other lanes). `W` = text/speech finalized. `V` = numbers verified. Figure status
lives in the §Figure manifest.

All 22 slides now drafted in the plain 7-part format (`slides/NN_*.md`). `Locked` = the
conductor + user have frozen the slide in review (not yet — content is drafted, review pending).
The **figure list is SET** regardless (Lane F may start). `Canva pg` = the matching page on
the existing deck (NEW = build fresh).

| # | Slide | Canva pg · % | Locked | W | V | Notes |
|---|-------|:-------------|:------:|:-:|:-:|-------|
| 01 | Title | pg 1 · ~70% | [ ] | [ ] | [ ] | drafted |
| 02 | TI puzzle (disc golf) | pg 2 · ~45% | [ ] | [ ] | [ ] | RESTRUCTURED → 7-beat disc-golf reveal; chain + 3 bracket rounds (Lane F), photo swap (Lane R) |
| 03 | sEEG dataset | pg 3 · ~35% | [ ] | [ ] | [ ] | + Reber precedent, + double-probe |
| 04 | Functional connectivity | pg 5 · ~55% | [x] | [x] | [x] | ✅ LOCKED M 07-11 — FC matrix MOSAIC 4pat×4phase (generic FCᵢⱼ cbar; sfdp matrix↔net DROPPED as too sparse; 3D brain→backup); text now anchors real net-neuro results (small-world/rich-club/communities) + functional≠structure (Honey'09); dropped struct-vs-func dichotomy. W: text+speech final (conductor). V: no cohort numbers — citations→Lane R (R10). |
| 05 | Higher-order & multiscale | pg 6 · ~30% | [ ] | [ ] | [ ] | drafted |
| 06 | Why this matters | NEW | [ ] | [x] | [x] | ✍ W+FIGS FINAL 07-11 — signature-detector stakes + nulls caveat; 3-axes split introduce=05·instantiate=07·reuse=24 (NO enumeration here). FIGS user-approved: img1 grid-cells-across-scales `_external/grid_cells_multiscale_giocomo_moser.jpg` (Giocomo/Moser 2011, CC BY-SA 4.0); img2 nested rhythms OURS `talk/slide06_nested_rhythms.pdf` (gen fig_talk_nested_rhythms.py). Both FUNCTIONAL. R2/R3 exemplars superseded. Only Canva page build left. |
| 07 | Diffusion / LRG | pg 7 · ~40% | [ ] | [ ] | [ ] | 3-axes CORRECTED (diffusion = topology) |
| 08 | The bet | pg 8 · ~55% | [x] | [x] | [x] | ✅ LOCKED M 07-11 (revised) — HERO = OUR 3D "form mirrors form" scene `pipeline_brain_tree_3d.png` (brain network at base → cophenetic tree rising above → diffusion streams; community-colour glue). Behrens/Park EXTERNAL figures DROPPED (user: "figure from the other papers completely wrong… not the one we're looking for"). On-slide text EXPANDED (7-line block). Optional small race-bracket callback + ≟. ⚠ REUSE: same fig is the slide-11 pipeline hero → bet=teaser/vision, 11=step-by-step build (or distinct camera); don't repeat unchanged. ⚠ black bg on light deck. Fixed 07→08 numbering lag. W: text+speech final. V: no cohort numbers. |
| 09 | Imaginary coherence | pg 9 · ~30% | [ ] | [ ] | [ ] | drafted |
| 10 | Frequency bands | pg 10 · ~35% | [ ] | [ ] | [ ] | drafted |
| 11 | Pipeline (ρ^coph) | pg 11 · ~55% | [x] | [x] | [x] | ✅ LOCKED M 07-11 — 9 compile-ready LaTeX formulas added (pipeline ladder A→L̂→K̂(τ)→D(τ)→D^coph + ρ^coph + ρ_sym trace schema + MS gate), MATCHED to code (LRG guide §2, cross_phase.py, rho-sym migration report). Protocol rule 9 convention. 5-panel SEQUENCE (timeseries→imcoh→propagator→dendrogram→ρ^coph tanglegram, Pat_05 β) building via subagent → talk/pipeline_seq_1..5. Fixed 10→11 numbering lag. ⚠ 3D hero shared w/ slide 08 → here lead with the step-by-step SEQUENCE not the same still. V: no cohort numbers; formulas code-verified. |
| 12 | Tanglegram categories | NEW | [ ] | [ ] | [ ] | trace/anchor/reset/emergent; FINALIZE 4 drafts |
| 13 | A lasting trace | pg 13 · ~30% | [ ] | [ ] | [ ] | reveal A |
| 14 | Who's right? — the nulls decide | pg 12 · ~30% | [ ] | [ ] | [ ] | MERGED old 14 (who's-right) + 15 (nulls) — reveal B + controls C; 3-method fig + null triangle |
| 15 | Who survives — it must be multiscale | NEW | [ ] | [ ] | [ ] | controls D; "it must be multiscale" (was 16) |
| 16 | β → OFC — held consolidation | pg 14–15 · ~20% | [ ] | [ ] | [ ] | MERGED old 17 (β→OFC) + 18 (held-not-replayed) — WHERE/HOW/WHO + persistent-across-traced-bands |
| 17 | Encoding vs inference (setup) | pg 16 · ~20% | [ ] | [ ] | [ ] | was 19 |
| 18 | Inference persists (β only) + localization | pg 17–18 · ~30% | [ ] | [ ] | [ ] | MERGED old 20 (climax) + 21 (localization) — enc→OFC / inf→cingulate |
| 19 | Per-band taxonomy | pg 16 · ~30% | [ ] | [ ] | [ ] | β flagship; low-γ cingulate explicit (was 22) |
| 20 | Epilepsy marker | pg 19 · ~30% | [ ] | [ ] | [ ] | cuttable (was 23) |
| 21 | Take-homes | pg 20 · ~30% | [ ] | [ ] | [ ] | + significance beat (was 24) |
| 22 | Outlook & thanks | pg 20 · ~20% | [ ] | [ ] | [ ] | drafted (was 25) |

*Structure locked 2026-07-11 → 22 slides. Restructure: 14+15+16 collapsed to the nulls (14) +
survivors table (15); 17+18 → β→OFC / held-not-replayed (16); 20+21 → inference-persists +
localization (18). Canva deck (`DAHO6jbUdrY`, 20 pp) still needs: the null moved into results, a new
tanglegram slide (12), the merged "who's right + nulls" (14) + survivors table (15), a "why it
matters" slide (06). 22 slides for 25–30 min ≈ 70–80 s each — comfortable; the real risk is
over-talking, so keep on-slide text minimal + speeches tight. Optional trim if needed: slide 20 (epilepsy coda).*

---

## Figure manifest — worklist for the Figure agent (Lane F + R)

**Target:** copy/produce every "ours" figure into `data/outputs/figures/talk/`; external → `talk/_external/`.
**The figure list is SET for all 22 slides — Lane F may start NOW** (gathering does not
wait on per-slide text lock). Record the generating script for each produced figure.

Action legend: **COPY** (final PDF exists → copy into talk/) · **FINALIZE** (a `_DRAFT` → finish → talk/)
· **BUILD** (new figure from data/code) · **RENDER** (interactive `.html` → static export) ·
**LOCATE** (find or build) · **CANVA** (presenter builds in Canva — NOT a figure-agent task)
· **EXTERNAL** (from a paper — Lane R).

### ⭐ Real figure-agent work — do these first
| slide | figure | action | source |
|---|---|---|---|
| 12 | 4 tanglegrams — trace / reset / anchor / reorganized | FINALIZE | ✅ **F 07-11** → 4 PDFs in `talk/` (dropped `_DRAFT`, kept fate/pat/band/N). Copied the user-approved drafts (seeded random pick → NOT re-run, would re-randomise). gen `scripts/01_compute/figures_embedded/fig_cophenetic_tanglegram_class.py` |
| 14 | raw vs Grassmann vs cophenetic, per band ("who's right") | BUILD | ✅ **F 07-11** → **`talk/fig_whos_right_three_methods.pdf`**. NEW gen `scripts/01_compute/figures_embedded/fig_whos_right_three_methods.py` — 3×6 method×band verdict matrix (disc size=−log10 p_MS, filled=clears null). ⚠ HONEST DATA: at the MS gate coph lights **α+β** (not β-only); β = the only band lit by all three ⇒ setup for "nulls decide". See Blockers. |
| 18 | inference "reasoning bloom" (β wedge glows alone) | FINALIZE | ✅ **F 07-11** → **`talk/fig_reasoning_bloom_inference_raw.pdf`**. gen `scripts/01_compute/figures_embedded/fig_reasoning_bloom.py --func T_infspec` |
| 19 | low-γ cingulate ENCODING — whole-brain zero masks focal | LOCATE/BUILD | ✅ **F 07-11** LOCATE→copy (no build) → **`talk/fig_arc_e_lowgamma_cingulate.pdf`** (already ρ_sym-current: 8/8, q=0.035, inference null). gen `scripts/01_compute/figures_embedded/fig_arc_e_lowgamma_cingulate.py` |
| 07 | diffusion τ-sweep video (ρ(τ) morphing on FC) | BUILD/RENDER | ✅ **F 07-11** assets fresh → **`talk/fig_fc_diffusion_tausweep.mp4`** (mp4 only; `.mov`+`_snapshots.pdf` DROPPED per user 07-11 — mov unreadable/unrequested). gen `scripts/01_compute/figures_embedded/fig_fc_diffusion_tausweep.py` |
| 04 | 3D cortex connectome | RENDER | ✅ **F 07-11** → **`talk/fig_brain_connectome_Pat_05_n10_rsPre_beta.{html,png}`** (Pat_05 β rsPre, n=10). static via NEW `scripts/07_figures/plotly_html_to_static.py --drop-title --drop-text-traces --zoom 0.62 --colorbar-title coupling`; html gen `scripts/07_figures/gen_brain_connectome_imcoh.py` · ⚠ NOW BACKUP for slide 04 (M 07-11) — slide 04 uses the FC-matrix mosaic instead; connectome held for later. |
| 16 | OFC pair-glow 3D | RENDER | ✅ **F 07-11** → **`talk/fig_beta_ofc_pairglow_3d.{html,png}`** (png already existed). gen `scripts/01_compute/figures_embedded/fig_beta_ofc_pairglow_3d.py` |
| 16 | laterality brain 3D (patients tinted by β net) | RENDER | ✅ **F 07-11** → **`talk/fig_trace_laterality_brain_beta.{html,png}`**. static via `scripts/07_figures/plotly_html_to_static.py`; gen `scripts/01_compute/figures_embedded/fig_trace_laterality_brain.py` |
| 02 | TI relation chain — 4 disc-golf player photo-boxes + ">" | BUILD | ✅ **F 07-11** → **`talk/ti_relation_chain.pdf`** (4 boxes FB·AG·DG·SM, faint initials + names, gold ">" = shown premises). NEW gen `scripts/01_compute/figures_embedded/fig_ti_relation_chain.py` |
| 02 | TI race bracket — 3 round-states (round0/1/2) + lightning | BUILD | ✅ **F 07-11** → **`talk/ti_race_bracket_round{0,1,2}.pdf`** (+ `ti_race_bracket.pdf` = round2). round0 both matchups unresolved "who wins these matches?" · round1 both semis→FB/AG amber lightning + final appears "who wins the final?" · round2 final gold lightning + champ FB ★ "seen (encoded): FB > AG". gen `scripts/01_compute/figures_embedded/fig_ti_race_bracket_demo.py` |
| 03 | 4-phase protocol timeline (rest_pre → task_learn → task_test → rest_post) | BUILD | ✅ **F 07-11** → **`talk/ti_protocol_timeline.pdf`** (time-arrow · 4 rounded blocks cool-rest/warm-task · chevron ordering · "transitive-inference task" bracket over task_learn+task_test · descriptors baseline/learn-adjacent-premises/test-non-adjacent-inferences/post-task · time axis). LIGHT slide / dark ink. **Durations OFF** but wire-ready: `DURATIONS` dict in the script — fill confirmed minutes (rest ≈ 10 min; task_learn/task_test still need Lane R's Ricci 2023 pin or the collaborator — cohort is heterogeneous, don't invent) and they render as "~X min" inside each block. NEW gen `scripts/01_compute/figures_embedded/fig_ti_protocol_timeline.py` |
| 03 | raw sEEG traces (wide banner) | BUILD | ✅ **F 07-11** → **`talk/seeg_traces_stacked.pdf`** (Pat_05 rest_pre · **wide horizontal banner** · 5 channels = 3 random non-SOZ dark + **2 random SOZ RED** · **8 s window** (long, many timepoints) · common gain · **continuation arrow off each trace's right end** ("recording goes on") · 1 s + amplitude scale bars · "2048 Hz" note). Vector PDF, light-slide/dark-ink. NEW gen `scripts/07_figures/gen_seeg_traces.py` (args `--patient/--phase/--start/--dur/--n-ctx/--n-soz/--seed/--uv`; amplitude a.u. — `--uv` if export confirmed µV) |
| 08 · 11 · 07 | pipeline HERO — network + brain + cophenetic tree in ONE 3D scene | BUILD | ✅ **F 07-11** → **`talk/pipeline_brain_tree_3d.{html,png}`** (Pat_05 β rest_pre). One plotly scene: translucent pial + implant spheres + top-120 **\|ImCoh\| network arcs** at the BOTTOM; **cophenetic dendrogram growing UP** above the head; **bundled "diffusion" root-streams** contact→leaf. GLUE = 6 mesoscale communities (distance cut, `--cut-frac 0.68`) tint spheres+roots+clade the same hue, rest grey (dense imcoh_abs tree is chained → k-way partition degenerate, so top-6 clades coloured + grey background). **Child-swap orders clades L→R by brain-x** so streams rise vertically, not crossing. Tree is a legible standing billboard (full-3D-embedded leaves = unreadable tangle). NEW gen `scripts/07_figures/gen_pipeline_brain_tree_3d.py` (`--patient/--band/--phase/--n-comms/--cut-frac/--min-size/--n-edges/--roots {comms,all,none}` + camera `--eye/--up/--center`). **Camera default = user's "hierarchy seen from BELOW" view** (camera lateral −x, slightly behind +y so orange sits left, and BELOW the tree −z looking UP → hierarchy underside floats up-left while the whole brain reads as a lateral silhouette, unobscured; `eye=-1.50,1.00,-0.45 center=0,0,0.28`) — baked as default so the HTML opens on it too. (Superseded from-ABOVE "idea rising" `eye=-0.72,1.36,0.67`; user wanted the camera dropped below to look UP at the tree, 07-11.) ✅ RESOLVED (M 07-11): user WANTED it on the BET (08) — now the slide-08 HERO ("form mirrors form"). Still ALSO the pipeline (11) hero → bet=teaser/vision, 11=step-by-step build (or distinct cameras); do NOT let the same still sit unchanged on both. |
| 11 | pipeline SEQUENCE — 5 single panels (timeseries → \|ImCoh\| → distance D(τ)=1/K → dendrogram → ρ^coph tanglegram) | BUILD | ✅ **M 07-11** → **`talk/pipeline_seq_{1..5}_*.png`** (high-res transparent PNG, user opt-in for Canva import 07-11; PDFs removed). Pat_05 β; compose L→R with arrows in Canva. Panel 3 = communication distance D(τ)=1/K (NOT K — at τ_min K≈FC; D is the distinct object UPGMA runs on). Panel 5 = two trees (rest_pre vs rest_post) with a RETAINED clade highlighted (5 parallel amber ribbons = kept its place; grey crossing = reorganised); subset ρ^coph=0.59 annotated, full-matrix ρ^coph=0.42 printed honestly beneath. House-style (use_lrg_style, generic FCᵢⱼ / D_ij(τ) cbars, no suptitle). NEW gen `scripts/07_figures/gen_pipeline_sequence.py`. |
| 06 | nested rhythms (θ–γ coupling) — "cognition across time-scales" | BUILD | ✅ **07-11** → **`talk/slide06_nested_rhythms.pdf`** (+`.png` for Canva). Synthetic cross-frequency-coupling cartoon; band-palette (θ warm / γ cool); zoom inset nests fast-in-slow. Pairs with the external grid-cells fig (space-scales). NEW gen `scripts/01_compute/figures_embedded/fig_talk_nested_rhythms.py` |

**TI cold-open build notes (slide 02 — disc golf; replaces the old static bracket AND the CANVA player-chain):**
Sport = DISC GOLF (a 1-v-1 sport; swapped from Ultimate, which is a team sport and breaks the
individual-ranking analogy). Two vector-PDF deliverables (a relation chain + a 3-state bracket).
Transparent background, DARK ink on a LIGHT slide (the deck is light), `use_lrg_style()`. Every box is a PHOTO PLACEHOLDER
(the four professor headshots are dropped in Canva). Palette (as Lane F built it): DARK `#1a1a1a`
(lines/borders/text), GRAY `#6f757c` (faint initial hints), AMBER `#e08d00` (inferred matches),
GOLD `#f2c14e` (seen final / champion). Order of
the four guests: FB > AG > DG > SM (Battiston, Gabrielli, Garlaschelli, Meloni).

These figures drive a 6-sub-slide reveal the presenter assembles in Canva: beat 1 = disc-golf photo →
beat 2 = relation chain → beat 3 = bracket round0 → beat 4 = bracket round1 → beat 5 = bracket round2
(champion) → beat 6 = closing text. Figures the agent must produce = the chain + the 3 bracket rounds.

- **`ti_relation_chain.pdf`** — the four players ranked by ability: 4 square photo-boxes in a row,
  faint initial hints FB · AG · DG · SM, full name under each, a large ">" between adjacent boxes.
  Each ">" is one SHOWN adjacent result (FB>AG, AG>DG, DG>SM); the non-adjacent matchups are never
  shown → must be inferred. Shown at beat 3 ("here's the ranking"), then vanishes as the bracket appears.
- **`ti_race_bracket_round{0,1,2}.pdf`** (+ keep `ti_race_bracket.pdf` = round2) — knockout hierarchy
  over leaves FB, DG, AG, SM. First-round matchups are the two INFERRED (non-adjacent) pairs: SF1 =
  FB vs DG, SF2 = AG vs SM. The FINAL FB vs AG is the one SEEN adjacent premise → champion FB. THREE
  round-states (round 1 resolves BOTH semifinals together):
    - round0 = both first-round matchups drawn, NO winners ("who wins these matches?"); matchup
      branches neutral/ghosted; tag each semifinal "inferred".
    - round1 = both semifinals resolved — FB and AG advance, their semifinal branches strike in as
      AMBER lightning; the final matchup FB vs AG now appears, still unresolved ("who wins the final?").
    - round2 = final resolved — champion FB (Battiston) on top; final branch strikes in as GOLD
      lightning + a ★; tag the final "seen: FB > AG".
  Ghost the latent skeleton faint so the shape is anticipated; winner nodes show the initials and are
  photo placeholders (a headshot can be dropped on top in Canva).
  Extend the existing generator `scripts/01_compute/figures_embedded/fig_ti_race_bracket_demo.py`.

### COPY — gather existing PDFs into `talk/` (mechanical)
| slide | figure | source | done |
|---|---|---|:--:|
| 02, 08 | tournament race bracket | ⚠ SUPERSEDED — now a 3-round BUILD (round0/1/2 + lightning, disc golf), see ⭐ table + TI cold-open build notes | [ ] |
| 03 | cohort implants overlay | `talk/cohort_implants_overlay_3d.{png,html}` | [x] **REBUILT 3D** — F 07-11 · flat glass brain → pooled 3D pial shell, contacts as mm-space **spheres**, one colour/patient, **3 angular views** (left·superior·right), 1169 contacts n=10. NEW gen `scripts/07_figures/gen_cohort_implants_3d.py` (uses `visuals.brain3d`); interactive html + static png. Old `implant_in_brain/cohort_implants_overlay.pdf` superseded |
| 04 | ~~FC matrix ↔ network~~ → **FC matrix MOSAIC** | ⚠ SUPERSEDED (M 07-11): sfdp matrix↔net too sparse/minimal. NOW → `data/outputs/figures/talk/n4_delta_imcoh_abs_generic_log_per_row.pdf` (4 pat × 4 phase, generic FCᵢⱼ cbar, no imCoh/band on-figure). gen `.agents/guides/05_plotting/fc_templates/mosaic_patient_phase.py --patients Pat_02,Pat_05,Pat_08,Pat_13 --band delta --fc-method imcoh_abs --tick-labels generic --colorbar-label '$\mathrm{FC}_{ij}$'` (added `colorbar_label` override to `plot_fc_adjacency_grid` + `use_lrg_style()` to the template) | [x] — M 07-11 |
| 05→19 | multiscale-brain result — RE-HOMED off 05 (don't spoil results on the gap slide) | `data/reports/rho_sym_band_map/fig_rho_sym_band_map.pdf` | [x] — F 07-11 · slide-6 is NOT a montage (grid cells + nested rhythms) → held for slide-18 (old 19, renumbered 2026-07-14) |
| 07 | dendrogram cut snapshots | `data/outputs/figures/lrg_diffusion_zoom/lrg_dendrogram_cut_snapshots.pdf` | [x] — F 07-11 |
| ~~08~~ (methods/results) | circular dendrogram network | `data/outputs/figures/network_templates/circular_dendrogram_network/Pat_05_beta_rest_pre_imcoh_abs_K7_b0.92_g2_with_dendro.pdf` · ⚠ DROPPED from slide 08 (M 07-11) — bet uses Behrens cognitive-map panel; retain for methods/results. NB has a busy suptitle (Pat_05·β·imcoh params) → needs a clean method-neutral regen if reused on a talk slide. | [x] — F 07-11 |
| 09 | spectral dist ImCoh vs MSC | `data/outputs/figures/section2/fig_C/fig_C2_spectral_distribution_{imcoh,msc}_rsPre.pdf` | [x] — F 07-11 |
| 09 | ~~ImCoh band average~~ → ImCoh spectral distribution | ⚠ user 07-11: `imcoh_band_average.pdf` was the WRONG kind — REMOVED; replaced by `fig_C2_spectral_distribution_imcoh_{rsPre,rsPost}.pdf` (rsPre+rsPost) | [x] — F 07-11 |
| 09, 11 | MSC vs ImCoh dendrograms | `data/outputs/figures/section3/fig_H/fig_H1_dendrograms_MSC_vs_ImCoh_3patients.pdf` | [x] — F 07-11 |
| 10 | per-band FC matrices | `data/outputs/figures/fc_templates/row_per_band/Pat_05_rest_pre_imcoh_abs_chnames_log_shared.pdf` | [x] — F 07-11 |
| 13 | raw vs multiscale | `data/preprint/figures/results_section1/fig_trace_c_raw_vs_multiscale.pdf` | [x] — F 07-11 |
| 13 | per-band cophenetic forest | `data/preprint/figures/results_section1/fig_trace_a_band_forest.pdf` | [x] — F 07-11 |
| 14 | cophenetic null triangle | `data/preprint/figures/all_bands/fig_bands_null_triangle_coph.pdf` | [x] — F 07-11 |
| 15 | coph ↔ Grassmann confirmation matrix | `data/outputs/figures/talk/fig_coph_grassmann_confirmation_matrix.pdf` (already in talk/) | [x] |
| 16 | OFC localization | `data/preprint/figures/results_section1/fig_trace_b_ofc_localization.pdf` | [x] — F 07-11 |
| 16 | β anatomy brain | `data/preprint/figures/beta/anatomy/fig_beta_anatomy_brain.pdf` | [x] — F 07-11 |
| 16 | reinstatement | `data/preprint/figures/results_section1/fig_trace_e_reinstatement.pdf` | [x] — F 07-11 |
| 16 | per-patient slopes | `data/outputs/figures/section_5_lrg_trace/headline/per_patient_slopes.pdf` | [x] — F 07-11 |
| 17 | decomposition forest | `data/preprint/figures/results_section2/fig_arc_a_decomposition_forest.pdf` | [x] — F 07-11 |
| 18 | encoding-vs-inference dissociation | `data/preprint/figures/results_section1/fig_encoding_inference_dissociation.pdf` | [x] — F 07-11 |
| 18 | inference → cingulate | `data/preprint/figures/results_section2/fig_arc_f_inference_cingulate.pdf` | [x] — F 07-11 |
| 19, 21 | per-band phenomenology | `data/reports/per_band_phenomenology/fig_per_band_phenomenology.pdf` | [x] — F 07-11 |
| 20 | SOZ divergence (β spares / α recruits) | `data/preprint/figures/results_section1/fig_trace_f_soz_divergence.pdf` | [x] — F 07-11 |
| 20 | SOZ diffusion community | `data/preprint/figures/results_section1/fig_soz_diffusion_community.pdf` | [x] — F 07-11 |
| 20 | epi marker + calibrated detector | `data/reports/results_section3/fig_epi_{a_relational_marker,b_calibrated_detector}.pdf` | [x] — F 07-11 |

### CANVA — presenter builds these schematics (figure agent: SKIP)
04 structure/function schematic ·
05 hypergraph sketch · 07 three-axes→one-lens (bands=time · sEEG=space · diffusion=topology — trio-completion visual; moved off 06 per the 3-axes split) ·
08 trace/anchor/reset/emergent strip (or reuse slide-12 tanglegrams) · 11 pipeline panel ·
14 null schematic (matched-strength + drift) · 15 survivors table (method × null) ·
17 e/f/p decomposition schematic · 18 OFC+cingulate on one brain (optional).

### EXTERNAL — worklist for the Reference agent (Lane R → `talk/_external/`)
Fetch each into `talk/_external/`; for each record the **full citation + figure number +
licence/attribution** (talks may show a cited figure with attribution — note the source on-slide).
For "propose" rows the agent suggests 2–3 candidates and the **user picks** — don't fetch blindly.

| # | slide | figure wanted | source | status |
|---|-------|---------------|--------|--------|
| R1 | 05, 06 | three-axis multiscale-brain schematic (time / space / topology) | Betzel & Bassett 2017, NeuroImage 160 | [x] ✅ FETCHED — `betzel_bassett_2017_fig1_multiscale_brain.png` (from arXiv:1608.08828) — R 07-11 · ⚠ nested 3rd-party sub-images, talk-only |
| R2 | ~~06~~ | clinical seizure-network exemplar | Khambhati et al. 2015 | [x] ✅ FETCHED — `khambhati2015_seizure_network_states_fig3.png` (Fig 3; PLoS Comput Biol 11:e1004608; CC BY 4.0) · ⚠ UNUSED on slide 06 after 07-11 refocus to functional cognition (grid cells + nested rhythms) — kept as backup — R 07-11 |
| R3 | ~~06~~ | ~~cognitive consolidation exemplar (Diekelmann)~~ | ⏹ SUPERSEDED | ⏹ SUPERSEDED 07-11 — slide 06 refocused to FUNCTIONAL cognition: img1 grid cells across scales (R14), img2 nested rhythms (ours). Diekelmann dropped; no PDF/redraw needed. |
| R4 | 02 / 08 | (optional) cognitive-map / grid-code schematic | Behrens 2018 (your pick) | [x] ✅ FETCHED — `behrens2018_fig6_transitive_inference.png` (Neuron 100:490; Fig 6 "TI"; via bioRxiv 365593; 🟡 © all-rights-reserved → talk-only; +panel-A alt) — R 07-11 · ⚠ slide 08 UN-picked it (M 07-11): user rejected external figures for the bet ("completely wrong") → bet now uses OUR 3D hero. Behrens/panel-A now only a candidate for slide 02 if wanted; else unused for the deck. |
| R5 | 02 | **disc-golf** photo (1-v-1 sport — SWAP from Ultimate) | stock / open-licence | [x] ✅ FETCHED — `disc_golf_scottish_open_2007_putt.jpg` (Conor Lawless, CC BY 2.0; Scottish Open; +2 alts) — R 07-11 · Ultimate retired |
| R6 | 03 | sEEG shaft / implant photo | clinical atlas / open-licence | [x] ✅ FETCHED — `seeg_electrode_implantation_herff2020_fig1.png` (Herff'20, CC BY 4.0; +panels A/C) — R 07-11 |
| R7 | 03 | pin exact citation: acquisition-protocol paper (collaborator lead) | Ricci et al. 2023 (Bambino Gesù group) | [x] 📌 PINNED — Ricci L. et al. 2023, Clin Neurophysiol 150:40–48, DOI 10.1016/j.clinph.2023.03.006 (Tamilia/Mercier/Specchio/de Palma) — R 07-11 · ⚠ title = PAC biomarker, not literally "acquisition protocol" — confirm this is the one. Paywalled → fig talk-only. |
| R8 | 03 / methods | pin exact citation: sEEG methods (collaborator lead) | Wang et al. 2023 | [x] 📌 CANDIDATE — Wang C. et al. 2023 = BrainBERT, ICLR 2023 / arXiv:2302.14367 — R 07-11 · ⚠ it's a representation-learning MODEL, not a referencing/acquisition-methods paper — confirm vs an iEEG-referencing paper if that's what you meant. |
| R9 | 03 | (Reber 2016 REMOVED by user) — Di Bello/Ferraina already confirmed | Di Bello et al. 2024 | [x] ✅ Commun Biol 7(1):1715, DOI 10.1038/s42003-024-07418-5 — CONFIRMED (macaque). No further action. |
| R10 | 04 | pin locators/DOIs for the 4 net-neuro RESULT citations now on slide 04 | Bassett&Bullmore'06 (small-world), vandenHeuvel&Sporns'11 (rich-club), Sporns&Betzel'16 (modules), Honey'09 (fn≠struct) | ⏳ FIND — all canonical, low-risk; pin DOI+locators. Author/year/journal already in slide 04 REFERENCES; conductor confident but flag any mismatch. Bullmore&Sporns'09 (umbrella) already correct. |
| R11 | 09 | volume-conduction / field-spread figure (motivates ImCoh) | Bastos & Schoffelen 2016, Front Syst Neurosci 9:175 | [x] ✅ EXTRACTED M 07-11 — `_external/bastos2016_fig6_field_spread_imcoh.png` (**Fig. 6** field-spread, 🟢 CC BY 4.0) + optional ImCoh exemplars `nolte2004_fig12_imcoh_beta_headmap.png` (β head-map) & `nolte2004_fig11_significant_imcoh.png` (🟡 © Elsevier, talk-only). All in CITATIONS.md (R11–R13). Lane R: verify Nolte DOI 10.1016/j.clinph.2004.04.029. |
| R14 | 06 | grid cells across scales (functional cognition × multiscale) — image 1 | Giocomo, Moser & Moser 2011, Neuron 71:589 | [x] ✅ FETCHED 07-11 — `_external/grid_cells_multiscale_giocomo_moser.jpg` (Wikimedia; CC BY-SA 4.0; 375×575 small → keep ≤ ½-slide). CITATIONS.md R14. Image 2 (nested rhythms) is OURS → see ⭐ table. |

**Lane R record →** `data/outputs/figures/talk/_external/CITATIONS.md` (full citations, licences,
ready-to-paste on-slide attribution lines, swap options, + the 3 PROPOSED shortlists).

Launch line: *"Read `.agents/talk/PLAN.md` §Figure manifest → EXTERNAL. You are Lane R. Fetch the
CLEAR/stock rows into `data/outputs/figures/talk/_external/` with full citation + attribution; for
PROPOSE/OPTIONAL rows, suggest 2–3 candidates and wait for the user to pick."*

**COPY note:** copy, don't move (originals stay where the manuscript expects them); record each
produced figure's generating script; `RENDER` = keep the `.html` AND a high-res static export in `talk/`.

---

## Verify lane (V) — load-bearing numbers to re-confirm fresh

Rule: **recompute / cross-check from the latest artifacts — never lift a number
from prior prose or a cached CSV** (two agents once disagreed 6/10 vs 10/10).
Confirm each against the ρ_sym-era SSOT + newest audit. Flag mismatches in Blockers.

| Slide | Numbers to verify |
|-------|-------------------|
| 14 · 15 | who-survives: raw ~half killed by MS + fails drift; Grassmann survives MS (k=27..55, 29=manuscript) but drift-vulnerable p≈0.116 (audit_169); coph β 23.7× MS + drift-clean C2 0.0137. Gate: β p=0.032/16.6×, α p=0.024/13.7×, γ/δ marginal p≈0.080 |
| 16 | β→OFC MS R=1000 BH q 0.009–0.013, all 4 conditions, LOO+shaft; ~59% within-OFC pro, 5/10 implanted |
| 16 | reinstatement 10/10 p=0.001, LOO 0.002; carriers L+0.34 / R+0.05; trait ρ=0.95 (audit_148) |
| 18 | β infspec +0.091 vs null +0.005 p=0.0098 (LO-P15 0.020, 6/10); τ≈2.6 p=0.0049; duration ρ=+0.10 p=0.78; other bands next-best p≈0.25; convergence raw 0.010 / eigcent 0.002 / closeness 0.024 |
| 18 | coph inference→cingulate MS BH q 0.03/0.04, LOO max 0.049/0.032 (audit_171); raw placeless; length-null p≈0.124 (contested/OPEN); low-γ encoding→cingulate 8/8 q 0.035 |
| 20 | co-diffusion AUC δ0.80 / γl0.74 / β0.69; off-shaft δ LOSO 0.72 8/10; detector AUC 0.81, p@5 ≈60%, 9/10 |
| 12 · 21 | 12 tanglegrams = illustrative examples (no cohort stat); 21 synthesis restates the above |

---

## Output & sync locations (decisions — veto if wrong)

- **Figures gathered →** `data/outputs/figures/talk/` (rides the travel-bundle rsync).
- **External figures →** `data/outputs/figures/talk/_external/`.
- **Interactive 3D →** keep `.html` + a static export side-by-side in `talk/`.
- **Generating code →** stays in `scripts/`/`src/`; **paths recorded per figure** in the
  manifest (rides git). No code copied into `.agents/`.
- **Naming:** keep the source filename when copying into `talk/` (traceable). Rename
  only if the conductor asks.

---

## Blockers / open questions

- **Structure locked at 22 slides (2026-07-11).** Content drafted; per-slide review/lock with the user still pending.
- **Length: 22 slides / 25–30 min ≈ 60–70 s each — fits.** Risk is over-talking → keep on-slide text minimal, speeches tight. Optional trims only if needed: slide 20 (epi), then merge 13 + 14.
- **Slide 18 duration control OPEN.** Inference→cingulate localizes under MS (q 0.03/0.04) but a length-matched null demotes it to p≈0.124; that null is likely invalid for a conditional. Present as a strong lead, not closed.
- **Slide 14 (who's right + nulls) — α/β honesty flag RESOLVED (M 07-11).** `talk/fig_whos_right_three_methods.pdf`. The merged slide 14 now states the cophenetic read-out lights **α (p=0.024) AND β (p=0.032)** at the matched-strength gate — NOT "β-selective"; β-selectivity is the *downstream* result (drift on 15 + localization on 16/18), not the raw coph verdict. Figure shows raw≈diffuse / Grassmann={δ,β,γ_low} / coph={α,β}, β = the only band lit by all three — the "who's right?" hook the nulls then resolve to β. The Careful line in `14_the-nulls.md` carries this guardrail.
- **Lane F gathering DONE (07-11):** all 8 original ⭐ rows + all 24 COPY rows → `data/outputs/figures/talk/` (41 files). 3 user corrections applied: dropped `fig_fc_diffusion_tausweep.mov` (unreadable ProRes/unrequested) + `_snapshots.pdf` (slide 07 keeps mp4 only); slide 09 `imcoh_band_average.pdf` was the WRONG kind → replaced with `fig_C2_spectral_distribution_imcoh_{rsPre,rsPost}.pdf`. RENDER statics reproducible via new `scripts/07_figures/plotly_html_to_static.py`.
- **Slide-02 TI cold-open BUILT (F 07-11):** `talk/ti_relation_chain.pdf` + `talk/ti_race_bracket_round{0,1,2}.pdf` (+ `ti_race_bracket.pdf`=round2). Disc-golf framing, both semis inferred (amber lightning), final seen (gold), champion FB ★, "seen (encoded): FB > AG". Photo-placeholder boxes/chips (drop headshots in Canva). gens: `fig_ti_relation_chain.py`, `fig_ti_race_bracket_demo.py`. (Old `stage{0-3}` naming removed — spec moved to 3 round-states.)
- Decisions in §Output above are the conductor's defaults — flag here if any is wrong.

---

## Next action (this chat / conductor)

22 slides drafted in plain format + figure list SET. Lane F can start gathering now
(⭐ work first: tanglegrams 12, 3-method panel 14, inference bloom 18, low-γ cingulate 19,
+ the four RENDER 3D/video). Conductor: begin per-slide review/lock with the user; as each
slide locks, flip its `Locked` box.
