> **⚠ QUANTITATIVE_STALE (2026-04-15)** — numbers in this report were
> computed under the pre-reset ImCoh mislabelling (stored `|ImCoh|²`
> under the name `imcoh`). Qualitative conclusions survive (rankings
> preserved under sqrt); regenerate absolute values against
> `fc_method="imcoh_abs"` before quoting. See
> `~/.claude/projects/.../memory/imcoh_taxonomy.md`.

# Section 2 Figures — Handoff to Writing Agent

**Target:** sections 2.1–2.3 of `notes_imcoh.tex` (the MSC→ImCoh rewrite of the
original `notes.tex` Section 2, which motivated the estimator switch at the
FC‑estimation level).

**Action requested:** read this document side‑by‑side with your original
`WRITING_AGENT_BRIEFING.md`, and **select ≤ 20 items** (figure files + report
`.md` files) that you actually need to draft Sections 2.1–2.3. Return the
selected list; we will then produce a short quantitative description per
selected figure.

---

## 1. What changed vs. your initial plan

Your initial briefing (`WRITING_AGENT_BRIEFING.md` §8 "Figures to reproduce
with ImCoh") proposed a generic list:

1. Unanimity maps (H1, H2a, H2b, H3)
2. Scalar contrast radar charts / heatmaps
3. VI(k) profiles by band
4. Adjacency matrix heatmaps (MSC vs ImCoh)
5. Entropy/susceptibility curves S(τ), C(τ)
6. Brain connectome plots (community structure)

This is the hypothesis‑testing layer (Sections 3+). **Section 2 sits
upstream** of that: it must *prove* the estimator switch at the FC level,
before any LRG hypothesis tests. We therefore produced a different, more
granular set of figures covering three logical blocks:

- **2.1 — same‑probe bias exists in MSC.** Adjacency heatmaps (A/B), weight
  distributions (C), enrichment vs. scale (D).
- **2.2 — ImCoh removes it (at the cost of sparsity/noise).** Same figures
  with MSC vs ImCoh contrasted; spectral‑distribution overlays (C2).
- **2.3 — graph‑level consequences.** Force‑directed networks on the LRG
  ultrametric layout (A‑network, E2), brain connectome (E1), nodewise graph
  metrics (F), Laplacian spectrum (G).

### Concrete deviations from your plan

| Your item | What we produced instead | Why |
|---|---|---|
| "Adjacency matrix heatmaps (MSC vs ImCoh side‑by‑side)" | **fig_A** per patient (phases × bands grid, per‑band vmax) **and** **fig_B** all patients at fixed (band, phase) — MSC vs ImCoh rows | The single side‑by‑side view loses either patient variety or band/phase variety; we split into two complementary views |
| "Brain connectome plots" | **fig_E1** nilearn glass brain (only `PATIENTS_WITH_COORDS = [Pat_02, Pat_03, Pat_05]`) **and** **fig_E2** force‑directed on LRG ultrametric for 3 patients × 2 methods | Brain coordinates exist only for 3 patients; E2 generalises to all patients and exposes modular structure MSC cannot produce |
| — (not in your plan) | **fig_C1** per‑band same/cross‑probe weight histograms, **fig_C2** Welch spectral‑distribution median + IQR + 5–95%, **fig_C3** KDE overlay across patients | These quantify the bias at the raw‑weight level, which is what Section 2.1 needs before any LRG talk |
| — (not in your plan) | **fig_D1** paired‑scatter w/ marginals (MSC vs ImCoh probe‑bias ratios per patient) and **fig_D2** bias‑reduction 3‑panel (heatmap + curves + strips) | Gives a single‑figure quantitative "ImCoh removes bias" summary |
| — (not in your plan) | **fig_A_network** (same (phases × bands) grid as fig_A but drawn as networks) | Shows that when the bias is removed, modular structure actually shows up in the graph |
| — (not in your plan) | **fig_F** 4‑panel metrics‑on‑graph (strength, clustering, participation, degree) per (patient, band, phase, method) | Demonstrates that nodewise metrics computed on MSC are driven by probe identity, on ImCoh by anatomy |
| — (not in your plan) | **fig_G** Laplacian eigenvalue spectrum, 4 phases overlaid, MSC vs ImCoh | Sanity check — spectral signature consistent with community collapse in MSC |
| **Fig_01** (the old `section2/fig_01/` set — heatmap + weight dists + enrichment vs scale, **beta only**) | **Kept**, but partially subsumed by the newer fig_D family which is multi‑band, multi‑phase, and has improved visual design | You may want these as the simplest/cleanest summary for a single‑page figure |

### Finalized visual recipe (applied everywhere)

All force‑directed network figures (A‑network, E2, F) use a **single recipe**
baked into `src/lrg_eegfc/config/const.py`:

- **Layout:** Kamada–Kawai with pairwise distances from the LRG ultrametric
  matrix (not `1/weight`). ImCoh's near‑uniform raw edges produce hairballs
  under any direct spring/KK; the LRG ultrametric carries the modular signal.
- **Same‑probe edges** drawn in the shaft's `tab20` color at full alpha.
- **Cross‑probe edges** drawn in black with rank‑based alpha fade.
- **Rank‑based γ edge scaling**: `t = ((rank+1)/N)**γ` with γ = 30 for ImCoh,
  γ = 2 for MSC. Only top‑ranked edges appear thick; weak edges stay
  hair‑thin. This is also what keeps the PDFs small.

Changing this recipe propagates to **every** Section 2 network figure — no
script edits needed.

---

## 2. Available figures (inventory)

All paths relative to `data/outputs/figures/section2/`. Phases covered:
`rsPre`, `taskLearn`, `rsPost`. Bands covered: `alpha`, `beta` (occasional
`theta`).

### fig_01/ — legacy single‑page probe‑bias summary (beta‑only)

- `fig_01A_probe_enrichment_heatmap.pdf` — MSC vs ImCoh enrichment (5 patients × scales)
- `fig_01B_weight_distributions_beta.pdf` — same/cross‑probe weight dists
- `fig_01C_enrichment_vs_scale_beta.pdf` — curves MSC vs ImCoh vs scale

### fig_A/ — adjacency & network, **per patient**, phases × bands grid

Heatmap view (6 files): `fig_A_{msc,imcoh}_{Pat_02,Pat_05,Pat_08}.pdf`
Network view (6 files): `fig_A_network_{msc,imcoh}_{Pat_02,Pat_05,Pat_08}.pdf`

### fig_B/ — all 5 patients, fixed (band, phase), MSC top row / ImCoh bottom

Band × phase combinations available (6 files):
`fig_B_{alpha,beta}_{rsPre,taskLearn,rsPost}_MSC_vs_ImCoh.pdf`

### fig_C/ — weight distributions and spectral profiles

- **fig_C1** same/cross‑probe weight histograms, per patient × phase ×
  method. 18 files: `fig_C1{a=msc,b=imcoh}_weight_dist_*_{Pat_02,Pat_05,Pat_08}_{rsPre,taskLearn,rsPost}.pdf`
- **fig_C2 spectral example** (single patient Welch CSD examples): 9 files,
  `fig_C2_spectral_example_{Pat_02,Pat_05,Pat_08}_{rsPre,taskLearn,rsPost}.pdf`
- **fig_C2 spectral distribution** (median / IQR / 5‑95 percentile across
  edges, per band): 6 files,
  `fig_C2_spectral_distribution_{msc,imcoh}_{rsPre,taskLearn,rsPost}.pdf`
- **fig_C3** KDE overlay all patients: 2 files (`alpha`, `beta`)

### fig_D/ — probe‑bias quantification (MSC vs ImCoh)

- **fig_D1 paired‑dots** (scatter + marginals, MSC vs ImCoh per patient): 3
  files (`rsPre`, `taskLearn`, `rsPost`)
- **fig_D1 enrichment‑heatmap** (single redesigned summary): 1 file
- **fig_D2 bias‑reduction** (3‑panel: log2 ratio heatmap + mean±SD curves +
  per‑scale strips): 4 files (`alpha_rsPre`, `alpha_taskLearn`, `beta_rsPre`,
  `beta_taskLearn`)

### fig_E/ — spatial network views

- **fig_E1 brain connectome (nilearn glass brain):** 6 files,
  `fig_E1_brain_connectome_{alpha,beta}_{rsPre,taskLearn,rsPost}.pdf` (only
  Pat_02, Pat_03, Pat_05 — others lack MNI coords)
- **fig_E2 force‑directed 3×2 grid** (3 patients × MSC/ImCoh): 6 files,
  `fig_E2_spring_Pat_02_Pat_05_Pat_08_{alpha,beta}_{rsPre,taskLearn,rsPost}.pdf`

### fig_F/ — nodewise graph metrics on the layout

4‑panel (strength, clustering, participation, degree) per
(patient × band × phase × method). 36 files:
`fig_F_{msc,imcoh}_{Pat_02,Pat_05,Pat_08}_{alpha,beta}_{rsPre,taskLearn,rsPost}.pdf`

### fig_G/ — Laplacian eigenvalue spectrum

2×6 grid, 4 phases overlaid, MSC vs ImCoh. 3 files (one per patient):
`fig_G1_eigenvalue_spectrum_{Pat_02,Pat_05,Pat_08}_MSC_vs_ImCoh.pdf`

### fig_helper/ — diagnostics, **not** for the paper

Layout sweeps, γ/k sweeps, layout galleries. Ignore unless you want a
diagnostic appendix.

---

## 3. Available supporting reports (`.agents/reports/`)

- `WRITING_AGENT_BRIEFING.md` — your initial brief (MSC→ImCoh narrative,
  tables, numbers). **Read this first; this handoff complements it.**
- `IMCOH_PROCESS_REPORT.md` — full discovery journey with technical details.
- `IMCOH_RESULTS_FOR_WRITING.md` — detailed MSC vs ImCoh comparison tables.
- `IMCOH_VERIFICATION_RESULTS.md` — verified numbers / cell counts.
- `IMCOH_PAT02_AND_CONTROLS.md` — Pat_02 and Pat_03 as controls.
- `IMCOH_GAP_ANALYSIS.md` — H2a beta k‑level gap analysis.
- `EPILEPTIC_IMCOH_FINAL.md` — paper‑ready framing of epileptic‑node clustering.

Method/guide references also available: `.agents/guides/02_methods/IMCOH_GUIDE.md`,
`PROBE_BIAS_GUIDE.md`, `MSC_METHOD_GUIDE.md`.

---

## 4. What we need back from you

1. **Pick at most 20 items total** (figures + `.md` reports combined) that
   you will actually use or reference when drafting Sections 2.1, 2.2, 2.3.
2. For each item, say **which subsection it goes in** and **whether it is a
   main‑text figure or supplementary**.
3. Flag any gap: if your narrative needs a figure we did not produce, say so
   explicitly and describe the intended panel layout.

Once we have your list, we will provide a ~5‑line quantitative description
per selected figure (what it shows, numbers to quote, sample size, caveats)
so you can write the captions and main text directly.
