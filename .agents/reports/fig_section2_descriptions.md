> **⚠ SUPERSEDED (2026-04-24)** — Figure descriptions in this document
> correspond to the n=5 era of Section 2 (MSC → ImCoh transition narrative)
> and predate the 2026-04-15 ImCoh reset. The current canonical handoff
> is [`MULTISCALE_TASK_TRACE_FOR_WRITING.md`](MULTISCALE_TASK_TRACE_FOR_WRITING.md).
> Methodological content (probe-bias figures, ImCoh vs MSC comparisons)
> is still valid conceptually; numerical values must be regenerated
> against `fc_method="imcoh_abs"` on the n=10 cohort before quoting.

# Section 2 figure descriptions (for writing agent)
## [SUPERSEDED — see MULTISCALE_TASK_TRACE_FOR_WRITING.md]

Quantitative handoff for drafting Sections 2.1–2.3 of `notes_imcoh.tex`
(MSC -> ImCoh rewrite). All numbers computed from cache
(`data/cache/{msc,imcoh,lrg,imcoh_lrg}`) on beta / alpha at `rsPre` unless
stated; same-/cross-probe means use |FC| with zero diagonal and
`re.match(r"([A-Za-z]+'?)", label)` probe extraction. Source reports:
`WRITING_AGENT_BRIEFING.md`, `IMCOH_PAT02_AND_CONTROLS.md`,
`.agents/guides/02_methods/PROBE_BIAS_GUIDE.md`,
`SECTION2_FIGURES_HANDOFF.md`.

Patients: 5 total (`Pat_02, 03, 05, 07, 08`); `Pat_06` excluded (missing
task phases). `Pat_03` recorded at 1024 Hz vs 2048 Hz others — flagged
as outlier / negative control wherever it is part of the sample.

---

## R-group: report summaries

### R1 — ImCoh H2a per-band (from `WRITING_AGENT_BRIEFING.md` §4.2 + §4.3)

"Task trace" = VI(rsPre,rsPost) > VI(taskTest,rsPost). Full k-range
(k = 2 … N-1, ~116 levels). "Unanimous cells" = (band, k) where all 5
patients share the sign. Scalar contrast = per-patient Δ_H2a averaged
over k (briefing §4.3 provides numbers for beta only).

| Band        | H2a unanimous cells | Mean scalar contrast (across patients) |
|-------------|:-------------------:|:---------------------------------------:|
| beta        | **82** (k = 16–108; all 5 patients positive) | **+0.217 ± 0.036** (group mean, briefing §4.3) |
| low_gamma   | 59 (k = 21–108; all 5 patients positive) | not given in briefing |
| alpha       | 32 (k = 17–61; Pat_02 = −0.05 breaks unanimity) | not given in briefing |
| delta       | 16 (k = 8–71) | not given in briefing |
| theta       | 5 (k = 45–49) | ≈ +0.061 (from `IMCOH_PAT02_AND_CONTROLS.md` Step 4c) |
| high_gamma  | 0  | not given in briefing |

Beta / low_gamma are the only "all-5-patients-positive" bands across a
wide k-range. Reference VI distances in beta (averaged across k, briefing
§4.3): VI(Pre,Post) ≈ 1.1, VI(TT,Post) ≈ 0.9, VI(TL,TT) ≈ 0.7 → Δ_H2a ≈
0.2 ≈ 20 % of VI(Pre,TT).

### R2 — Pat_02 H2a rank by band (from `IMCOH_PAT02_AND_CONTROLS.md` Step 3)

Pat_02 is NOT globally weak — only H2a-specific and band-dependent.

| Hypothesis          | Pat_02 typical rank (1=weakest, 5=strongest) | Interpretation                  |
|---------------------|:--------------------------------------------:|---------------------------------|
| H1                  | 2–4 / 5                                      | Normal to strong                |
| H2a (delta/theta/alpha) | **1–2 / 5**                              | Weak trace in slow bands        |
| **H2a beta**        | **5 / 5 (strongest)** — scalar +0.262        | Strong trace in beta            |
| H2b                 | 4–5 / 5                                      | Strong approach                 |
| H3                  | 4–5 / 5                                      | Strongest in most bands         |

Since the paper headline is beta H2a, Pat_02 is the *strongest* patient
for that claim; its slow-band weakness is a *separate* band-specific
pattern.

### R3 — Probe-bias enrichment by scale (from `PROBE_BIAS_GUIDE.md`)

Same-community / same-probe enrichment ratio under MSC (beta) — the
inflation that motivates the estimator switch.

| n_communities | Pat_02 | Pat_03 | Pat_05 | Pat_07 | Pat_08 |
|:-------------:|:------:|:------:|:------:|:------:|:------:|
| 3             | 0.9×   | 1.0×   | 1.0×   | 1.1×   | 1.0×   |
| 5             | 1.6×   | 1.0×   | 1.1×   | 1.8×   | 1.3×   |
| 10            | 3.9×   | 1.4×   | 1.4×   | 2.0×   | 1.5×   |
| 15            | 5.6×   | 1.7×   | 3.7×   | 2.7×   | 2.1×   |
| 20            | 6.2×   | 4.0×   | 3.7×   | 2.8×   | 2.4×   |
| 30            | 8.2×   | 5.8×   | 6.5×   | 4.1×   | 4.0×   |

Under ImCoh (briefing §6, same metric) enrichment at n = 10 drops to
1.1–2.4× and at n = 30 to 1.3–3.8×. The bias is concentrated at coarse
scales and effectively eliminated for n ≤ 10 under ImCoh.

---

## F-group: figure descriptions

### F1 — MSC adjacency heatmap grid, per patient (phases × bands)
Path: `/home/giulio/Documents/research/neural_networks/lrgeegfc/data/outputs/figures/section2/fig_A/fig_A_msc_Pat_02.pdf` (Pat_05, Pat_08 companions alongside).
1. **What it shows:** Heatmap of |MSC| (zero diagonal, probe-sorted) in a 4-phase × 6-band grid; probe blocks outlined. The conspicuous diagonal blocks = same-probe inflation — the visual smoking gun for Section 2.1.
2. **Key numbers (Pat_02, rsPre, beta):** |MSC| mean = 0.066, median = 0.025, q95 = 0.269; density(|MSC|>0.05) = 34 %. Same-probe / cross-probe mean ratio = **5.33×** (sp = 0.254, cp = 0.048).
3. **Sample:** 3 patients (Pat_02 / 05 / 08); 6 bands × 4 phases each.
4. **Caveats:** Pat_02 / 05 / 08 are the "with-coords" triad; Pat_03 (1024 Hz, MSC 3× higher than any other patient) and Pat_07 are not in this grid but are in R3. Only a visual/qualitative panel — quote numbers from F2.
5. **Paired with:** F2 (same grid rendered as a force-directed network) — "adjacency heatmap vs network" subfloats a/b.

### F2 — MSC network-view grid, per patient
Path: `/home/giulio/Documents/research/neural_networks/lrgeegfc/data/outputs/figures/section2/fig_A/fig_A_network_msc_Pat_02.pdf`
1. **What it shows:** Same (phases × bands) grid as F1, drawn as Kamada-Kawai on LRG ultrametric. Same-probe edges in probe-colour; cross-probe black, rank-gamma fade (γ_MSC = 2). Collapse to probe-coloured "hairballs" visualises the coarse-scale probe dominance (R3).
2. **Key numbers:** See R3 — at n=10 Pat_02 MSC enrichment = 3.9×; at n=30 = 8.2×.
3. **Sample:** 3 patients × 4 phases × 6 bands.
4. **Caveats:** Force-directed layout is cosmetic; quantitative claims must rely on R3.
5. **Paired with:** F1 (adjacency a / network b).

### F3 — MSC weight distribution, same-probe vs cross-probe
Path: `/home/giulio/Documents/research/neural_networks/lrgeegfc/data/outputs/figures/section2/fig_C/fig_C1a_weight_dist_msc_Pat_02_rsPre.pdf`
1. **What it shows:** Histogram / KDE of |MSC| split by same-probe (orange) and cross-probe (blue) edges, per patient × phase.
2. **Key numbers (Pat_02 rsPre beta):** same-probe mean 0.254 vs cross-probe 0.048 (ratio **5.33×**); alpha 0.251 vs 0.046 (**5.49×**). Pat_05 beta 4.25×, Pat_08 beta 2.62×, Pat_07 beta 2.24×, Pat_03 beta 3.38×.
3. **Sample:** 3 patients × 3 phases (9 files).
4. **Caveats:** Pat_03 is in R3 / WRITING_AGENT_BRIEFING §6 but is NOT in this fig_C1 family (only Pat_02/05/08 rendered). Same-probe pairs ≈ 9–10 % of edges (PROBE_BIAS_GUIDE).
5. **Paired with:** F4 (ImCoh equivalent) as subfloats a/b in a single "MSC vs ImCoh weight-distribution" figure.

### F4 — ImCoh weight distribution, same-probe vs cross-probe
Path: `/home/giulio/Documents/research/neural_networks/lrgeegfc/data/outputs/figures/section2/fig_C/fig_C1b_weight_dist_imcoh_Pat_02_rsPre.pdf`
1. **What it shows:** Same split as F3 but on |ImCoh|. Distributions nearly overlap — the raw-weight signature that volume conduction is gone.
2. **Key numbers (Pat_02 rsPre beta):** same-probe mean 0.0071 vs cross-probe 0.0039 (ratio **1.81×**); alpha 0.0090 vs 0.0067 (**1.35×**). Pat_03 beta 0.99×; Pat_05 beta 1.59×, Pat_07 beta 1.68×, Pat_08 beta 2.43×. MSC→ImCoh drops the ratio from 3–5× (typical) to 0.6–2× (briefing §3).
3. **Sample:** 3 patients × 3 phases (9 files).
4. **Caveats:** ImCoh mean is ~0.004 (near noise floor ~0.001, briefing §7 ·5); the residual 2× (e.g. Pat_08 beta) is genuine lagged coupling, not artifact.
5. **Paired with:** F3.

### F5 — Weight KDE overlay, all 5 patients (beta)
Path: `/home/giulio/Documents/research/neural_networks/lrgeegfc/data/outputs/figures/section2/fig_C/fig_C3_weight_overlay_all_patients_beta.pdf`
1. **What it shows:** KDE of same-probe vs cross-probe |FC| overlaid across all 5 patients, MSC on top / ImCoh on bottom (beta). Distributions cluster in MSC by patient; in ImCoh they collapse to a tight near-identical density.
2. **Key numbers:** Group MSC same-probe mean range 0.24–0.39, cross-probe 0.05–0.21; ImCoh same-probe 0.003–0.009, cross-probe 0.002–0.005. Same-probe / cross-probe ratio range shrinks from MSC [2.24–5.33]× to ImCoh [0.99–2.43]× (beta, rsPre, all 5 patients; see F3/F4).
3. **Sample:** All 5 patients, rsPre, beta (also `_alpha.pdf` available).
4. **Caveats:** Pat_03 (1024 Hz) is one of the overlaid curves — flag distinctly, its MSC is 3× higher than any other patient's.
5. **Paired with:** standalone summary panel for 2.1 → 2.2 bridge; or subfloat with alpha companion `fig_C3_weight_overlay_all_patients_alpha.pdf`.

### F6 — MSC adjacency, all 5 patients, fixed (beta, rsPre)
Path: `/home/giulio/Documents/research/neural_networks/lrgeegfc/data/outputs/figures/section2/fig_B/fig_B_beta_rsPre_MSC_vs_ImCoh.pdf`
1. **What it shows:** MSC row (top) of the 5-patient matrix grid at fixed (beta, rsPre). Checkerboard from same-probe blocks is visible in every patient.
2. **Key numbers:** Same/cross ratio at (beta, rsPre) across patients: Pat_02 5.33×, Pat_03 3.38×, Pat_05 4.25×, Pat_07 2.24×, Pat_08 2.62× (median ≈ 3.4×).
3. **Sample:** 5 patients × 1 band × 1 phase (MSC row); file also carries the ImCoh bottom row (F7).
4. **Caveats:** Vmax set per-band in the script; visual comparisons should not infer magnitude from the colour scale alone.
5. **Paired with:** F7 (same file, bottom row) — natural two-row single figure.

### F7 — ImCoh adjacency, all 5 patients, fixed (beta, rsPre)
Path: `/home/giulio/Documents/research/neural_networks/lrgeegfc/data/outputs/figures/section2/fig_B/fig_B_beta_rsPre_MSC_vs_ImCoh.pdf` (bottom row)
1. **What it shows:** ImCoh row at the same (beta, rsPre). Diagonal probe blocks are strongly reduced; residual structure is sparse, near-uniform, closer to cross-probe.
2. **Key numbers:** Same/cross ratio collapses to Pat_02 1.81×, Pat_03 0.99×, Pat_05 1.59×, Pat_07 1.68×, Pat_08 2.43× (median ≈ 1.68×). Mean |ImCoh| ≈ 0.003 (near the ~0.001 noise floor, briefing §7·5).
3. **Sample:** 5 patients × 1 band × 1 phase.
4. **Caveats:** Noise-floor proximity means dense ImCoh (unsparsified) is required — circular-shift surrogates are invalid for ImCoh (briefing §3).
5. **Paired with:** F6 (one figure).

### F8 — Spectral distribution MSC vs ImCoh (median + IQR + 5–95 %)
Path: `/home/giulio/Documents/research/neural_networks/lrgeegfc/data/outputs/figures/section2/fig_C/fig_C2_spectral_distribution_msc_rsPre.pdf` + `fig_C2_spectral_distribution_imcoh_rsPre.pdf`
1. **What it shows:** Per-frequency median / IQR / 5–95 % envelope of Welch coherence values over all edges, per band panel. MSC shows a systematic baseline offset (volume conduction raises the floor across all f); ImCoh is centred near 0 with symmetric tails.
2. **Key numbers:** Not directly extractable from cache without re-running Welch; the briefing §3 notes "computation time identical (~35 s per patient/phase)" and the MSC same/cross contrast of 3–5× (vs ImCoh ~1×) is the take-home.
3. **Sample:** 3 patients aggregated × 3 phases, 6 bands per panel.
4. **Caveats:** `value not computable from cache` for per-frequency percentiles — numbers need a rerun of `fig_C2_spectral_distribution.py`. Use the figure itself qualitatively.
5. **Paired with:** Standalone mid-2.2 figure, or subfloat a/b (msc/imcoh) in one panel.

### F9 — Enrichment-heatmap, MSC vs ImCoh (redesigned master summary)
Path: `/home/giulio/Documents/research/neural_networks/lrgeegfc/data/outputs/figures/section2/fig_D/fig_D1_enrichment_heatmap_redesigned.pdf`
1. **What it shows:** Heatmap of same-community / same-probe enrichment across (patient × n_communities), MSC vs ImCoh panel-by-panel. This is the single-page quantitative summary of Section 2.
2. **Key numbers:** See R3 (MSC) — Pat_02 rises 0.9→8.2×, Pat_08 1.0→4.0× across n = 3→30. Under ImCoh (briefing §6), Pat_02 rises only 1.0→2.5×, Pat_08 1.0→3.8×. Coarse-scale (n≤10) bias reduction is 3–8× → 1–2×.
3. **Sample:** All 5 patients × 6 scales × 2 methods = 60 cells.
4. **Caveats:** Enrichment is computed on LRG community labels, so it inherits the LRG random-seed choices; Pat_03 is included but flag separately.
5. **Paired with:** Standalone main-text figure (Section 2.2 closer); optionally pair with F10.

### F10 — Bias-reduction 3-panel (log2 ratio + curves + strips)
Path: `/home/giulio/Documents/research/neural_networks/lrgeegfc/data/outputs/figures/section2/fig_D/fig_D2_bias_reduction_beta_rsPre.pdf`
1. **What it shows:** (i) log2(same/cross ratio) heatmap per patient × band; (ii) mean ± SD curves of that ratio vs scale; (iii) per-scale strips of individual patients. Shows the directionality + magnitude of bias change per patient.
2. **Key numbers (beta, rsPre):** log2(MSC ratio) across patients = {2.41, 1.76, 2.09, 1.16, 1.39}; log2(ImCoh ratio) = {0.86, −0.01, 0.67, 0.75, 1.28}. Mean reduction ≈ −1.1 log2 units (~2.1× reduction).
3. **Sample:** 5 patients, 1 band, 1 phase per file (4 files cover {alpha, beta} × {rsPre, taskLearn}).
4. **Caveats:** Pat_03 log2 = −0.01 (ImCoh) — essentially bias-free after the switch, strongest single-patient confirmation; also most extreme MSC outlier (briefing §7·7). Pat_06 absent (no task phases).
5. **Paired with:** F9 (2-subfigure summary of probe bias and its removal).

### F11 — Brain connectome, nilearn glass brain
Path: `/home/giulio/Documents/research/neural_networks/lrgeegfc/data/outputs/figures/section2/fig_E/fig_E1_brain_connectome_beta_rsPre.pdf`
1. **What it shows:** Top-N edges rendered on nilearn glass brain at (beta, rsPre), MSC vs ImCoh side-by-side. Under MSC edges align with probe trajectories; under ImCoh edges jump across probes / hemispheres.
2. **Key numbers:** Only patients with MNI coordinates are rendered — `PATIENTS_WITH_COORDS = {Pat_02, Pat_03, Pat_05}` (from `_shared.py`). Pat_07, Pat_08 omitted for lack of coords.
3. **Sample:** 3 patients × (alpha,beta) × (rsPre, taskLearn, rsPost).
4. **Caveats:** **Missing MNI coords for Pat_07 and Pat_08.** Thresholding is top-N edges, so absolute magnitudes differ between MSC and ImCoh panels. Pat_03 (1024 Hz, MSC outlier) present.
5. **Paired with:** F12 (force-directed network triad) to complement anatomical vs topological view.

### F12 — Force-directed 3-patient × 2-method network grid
Path: `/home/giulio/Documents/research/neural_networks/lrgeegfc/data/outputs/figures/section2/fig_E/fig_E2_spring_Pat_02_Pat_05_Pat_08_beta_rsPre.pdf`
1. **What it shows:** 3 patients × 2 methods grid of Kamada-Kawai layouts on the LRG ultrametric. MSC column = probe-coloured clumps; ImCoh column = modular structure crossing probes.
2. **Key numbers:** Graph sizes: Pat_02 N = 117, Pat_05 = 118, Pat_08 = 120 (PROBE_BIAS_GUIDE §"Probe Structure"). Under ImCoh at n = 10, same-probe enrichment drops to 1.3× (Pat_02), 1.2× (Pat_05), 2.4× (Pat_08) vs 3.9× / 1.4× / 1.5× under MSC.
3. **Sample:** 3 patients × 2 methods × 1 (band, phase) per file — 6 files cover {alpha, beta} × {rsPre, taskLearn, rsPost}.
4. **Caveats:** Layout depends on LRG seed and γ = 30 / 2 recipe (handoff §1 "Finalized visual recipe") — not a quantitative visualisation.
5. **Paired with:** F11.

### F13 — Nodewise graph-metrics-on-graph (strength / clustering / participation / degree)
Path: `/home/giulio/Documents/research/neural_networks/lrgeegfc/data/outputs/figures/section2/fig_F/fig_F_msc_Pat_02_beta_rsPre.pdf` (+ `fig_F_imcoh_Pat_02_beta_rsPre.pdf`)
1. **What it shows:** 4-panel per (method, patient, band, phase): node colour = one of {strength, weighted clustering, participation coefficient wrt probes, degree}. Participation is the key: near 0 under MSC (edges trapped within probe) → closer to 1-1/n_probe under ImCoh (edges spread across probes).
2. **Key numbers (Pat_02 rsPre beta, computed):** network mean |FC| MSC 0.066 vs ImCoh 0.0042 → strength scales linearly with mean. Density(|FC|>0.05) MSC 0.34 vs ImCoh 0.01 — degree definitions must be thresholded carefully.
3. **Sample:** 3 patients × 2 methods × 2 bands × 3 phases = 36 files.
4. **Caveats:** Weighted vs thresholded metrics behave very differently for ImCoh near the noise floor; the same colour-scale is NOT shared between MSC and ImCoh panels.
5. **Paired with:** F14 (ImCoh counterpart of same patient/band/phase) as subfloats a/b of a single "metrics on the graph" figure.

### F14 — Nodewise graph-metrics, ImCoh counterpart
Path: `/home/giulio/Documents/research/neural_networks/lrgeegfc/data/outputs/figures/section2/fig_F/fig_F_imcoh_Pat_02_beta_rsPre.pdf`
1. **What it shows:** Same four panels as F13 under ImCoh. Participation map looks much more uniform across the brain than under MSC.
2. **Key numbers:** See F13 caveats; participation coefficient is the panel that swings most (definition in `_shared.py:compute_participation_coefficient` is 1 − Σ_s (k_{i,s}/k_i)²).
3. **Sample:** Same as F13.
4. **Caveats:** Same as F13 — node colour scales not shared with MSC panel.
5. **Paired with:** F13.

### F15 — Legacy single-page probe-bias summary (beta only)
Path: `/home/giulio/Documents/research/neural_networks/lrgeegfc/data/outputs/figures/section2/fig_01/fig_01A_probe_enrichment_heatmap.pdf`, `fig_01B_weight_distributions_beta.pdf`, `fig_01C_enrichment_vs_scale_beta.pdf`
1. **What it shows:** Three-panel legacy summary — (A) MSC vs ImCoh enrichment heatmap across 5 patients × 6 scales; (B) same-/cross-probe weight dists beta; (C) enrichment curves vs n_communities. Superseded by F9/F10 but cleaner for a single-page figure.
2. **Key numbers:** See R3 (MSC panel) and briefing §6 (ImCoh panel). At n = 10, MSC mean 2.1×, ImCoh mean 1.5× (across 5 patients, beta).
3. **Sample:** All 5 patients, beta-only.
4. **Caveats:** Pat_03 in panel A is an MSC outlier — its row is visibly lighter (lower enrichment) because its flat MSC hides structure; do NOT interpret as "less bias".
5. **Paired with:** Optional supplementary single-page summary; redundant with F9 + F10 for the main text.

### F16 — Laplacian eigenvalue spectrum, MSC vs ImCoh, 4 phases overlaid
Path: `/home/giulio/Documents/research/neural_networks/lrgeegfc/data/outputs/figures/section2/fig_G/fig_G1_eigenvalue_spectrum_Pat_02_MSC_vs_ImCoh.pdf`
1. **What it shows:** Sorted Laplacian eigenvalues for (Pat_02 / 05 / 08) in a 2×(rows for MSC/ImCoh) panel with 4 phases overlaid. Huge MSC λ_max dominates; ImCoh spectrum is 10× more compact with a cleaner Fiedler gap.
2. **Key numbers (rsPre, beta):** λ_max / λ_2:  Pat_02 MSC 18.1 / 0.53 (ratio 34.2) → ImCoh 2.44 / 0.109 (ratio 22.4); Pat_05 19.1 / 2.24 (8.6) → 1.07 / 0.156 (6.9); Pat_08 29.8 / 0.77 (38.6) → 1.94 / 0.134 (14.5). Condition number λ_max/λ_2 drops ~1.5–2.6× under ImCoh.
3. **Sample:** 3 patients × 1 file each (all 4 phases inside).
4. **Caveats:** Not a surrogate-tested statistic; qualitative structural check.
5. **Paired with:** fig_G2 (below) as subfigures a (spectrum) and b (susceptibility) of one Laplacian / LRG-thermodynamics figure.

### F17 — Susceptibility C(τ) MSC vs ImCoh (Pat_02, beta, rsPre) — **new**
Path: `/home/giulio/Documents/research/neural_networks/lrgeegfc/data/outputs/figures/section2/fig_G/fig_G2_susceptibility_Pat02_beta_MSC_vs_ImCoh.pdf`
1. **What it shows:** Two panels (MSC / ImCoh), each plotting C(τ) on log-x, with vertical dashed lines at 1/λ_max, 1/λ_2 and the empirical τ* = argmax_{τ>10⁻⁴} C(τ). Demonstrates that the LRG specific-heat curve sits on a completely different τ scale under ImCoh (pushed 10× later to the right).
2. **Key numbers (Pat_02 rsPre beta):** MSC: 1/λ_max = 0.0554, 1/λ_2 = 1.89, τ* ≈ 1.05e-4, C_max ≈ 17.0. ImCoh: 1/λ_max = 0.409, 1/λ_2 = 9.18, τ* ≈ 1.08e-4, C_max ≈ 0.117. MSC C is ~150× larger than ImCoh C — the Laplacian is ~10× stiffer.
3. **Sample:** 1 patient / 1 band / 1 phase (target headline case).
4. **Caveats:** C near τ = 0 is numerically noisy (plateau); the τ > 1e-4 mask is a pragmatic filter. `entropy_C` array in cache has length 400 (same as `entropy_tau`), plotted directly.
5. **Paired with:** F16 — one "LRG thermodynamics" figure with spectrum (a) and susceptibility (b).

---

## Quick suggested figure layout for 2.1 / 2.2 / 2.3

- **Section 2.1 — the MSC probe bias exists.** Main text: F1+F2 (single figure, a/b), F3+F4 (single figure, a/b), F5 (single figure, 2 bands). Supplementary: F15.
- **Section 2.2 — ImCoh removes it.** Main text: F6+F7 (single figure, stacked rows), F9 (single figure), F10 (single figure). Supplementary: F8.
- **Section 2.3 — graph-level consequences.** Main text: F11+F12 (single figure, a/b), F13+F14 (single figure, a/b), F16+F17 (single figure, a/b).
