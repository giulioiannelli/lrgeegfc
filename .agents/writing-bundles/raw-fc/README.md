---
name: README
type: writing-bundle
era: IMCOH_ABS × COHORT_N10
status: historical
created: 2026-04-30
updated: 2026-07-04
---

# Raw-FC phase-distance audit — writing bundle

> **⚠ Historical (banner added 2026-07-04, `.agents/` reorg).** This bundle is
> **substrate context, not a standalone result.** Three things below are
> superseded:
> 1. **Framing.** "Raw FC as headline Result 1" is retired — since 2026-06-18
>    (`feedback_results_only_in_laplacian_framework`) raw FC is the **comparison
>    baseline, never a result**. The trace result lives in the LRG/Laplacian
>    framework (cophenetic `ρ^coph` + Grassmann); see `.agents/preprint/headlines/`.
> 2. **`imcoh_sq` substrate.** The "two FC substrates" section is unsupported —
>    no `*_imcoh_sq.pdf` or `data/audit/raw_fc_phase_distance_imcoh_sq/` exists.
>    Only `imcoh_abs` remains.
> 3. **Figure names + Pat_03.** The canonical renders are the **`_imcoh_abs`**
>    PDFs in this folder; the 7 old unsuffixed duplicates were archived to
>    `archive/superseded-renders/`. Body references to `figN_*.pdf` (unsuffixed)
>    and to "Pat_03 = 1024 Hz outlier" are stale (Pat_03 is a full cohort member).
>
> The `imcoh_abs` figures + numbers are kept as accurate substrate context.

**Self-contained handoff for the writing agent.** One README +
parallel sets of figures for **two FC substrates** as a
sensitivity-check pair:

- `*_imcoh_abs.pdf` (default) — the canonical
  `imcoh_abs = mean_f |Im(S_ij)/√(S_ii·S_jj)|` substrate. Source
  CSVs at `data/audit/raw_fc_phase_distance/`.
- `*_imcoh_sq.pdf` — the squared variant
  `imcoh_sq = mean_f (Im(S_ij)/√(S_ii·S_jj))²` (band-averaged power-
  weighted, see §11 caveat). Source CSVs at
  `data/audit/raw_fc_phase_distance_imcoh_sq/`. **figS1 and figS3
  are not regenerated for `imcoh_sq`** because they rely on the
  expensive within-RPre split-half null (audit_25), which was not
  re-run for this substrate.

**Verdict robustness across substrates** (`final_verdict_table.csv`
on each substrate):

| band | imcoh_abs verdict | imcoh_abs n_S | imcoh_sq verdict | imcoh_sq n_S |
|---|---|:---:|---|:---:|
| δ      | trace        | 6/10 | borderline   | 5/10 |
| θ      | drift-only   | 3/10 | drift-only   | 3/10 |
| α      | trace        | 8/10 | trace        | 8/10 |
| β      | trace        | 7/10 | **trace**    | **8/10** |
| low_γ  | trace        | 7/10 | trace        | 7/10 |
| high_γ | borderline   | 5/10 | borderline   | 5/10 |

5 of 6 bands give the same verdict; β actually strengthens to 8/10
under squaring; only δ flips from `trace` to `borderline` (drops one
patient: 6/10 → 5/10). The headline finding is robust to the
abs-vs-squared choice.

> **Terminology.** Everything below uses the **trace / anchor / reset
> / emergent** taxonomy from
> `.agents/guides/01_project/terminology.md`. Our `T_d^(d_S) < 0`
> finding is a **trace** (task reorganized AND the change persists
> into `RPost`) — never the bare word "persistence", which silently
> flips to "anchor" (module unchanged) on a careful reader.

---

## 1. Headline

### Head (one paragraph)

At the bare-substrate level — `imcoh_abs` FC matrices on the
per-patient, per-band giant-component intersection `V*`, no diffusion,
no hierarchy — the cohort (n = 10) shows a **task-shaped trace** on
the rank distance `d_S`. The effect is small, does not survive
multiple-comparison correction, and is fragile to one outlier
(Pat_06). The geometric positive control (within-task < within-rest)
holds robustly and is the strongest single piece of cohort evidence.

### Claim (one sentence)

In 4 of 6 frequency bands (δ, α, β, low_γ), `task_test` leaves a
band-specific, structurally-coherent trace in `rest_post`: the rank
ordering of FC edges in `rest_post` is closer to `task_test` than to
`rest_pre`, in 6–8 of 10 patients per band.

### Key numbers (drop-in)

- **Cohort:** n = 10 (`PATIENTS_4PHASE`). Pat_03 included (1024 Hz;
  Z and ratios are dimensionless w.r.t. fs), marked orange triangle.
- **Bands with trace:** δ, α, β, low_γ (cohort-median `T_d^(d_S) < 0`,
  trace count ≥ 6/10).
- **Drift-only band:** θ (trace count 3/10 on `d_S`).
- **Borderline band:** high_γ (cohort-median `T_d^(d_S) = +0.013`).
- **Strongest single rank-distance cell:** α × `d_S` = 8/10.
- **Strongest cohort-median asymmetry:** δ × `d_S` = −0.131.
- **Best uncorrected p:** 0.042 (β × `d_F`).
- **Smallest BH-FDR q:** 0.29 — **zero cells survive q = 0.05.**
- **Geometric positive control:** cohort-median `d_S(TL, TT)` ≈
  0.5 × `d_S(RPre, RPost)` in every band.
- **Structural-vs-drift:** all 6 bands flag "structural-shift";
  per-patient Spearman ρ between `T_d^(d_S)` and `T_d^(d_F)` is
  0.78–0.92.

---

## 2. Operational definitions

### Cohort

`PATIENTS_4PHASE = {Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15}`,
n = 10. Pat_03 = 1024 Hz outlier (others 2048 Hz). Pat_06 =
load-bearing patient (largest `T_d^(d_F)`; motivated the audit).

### FC substrate

`imcoh_abs(i, j) = mean_f |Im(S_ij(f)) / √(S_ii(f)·S_jj(f))|`
band-averaged from Welch (`nperseg = 4096` at 2048 Hz, `2048` at 1024
Hz). Band-averaged absolute imaginary coherence (Nolte 2004 / Ewald
2012). Range [0, 1].

### V* — per-patient, per-band substrate

`V*(p, b, comparison) = ⋂_{φ ∈ comparison} GC(|A^φ(p, b)|)` —
intersection across the involved phases of giant-component node sets.
Adjacency matrices are restricted to `V*` before any distance
computation.

### Three distances on `triu(A, k=1)` — DEFENSIBLE FRAMING

| distance | formula | what it measures | drift exposure | role |
|---|---|---|---|---|
| **`d_S`** | `1 − ρ_S(triu A_a, triu A_b)` | **rank-only** — ordering of edges | LOW (rank-invariant) | **load-bearing** |
| `d_P`     | `1 − ρ_P(triu A_a, triu A_b)` | magnitude-weighted complement to `d_S` (mixes rank and amplitude) | partial | convergence check |
| `d_F`     | `‖A_a − A_b‖_F / √(‖A_a‖_F·‖A_b‖_F)` | **amplitude-only** — Euclidean on raw values | HIGH | drift diagnostic |

**Critical framing notes** (do NOT overreach):

- `d_P` is NOT "volume + topology". Pearson on `triu(A)` is a
  magnitude-weighted linear correlation; it does NOT decompose into
  "topology" plus "volume".
- `d_S` and `d_P` are NOT "orthogonal". Cohort Spearman ρ between
  `T_d^(d_S)` and `T_d^(d_P)` is 0.85–0.95 per band — they are
  strongly correlated. The β cell is the load-bearing **convergence**
  claim (`d_S` 7/10, `d_P` 7/10), not "orthogonal agreement".
- `d_F` triangle pattern (5/6 bands cohort-wide) is the **predicted
  footprint of monotone session drift**, not a structural finding.
  Report as diagnostic only.

### Trace scalar

`T_d = d(task_test, rest_post) − d(rest_pre, task_test)`. Negative ⇒
`rest_post` closer to `task_test` than `rest_pre` is — the **trace**
pattern. Cohort trace count `n_trace = #{p : T_d(p) < 0}` out of 10.

(Existing CSV column `n_persist` is the same quantity; rename pending
on next artefact rebuild.)

### Null models + positive control

- **Within-`rest_pre` split-half null** (audit_25): 50 non-overlapping
  Welch segments → distribution `N`; `Z = (d_obs − median(N))/MAD(N)`.
  Significance axis only — does NOT control for monotone drift.
- **4-phase geometry positive control** (audit_26): all 6 inter-phase
  distances on cached FCs; cohort medians + Q1/Q3.
- **Drift-triangle null** (audit_28, **Pat_06 only**): K=4 chunk
  decomposition + linear fit `D ≈ α + β·gap`. R² ≤ 0.05 → no monotone
  within-baseline drift. Cohort sweep deferred.

---

## 3. Per-band verdict (cohort median, n = 10)

From `final_verdict_table.csv` (column `n_persist` = trace count):

| band | d_S(RPre, TT) | d_S(TT, RPost) | d_S(RPre, RPost) | T_d^(d_S) | n_trace d_S | n_trace d_P | n_trace d_F | verdict |
|------|---:|---:|---:|---:|:---:|:---:|:---:|---|
| δ          | 0.464 | 0.333 | 0.538 | −0.131 | 6/10 | 7/10 | 7/10 | **trace** |
| θ          | 0.427 | 0.375 | 0.384 | −0.052 | 3/10 | 5/10 | 6/10 | drift-only |
| α          | 0.430 | 0.395 | 0.576 | −0.035 | 8/10 | 7/10 | 7/10 | **trace** |
| β          | 0.430 | 0.389 | 0.419 | −0.040 | 7/10 | 7/10 | 8/10 | **trace** |
| low_γ      | 0.584 | 0.479 | 0.612 | −0.105 | 7/10 | 5/10 | 5/10 | **trace** |
| high_γ     | 0.420 | 0.433 | 0.445 | +0.013 | 5/10 | 7/10 | 6/10 | borderline |

**Per-band one-liner.**
- **δ.** Strongest cohort-median asymmetry on `d_S` (−0.131), 6/10.
- **θ.** Drift-only. Trace counts on rank distances below majority.
- **α.** Strongest single rank cell: `d_S` 8/10.
- **β.** Convergence cell — both `d_S` and `d_P` at 7/10. Best
  uncorrected significance lives here too.
- **low_γ.** Strong cohort-median asymmetry, 7/10 on `d_S`.
- **high_γ.** Borderline. Cohort-median `T_d^(d_S)` marginally
  positive but `d_P` 7/10.

**Headline sentence for the paper:**
> "4 of 6 bands (δ, α, β, low_γ) show a cohort-directional task-trace
> on `d_S`, with 6–8 of 10 patients per band."

---

## 4. Figures — what to look at, what to claim

### Fig 1 — Trace scatter on `d_S`  →  `fig1_trace_scatter_dS.pdf`

**The headline figure.** Per band (6 panels): x = `d_S(RPre, TT)`,
y = `d_S(TT, RPost)`, **one point per patient** (white circles, edge
blue; Pat_03 = orange triangle). Visual elements:

- **Red-shaded triangle** below the dashed identity = **trace zone**
  (points here mean `RPost` is closer to `TT` than `RPre` is).
- **Gray ellipse** = cohort 1σ covariance ellipse around the cohort
  mean.
- **Green star** = cohort mean (centroid of all 10 patient points).
- **Inset** = `trace: N/10` count.

**Drop-in caption:**
> "Fig. 1. Per-band trace scatter on the rank distance `d_S`. Each
> point is one patient (n = 10; Pat_03 marked as orange triangle —
> 1024 Hz methodological outlier whose Z-score is dimensionless w.r.t.
> sampling rate). Points in the red-shaded trace zone (below the
> dashed identity) indicate `rest_post` is closer to `task_test` than
> `rest_pre` is. The gray ellipse is the cohort 1σ covariance ellipse
> around the cohort mean (green star). 6–8 of 10 patients fall in
> the trace zone in 4 of 6 bands (δ, α, β, low_γ)."

### Fig 2 — Phase geometry  →  `fig2_phase_geometry.pdf`

**The strongest single piece of cohort evidence.** Per band (6
panels): the four reference phase pairs, **sorted per band** from
smallest to largest cohort median (so x-position 1 is always the
tightest pair within that band, x-position 4 is the widest). The
phase-pair colour is fixed across panels:

| pair | colour | role |
|---|---|---|
| TL↔TT | blue | within-task cohesion |
| TT↔RPost | red | the trace pair |
| RPre↔TT | orange | pre→task baseline shift |
| RPre↔RPost | gray | within-rest drift floor |

Per pair: translucent IQR box (Q1–Q3), bold median line in the pair
colour, jittered per-patient strip dots. Thin gray lines connect each
patient's four points across pair-positions (per-band-sorted), so a
monotone gray trace says "this patient agrees with the cohort sort".

**Cohort-wide reading:** in every band, the leftmost pair is **TL↔TT**
(within-task is tightest), and **RPre↔RPost** sits at or near the
right (within-rest is widest). The trace pair (TT↔RPost) typically
sits second (just above within-task), placing `RPost` closer to `TT`
than `RPre` is. The sequence "task-tight → trace pulls post toward
task → pre-to-task baseline shift → rest-rest drift floor" is the
geometric statement of the trace.

**Drop-in caption:**
> "Fig. 2. Four-phase geometry on `d_S` per band, with phase pairs
> sorted **per panel** from smallest to largest cohort median. Each
> per-pair element shows an IQR box, the cohort median (bold
> horizontal bar in the pair colour), and per-patient jittered dots
> (n = 10); thin gray lines connect each patient's four points across
> the band-specific sort positions. Pair colours are fixed across
> bands: TL↔TT blue (within-task cohesion), TT↔RPost red (trace pair),
> RPre↔TT orange (pre→task baseline shift), RPre↔RPost gray (within-rest
> drift floor). In every band the leftmost (smallest) pair is TL↔TT;
> the trace pair TT↔RPost is typically next, indicating `rest_post`
> sits closer to `task_test` than `rest_pre` does. The within-task
> bar being tighter than the within-rest bar (cohort-median
> `d_S(TL, TT) ≈ 0.5 × d_S(RPre, RPost)` in every band) rejects the
> alternative that the trace triangle could be produced by monotone
> session drift alone, which would inflate both equally."

### Fig 3 — `d_S` × `d_P` convergence  →  `fig3_dS_dP_convergence.pdf`

**The convergence story (β).** Per band: x = `T_d^(d_S)` (rank-only),
y = `T_d^(d_P)` (magnitude-weighted), one **numbered marker** per
patient (the patient ID is written inside the circle). Visual:

- **Green-shaded lower-left quadrant** = trace zone
  (`T_d^(d_S) < 0` AND `T_d^(d_P) < 0`).
- **Dashed crosshair** at zero, **dotted identity** line.
- **Blue star** = cohort mean (centroid).
- **Numbered circles**: green-edged when in trace zone, gray-edged
  otherwise; orange-edged for Pat_03.
- **Inset** = `trace: N/10` count for the trace zone.

The cohort `T_d^(d_S)` and `T_d^(d_P)` are strongly correlated
(Spearman ρ 0.85–0.95 per band) — not orthogonal axes; the
informative claim is **convergence** (β: 7/10 each).

**Drop-in caption:**
> "Fig. 3. Convergence of the rank distance `d_S` (Spearman) and the
> magnitude-weighted distance `d_P` (Pearson) on the trace scalar
> `T_d`. Each numbered circle is one patient (last two digits of the
> patient ID); green-edged in the trace zone (lower-left, both
> `T_d < 0`), gray-edged otherwise; orange-edged for Pat_03. The
> blue star marks the cohort mean. Cohort-wise Spearman ρ between
> the two axes is 0.85–0.95 per band — `d_S` and `d_P` are not
> orthogonal but co-confirming. Bands β, α and low_γ show the
> strongest convergence; θ is null on both."

### Fig 4 — Structural vs drift  →  `fig4_structural_vs_drift.pdf`

**The single strongest "this is not drift" argument.** Per band:
x = `T_d^(d_S)` (rank-only), y = `T_d^(d_F)` (amplitude-only). Each
patient is a **numbered circle** (last 2 digits of patient ID
inside). Visual aids:

- **Green-shaded lower-left quadrant** = trace zone (`T_d < 0` on
  both axes).
- **Dotted identity line** = perfect rank/amplitude agreement.
- **Dashed crosshair** at zero (separates the four sign quadrants).
- **Gray ellipse** = cohort 1σ covariance ellipse — its **elongation
  along the identity directly visualises the strong correlation**
  that defeats the drift hypothesis.
- **Blue star** = cohort mean.
- **Numbered circles**: green-edged in trace zone, gray-edged
  otherwise; orange-edged for Pat_03.
- **Inset** per band: `trace: N/10`, `ρ_S = X.XX` (per-band Spearman
  ρ between the two axes), `sign: N/10` (sign-agreement count).

The argument: under monotone session drift, `T_d^(d_F)` would be
uncorrelated with `T_d^(d_S)` (rank-invariance shields `d_S` from
uniform amplitude shifts) — the cohort cloud would be a circular
blob, not an ellipse. Instead the ellipse is elongated along the
identity in every band: per-patient Spearman ρ = 0.78–0.92, sign
agreement 7–9 of 10. **All 6 bands flag "structural-shift"** in
`Td_dS_vs_dF_band_contrast.csv`.

**Drop-in caption:**
> "Fig. 4. Structural-shift vs amplitude-drift contrast. Per band,
> one numbered circle per patient (last two digits of the patient ID,
> n = 10): x = `T_d^(d_S)` (rank-only trace scalar), y = `T_d^(d_F)`
> (amplitude-only). Green-edged circles fall in the trace zone
> (lower-left, both `T_d < 0`); gray-edged circles fall off the
> trace zone; Pat_03 is marked orange. The blue star marks the
> cohort mean; the gray ellipse is the cohort 1σ covariance ellipse
> and its elongation along the dotted identity line directly
> visualises the rank/amplitude correlation. Under monotone session
> drift, `T_d^(d_F)` would be uncorrelated with `T_d^(d_S)` because
> rank-invariance shields `d_S` from uniform amplitude shifts;
> instead the per-patient Spearman ρ between the two axes (inset per
> panel) is 0.78–0.92 across all six bands, with sign agreement 7–9
> of 10. The trace moves coherently on rank and amplitude axes —
> incompatible with pure amplitude drift."

### Fig 5 — Distance-class illustration  →  `fig5_distance_class_examples.pdf`

**The "why we report three distances" figure.** Two real (patient,
band, phase pair) cells, each as one row. Each row has three
super-columns separated by a controlled `wspace`:

- **Super-col A — heatmap pair (cols 1–2):** the two adjacency
  matrices side-by-side with a tight inner `wspace = 0.06`. Rows /
  columns of both heatmaps are **probe-sorted** so same-probe
  contacts form contiguous diagonal blocks. Tick labels are the sEEG
  probe IDs at probe-block centres (rotated 90° on x). White
  block outlines via `draw_probe_outlines`. Both heatmaps share a
  log-scale range (`LogNorm`) and a per-row colorbar attached to the
  rightmost heatmap via `imshow_colorbar_caxdivider` (5% size).
  **Colorbar ticks** sit only at integer powers of 10 strictly inside
  `[lo, hi]`; `cbar.minorticks_off()` removes the wide
  `2×10^k`/`6×10^k` minor labels.
- **Super-col B — log-log scatter (col 3):** `pcolormesh` density on
  a 60-bin geometric grid (cividis), white-dashed identity, red
  power-law fit (slope annotated), `d_S / d_P / d_F` inset. Tick
  positions are the **same `decade_ticks` set as the heatmap
  colorbar** — visually synchronised, no minor labels.
- **Super-col C — network pair (cols 4–5):** classical-MDS layout on
  the LRG heat-kernel transition distance (`mds_lrg_continuous`, the
  configured default for `imcoh_abs`) computed **once** on
  `0.5·(A_a + A_b)` so both panels share node positions; nodes
  coloured by sEEG probe (tab20), edges via `draw_network_edges` with
  the **canonical defaults** (rank-based scaling, per-fc_method
  `EDGE_RANK_GAMMA / WIDTH / ALPHA` from `config.const`) — same
  treatment as the `fig_E2_*` reference figures. Inner pair
  `wspace = 0.04`.

Per-row title sits in the row gap (`gs_outer.hspace = 0.55`), centred
across all three super-cols. Bottom-of-figure legend explains the
network and scatter conventions.

**Row 1 — Pat_14, low_γ, TL ↔ RPost** (`d_P = 0.08`, `d_S = 0.67`):
scatter monotone but compressed in log-log → topological
reorganization without amplitude redistribution.

**Row 2 — Pat_02, high_γ, TT ↔ RPost** (`d_P = 0.26`, `d_S = 0.13`):
scatter rank-preserving but non-linear → magnitude redistribution
without topological reorganization.

**Drop-in caption:**
> "Fig. 5. Two example phase-pair comparisons illustrating that `d_S`
> and `d_P` capture different aspects of FC reorganization. Each row
> shows two phases of `imcoh_abs` for one (patient, band) cell: row 1
> Pat_14 / low_γ / `task_learn ↔ rest_post`, row 2 Pat_02 / high_γ /
> `task_test ↔ rest_post`. Cols 1–2: adjacency heatmaps on per-pair
> `V*` with shared `vmax` (per-row colorbar). Col 3: log-log
> density-shaded scatter of upper-triangle entries; white dashed
> identity, red line is the log-log power-law fit (slope inset),
> distance triplet inset. Cols 4–5: spring-layout network drawings
> of the same `V*`-restricted adjacency, with shared node positions
> across the two phases for direct visual comparison; nodes coloured
> by sEEG probe (tab20), same-probe edges drawn in the probe colour
> at full alpha, cross-probe edges in dark gray with weight-scaled
> width and alpha. Pat_14 / low_γ shows compressed-but-monotone
> scatter (low `d_P = 0.08`, high `d_S = 0.67`) — topological
> reorganization without amplitude redistribution. Pat_02 / high_γ
> shows rank-preserving but non-linear scatter (low `d_S = 0.13`,
> high `d_P = 0.26`) — magnitude redistribution without topological
> reorganization."

### Fig 6 — 4-phase geometry (chord diagrams)  →  `fig6_4phase_geometry.pdf`

**The full inter-phase geometry, no redundant cells.** A 4×4
symmetric distance matrix would have 16 cells but only 6 unique
values (4 choose 2); the diagonal is trivially 0 and the lower
triangle is symmetric to the upper. We replace the matrix with a
**chord/arc diagram** per `(band × distance class)` cell:

- 4 phase markers on a baseline. Their **horizontal positions are
  computed via a 1-D spring layout** (weighted MDS): each phase pair
  has a spring with **stiffness `w_ij = 1/d_ij²`** (closer pairs pull
  much harder than far pairs — the squared denominator amplifies
  small-distance differences) and rest length `d_ij`. The equilibrium
  positions minimise `Σ w_ij · (|p_i − p_j| − d_ij)²`. Temporal order
  RPre < TL < TT < RPost is enforced via positive-delta
  parameterisation. Positions are rescaled to **[0, 1]** per panel.
- **Reference dashed grid** at x = 0, 0.25, 0.5, 0.75, 1 — the eye
  reads how far each phase has moved from its "evenly spaced"
  position.
- 6 arcs above the baseline connect each phase pair (no
  self-distances, no duplicated symmetric values).
- **Arc COLOUR** = within-panel distance value on a custom
  **2-colour green → gray → red** map (no yellow waypoint —
  semantically binary "close/far" with a clean neutral middle).
  Each panel has its **own colorbar** showing local `[vmin, vmax]`
  so the colour gradient is meaningful within the panel; distances
  are NOT directly comparable across panels through colour alone.
- **Arc height** is purely proportional to the horizontal span
  (semicircular). Distance is encoded by colour AND by the spring
  positions; height carries no extra information.
- The cohort-median value is annotated in small text at each arc
  apex.

Three rows × six columns = 18 panels; per-row colorbar on the right.

**Drop-in caption:**
> "Fig. 6. Four-phase distance geometry per band × distance class
> shown as chord diagrams. The four phase nodes (RPre, TL, TT,
> RPost) sit on a baseline in temporal order; the six unique pair
> distances are shown as arcs above. Arc colour encodes the
> cohort-median distance value (n = 10; viridis, per-row shared
> scale, see colorbar at right); arc height encodes the temporal
> gap of the pair (adjacent / skip-1 / skip-2). The within-task arc
> (TL↔TT) is consistently the lowest-distance arc (smallest gap and
> smallest distance), and the trace pair (TT↔RPost) is consistently
> shorter than (RPre↔TT) on `d_S`."

---

## 5. Supplementary figures

### Fig S1 — Z inflation diagnostic  →  `figS1_z_inflation.pdf`

Per `(band × distance)` (3 rows × 6 cols, 18 panels): for each of
the 10 patients (rows within a panel) we draw the within-`rest_pre`
null IQR as a gray ribbon (Q1–Q3 from 50 split-halves), a black
tick at the null median, and the three observed cross-phase
distances as coloured dots (red = RPre→TT, blue = RPre→RPost,
green = TT→RPost). Footer legend explains every element.

The gray ribbons are extremely narrow (1–4 % of their median), so
`Z = (d_obs − median)/MAD` blows up to 12–685 even when
`d_obs / null_median` is only ~1.5. **Use the ratio (or absolute
distance), never Z, as an effect size.**

### Fig S2 — `T_d` swarm per band, all three distances  →  `figS2_Td_swarm_all_distances.pdf`

Per band (6 panels, **independent y-axes** because per-band
variability differs by an order of magnitude — Pat_06's contributions
to low_γ / high_γ `T_d^(d_F)` reach −1.5, dwarfing the rest of the
cohort): jittered per-patient `T_d` for all three distances
side-by-side. Each dot is one patient, **green if `T_d < 0` (trace),
red if `T_d > 0`**; Pat_03 = orange triangle. Black bar = cohort
median; `n_neg / n_tot` annotated in a small badge at the top of
each column. Dashed gray line at zero. Footer legend explains all
markers.

### Fig S3 — Significance counts (stoplight heatmap)  →  `figS3_significance_counts.pdf`

A single 9 × 6 stoplight heatmap that reads at a glance. Rows = the
9 `(distance × phase-pair)` combinations grouped into three blocks
by phase pair (RPre→TT, RPre→RPost, TT→RPost), separated by black
lines. Columns = the 6 bands. Each cell shows `n_+/10` (number of
patients with `Z > 2` against the within-`rest_pre` split-half null)
with hard-threshold coloring:

- **green** (`n_+ ≥ 8/10`): cohort screen passes;
- **yellow** (`5 ≤ n_+ ≤ 7/10`): marginal;
- **red** (`n_+ ≤ 4/10`): cohort screen fails.

Pattern at sight: `d_F` is mostly green from α onward (the predicted
drift footprint); `d_S` is largely red except for high_γ; rank
distances on TT→RPost specifically are nearly all red — confirming
that this within-baseline screening is NOT the trace test (read
together with the verdict table in §3 and Fig 1).

---

## 6. Significance + Pat_06 sensitivity

Source: `Td_significance_tests.csv`. Per-cell Wilcoxon signed-rank,
one-sided, alternative `T_d < 0`.

### Top cells (sorted by uncorrected p)

| band × dist   | n_neg/10 | median T_d | Wilcoxon p | sign p | Wilcoxon p (no Pat_06) | Wilcoxon q (BH) |
|---|---:|---:|---:|---:|---:|---:|
| β × `d_F`     | 8/10 | −0.0715 | **0.042** | 0.055 | 0.082 | 0.29 |
| β × `d_S`     | 7/10 | −0.0377 | 0.053 | 0.172 | 0.102 | 0.29 |
| β × `d_P`     | 7/10 | −0.0130 | 0.065 | 0.172 | 0.125 | 0.29 |
| low_γ × `d_S` | 7/10 | −0.0744 | 0.065 | 0.172 | 0.102 | 0.29 |
| α × `d_S`     | 8/10 | −0.0420 | 0.116 | 0.055 | 0.180 | 0.35 |
| α × `d_F`     | 7/10 | −0.0670 | 0.116 | 0.172 | 0.180 | 0.35 |

**Zero cells survive BH-FDR at q = 0.05** across the 18 (band ×
distance) cells; smallest q ≈ 0.29.

**Pat_06 contributions:** Pat_06 contributes the largest `T_d^(d_F)`
values — low_γ = −1.48, high_γ = −1.86, β = −0.57 — order of
magnitude above the rest of the cohort. Removing Pat_06 pushes every
borderline cell out of even nominal significance.

**Methods sentence:**
> "Per-cell Wilcoxon signed-rank tests on the cohort `T_d` vector
> (one-sided, alternative `T_d < 0`) yielded a smallest uncorrected
> p ≈ 0.042 (β × `d_F`); zero cells survived Benjamini-Hochberg
> correction across the 18 (band × distance) cells (smallest q ≈
> 0.29). Removing the load-bearing patient (Pat_06) pushed every
> borderline cell out of nominal significance, so the headline is
> reported as a directional-consistency claim, not as a hypothesis
> test at α = 0.05."

---

## 7. Caveats — what NOT to claim

The honest list. Each item is a real reviewer concern.

1. **Not multiple-comparison significant.** Smallest BH-FDR q ≈ 0.29.
   Frame as "directionally consistent task-trace in 4 of 6 bands
   cohort-wide", not "significant trace in 4 of 6 bands".
2. **Pat_06 load-bearing.** Removing Pat_06 pushes every borderline
   cell out of nominal significance. Pat_06 also motivated the audit
   originally — disclose the circularity.
3. **n = 10 underpowered.** Median `|T_d^(d_S)|` per band is 0.01–0.07
   against cohort SD ~0.10–0.20; SNR ≈ 0.3–0.6 σ. Need n ≈ 30–50 for
   80 % power on a one-sided Wilcoxon.
4. **Three causes confusable** for `T_d < 0`: (a) genuine task-induced
   trace (what we want to claim), (b) anomalous `RPre` (warm-up /
   vigilance baseline) that both task and `RPost` differ from,
   (c) coupled monotone drift. The geometric positive control (Fig 2)
   **rejects (c)**; (a) vs (b) requires a within-`rest_post` split-half
   null (deferred).
5. **`d_F` triangle = drift fingerprint, not structure.** Report only
   as a diagnostic; can co-confirm `d_S`/`d_P` when those agree.
6. **Within-`rest_pre` null does NOT control for drift** — only for
   sampling jitter inside `rest_pre`. The owed controls
   (within-`rest_post` split-half null, cross-baseline distance,
   cohort drift-triangle sweep) are deferred.
7. **`V*` is band-dependent.** Different bands compute distances on
   different node sets per patient; `|V*|` for high_γ can drop to
   60–80 vs 100+ for δ/θ.
8. **Pat_03 inclusion is methodological.** Justified because Z and
   ratios are dimensionless w.r.t. fs; absolute distances on Pat_03
   are not directly comparable. Pat_03 turns out consistent with the
   cohort (negative `T_d^(d_S)` in 5/6 bands).
9. **Do NOT call `d_P` "volume + topology" or "orthogonal" to `d_S`.**
   Pearson on `triu(A)` is a magnitude-weighted linear correlation
   that mixes rank and amplitude — it does NOT decompose into
   "topology" + "volume", and cohort Spearman ρ between `T_d^(d_S)`
   and `T_d^(d_P)` is 0.85–0.95 per band. Use **"convergence"**.
10. **Do NOT call our finding the bare "persistence".** It is a
    **trace** — task reorganized AND change persists. The bare word
    flips the reading to "anchor" (module unchanged) on a careful
    reader. See `.agents/guides/01_project/terminology.md`.
11. **Band-averaging is conventional, not principled.** The substrate
    of the entire audit is `imcoh_abs(i, j) = mean_f |Im(S_ij(f)) /
    √(S_ii(f)·S_jj(f))|` — band-averaged absolute imaginary
    coherence. This is the standard EEG/MEG-FC summary, but it is a
    pragmatic compromise, not a first-principles estimator:
      - Information loss within the band: a sharp narrow-band peak
        (e.g. one bin in α) is diluted in proportion to the band's
        bandwidth. Band-averaging can only *hide* signal; it cannot
        fabricate it. So a null under band-averaging may hide a
        frequency-specific effect, but a positive finding is at
        least as strong as the band average suggests.
      - Jensen's inequality: `mean_f |x_f| ≠ |mean_f x_f|`, and the
        order of `|·|` and `mean_f` is a convention. We take `|·|`
        per-frequency-bin first then band-average, which prevents
        sign-flips of `Im(S_ij)` across the band from cancelling
        genuine phase-locking. Different conventions (e.g. `imcoh_sq
        = mean_f |x_f|²`) would give numerically different summaries
        of the same underlying coupling.
      - Treats frequencies inside a band as interchangeable, which
        physiologically they are not (e.g. low-α vs high-α generators).
      - Downstream LRG implication: the dendrogram, heat kernel and
        ultrametric all operate on the band-averaged adjacency, not
        on a specific frequency's coupling. Renormalization on a
        summary graph has no clean physical interpretation.

    **All trace findings in this report are statements about
    band-level FC**, not about coupling at any specific frequency
    inside the band. This substrate is what the entire pipeline,
    cache, and cohort definition is built on; replacing it would
    require rebuilding everything. We therefore disclose this as a
    bounded-scope caveat rather than a fixable defect, and treat
    frequency-resolved LRG as an optional follow-up if a result needs
    sub-band localisation.

---

## 8. Writing-agent guidance

### Lead with `d_S`, frame `d_F` as drift diagnostic

`d_S` is the load-bearing claim. Reference `d_P` as a **convergence
check** (β: 7/10 each), and reference `d_F` only as a drift
diagnostic — its triangle pattern is the predicted footprint of
monotone drift, not structural evidence on its own. Pre-empt the
reviewer question: "task↔post are temporally adjacent and pre↔post are
far apart, so monotone drift produces `T_d^(d_F) < 0` for free; the
rank distance `d_S` is invariant to this and is the load-bearing axis."

### Foreground the structural-vs-drift result (Fig 4)

The single strongest "this is not a drift artefact" argument is the
per-patient Spearman ρ between `T_d^(d_S)` and `T_d^(d_F)` = 0.78–0.92
per band, sign agreement 7–9 of 10. Give it a dedicated paragraph.

### Geometric positive control deserves prominence (Fig 2)

`d_S(TL, TT) ≈ 0.5 × d_S(RPre, RPost)` cohort-wide is the most
cohort-robust single fact in the audit. It rejects "everything is
drift" without needing a hypothesis test (it's about cluster widths,
not signed effect direction).

### Use ratios + absolute distances, NOT Z

Z = 12–685 is a tiny-MAD inflation artefact (Fig S1). Report
`d_obs / median(N)` (`cohort_scale_summary.csv`) or absolute distances
(`cohort_geometry_4phase_summary.csv`).

### Frame the headline as exploratory

Smallest BH-FDR q ≈ 0.29. "Directionally consistent task-trace in 4
of 6 bands cohort-wide", not "significant in 4 of 6 bands". Disclose
Pat_06 sensitivity and the deferred controls.

### What the LRG step is for

Frame the raw-FC trace as a hypothesis to be tested at higher
resolution by the LRG hierarchy, not as an established fact. The LRG
step should (i) confirm β/α/δ/low_γ on a multiscale axis that pools
edge evidence, and (ii) ideally surface θ or high_γ trace where the
raw-FC test was null/borderline.

---

## 9. Glossary

**Phases.** RPre = `rest_pre`, TL = `task_learn`, TT = `task_test`,
RPost = `rest_post`. Order in time: RPre → TL → TT → RPost.

**Cohort.** `PATIENTS_4PHASE = {Pat_02, 03, 05, 06, 07, 08, 10, 13,
14, 15}`, n = 10. Pat_03 = 1024 Hz outlier (others 2048 Hz).

**FC method.** `imcoh_abs` = band-averaged absolute imaginary
coherence (Nolte 2004 / Ewald 2012). Range [0, 1].

**Substrate.** `V*(p, b)` = per-patient, per-band giant-component
intersection of `|A^φ|` across the involved phases.
`triu(A, k=1)` = upper-triangular off-diagonal vector of A.

**Distances.**
- `d_S = 1 − ρ_S(triu A_a, triu A_b)` — rank-only.
- `d_P = 1 − ρ_P(triu A_a, triu A_b)` — magnitude-weighted complement.
- `d_F = ‖A_a − A_b‖_F / √(‖A_a‖_F · ‖A_b‖_F)` — amplitude-only.

**Trace scalar.** `T_d = d(TT, RPost) − d(RPre, TT)`. Negative ⇒
trace (RPost closer to TT than RPre is). `n_trace = #{p : T_d(p) <
0}` out of 10. (CSV column `n_persist` == this.)

**Statistics.** `Z = (d_obs − median(N)) / MAD(N)` — robust z on the
within-`rest_pre` split-half null `N`. **MAD** = median absolute
deviation. **`n_+`** = patients with `Z > 2`. **BH** =
Benjamini-Hochberg FDR. **Wilcoxon** = signed-rank test on the
cohort `T_d` vector, one-sided.

**Bands.** δ 1–4 Hz, θ 4–8, α 8–13, β 13–30, low_γ 30–80, high_γ
80–150.

**Module behaviour taxonomy** (see
`.agents/guides/01_project/terminology.md`): **trace** (changed in TT,
stays in RPost — our finding), **anchor** (unchanged across all
phases), **reset** (changed in TT, reverts in RPost), **emergent**
(did not exist in RPre).

---

## 10. Canonical artefacts (for source-tracing)

All files are at `data/audit/raw_fc_phase_distance/`.

- `final_verdict_table.csv` — one row per band: cohort medians, trace
  counts (column name: `n_persist`), verdict.
- `Td_per_patient_per_band.csv` — per-patient `T_d` for each distance.
- `Td_significance_tests.csv` — Wilcoxon + BH per (band, distance).
- `Td_dS_vs_dF_band_contrast.csv` — structural-vs-drift flag per band.
- `cohort_geometry_4phase_summary.csv` — full 4-phase geometry.
- `cohort_scale_summary.csv` — scale + ratio + Z per (band, distance,
  pair).
- `per_patient_with_scale.csv` — long-format per-patient scale CSV.
- `cohort_summary.csv` — audit_25 within-baseline `n_+` counts.
- `session_timing.csv` — per-patient phase durations + sampling rate.
- `Pat_NN/{report.pdf, audit.csv, distance_4phase.csv, null.npz}` —
  per-patient artefacts (n = 10).

**Canonical prose report:**
`.agents/reports/2026-04-29_result-1-raw-fc-phase-trace.md`.

**Bundle figure-rebuild script:**
`scripts/01_compute/audit/queries/q_bundle_figures.py` — regenerates
fig1/2/3 with the polished visuals shown here.
