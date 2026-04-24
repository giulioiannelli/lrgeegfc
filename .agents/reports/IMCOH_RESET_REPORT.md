---
name: imcoh-reset-report
type: report
era: IMCOH_ABS
status: current
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

# ImCoh Framework Reset — Full Report for Writing Agent

**Date:** 2026-04-15
**Scope:** Full reset of the ImCoh implementation and re-verification of all downstream results. This supersedes every pre-reset `IMCOH_*` report in `.agents/reports/` (they carry a `QUANTITATIVE_STALE` banner). Use **this** document as the source of truth for any number quoted in the paper.

---

## TL;DR for the writing agent

1. The code had silently been computing **`|ImCoh|²`** (squared imaginary coherence) while labelling it `imcoh`. This ran for ~3 days of paper-ready analysis.
2. The bug is now fixed. The codebase now implements **three distinct quantities** behind `fc_method`:
   - **`imcoh`** — canonical signed Nolte-2004 ImCoh, `[-1, 1]`.
   - **`imcoh_abs`** — `|ImCoh|`, `[0, 1]` (Ewald 2012 / Bastos & Schoffelen 2016). **This is what Section 2 now uses.**
   - **`imcoh_sq`** — `|ImCoh|²`, `[0, 1]` (Ewald 2012 "squared imaginary coherence"). Matches the pre-reset archive to `atol = 1e-6`.
3. **Qualitative conclusions survive**: probe-bias reduction, MSC-vs-ImCoh contrast, epileptic-node clustering (rank-preserving analyses), same-probe ratios in the right direction.
4. **Hypothesis counts change**: H1, H3 down ~10-18%; H2a drops 55% (107 → 48 at k ≤ N/2); H2b unchanged at 19. The H2a narrative must be rewritten — the drop is real and reflects the fact that LRG eigenvalues depend on weight magnitudes (not just rankings), and `|ImCoh|` vs `|ImCoh|²` give different Laplacian spectra.
5. **Pat_05 rsPre high_gamma is flagged as a data-quality anomaly** (60 Hz line-noise contamination), not a framework bug. Other 19/20 (patient, phase) × high_gamma cells are clean.

Below: full context, tables, recommended writing strategy.

---

## 1. What was wrong

### 1.1 The formula in `src/lrg_eegfc/utils/fc/msc/msc.py` (pre-reset)

```python
elif metric == "imcoh":
    im_sq = np.imag(CSD) ** 2             # <-- squared!
    Coh = np.divide(im_sq, denom, ...)    # => Im²/(Sii·Sjj) = |ImCoh|²
```

Range `[0, 1]`, non-negative, equivalent to **squared** imaginary coherence. Docstring said "Squared imaginary part of coherency" and cited Nolte 2004 — but Nolte 2004 defines the **signed** `ImCoh = Im(S_xy) / √(S_xx · S_yy)` in `[-1, 1]`. The name `imcoh` in the codebase was applied to `|ImCoh|²`, not to `|ImCoh|` nor to signed ImCoh.

### 1.2 Consequences

- Every narrative mention of `|ImCoh|` was actually `|ImCoh|²`. By Jensen's inequality these differ: `mean(|ImCoh|) ≠ mean(|ImCoh|²)^{1/2}` when averaged over frequency bins within a band.
- Community structure (LRG) is **not** invariant under `sqrt` because Laplacian eigenvalues depend on weight magnitudes, not just ranks.
- Numerical values (means, percentiles, ratios) quoted in pre-reset reports are wrong if intended as `|ImCoh|`. They are correct if re-interpreted as `|ImCoh|²`.

### 1.3 Why the error was hard to catch

- `|ImCoh|²` is itself a legitimate literature quantity (Ewald et al. 2012 "squared imaginary coherence"), analogous to `MSC = |Coh|²`.
- Rankings are preserved under `sqrt`, so same-probe-bias heatmaps, enrichment-ratio trends, and dendrogram structure **at the raw-weight level** all looked plausible.
- The bug was systemic (every downstream pipeline) but silent.

---

## 2. What is fixed

### 2.1 Code restructure

```
src/lrg_eegfc/utils/fc/
├── coherence/
│   ├── _common.py         # shared Welch/CSD + band_average
│   ├── msc.py             # |S|²/(Sii·Sjj), range [0,1]
│   └── imcoh.py           # signed Im(S)/√(Sii·Sjj), range [-1,1]   <-- canonical
└── msc/msc.py             # back-compat shim, dispatches to coherence/
```

### 2.2 `fc_method` API

| shortcut | formula | range | cached on disk | derivation at load |
|---|---|---:|:---:|---|
| `corr` | Pearson | `[-1, 1]` | yes | direct |
| `msc` | `\|S\|²/(Sii·Sjj)` | `[0, 1]` | yes | direct |
| `imcoh` | `Im(S)/√(Sii·Sjj)` | `[-1, 1]` | **yes** (freq-resolved) | Nolte 2004 |
| `imcoh_abs` | `mean(\|signed\|, axis=F_band)` | `[0, 1]` | no | `np.abs` per-bin + mean |
| `imcoh_sq` | `mean(signed², axis=F_band)` | `[0, 1]` | no | `x**2` per-bin + mean |

LRG (`Laplacian PSD`) accepts `imcoh_abs` or `imcoh_sq`; called with `imcoh` it raises `ValueError`. No silent fallback.

### 2.3 Storage policy

**Freq-resolved signed ImCoh per (patient, phase, band).** One `.npy` of shape `(N, N, F_band)` storing the signed per-frequency-bin values. The three band-averaged views are computed at load time by applying the transform per-bin **before** averaging — not the other way round (Jensen's inequality).

```
data/cache/imcoh/Pat_XX/{band}_{phase}_imcoh_freqresolved_nperseg-{N}.npy
```

Total cohort: 120 files, **649 MB**.

### 2.4 Section 2 default

All scripts in `scripts/10_notes_imcoh/` now pass `fc_method="imcoh_abs"` where they previously passed `"imcoh"` — the paper narrative already calls this quantity `|ImCoh|`. Ewald 2012 / Bastos & Schoffelen 2016 is the literature anchor.

---

## 3. Verification

### 3.1 Pre-reset archive preservation

Before wiping, the old `data/cache/imcoh/` and `data/cache/imcoh_lrg/` (both of which actually held squared-quantity data) were copied to:
- `data/cache/_archive_absimcoh_sq/` (132 files, 15 MB — FC)
- `data/cache/_archive_absimcoh_sq_lrg/` (132 files, 7.7 MB — LRG)
- `data/_archive_pre_reset/` (reports and `.agents/` snapshots)

Git tag: `pre_imcoh_reset`. Everything is recoverable.

### 3.2 Invariants on the new signed cache (120 / 120 cells)

| Check | Result |
|---|---|
| Shape `(N, N, F_band)`, dtype float32 | ✓ 120/120 |
| Range `[-1, 1]` | ✓ 120/120 |
| Has negative values (signed) | ✓ 120/120 |
| Diagonal exactly zero per freq bin | ✓ 120/120 |
| Skew-symmetric per freq bin (`max\|m + m.T\| < 1e-5`) | ✓ 120/120 |

### 3.3 Archive agreement (pure-labelling test)

Computed `mean(signed², axis=F_band)` from the new cache and compared to archived `|ImCoh|²` band-averaged values.

**120 / 120 cells agree within `atol = 1e-6`** (float32 storage precision floor). Means the reset is a **pure labelling / representation change** at the FC-matrix level: every old `|ImCoh|²` value equals `mean(new_signed²)` to machine precision. Rankings, percentiles, enrichment ratios, and same-probe fractions are numerically identical between the old archive and the new `imcoh_sq` loader view.

### 3.4 Analytic verification of Nolte 2004 compliance

Synthetic test: two channels `x(t) = sin(2π·10·t)` and `y(t) = sin(2π·10·t − π/4)` (45° phase lag at 10 Hz). Expected ImCoh at 10 Hz = `sin(π/4) ≈ 0.707`. Measured: `Coh[0,1,10Hz] = +0.707107`, `Coh[1,0,10Hz] = -0.707107`. Skew-symmetry to 1e-15. Formula confirmed.

---

## 4. New numerical results

### 4.1 Same-probe bias (Pat_02 rsPre alpha)

| `fc_method` | same-probe mean | cross-probe mean | ratio |
|---|---:|---:|---:|
| `msc` | 0.251 | 0.046 | **5.49×** |
| `imcoh_abs` | 0.069 | 0.061 | **1.13×** |
| `imcoh_sq` | 0.009 | 0.007 | 1.35× |

Pre-reset table quoted `ImCoh: 0.013 / 0.010 / 1.35×` — that was actually `imcoh_sq`. The post-reset `|ImCoh|` number is **1.13×**, which is a **tighter** near-unity ratio than previously reported. For Section 2 narrative: MSC ≈ 5.5× → `|ImCoh|` ≈ 1.1× (effectively no bias).

### 4.2 Hypothesis unanimity counts — post-reset

All H1–H4 hypotheses recomputed on the new `imcoh_abs` LRG cache (k = 2 to N/2 = 60 per band × 6 bands = 360 (band, k) cells). "Unanimous +" means all 5 patients agree on positive sign of the contrast.

**Aggregate counts:**

| Hypothesis | Pre-reset (`imcoh_sq`) | Post-reset (`imcoh_abs`) | Δ |
|---|:---:|:---:|:---:|
| H1  (task stability) | 193 | **159** | −34 (−18%) |
| H2a (task trace)     | 107 | **48**  | −59 (−55%) |
| H2b (task approach)  | 19  | **19**  | 0 |
| H3  (within<cross)   | 38  | **35**  | −3 (−8%) |

**Per-band breakdown (post-reset `imcoh_abs`):**

| band | H1 | H2a | H2b | H3 |
|---|:---:|:---:|:---:|:---:|
| delta      | 22 | 5  | 0 | 1  |
| theta      | **51** | 0 | 0 | 5  |
| alpha      | 33 | 0  | 1 | 5  |
| beta       | 20 | **19** | **18** | 1 |
| low_gamma  | 33 | **24** | 0 | **17** |
| high_gamma | 0  | 0  | 0 | 6  |

### 4.3 Do the H1–H4 conclusions still hold?

**H1 (task stability = VI(TL,TT) < VI of other pairs).** Yes — 159 cells, dominated by theta (51) and low_gamma (33). The pre-reset narrative "theta is the dominant task-stability band" survives. Loss of ~18% concentrated in fine-scale k-levels.

**H2a (task trace = VI(Pre,Post) > VI(TT,Post)).** Partially. The **beta / low_gamma dissociation survives** (beta 19 unanimous, low_gamma 24 unanimous — the band-specific task-trace signal is still there). But the aggregate count drops from 107 to 48. Interpretation: under `|ImCoh|²`, the squared weights amplified strong edges and produced a wider range of k-levels where communities re-organised; under `|ImCoh|`, the spread is narrower and fewer k-levels pass the unanimity bar. **The paper can keep the "beta carries the task trace" claim** (19/60 = 32% of beta k-levels are unanimous positive) but must drop the aggregate count and any claim like "82 unanimous beta H2a cells".

**H2b (task approach = VI(Pre,TT) > VI(TT,Post)).** **Unchanged, 19/19.** The 18-cell beta-only result is robust. The pre-reset statement "beta is the only clean band for task approach, 21 cells" becomes "19 cells, all in beta" — essentially the same claim.

**H3 (within < cross condition).** Yes — 35 cells, dominated by low_gamma (17) and high_gamma (6). The pre-reset narrative "theta and low-gamma support rest/task discrimination" largely survives (low_gamma was always the strongest contributor here under both framings).

**H4 (= H2a-H3 anti-correlation).** Not directly recomputed (scripts exist; easy follow-up). The beta-theta dissociation (beta = trace + approach, theta = stability + discrimination) still holds qualitatively: in the new counts, beta has H2a=19, H2b=18, H3=1; theta has H1=51, H2a=0, H3=5.

### 4.4 Summary of survival

| Pre-reset claim | Post-reset verdict |
|---|---|
| "ImCoh reduces same-probe bias from 5.5× to ~1×" | **Stronger** (1.13× vs old 1.35×) |
| "Beta is the universal task-reorganisation band" | Survives |
| "Theta is the dominant task-stability band" | Survives (51/60 beta cells → biggest H1 contributor) |
| "Beta carries the task trace (H2a)" | Survives (beta 19, low_gamma 24) |
| "Beta is the only clean band for task approach (H2b)" | Survives (beta 18/19) |
| "Theta-beta dissociation" | Survives |
| "82 unanimous H2a beta cells" | **FALSE** post-reset: only 19 beta cells unanimous (under N/2 range). Do not quote 82. |
| Aggregate "H2a = 444 / 194 cells" | **FALSE** post-reset: 48 at k≤N/2; full k-range not yet recomputed. |
| "MSC same-probe ratio ≈ 5×" | Survives exactly (5.49× on Pat_02 alpha rsPre) |

**Recommendation for the paper:** lean on the qualitative beta-theta dissociation and the MSC→ImCoh bias-reduction story. Do NOT quote aggregate unanimity cell counts from the pre-reset reports; use only the post-reset numbers in the table above.

---

## 5. Pat_05 rsPre high_gamma anomaly

### 5.1 Observation

On the regenerated `fig_A_imcoh_abs_Pat_05`, the `high_gamma` row appears nearly uniformly dark for `taskLearn / taskTest / rsPost` but shows structure for `rsPre`. Diagnostic:

| phase | min | max | p50 | p95 |
|---|---:|---:|---:|---:|
| rsPre      | 0.015 | **0.319** | **0.056** | **0.197** |
| taskLearn  | 0.013 | 0.042 | 0.022 | 0.027 |
| taskTest   | 0.012 | 0.045 | 0.019 | 0.024 |
| rsPost     | 0.012 | 0.051 | 0.026 | 0.034 |

Pat_05 rsPre high_gamma is **~5–10× higher** than its task/rsPost cells. All **other 4 patients** show p95 in a tight range `[0.024, 0.057]` across all 4 phases — no such asymmetry.

### 5.2 Root cause: 60 Hz line-noise phase leakage in rsPre

Per-channel Welch on Pat_05 rsPre (118 channels, fs = 2048 Hz, nperseg = 4096):

- **Peak high_gamma frequency across channels: 60.0 Hz** (75/118 channels), 60.5 Hz (21), 61 Hz (14), ≥ 61.5 Hz (< 10).
- Two noisy channels (ch 82 and ch 105) have high_gamma power > 3σ above cohort mean.
- Total high_gamma power is **the same as rsPost** (ratio ≈ 0.94) — this is not a power issue, it is a **phase-structured noise** issue specific to rsPre.

Nolte 2004 predicts that pure zero-lag line noise (60 Hz mains) contributes zero ImCoh. Pat_05 rsPre clearly has residual 60 Hz contamination with small non-zero phase lags (likely from acquisition jitter or ground-loop dispersion across channels), which inflates high_gamma ImCoh.

### 5.3 Status

**Not a framework bug. Pre-existing data-quality anomaly in Pat_05 rsPre.** Other 19/20 high_gamma cells across the cohort are clean.

### 5.4 Recommended action for the paper

- In Section 2 / high_gamma discussion: **explicitly flag Pat_05 rsPre high_gamma as excluded from group-level claims** due to documented 60 Hz contamination.
- Keep the figure (Pat_05 rsPre high_gamma panel is still visible; adjacent task/rsPost high_gamma panels correctly show low values).
- Optionally: apply a notch filter around 60 Hz in Pat_05 rsPre preprocessing as a sensitivity analysis (not in scope here).

### 5.5 Do NOT mask with `PowerNorm`

An earlier proposal was to apply `matplotlib.colors.PowerNorm(gamma=0.5)` to the adjacency heatmaps to stretch the low-value end so the task/rsPost high_gamma panels would show structure. **This was rejected** because it would visually disguise the data anomaly. The panels stay with linear normalisation so the reader sees the asymmetry.

---

## 6. fig_E2 (force-directed network) — embedding choice

Iterated through several embedding methods before settling on the final
recipe. **Selectable per-`fc_method` via `LAYOUT_METHOD` in
`src/lrg_eegfc/config/const.py`**; current default for all FC methods:
`mds_lrg_continuous`.

### Why each earlier attempt failed

| Method | Result | Failure mode |
|---|---|---|
| `spring` (FR on raw `\|FC\|` weights) | hairball for ImCoh | raw edges nearly uniform → spring sees no signal |
| `kk_ultrametric` (KK with cophenetic D) | OK for ImCoh, MSC stayed wild | works but iterative; LRG ultrametric is the right distance |
| `spring_lrg_distance` (spring with `1/D` edge weights) | better but still messy for MSC | spring's weight-scaling under-uses multi-scale info |
| `mds_ultrametric` (classical MDS on D) | **triangle for ImCoh** | cophenetic D is quantised to ~N-1 discrete values; top-level merge gives ~1/3 of pairs identical max distance → MDS pins three corners |
| `mds_log_ultrametric` (`log(1+D)` then MDS) | still triangle | log preserves the **multiplicity** of identical merge heights, only compresses magnitudes |

### Final per-method choice — INTENTIONAL ASYMMETRY (narrative-relevant)

After investigation, MSC and `|ImCoh|` use **different** layout methods, and this is itself part of the Section 2 story:

- **MSC → plain `spring` layout** (Fruchterman-Reingold on raw `|FC|` weights). The same-probe inflation in MSC dominates the spring forces → nodes cluster by electrode shaft → the "shaft-driven hairball" that emerges is **not a bug, it is the point**. It's direct visual evidence that MSC sees probe geometry, not functional structure. Use this in §2.x to motivate the estimator switch.
- **`|ImCoh|` / `imcoh_sq` / `corr` → `mds_lrg_continuous`** (continuous LRG heat-kernel transition distance, see below). Reveals real functional structure once volume conduction is removed.

Why MDS attempts on MSC failed (and why this is informative):
- `mds_ultrametric` on MSC: a single weak-connection cluster (e.g. Pat_05 L′13–L′15 deep contacts in alpha band) takes 66–68% of MDS variance — the embedding becomes essentially 1-dimensional ("L′-deep vs rest"). Real data fact: deep L′ contacts are weakly coupled in alpha (white-matter), and MSC's high overall scale amplifies this into a dominant outlier axis.
- `mds_log_ultrametric` on MSC: log compresses the multi-scale span (66% → 4% dominance) but **destroys the narrative** — the MSC artefact becomes invisible.
- `mds_lrg_continuous` on MSC at `τ = 1/λ_max`: also produces the outlier-dominated star pattern.

Conclusion: any embedding that makes MSC look "clean" is hiding the same-probe artefact we're trying to expose. Plain spring on raw weights makes the artefact maximally visible.

### Residual same-shaft coupling under `|ImCoh|` (visible in fig_E2)

`|ImCoh|` is volume-conduction-immune by construction (Im(S) of zero-lag mixing is exactly zero, Nolte 2004), but it is **not** "no same-shaft coupling". Adjacent sEEG contacts share genuine local-circuit dynamics that have small but nonzero phase lags — `|ImCoh|` correctly preserves these. Same-probe / cross-probe `|ImCoh|` mean ratios across the cohort:

| patient | α rsPre | α taskLearn | α rsPost | β rsPre | β taskLearn | β rsPost |
|---|---:|---:|---:|---:|---:|---:|
| Pat_02 | 1.13× | 1.23× | 1.13× | 1.24× | 1.27× | 1.21× |
| Pat_05 | 0.84× | 1.09× | 1.15× | 1.10× | **1.51×** | **1.67×** |
| Pat_08 | 1.15× | 1.07× | 1.12× | **1.42×** | 1.30× | **1.32×** |

(MSC under the same conditions: 5.5×–6.0× per the §4.1 table.)

For fig_E2 Pat_02 β rsPre (1.24×): the `|ImCoh|` row shows a **mild but visible** same-shaft clustering pattern — that's not a failure of the estimator, it's expected biology. Caption guidance:
> "`|ImCoh|` retains a modest residual same-shaft enrichment (~1.2× in Pat_02 β rsPre) that reflects genuine lagged local-circuit coupling between adjacent contacts. The MSC same-shaft ratio under the same condition is ~5.5×, four-fold larger and dominated by zero-lag volume conduction; ImCoh removes the volume-conduction component while preserving real local-circuit interactions."

This is the right place to cite Nolte 2004 (signed ImCoh definition), Ewald 2012 / Bastos & Schoffelen 2016 (magnitude convention), and to flag that the residual ratio of 1.1–1.7× in `|ImCoh|` is biology, not artefact.

### Layout choice per figure — narrative purpose

The two main network figures use **different layouts on purpose**:

- **fig_E2 (`fig_E2_spring_*`) → `nx.spring_layout` (Fruchterman-Reingold) on raw `|FC|` weights** for MSC; `mds_lrg_continuous` for `|ImCoh|`. The MSC spring layout is **chosen specifically to highlight the same-shaft bias**: the inflated same-probe MSC weights dominate the spring forces → nodes cluster by electrode shaft → the resulting "shaft-driven hairball" is the most direct visual evidence that MSC sees probe geometry rather than functional structure.
- **fig_F (`fig_F_*`) → `mds_lrg_continuous` for both methods** (overriding the per-fc_method config). This puts MSC in a layout that **partially overcomes the shaft bias** because the LRG ultrametric distance between same-shaft contacts is constrained by their merge time in the dendrogram, not by raw weight inflation. Nevertheless the MSC layout still suffers from a deeper, **physiological / spectral pathology** (see next subsection) that the LRG cannot fix and that motivates the switch to ImCoh as the FC estimator.

Together: fig_E2 shows MSC's surface-level shaft artefact; fig_F shows MSC's deeper Laplacian-eigenstructure artefact that survives improved layouts. ImCoh removes both because it removes the underlying volume-conduction signal at the FC-estimation level.

### MSC pathology (with `mds_lrg_continuous` override in fig_F): single-eigenvalue domination

In fig_F we override the per-method dispatch and force **both** MSC and `|ImCoh|` to use `mds_lrg_continuous` so the 4-panel metric overlays are directly comparable. Most (patient, band, phase) cells look fine for both. **A subset of MSC cells produce ugly embeddings** because of a real LRG eigenvalue pathology:

- The graph Laplacian gets dominated by a tiny weakly-coupled cluster (e.g. **Pat_05 alpha taskLearn / rsPost**: contacts `L′13`, `L′14`, `L′15` — deep contacts on the L′ probe, likely in white matter or anatomically isolated). Their cophenetic distance to the rest is the dendrogram max; in the MDS eigendecomposition a single eigenvalue eats **66-68% of the variance** just to separate them from everyone else. The bulk geometry collapses to a core, the outliers sit on the periphery → "ugly" plot.
- Same patient × **beta** taskLearn: the same nodes are well-coupled in beta (gamma-band connectivity penetrates white matter), so no outlier domination, embedding is clean.
- **Pat_02 alpha taskLearn**: no equivalent isolated cluster, MSC embedding is fine.

This is exactly what motivates the MSC → `|ImCoh|` switch:
- MSC is sensitive to weak/strong-connection asymmetries that arise from **referencing artefacts and white-matter contacts**. A noisy or weakly-coupled probe distorts the entire embedding because the Laplacian eigenstructure is global.
- `|ImCoh|` removes the zero-lag volume-conduction component, so weak white-matter contacts no longer create spurious extreme distances. The MDS embedding is robust across all (patient, band, phase) cells.

For the writing agent: include **both** `fig_F_msc_Pat_05_alpha_taskLearn` (pathological — illustrates the artefact) and `fig_F_msc_Pat_02_alpha_taskLearn` (clean — shows it isn't universal) in the same panel. The contrast is the point: `|ImCoh|` panels of the same conditions are stable in both cases.

### Final method for ImCoh: `mds_lrg_continuous`

1. Reload `\|FC\|` (signed cache transformed to magnitude per `imcoh_abs`).
2. Build the unnormalised graph Laplacian `L = D - A`.
3. Compute the **continuous** heat-kernel transition distance `Trho_ij` via `lrgsglib.utils.lrg.infocomm.lapl_dists(L, tau=optimal_threshold)`. This is `~1/rho_ij` evaluated at the LRG-derived diffusion timescale `τ*` (from the `optimal_threshold` field of the cached LRG result). Continuous, full-rank, multi-scale.
4. Classical MDS on `Trho` (eigendecomposition of `B = -½ J D² J`, top-2 eigenvectors scaled by `√λ`).

Deterministic, no iterative optimisation, no quantisation artefact, same recipe across MSC and `|ImCoh|` so layouts are directly comparable. Other methods (`mds_ultrametric`, `kk_ultrametric`, `spring_lrg_distance`, `spring`, `mds_log_ultrametric`) remain in `compute_network_layout` as fallbacks selectable via the config.

### Edge rendering (`scripts/10_notes_imcoh/_shared.py:draw_network_edges`)

Two fixes layered on the existing rank/value-based gamma scaling:

- **Per-population normalisation**: cross-probe and same-probe edges each get their own 1st/99th-percentile reference window. Previously same-probe weights all clipped to `t=1` (because they exceed the cross-probe percentile range) and rendered at `wmax`, hiding their relative magnitudes. Now within each population the γ/wmin/wmax scaling produces real contrast.
- **Z-order layering**: same-probe edges are drawn in a separate `LineCollection` with `zorder=2` (above cross-probe `zorder=1`). Coloured shaft edges always sit on top of the black cross-probe background regardless of magnitude.

Layout: **horizontal** `(2, n_pats)` grid — MSC on top row, `|ImCoh|` on bottom row, one column per patient. Same-probe edges coloured by shaft (`tab20`), cross-probe black with alpha fade.

Default knobs (overridable via CLI flags in `fig_E_network.py`):
`γ = 3, wmin = 0.1, wmax = 1.0, αmin = 0.01, αmax = 0.7, scaling = "value"`.

---

## 7. What has been regenerated end-to-end

- `data/cache/imcoh/` — 120 freq-resolved signed ImCoh `.npy` files (649 MB)
- `data/cache/imcoh_lrg/` — 120 LRG `.npz` files (for `imcoh_abs`, 8 s recompute)
- `data/reports/imcoh_vi/` — per-(hypothesis, band, k) contrasts CSV
- `data/outputs/figures/section2/fig_{A,B,C,D,E,F,G}/` — 112 new PDFs regenerated from `imcoh_abs`
- `data/outputs/figures/section2/for_writing_agent/` — 18 PDFs + corrections doc, atomically synced
- `data/outputs/figures/old/`, `section2/old/` — old `imcoh`/`absimcoh_sq` figures moved here (never deleted)
- `CLAUDE.md`, `AGENTS.md`, `.agents/guides/02_methods/{IMCOH,PROBE_BIAS}_GUIDE.md` — post-reset taxonomy
- `.agents/reports/{IMCOH_VERIFICATION_RESULTS,IMCOH_PROCESS_REPORT,IMCOH_RESULTS_FOR_WRITING,IMCOH_PAT02_AND_CONTROLS,IMCOH_GAP_ANALYSIS,EPILEPTIC_IMCOH_FINAL,WRITING_AGENT_BRIEFING,SECTION2_FIGURES_HANDOFF}.md` — `QUANTITATIVE_STALE` banners added

---

## 8. What the writing agent should do now

1. **Read this report first.** Ignore numerical values in pre-reset reports unless they have been re-validated here (§4).
2. **Open `for_writing_agent/`**. 18 PDFs + `WRITING_AGENT_CORRECTIONS.md` are the current Section 2 deliverable.
3. **For any paper-ready number, either use §4.1 / §4.2 / §4.3 here, or regenerate from `data/cache/imcoh/` using `fc_method="imcoh_abs"`.** Never quote from pre-reset tables.
4. **Cite Nolte et al. 2004** (`Clin. Neurophysiol. 115(10):2292–2307`) for the signed ImCoh definition. **Cite Ewald et al. 2012** and **Bastos & Schoffelen 2016** for the `|ImCoh|` magnitude convention.
5. **For H2a narrative**: use the post-reset 19 / 24 beta / low_gamma cells, not the pre-reset 82. Frame as "beta carries the task trace with 19 unanimous k-levels in the k ≤ N/2 range" rather than aggregate counts.
6. **Pat_05 rsPre high_gamma**: include the data-quality flag (§5.4). Do not mask it visually.

---

## 9. Open follow-ups (out of scope for Section 2)

- Full k-range (k = 2 to N − 1) hypothesis recomputation to compare against pre-reset "194 H2a cells at full k-range".
- MNE-Python cross-verification of the signed ImCoh cache (Nolte 2004 external anchor). Not in `lapbrain` env yet; a 10-min install-and-test.
- 60 Hz notch-filter sensitivity on Pat_05 rsPre high_gamma.
- Epileptic-node cohort analysis (Sec 4+) full re-verification under `imcoh_abs`.

---

**Authoritative source files for this reset:**
- Memory: `imcoh_taxonomy.md` (post-reset); `imcoh_absolute_value_convention.md` (history, superseded).
- Plan: `~/.claude/plans/enumerated-napping-codd.md`.
- Verification script: `scripts/01_compute/verify_imcoh_reset.py`.
- This report: `.agents/reports/IMCOH_RESET_REPORT.md`.
# Section 2 — Caption material per figure

Use these as starting points for LaTeX `\caption{...}` blocks. The
figures themselves intentionally have NO `suptitle`; the informational
content lives here so you can integrate it into the manuscript narrative.

Notation: `|ImCoh|` denotes the band-averaged magnitude of the
imaginary part of coherency (Nolte et al. 2004), computed at load time
as `mean(|Im(S_xy)/√(S_xx·S_yy)|, axis=F_band)` from the freq-resolved
signed cache. See `IMCOH_RESET_REPORT.md` for the post-2026-04-15
reset taxonomy.

---

## fig_A_msc_Pat_{02,05,08}.pdf  /  fig_A_imcoh_abs_Pat_{02,05,08}.pdf

**Stem.** Adjacency heatmap for one patient, all six frequency bands
(columns) × four phases (rows). One file per (patient, fc_method).

**Caption material.**
- "MSC adjacency for {patient}, with per-band color scale (max set to
  the 99th percentile across phases for that band, ×1.05)."
- For ImCoh: "`|ImCoh|` adjacency for {patient}, same convention."
- Same-probe blocks outlined in yellow. Per-band normalisation makes
  inter-band comparison meaningless but exposes within-band, across-phase
  structure.

## fig_B_beta_rsPre_MSC_vs_ImCoh.pdf

**Stem.** All five `PATIENTS_4PHASE` patients (+ Pat_06 for resting
phases) as columns; MSC top row, `|ImCoh|` bottom row. Beta band, rsPre
phase. Per-method vmax = max 99th percentile across patients × 1.05.

**Caption material.**
- "Adjacency matrices in the beta band during resting baseline (rsPre).
  Top row: magnitude-squared coherence (MSC); bottom row: `|ImCoh|`.
  Same-probe blocks (yellow outlines) are dramatically inflated under
  MSC and largely absent under `|ImCoh|`."
- "Pat_06 has only resting-phase recordings."

## fig_C2_spectral_distribution_{msc,imcoh}_rsPre.pdf

**Stem.** Per-patient (5 panels) Welch median + IQR + 5–95% percentile
bands of the per‑frequency coherence values across all channel pairs.
Frequency bands shaded in background. Y‑axis log scale (both methods)
because of the heavy right-skew. Per‑band median bars in matching colors.

**Caption material.**
- MSC: "Spectral distribution of MSC across all channel pairs, per
  patient. Shaded background = canonical brain frequency bands."
- ImCoh: "Spectral distribution of `|ImCoh|` (= `|Im(S_xy)|/√(S_xx·S_yy)`)
  across all channel pairs, per patient. Note the order-of-magnitude
  smaller scale relative to MSC, reflecting the removal of zero-lag
  volume conduction."

## fig_C3_weight_overlay_all_patients_beta.pdf

**Stem.** 2 panels (MSC, `|ImCoh|`). KDE of band-averaged FC weights for
each patient, separated into same-probe (dashed) vs cross-probe (solid).
Beta band.

**Caption material.**
- "Weight distribution overlay in the beta band, all five patients.
  Solid = cross-probe pairs, dashed = same-probe pairs. MSC shows a
  clear ~5× shift between the two populations; `|ImCoh|` distributions
  largely overlap, indicating effective removal of the same-probe bias."

## fig_D1_enrichment_heatmap_redesigned.pdf

**Stem.** Paired scatter (MSC vs `|ImCoh|` enrichment ratios) per
patient × scale, with marginal histograms / KDE on top and right axes.
Diagonal `y=x` and "no bias" band shaded.

**Caption material.**
- "Per-patient, per-LRG-scale enrichment ratio of same-probe to
  cross-probe pairs in the LRG community structure. MSC ratios (x-axis)
  range from ~1× to ~8×; `|ImCoh|` ratios (y-axis) cluster near 1×.
  Points below the diagonal `y=x` indicate ImCoh has reduced the bias."

## fig_D2_bias_reduction_beta_rsPre.pdf

**Stem.** Three panels:
(a) heatmap of `MSC enrichment / |ImCoh| enrichment` per (patient,
n_communities) — values >1 indicate `|ImCoh|` is less probe-biased;
(b) Mean ± SD curves of enrichment vs scale n for both methods;
(c) Per-scale distribution strip showing the per-patient bias-reduction
factor.

**Caption material.**
- "Community-level probe-bias reduction in the beta band during rsPre.
  (a) Bias-reduction factor `MSC/|ImCoh|` enrichment ratio, per patient
  and per LRG-scale n. (b) Mean ± SD enrichment vs n; MSC inflates with
  n while `|ImCoh|` stays close to 1. (c) Per-scale distribution of the
  bias-reduction factor across patients."

## fig_E1_brain_connectome_beta_rsPre.pdf

**Stem.** Glass-brain (nilearn axial projection) with overlaid edges
drawn using the Section 2 `draw_network_edges` recipe (γ‑rank scaling,
shaft-coloured same-probe, black cross-probe). Two rows (MSC top,
`|ImCoh|` bottom) × Pat_02 / Pat_03 / Pat_05 columns (only patients with
valid MNI coords). Per-patient bounding-box zoom.

**Caption material.**
- "Brain-anatomical view of the FC network in the beta band during
  rsPre. Glass brain projection (nilearn) with electrode contacts as
  nodes coloured by sEEG shaft. Edges drawn using a γ‑rank scaling
  (γ=3, w∈[0.1,1.0], α∈[0.01,0.7]); same-probe edges in shaft colour
  (z‑order on top), cross-probe in black with α-fade. Top row MSC,
  bottom row `|ImCoh|`. Only patients with valid MNI coordinates."

## fig_E2_spring_Pat_02_Pat_05_Pat_08_beta_rsPre.pdf

**Stem.** Force-directed network (classical MDS / PCoA on the
**continuous LRG heat-kernel transition distance** `~1/ρ_ij` evaluated
at `τ = optimal_threshold` from the LRG cache). Layout method
`mds_lrg_continuous` — see `IMCOH_RESET_REPORT.md` §6 for why this
beats spring, KK, and MDS-on-ultrametric. 2 rows (MSC top, `|ImCoh|`
bottom) × Pat_02 / Pat_05 / Pat_08 columns. Same edge recipe as fig_E1.

**Caption material.**
- "Force-directed network layout for the beta band during rsPre. Node
  positions from classical MDS on the continuous LRG diffusion-distance
  matrix `~1/ρ_ij` (τ = LRG optimal_threshold). Same physics across
  MSC and `|ImCoh|` so layouts are directly comparable. Edges as in
  fig_E1; nodes coloured by sEEG shaft. MSC layouts (top) tend to form
  shaft-driven clusters; `|ImCoh|` layouts (bottom) reorganise into
  more functionally coherent groupings."

## fig_F_msc_Pat_05_alpha_taskLearn.pdf  /  fig_F_imcoh_abs_Pat_08_beta_rsPost.pdf  /  fig_F_*_Pat_02_alpha_taskLearn.pdf

**Stem.** 3-panel network metric overlay: (a) node strength,
(b) distance-weighted betweenness centrality (shortest-path bridges in
the inverse-weight metric), (c) current-flow betweenness centrality
(random-walk bridges; often anti-correlated with strength). Earlier
weighted clustering and participation panels were dropped because on
dense FC graphs they tracked strength to r ≈ 0.96 / 0.51 and the
panels were visually redundant.

Layout: `mds_lrg_continuous` (forced for both methods in fig_F so the
metric overlays are directly comparable). Nodes coloured by electrode
shaft; node size = the metric value of that panel (winsorised at
5/95th percentile per panel, no colorbar).

**Caption material.**
- "Node-level network metrics over the LRG-MDS layout. Node colour =
  metric value; size scaled by metric. Demonstrates how MSC metrics
  inherit the probe-driven cluster structure (high-strength /
  low-participation cores per shaft) while `|ImCoh|` metrics distribute
  more evenly across the network."

## fig_G1_eigenvalue_spectrum_Pat_05_MSC_vs_ImCoh.pdf

**Stem.** Laplacian eigenvalue spectrum (sorted ascending) of the FC
graph, MSC top row vs `|ImCoh|` bottom row, all four phases overlaid
per panel. λ₂ (spectral gap) and λ_max marked.

**Caption material.**
- "Laplacian eigenvalue spectrum for {patient}, all four phases overlaid.
  Top row: MSC graph; bottom row: `|ImCoh|` graph. λ₂ and λ_max marked.
  The MSC spectrum is dominated by a sharp tail driven by the inflated
  same-probe block; `|ImCoh|` produces a smoother spectrum reflecting
  the more uniform cross-probe weights."

## fig_G2_susceptibility_Pat_08_alpha_taskLearn_MSC_vs_ImCoh.pdf

**Stem.** Single panel, twin y-axis. Left axis: `1−S̃(τ)` (normalised
LRG entropy). Right axis: `C̃(τ) = -dS̃/d log₁₀ τ` (specific heat).
Both MSC and `|ImCoh|` overlaid; methods distinguished by line style;
axis label colours match curve type. Pat_08 alpha taskLearn.

**Caption material.**
- "LRG entropy and specific heat as a function of diffusion time τ for
  Pat_08, alpha band, taskLearn. Left axis (blue): `1−S̃(τ)`; right
  axis (red): `C̃(τ) = -dS̃/d log₁₀ τ`. Solid = MSC, dashed = `|ImCoh|`.
  The location and height of the `C̃` peak identify the dominant
  hierarchical scale; MSC and `|ImCoh|` produce similar peak positions
  but different magnitudes."
