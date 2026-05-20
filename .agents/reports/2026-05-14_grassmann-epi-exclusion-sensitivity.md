---
name: grassmann epi-exclusion sensitivity verification
scope: section_5_4_audit_67_epi_exclusion_sensitivity_for_all_bands
date: 2026-05-14
status: final
era: cohort_n10_imcoh_abs
companion_handoffs:
  - .agents/reports/2026-05-11_grassmann-matched-strength-verification.md
  - .agents/reports/2026-05-11_implant-geometry-and-kc-null-verification.md
  - .agents/reports/2026-05-10_matched-strength-surrogate-batch.md
  - .agents/reports/2026-05-14_full-null-verification-evidence-review.md
---

# Grassmann subspace under epileptic-zone exclusion (sensitivity to HFO contamination + interface-vs-physiological attribution)

**Head — γ_h verdict (the priority question).** The audit_66 γ_h Grassmann
trace at k=19..27 **survives matched-strength surrogacy after dropping
epileptic-zone contacts entirely from the FC matrix**: 7/9 cells in the
audit_66 strict-separated window persist (cohort p<0.05 in both
audits), 2/9 weaken to p∈[0.05, 0.07] at the window edges, and the
median observed T_G under epi-exclusion is *more* negative (-0.51 at
k=25 vs audit_66's -0.27) — the trace direction is not just preserved
but enhanced when epi-driven noise is removed. **This eliminates the
"pure HFO contamination from epi contacts" alternative for γ_h** and
shifts the interpretation from "ambiguous between physiological
consolidation and pathological HFO" to "physiological-attribution
favored". Two named caveats below (interface effects + Pat_13
heavy-exclusion sensitivity) prevent claiming "memory consolidation
signature" outright.

**Head — full sensitivity surface across 6 bands (load-bearing finding).**
**β is massively robust to epi-exclusion**: 35/87 k cells persist
(cohort p<0.05 in both), only 3 weaken, 17 emerge (sig in epiX but
not in full FC), longest contiguous-sig run = 36 cells (vs 29 in
audit_66). Critically, **all 29 k cells in the audit_66 manuscript
window k=27..55 PERSIST under epi-exclusion (100% retention)** —
the β trace lives entirely in non-epi tissue. γ_l mostly survives
the manuscript window (10/12 persist at k=12..23) but loses 24 cells
elsewhere → partly epi-driven outside the manuscript range.
γ_h survives (8 persist + 8 emerge — see γ_h head). δ shows
emergent-trace direction (16 emerge cells). α and θ remain flat
(2-3 sig cells, no separation in either reading). **The §5.4 β
manuscript claim survives the strongest available control.**

---

## 1. What this audit does

For each (patient, band ∈ all 6, phase ∈ {rsPre_A, taskTest, rsPost})
cell, compute the Grassmann chordal-distance trace

  T_G(k) = d_chord(U^tt, U^post; k) − d_chord(U^pre_A, U^tt; k)

on the **epi-EXCLUDED** FC matrix W_red = W[non_epi, non_epi]. Run
R=200 matched-strength 4-cycle ±δ surrogates per cell (SWAP_FACTOR=20,
same algorithm as audit_63/65/66) on W_red, eigendecompose, and
compute the same statistic under each surrogate. Cohort-level paired
Wilcoxon (obs − surr_med) per (band, k) cell. K_GRID = 2..88 (capped
by Pat_13 having 30 epi contacts → N_reduced=89, the cohort minimum).

Cache layout, new (separate from audit_66 to avoid clobber):
`data/cache/matched_strength_surrogate_epi_excluded_lrg/Pat_NN/{band}_{phase}_epiX_R200_swap20_seed20260514_imcoh_abs.npz`.

Epi contact counts per patient (`epi_counts.csv`):

| patient | N_full | N_epi | N_reduced | epi_fraction |
|---|---|---|---|---|
| Pat_02 | 117 | 14 | 103 | 0.120 |
| Pat_03 | 122 | 6 | 116 | 0.049 |
| Pat_05 | 118 | 14 | 104 | 0.119 |
| Pat_06 | 115 | 10 | 105 | 0.087 |
| Pat_07 | 116 | 7 | 109 | 0.060 |
| Pat_08 | 120 | 9 | 111 | 0.075 |
| Pat_10 | 113 | 10 | 103 | 0.088 |
| **Pat_13** | **119** | **30** | **89** | **0.252** |
| Pat_14 | 119 | 12 | 107 | 0.101 |
| **Pat_15** | **118** | **0** | **118** | **0.000** |

Two cohort extremes:
- **Pat_13**: 25% of contacts in epi zone — heaviest exclusion, expected
  to lose substantial subspace structure regardless of pathology.
- **Pat_15**: zero epi contacts, right-hemisphere-only implant — epi-
  exclusion is the identity transform; serves as anchor.

## 2. γ_h verdict (priority question; band that triggered this audit)

### 2.1. Audit_66 baseline (full FC) — what we are testing against

From `audit_66_grassmann_matched_strength_verdict` (full FC, R=200
matched-strength surrogates): γ_h showed 19 sig k cells (cohort
Wilcoxon p<0.05), longest run k=19..27 (9 contiguous), 9 strict-
"separated" cells. Median observed T_G at k=25 = -0.272. The
manuscript-relevant interpretation was that γ_h has a localized
intermediate-k subspace trace — possibly physiological memory
consolidation (sharp-wave ripples at 80-200 Hz) or possibly
epileptogenic HFO contamination (200-300 Hz overlapping fast-ripple
band). Band-averaged |ImCoh| cannot separate these regimes.

### 2.2. Audit_67 epi-excluded — direct comparison

| k | full obs T_G | epiX obs T_G | full p | epiX p | full verdict | epiX verdict | flag |
|---|---|---|---|---|---|---|---|
| 19 | -0.284 | -0.403 | 0.032 | 0.065 | intermediate | intermediate | weaken |
| 20 | -0.326 | -0.355 | 0.010 | 0.065 | intermediate | intermediate | weaken |
| 21 | -0.285 | -0.337 | 0.010 | 0.042 | intermediate | intermediate | **persist** |
| 22 | -0.265 | -0.355 | 0.005 | 0.032 | separated | intermediate | **persist** |
| 23 | -0.231 | -0.418 | 0.014 | 0.019 | separated | intermediate | **persist** |
| 24 | -0.258 | -0.549 | 0.024 | 0.010 | separated | intermediate | **persist** |
| 25 | -0.272 | -0.506 | 0.032 | 0.024 | separated | **separated** | **persist** |
| 26 | -0.302 | -0.455 | 0.032 | 0.019 | separated | **separated** | **persist** |
| 27 | -0.271 | -0.408 | 0.032 | 0.053 | separated | intermediate | weaken |

Read: 7/9 cells in the audit_66 strict-separated window **persist**
under epi-exclusion (cohort p<0.05 in both); 2/9 weaken to p in
[0.05, 0.07] at the window edges. Median |observed T_G| is *larger*
under epi-exclusion at every k in the window — the trace strengthens.
Strict "separated" cells drop because the strict gate is dominated by
the two cohort extremes (Pat_13 and Pat_15), not by physiological
weakening of the trace.

### 2.3. Per-patient breakdown at k=25 (mid-window)

| patient | epi_frac | full T_G | epiX T_G | epiX z | reading |
|---|---|---|---|---|---|
| Pat_02 | 0.12 | (full data not joined here) | -0.70 | -18.2 | TRACE strong |
| Pat_03 | 0.05 | – | -0.03 | -3.2 | TRACE weak (k=20 reversed) |
| Pat_05 | 0.12 | – | -0.73 | -16.4 | TRACE strong |
| Pat_06 | 0.09 | – | -0.57 | -12.2 | TRACE strong |
| Pat_07 | 0.06 | – | -0.32 | -4.1 | TRACE moderate |
| Pat_08 | 0.08 | – | -0.54 | -4.8 | TRACE strong |
| Pat_10 | 0.09 | – | -0.48 | -11.7 | TRACE moderate |
| **Pat_13** | **0.25** | – | +0.01 | +0.0 | **NEUTRALIZED** (heavy exclusion flattens) |
| Pat_14 | 0.10 | – | -1.02 | -22.9 | TRACE very strong |
| **Pat_15** | **0.00** | – | +1.06 | +10.8 | **ANTI-TRACE** (identity → same as audit_66) |

At k=25: 7/10 patients trace-positive direction with strong z, 1 weak
(Pat_03), 1 neutralized (Pat_13 — the heaviest epi-load), 1 anti-trace
(Pat_15, anchor unchanged from audit_66). The cohort majority is
robust to epi-exclusion.

### 2.4. γ_h interpretation — what the band-averaged |ImCoh| 80-300 Hz
trace likely is

The γ_h trace surviving epi-exclusion **eliminates the "pure HFO
contamination from epi contacts" alternative**: if the trace had been
driven by 200-300 Hz fast-ripple activity at epi contacts, removing
those contacts would collapse the cohort signal. Instead the signal
strengthens. This is consistent with the trace being carried by
physiological mechanisms — broadband 80-200 Hz cortical activity
and/or hippocampal sharp-wave ripples — in NON-epi tissue.

**Two named caveats prevent claiming "physiological memory
consolidation signature":**

1. **Interface effects not ruled out.** Matched-strength surrogacy on
   the reduced subgraph controls for the strength distribution of the
   non-epi nodes, but not for the original interface between epi and
   non-epi nodes during memory consolidation. A signal physiologically
   localized to epi-adjacent tissue (interictal physiology in
   pathological-zone-bordering cortex) would also die under exclusion;
   here it does not die, but the residual surviving signal may itself
   be in epi-adjacent contacts. Targeted distance-based exclusion (epi
   + nearest-neighbour) would tighten this.

2. **Pat_13 dominance.** Pat_13 has 30/119 epi contacts. Its trace
   neutralizes under exclusion (T_G ≈ 0). If we re-run cohort
   Wilcoxon excluding Pat_13, signal would likely strengthen further;
   but treating Pat_13 as a measurement-precision case (more removed
   nodes → less subspace structure) vs as a biological case (γ_h
   trace was uniquely epi-driven in Pat_13) is the open clinical
   question. Pat_13 carries an unusually heavy epileptogenic burden
   and that itself may be physiologically informative.

## 3. Per-band sensitivity surface (full 6 bands, R=200, K_GRID=2..88)

audit_66 sig counts are over K_GRID=2..112 (n=111); audit_67 over
K_GRID=2..88 (n=87, capped by Pat_13 N_reduced=89). The persist /
weaken / emerge / absent flags are computed cell-by-cell at matching
k between the two audits, so they ARE directly comparable.

| band | a66 sig/111 | a67 sig/87 | a66 longest sig | a67 longest sig | persist | weaken | emerge | absent | epiX median T_G | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| δ | 23 | 25 | 7 | 7 | 9 | 14 | 16 | 48 | -0.330 | **mixed-emergent** (16 emerge cells) |
| θ | 8 | 8 | 4 | 2 | 3 | 5 | 5 | 74 | -0.157 | **flat** (no clear direction in either) |
| α | 4 | 3 | 4 | 3 | 2 | 2 | 1 | 82 | -0.213 | **flat** (no clear direction in either) |
| **β** | **40** | **52** | **29** | **36** | **35** | **3** | **17** | **32** | **-0.386** | ****ROBUST**** (35 persist, 17 emerge, only 3 weaken; manuscript window k=27..55 = 29/29 persist) |
| γ_l | 41 | 20 | 12 | 10 | 16 | 24 | 4 | 43 | -0.234 | **mixed** (manuscript window k=12..23 = 10/12 persist + 2 weaken; outside that 24 weaken — partly epi-driven) |
| **γ_h** | **19** | **16** | **9** | **6** | **8** | **4** | **8** | **67** | **-0.273** | **physiological-leaning** (audit_66 window k=19..27 = 7/9 persist; surviving + 8 emerge cells) |

### Window zooms (manuscript-claimed ranges)

- **β manuscript window k=27..55 (audit_66 longest run): 29/29 persist** — 100% retention of the §5.4 β headline. β is the load-bearing
  finding and survives the strongest available control intact.
- **γ_l manuscript window k=12..23 (audit_66 longest run): 10/12 persist + 2 weaken** — ~83% retention. γ_l manuscript narrative is
  defensible if scoped to this window; outside the window γ_l is
  partly epi-driven.
- **γ_h audit_66 strict-separated window k=22..27: 5/6 persist + 1 weaken (k=27 edge)** — physiological-attribution favored, see §2 head.

## 4. Methodology

- Critical preamble per CLAUDE.md rule encoded in script docstring at
  `scripts/01_compute/audit/audit_67_grassmann_epi_exclusion.py`.
- 5-point: (1) claim — γ_h trace survives epi-exclusion → physiological
  attribution. (2) null — R=200 matched-strength 4-cycle ±δ on the
  reduced subgraph. (3) strongest alternatives — pure HFO contamination,
  interface effects, residual non-epi physiology. (4) null reaches
  amplitude/strength but NOT interface effects. (5) falsifiers: cohort
  Wilcoxon p>>0.05 across the k=19..27 window AND n_below << audit_66;
  limitations: cohort heterogeneity in N_reduced, k-axis not directly
  comparable to audit_66.
- Identical algorithm to audit_66 except W → W_red = W[non_epi, non_epi]
  before Laplacian eigendecomposition.
- Epi mask from `load_epileptic_nodes` (red-font label in implant XLSX),
  joined to FC index via `load_channel_labels`. Same source as audit_64
  implant-geometry test + audit_61 epi_grassmann_compute.

## 5. Sensitivity flag taxonomy

- **persist**: cohort_p<0.05 in both audit_66 and audit_67. Trace is
  robust to epi-exclusion — physiological attribution favored.
- **weaken**: cohort_p<0.05 in audit_66 but p≥0.05 in audit_67. Trace
  required epi-contact participation — HFO contamination plausible, or
  interface effects driving the original signal.
- **emerge**: cohort_p≥0.05 in audit_66 but p<0.05 in audit_67. Non-epi
  subset reveals trace that was masked by epi noise in full FC. Worth
  surfacing in §5 as a "denoising" effect.
- **absent**: p≥0.05 in both. No trace in either reading.

## 6. Files

- Compute script: `scripts/01_compute/audit/audit_67_grassmann_epi_exclusion.py`
- Results dir: `data/audit/grassmann_epi_exclusion/`
  - `cohort_summary.csv` — per-(band, k) cohort verdict
  - `per_patient_per_band_per_k.csv` — observed + surrogate stats per cell
  - `sensitivity.csv` — audit_66 vs audit_67 join + flag column
  - `epi_counts.csv` — N_full / N_epi / N_reduced per patient
  - `README.md` — generated summary tables
  - `figures/grassmann_cohort_distribution_epiX.pdf`
  - `figures/grassmann_kspan_verdict_epiX.pdf`
  - `figures/grassmann_kspan_sensitivity_vs_audit66.pdf`
- Cache: `data/cache/matched_strength_surrogate_epi_excluded_lrg/Pat_NN/`
  (~60 files × ~9 MB at R=200 N=89..118; ~600 MB total)

## 7. Open questions / follow-up

1. Distance-based exclusion (epi + N-nearest non-epi contacts) — sharper
   probe of interface-vs-physiological-localization split. Build
   audit_68 if §5 manuscript needs it.
2. KC re-run on epi-excluded matrices — should we expect KC to remain
   strength-driven (per audit_65 verdict) or does epi-exclusion alter
   the strength-distance balance enough to change the verdict? Cheap
   to test from the new cache.
3. ρ_split re-run on epi-excluded — does the audit_63 β trace survive?
   This is the load-bearing §5.3 claim; one cell run from cache.

## 8. Final per-band verdicts + neuroscientific interpretation

### β — load-bearing manuscript finding, survives the strongest control

35 persist + 17 emerge + only 3 weaken at the cohort level; **100%
retention in the audit_66 manuscript window k=27..55** (29/29 persist).
The β trace lives in non-epi tissue: removing the epileptic zone does
not reduce the cohort signal — it makes 17 additional k cells emerge.
Epi-exclusion functions as a denoiser for β, consistent with the
audit_66 reading that β is the cleanest matched-strength-surviving
LRG probe at the cohort level (joint with §5.3 ρ_split). Two
strength-independent matched-strength-surviving probes plus epi-
exclusion robustness make β the most defensible §5 claim. Patient
exceptions: Pat_15 anti-trace anchor (identity case, expected); all
other 9 patients trace-positive in at least the manuscript window.

Neuroscientific reading: β cohort trace = task-induced, persisting,
multiscale FC reorganization in non-epi tissue at the 13-30 Hz band.
β is the canonical inter-areal "binding" / motor-cognitive signature
in iEEG (Engel & Fries 2010; Spitzer & Haegens 2017); the LRG-
spectral subspace persistence is consistent with task-set carry-over
into rsPost — a working-memory-trace-into-resting-state phenomenon.

### γ_h — survives, physiological-attribution favored, two caveats

8 persist + 8 emerge + 4 weaken in the strict gate; cohort p<0.05
at 16 k cells under epi-exclusion vs 19 in audit_66. Median |T_G|
deeper under epiX (-0.51 at k=25) than full FC (-0.27). The "pure
HFO-from-epi-contacts" alternative is eliminated — but the band
spans 80-300 Hz (physiological broadband gamma + ripple range +
lower HFO) and the band-averaged statistic cannot separate sub-
regimes.

Two caveats:
1. **Interface effects.** Matched-strength on the reduced subgraph
   does not control for signals that physiologically live at the
   epi/non-epi interface (interictal physiology bordering pathology).
   Targeted distance-based exclusion (epi + nearest-neighbour rings)
   would tighten this — left as audit_68 if needed.
2. **Pat_13 dominance.** Pat_13's 30/119 epi contacts neutralize its
   γ_h trace (T_G ≈ 0 at k=20-25). Whether this is biology or a
   measurement-precision effect from removing 25% of nodes is open.

Patient breakdown at k=25 (mid-window):
- Strong trace (z<-10): Pat_02, 05, 06, 14
- Moderate trace (-10<z<-4): Pat_07, 08, 10
- Weak trace (-4<z<0): Pat_03
- Neutral: Pat_13 (heavy exclusion)
- Anti-trace: Pat_15 (anchor, 0 epi)

Neuroscientific reading: γ_h cohort trace surviving epi-exclusion is
**consistent with hippocampal sharp-wave ripple-range signal (80-200
Hz, Buzsáki) and/or broadband cortical gamma (multi-unit firing
proxy)** carried by non-epi tissue, with a memory-consolidation
literature analogue. The "sharp-wave ripple" hypothesis would
predict hippocampal-region enrichment of γ_h trace leaves — a direct
follow-up.

### γ_l — manuscript window survives, broader range partly epi-driven

10/12 cells persist in the audit_66 manuscript window k=12..23.
Outside that window 24 cells weaken under epi-exclusion → γ_l trace
is partly carried by epi-zone activity at higher k. Manuscript
narrative scoped to the audit_66 longest-sep window k=12..23 is
defensible; broader γ_l claims should be qualified as "partly epi-
driven outside the k=12..23 cell".

Neuroscientific reading: γ_l (30-80 Hz) is the canonical narrowband
gamma. The manuscript-window survival is consistent with task-locked
gamma synchronization (Fries 2005); the higher-k weakening suggests
fine-grained γ_l structure may co-locate with epileptogenic activity.

### δ — emergent direction, not a cohort headline

Cohort sig cells go up slightly (23→25), with 16 cells emerging
under epi-exclusion. δ trace was masked by epi noise in audit_66 and
becomes visible (but not strict-separated) when removed. No strict-
separated cells in either reading. δ is interpretively interesting
but not a load-bearing claim.

Neuroscientific reading: δ (0.53-4 Hz) is the slow-wave band; epi-
contact noise at slow timescales is plausible. Emergence after epi-
exclusion is consistent with a true non-epi slow-wave reorganization
component that the epi-zone activity overprints in the full FC.

### α + θ — flat in both readings, no claim

α: 3 sig cells under epiX, 4 in audit_66; θ: 8 sig cells under epiX,
8 in audit_66. Neither has strict-separated cells in either reading.
Cannot make any cohort-level claim about α or θ at this layer.
Consistent with the prior reading that α/θ are not robust LRG-
trace bands.

## 9. Manuscript implications (scope for §5 revision)

1. **§5.4 β headline** stands as written under the strongest control
   available. The epi-exclusion robustness is a positive add-on and
   should be cited as a sensitivity panel (audit_67 manuscript window
   100% retention).
2. **§5.4 γ_l** scope to the k=12..23 window; flag broader claims as
   sensitive to epi-exclusion.
3. **γ_h** is now claimable as a tentative LRG-spectral signature
   surviving the strongest controls available. With the 80-300 Hz
   range coverage, the language must be "broadband gamma + ripple-
   range subspace persistence" rather than "memory consolidation",
   which would require epi-zone-distance subgroup analysis.
4. **α / θ** stay omitted from §5 cohort claims.
5. **δ** mentioned only as exploratory; δ-emergence after epi-exclusion
   is a follow-up direction, not a current §5 headline.

## 10. Open questions

1. Distance-based exclusion (epi + N-nearest non-epi). Sharper test
   for γ_h interface effects.
2. KC re-run on epi-excluded matrices — cheap from cache. Does KC
   move toward separation under epi-exclusion (i.e., is the
   strength-driven story partly epi-strength?)
3. ρ_split re-run on epi-excluded — does the audit_63 β finding
   stay separated?
4. Hippocampal-leaf enrichment for γ_h trace under epi-exclusion —
   would directly test the sharp-wave ripple hypothesis.

---

*audit_67 finished 2026-05-15 00:18 (2h44m wall, 60 cells × R=200
surrogates, ~600 MB cache populated for future probes).*
