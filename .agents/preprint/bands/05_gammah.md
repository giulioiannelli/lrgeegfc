---
name: preprint-gammah-band
era: IMCOH_ABS_COHORT_N10
status: current
kind: preprint-result-report
band: high_gamma
range_hz: [80, 300]
range_note: "γ_h covers 80-300 Hz spanning broadband gamma + hippocampal ripple + lower HFO — NOT textbook 60-120 Hz (per memory brain_bands_definition.md)"
cohort: [Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15]
n_patients: 10
fc_method: imcoh_abs
tau: 1/lambda_max
canonical_lrg_object: "D_coph = cophenet(UPGMA(D(tau_max)))"
canonical_subspace_object: "U_k = span{phi_2, ..., phi_{k+1}}"
verdict_tag: "no trace"
verdict_source: VERDICT_LEDGER.md (locked 2026-05-18, REVISED 2026-05-19 — demoted from weak under cluster-extent permutation)
verdict_layers:
  substrate_rank: borderline_positive_direction (raw FC ρ_split^raw 7/10 p=0.188 ratio 88×; cohort-median is positive but does not clear matched-strength)
  rho_split_coph: no_trace_near_null (C3 fails decisively: paired Wilcoxon p=0.246, n_above_surrogate 4/10, obs_median +0.0003 with ratio 0.08× — observation effectively equal to surrogate median)
  grassmann: no_trace_borderline (cluster-extent permutation p=0.0547, audit_70; 9-cell observed run at the null 95th percentile = 8.05 — just outside the 0.05 gate)
  grassmann_epi_excluded: shifted_smaller_window (audit_67; 6-cell longest run at k=21..26 plus emergent cells at k=15..18 — descriptive physiological-attribution hint, does NOT change the verdict)
  cluster_extent_revision: "8-cell hardcoded threshold (previous gate) had γ_h passing → revised gate (audit_70 cluster-extent permutation) demotes to no trace (Decision 6, 2026-05-19)"
  anatomy: not_yet_audited (no γ_h-specific Desikan-Killiany enrichment run; not primary under no-trace verdict)
sources:
  - .agents/preprint/locked/CONTROLS.md (locked 5-control battery)
  - .agents/preprint/locked/VERDICT_LEDGER.md (locked verdicts; Decision 6 for γ_h)
  - .agents/preprint/methods/methods_revision_2026-05-18_cophenet.md (binding methods directive)
  - data/audit/ctm_triangle/cohort_summary.csv (γ_h row: rho_split_median, n_above_drift, rho_xprobe, Wilcoxon)
  - data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv (γ_h row: paired Wilcoxon vs surrogate)
  - data/audit/grassmann_cluster_extent/cohort_summary.csv (γ_h cluster-extent p=0.0547, audit_70)
  - data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv (γ_h per-k T_G + surr stats, audit_66)
  - data/audit/grassmann_epi_exclusion/sensitivity.csv (γ_h per-k under C5 epi-X, audit_67)
  - data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv (γ_h substrate row)
  - data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv (γ_h raw D sensitivity)
revision_history:
  - 2026-05-19: initial brief produced from VERDICT_LEDGER.md lockdown — γ_h Grassmann demoted weak → no trace via cluster-extent permutation (Decision 6)
---

# γ_high band (80–300 Hz) — preprint result report

## Head

γ_h carries **no trace on either probe** under the 2026-05-19 cluster-extent permutation Grassmann gate, demoted from the previous `weak trace` reading that depended on a hardcoded 8-cell contiguous-significant threshold (`locked/VERDICT_LEDGER.md` Decision 6, audit_70). The cophenet `ρ_split^coph` C3 paired Wilcoxon p = 0.246 with `obs_median = +0.0003` (effectively zero, ratio 0.08× to surrogate median), 4/10 patients above their own surrogate — well into the no-trace regime. The Grassmann probe shows a 9-cell observed longest contiguous-significant run at k = 19..27 with cluster_p = 0.0547 — **at the null 95th percentile (8.05)** — just outside the 0.05 gate. Under C5 epi-X (audit_67), the Grassmann window shrinks to a 6-cell run at k = 21..26 with 4 additional cells emerging at k = 15..18 ("physiological-attribution hint"), but the cluster_p does not improve into the strong/weak regime. **Note on band range**: γ_h here is **80–300 Hz** spanning broadband gamma, hippocampal ripple, and the lower HFO range (not the textbook 60–120 Hz γ); this matters when interpreting the band physiologically. Verdict from `locked/VERDICT_LEDGER.md`: **`no trace`** — useful as the demoted band that the cluster-extent revision corrected.

## Headline three-layer cohort table (γ_h-band row, copied from `locked/VERDICT_LEDGER.md`)

```
Layer                                        γ_h (full n=10)
-------------------------------------       ---------------------
Raw FC          ρ_split^raw  (raw_fc_ms)    7/10 p=.188  ratio 88×    ← borderline +
Raw D(τ_max)    ρ_split      (preprint_05)  7/10 p=.161  ratio 115×   ← borderline + (tiny denom)
Cophenet D_coph ρ_split^coph (audit_63)     4/10 p=.246  ratio 0.08×  ← FAILS (obs ≈ surrogate)
Grassmann       d_G(k)       (audit_70)     9-cell run k=19..27, cluster_p = 0.0547  ← no trace (just outside)
```

Sources: `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv`, `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv`, `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv`, `data/audit/grassmann_cluster_extent/cohort_summary.csv`.

**Key reading**: γ_h is the **band where the locked Grassmann gate revision changed the verdict**. Under the previous 8-cell hardcoded threshold (audit_66 era), the 9-cell run at k=19..27 would have passed and γ_h would have been `weak trace, only Grassmann` — the same coverage tag as γ_l. Under the cluster-extent permutation gate (audit_70, 2026-05-19), the 9-cell run is at the null 95th percentile (matched-strength surrogates produce LR=8 at the 95th percentile, so a 9-cell observed run yields cluster_p = 0.0547 — just outside the 0.05 gate). γ_h Grassmann is therefore **no trace** under the locked rule. The cophenet probe fails decisively (obs_median ≈ surrogate_median, paired Wilcoxon p = 0.246 with 4/10 patients above own surrogate). The substrate-level borderline positive cohort medians at raw FC (7/10 p = 0.188) and raw D (7/10 p = 0.161) do not clear the matched-strength gate either. γ_h is `no trace` across the entire 5-control battery.

## 1. Scientific claim

**Cohort-level question.** Does the γ_h-band (80–300 Hz) post-task resting state at the LRG-CTM layer carry task-related reorganization on either of the two locked probes?

**Refined biological claim (null direction).** γ_h coupling does **not** carry a statistically supported task-related reorganization at either the per-pair cophenet multiscale geometry or the leading-mode Laplacian subspace, under the locked 5-control battery. The borderline cohort medians at the substrate and raw-D layers (positive +0.111 and +0.095 respectively at 7/10 patients) do not clear matched-strength, and the LRG-level probes both fail their primary gates. The Grassmann probe sits exactly at the 0.05 boundary (cluster_p = 0.0547) — close enough that the cluster-extent revision was the deciding factor between `weak` and `no trace` (`locked/VERDICT_LEDGER.md` Decision 6).

**Falsification budget.** This `no trace` verdict would flip back to `weak trace` if any of: (a) a re-run with R > 200 surrogates produces an empirical null 95th percentile below 8.0 (currently exactly at 8.05; with a tighter null the 9-cell observation would cross the gate); (b) a different cluster-extent statistic (e.g., cluster mass at a different p-threshold) produces cluster_p < 0.05 — currently the cluster_mass p = 0.065, also outside; (c) the cophenet probe is re-audited and clears C3.

## 2. Cohort, substrate, and library entry points

Same as other bands: n=10 patients (Pat_02/03/05/06/07/08/10/13/14/15); `imcoh_abs` band-averaged magnitude of imaginary coherence; γ_h band `B = [80, 300] Hz`. **Important**: γ_h here is the broadband-gamma + ripple + lower-HFO range, *not* the textbook narrowband γ (60–120 Hz). See `memory/brain_bands_definition.md`. Pat_03 acquired at 1024 Hz — Nyquist 512 Hz, in-band content up to 300 Hz is extracted by the filterbank with anti-aliasing; this is a config-layer fact, not a dropout or outlier flag (policy updated 2026-05-18). Pat_03 is a full cohort member at the analysis layer.

```python
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
W = load_fc_matrix(patient="Pat_02", phase="rest_pre", band="high_gamma", fc_method="imcoh_abs")
lrg = load_lrg_result(patient="Pat_02", phase="rest_pre", band="high_gamma", fc_method="imcoh_abs")
```

## 3. The two LRG probes — methodology, results, and provenance

### 3.1 Per-pair multiscale correlation on `D_coph` — **no trace, near null**

#### Critical preamble (5-point)
1. **Claim.** γ_h `ρ_split^coph > 0` at cohort level.
2. **Null (C3, mandatory).** Strength-preserving 4-cycle ±δ rewiring, R=200, seed 20260511. Paired Wilcoxon one-sided greater across patients: obs > surrogate.
3. **Strongest plausible alternative.** Within-session drift; same-strength edge re-arrangement producing a small positive cohort-median by chance.
4. **Does C3 cover (3).** Yes — matched-strength is the mandatory FC-level null. At γ_h, obs_median and surr_median are essentially equal, so the test directly rejects the trace alternative.
5. **Falsification of the verdict.** Would require C3 paired Wilcoxon p < 0.05. Currently 0.246, far from threshold.

#### Results (γ_h, full cohort)

| Statistic | Value | Source |
|---|---|---|
| Cohort-median `ρ_split^coph` | **−0.0135** | `ctm_triangle/cohort_summary.csv` high_gamma row |
| n_trace_split (`ρ > 0`) | 4/10 | same |
| n_above_drift | 5/10 | same |
| C1 Wilcoxon p (`ρ_split > 0`) | **0.3477** | same col `wilcoxon_split_gt_0_p` — ✗ |
| C2 Wilcoxon p (`ρ_split > ρ_drift`) | **0.4229** | same col `wilcoxon_split_gt_drift_p` — ✗ |
| C3 matched-strength: obs_median | +0.00034 | `matched_strength_surrogate_split_baseline/cohort_summary.csv` high_gamma |
| C3 surr median (median over surrogate medians) | +0.00418 | same |
| C3 ratio obs/surr | **0.08×** (obs effectively at surrogate level) | derived |
| C3 paired Wilcoxon z / p | 35.0 / **0.2461** | same — ✗ |
| C3 n_above_surrogate (per-patient) | **4/10** | same |
| C4 cross-probe rho_xprobe | −0.0161 | `ctm_triangle/cohort_summary.csv` high_gamma |
| C4 n_trace_xprobe | 4/10 (−sign matches rho_split −0.014) | same — sign matches but n=4 < 6 → ✗ |
| C5 epi-X cophenet | **not run** for γ_h (audit_68 covered α only) | — |

**Every primary control fails on the cophenet probe.** The obs_median is essentially identical to the surrogate median (+0.0003 vs +0.0042 — both near zero), the paired Wilcoxon p = 0.246 is well above 0.05, and only 4/10 patients are above their own surrogate. Note that the rho_split cohort median (−0.014) and the C3 obs_median (+0.0003) differ — the former is the median over patients of per-patient ρ_split^coph values, the latter is the median over patients of per-patient median obs_ρ from the surrogate comparison. Both are essentially zero; the small numerical difference reflects definitional plumbing, not a meaningful signal.

#### Reading
γ_h cophenet has **no trace**. The observation is statistically indistinguishable from matched-strength surrogates at every primary control. The substrate-level positive cohort median at γ_h (raw FC +0.111, raw D +0.095) does not survive into the cophenet representation — the multiscale dendrogram aggregation extinguishes whatever weak per-pair structure existed at the substrate level.

#### Cache + script provenance
- Full cohort: `audit_63_split_baseline_surrogate.py` → `data/audit/matched_strength_surrogate_split_baseline/`
- Within-baseline triangle: `audit_*_ctm_triangle.py` → `data/audit/ctm_triangle/`
- Object built from `lrg.ultrametric_matrix` field in `data/cache/imcoh_lrg/Pat_NN/high_gamma_{phase}_lrg_imcoh-abs.npz`

### 3.2 Global subspace rotation: Grassmann chordal distance `d_G(k)` — **no trace (cluster-extent demotion)**

#### Critical preamble (5-point)
1. **Claim.** `T_G(k) < 0` at cohort level for some contiguous window of `k`.
2. **Null (C3, audit_70 cluster-extent permutation).** R=200 matched-strength surrogates; cluster_p = (1 + #(null_LR ≥ obs_LR)) / (R + 1). **Locked gate as of 2026-05-19, replacing the previous 8-cell hardcoded threshold.**
3. **Strongest plausible alternative.** Same-strength edge rewiring producing a similar-length contiguous-significant cluster by chance. At γ_h this is the *exact* failure mode that flipped the verdict: under R=200 surrogates, the empirical null distribution has 95th percentile LR = 8.05, so a 9-cell observed run yields cluster_p ≈ 0.055.
4. **Does C3 cover (3).** Yes by construction — the cluster-extent null IS the matched-strength chance-cluster distribution.
5. **Falsification of the verdict.** Would require cluster_p < 0.05. Currently 0.0547 (longest run) or 0.0647 (cluster mass), both just outside.

#### Results (γ_h)

| Statistic | Value | Source |
|---|---|---|
| Observed longest contiguous-significant run | **9 cells** at k = 19..27 | `grassmann_cluster_extent/cohort_summary.csv` high_gamma + per-k inspection |
| Cluster-extent null mean LR | 2.50 | `grassmann_cluster_extent/cohort_summary.csv` high_gamma |
| Cluster-extent null median LR | 2.0 | same |
| Cluster-extent null **95th LR** | **8.05** | same — directly determines γ_h's borderline status |
| Cluster-extent null max LR | 19 | same |
| Observed cluster mass (Σ −log10p) | 15.75 | same |
| Null mean cluster mass | 4.52 | same |
| Null 95th cluster mass | 16.12 | same |
| **cluster_p_longest_run** | **0.0547** | same — ✗ (just outside 0.05) |
| **cluster_p_cluster_mass** | 0.0647 | same — ✗ |
| C5 epi-X (audit_67): obs longest run | **6 cells at k = 21..26** | `grassmann_epi_exclusion/sensitivity.csv` high_gamma rows |
| C5 epi-X: additional emergent cells | k = 15..18 (4 cells, "emerge" flag) | same — physiological-attribution hint |
| C5 epi-X: cluster_p (re-estimate) | does not cross the 0.05 gate | same |

The observed 9-cell run at k=19..27 sits **exactly at the null 95th percentile (8.05)**. Both cluster-p statistics (longest-run 0.055; cluster-mass 0.065) are just outside the 0.05 gate. By the locked cluster-extent rule, γ_h Grassmann is **no trace**.

This is the band where the 2026-05-19 cluster-extent revision (`locked/VERDICT_LEDGER.md` Decision 6) changed the verdict. Under the previous 8-cell hardcoded threshold, the 9-cell run would have passed and γ_h would have been `weak trace, only Grassmann`. Under cluster-extent, the empirical null shows that 8-cell runs are produced by 5% of matched-strength surrogates, so a 9-cell observed run is at the boundary of chance.

#### C5 epi-X sensitivity (audit_67)
Under epi-zone exclusion:
- Original full-cohort: 9-cell run at k=19..27.
- Epi-X: 6-cell run at k=21..26 (cell retention 67%, within-window overlap k=21..26 = 6/9 ≈ 67%).
- Additional 4 cells **emerge** under epi-X at k=15..18 (significant under epi-X, not significant in full cohort).

The "emergent cells" pattern at k=15..18 is the **physiological-attribution hint**: when epi-zone contacts are removed, a new contiguous-significant region appears at slightly lower k. Mechanistically this is consistent with the epi-zone *masking* a physiological subspace rotation in the lower-k modes — when removed, the underlying physiological signal becomes visible. However:
- The emergent k=15..18 window (4 cells) does NOT form a long-enough contiguous cluster on its own to clear the cluster-extent gate.
- The retained k=21..26 window (6 cells) is shorter than the original 9-cell window.
- A re-estimate of cluster_p under epi-X (treating epi-X as a separate sweep) does not cross the 0.05 gate.

This is descriptive evidence that γ_h *might* carry a physiologically-attributable signal masked by the epi zone — but under the locked rule, this does not change the verdict. **γ_h is `no trace`** under the locked battery; the epi-X observation is a flagged open question, not a positive finding.

#### Reading
γ_h Grassmann fails the cluster-extent permutation gate at cluster_p = 0.055, just outside the 0.05 boundary. The 9-cell observed run at the null 95th percentile is the textbook borderline case where a stronger null discriminates between chance clustering and a real subspace trace. The C5 epi-X "emergent" cells at k=15..18 are interesting physiologically but cannot rescue the verdict — they form a 4-cell cluster that is itself within the null.

#### Cache + script provenance
- `audit_66_grassmann_matched_strength_surrogate.py` → `data/audit/grassmann_matched_strength_surrogate/`
- `audit_67_grassmann_epi_exclusion.py` → `data/audit/grassmann_epi_exclusion/`
- `audit_70_grassmann_cluster_extent.py` → `data/audit/grassmann_cluster_extent/`
- Object built from `eigvecs` field in `data/cache/imcoh_lrg/Pat_NN/high_gamma_{phase}_lrg_imcoh-abs.npz` (obs) and `data/cache/matched_strength_surrogate_lrg/Pat_NN/high_gamma_{phase}_R200_swap20_seed20260511_imcoh_abs.npz` (surrogate)

## 4. Multiscale structure of γ_h (null at LRG layer)

### Per-pair multiscale via `D_coph` (extinguished)

At γ_h the cophenet step **extinguishes** any substrate-level positive direction. Raw FC has 7/10 patients in trace direction with cohort-median +0.111; raw D has 7/10 with cohort-median +0.095; cophenet D_coph has 4/10 with cohort-median −0.014 and a paired-Wilcoxon p = 0.246 against matched-strength. The substrate-to-cophenet trajectory at γ_h is **complete extinction** rather than amplification (β), demotion (γ_l), or sign inversion (θ).

### Subspace multiscale via Grassmann `d_G(k)` (borderline-null)

The 9-cell run at k=19..27 is in the same `k`-range as β's signature (k=27..55) and γ_l's (k=12..23) — these are the mid-`k` modes that carry the largest spatial-scale information. But unlike β and γ_l, γ_h's contiguous-significant window does not extend beyond the empirical null 95th percentile. The probe is "trying to see" something at mid-`k` but the signal is at chance.

The audit_67 epi-X result — a 6-cell core at k=21..26 plus 4 emergent cells at k=15..18 — suggests that under a larger-cohort version of this study, the γ_h signal might cross the cluster-extent gate. With n=10 patients and R=200 surrogates, the boundary at cluster_p = 0.055 is at the resolution limit; a study with n=20 or R=500 might disambiguate this.

## 5. Anatomical distribution — not relevant under no-trace verdict

No γ_h-specific Desikan-Killiany hypergeometric or implant-geometry analysis exists, and none is required under the `no trace` verdict. The audit_67 emergent-cell pattern at k=15..18 under epi-X is the only descriptive hint that γ_h *could* have a physiological signature; localizing those modes anatomically would require a mode-decomposition analysis but is not part of the locked package.

## 6. Patient-by-patient reading

### γ_h `ρ_split^coph` per-patient signs (descriptive — probe is null under C3)
From `ctm_triangle/cohort_summary.csv`: cohort median `ρ_split^coph` = −0.014 with 4/10 in + direction. The C3 paired Wilcoxon p = 0.246 has only 4 patients above their own matched-strength surrogate. No single patient drives the verdict; the cohort is uniformly null at this probe.

### γ_h Grassmann subspace trace per-patient (descriptive)
Per-patient `n_patients_below_own_surrogate` ranges from 4–9 across the k = 19..27 window (from `grassmann_matched_strength_surrogate/cohort_summary.csv` high_gamma rows). The 8/10 and 9/10 cells at k=22..26 are individually significant but do not aggregate to a cluster above the empirical null.

### Pat_03 γ_h Nyquist constraint
Pat_03 is acquired at 1024 Hz (Nyquist 512 Hz vs the 80–300 Hz band); the filterbank extracts in-band content up to 300 Hz with anti-aliasing. This is a config-layer fact about Pat_03's signal coverage at γ_h, not a dropout or outlier framing — Pat_03 remains a full cohort member at the analysis layer.

## 7. Robustness panel (under the locked battery)

| Check | γ_h verdict | Source |
|---|---|---|
| C1 within-rsPre split-half null | ✗ 4/10 cohort, p = 0.348 | `ctm_triangle/cohort_summary.csv` high_gamma |
| C2 drift-floor `ρ_split > ρ_drift` | ✗ 5/10 above drift, p = 0.423 | same |
| C3 matched-strength surrogate `ρ_split^coph` | ✗✗ ratio 0.08× (obs ≈ surrogate), 4/10, p = 0.246 | `matched_strength_surrogate_split_baseline` high_gamma |
| C3 matched-strength surrogate Grassmann (cluster-extent) | ✗ cluster_p = 0.0547, 9-cell run k=19..27 — just outside | `grassmann_cluster_extent` high_gamma |
| C4 cross-probe restriction | ✗ −0.016 sign-matches rho_split but 4/10 < 6/10 | `ctm_triangle` high_gamma (rho_xprobe column) |
| C5 epi-X `ρ_split^coph` | **not run** | (audit_68 was α-only) |
| C5 epi-X Grassmann | 6-cell at k=21..26 + emergent 4 cells at k=15..18; cluster_p not improved | `grassmann_epi_exclusion/sensitivity.csv` high_gamma |
| Three-layer cohort table | substrate borderline + → raw D borderline + → cophenet **extinguishes** to near-zero | headline section above |

## 8. Sensitivity panel under the locked 5-control battery

| Control | γ_h status (D_coph) | γ_h status (Grassmann) | Notes |
|---|---|---|---|
| C1 split | ✗ p = 0.348 | n/a | C1 specific to cophenet |
| C2 drift | ✗ p = 0.423 | n/a | C2 specific to cophenet |
| C3 matched-strength | ✗✗ p = 0.246 (obs ≈ surr, ratio 0.08×) | ✗ cluster_p = 0.0547 (just outside; **demoted from weak under cluster-extent revision**) | mandatory; both probes fail |
| C4 cross-probe | ✗ −sign at 4/10 (< 6/10) | n/a | C4 specific to cophenet |
| C5 epi-X | **not run** (no audit_68 for γ_h cophenet) | partial: 6-cell retention + 4 emergent cells; descriptive only | sensitivity layer |

**Verdict (locked)**: `no trace` per `locked/VERDICT_LEDGER.md` (revised 2026-05-19, Decision 6).

## 9. Interpretation

**Neuroscientific framing.** γ_h here (80–300 Hz) spans **broadband gamma** (80–150 Hz), the **hippocampal ripple range** (~150–250 Hz), and the **lower HFO range** (≥ 80 Hz). This is a clinically and physiologically heterogeneous band: broadband gamma reflects unstructured cortical computation (correlated with multi-unit firing rate; Manning et al. 2009; Ray & Maunsell 2011); ripples are canonical memory-consolidation signatures in hippocampus (Buzsáki 2015); HFOs are clinically associated with epileptogenesis (Jacobs et al. 2012). The negative LRG-CTM verdict at γ_h is therefore **not** a claim about HFOs being absent — it is a claim that whatever γ_h-band coupling reorganization the task induces does not leave a multiscale-LRG signature in either the per-pair cophenet or the leading subspace.

**Why γ_h was demoted.** The 9-cell Grassmann window at k=19..27 was a "promising" signal under the old 8-cell threshold, but the cluster-extent permutation (audit_70) showed that the matched-strength null produces 8-cell runs at the 95th percentile of R=200 surrogates. The 9-cell observation is therefore at the exact boundary of chance clustering. The cluster_p = 0.0547 is a clean expression of "this is one cell above the 95th-percentile chance threshold" — not significant by the 0.05 convention. Decision 6 (`locked/VERDICT_LEDGER.md`) explicitly logs this demotion and is the canonical reference for the verdict change.

**The audit_67 emergent-cell hint.** Under C5 epi-X, the 4-cell window at k=15..18 emerges as significant where it was not in the full cohort. Mechanistically this is consistent with **epi-zone γ_h coupling masking a physiological signal**: the high-amplitude HFO-like coupling near the seizure focus dominates the leading-mode subspace; removing those contacts uncovers a weaker but more physiologically-attributable signal in the slow modes. However, this emergent cluster is itself only 4 cells and does not pass cluster-extent. **The locked verdict is `no trace`** — the audit_67 hint is filed as a descriptive open question for a larger-cohort follow-up, not as a positive finding for this manuscript.

**Why γ_h is informative even as a null.** γ_h's demotion from weak → no trace demonstrates that **the locked cluster-extent gate is more conservative than the previous 8-cell hardcoded threshold** for borderline contiguous-significant runs. Without the cluster-extent revision, γ_h would have appeared as a `weak trace, only Grassmann` band — promoting a chance-clustering pattern to verdict status. The fact that the same revision *promoted* δ from no → weak (cluster_p = 0.0249, well above null mean and 95th) AND *demoted* γ_h from weak → no (cluster_p = 0.0547, exactly at 95th) is evidence that the locked gate is correctly discriminating real subspace structure from chance clustering.

## 10. Source-of-truth references

### Lockdown documents
- `.agents/preprint/locked/CONTROLS.md` — 5-control battery (Grassmann gate locked at cluster-extent permutation)
- `.agents/preprint/locked/VERDICT_LEDGER.md` — locked verdict tag `no trace` for γ_h; Decision 6 logs the cluster-extent revision
- `.agents/preprint/methods/methods_revision_2026-05-18_cophenet.md` — binding methods directive
- `.agents/preprint/bands/01_beta.md`, `.agents/preprint/bands/02_alpha.md`, `.agents/preprint/bands/03_gammalow.md`, `.agents/preprint/bands/04_theta.md` — companion band briefs

### Audit data (rows the numbers come from)
- `data/audit/ctm_triangle/cohort_summary.csv` — γ_h row: C1, C2, C4 statistics
- `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` — γ_h row: C3 cophenet (fails)
- `data/audit/grassmann_cluster_extent/cohort_summary.csv` — γ_h cluster-extent permutation (audit_70)
- `data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv` — γ_h per-k T_G + surrogate stats (audit_66)
- `data/audit/grassmann_epi_exclusion/sensitivity.csv` — γ_h per-k under C5 epi-X (audit_67)
- `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv` — γ_h substrate row
- `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv` — γ_h raw D sensitivity

### Audit scripts (provenance)
- `audit_63_split_baseline_surrogate.py` — C3 ρ_split^coph matched-strength
- `audit_66_grassmann_matched_strength_surrogate.py` — Grassmann C3 per-k
- `audit_67_grassmann_epi_exclusion.py` — Grassmann C5 per-k under epi-X (emergent-cell hint at γ_h)
- `audit_70_grassmann_cluster_extent.py` — Grassmann cluster-extent permutation (2026-05-19; the verdict-flipping audit at γ_h)

### Companion guides
- `.agents/guides/02_methods/lrg-framework-guide.md` — LRG primitive
- `.agents/guides/02_methods/imcoh-guide.md` — imcoh_abs convention
- `.agents/guides/04_rules/never-always-list.md` — coding rules + no-hardcoded-thresholds feedback
- `memory/brain_bands_definition.md` — γ_h is 80–300 Hz (not textbook 60–120)

## 11. What follows

This brief freezes the γ_h `no trace` verdict for the preprint manuscript under the 2026-05-18 / 2026-05-19 lockdown. Next steps:

1. **γ_h figures (optional / contextual)** — if a γ_h panel is wanted to show the cluster-extent demotion explicitly (i.e., why a 9-cell window doesn't pass), it can show: (i) per-k cohort `T_G(k)` curve with the k=19..27 window highlighted, (ii) the audit_70 empirical null distribution of longest-run lengths with the 9-cell observation marked relative to the 95th percentile (8.05) and the cluster_p = 0.0547 annotated. This would be useful in cross-band synthesis as the methodological-control panel. No anatomy figure is required.
2. **No anatomy audit owed.**
3. **No further sensitivity tests owed under the locked battery.** The audit_67 emergent-cell hint at k=15..18 under epi-X is filed as a descriptive open question; resolving it would require either a larger cohort (n ≥ 20) or a tighter null (R ≥ 500 surrogates) to push the 8.05 95th-percentile of the LR null below 9. Out of scope for this manuscript.

Verdict ready for writing-agent handoff: **`no trace`**. The headline story is the **cluster-extent demotion** — γ_h was `weak trace, only Grassmann` under the previous 8-cell hardcoded threshold and is `no trace` under the locked cluster-extent permutation gate (audit_70, 2026-05-19, Decision 6). The audit_67 emergent-cell pattern at k=15..18 under epi-X is a descriptive open question for follow-up, not a positive finding.
