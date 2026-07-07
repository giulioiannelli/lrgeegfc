---
name: preprint-delta-band
era: IMCOH_ABS_COHORT_N10
status: current
kind: preprint-result-report
band: delta
range_hz: [0.53, 4]
cohort: [Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15]
n_patients: 10
fc_method: imcoh_abs
tau: 1/lambda_max
canonical_lrg_object: "D_coph = cophenet(UPGMA(D(tau_max)))"
canonical_subspace_object: "U_k = span{phi_2, ..., phi_{k+1}}"
verdict_tag: "weak trace, only Grassmann"
verdict_source: VERDICT_LEDGER.md (locked 2026-05-18, REVISED 2026-05-19 — promoted from no trace under cluster-extent permutation)
verdict_layers:
  substrate_rank: trace_direction_passing (raw FC ρ_split^raw matched-strength p=0.042 at 7/10 cohort, ratio 24.6× — the substrate at δ passes matched-strength, one of only two bands that do so along with α)
  rho_split_coph: no_trace_obs_equals_surrogate (C3 fails: paired Wilcoxon p=0.278, n_above_surrogate 4/10, obs_median +0.0076 vs surr_median +0.0078 — ratio 0.98×, observation effectively indistinguishable from surrogate)
  grassmann: weak_trace (cluster-extent permutation cluster_p_mass=0.005 Decision-8 gate at floor, cluster_p_LR=0.025 descriptive co-statistic, audit_70 all-clusters; obs cluster mass 38.07 raw / 0.149 normalized; 7-cell observed run at k=57..63 — above null 95th percentile of 18.37; **Decision-12 LOO max p_mass=0.055 Pat_08 FAILS < 0.05 precondition — this is why the verdict stays "weak" rather than promoting to "strong"**)
  grassmann_epi_excluded: secondary mechanistic observation strengthens decisively (audit_72; raw mass strengthens 38.07 → 43.99, p_mass^epi-X = 0.005 at floor; LOO under epi-X max = 0.005 Pat_02 fully robust — the full-data Pat_08 leverage is attributable to epi-zone interactions rather than the true trace; reported as secondary, not as verdict-promoter per Decision 10)
  rho_split_coph_wm_excluded: SECONDARY_absent_under_wm_exclusion (C6, audit_83 2026-06-08 / amended audit_85 2026-06-12 → data/audit/wm_stratified/; cophenet null with or without white matter — full p=0.216, exclude_wm p=0.500, obs ≈ 0; δ has no per-pair cophenet trace, unchanged. The δ RAW substrate "emerges" under exclude_wm (matched-strength p=0.019) but this is NOT WM-specific — the random-node-decimation control (audit_85) gives raw-δ cohort p_dec=0.355, generic_nodecount; do not cite raw-δ WM-emergence as a tissue effect)
  grassmann_wm_excluded: SECONDARY_WEAKENS_on_C3_gate_under_wm_exclusion (C6, audit_86 2026-06-12; on the LOCKED C3 cluster-extent mass gate δ Grassmann FAILS under WM exclusion — cluster_p_mass^wmX=0.105 (no_trace), down from the full-graph C3 mass at floor; confirming the per-k weakening (23→7) ON THE ACTUAL GATE, not just the lighter per-k count. δ is the ONE trace component that depends on white matter, and it is precisely the cross-probe epileptic-biology δ channel (Decision 5 anchor anatomy), mechanistically separate from the α/β task trace. Dissociation: δ Grassmann STRENGTHENS under epi-exclusion (audit_72) yet WEAKENS under WM-exclusion → the δ subspace is WM-supported and epi-independent. SECONDARY/mechanistic per CONTROLS §C6; δ verdict tag UNCHANGED)
  anchor_anatomy_known_biology: descriptive_only (C4 passes formally at +0.032/6/10 +sign — but this is the published δ "anchor" cross-probe ratio 1.55× known epileptogenesis pattern, NOT a positive trace claim — LEDGER Decision 5)
  anchor_anatomy_descriptive: documented (δ anchor anatomy is the n=10 verdict in memory/epileptic_imcoh_universal.md — cross-probe ratio 1.55× is known biology, not a discovery)
  anatomy_grassmann_full: RETRACTED_0_of_4_cohort_supported (2026-05-30 signed sampling-aware audit, locked in ANATOMY_LEDGER: the 4-DK-region "temporal + parietal + frontal" list does NOT survive — only ctx-lh-inferiorparietal passes and it is n=1, while the well-sampled regions [superiortemporal n=4, inferiortemporal n=5] are ANTI-localized [−0.10, −0.03]; per-patient null. The Grassmann subspace probe does not localize in any band [2026-06-10 atlas]. β is the only band that localizes [→OFC, cophenet]. Historical list retained in §5/body as record. Source data/audit/anatomy_localization_wilcoxon/README.md)
  anatomy_grassmann_epiX: RETRACTED_0_of_3_and_noop_mask_bug (2026-05-30, locked in ANATOMY_LEDGER: the "fully disjoint network" framing is MOOT — 0/3 epi-X regions are cohort-localized AND the δ epi-X mask was a no-op `np.isin(int,str)` bug that never excluded epileptic nodes [fixed 2026-05-30], so the "two disjoint networks" contrast was partly an artifact. ctx-rh-rostralmiddlefrontal is anti-localized. Neither network localizes; the dissociation claim is withdrawn. Historical detail retained as record.)
sources:
  - .agents/preprint/locked/CONTROLS.md (locked 5-control battery)
  - .agents/preprint/locked/VERDICT_LEDGER.md (locked verdicts; Decision 5 + Decision 6 for δ)
  - .agents/preprint/methods/methods_revision_2026-05-18_cophenet.md (binding methods directive)
  - data/audit/ctm_triangle/cohort_summary.csv (δ row: rho_split_median, n_above_drift, rho_xprobe, Wilcoxon)
  - data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv (δ row: paired Wilcoxon vs surrogate)
  - data/audit/grassmann_cluster_extent/cohort_summary.csv (δ cluster-extent p=0.0249, audit_70)
  - data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv (δ per-k T_G + surr stats, audit_66)
  - data/audit/grassmann_epi_exclusion/sensitivity.csv (δ per-k under C5 epi-X, audit_67)
  - data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv (δ substrate row — passes matched-strength)
  - data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv (δ raw D sensitivity)
  - memory/epileptic_imcoh_universal.md (δ anchor anatomy n=10 verdict)
revision_history:
  - 2026-05-19: initial brief produced from VERDICT_LEDGER.md lockdown — δ Grassmann promoted no trace → weak trace via cluster-extent permutation (Decision 6); cophenet remains no trace; C4 descriptive anchor-anatomy noted (Decision 5)
  - 2026-05-26: Decision-12 cascade — δ Grassmann verdict STAYS "weak" but the gate citation switches from `cluster_p_LR = 0.0249` to `cluster_p_mass = 0.005` (Decision-8 mass-only) with LOO max `p_mass = 0.055 Pat_08` as the actual demoter (fails Decision-12 < 0.05 LOO precondition). Pre-fix raw values 12.78 / 4.01 / 9.81 updated to post-fix all-clusters values 38.07 / 8.06 / 18.37. C5 epi-X reframed as secondary mechanistic observation per Decision 10 (strengthens rather than shifts; LOO under epi-X = 0.005 Pat_02 fully robust; the previous "0% overlap window shift" framing is descriptive only, not verdict-bearing).
  - 2026-06-08: added **C6 white-matter-exclusion** sensitivity (frontmatter `verdict_layers.rho_split_coph_wm_excluded` + `grassmann_wm_excluded`) from audit_83/84 (`data/audit/wm_stratified/`). δ is the ONE band where WM removal WEAKENS a component: Grassmann per-k significant cells 23→7 (cophenet stays absent, full p=0.216 → exclude_wm p=0.500; raw-δ exclude_wm emerges p=0.019). The WM-dependent δ subspace is precisely the cross-probe epileptic-biology δ channel (Decision 5), mechanistically separate from the α/β task trace — δ Grassmann STRENGTHENS under epi-X (audit_72) yet WEAKENS under WM-X. SECONDARY/mechanistic per CONTROLS §C6, NOT a gate (per-k Wilcoxon basis, not C3 cluster-extent); **δ verdict tag UNCHANGED** ("weak trace, only Grassmann"). See VERDICT_LEDGER.md 2026-06-08 revision + `responses/2026-06-08_wm-exclusion-cascade.md`.
  - 2026-06-12: **resolved the two C6 honesty flags** (audit_85 + audit_86). Flag B (Grassmann gate): the δ Grassmann weakening is now CONFIRMED on the locked C3 cluster-extent mass gate — `cluster_p_mass^wmX=0.105` (no_trace) under WM exclusion, confirming the per-k 23→7 weakening on the actual gate (δ is WM-supported, epi-independent; the cross-probe epi-biology channel). Flag A (node-count): the δ RAW "emergence" under exclude_wm is NOT WM-specific either (decimation raw-δ p_dec=0.355, generic_nodecount) — do not cite it as a tissue effect. Updated both `verdict_layers.*_wm_excluded` entries. **δ verdict tag UNCHANGED** ("weak trace, only Grassmann"). See VERDICT_LEDGER.md 2026-06-12 amendment.
---

# δ band (0.53–4 Hz) — preprint result report

## Head

δ carries a **whole-network subspace trace** detected by the Grassmann probe under cluster-extent permutation (audit_70 all-clusters `cluster_p_mass = 0.005` Decision-8 gate at empirical floor; `cluster_p_LR = 0.025` descriptive co-statistic) at a 7-cell contiguous-significant window k=57..63. **The verdict is "weak" rather than "strong" because full-data LOO max `p_mass = 0.055 (Pat_08)` fails the Decision-12 < 0.05 LOO robustness precondition** — a single patient leverages the cohort verdict over the gate. **No per-pair cophenet trace**: `ρ_split^coph` C3 paired Wilcoxon p = 0.278, with `obs_median = +0.0076` essentially equal to `surr_median = +0.0078` (ratio 0.98×) and only 4/10 patients above their own surrogate. The substrate-level raw FC at δ **passes** matched-strength (cohort-median +0.111, ratio 24.6×, p=0.042 at 7/10) — one of only two bands where the substrate alone clears the matched-strength gate (α is the other). Under C5 epi-X (audit_72, secondary mechanistic observation per Decision 10): the trace **strengthens** (raw mass 38.07 → 43.99, `p_mass^epi-X = 0.005`, LOO under epi-X max = 0.005 (Pat_02) fully robust); the full-data Pat_08 leverage is attributable to epi-zone interactions rather than the true biological trace. C5 epi-X is reported as interesting mechanistic biology, not as a verdict-promoter — the full-data verdict remains "weak". The δ C4 cross-probe trace (+0.032 at 6/10 +sign) **formally passes the C4 gate** but is the published δ "anchor anatomy" 1.55× cross-probe ratio (known epileptogenesis pattern per `memory/epileptic_imcoh_universal.md`), **not** a trace claim — LEDGER Decision 5. Verdict from `locked/VERDICT_LEDGER.md` (Decision 12 cascade 2026-05-26): **`weak trace, only Grassmann`** — cohort gate passes at floor but full-data LOO Pat_08 fragility prevents strong-tier promotion.

## Headline three-layer cohort table (δ-band row, copied from `locked/VERDICT_LEDGER.md`)

```
Layer                                       δ (full n=10)
-------------------------------------      ---------------------
Raw FC          ρ_split^raw  (raw_fc_ms)   7/10 p=.042  ratio 24.6×  ← PASSES substrate
Raw D(τ_max)    ρ_split      (preprint_05) 7/10 p=.116  ratio 11.5×  ← borderline (just outside)
Cophenet D_coph ρ_split^coph (audit_63)    4/10 p=.278  ratio 0.98×  ← FAILS (obs ≈ surrogate)
Grassmann       d_G(k)       (audit_70)    7-cell run k=57..63, cluster_p_mass = 0.005  ← weak trace (LOO Pat_08=0.055 fails Decision-12 precondition)
```

Sources: `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv`, `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv`, `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv`, `data/audit/grassmann_cluster_extent/cohort_summary.csv`.

**Key reading**: δ is the **band where the substrate signal exists but does not propagate into per-pair LRG**, and the Grassmann probe rescues a subspace signature at high `k` (60-ish modes). The substrate-to-cophenet trajectory at δ shows substrate passing → raw D borderline → cophenet fully extinguished (obs and surrogate essentially identical). The LRG cophenet step does *not* preserve the substrate signal at δ. The Grassmann probe rescues a subspace signature at k=57..63 (`cluster_p_mass = 0.005` at floor; `T_G^* = 38.07` raw / 0.149 normalized; null mean 8.06, null 95th 18.37). The previous 8-cell hardcoded threshold had δ as `no trace`; the cluster-extent permutation (audit_70, 2026-05-19) promoted it under the Decision-8 mass gate. **Under Decision 12 (2026-05-26) the verdict stays "weak" because full-data LOO max p_mass = 0.055 (Pat_08) fails the < 0.05 robustness precondition — a single patient leverages the cohort verdict.** C5 epi-X (audit_72) strengthens the trace (mass 38.07 → 43.99, LOO under epi-X = 0.005 fully robust), interpreted as the full-data Pat_08 leverage being attributable to epi-zone interactions rather than the true biological signal — reported as secondary mechanistic observation per Decision 10, not as a verdict-promoter.

## 1. Scientific claim

**Cohort-level question.** Does the δ-band (0.53–4 Hz) post-task resting state at the LRG-CTM layer carry task-related reorganization on either of the two locked probes? Separately: does the δ-band cross-probe anchor anatomy reproduce the published n=10 verdict (1.55× cross-probe ratio, known epileptogenesis biology)?

**Refined biological claim.** δ coupling reorganizes during the task at a level visible in the raw FC substrate (cohort-median +0.111 at 7/10 patients, matched-strength p=0.042), but the per-pair LRG cophenet representation **does not retain** this signal (obs_median +0.0076 essentially equal to surrogate +0.0078). The Grassmann probe detects a subspace trace at high `k` (k=57..63, 7-cell window, `cluster_p_mass = 0.005` Decision-8 gate at floor; `T_G^* = 38.07` raw / 0.149 normalized). The verdict is **weak** because full-data LOO max `p_mass = 0.055 (Pat_08)` fails the Decision-12 < 0.05 LOO precondition — Pat_08 single-handedly drags the cohort verdict over the gate at full data. Under C5 epi-X (audit_72) the trace strengthens decisively (mass 38.07 → 43.99, p_mass^epi-X = 0.005, LOO under epi-X = 0.005 Pat_02 fully robust), interpreted as the full-data Pat_08 leverage being attributable to epi-zone interactions rather than the true biological trace — secondary mechanistic observation per Decision 10, not a verdict-promoter. The separate δ cross-probe anchor anatomy (1.55× cross-probe ratio, +sign at 6/10) is a **known epileptogenesis biology** (per memory `epileptic_imcoh_universal.md`) and is not a positive trace claim — its inclusion in the per-band brief is *descriptive*, confirming the pipeline reads the expected anchor structure (LEDGER Decision 5).

**Falsification budget.** Under Decision 12 the verdict would flip to `no trace` if `cluster_p_mass ≥ 0.05` (currently 0.005 at floor) — far from threshold. Promotion to `strong` would require full-data LOO max `p_mass < 0.05` (currently 0.055, just above gate). The cophenet verdict (`no trace`) would flip to `weak trace` only if a re-audit shows obs_median substantially above surrogate (currently they are within 3% of each other).

## 2. Cohort, substrate, and library entry points

Same as other bands: n=10 patients (Pat_02/03/05/06/07/08/10/13/14/15); `imcoh_abs` band-averaged magnitude of imaginary coherence; δ band `B = [0.53, 4] Hz`. The lower bound 0.53 Hz is the high-pass cutoff of the iEEG amplifier. Pat_03 1024 Hz outlier handled by rank-/ratio-based statistics. Pat_10 task rows [53, 54, 55] dropped at load. Pat_14 vendor-replaced 2026-04-25.

```python
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
W = load_fc_matrix(patient="Pat_02", phase="rest_pre", band="delta", fc_method="imcoh_abs")
lrg = load_lrg_result(patient="Pat_02", phase="rest_pre", band="delta", fc_method="imcoh_abs")
```

## 3. The two LRG probes — methodology, results, and provenance

### 3.1 Per-pair multiscale correlation on `D_coph` — **no trace, obs ≈ surrogate**

#### Critical preamble (5-point)
1. **Claim.** δ `ρ_split^coph > 0` at cohort level.
2. **Null (C3, mandatory).** Strength-preserving 4-cycle ±δ rewiring, R=200, seed 20260511. Paired Wilcoxon one-sided greater across patients: obs > surrogate.
3. **Strongest plausible alternative.** At δ specifically: the published "anchor anatomy" cross-probe pattern (1.55× ratio) is a strong known-biology baseline; the matched-strength test must reject this as the explanation for any positive cohort-median signal.
4. **Does C3 cover (3).** Yes — matched-strength preserves node strengths and edge magnitudes, only topology shuffles. The anchor-anatomy pattern is anatomically distributed (driven by spatial proximity of probe contacts to epi-zone and to each other), so a per-pair matched-strength test correctly rejects it as a trace alternative. The cohort-median ρ_split^coph at δ is +0.031 — small but positive — and the matched-strength test directly evaluates whether this exceeds chance under strength-preserving topology shuffling.
5. **Falsification of the verdict.** Would require C3 paired Wilcoxon p < 0.05. Currently 0.278, far from threshold.

#### Results (δ, full cohort)

| Statistic | Value | Source |
|---|---|---|
| Cohort-median `ρ_split^coph` | +0.0313 | `ctm_triangle/cohort_summary.csv` delta row |
| n_trace_split | 6/10 | same |
| n_above_drift | 7/10 | same |
| C1 Wilcoxon p (`ρ_split > 0`) | **0.2158** | same col `wilcoxon_split_gt_0_p` — ✗ |
| C2 Wilcoxon p (`ρ_split > ρ_drift`) | 0.2461 | same col `wilcoxon_split_gt_drift_p` — ✗ |
| C3 matched-strength: obs_median | +0.00765 | `matched_strength_surrogate_split_baseline/cohort_summary.csv` delta |
| C3 surr median (median over surrogate medians) | +0.00782 | same |
| C3 ratio obs/surr | **0.98×** (obs essentially equal to surrogate) | derived |
| C3 paired Wilcoxon z / p | 34.0 / **0.2783** | same — ✗ |
| C3 n_above_surrogate (per-patient) | 4/10 | same |
| C4 cross-probe rho_xprobe | +0.0323 | `ctm_triangle/cohort_summary.csv` delta |
| C4 n_trace_xprobe | **6/10 (+sign matches +0.031)** | same — formally passes C4 gate |
| C5 epi-X cophenet | **not run** for δ (audit_68 covered α only) | — |

**C1, C2, C3 all fail.** C4 formally passes (6/10 +sign, sign-match with `rho_split`), but per LEDGER Decision 5 the C4 reading at δ is the **published δ cross-probe anchor anatomy** (1.55× cross-probe ratio, known epileptogenesis pattern from `memory/epileptic_imcoh_universal.md`) — *not* a positive trace claim. Under the locked rule, C3 mandatory failure means **no trace at cophenet** regardless of C4.

#### Reading
δ cophenet has **no trace**. The obs_median and surr_median are essentially identical (+0.0076 vs +0.0078), so the matched-strength test directly confirms that the cophenet representation of δ does not retain task-related per-pair reorganization. The substrate-level positive cohort median (+0.111 at raw FC, passing matched-strength at p=0.042) does not survive into the cophenetic dendrogram aggregation — at δ, the cophenet step **fully extinguishes** the substrate signal.

The C4 cross-probe formal pass is recorded for completeness and tagged in the brief as the published δ anchor-anatomy known biology, not as a verdict-changing observation. The δ cross-probe pattern reproduces the n=10 verdict from `memory/epileptic_imcoh_universal.md` (1.55× cross-probe ratio at δ, known biology of low-frequency synchronization near epileptogenic zones) — this is a **pipeline-validation observation** (the pipeline reads what it should read), not a discovery.

#### Cache + script provenance
- Full cohort: `audit_63_split_baseline_surrogate.py` → `data/audit/matched_strength_surrogate_split_baseline/`
- Within-baseline triangle: `audit_*_ctm_triangle.py` → `data/audit/ctm_triangle/`
- Object built from `lrg.ultrametric_matrix` field in `data/cache/imcoh_lrg/Pat_NN/delta_{phase}_lrg_imcoh-abs.npz`
- Anchor anatomy reference: `.agents/reports/2026-05-07_epileptic-n10-revisit.md` + `memory/epileptic_imcoh_universal.md`

### 3.2 Global subspace rotation: Grassmann chordal distance `d_G(k)` — **weak trace (LOO Pat_08 prevents strong-tier promotion)**

#### Critical preamble (5-point)
1. **Claim.** `T_G(k) < 0` at cohort level for some contiguous window of `k`.
2. **Null (C3, audit_70 cluster-extent permutation).** R=200 matched-strength surrogates; cluster_p = (1 + #(null_LR ≥ obs_LR)) / (R + 1). **Locked gate as of 2026-05-19**, replacing the previous 8-cell hardcoded threshold (under which δ was `no trace`).
3. **Strongest plausible alternative.** Same-strength edge rewiring producing a 7-cell contiguous cluster by chance. Under R=200 surrogates, the empirical null produces 5-cell runs at the 95th percentile, so a 7-cell observed run is above the 95th-percentile chance threshold.
4. **Does C3 cover (3).** Yes by construction — the cluster-extent null IS the matched-strength chance-cluster distribution.
5. **Falsification of the verdict.** Would require `cluster_p_mass ≥ 0.05` (currently 0.005 at floor). Promotion to "strong" tier would require full-data LOO max `p_mass < 0.05` (currently 0.055 Pat_08, just above gate — fails Decision-12 precondition).

#### Results (δ)

| Statistic | Value | Source |
|---|---|---|
| Observed longest contiguous-significant run | **7 cells** at k = 57..63 | `grassmann_cluster_extent/cohort_summary.csv` delta + per-k inspection |
| Cluster-extent null mean LR | 2.17 | `grassmann_cluster_extent/cohort_summary.csv` delta |
| Cluster-extent null median LR | 2.0 | same |
| Cluster-extent null 95th LR | **5.0** | same |
| Cluster-extent null max LR | 10 | same |
| Observed cluster mass `T_G^*` (raw Σ −log10p) | **38.07** (all-clusters formula, post-2026-05-19-pm fix) | same |
| Observed `T_G^*` (normalized per C1) | **0.149** (raw / 255.65, n_k=111, R=200) | same |
| Null mean cluster mass | 8.06 | same |
| Null 95th cluster mass | 18.37 | same |
| **cluster_p_mass (Decision-8 gate)** | **0.005** (at empirical floor `1/(R+1)`) | same |
| cluster_p_longest_run (descriptive co-statistic) | 0.025 | same |
| **LOO max p_mass (Decision-12 precondition)** | **0.055 (Pat_08) — FAILS < 0.05** | same — this is why the verdict stays "weak" |
| C5 epi-X (audit_72; secondary): raw mass | 38.07 → **43.99** (strengthens) | `grassmann_epi_exclusion/c5_wilcoxon_cohort.csv` delta |
| C5 epi-X cohort gate `p_mass^epi-X` | **0.005** (at floor) | same |
| C5 epi-X LOO max under epi-X | **0.005 (Pat_02)** — fully robust under epi-X | same |
| C5 epi-X obs longest run | **7 cells at k = 33..39** | `grassmann_epi_exclusion/sensitivity.csv` |
| C5 epi-X overlap with full-cohort window | 0% (k=57..63 vs k=33..39) — descriptive only | per-k inspection |

The 7-cell observed run sits clearly above the null 95th percentile (5.0) and well above the null mean (2.17). Under the Decision-8 mass gate, `cluster_p_mass = 0.005` is at the empirical floor — δ clears the cohort gate at the same nominal strength as β and γ_l. However, **full-data LOO max p_mass = 0.055 (Pat_08) fails the Decision-12 < 0.05 LOO precondition**: Pat_08 single-handedly leverages the cohort verdict over the gate at full data. By the locked Decision-12 rule, δ Grassmann is therefore **weak trace** (cohort gate passes; LOO robustness precondition fails). The cluster-extent revision (Decision 6) had promoted δ no trace → weak; Decision 12 confirms "weak" rather than further promoting to "strong".

#### C5 epi-X — secondary mechanistic observation (audit_72)
Under epi-zone exclusion (treated as a secondary mechanistic observation per Decision 10, **not** a verdict-promoter):
- Cohort gate strengthens: `p_mass^epi-X = 0.005` (floor); raw mass 38.07 → **43.99** (strengthens, opposite of γ_l which contracts).
- LOO under epi-X max: `p_mass^epi-X = 0.005 (Pat_02)` — **fully robust** under epi exclusion.
- Contiguous window shifts: full-cohort k=57..63 (7 cells) → epi-X k=33..39 (7 cells), 0% within-window overlap.

**Mechanistic interpretation (informative but not verdict-bearing).** The full-data Pat_08 leverage (LOO 0.055) resolves cleanly under epi-X (LOO 0.005, fully robust) — this strongly suggests that Pat_08's contribution to the full-cohort verdict was driven by epi-zone-coupled interactions rather than the true biological trace. Removing the epi-zone contacts unmasks a more robust subspace signature at a different mode range (k=33..39 vs k=57..63). The cohort-level signal not only survives epi exclusion but strengthens (mass 38.07 → 43.99) and resolves the single-patient leverage. **This is interesting biology to discuss but does not promote the full-data verdict to "strong"** — per Decision 10, C5 epi-X is a secondary mechanistic observation, never a verdict-driver. The full-data verdict remains "weak" because the *full-data* LOO at full data fails the Decision-12 precondition; C5 resolution is descriptive of *why* (epi-zone interactions), not a recipe for promotion.

#### Reading
δ Grassmann passes the Decision-8 cluster-mass gate at `cluster_p_mass = 0.005` (floor) with cluster-extent co-statistic `cluster_p_LR = 0.025`. The cohort gate is held decisively; what stops the verdict from being "strong" is the Decision-12 LOO precondition: full-data LOO max `p_mass = 0.055 (Pat_08)` fails < 0.05 — Pat_08 single-handedly leverages the cohort verdict. Under C5 epi-X (secondary, not verdict-driver) the trace strengthens (mass 38.07 → 43.99) and LOO resolves to 0.005 (Pat_02) fully robust, interpretable as Pat_08's leverage at full data being driven by epi-zone interactions rather than the true trace. **Final verdict per Decision 12: weak trace, only Grassmann.**

#### Cache + script provenance
- `audit_66_grassmann_matched_strength_surrogate.py` → `data/audit/grassmann_matched_strength_surrogate/`
- `audit_67_grassmann_epi_exclusion.py` → `data/audit/grassmann_epi_exclusion/`
- `audit_70_grassmann_cluster_extent.py` → `data/audit/grassmann_cluster_extent/`
- Object built from `eigvecs` field in `data/cache/imcoh_lrg/Pat_NN/delta_{phase}_lrg_imcoh-abs.npz` (obs) and `data/cache/matched_strength_surrogate_lrg/Pat_NN/delta_{phase}_R200_swap20_seed20260511_imcoh_abs.npz` (surrogate)

## 4. Multiscale structure of the δ trace

### Per-pair multiscale via `D_coph` (extinguished)

The δ cophenet probe is silent under C3 with obs ≈ surrogate. Substrate-level positive cohort signal (raw FC matched-strength passes at p=0.042, ratio 24.6×) does not propagate through the cophenetic dendrogram aggregation. This is the most complete "extinction" pattern in the panel — at γ_l the substrate signal at least retains positive sign through to cophenet (cohort-median +0.083 with ratio 30.2× to surrogate, even though paired test fails); at δ, obs and surrogate medians are within 3% of each other (+0.0076 vs +0.0078).

Interpretation: δ-band per-pair coupling is dominated by slow ambient activity (e.g., the slow waves of NREM-like states; the 1/f background; near-DC drift). The substrate-level cohort signal that passes matched-strength at the raw FC layer captures patient-level coupling strength differences, but the cophenetic representation re-distributes this across `N−1` merge heights without preserving the per-pair task-related component.

### Subspace multiscale via Grassmann `d_G(k)` (weak trace at high k; LOO Pat_08 binding)

The δ subspace trace at k=57..63 is at **the highest k of any band's contiguous-significant window**: β (k=27..55), γ_l (k=12..23), γ_h (k=19..27), δ (k=57..63). High-k modes are faster diffusion modes with smaller spatial coherence length — locally-extending rather than network-wide patterns. The δ Grassmann signal at k=57..63 is therefore a **mid-spectrum subspace rotation** rather than a slow-mode (low-k) rotation.

The C5 epi-X shift to k=33..39 brings the trace down to slower modes — these are more macroscopic patterns that emerge when the local high-k epi-zone-coupled modes are removed. Mechanistically, this is consistent with the epi-zone forming local coupling clusters that contribute to mid-spectrum modes; removing those clusters leaves the remaining network with a macroscopic δ-band coupling reorganization at the slow-mode level.

## 5. Anatomical distribution — ~~two FULLY DISJOINT networks~~ **RETRACTED 2026-05-30 / DELOCALIZED 2026-06-01**

> ⚠️ **RETRACTED — the δ Grassmann "two disjoint networks" claim does not survive.** Neither
> network localizes: **0/4 (full) and 0/3 (epi-X) cohort-supported** (the well-sampled locked
> regions are *anti*-localized: superiortemporal n=4 = −0.10, inferiortemporal n=5 = −0.03),
> and **per-patient localization is null** too. Two extra problems compound it: (a) the
> δ-epi-X "network" was computed with a **no-op epi-X mask bug** (`np.isin(int,str)` excluded
> 0 nodes — it never removed epileptic contacts; fixed 2026-05-30), so the "fully disjoint
> two networks" contrast was partly an artifact; (b) δ is only a **weak** trace to begin with.
> The region lists were top-decile per-node-participation + A3-on-counts on the phase-AVERAGED
> anchor quantity (not a cross-phase trace). Tables retained as historical record. Sources:
> `data/audit/anatomy_localization_wilcoxon/README.md`,
> `data/audit/per_patient_localization/README.md`, `locked/ANATOMY_LEDGER.md`.

~~The δ Grassmann trace is **localized to two anatomically distinct cortical networks** ... The dissociation is **fully disjoint** under the cluster-extent paradigm.~~ (retracted — see banner)

Audited under `locked/ANATOMY_CONTROLS.md`. **Verdict source NOW: `data/audit/anatomy_localization_wilcoxon/` + `data/audit/per_patient_localization/` (2026-05-30/06-01 retraction); `ANATOMY_LEDGER.md` 2026-05-19 pm entry superseded.**

### δ Grassmann anatomy, full cohort (over `S(δ)`)

Aggregation is the **unweighted average of per-node participation across `S(δ) = {k : p_k(δ) < α_k}`** — the support of the cluster-mass statistic `T_G^*`, 23 cells. The retired contiguous-significant window `K*(δ) = [57, 63]` (7 cells) is **not used** for anatomy under the locked paradigm.

Four named DK regions pass A3 alone (A1 sparse for Grassmann probe):

| Region | k_trace / K_region | A1 enrichment | A3 obs_z | A3 p_emp |
|---|---|---|---|---|
| ctx-lh-inferiortemporal | 4 / 29 | 1.34× | very large | 0.005 |
| ctx-lh-inferiorparietal | 2 / 5 | 3.89× | very large | 0.005 |
| ctx-rh-parstriangularis | 1 / 13 | 0.65× | very large | 0.005 |
| ctx-lh-superiortemporal | 2 / 32 | 0.61× | 14.11 | 0.010 |

~~Network anatomy: **temporal + parietal + frontal** — left inferior temporal + left superior temporal + left inferior parietal + right pars triangularis.~~ **RETRACTED (see §5 banner): 0/4 cohort-supported; well-sampled superiortemporal (n=4) and inferiortemporal (n=5) are *anti*-localized; per-patient null. Only passing locked region (inferiorparietal) is n=1.**

Source: `data/audit/anatomy_delta_grassmann_clusterext/cohort_summary.csv` (audit_72 --cluster-extent, 2026-05-19 pm).

### δ Grassmann anatomy, C5 epi-X (over `S^epiX(δ)`)

`S^epiX(δ)` is derived from per-k cohort-paired Wilcoxon on `data/audit/grassmann_epi_exclusion/per_patient_per_band_per_k.csv` (one-sided "less", α_k = 0.05), 22 cells with span `k = [2, 34–52, 87, 88]`. The retired `K*^epiX(δ) = [33, 39]` (7 cells) is **not used**.

Three named DK regions pass A3 alone:

| Region | k_trace / K_region | A1 enrichment | A3 obs_z | A3 p_emp |
|---|---|---|---|---|
| ctx-lh-superiorparietal | 1 / 1 | 9.73× | very large | 0.005 |
| ctx-rh-rostralmiddlefrontal | 5 / 33 | 1.47× | 2.71 | 0.005 |
| ctx-lh-superiorfrontal | 2 / 32 | 0.88× | 14.11 | 0.010 |

~~Network anatomy: **left superior parietal + right rostral middle frontal + left superior frontal**.~~ **RETRACTED (see §5 banner): 0/3 cohort-supported (rostralmiddlefrontal anti-localized); per-patient null; AND this network was computed with the no-op epi-X mask bug (never excluded epileptic nodes; fixed 2026-05-30).**

~~**Region-set overlap with full cohort**: **0/3 regions are shared**. The two networks are **fully disjoint** ...~~ **RETRACTED — the "fully disjoint two networks" contrast does not survive (neither network localizes; the epi-X side was a masking-bug artifact).**

Source: `data/audit/anatomy_delta_grassmann_epiX_clusterext/cohort_summary.csv` (audit_72 --epi-x --k-list ..., 2026-05-19 pm).

### Reading: ~~δ Grassmann is a mixture of two distinct phenomena~~ **RETRACTED — neither "phenomenon" localizes**

> ⚠️ **RETRACTED 2026-05-30/06-01.** The "two distinct, anatomically localized phenomena"
> reading does not survive: neither network is cohort-supported (0/4 full, 0/3 epi-X;
> well-sampled regions anti-localized) or per-patient-supported, and the epi-X "network"
> rested on a no-op masking bug (fixed). The correct reading is that the weak δ Grassmann
> trace, to the extent it exists, is **anatomically delocalized** — there is no "two-network"
> structure to report. The bullets below are retained only as historical record.

~~The two fully disjoint anatomies at fully disjoint significance-thresholded sets mean:~~
- ~~**Full-cohort signal** (over `S(δ)`): temporal + parietal + frontal.~~ (retracted)
- ~~**C5 epi-X signal** (over `S^epiX(δ)`): superior parietal + superior frontal + rostral middle frontal.~~ (retracted — masking-bug artifact)

~~This is **not** the same trace seen through two windows...~~ **RETRACTED — see banner.**

~~The δ Grassmann verdict ... the anatomy lockdown reveals that this weak trace is actually **two phenomena, both anatomically localized, both physiologically interpretable, with zero anatomical overlap**. The preprint should report both networks and note the strengthened dissociation.~~ **RETRACTED 2026-05-30/06-01 — neither "phenomenon" is anatomically localized (0/4 and 0/3 cohort-supported, per-patient null), and the epi-X side rested on a no-op masking bug. The δ Grassmann verdict stays "weak trace"; do NOT report a "two-network" anatomy. The trace, such as it is, is delocalized.**

### δ anchor anatomy (descriptive only, substrate-layer)

Separately, the δ anchor anatomy at the C4 cross-probe level is **already documented**:
- δ cross-probe ratio 1.55× (confirmed at n=10, `memory/epileptic_imcoh_universal.md`)
- Known biology: low-frequency synchronization near epileptogenic zones
- **NOT a discovery** — confirms the pipeline reads expected anchor structure
- **Substrate-layer only** (raw FC at C4 cross-probe restriction) — independent of the Grassmann anatomy result above

~~Under the locked cluster-extent paradigm, the Grassmann anatomy at δ full-cohort is **temporal + parietal + frontal**, NOT Amy + cingulate + OFC.~~ **RETRACTED 2026-05-30 (signed sampling-aware audit, locked in `ANATOMY_LEDGER.md`): the δ Grassmann DK-region localization does not survive — 0/4 cohort-supported, and the well-sampled superiortemporal (n=4) / inferiortemporal (n=5) are *anti*-localized — so there is no δ "temporal + parietal + frontal" localization to report (the Grassmann subspace probe localizes in no band).** What stands: the δ cross-probe anchor anatomy (1.55× ratio, known biology per `memory/epileptic_imcoh_universal.md`) remains a separate, descriptive substrate-layer observation (LEDGER Decision 5), independent of the (retracted) Grassmann anatomy.

### Caveats

- A1 hypergeometric is sparse for Grassmann (top-decile-per-patient endpoints). A3 is the gate.
- "Very large" A3 obs_z indicates all 200 matched-strength surrogate enrichments fell strictly below the observation.
- The two anatomies in this section come from two separate audits at two non-overlapping k-windows; neither directly establishes that the **same** patients drive both signals. A per-patient overlay (which patients contribute the most to the full-cohort window vs the epi-X window) would tighten this — descriptive supplement, not done.

### Cache + script provenance
| Artifact | Path |
|---|---|
| δ Grassmann anatomy (locked, `S(δ)` aggregation) | `data/audit/anatomy_delta_grassmann_clusterext/cohort_summary.csv` |
| δ Grassmann anatomy (locked, `S^epiX(δ)` aggregation) | `data/audit/anatomy_delta_grassmann_epiX_clusterext/cohort_summary.csv` |
| δ Grassmann anatomy (retired, `K*(δ) = [57,63]`) | `data/audit/anatomy_delta_grassmann/cohort_summary.csv` (do not cite) |
| δ Grassmann anatomy (retired, `K*^epiX(δ) = [33,39]`) | `data/audit/anatomy_delta_grassmann_epiX/cohort_summary.csv` (do not cite) |
| Producing script | `scripts/01_compute/audit/audit_72_anatomy_grassmann.py --band delta --cluster-extent [--epi-x --k-list ...]` |
| Anatomy battery doc | `.agents/preprint/locked/ANATOMY_CONTROLS.md` |
| Anatomy verdict ledger | `.agents/preprint/locked/ANATOMY_LEDGER.md` (2026-05-19 pm entry) |
| δ anchor anatomy memory (substrate-layer only) | `memory/epileptic_imcoh_universal.md` |

## 6. Patient-by-patient reading

### δ `ρ_split^coph` per-patient signs (descriptive — probe is null under C3)
From `ctm_triangle/cohort_summary.csv`: cohort median `ρ_split^coph` = +0.031 with 6/10 in + direction. The C3 paired Wilcoxon p = 0.278 has 4 patients above their own matched-strength surrogate. The 6/10 +sign agreement combined with the obs≈surrogate finding means the +sign cohort consensus does not represent task-related per-pair reorganization at the LRG layer.

### δ Grassmann subspace trace per-patient (descriptive)
Per-patient `n_patients_below_own_surrogate` at k=57..63 ranges from 6–9 patients (e.g., k=59: 9/10; k=57: 8/10; k=63: 7/10 — from `grassmann_matched_strength_surrogate/cohort_summary.csv`). Cohort agreement at the Grassmann window is therefore strong. Under epi-X at k=33..39 the cohort agreement reaches 8–9/10 at most cells.

### Pat_03 (1024 Hz) note
Pat_03 sits within the cohort-typical distribution at both probes. The δ band is well-resolved at 1024 Hz (Nyquist 512 Hz comfortably above the 4 Hz upper bound), so Pat_03's sampling-rate outlier status does not constrain the δ analysis.

## 7. Robustness panel (under the locked battery)

| Check | δ verdict | Source |
|---|---|---|
| C1 within-rsPre split-half null | ✗ 6/10 cohort, p = 0.216 | `ctm_triangle/cohort_summary.csv` delta |
| C2 drift-floor `ρ_split > ρ_drift` | ✗ 7/10 above drift, p = 0.246 | same |
| C3 matched-strength surrogate `ρ_split^coph` | ✗ ratio 0.98× (obs ≈ surr), 4/10, p = 0.278 | `matched_strength_surrogate_split_baseline` delta |
| C3 matched-strength surrogate Grassmann (cluster-extent) | ✓ cluster_p_mass = 0.005 (Decision-8 gate at floor); cluster_p_LR = 0.025 descriptive; **LOO max p_mass = 0.055 Pat_08 fails Decision-12 precondition** → verdict stays "weak" rather than promoting to "strong" | `grassmann_cluster_extent` delta |
| C4 cross-probe restriction | (✓) +0.032, 6/10 +sign matches (but reading is anchor-anatomy known biology, LEDGER Decision 5) | `ctm_triangle` delta (rho_xprobe column) |
| C5 epi-X `ρ_split^coph` | **not run** | (audit_68 was α-only) |
| C5 epi-X Grassmann (secondary observation per Decision 10) | strengthens: mass 38.07 → 43.99, p_mass^epi-X = 0.005 (floor), LOO under epi-X = 0.005 Pat_02 fully robust — full-data Pat_08 leverage attributable to epi-zone interactions, NOT a verdict-promoter | `grassmann_epi_exclusion/c5_wilcoxon_cohort.csv` delta |
| Pat_03 dropout | not decisive | inspection |
| Three-layer cohort table | substrate **passes** matched-strength → raw D borderline → cophenet **extinguishes** (obs ≈ surrogate); Grassmann rescues at high k | headline section above |
| Anchor anatomy (memory n=10 cross-probe 1.55×) | reproduced at C4; descriptive only per LEDGER Decision 5 | `memory/epileptic_imcoh_universal.md` |

## 8. Sensitivity panel under the locked 5-control battery

| Control | δ status (D_coph) | δ status (Grassmann) | Notes |
|---|---|---|---|
| C1 split | ✗ p = 0.216 | n/a | C1 specific to cophenet |
| C2 drift | ✗ p = 0.246 | n/a | C2 specific to cophenet |
| C3 matched-strength | ✗ p = 0.278 (obs ≈ surr, ratio 0.98×) | ✓ cluster_p_mass = 0.005 cohort gate (Decision-8 floor); LOO max p_mass = 0.055 Pat_08 fails Decision-12 precondition → "weak" not "strong" | mandatory |
| C4 cross-probe | (✓) +sign 6/10 — passes formally; reads as anchor-anatomy known biology (LEDGER Decision 5) | n/a | C4 specific to cophenet |
| C5 epi-X | **not run** (no audit_68 for δ cophenet) | strengthens: mass 38.07 → 43.99, LOO under epi-X 0.005 Pat_02 fully robust — secondary observation per Decision 10 | sensitivity layer; not verdict-driver |

**Verdict (locked, Decision 12)**: `weak trace, only Grassmann` per `locked/VERDICT_LEDGER.md` (revised 2026-05-19 Decision 6; Decision-12 cascade 2026-05-26 — δ cohort gate passes at floor but full-data LOO Pat_08 = 0.055 fails the Decision-12 < 0.05 LOO precondition, blocking promotion to "strong"). C4 anchor-anatomy reading is descriptive (LEDGER Decision 5).

## 9. Interpretation

**Neuroscientific framing.** δ (0.53–4 Hz) — the canonical slow-wave / deep-sleep oscillation band, also implicated in resting-state large-scale network organization (Steriade et al. 1993; He et al. 2008) and known to be elevated near epileptogenic zones in interictal recordings (Imamura et al. 2011; Ren et al. 2015). The known anchor-anatomy biology (cross-probe ratio 1.55× confirmed at n=10, memory `epileptic_imcoh_universal.md`) reflects the epileptogenesis-related δ synchronization pattern that ImCoh recovers under volume-conduction-immune assumptions.

**Two separate findings to report at δ.**

1. **Anchor anatomy (descriptive, confirmation of known biology).** The C4 cross-probe +sign at 6/10 reproduces the published δ "anchor" pattern. This is a *pipeline-validation observation*: the pipeline reads what it should read. Not a trace claim.

2. **Whole-network subspace trace at high `k` (weak, novel).** The Grassmann probe detects a 7-cell contiguous-significant window at k=57..63 (`cluster_p_mass = 0.005` at empirical floor, all-clusters formula, `T_G^* = 38.07` raw / 0.149 normalized). The cohort gate clears decisively, but full-data LOO max `p_mass = 0.055 (Pat_08)` fails the Decision-12 < 0.05 LOO precondition — preventing promotion to "strong" tier. Under C5 epi-X (secondary mechanistic observation per Decision 10) the trace strengthens (mass → 43.99, LOO → 0.005 fully robust), interpretable as the full-data Pat_08 leverage being driven by epi-zone interactions rather than the true biological signal.

**Why the substrate passes but cophenet doesn't.** This is the structurally interesting feature of δ: the raw FC matched-strength gate clears at p=0.042 (ratio 24.6×, 7/10 cohort), one of only two bands where the substrate passes (α also does). But the cophenet representation extinguishes the signal entirely (obs +0.0076 ≈ surrogate +0.0078). At α, the cophenet step *demotes* the substrate signal but C5 epi-X *promotes* it back (full cohort 5/10 → epi-X 7/10). At δ, the substrate signal does **not** re-emerge under any cophenet treatment — the dendrogram aggregation simply does not preserve δ per-pair task-related structure.

Mechanistically: at δ frequencies, per-pair coupling is dominated by spatially-extensive slow rhythms (slow waves; 1/f background). The cohort-median substrate-level reorganization captured by raw FC ρ_split is a **coupling-strength-distribution** difference between phases (which patient/pair is more strongly coupled). When ρ_split^coph is computed on the cophenetic image, this magnitude information is removed (cophenet records only merge heights, which depend on the relative ordering of pair distances), and the residual δ pattern is uniformly redistributed across merge scales. The full-cohort δ Grassmann trace lives in mid-spectrum modes (k=57..63) where local subspace rotation is detectable; this is precisely where epi-zone removal shifts the signal — the local coupling structure near epi-zones contributes to those modes, and removing them moves the trace to slower physiological modes.

**Why δ is `weak trace` rather than `strong trace` or `no trace` (Decision-12 cascade 2026-05-26).** The cohort gate is held: `cluster_p_mass = 0.005` at the empirical floor (Decision-8 mass gate) clears the < 0.01 strong threshold. **What prevents promotion to "strong" is the Decision-12 LOO precondition: full-data LOO max p_mass = 0.055 (Pat_08) > 0.05** — Pat_08 single-handedly leverages the cohort verdict over the gate at full data. The cluster-extent revision (Decision 6) promoted no trace → weak; Decision 12 confirms weak without further promotion. C5 epi-X strengthens decisively (mass 38.07 → 43.99, LOO under epi-X = 0.005 Pat_02 fully robust), reported as secondary mechanistic observation per Decision 10 — *not* a verdict-promoter. Cophenet fails fully (obs ≈ surrogate), so the trace exists only at the subspace level and the per-pair multiscale geometry is silent.

## 10. Source-of-truth references

### Lockdown documents
- `.agents/preprint/locked/CONTROLS.md` — 5-control battery (Grassmann gate at cluster-extent permutation)
- `.agents/preprint/locked/VERDICT_LEDGER.md` — locked verdict tag `weak trace, only Grassmann` for δ; Decision 5 (anchor anatomy descriptive) + Decision 6 (cluster-extent promotion) + Decision 12 (LOO+extent preconditions for strong tier; δ fails LOO at full data, stays weak)
- `.agents/preprint/methods/methods_revision_2026-05-18_cophenet.md` — binding methods directive
- `.agents/preprint/bands/01_beta.md`, `02_alpha.md`, `03_gammalow.md`, `04_theta.md`, `05_gammah.md` — companion band briefs

### Audit data (rows the numbers come from)
- `data/audit/ctm_triangle/cohort_summary.csv` — δ row: C1, C2, C4 statistics
- `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` — δ row: C3 cophenet (fails, obs ≈ surr)
- `data/audit/grassmann_cluster_extent/cohort_summary.csv` — δ cluster-extent permutation (audit_70, the promotion-determining audit)
- `data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv` — δ per-k T_G + surrogate stats (audit_66)
- `data/audit/grassmann_epi_exclusion/sensitivity.csv` — δ per-k under C5 epi-X (audit_67) — window shift k=57..63 → k=33..39
- `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv` — δ substrate row (passes matched-strength)
- `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv` — δ raw D sensitivity

### Audit scripts (provenance)
- `audit_63_split_baseline_surrogate.py` — C3 ρ_split^coph matched-strength (full cohort)
- `audit_66_grassmann_matched_strength_surrogate.py` — Grassmann C3 per-k
- `audit_67_grassmann_epi_exclusion.py` — Grassmann C5 per-k under epi-X (window shift at δ)
- `audit_70_grassmann_cluster_extent.py` — Grassmann cluster-extent permutation (2026-05-19; the promotion-determining audit for δ)

### Companion guides + memory
- `.agents/guides/02_methods/lrg-framework-guide.md` — LRG primitive
- `.agents/guides/02_methods/imcoh-guide.md` — imcoh_abs convention
- `.agents/guides/04_rules/never-always-list.md` — coding rules + no-hardcoded-thresholds feedback
- `.agents/reports/2026-05-07_epileptic-n10-revisit.md` — δ anchor anatomy n=10 verdict (1.55× cross-probe ratio, known biology)
- `memory/epileptic_imcoh_universal.md` — δ anchor anatomy known-biology framing

## 11. What follows

This brief freezes the δ `weak trace, only Grassmann` verdict for the preprint manuscript under the 2026-05-18 / 2026-05-19 lockdown. Next steps:

1. **δ figures** — produce F1 (4-control panel for `ρ_split^coph` showing all primary controls fail at cophenet; C4 +sign annotated as anchor-anatomy known biology), F2 (Grassmann `T_G(k)` cohort cells highlighting the 7-cell window k=57..63 and the audit_70 empirical null), F3 (per-patient subspace-rotation visualization at k=59 — the within-window k with strongest cohort agreement 9/10), and F4 (audit_72 c5_wilcoxon overlay showing the secondary epi-X analysis: mass strengthens 38.07 → 43.99, LOO under epi-X = 0.005 Pat_02 fully robust, window shifts k=57..63 → k=33..39). The F4 epi-X panel illustrates *why* the full-data LOO Pat_08 leverage fails — epi-zone interactions drive Pat_08's contribution at full data; removing them resolves the leverage. Output at `data/preprint/figures/delta/`.
2. **No anatomy audit owed** for the trace verdict. The δ anchor anatomy is already documented in `memory/epileptic_imcoh_universal.md` as known biology; a mode-decomposition anatomy of the k=57..63 window would be a descriptive supplement.
3. **No further sensitivity tests owed under the locked battery.** The C5 epi-X Grassmann audit (audit_72 c5_wilcoxon) is the secondary mechanistic observation per Decision 10 for δ; cophenet C5 epi-X was not run (audit_68 was α-only). The verdict is determined by cophenet C3 failure + Grassmann Decision-8 cohort gate pass + Decision-12 LOO precondition failure at full data.

Verdict ready for writing-agent handoff: **`weak trace, only Grassmann`** (Decision-12 cascade 2026-05-26 — cohort gate passes at floor but full-data LOO Pat_08 = 0.055 fails the < 0.05 LOO precondition). The headline story has three components: (i) **substrate passes but cophenet extinguishes** (one of two bands where substrate clears matched-strength along with α; cophenet aggregation does not preserve δ per-pair signal); (ii) **Grassmann cohort gate clears at floor but LOO Pat_08 prevents strong-tier promotion** (`cluster_p_mass = 0.005`, `T_G^* = 38.07` raw / 0.149 normalized at k=57..63; full-data LOO max 0.055 Pat_08 fails Decision-12 precondition); (iii) **C5 epi-X strengthens and resolves LOO** (mass 38.07 → 43.99, LOO under epi-X = 0.005 Pat_02 fully robust — full-data Pat_08 leverage attributable to epi-zone interactions, not the true trace; secondary mechanistic observation per Decision 10, not a verdict-promoter). The δ anchor-anatomy cross-probe 1.55× reading is reported separately as descriptive confirmation of known epileptogenesis biology (LEDGER Decision 5), not as a trace claim.
