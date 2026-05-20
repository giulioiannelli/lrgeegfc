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
  grassmann: weak_trace (cluster-extent permutation p=0.0249, audit_70; 7-cell observed run at k=57..63 — above null 95th percentile of 5.0; promoted from no trace under previous 8-cell hardcoded threshold)
  grassmann_epi_excluded: shifts_window_decisively (audit_67; 7-cell longest run at k=33..39 vs full-cohort k=57..63 — 0% within-window overlap, the epi-X signal is at a completely different mode range than the full-cohort signal)
  anchor_anatomy_known_biology: descriptive_only (C4 passes formally at +0.032/6/10 +sign — but this is the published δ "anchor" cross-probe ratio 1.55× known epileptogenesis pattern, NOT a positive trace claim — LEDGER Decision 5)
  anchor_anatomy_descriptive: documented (δ anchor anatomy is the n=10 verdict in memory/epileptic_imcoh_universal.md — cross-probe ratio 1.55× is known biology, not a discovery)
  anatomy_grassmann_full: strong_localized (4 named DK regions under `S(δ)`; temporal + parietal + frontal: left inferior temporal + left inferior parietal + right pars triangularis + left superior temporal; audit_72 --cluster-extent under `S(δ)`, 2026-05-19 pm — supersedes the retired `K*(δ)` audit which had picked up Amy + cingulate + fusiform + OFC + bankssts)
  anatomy_grassmann_epiX: strong_localized_FULLY_DISJOINT_NETWORK (3 named regions under `S^epiX(δ)` with **0 shared with full cohort under `S(δ)`**; left superior parietal + right rostral middle frontal + left superior frontal; audit_72 --cluster-extent --epi-x --k-list, 2026-05-19 pm — supersedes the retired `K*^epiX(δ)` audit; dissociation strengthens from 2/6 shared to fully disjoint)
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
---

# δ band (0.53–4 Hz) — preprint result report

## Head

δ carries a **whole-network subspace trace** detected by the Grassmann probe under cluster-extent permutation (audit_70 cluster_p = 0.0249) at a 7-cell contiguous-significant window k=57..63 — well above the empirical null 95th percentile of 5.0. **No per-pair cophenet trace**: `ρ_split^coph` C3 paired Wilcoxon p = 0.278, with `obs_median = +0.0076` essentially equal to `surr_median = +0.0078` (ratio 0.98×) and only 4/10 patients above their own surrogate. The substrate-level raw FC at δ **passes** matched-strength (cohort-median +0.111, ratio 24.6×, p=0.042 at 7/10) — one of only two bands where the substrate alone clears the matched-strength gate (α is the other). Under C5 epi-X (audit_67), the Grassmann window **shifts decisively** to k=33..39 (7-cell run) with **zero within-window overlap** vs the full-cohort k=57..63 window — the epi-X signal is at a completely different mode range than the full-cohort signal. This is descriptively informative (suggests the full-cohort signal is at least partially epi-zone-driven; the epi-X signal is a separate physiological subspace pattern) but does not change the locked verdict. The δ C4 cross-probe trace (+0.032 at 6/10 +sign) **formally passes the C4 gate** but is the published δ "anchor anatomy" 1.55× cross-probe ratio (known epileptogenesis pattern per `memory/epileptic_imcoh_universal.md`), **not** a trace claim — LEDGER Decision 5. Verdict from `locked/VERDICT_LEDGER.md`: **`weak trace, only Grassmann`** (promoted from no trace under the previous 8-cell hardcoded threshold via cluster-extent permutation, Decision 6).

## Headline three-layer cohort table (δ-band row, copied from `locked/VERDICT_LEDGER.md`)

```
Layer                                       δ (full n=10)
-------------------------------------      ---------------------
Raw FC          ρ_split^raw  (raw_fc_ms)   7/10 p=.042  ratio 24.6×  ← PASSES substrate
Raw D(τ_max)    ρ_split      (preprint_05) 7/10 p=.116  ratio 11.5×  ← borderline (just outside)
Cophenet D_coph ρ_split^coph (audit_63)    4/10 p=.278  ratio 0.98×  ← FAILS (obs ≈ surrogate)
Grassmann       d_G(k)       (audit_70)    7-cell run k=57..63, cluster_p = 0.0249  ← weak trace
```

Sources: `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv`, `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv`, `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv`, `data/audit/grassmann_cluster_extent/cohort_summary.csv`.

**Key reading**: δ is the **band where the substrate signal exists but does not propagate into per-pair LRG**, and the Grassmann probe rescues a subspace signature at high `k` (60-ish modes). The substrate-to-cophenet trajectory at δ shows substrate passing → raw D borderline → cophenet fully extinguished (obs and surrogate essentially identical). The LRG cophenet step does *not* preserve the substrate signal at δ. The Grassmann probe rescues a subspace signature at k=57..63, but the C5 epi-X shift to k=33..39 (0% within-window overlap) suggests the full-cohort high-k signal is at least partially driven by epi-zone contacts, while the epi-X-revealed lower-k signal is a separate physiological pattern. The previous 8-cell hardcoded threshold had δ as `no trace`; the cluster-extent permutation (audit_70, 2026-05-19) promoted it to `weak trace` because the 7-cell observed run is well above the null 95th percentile (5.0).

## 1. Scientific claim

**Cohort-level question.** Does the δ-band (0.53–4 Hz) post-task resting state at the LRG-CTM layer carry task-related reorganization on either of the two locked probes? Separately: does the δ-band cross-probe anchor anatomy reproduce the published n=10 verdict (1.55× cross-probe ratio, known epileptogenesis biology)?

**Refined biological claim.** δ coupling reorganizes during the task at a level visible in the raw FC substrate (cohort-median +0.111 at 7/10 patients, matched-strength p=0.042), but the per-pair LRG cophenet representation **does not retain** this signal (obs_median +0.0076 essentially equal to surrogate +0.0078). The Grassmann probe detects a weak subspace trace at high `k` (k=57..63, 7-cell window, cluster_p = 0.0249) that is **partially epi-zone-driven** — under epi-X, the window shifts to k=33..39 with 0% overlap, suggesting the slow-mode subspace rotation at the full cohort scale incorporates epi-zone contributions that mask a separate physiological lower-k rotation. The separate δ cross-probe anchor anatomy (1.55× cross-probe ratio, +sign at 6/10) is a **known epileptogenesis biology** (per memory `epileptic_imcoh_universal.md`) and is not a positive trace claim — its inclusion in the per-band brief is *descriptive*, confirming the pipeline reads the expected anchor structure (LEDGER Decision 5).

**Falsification budget.** The Grassmann verdict (`weak trace`) would flip to `no trace` if a re-run with a tighter null (R ≥ 500) raises the null 95th percentile LR above 7, or if the 7-cell window is shown to be entirely epi-zone-driven (audit_67 partial: the window shifts but the count remains 7). The cophenet verdict (`no trace`) would flip to `weak trace` only if a re-audit shows obs_median substantially above surrogate (currently they are within 3% of each other).

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

### 3.2 Global subspace rotation: Grassmann chordal distance `d_G(k)` — **weak trace (cluster-extent promotion)**

#### Critical preamble (5-point)
1. **Claim.** `T_G(k) < 0` at cohort level for some contiguous window of `k`.
2. **Null (C3, audit_70 cluster-extent permutation).** R=200 matched-strength surrogates; cluster_p = (1 + #(null_LR ≥ obs_LR)) / (R + 1). **Locked gate as of 2026-05-19**, replacing the previous 8-cell hardcoded threshold (under which δ was `no trace`).
3. **Strongest plausible alternative.** Same-strength edge rewiring producing a 7-cell contiguous cluster by chance. Under R=200 surrogates, the empirical null produces 5-cell runs at the 95th percentile, so a 7-cell observed run is above the 95th-percentile chance threshold.
4. **Does C3 cover (3).** Yes by construction — the cluster-extent null IS the matched-strength chance-cluster distribution.
5. **Falsification of the verdict.** Would require cluster_p ≥ 0.05. Currently 0.0249 (longest run) and 0.0249 (cluster mass), both clearly inside the gate.

#### Results (δ)

| Statistic | Value | Source |
|---|---|---|
| Observed longest contiguous-significant run | **7 cells** at k = 57..63 | `grassmann_cluster_extent/cohort_summary.csv` delta + per-k inspection |
| Cluster-extent null mean LR | 2.17 | `grassmann_cluster_extent/cohort_summary.csv` delta |
| Cluster-extent null median LR | 2.0 | same |
| Cluster-extent null 95th LR | **5.0** | same |
| Cluster-extent null max LR | 10 | same |
| Observed cluster mass (Σ −log10p) | 12.78 | same |
| Null mean cluster mass | 4.01 | same |
| Null 95th cluster mass | 9.81 | same |
| **cluster_p_longest_run** | **0.0249** | same — ✓ (weak gate) |
| **cluster_p_cluster_mass** | **0.0249** | same — ✓ |
| C5 epi-X (audit_67): obs longest run | **7 cells at k = 33..39** | `grassmann_epi_exclusion/sensitivity.csv` delta rows |
| C5 epi-X: within-window overlap with full-cohort | **0% (no overlap between k=57..63 and k=33..39)** | per-k inspection |

The 7-cell observed run sits clearly above the null 95th percentile (5.0) and well above the null mean (2.17). Both cluster-p statistics (longest-run and cluster-mass) are at 0.0249 — within the `weak` gate [0.01, 0.05). By the locked rule, δ Grassmann is **weak trace** (promoted from `no trace` under the cluster-extent revision, Decision 6).

#### C5 epi-X sensitivity (audit_67) — decisive window shift
Under epi-zone exclusion, the contiguous-significant window shifts:
- Original full-cohort: k=57..63 (7-cell run at p<0.05 gate).
- Epi-X: k=33..39 (7-cell run at p<0.05 gate, with all cells flagged "persist" or "emerge").

**Within-window overlap: 0%** — the full-cohort window k=57..63 and the epi-X window k=33..39 do not intersect at any `k`. This is the **most extreme window shift** in the panel (γ_l: partial overlap k=19..23, 42% overlap; γ_h: partial overlap, 67% overlap; δ: zero overlap).

Mechanistic reading: the full-cohort signal at k=57..63 is at least partially driven by epi-zone contacts contributing to high-`k` (faster) modes. Removing those contacts unmasks a **separate physiological signal at lower k=33..39** (slower modes, more macroscopic). These are two qualitatively different subspace rotations:
- **Full-cohort k=57..63**: roughly modes φ_58..φ_64 — these are mid-spectrum modes with intermediate spatial scale; the presence of epi-zone contacts in the network biases these modes toward epi-zone-related coupling.
- **Epi-X k=33..39**: roughly modes φ_34..φ_40 — these are slower (more macroscopic) modes that emerge when the epi-zone contributions are removed.

Both observations are 7-cell windows above their respective empirical 95th percentiles — but they describe **different subspace patterns**.

#### Reading
δ Grassmann passes cluster-extent at `weak` strength (cluster_p = 0.0249). The C5 epi-X window shift to k=33..39 with 0% overlap is descriptively the most interesting feature of the δ verdict — it suggests the δ subspace trace has two superposed components (an epi-zone-driven high-k component visible in the full cohort, and a physiological lower-k component revealed under epi-X) that are independent of each other. Both components individually satisfy the 7-cell-above-95th-percentile criterion, but the locked verdict is based on the full-cohort window (k=57..63, cluster_p = 0.0249).

#### Cache + script provenance
- `audit_66_grassmann_matched_strength_surrogate.py` → `data/audit/grassmann_matched_strength_surrogate/`
- `audit_67_grassmann_epi_exclusion.py` → `data/audit/grassmann_epi_exclusion/`
- `audit_70_grassmann_cluster_extent.py` → `data/audit/grassmann_cluster_extent/`
- Object built from `eigvecs` field in `data/cache/imcoh_lrg/Pat_NN/delta_{phase}_lrg_imcoh-abs.npz` (obs) and `data/cache/matched_strength_surrogate_lrg/Pat_NN/delta_{phase}_R200_swap20_seed20260511_imcoh_abs.npz` (surrogate)

## 4. Multiscale structure of the δ trace

### Per-pair multiscale via `D_coph` (extinguished)

The δ cophenet probe is silent under C3 with obs ≈ surrogate. Substrate-level positive cohort signal (raw FC matched-strength passes at p=0.042, ratio 24.6×) does not propagate through the cophenetic dendrogram aggregation. This is the most complete "extinction" pattern in the panel — at γ_l the substrate signal at least retains positive sign through to cophenet (cohort-median +0.083 with ratio 30.2× to surrogate, even though paired test fails); at δ, obs and surrogate medians are within 3% of each other (+0.0076 vs +0.0078).

Interpretation: δ-band per-pair coupling is dominated by slow ambient activity (e.g., the slow waves of NREM-like states; the 1/f background; near-DC drift). The substrate-level cohort signal that passes matched-strength at the raw FC layer captures patient-level coupling strength differences, but the cophenetic representation re-distributes this across `N−1` merge heights without preserving the per-pair task-related component.

### Subspace multiscale via Grassmann `d_G(k)` (weak trace at high k)

The δ subspace trace at k=57..63 is at **the highest k of any band's contiguous-significant window**: β (k=27..55), γ_l (k=12..23), γ_h (k=19..27), δ (k=57..63). High-k modes are faster diffusion modes with smaller spatial coherence length — locally-extending rather than network-wide patterns. The δ Grassmann signal at k=57..63 is therefore a **mid-spectrum subspace rotation** rather than a slow-mode (low-k) rotation.

The C5 epi-X shift to k=33..39 brings the trace down to slower modes — these are more macroscopic patterns that emerge when the local high-k epi-zone-coupled modes are removed. Mechanistically, this is consistent with the epi-zone forming local coupling clusters that contribute to mid-spectrum modes; removing those clusters leaves the remaining network with a macroscopic δ-band coupling reorganization at the slow-mode level.

## 5. Anatomical distribution — δ Grassmann trace localizes to two FULLY DISJOINT networks (full cohort vs C5 epi-X)

The δ Grassmann trace is **localized to two anatomically distinct cortical networks** depending on whether epi-zone contacts are included. Under the locked all-clusters paradigm (`S(δ)` aggregation, `data/audit/anatomy_delta_grassmann{_epiX}_clusterext/`), the full-cohort and C5 epi-X networks **share 0 regions** (was 2/6 under the retired `K*(δ)`/`K*^epiX(δ)` windows). The dissociation is **fully disjoint** under the cluster-extent paradigm.

Audited under `locked/ANATOMY_CONTROLS.md`. Verdict source: `locked/ANATOMY_LEDGER.md` 2026-05-19 pm cluster-extent revision.

### δ Grassmann anatomy, full cohort (over `S(δ)`)

Aggregation is the **unweighted average of per-node participation across `S(δ) = {k : p_k(δ) < α_k}`** — the support of the cluster-mass statistic `T_G^*`, 23 cells. The retired contiguous-significant window `K*(δ) = [57, 63]` (7 cells) is **not used** for anatomy under the locked paradigm.

Four named DK regions pass A3 alone (A1 sparse for Grassmann probe):

| Region | k_trace / K_region | A1 enrichment | A3 obs_z | A3 p_emp |
|---|---|---|---|---|
| ctx-lh-inferiortemporal | 4 / 29 | 1.34× | very large | 0.005 |
| ctx-lh-inferiorparietal | 2 / 5 | 3.89× | very large | 0.005 |
| ctx-rh-parstriangularis | 1 / 13 | 0.65× | very large | 0.005 |
| ctx-lh-superiortemporal | 2 / 32 | 0.61× | 14.11 | 0.010 |

Network anatomy: **temporal + parietal + frontal** — left inferior temporal + left superior temporal + left inferior parietal + right pars triangularis. The Amygdala, medial-OFC, caudal anterior cingulate, fusiform, and bankssts that appeared in the retired `K*(δ)` analysis are **not present under `S(δ)`** — the all-clusters aggregation reads a different network. The "anchor anatomy" interpretation (Amy + cingulate + fusiform overlap with published δ epi-zone synchronization) is **not supported** at the Grassmann probe under the locked paradigm.

Source: `data/audit/anatomy_delta_grassmann_clusterext/cohort_summary.csv` (audit_72 --cluster-extent, 2026-05-19 pm).

### δ Grassmann anatomy, C5 epi-X (over `S^epiX(δ)`)

`S^epiX(δ)` is derived from per-k cohort-paired Wilcoxon on `data/audit/grassmann_epi_exclusion/per_patient_per_band_per_k.csv` (one-sided "less", α_k = 0.05), 22 cells with span `k = [2, 34–52, 87, 88]`. The retired `K*^epiX(δ) = [33, 39]` (7 cells) is **not used**.

Three named DK regions pass A3 alone:

| Region | k_trace / K_region | A1 enrichment | A3 obs_z | A3 p_emp |
|---|---|---|---|---|
| ctx-lh-superiorparietal | 1 / 1 | 9.73× | very large | 0.005 |
| ctx-rh-rostralmiddlefrontal | 5 / 33 | 1.47× | 2.71 | 0.005 |
| ctx-lh-superiorfrontal | 2 / 32 | 0.88× | 14.11 | 0.010 |

Network anatomy: **left superior parietal + right rostral middle frontal + left superior frontal**. Fusiform, inferior parietal, postcentral, and inferior temporal that appeared in the retired `K*^epiX(δ)` analysis are **not present under `S^epiX(δ)`**; left superior frontal emerges as new.

**Region-set overlap with full cohort**: **0/3 regions are shared**. The two networks are **fully disjoint** under the cluster-extent paradigm (was 2/6 shared under the retired windows; the dissociation strengthens to "completely non-overlapping").

Source: `data/audit/anatomy_delta_grassmann_epiX_clusterext/cohort_summary.csv` (audit_72 --epi-x --k-list ..., 2026-05-19 pm).

### Reading: δ Grassmann is a mixture of two distinct phenomena (strengthened under `S(b)`)

The two fully disjoint anatomies at fully disjoint significance-thresholded sets mean:
- **Full-cohort signal** (over `S(δ)`): temporal + parietal + frontal — a cortical pattern that does NOT overlap with the published δ epi-zone anchor anatomy.
- **C5 epi-X signal** (over `S^epiX(δ)`): superior parietal + superior frontal + rostral middle frontal — a tighter parietal-frontal cortical signal revealed only when epi contacts are removed.

This is **not** the same trace seen through two windows. The full-cohort δ Grassmann probe and the C5 epi-X probe each read a distinct cortical reorganization pattern; the two patterns share no DK regions under the locked methodology.

The δ Grassmann verdict from `locked/VERDICT_LEDGER.md` is "weak trace" (cluster_p=0.0249); the anatomy lockdown reveals that this weak trace is actually **two phenomena, both anatomically localized, both physiologically interpretable, with zero anatomical overlap**. The preprint should report both networks and note the strengthened dissociation.

### δ anchor anatomy (descriptive only, substrate-layer)

Separately, the δ anchor anatomy at the C4 cross-probe level is **already documented**:
- δ cross-probe ratio 1.55× (confirmed at n=10, `memory/epileptic_imcoh_universal.md`)
- Known biology: low-frequency synchronization near epileptogenic zones
- **NOT a discovery** — confirms the pipeline reads expected anchor structure
- **Substrate-layer only** (raw FC at C4 cross-probe restriction) — independent of the Grassmann anatomy result above

Under the locked cluster-extent paradigm, the Grassmann anatomy at δ full-cohort is **temporal + parietal + frontal**, NOT Amy + cingulate + OFC. The "anchor anatomy" interpretation belongs to the substrate / C4 cross-probe discussion only and is not supported at the LRG-subspace anatomy layer. LEDGER Decision 5 logs the C4 cross-probe known biology; the cluster-extent revision (2026-05-19 pm) decouples the δ Grassmann anatomy from this anchor framing.

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
| C3 matched-strength surrogate Grassmann (cluster-extent) | ✓ cluster_p = 0.0249, 7-cell run k=57..63 (**promoted from no trace** under cluster-extent revision) | `grassmann_cluster_extent` delta |
| C4 cross-probe restriction | (✓) +0.032, 6/10 +sign matches (but reading is anchor-anatomy known biology, LEDGER Decision 5) | `ctm_triangle` delta (rho_xprobe column) |
| C5 epi-X `ρ_split^coph` | **not run** | (audit_68 was α-only) |
| C5 epi-X Grassmann | 7-cell run at k=33..39 (0% within-window overlap with full-cohort k=57..63) | `grassmann_epi_exclusion/sensitivity.csv` delta |
| Pat_03 dropout | not load-bearing | inspection |
| Three-layer cohort table | substrate **passes** matched-strength → raw D borderline → cophenet **extinguishes** (obs ≈ surrogate); Grassmann rescues at high k | headline section above |
| Anchor anatomy (memory n=10 cross-probe 1.55×) | reproduced at C4; descriptive only per LEDGER Decision 5 | `memory/epileptic_imcoh_universal.md` |

## 8. Sensitivity panel under the locked 5-control battery

| Control | δ status (D_coph) | δ status (Grassmann) | Notes |
|---|---|---|---|
| C1 split | ✗ p = 0.216 | n/a | C1 specific to cophenet |
| C2 drift | ✗ p = 0.246 | n/a | C2 specific to cophenet |
| C3 matched-strength | ✗ p = 0.278 (obs ≈ surr, ratio 0.98×) | ✓ cluster_p = 0.0249 (**promoted from no trace** under cluster-extent revision) | mandatory |
| C4 cross-probe | (✓) +sign 6/10 — passes formally; reads as anchor-anatomy known biology (LEDGER Decision 5) | n/a | C4 specific to cophenet |
| C5 epi-X | **not run** (no audit_68 for δ cophenet) | window shifts decisively k=57..63 → k=33..39 (0% overlap) | sensitivity layer; descriptive |

**Verdict (locked)**: `weak trace, only Grassmann` per `locked/VERDICT_LEDGER.md` (revised 2026-05-19, Decision 6). C4 anchor-anatomy reading is descriptive (LEDGER Decision 5).

## 9. Interpretation

**Neuroscientific framing.** δ (0.53–4 Hz) — the canonical slow-wave / deep-sleep oscillation band, also implicated in resting-state large-scale network organization (Steriade et al. 1993; He et al. 2008) and known to be elevated near epileptogenic zones in interictal recordings (Imamura et al. 2011; Ren et al. 2015). The known anchor-anatomy biology (cross-probe ratio 1.55× confirmed at n=10, memory `epileptic_imcoh_universal.md`) reflects the epileptogenesis-related δ synchronization pattern that ImCoh recovers under volume-conduction-immune assumptions.

**Two separate findings to report at δ.**

1. **Anchor anatomy (descriptive, confirmation of known biology).** The C4 cross-probe +sign at 6/10 reproduces the published δ "anchor" pattern. This is a *pipeline-validation observation*: the pipeline reads what it should read. Not a trace claim.

2. **Whole-network subspace trace at high `k` (weak, novel).** The Grassmann probe detects a 7-cell contiguous-significant window at k=57..63 (cluster_p = 0.0249) — promoted from `no trace` under the cluster-extent permutation revision. The epi-X window shift to k=33..39 with 0% overlap suggests a **two-component** δ subspace signature: an epi-zone-driven mid-spectrum (k=57..63) component visible in the full cohort, and a physiological lower-k (k=33..39) component revealed under epi-X. Both are 7-cell windows above their respective empirical 95th percentiles.

**Why the substrate passes but cophenet doesn't.** This is the structurally interesting feature of δ: the raw FC matched-strength gate clears at p=0.042 (ratio 24.6×, 7/10 cohort), one of only two bands where the substrate passes (α also does). But the cophenet representation extinguishes the signal entirely (obs +0.0076 ≈ surrogate +0.0078). At α, the cophenet step *demotes* the substrate signal but C5 epi-X *promotes* it back (full cohort 5/10 → epi-X 7/10). At δ, the substrate signal does **not** re-emerge under any cophenet treatment — the dendrogram aggregation simply does not preserve δ per-pair task-related structure.

Mechanistically: at δ frequencies, per-pair coupling is dominated by spatially-extensive slow rhythms (slow waves; 1/f background). The cohort-median substrate-level reorganization captured by raw FC ρ_split is a **coupling-strength-distribution** difference between phases (which patient/pair is more strongly coupled). When ρ_split^coph is computed on the cophenetic image, this magnitude information is removed (cophenet records only merge heights, which depend on the relative ordering of pair distances), and the residual δ pattern is uniformly redistributed across merge scales. The full-cohort δ Grassmann trace lives in mid-spectrum modes (k=57..63) where local subspace rotation is detectable; this is precisely where epi-zone removal shifts the signal — the local coupling structure near epi-zones contributes to those modes, and removing them moves the trace to slower physiological modes.

**Why δ is `weak trace` rather than `strong trace` or `no trace`.** Three reasons: (i) cluster_p = 0.0249 is well within [0.01, 0.05), the `weak` band; (ii) the epi-X 0% within-window overlap argues against a strong physiologically-coherent trace at full cohort — the signal mode-identity is unstable to epi-zone removal; (iii) cophenet fails fully (obs ≈ surrogate), so the trace exists only at the subspace level and the per-pair multiscale geometry is silent.

## 10. Source-of-truth references

### Lockdown documents
- `.agents/preprint/locked/CONTROLS.md` — 5-control battery (Grassmann gate at cluster-extent permutation)
- `.agents/preprint/locked/VERDICT_LEDGER.md` — locked verdict tag `weak trace, only Grassmann` for δ; Decision 5 (anchor anatomy descriptive) + Decision 6 (cluster-extent promotion)
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

1. **δ figures** — produce F1 (4-control panel for `ρ_split^coph` showing all primary controls fail at cophenet; C4 +sign annotated as anchor-anatomy known biology), F2 (Grassmann `T_G(k)` cohort cells highlighting the 7-cell window k=57..63 and the audit_70 empirical null), F3 (per-patient subspace-rotation visualization at k=59 — the within-window k with strongest cohort agreement 9/10), and F4 (audit_67 epi-X overlay showing the window shift k=57..63 → k=33..39 with 0% overlap). The F4 epi-X panel is **especially informative** at δ because of the 0% overlap. Output at `data/preprint/figures/delta/`.
2. **No anatomy audit owed** for the trace verdict. The δ anchor anatomy is already documented in `memory/epileptic_imcoh_universal.md` as known biology; a mode-decomposition anatomy of the k=57..63 window would be a descriptive supplement.
3. **No further sensitivity tests owed under the locked battery.** The C5 epi-X Grassmann audit (audit_67) is the sensitivity layer for δ; cophenet C5 epi-X was not run (audit_68 was α-only). The verdict is determined by cophenet C3 failure + Grassmann cluster-extent + descriptive epi-X window shift.

Verdict ready for writing-agent handoff: **`weak trace, only Grassmann`**. The headline story has three components: (i) **substrate passes but cophenet extinguishes** (one of two bands where substrate clears matched-strength along with α; cophenet aggregation does not preserve δ per-pair signal); (ii) **cluster-extent promotes Grassmann no trace → weak trace** (7-cell run at k=57..63, cluster_p = 0.0249, audit_70 Decision 6); (iii) **C5 epi-X window shift 0% overlap** (full-cohort k=57..63 vs epi-X k=33..39 — two-component subspace signature with epi-zone-driven high-k and physiological lower-k components). The δ anchor-anatomy cross-probe 1.55× reading is reported separately as descriptive confirmation of known epileptogenesis biology (LEDGER Decision 5), not as a trace claim.
