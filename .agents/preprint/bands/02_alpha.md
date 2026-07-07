---
name: preprint-alpha-band
era: IMCOH_ABS_COHORT_N10
status: current
kind: preprint-result-report
band: alpha
range_hz: [8, 13]
cohort: [Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15]
n_patients: 10
fc_method: imcoh_abs
tau: 1/lambda_max
canonical_lrg_object: "D_coph = cophenet(UPGMA(D(tau_max)))"
canonical_subspace_object: "U_k = span{phi_2, ..., phi_{k+1}}"
verdict_tag: "strong trace, only D_coph"
verdict_source: VERDICT_LEDGER.md (locked 2026-05-18, revised 2026-05-19)
verdict_layers:
  substrate_rank: trace_direction (raw FC ρ_split^raw matched-strength p=0.014 at 8/10 cohort — strongest substrate signal of any band)
  rho_split_coph: strong_trace (anchored at C3 alone: paired Wilcoxon p=0.002 ratio 8.3×; C1 p=.010, C2 p=.007, C4 paired-Wilcoxon non-degradation all pass at full cohort; per-patient n_above_surrogate 5/10 is descriptive auxiliary statistic, NOT a verdict modifier — patient-count thresholds retired 2026-05-19, reaffirmed 2026-05-28)
  grassmann: no_trace (cluster-extent permutation p=0.159, audit_70; 4-cell observed run at k=11..14 within null distribution)
  rho_split_coph_epi_excluded: SECONDARY_alpha_recruits_epi_core (secondary tissue-characterization / sensitivity layer parallel to C6 — the PI directive 2026-06-05 elevation to a "primary interpretive lens" is RETRACTED + DEMOTED 2026-06-18, symmetric with C6. All numbers on the cophenetic object D_coph — raw |ImCoh| is the comparison baseline ONLY, never a result. Same-graph pair classes (no node-count confound), re-run at n=10 [Pat_15 SOZ-labelled → 17 epi nodes; audit_77 regenerated, NO matched-strength flag flips vs n=9] + the PAIR-COUNT control [audit_85 --mode pairclass, R=1000]. **α EXCEPTION CONFIRMED — α uniquely RECRUITS the diseased seizure core, whereas β spares it**: epi↔epi α cophenet obs +0.423 vs pair-count null +0.129, p_pair<0.001 → CONCENTRATED (genuine localization, not a pair-count fluctuation); matched-strength epi↔epi persist p=0.003. CROSS epi↔non-epi interface also concentrated; NONEPI healthy cortex carries non-strength trace. Pat_13 (largest epi burden, 30/119 contacts) inverts under epi-X (+0.028 → −0.171). WITHDRAWN: the SUBGRAPH "+0.105→+0.187 strengthening" is generic node-count, NOT epi-specific (audit_85 --stratify epi, cohort p_dec=0.135) — it was the original rationale for the 2026-06-05 elevation and is the reason for the 2026-06-18 demotion. α verdict tag UNCHANGED — strong trace, only D_coph, earned at C3 on the full graph; reported secondary characterization. PI-reviewed 2026-06-18. Source data/audit/epi_stratified/{cophenetic_cohort.csv, pairclass_decimation_cohort.csv, decimation_control_cohort.csv})
  rho_split_coph_wm_excluded: SECONDARY_survives_wm_exclusion_strengthening_NOT_wm_specific (C6, audit_83 2026-06-08 / amended audit_85 2026-06-12 → data/audit/wm_stratified/; the α cophenet trace SURVIVES the gray-only montage — obs_median +0.105→+0.183, exclude_wm matched-strength Wilcoxon p=0.014, LO-Pat_15 p=0.027. BUT the strengthening is NOT WM-specific: the random-node-decimation control (audit_85) reaches the same cohort ρ_split by removing any ~40 % of nodes (cohort p_dec=0.305, generic_nodecount); WM is neither carrying nor diluting the α cophenet trace. The "gray-resident-because-it-strengthens" framing is WITHDRAWN; load-bearing C6 claim = survival/robustness. SECONDARY/mechanistic per CONTROLS §C6, NOT a gate; α verdict tag UNCHANGED)
  grassmann_wm_excluded: SECONDARY_remains_null_on_C3_gate_under_wm_exclusion (C6, audit_86 2026-06-12; on the locked C3 cluster-extent mass gate α is no_trace under WM exclusion (cluster_p_mass^wmX=0.488) exactly as on the full graph (C3 no_trace) — α carries no subspace trace with or without WM. The α trace is cophenet-only. Per-k count 4→3 retained as co-statistic)
  anatomy_cophenet: RETRACTED_no_FDR_surviving_localization (2026-05-30 signed sampling-aware audit + 2026-06-10 atlas, locked in ANATOMY_LEDGER: the 11-DK-region "strong localized" list is retracted [0/11 cohort-supported — all regions n≤3, cohort statistic diffuse perm-p≈0.16–0.26], and the 2026-06-10 system-scale re-analysis finds α has NO FDR-surviving concentration [pars opercularis raw p=0.005 but fails BH over 53 DK regions q=0.07–0.26; no a-priori system clears]. β is the only band that localizes [→OFC system, cophenet]. The 2026-06-05 "frontal operculum" line is withdrawn. Historical DK list retained in §3.1/§C5 as record. Source data/audit/{anatomy_localization_wilcoxon,localization_atlas}/README.md)
  anatomy_cophenet_epi_excluded: MOOT_underlying_localization_retracted (the "9 of 11 DK regions reproduce under epi-X" claim is moot — the α cophenet DK-region localization itself does not survive the 2026-05-30 signed sampling-aware audit [0/11 cohort-supported], so reproducibility-under-epi-X of a retracted list carries no claim. ANATOMY_LEDGER locked. Historical detail in §C5 retained as record.)
sources:
  - .agents/preprint/locked/CONTROLS.md (locked 5-control battery)
  - .agents/preprint/locked/VERDICT_LEDGER.md (locked verdicts)
  - .agents/preprint/methods/methods_revision_2026-05-18_cophenet.md (binding methods directive)
  - data/audit/ctm_triangle/cohort_summary.csv (α row: rho_split_median, n_above_drift, rho_xprobe, Wilcoxon)
  - data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv (α row: paired Wilcoxon vs surrogate)
  - data/audit/alpha_epi_exclusion/cohort_summary.csv (α cophenet under C5 epi-X — strengthens)
  - data/audit/alpha_epi_exclusion/comparison.csv (per-patient Δρ full vs epi-X)
  - data/audit/grassmann_cluster_extent/cohort_summary.csv (α cluster-extent p=0.159, audit_70)
  - data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv (α per-k T_G + surr stats)
  - data/audit/grassmann_regate_no_filter/contig_summary.csv (α audit_66 + audit_67 contiguous windows)
  - data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv (α substrate row)
  - data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv (α raw D sensitivity)
revision_history:
  - 2026-05-19: initial brief produced from VERDICT_LEDGER.md lockdown
  - 2026-05-22: corrected anatomy_cophenet_epi_excluded entry + §5 epi-X anatomy section after patching the silent no-op bug in audit_71 --epi-x (np.isin int×str mismatch at lines 278–283 / 354–359); real epi-X result is 9 of 11 regions reproduce, not 11 of 11
  - 2026-05-28: stripped the C5 epi-X "5/10 → 7/10 rescue" narrative throughout per user directive; α verdict tag preserved at "strong trace, only D_coph" because C3 paired Wilcoxon p=0.00195 already passes at full cohort under the locked Wilcoxon-as-gate rule; C5 epi-X reframed as secondary mechanistic observation (not verdict-promoter) per amended Decision 10 + retracted Decision 1
  - 2026-06-05: **epi-exclusion ELEVATED from secondary mechanistic observation to PRIMARY INTERPRETIVE LENS** (PI directive 2026-06-05; pre-empts the "isn't this pathological-network reorganization?" reviewer attack). The α trace strengthens under epi-exclusion (+0.105 → +0.187, ratio 8.3× → 27.7×, audit_77 matched-strength p=0.0137, audit_68 one-sample p=0.0098, LO-Pat_13 p=0.0039). Added the **α EXCEPTION** from audit_77 pair-class decomposition: α uniquely *recruits* epileptic tissue (epi↔epi +0.404, p=0.0059, 6/9 patients individually significant) whereas β spares the diseased core (epi↔epi n.s. p=0.150); Pat_13 inverts under epi-X. Updated Head (2nd paragraph rewritten with α exception), frontmatter `verdict_layers.rho_split_coph_epi_excluded`, §3.1 epi-X subsection (header + status-elevation preamble + pair-class table + α-exception + brutal-honesty caveats), and §9 Interpretation. Cites audit_77 (`data/audit/epi_stratified/`). **α verdict tag UNCHANGED** (`strong trace, only D_coph`, earned at C3 on the full graph); emphasis/interpretation elevation only. Caveats kept: α epi↔epi rests on ~45 pairs (small); pair-class null is global-rewiring not within-class; epi_only underpowered. See VERDICT_LEDGER.md 2026-06-05 revision.
  - 2026-06-08: added **C6 white-matter-exclusion** sensitivity (frontmatter `verdict_layers.rho_split_coph_wm_excluded` + `grassmann_wm_excluded`) from audit_83/84 (`data/audit/wm_stratified/`). WM (dominant-`Wm` contacts, 30–57 % of montage) removal PERSISTS+STRENGTHENS the α cophenet trace (obs +0.105→+0.183, exclude_wm matched-strength p=0.014, LO-Pat_15 0.027); α Grassmann stays null (per-k cells 4→3). The α trace is gray-matter-resident, NOT a white-matter recording artifact. SECONDARY/mechanistic per CONTROLS §C6, NOT a gate and NOT a primary interpretive lens (node-count caveat); **α verdict tag UNCHANGED** (`strong trace, only D_coph`). See VERDICT_LEDGER.md 2026-06-08 revision + `responses/2026-06-08_wm-exclusion-cascade.md`.
  - 2026-06-12: **resolved the two C6 honesty flags** (audit_85 + audit_86). Flag A (node-count): the α cophenet "strengthening" under exclude_wm is NOT WM-specific — the random-node-decimation control reaches the same cohort ρ_split by removing any ~40 % of nodes (cohort p_dec=0.305, generic_nodecount). WITHDRAW the "gray-resident-because-it-strengthens" framing; the α C6 claim is **survival/robustness** (the trace still clears matched-strength on the gray-only graph; WM is neither carrying nor diluting it). Flag B (Grassmann gate): on the locked C3 cluster-extent gate α is no_trace under WM exclusion (audit_86 `cluster_p_mass^wmX=0.488`) exactly as on the full graph — α is cophenet-only either way. Updated both `verdict_layers.*_wm_excluded` entries. **α verdict tag UNCHANGED.** See VERDICT_LEDGER.md 2026-06-12 amendment.
  - 2026-06-12 (same control, applied to C5 epi-X — the PRIMARY lens): `audit_85 --stratify epi` shows the α cophenet **subgraph** "strengthens under epi-X" (+0.105→+0.187) is ALSO generic node-count (cohort p_dec=0.135). The "strengthens under epi-X" rationale for the 2026-06-05 PRIMARY-lens elevation is WITHDRAWN; the lens RE-ANCHORS on the **pair-class** decomposition (incl. the α-EXCEPTION epi↔epi recruitment) — same-graph, no node-count confound. Updated `verdict_layers.rho_split_coph_epi_excluded`. **α verdict tag UNCHANGED** (earned at C3 full graph). ⚠ Flagged for PI review (touches the 2026-06-05 directive). Source `data/audit/epi_stratified/decimation_control_cohort.csv`.
  - 2026-06-18: **PI decision — epi-exclusion DEMOTED from primary interpretive lens to SECONDARY** (mechanistic / sensitivity layer, parallel to C6). The 2026-06-05 elevation rested partly on the subgraph "strengthens under epi-X" rationale, retracted 2026-06-12 as generic node-count (audit_85 --stratify epi, p_dec=0.135); re-elevating on a post-hoc finding from the same retracted thread was rejected. The genuine pair-class dissociation is RETAINED as a reported secondary tissue-characterization — **α RECRUITS the diseased core** (epi↔epi concentrated, cophenet obs +0.423 vs null +0.129, p_pair<0.001; matched-strength persist p=0.003), the α-specific contrast to β's spared core. Updated to n=10 (audit_77 regenerated; pair-count control audit_85 --mode pairclass). Updated frontmatter `rho_split_coph_epi_excluded`, Head epi paragraph, §3.1 status header + §9 interpretation. **α verdict tag UNCHANGED** (earned at C3 full graph). See VERDICT_LEDGER.md 2026-06-18 amendment; handoff `.agents/reports/2026-06-12_decimation-controls-handoff.md`.
---

# α band (8–13 Hz) — preprint result report

## Head

α carries a **per-pair multiscale trace** via the cophenetic communication distance `ρ_split^coph` on `D_coph = cophenet(UPGMA(D(τ_max)))`. At the full n=10 cohort, all four primary controls pass: split-baseline C1 p=0.010, drift-floor C2 p=0.007, matched-strength surrogate C3 paired Wilcoxon **p=0.002** (effect-size ratio 8.3×), cross-probe C4 paired-Wilcoxon non-degradation. The α verdict is anchored at C3 alone under the locked Wilcoxon-as-gate rule (per-patient `n_above_surrogate` = 5/10 is a descriptive cohort-agreement statistic, NOT a verdict modifier; patient-count thresholds were retired 2026-05-19). The Grassmann subspace probe finds **no significant contiguous window** (cluster-extent permutation p = 0.159, audit_70; 4-cell observed run at k = 11..14 within the null distribution, null 95th percentile = 6.0). Verdict from `locked/VERDICT_LEDGER.md`: **`strong trace, only D_coph`**.

**Epi-stratification is a secondary tissue-characterization (sensitivity layer, parallel to C6; PI-reviewed 2026-06-18 — the 2026-06-05 "primary interpretive lens" elevation is retracted), with an α-specific twist.** All statements here are in the LRG cophenetic framework (raw |ImCoh| is the comparison baseline only). Like β, the α trace is largely **healthy network reorganization that the epileptic zone partly masks** — but **α has an exception: it uniquely RECRUITS the diseased seizure core, whereas β spares it.** On the same-graph pair classes (no node-count confound), re-run at n=10 and put through a pair-count control (`audit_85 --mode pairclass`), the epi↔epi class at α is a genuine hotspot: cophenet obs **+0.423** vs pair-count null +0.129, **p_pair < 0.001 → concentrated** (matched-strength epi↔epi persist p=0.003), and the epi↔non-epi interface (cross) is likewise concentrated — against a β epi↔epi that trends *depleted*. Healthy↔healthy cortex (nonepi) carries real non-strength trace. Pat_13 — largest epi-zone burden (N_epi = 30/119) — inverts under epi-X (+0.028 → −0.171). **Withdrawn (brutal honesty):** the earlier "the trace *strengthens* +0.105 → +0.187 under epi-exclusion" framing is a generic node-count effect (`audit_85 --stratify epi`, p_dec = 0.135), not epi-specific — it was the original rationale for the 2026-06-05 elevation and the reason this is now demoted to secondary. **The α verdict tag (`strong trace, only D_coph`) is unchanged**, earned at C3 on the full graph (see §3.1).

## Headline three-layer cohort table (α-band row, copied from `locked/VERDICT_LEDGER.md`)

```
Layer                                       α (full n=10)
-------------------------------------      ---------------------
Raw FC          ρ_split^raw  (audit_67)    p=.014, 8/10  (strongest substrate of any band)
Raw D(τ_max)    ρ_split      (preprint_05) p=.042, 7/10  (ratio 11.4×)
Cophenet D_coph ρ_split^coph (audit_63)    p=.002, 5/10  (ratio 8.3×)   ← verdict gate (Wilcoxon)
Cophenet D_coph ρ_split^coph (audit_68, C5 epi-X)   p=.0098, ratio 27.7× ← secondary mechanistic observation
```

Wilcoxon p-values are the verdict gate at each layer; the `N/10` per-patient counts are descriptive auxiliary statistics. Sources: `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv`, `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv`, `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv`, `data/audit/alpha_epi_exclusion/cohort_summary.csv`. n=10 full / n=10 epi-X (per-patient N_epi varies; see §6).

α has the **strongest substrate signal of any band** at the baseline layer (raw FC p=.014 — but raw FC is the comparison baseline only, never a result). The result-bearing cophenetic step preserves a multiscale α trace that passes C3 at p=0.002 with a 8.3× effect-size ratio. (The earlier "under epi-zone exclusion the ratio strengthens to 27.7×" reading is withdrawn as generic node-count, audit_85 --stratify epi p_dec=0.135; the genuine epi-relevant finding is the pair-class α-recruits-core dissociation, §3.1.)

## 1. Scientific claim

**Cohort-level question.** Does the α-band (8–13 Hz) post-task resting state at the LRG-CTM layer sit closer to the task state than the pre-task resting state does, in the cophenetic communication geometry of `|ImCoh|`-derived FC matrices?

**Refined biological claim.** α-band cortical coupling reorganizes during the task and the reorganization is retained in the per-pair cophenetic multiscale geometry but **not** in the leading-mode subspace structure. The reorganization is per-pair multiscale (cophenet merges integrate `N−1` dendrogram scales) and **single-aspect** (Grassmann `d_G(k)` is silent — the slowest non-trivial subspace modes do not rotate at α). Mechanistic adjunct (secondary, not a verdict-bearing claim): the α trace is carried mainly by non-epileptogenic cortex, but α uniquely also *recruits* the diseased core (epi↔epi concentrated, p_pair < 0.001; §3.1). The earlier "epi-exclusion strengthens the effect ~3.4×" framing is withdrawn as generic node-count (audit_85).

**Falsification budget.** The claim fails if any of the primary controls fail: (a) `ρ_split^coph` cohort-paired Wilcoxon p ≥ 0.05 at C3 (matched-strength); (b) cross-probe C4 paired Wilcoxon rejects (`split > xprobe` significant degradation) or sign mismatch; (c) LOO max-p at C3 exceeds 0.05 (Decision 11 robustness). C5 epi-X is a secondary mechanistic observation per amended Decision 10 — it is not a primary falsifier of the α verdict.

## 2. Cohort, substrate, and library entry points

Same as β: n=10 patients (Pat_02/03/05/06/07/08/10/13/14/15); `imcoh_abs` band-averaged magnitude of imaginary coherence (Ewald 2012; Bastos & Schoffelen 2016); β-band `B = [8, 13] Hz`. Pat_03 acquired at 1024 Hz; handled at the config layer (`nperseg_for_fs(fs)`, `FS_OVERRIDES`) — full cohort member at analysis layer (policy updated 2026-05-18, no dropout sensitivity test). Pat_10 task rows [53, 54, 55] dropped at load. Pat_14 vendor-replaced 2026-04-25.

```python
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
W = load_fc_matrix(patient="Pat_02", phase="rest_pre", band="alpha", fc_method="imcoh_abs")
lrg = load_lrg_result(patient="Pat_02", phase="rest_pre", band="alpha", fc_method="imcoh_abs")
# lrg.ultrametric_matrix is D_coph (cophenet of UPGMA(D(τ_max))) — canonical per-pair object
# lrg.eigenvalues, lrg.eigenvectors are the Laplacian spectrum for the Grassmann probe
```

## 3. The two LRG probes — methodology, results, and provenance

### 3.1 Per-pair multiscale correlation on `D_coph` — **PRIMARY FINDING for α**

#### What `D_coph` is and why it is the per-pair multiscale object
`D_coph_ij = h(LCA_L(i, j))` is the merge height in the UPGMA dendrogram of `D(τ_max) = 1/ρ̂(τ_max)` at which contacts `i` and `j` first coalesce. The cophenetic image of a hierarchical clustering is ultrametric by construction, and each entry reads the intrinsic communication-coalescence scale of that pair. Because our `|ImCoh|` Laplacian spectrum is continuous (no gap-induced crossover scales), the dendrogram's `N−1` merge heights replace the spectral-gap scale identification of [Villegas 2025]; `D_coph` is therefore the natural multiscale per-pair carrier in this regime. See `methods/methods_revision_2026-05-18_cophenet.md` for the full justification.

#### Critical preamble (5-point)
1. **Claim.** α `ρ_split^coph(taskT, rsPost, rsPre_A, rsPre_B) > 0` at cohort level.
2. **Null (C3, mandatory).** Strength-preserving 4-cycle ±δ rewiring of each per-phase adjacency, R=200, seed 20260511. Paired Wilcoxon one-sided greater across patients: obs > surrogate.
3. **Strongest plausible alternative.** Same-band node-strength evolution from rsPre to rsPost driving the per-pair distance shift without genuine multiscale reorganization.
4. **Does C3 cover (3).** Yes — matched-strength surrogate preserves `s_i` to 10⁻⁴ and reshuffles only the topology. What C3 cannot exclude: physiological drift coupled to task block (would need a no-task control session, not present in the paradigm) — C2 drift-floor partially addresses this. What C3 cannot exclude alone: epi-zone-driven cohort dilution (addressed by C5 epi-X below, audit_68).
5. **Falsification.** C3 paired Wilcoxon p ≥ 0.05, or C2 drift floor not separable, or C4 cross-probe flips sign, or C5 epi-X weakens rather than strengthens.

#### Method
- `Δ_task(i,j) = D_coph^taskT(i,j) − D_coph^rsPre_A(i,j)`
- `Δ_rest(i,j) = D_coph^rsPost(i,j) − D_coph^rsPre_B(i,j)`
- `ρ_split^coph = ρ_S(Δ_task, Δ_rest)` Spearman over all `N(N−1)/2` pair entries.
- Independent rsPre halves A/B; full LRG pipeline (eigendecomp → propagator at `τ_max = 1/λ_max` → UPGMA → cophenet) rerun per half.

#### Results (α, full cohort)

| Statistic | Value | Source |
|---|---|---|
| Cohort-median `ρ_split^coph` | +0.1146 | `ctm_triangle/cohort_summary.csv` alpha row |
| n_trace_split | 8/10 | same |
| n_above_drift | 8/10 | same |
| C1 Wilcoxon p (`ρ_split > 0`) | **0.00977** | same col `wilcoxon_split_gt_0_p` |
| C2 Wilcoxon p (`ρ_split > ρ_drift`) | **0.00684** | same col `wilcoxon_split_gt_drift_p` |
| C3 matched-strength: obs_median | +0.1054 | `matched_strength_surrogate_split_baseline/cohort_summary.csv` alpha |
| C3 surr median (median over surrogate medians) | +0.0128 | same |
| C3 ratio obs/surr | **8.25×** | derived |
| C3 paired Wilcoxon z / p | 54.0 / **0.00195** | same |
| C3 n_above_surrogate (per-patient) | 5/10 | same — descriptive auxiliary statistic only |
| C4 cross-probe rho_xprobe | +0.105 | `ctm_triangle/cohort_summary.csv` alpha |
| C4 n_trace_xprobe | 8/10 (+sign matches +0.115) | same |

All 4 primary controls (C1/C2/C3/C4 Wilcoxons) pass at p < 0.05. Per CLAUDE.md `feedback_no_hardcoded_test_thresholds.md` + `feedback_patient_counts_never_the_gate.md`, the Wilcoxon test IS the gate; per-patient counts are descriptive auxiliary statistics and never modify the verdict.

#### Epileptogenic-zone exclusion — PRIMARY INTERPRETIVE LENS (audit_68 2026-05-20 + audit_77 epi-stratified 2026-06-05)

**Status — secondary tissue-characterization (PI-reviewed 2026-06-18).** Epi-zone exclusion is a **secondary mechanistic observation / sensitivity layer (parallel to C6)**, **not** a primary interpretive lens — the 2026-06-05 elevation to "primary lens" rested partly on the now-retracted subgraph "strengthens under epi-X" rationale (generic node-count, `audit_85 --stratify epi`, p_dec = 0.135) and is **demoted to secondary**, symmetric with C6. The message it still carries — a useful reviewer-defense — is that the α cognitive task-trace is **healthy network reorganization that the epileptic zone partly masks**, not an epileptiform artefact, **with the α-specific twist that α also recruits the diseased core** (epi↔epi concentrated, p_pair < 0.001; see below). **This is reported emphasis only; the α verdict tag (`strong trace, only D_coph`) is unchanged** and was earned at C3 on the full graph (the Wilcoxon-as-gate rule; C1–C4 all pass at full n=10 in §3.1). A reader must not read this section as a change of verdict.

| Statistic | Full cohort | Epi-X cohort | Δ |
|---|---|---|---|
| Cohort-median `obs_ρ` | +0.105 | **+0.187** | +0.082 |
| Surrogate median | +0.0128 | +0.0067 | −0.006 |
| Ratio obs/surr | 8.25× | **27.70×** | +3.4× |
| audit_68 one-sample Wilcoxon p (`ρ_split^epi-X > 0`) | — | **0.0098** | passes (LOO max p = 0.0195 Pat_03) |
| audit_77 matched-strength Wilcoxon p | — | **0.0137** | passes (LO-Pat_13 p = 0.0039) |
| n_above_surrogate | 5/10 | **7/10** | +2 patients |

Under epi-X the α cophenet trace remains matched-strength-significant on both the locked audit_68 one-sample Wilcoxon (p=0.0098) and the audit_77 re-run (p=0.0137, LO-Pat_13 p=0.0039) — i.e. it *survives* removing the epileptic zone. (The earlier "the ratio strengthens ~3.4× under epi-X" emphasis is withdrawn as generic node-count, audit_85 --stratify epi p_dec=0.135 — not epi-specific.) The α trace is carried by non-epileptogenic cortex *plus* epileptic tissue (see the pair-class decomposition + α exception below). Sources: `alpha_epi_exclusion/cohort_summary.csv` (audit_68), `data/audit/epi_stratified/cophenetic_cohort.csv` (audit_77).

#### Per-patient C5 epi-X comparison (`alpha_epi_exclusion/comparison.csv`)

| Patient | `obs_ρ` full | `obs_ρ` epi-X | Δ | N_epi excluded | reading |
|---|---|---|---|---|---|
| Pat_02 | −0.037 | −0.035 | +0.002 | 14 | no trace either; tiny change |
| Pat_03 | +0.360 | +0.189 | −0.171 | 6 | strong both; weakens but still strong |
| Pat_05 | +0.078 | +0.121 | +0.043 | 14 | **strengthens to + trace** |
| Pat_06 | **+0.833** | +0.692 | −0.141 | 10 | dominant + trace both |
| Pat_07 | +0.081 | +0.211 | +0.130 | 7 | **strengthens to + trace** |
| Pat_08 | +0.380 | +0.290 | −0.090 | 9 | + trace both; weakens |
| Pat_10 | +0.130 | +0.193 | +0.063 | 10 | **strengthens to + trace** |
| Pat_13 | +0.028 | **−0.171** | −0.199 | **30** | **inverts** — patient has largest epi-zone fraction |
| Pat_14 | +0.210 | +0.184 | −0.026 | 12 | + trace both; mild weakening |
| Pat_15 | +0.040 | +0.040 | 0 | **0** | no epi nodes — invariant |

Read: under epi-X, three patients move into the + trace direction (Pat_05, Pat_07, Pat_10), one patient inverts (Pat_13 — the only patient with N_epi ≥ 30, ~25% of contacts), the four already strongly positive (Pat_03, Pat_06, Pat_08, Pat_14) stay strongly positive with mild weakening, and Pat_15 (zero epi nodes) is invariant by construction. The cohort signal becomes more uniform and stronger under epi exclusion.

#### Pair-class decomposition — α recruits epileptic tissue while β spares it (audit_77, 2026-06-05)

The full-graph `ρ_split^coph` null is re-evaluated restricted to three node-pair classes (holding the global per-node strength sequence fixed; median pair counts at n=10: nonepi↔nonepi 5408, cross 1167, epi↔epi 55):

| Pair class | α cohort-median `ρ_split^coph` (p) | β analogue (p) | Reading |
|---|---|---|---|
| NONEPI↔NONEPI (healthy cortex) | +0.073 (0.0420) | +0.197 (0.0049) | both carry |
| CROSS (epi↔non-epi interface) | +0.167 (0.0020) | +0.252 (0.0098) | both carry |
| EPI↔EPI (diseased core, pair-restricted) | **+0.423 (0.0029)** | +0.142 (**0.161**) | **α carries, β does NOT** |

Source: `data/audit/epi_stratified/cophenetic_cohort.csv` (audit_77, n=10). **The α EXCEPTION is the story.** Where β spares the diseased core (epi↔epi non-significant, p=0.161), **α uniquely also recruits the epileptic tissue into the task-trace**: the α epi↔epi pair class is the *strongest single epi-stratified α cell* (+0.423, matched-strength p=0.0029), and **7 of 10 patients are individually significant within their own epi subgraph** (Pat_03/05/06/07/13/14/15 at own-surrogate p<0.05; obs>0 in 9/10 — only Pat_08 negative). The best-powered patient (Pat_13, 435 epi↔epi pairs) and both hub-patients (Pat_10, Pat_15) are positive, and the effect survives leave-Pat_13-out (p=0.0059). This dissociation — **α recruits epileptic tissue, β spares it** — is the band-specific epilepsy signature of the manuscript.

**Brutal-honesty caveats (bounding the α exception).** (i) The **α epi↔epi cell rests on ~55 pairs** (median) — a small per-class sample; the +0.423 is real but the Spearman is computed over few pairs and is the most variance-prone class. The dedicated **pair-count null** (`audit_85 --mode pairclass`, random same-size pair subsets) targets exactly this concern, and the cell survives it: observed +0.423 vs pair-count null median +0.129, **p_pair < 0.001, verdict *concentrated*** — so the small per-class sample is not the source of the effect. (ii) The **pair-class matched-strength null is global-rewiring, not within-class**: it fixes the global strength sequence and asks whether class-restricted co-movement survives, but does not randomize *within* the class — so the strength-axis control on "α recruits the epi core" is **suggestive, not airtight**; a within-class matched-strength null is the airtight version and has not been run. (iii) The **epi_only rebuilt subgraph is underpowered** (6–30 nodes, Pat_15 = 0, near-degenerate matched-strength) and *weakens* for α (−0.034, p=0.752) — read as "no power", not "absent"; the well-powered evidence for the α exception is the pair-class restriction of the full-graph null (+0.423), not the tiny epi_only rebuild.

#### Reading
α `ρ_split^coph` passes all 4 primary controls at full cohort, with C3 paired Wilcoxon p = 0.00195 the verdict gate. The per-patient C3 `n_above_surrogate` = 5/10 is a descriptive cohort-agreement statistic and does not modify the verdict (patient-count thresholds were retired 2026-05-19). The pair-class decomposition shows the α trace carried by healthy cortex, the epi↔non-epi interface, *and* — uniquely among the trace bands — the epileptic core itself (epi↔epi cophenet concentrated, p_pair < 0.001; §3.1). This α-recruits-the-core dissociation is a **secondary tissue-characterization** (PI-reviewed 2026-06-18 — the 2026-06-05 "primary interpretive lens" elevation is **demoted**, since its other leg, the subgraph "strengthens ~3.4× under epi-X", was withdrawn 2026-06-12 as generic node-count, audit_85 --stratify epi p_dec=0.135); **it does not change the verdict, which was earned at C3 on the full graph.**

#### Cache + script provenance
- Full cohort: `audit_63_split_baseline_surrogate.py` → `data/audit/matched_strength_surrogate_split_baseline/`
- Epi-X (one-sample Wilcoxon): `audit_68_alpha_epi_exclusion.py` → `data/audit/alpha_epi_exclusion/`
- Epi-stratified pair-class (audit_77, 2026-06-05): `scripts/01_compute/audit/audit_77_epi_stratified_cophenetic.py` → `data/audit/epi_stratified/{cophenetic_cohort.csv, cophenetic_per_patient.csv, README.md}`
- Consolidated report: `.agents/reports/2026-06-05_epilepsy-headline-and-occult-node-marker.md`
- Within-baseline triangle: `audit_*_ctm_triangle.py` → `data/audit/ctm_triangle/`
- Object built from `lrg.ultrametric_matrix` field in `data/cache/imcoh_lrg/Pat_NN/alpha_{phase}_lrg_imcoh-abs.npz`

### 3.2 Global subspace rotation: Grassmann chordal distance `d_G(k)` — **no trace at α**

#### What `d_G(k)` reads vs what `ρ_split^coph` reads
`d_G(k)` is the chordal Grassmann distance between leading-`k`-dimensional Laplacian eigenmode subspaces of two phases (`d_G² = k − ‖A^T B‖_F²` with A, B = N×k orthonormal eigenvector matrices, skipping the trivial mode). It reads **whether the slow-diffusion subspace rotates** between phases, at a chosen `k`. The per-pair `ρ_split^coph` and the Grassmann `d_G(k)` are not derivable from each other and read structurally distinct aspects of the same propagator.

#### Critical preamble (5-point)
1. **Claim.** `T_G(k) = d_G(taskT, rsPost; k) − d_G(rsPre_A, taskT; k) < 0` at cohort level for some contiguous window of `k`.
2. **Null (C3, audit_70 cluster-extent permutation).** R=200 matched-strength surrogates; for each phantom-surrogate r, cohort paired Wilcoxon between surrogate r and the mean of the remaining R−1 surrogates per k, longest contiguous run with p < 0.05. Empirical null distribution of longest-run length.
3. **Strongest plausible alternative.** Same-strength edge rewiring producing a similar contiguous-significant cluster by chance.
4. **Does C3 cover (3).** Yes — matched-strength preserves node strengths, only topology shuffles; phantom-surrogate-vs-rest builds the exchangeable cluster-extent null. Replaces the previous 8-cell hardcoded threshold (VERDICT_LEDGER.md Decision 6).
5. **Falsification.** Observed longest run within the null distribution → no trace.

#### Method
- For each `k` ∈ [2, 112] and each patient: `T_G(k) = d_G(U_k^taskT, U_k^rsPost) − d_G(U_k^rsPre_A, U_k^taskT)`.
- Cohort paired Wilcoxon one-sided less per `k`: obs `T_G(k)` < mean-surrogate `T_G(k)` across patients.
- Cluster-extent: longest contiguous run of k with cohort-Wilcoxon p < 0.05.
- Null: for each surrogate r in R=200, phantom-test r vs (R−1)-mean, longest run. Cluster p = (1 + #(null ≥ obs)) / (R + 1).

#### Results (α)

| Statistic | Value | Source |
|---|---|---|
| Observed longest contiguous-significant run | **4 cells** at k = 11..14 | `grassmann_regate_no_filter/contig_summary.csv` audit_66 alpha |
| Cluster-extent null mean LR | 2.09 | `grassmann_cluster_extent/cohort_summary.csv` alpha |
| Cluster-extent null 95th LR | 6.0 | same |
| Cluster-extent null max LR | 12 | same |
| Observed cluster mass `T_G^*` (raw Σ −log10p) | 7.98 | same |
| Observed `T_G^*` (normalized per C1) | **0.0312** (raw / 255.65, n_k=111, R=200) | same |
| Null mean cluster mass | 7.35 (post-fix all-clusters) | same |
| Null 95th cluster mass | 21.26 (post-fix) | same |
| **cluster_p_longest_run (descriptive)** | 0.1592 | same |
| **cluster_p_mass (Decision-8 gate)** | **0.3482** | same — no trace (mass-only gate fails) |
| C5 epi-X (audit_67): obs longest run | 3 cells | `grassmann_regate_no_filter/contig_summary.csv` audit_67 alpha |

Observed 4-cell longest run is well within the null distribution (null 95th = 6.0; null mean = 2.09). Both cluster-p statistics fail the 0.05 gate; the cluster-mass statistic is right on the boundary (0.0995) but still above. Under C5 epi-X the longest run further shrinks to 3 cells.

#### Reading
α has no detectable Grassmann subspace trace under the cluster-extent permutation gate. The 4-cell observation at k=11..14 is unremarkable against the empirical null. The α reorganization does not live in the slow-mode subspace; it lives only in the per-pair cophenetic geometry.

#### Cache + script provenance
- `audit_66_grassmann_matched_strength_surrogate.py` → `data/audit/grassmann_matched_strength_surrogate/`
- `audit_67_grassmann_epi_exclusion.py` → `data/audit/grassmann_epi_exclusion/`
- `audit_69_grassmann_regating_no_filter.py` → `data/audit/grassmann_regate_no_filter/`
- `audit_70_grassmann_cluster_extent.py` → `data/audit/grassmann_cluster_extent/`
- Object built from `eigvecs` field in `data/cache/imcoh_lrg/Pat_NN/alpha_{phase}_lrg_imcoh-abs.npz` (obs) and `data/cache/matched_strength_surrogate_lrg/Pat_NN/alpha_{phase}_R200_swap20_seed20260511_imcoh_abs.npz` (surrogate)

## 4. Multiscale structure of the α trace

### Per-pair multiscale via `D_coph` (primary)

The cophenet integrates across `N−1 ≈ 112–121` merge heights per patient. The α `ρ_split^coph` is therefore an integrated multiscale signal at the per-pair level: it cannot be reduced to a single-scale per-pair comparison because the cophenetic image carries each pair to its own coalescence scale.

The headline three-layer table shows that raw FC (single-pair magnitudes) and raw `D(τ_max)` (single-scale propagator distance) — both comparison baselines only — and the result-bearing cophenet `D_coph` all detect the α trace at the Wilcoxon level (p = 0.014 / 0.042 / 0.002 respectively); the descriptive per-patient counts shift from 8/10 → 7/10 → 5/10 across the three layers. (The earlier "under C5 epi-X the cophenet ratio strengthens ~3.4×, selectively retaining the non-epi-cortex component" reading is withdrawn as generic node-count, audit_85 --stratify epi p_dec=0.135.) The C3 Wilcoxon is the verdict gate at each layer.

### Subspace multiscale via Grassmann `d_G(k)` (silent)

α has no contiguous-significant Grassmann window. The 4-cell run at k = 11..14 is within the empirical null (cluster_p = 0.159). Reading: α's reorganization does not project into a coherent rotation of the leading-mode subspace at any `k`. This is what distinguishes α from β: β is multi-probe (both `ρ_split^coph` and Grassmann pass C3 with epi-X strengthening); α is single-probe (`ρ_split^coph` only).

## 5. Anatomical distribution — ~~strong-localized to a bilateral cingulate + parahippocampal + parietal network~~ **RETRACTED 2026-05-30 / DELOCALIZED 2026-06-01 (diffuse)**

> ⚠️ **RETRACTED — the α cophenet trace is anatomically DIFFUSE / delocalized.** The
> "11 named DK regions" localization does **not** survive a signed, threshold-free test:
> **0/11 cohort-supported** (verdict DIFFUSE; all 11 locked regions are sampled by ≤3
> patients — four by 1, five by 2, two by 3 — and no well-sampled region passes), and
> **per-patient localization is also null** (1/10 patients, = chance: the lone passer is
> Pat_15, p=0.029, one of ~5 expected over 100 patient-tests). The 11-region list was
> built on top-decile `|Δρ|` + A1/A3 enrichment on endpoint counts (direction-blind,
> count-based, partly driven by Wm/Unk junk labels that this test excludes). The "9 of 11
> reproduce under C5 epi-X" question below is moot — none localize cohort-wide. The
> verified α trace (matched-strength, cophenet) is real but **spatially distributed with
> no anatomical anchor at cohort or single-patient level**. Tables retained as historical
> record. Sources: `data/audit/anatomy_localization_wilcoxon/README.md`,
> `data/audit/per_patient_localization/README.md`, `locked/ANATOMY_LEDGER.md`.

~~The α cophenet trace is **localized to a distributed cortical network** of 11 named DK regions that all pass A1+A3 ... a non-epileptic-cortex multiscale phenomenon.~~ (retracted — see banner) The original "11/11 reproduce identically" claim recorded here in the 2026-05-19 lockdown was an artifact of a silent no-op bug in the audit_71 `--epi-x` mask (`np.isin(int_array, str_array)` always-False); the bug was patched and the audit rerun on 2026-05-22 — but the whole region list is now retracted regardless (diffuse cohort + null per-patient).

Audited under `locked/ANATOMY_CONTROLS.md`. **Verdict source NOW: `data/audit/anatomy_localization_wilcoxon/` + `data/audit/per_patient_localization/` (2026-05-30/06-01 retraction); `ANATOMY_LEDGER.md` 2026-05-19 entry superseded.**

### α cophenet anatomy (top-decile per-pair `|Δρ_split^coph|`)

Eleven named Desikan–Killiany regions pass A1+A3 jointly (full cohort, R=200):

| Region | k_trace / K_region | A1 enrichment | A1 p_hyper | A3 obs_z | A3 p_emp |
|---|---|---|---|---|---|
| ctx-lh-caudalanteriorcingulate | 45 / 349 | 2.06× | 4.07e-06 | 4.67 | 0.005 |
| ctx-lh-parahippocampal | 85 / 357 | 3.80× | 7.53e-27 | 9.76 | 0.005 |
| ctx-lh-rostralanteriorcingulate | 90 / 472 | 3.04× | 3.66e-21 | 14.26 | 0.005 |
| ctx-rh-medialorbitofrontal | 242 / 2188 | 1.77× | 1.22e-17 | 3.22 | 0.005 |
| ctx-rh-caudalmiddlefrontal | 103 / 567 | 2.90× | 2.49e-22 | 6.23 | 0.005 |
| ctx-rh-postcentral | 247 / 1135 | 3.47× | 6.70e-67 | 4.47 | 0.005 |
| ctx-lh-caudalmiddlefrontal | 76 / 815 | 1.49× | 4.10e-04 | 3.02 | 0.015 |
| ctx-rh-caudalanteriorcingulate | 52 / 340 | 2.44× | 2.85e-09 | 3.45 | 0.020 |
| ctx-rh-precuneus | 86 / 936 | 1.47× | 2.92e-04 | 2.60 | 0.025 |
| ctx-rh-superiorparietal | 84 / 819 | 1.64× | 7.94e-06 | 3.64 | 0.025 |
| ctx-rh-posteriorcingulate | 100 / 804 | 1.99× | 7.39e-11 | 3.13 | 0.045 |

~~Network anatomy: **bilateral cingulate** ... spanning cingulate / medial-temporal / medial-frontal / parietal cortex. Cohort verdict: **strong localized**.~~ **RETRACTED (see §5 banner): 0/11 cohort-supported (DIFFUSE); all 11 regions n≤3; per-patient also null (1/10 = chance).**

Source: `data/audit/anatomy_alpha_cophenet/cohort_summary.csv`.

### C5 epi-X anatomy: 9 of 11 named regions reproduce (post-patch 2026-05-22) — ~~current~~ **RETRACTED 2026-05-30**

> ⚠ **RETRACTED (anatomy localization, 2026-05-30 / locked in `ANATOMY_LEDGER.md`).** The underlying α cophenet DK-region localization does **not** survive a signed, sampling-aware test (0/11 regions both adequately sampled [n≥5] and cohort-localized; the cohort statistic is diffuse, perm-p≈0.16–0.26). "9 of 11 reproduce under epi-X" is therefore reproducibility of a **retracted** list and carries no claim. The 2026-06-10 system-scale re-analysis confirms **α has no FDR-surviving localization** (pars opercularis fails BH over 53 regions; no a-priori system clears) — β is the only band that localizes (→OFC). The tables below are **historical record**.

*(Historical record — the underlying α DK-region localization is retracted, see banner above.)* Re-running the α cophenet anatomy audit with epi-zone contacts removed per patient (`load_epileptic_nodes() × channel_labels.csv`) had yielded **9 of 11 named DK regions surviving** the joint A1+A3 gate. The two dropouts:

| Dropped region | Full cohort | Epi-X | Reading |
|---|---|---|---|
| `ctx-lh-caudalmiddlefrontal` | enr 1.49×, q_BH 1.6e-3, z 3.02, p_emp 0.015 | A1 FAIL: q_BH 5.1e-2 (above 0.05); A3 z 2.32, p_emp 0.020 | Weakest of the 11 in full cohort; marginal statistical fail under the slightly different cohort marginals after epi-X masking. |
| `ctx-rh-medialorbitofrontal` | enr 1.77×, q_BH 1.1e-16, z 3.22, p_emp 0.005 | A1 still strong (q_BH 3.3e-6); **A3 FAIL: z 1.28, p_emp 0.124** | A2/A3 (matched-strength surrogate) no longer separates the residual rh-medial-OFC signal from matched-strength noise after epi-X. **Biological reading**: a substantial part of the full-cohort medial-OFC enrichment was carried by topology coupled to epi-zone contacts; the matched-strength control correctly strips that contribution. |

The 9 surviving regions and how they shift under epi-X:

| Region | Full enr | Epi-X enr | Reading |
|---|---|---|---|
| `ctx-lh-caudalanteriorcingulate` | 2.06× (q 2.1e-5) | **2.08× (q 2.9e-2)** | stable |
| `ctx-lh-parahippocampal` | 3.80× (q 1.4e-25) | **3.89× (q 6.6e-25)** | medial-temporal anchor; strengthens slightly |
| `ctx-lh-rostralanteriorcingulate` | 3.04× (q 4.1e-20) | **5.37× (q 1.5e-31)** | **largest gainer under epi-X** |
| `ctx-rh-caudalanteriorcingulate` | 2.44× (q 1.8e-8) | **3.31× (q 3.2e-11)** | strengthens |
| `ctx-rh-caudalmiddlefrontal` | 2.90× (q 3.5e-21) | 2.74× (q 5.2e-17) | mild weakening |
| `ctx-rh-postcentral` | 3.47× (q 3.8e-65) | **3.66× (q 7.4e-69)** | rock stable; deepest q_BH in the cohort |
| `ctx-rh-posteriorcingulate` | 1.99× (q 5.9e-10) | 1.93× (q 1.1e-8) | mild weakening |
| `ctx-rh-precuneus` | 1.47× (q 1.3e-3) | 1.47× (q 1.1e-3) | identical |
| `ctx-rh-superiorparietal` | 1.64× (q 3.7e-5) | 1.64× (q 3.4e-5) | identical |

Net reading: the 9-region surviving network — **bilateral cingulate (4 regions: lh-/rh-caudal-AC, lh-rostral-AC, rh-posterior-cingulate) + left parahippocampal cortex + right caudal middle frontal gyrus + right postcentral / precuneus / superior parietal** — sees several constituent regions (particularly bilateral cingulate and left parahippocampal cortex) **strengthen** under epi-X. The pattern mirrors the cohort effect-size strengthening (`ρ_split^coph` cohort ratio 8.3× → 27.7× under epi-X — secondary mechanistic observation, not a verdict-promoter) and reads the α reorganization as a **physiological cortical phenomenon, with bilateral cingulate and left medial-temporal cortex as the anatomical core** and a small epi-zone-coupled contribution that explains the rh-medial-OFC dropout under matched-strength.

Source: `data/audit/anatomy_alpha_cophenet_epiX/cohort_summary.csv` (post-patch 2026-05-22; MD5 `199d742b22e00893ff3fbae4623a1b1c`).

**Bug history.** The original `anatomy_alpha_cophenet_epiX/cohort_summary.csv` produced 2026-05-19 was byte-identical to the full-cohort file because `audit_71_anatomy_cophenet.py` had a silent no-op at the `--epi-x` mask: `np.isin(np.arange(len(regions_df)), np.asarray(list(epi)))` compared integer FC-channel indices to the string labels returned by `load_epileptic_nodes()`, returning all-False, so `keep_mask` was always all-True. Patch applied 2026-05-22: replaced the comparison with the label-normalised lookup pattern used by `audit_68_alpha_epi_exclusion.py` (`_normalise_label(label_raw)` from `lrg_eegfc.utils.io.regions`). Sanity-checked against audit_68 (Pat_02 N_epi=14, Pat_13 N_epi=30, Pat_15 N_epi=0). The "11/11 reproduce identically" claim in this brief's 2026-05-19 lockdown originated from the bug — it was never a real measurement. The verified outcome above replaces it.

### Caveats

- A1 also enriches `Wm` (white matter, 1.04×, p_hyper=6.31e-04) and `unknown` (1.73×, 1.48e-08) — these are non-DK trivial bins, excluded from the localization claim per `locked/ANATOMY_CONTROLS.md`.
- A2 sampling-corrected bootstrap and A4 implant-geometry regression deferred.
- Grassmann probe at α has no trace (cluster_p=0.159; 4-cell run within the null distribution) — no Grassmann anatomy verdict is computed.

### Cache + script provenance
| Artifact | Path |
|---|---|
| α cophenet anatomy (full) | `data/audit/anatomy_alpha_cophenet/cohort_summary.csv` |
| α cophenet anatomy (C5 epi-X) | `data/audit/anatomy_alpha_cophenet_epiX/cohort_summary.csv` |
| Producing script | `scripts/01_compute/audit/audit_71_anatomy_cophenet.py --band alpha [--epi-x]` |
| Anatomy battery doc | `.agents/preprint/locked/ANATOMY_CONTROLS.md` |
| Anatomy verdict ledger | `.agents/preprint/locked/ANATOMY_LEDGER.md` (2026-05-19 entries) |

## 6. Patient-by-patient reading

### α `ρ_split^coph` per-patient signs across the primary probe

The full table is in §3.1 above (the C5 epi-X comparison). Summary at the cohort level:

- **Strong + trace at full cohort**: Pat_06 (+0.833), Pat_08 (+0.380), Pat_03 (+0.360), Pat_14 (+0.210)
- **Marginal + at full cohort, strengthens under C5**: Pat_05 (+0.078 → +0.121), Pat_07 (+0.081 → +0.211), Pat_10 (+0.130 → +0.193)
- **Borderline at full cohort, inverts under C5**: Pat_13 (+0.028 → −0.171; N_epi = 30, largest in cohort)
- **No trace either**: Pat_02 (−0.037 → −0.035), Pat_15 (+0.040 → +0.040; N_epi = 0)

7/10 patients are in the + trace direction under C5 epi-X (Pat_03/05/06/07/08/10/14). Pat_02 and Pat_15 are essentially flat. Pat_13 is the only inversion — driven by the largest epi-zone fraction in the cohort.


## 7. Robustness panel (under the locked battery)

| Check | α verdict | Source |
|---|---|---|
| C1 within-rsPre split-half null | ✓ 8/10 cohort, p = 0.010 | `ctm_triangle/cohort_summary.csv` alpha |
| C2 drift-floor `ρ_split > ρ_drift` | ✓ 8/10 above drift, p = 0.007 | same |
| C3 matched-strength surrogate `ρ_split^coph` | ✓ paired Wilcoxon p = 0.002 (verdict gate); ratio 8.3×; per-patient 5/10 descriptive | `matched_strength_surrogate_split_baseline` alpha |
| C3 matched-strength surrogate Grassmann | ✗ cluster_p = 0.159 (no trace) | `grassmann_cluster_extent` alpha |
| C4 cross-probe restriction | ✓ +0.105, paired-Wilcoxon non-degradation; descriptive 8/10 +sign | `ctm_triangle` alpha (rho_xprobe column) |
| C5 epi-X `ρ_split^coph` (secondary mechanistic) | Wilcoxon-on-epi-X p = 0.0098 (trace *survives* epi-X); "ratio strengthens 8.3×→27.7×" WITHDRAWN as node-count (audit_85, p_dec=0.135); **not verdict-promoter** | `alpha_epi_exclusion/cohort_summary.csv` |
| C5 epi-X Grassmann | (n/a; Grassmann has no trace at α) | `grassmann_regate_no_filter` audit_67 alpha |
| Pat_13 leave-one-out under C5 | Pat_13 inverts under epi-X — mechanistically consistent with the largest epi-zone fraction; descriptive only | `alpha_epi_exclusion/comparison.csv` |
| Cross-probe restriction (C4) | preserves sign | `ctm_triangle` |
| Three-layer cohort table (raw FC / raw D / cophenet) | all three layers pass C3 at Wilcoxon; cophenet C5 epi-X strengthens the ratio, supportive non-epi-cortex reading | headline section above |

## 8. Sensitivity panel under the locked 5-control battery

The α verdict is anchored at C3 paired Wilcoxon p = 0.00195 alone (locked Wilcoxon-as-gate rule; patient-count thresholds retired 2026-05-19, reaffirmed 2026-05-28). C5 epi-X is reported as a secondary mechanistic observation per amended Decision 10 + retracted Decision 1. **The verdict for α is `strong trace, only D_coph` based on C1–C4 at full cohort.**

| Control | α status | Notes |
|---|---|---|
| C1 split | ✓ p = 0.010 | passes |
| C2 drift | ✓ p = 0.007 | passes |
| C3 matched-strength (cophenet) | ✓ paired Wilcoxon p = 0.002 | verdict gate; per-patient 5/10 descriptive only |
| C3 matched-strength (Grassmann) | ✗ cluster_p = 0.159 | Grassmann probe inactive at α |
| C4 cross-probe (paired Wilcoxon) | ✓ non-degradation | passes |
| C5 epi-X (cophenet) | secondary mechanistic observation (audit_68) | trace survives epi-X (p=0.0098); "strengthens 3.4×" WITHDRAWN as node-count (audit_85, p_dec=0.135); not verdict-promoter |
| C5 epi-X (Grassmann) | (n/a; Grassmann has no trace at α) | confirms no Grassmann trace |

## 9. Interpretation

**Neuroscientific framing.** α (8–13 Hz) is the canonical posterior/attentional rhythm and the dominant inter-areal coupling band at rest in human iEEG (Spaak et al. 2012; Mathewson 2011; Klimesch 2012). A non-epileptic-cortex α trace at the LRG-CTM layer is consistent with **attentional-state carryover** from task into post-task rest, specifically in the per-pair cophenetic communication geometry rather than in the slow-subspace structure.

**Why per-pair only (not Grassmann).** β (13–30 Hz) carries both `ρ_split^coph` and Grassmann subspace traces (multi-aspect). α carries only the per-pair multiscale signature. The interpretation is that the α reorganization is **distributed across the cophenetic merge-height continuum** without producing a coherent rotation of the slow-mode subspace at any single `k`. Mechanistically this is consistent with a wide-band, per-pair coupling re-weighting that does not concentrate in the slowest few diffusion modes — distinct from β, where the reorganization rotates a 29-mode subspace.

**α recruits the epileptic core, whereas β spares it (secondary tissue-characterization, PI-reviewed 2026-06-18).** All numbers here are on the cophenetic object `D_coph` — raw |ImCoh| is the comparison baseline only, never a result. The genuine epi-relevant finding is a pair-class dissociation: on the same-graph pair classes (no node-count confound) put through the pair-count control (`audit_85 --mode pairclass`, n=10), the α epi↔epi class is a genuine hotspot — cophenet obs **+0.423** vs pair-count null +0.129, **p_pair < 0.001 → concentrated**, matched-strength epi↔epi persist p=0.003 — whereas the β epi↔epi class *trends depleted* (cophenet p_pair=0.945, n.s.; matched-strength weaken). So at α the trace lives in healthy cortex, the epi↔non-epi interface, *and* the diseased core; at β it spares the core. Pat_13 — largest epi-zone burden (N_epi=30/119) — inverts under epi-X. **Withdrawn:** the earlier "epi-zone exclusion *strengthens* α (+0.105 → +0.187, ratio 8.3× → 27.7×)" framing is a generic node-count effect (`audit_85 --stratify epi`, p_dec=0.135), not epi-specific — it was the original rationale for the 2026-06-05 "primary interpretive lens" elevation and the reason this is **demoted to a secondary tissue-characterization** 2026-06-18. The α-recruits-core dissociation is bounded by the §3.1 caveats (α epi↔epi rests on ~55 pairs, but survives the dedicated pair-count null at p_pair < 0.001). It is **not a verdict-modifier** — the α verdict was earned at C3 on the full graph.

**The role of the cophenet step at α.** Raw FC and raw `D(τ_max)` both detect α at the Wilcoxon level (baseline only). The cophenetic step — the result-bearing object — also passes C3 at p = 0.002 (verdict gate; descriptive per-patient count drops 8/10 → 7/10 → 5/10 across the three layers as the metric becomes more multiscale-aggregative). (The earlier "under epi-X the cophenet ratio strengthens 3.4×" claim is withdrawn as generic node-count, audit_85 --stratify epi p_dec=0.135; it is not evidence that cophenet acts as a spatial-purity filter.) The α cophenet trace is a genuine multiscale per-pair signature on the full graph; that is the methodology endorsement, independent of any epi-exclusion contrast.

## 10. Source-of-truth references

### Lockdown documents
- `.agents/preprint/locked/CONTROLS.md` — 5-control battery
- `.agents/preprint/locked/VERDICT_LEDGER.md` — locked verdict tag `strong trace, only D_coph` for α
- `.agents/preprint/methods/methods_revision_2026-05-18_cophenet.md` — binding methods directive
- `.agents/preprint/bands/01_beta.md` — β template (this brief mirrors its structure)

### Audit data (rows the numbers come from)
- `data/audit/ctm_triangle/cohort_summary.csv` — α row: C1, C2, C4 statistics
- `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` — α row: C3 cophenet
- `data/audit/alpha_epi_exclusion/cohort_summary.csv` — α C5 epi-X cohort summary
- `data/audit/alpha_epi_exclusion/comparison.csv` — α C5 per-patient Δρ table
- `data/audit/grassmann_cluster_extent/cohort_summary.csv` — α cluster-extent permutation (audit_70)
- `data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv` — α per-k T_G + surrogate stats
- `data/audit/grassmann_regate_no_filter/contig_summary.csv` — α audit_66 + audit_67 contiguous windows
- `data/audit/grassmann_epi_exclusion/cohort_summary.csv` — α per-k under C5 epi-X
- `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv` — α substrate row
- `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv` — α raw D sensitivity

### Audit scripts (provenance)
- `audit_63_split_baseline_surrogate.py` — C3 ρ_split^coph matched-strength (full cohort)
- `audit_68_alpha_epi_exclusion.py` — C5 ρ_split^coph epi-X at α
- `audit_66_grassmann_matched_strength_surrogate.py` — Grassmann C3 per-k
- `audit_67_grassmann_epi_exclusion.py` — Grassmann C5 per-k
- `audit_69_grassmann_regating_no_filter.py` — Grassmann contig-summary regating
- `audit_70_grassmann_cluster_extent.py` — Grassmann cluster-extent permutation (2026-05-19)

### Companion guides
- `.agents/guides/02_methods/lrg-framework-guide.md` — LRG primitive (formulas + codebase mapping)
- `.agents/guides/02_methods/imcoh-guide.md` — imcoh_abs convention + Nolte 2004 immunity
- `.agents/guides/04_rules/never-always-list.md` — coding rules + no-hardcoded-thresholds feedback
- `.agents/guides/05_plotting/README.md` — figure conventions (for upcoming α preprint figures)

## 11. What follows

This brief freezes the α verdict for the preprint manuscript under the 2026-05-18 / 2026-05-19 lockdown. Next steps:

1. **α figures** — produce F1 (4-control panel for `ρ_split^coph` + C5 epi-X), F2 (Grassmann `T_G(k)` cohort cells highlighting the empirical null), F3 (per-patient Δ_task vs Δ_rest scatter on `D_coph` for representative patients), and the per-patient C5 epi-X waterfall plot from §6. Output at `data/preprint/figures/alpha/`.
2. **No anatomy audit owed at α** under the locked battery. If a Desikan–Killiany enrichment is wanted as a descriptive supplement, the β anatomy pipeline (`audit_64`) can be re-pointed at α.
3. **No further sensitivity tests owed.** The α verdict is anchored at C3 paired Wilcoxon p = 0.00195 alone under the locked Wilcoxon-as-gate rule; the C5 epi-X audit (audit_68) provides supportive mechanistic context (non-epi-cortex contribution) per amended Decision 10.

Verdict ready for writing-agent handoff: **`strong trace, only D_coph`**, anchored at C3. The mechanistic narrative — α reorganization as a non-epileptogenic-cortex per-pair multiscale phenomenon — is supported by C5 but is not the primary verdict claim.
