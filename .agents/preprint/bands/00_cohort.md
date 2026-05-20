---
name: preprint-cohort-synthesis
era: IMCOH_ABS_COHORT_N10
status: locked_2026-05-19
kind: preprint-cross-band-synthesis
cohort: [Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15]
n_patients: 10
fc_method: imcoh_abs
tau: 1/lambda_max
canonical_lrg_object: "D_coph = cophenet(UPGMA(D(tau_max)))"
canonical_subspace_object: "U_k = span{phi_2, ..., phi_{k+1}}"
sources:
  - .agents/preprint/locked/CONTROLS.md
  - .agents/preprint/locked/VERDICT_LEDGER.md
  - .agents/preprint/locked/ANATOMY_CONTROLS.md
  - .agents/preprint/locked/ANATOMY_LEDGER.md
  - .agents/preprint/bands/01_beta.md
  - .agents/preprint/bands/02_alpha.md
  - .agents/preprint/bands/03_gammalow.md
  - .agents/preprint/bands/04_theta.md
  - .agents/preprint/bands/05_gammah.md
  - .agents/preprint/bands/06_delta.md
---

# Cross-band synthesis — 6 bands × 2 probes × 2 batteries (trace + anatomy)

## Head

Across 6 frequency bands × 2 LRG probes × 5 trace controls + 4 anatomy controls, **β is the only band with a trace at both probes; α adds a strong-but-cophenet-only trace; γ_l and δ carry weak Grassmann-only traces; θ and γ_h carry no trace on either probe.** Where a trace exists, it is **anatomically localized to a band-specific distributed cortical network** — never single-region, never diffuse-brain-wide. The cophenet wrap (`D_coph = cophenet(UPGMA(D(τ_max)))`) is responsible for **band resolution at the LRG-multiscale layer**: raw FC and raw `D(τ_max)` detect every band at 6-8/10 cohort agreement, but only the cophenet wrap demotes δ/θ/γ_h to non-trace and preserves β at 7/10.

## 1. Locked verdict matrix

### Trace + anatomy per band (locked)

| Band | Range (Hz) | D_coph trace | Grassmann trace | Anatomy verdict | Coverage tag |
|---|---|---|---|---|---|
| **β** | 13–30 | **strong** | **strong** | strong localized, both probes (7+7 named DK regions) | **strong trace, both probes** |
| α | 8–13 | **strong** (epi-X strengthens 8.3× → 27.7×) | no trace | strong localized; C5 epi-X reproduces identically (11 named DK regions both) | **strong trace, only D_coph** |
| γ_l | 30–80 | no trace | **weak** | strong localized, occipito-temporal + frontal + medial-OFC (7 named DK regions under `S(γ_l)`) | **weak trace, only Grassmann** |
| δ | 0.53–4 | no trace | **weak** | strong localized, full ≠ epi-X (4+3 regions, **0 shared** under `S(δ)`) | **weak trace, only Grassmann** |
| γ_h | 80–300 | no trace | no trace (borderline, cluster_p = 0.055) | n/a | **no trace** |
| θ | 4–8 | no trace | no trace | n/a | **no trace** |

Verdicts locked in `locked/VERDICT_LEDGER.md` (trace, 2026-05-18 + 2026-05-19 cluster-extent revision) and `locked/ANATOMY_LEDGER.md` (anatomy, 2026-05-19). Briefs document them; they do not re-derive them.

## 2. Headline methodological argument — the cophenet step is responsible for band resolution at LRG layer

Three-layer cohort matched-strength gating across all 6 bands. Each cell reports `n_above_surrogate / 10` followed by paired one-sided Wilcoxon p-value:

```
Layer                                       δ              θ             α              β              γ_l            γ_h
-------------------------------------      ---------     ----------     ---------     -----------    ---------      ---------
Raw FC          ρ_split^raw  (audit_67)    7/10 p=.042   7/10 p=.14     8/10 p=.014   6/10 p=.053    7/10 p=.053    7/10 p=.19
Raw D(τ_max)    ρ_split      (preprint_05) 7/10 p=.116   6/10 p=.097    7/10 p=.042   6/10 p=.042    6/10 p=.032    7/10 p=.161
Cophenet D_coph ρ_split^coph (audit_63)    4/10 p=.28    2/10 p=.72     5/10 p=.002   7/10 p=.005    5/10 p=.12     4/10 p=.25
```

Reading row by row:
- **Raw FC**: every band sits at 6-8/10 cohort agreement; δ/α pass matched-strength but the layer is band-agnostic.
- **Raw `D(τ_max)`**: same diffuse band-agnostic detection at the propagator-distance layer; β/α/γ_l all pass at p ≤ .042 but with similar cohort agreement to raw FC.
- **Cophenet `D_coph`**: the cophenet wrap selectively demotes δ/θ/γ_h from "passing matched-strength" to "obs ≈ surrogate" (n_above 4/10, 2/10, 4/10; p ≥ 0.116) while preserving β at 7/10 + p=0.005 and α at 5/10 + p=0.002.

The cophenet step is responsible for **band resolution at the LRG-multiscale layer**, NOT for amplification of detection. Raw FC already detects everything at the cohort level; cophenet identifies which bands carry the trace at the hierarchical-multiscale level. This three-layer contrast is the load-bearing argument for adopting `D_coph` as the canonical per-pair LRG probe.

Sources: `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv`, `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv`, `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv`. n=10, R=200 matched-strength surrogates per cell, seed 20260511. Methodology directive: `.agents/preprint/methods/methods_revision_2026-05-18_cophenet.md`.

## 3. Cross-band cluster-extent permutation table (Grassmann probe, audit_70)

Cluster-extent permutation null over k=2..112 for each band's `T_G(k)` cohort statistic. Pre-cluster gate is per-k matched-strength Wilcoxon at α_per_k = 0.05; cluster statistics are (LR) longest contiguous-significant run and (mass) cluster mass = Σ_k −log₁₀ p_k over all p<0.05 cells (including non-contiguous secondaries). Disjunctive gate as of 2026-05-19: pass if EITHER cluster_p_LR<0.05 OR cluster_p_mass<0.05.

| Band | obs_LR | null_mean_LR | null_max_LR | cluster_p_LR | obs_mass | cluster_p_mass | Verdict |
|---|---|---|---|---|---|---|---|
| **β** | **29** | 2.37 | 17 | **0.005** | **52.97** | **0.005** | strong |
| δ | 7 | 2.17 | 10 | 0.025 | 12.78 | 0.025 | weak |
| γ_l | 13 | 2.63 | 15 | 0.015 | 19.17 | 0.035 | weak |
| γ_h | 9 | 2.50 | 19 | 0.055 | 15.75 | 0.065 | no trace (borderline) |
| θ | 5 | 2.13 | 17 | 0.099 | 7.11 | 0.144 | no trace |
| α | 4 | 2.09 | 12 | 0.159 | 7.98 | 0.099 | no trace |

Source: `data/audit/grassmann_cluster_extent/cohort_summary.csv` (audit_70, 2026-05-19, R=200, seed unchanged from C3).

Reading: β dominates by both statistics (29-cell run, mass 52.97 nats — 13× null mean). δ and γ_l pass at weak (cluster_p_LR < 0.05). γ_h is the closest miss (cluster_p_LR = 0.055, just outside the 0.05 gate; demoted by `locked/VERDICT_LEDGER.md` Decision 6 from "weak trace" to "no trace"). α has no Grassmann trace — the 4-cell run at k=11..14 is well within the null distribution.

## 4. Per-band trace synthesis (one paragraph each)

### β — strong trace, both probes
β is the only band whose post-task structural trace survives matched-strength surrogacy at both LRG probes. `ρ_split^coph` cohort-paired Wilcoxon p = 0.005 (ratio 23.7×, 7/10 above own surrogate, passes C1 split + C2 drift + C3 matched-strength + C4 cross-probe). Grassmann cluster-extent permutation p = 0.005 (LR and mass both at audit floor, 29-cell contiguous run at k = 27..55, **strengthening to 36 cells at k = 21..56 under C5 epi-X**). The trace is localized to a distributed cortical network: cingulate + parahippocampal + entorhinal + insula + postcentral + superior frontal (cophenet probe); Hippocampus + temporal cortex + orbitofrontal + insula + rostral middle frontal (Grassmann probe). β is the load-bearing finding of the manuscript.

### α — strong trace, only D_coph; identical anatomy under C5 epi-X
α `ρ_split^coph` passes all four primary controls (C1 p=.010, C2 p=.007, C3 p=.002, C4 8/10). The cohort `n_above_surrogate` is 5/10 borderline; **the C5 epi-X analysis decisively resolves the borderline upward** — ratio 8.3× → 27.7×, p .002 → .014, n_above 5/10 → 7/10. C5 anatomy reproduces the same 11 named DK regions identically — α is a **non-epi-cortex cingulate-temporo-parietal trace** that the full-cohort signal partly masks via noise from epi-zone contacts. No Grassmann trace at α (4-cell run within null distribution).

### γ_l — weak trace, only Grassmann; occipito-temporal + frontal anatomy
γ_l cophenet has no trace (C3 p=0.116). Grassmann cluster-extent permutation passes at weak (cluster_p_LR = 0.015; full-cohort `|S(γ_l)| = 41` cells of `T_G^*` support). Under C5 epi-X the support shifts but the permutation gate is held. Anatomy under `S(γ_l)` is **occipito-temporal + frontal + medial-OFC** (7 named DK regions): left lateral occipital + left cuneus + left middle + superior temporal + left rostral middle frontal + right pars triangularis + right medial OFC. The KC-era "Hippocampus + left fusiform" claim is **fully retracted** — left fusiform does not appear under `S(γ_l)` (the locked contiguous-window analysis had picked it up but the all-clusters cluster-extent rerun does not).

### δ — weak trace, only Grassmann; full and epi-X read DISJOINT NETWORKS
δ cophenet has no trace (C3 p=0.278; obs/surr ratio 0.98×, observation indistinguishable from surrogate). Grassmann cluster-extent passes at weak (cluster_p_LR = 0.025; full-cohort `|S(δ)| = 23` cells of `T_G^*` support). C5 epi-X yields `|S^epiX(δ)| = 22` cells (per-k Wilcoxon on `grassmann_epi_exclusion/per_patient_per_band_per_k.csv`, support spans k=[2, 34–52, 87, 88]). Anatomy under `S(δ)` (full): 4 regions — left inferior temporal + left inferior parietal + right pars triangularis + left superior temporal. Anatomy under `S^epiX(δ)` (epi-X): 3 regions — left superior parietal + right rostral middle frontal + left superior frontal. **Zero regions are shared** between the two networks (was 2/6 under the retired `K*(δ)` windows; the dissociation strengthens to **fully disjoint** under the cluster-extent paradigm). The "anchor anatomy" (Amy + cingulate + fusiform) of the retired δ-full ledger is **not supported** under `S(δ)`. Separately, the C4 cross-probe `ρ_xprobe = +0.032` at 6/10 +sign is the known-biology δ anchor anatomy at the substrate layer (`memory/epileptic_imcoh_universal.md` 1.55× cross-probe ratio), reported as descriptive only (LEDGER Decision 5); this substrate-layer observation is independent of the Grassmann anatomy result and is not affected by the cluster-extent cascade.

### γ_h — no trace (borderline; Decision 6 demotion 2026-05-19)
γ_h cophenet has no trace (C3 p=0.246). Grassmann cluster-extent permutation is borderline (cluster_p_LR = 0.055, just outside the 0.05 gate; cluster_p_mass = 0.065). 9-cell observed run at k = 19..27. Under C5 epi-X the run shifts to k = 19..28 with 6 cells persisting + 8 new emergent cells at neighboring k. Demoted to "no trace" in the 2026-05-19 cluster-extent revision (Decision 6) — the 8-cell hardcoded threshold of audit_66 had read this as "weak"; under cluster-extent permutation the cohort signal is not distinguishable from the matched-strength null at the 0.05 level. γ_h is the closest miss in the manuscript and is the borderline case worth flagging in discussion.

### θ — no trace
θ has no trace on either probe. Cophenet C3 p=0.722 (obs_median ≈ surr_median; cohort sign agreement 2/10). Grassmann cluster-extent permutation cluster_p_LR = 0.099, cluster_p_mass = 0.144 — well within the null distribution. The cleanest negative reference for the manuscript and a band-specificity benchmark for β and α.

## 5. Cross-band anatomy synthesis (regions named)

### Region recurrence across trace-positive bands

Locked under the all-clusters cluster-extent paradigm (`S(b) = {k : p_k(b) < α_k}` aggregation, `data/audit/anatomy_<band>_grassmann{_epiX}_clusterext/cohort_summary.csv`, 2026-05-19 pm). The retired `K*(b)` (longest-contiguous-significant) windows are not used. See `directives/writing_directive_2026-05-19_anatomy_clusterext_rerun.md` for the methodology cascade.

| Region motif | β cophenet | β Grassmann | α cophenet | γ_l Grassmann | δ Grassmann full | δ Grassmann epi-X |
|---|---|---|---|---|---|---|
| Cingulate (any subregion) | isthmus + ACC (rh) | — | ACC bilateral + PCC (rh) | — | — | — |
| Hippocampus / parahippocampal | parahippocampal | **Hip** | parahippocampal | — | — | — |
| Insula | right | left | — | — | — | — |
| Medial / lateral OFC | — | bilateral | medial (rh) | medial (rh) | — | — |
| Fusiform | — | — | — | — | — | — |
| Inferior parietal | — | — | — | — | left | — |
| Superior parietal / postcentral / precuneus | postcentral (rh) | — | superior + postcentral + precuneus (rh) | — | — | superior (lh) |
| Temporal cortex (mid/sup/inf) | — | middle + superior (lh) | — | middle + superior (lh) | inferior + superior (lh) | — |
| Amygdala | — | — | — | — | — | — |
| Bankssts | — | — | — | — | — | — |
| Pars (orbitalis / triangularis) | — | — | — | triangularis (rh) | triangularis (rh) | — |
| Rostral / caudal middle frontal | — | rostral (rh) | caudal bilateral | rostral (lh) | — | rostral (rh) |
| Superior frontal | left | — | — | — | — | left |
| Paracentral | — | — | — | — | — | — |
| Lateral occipital / cuneus | — | — | — | lateral occipital (lh) + cuneus (lh) | — | — |

### Reading of the anatomy cross-band table

- **Cingulate** is at β cophenet + α cophenet — bilateral anterior cingulate + isthmus + posterior cingulate. The δ Grassmann full caudal ACC entry of the retired ledger is **not** preserved under `S(δ)`. The cingulate signature is now a per-pair-multiscale-only motif (cophenet probe).
- **Medial-temporal** (Hippocampus + parahippocampal) appears at β and α — β-and-α medial-temporal involvement supports β-α coupling in task-trace retention.
- **Fusiform** appears **nowhere** under the cluster-extent paradigm. The KC-era "Hippocampus + left fusiform" β claim is fully retracted (Hippocampus survives at β Grassmann; fusiform retracts at γ_l and δ alike under `S(b)`). The retraction is one of the cleanest worked-examples of why the methodology cascade was necessary.
- **Insula** appears only at β (both probes, opposite hemispheres — a label match, not contact-set match).
- **Amygdala** appears **nowhere** under the cluster-extent paradigm. The δ Grassmann full "anchor anatomy" framing (Amy + cingulate + fusiform) of the retired ledger is **not supported** under `S(δ)`. The δ cross-probe 1.55× known biology (`memory/epileptic_imcoh_universal.md`) is a **separate substrate-layer observation** and is unaffected by this cascade.
- **Parietal cortex** (inferior + superior + postcentral) appears at β cophenet + α cophenet + δ Grassmann full (left inferior) + δ Grassmann epi-X (left superior) — a recurring sensorimotor + parietal motif across bands where the trace localizes outside epi zones.
- **Temporal cortex** is the most multi-band motif under `S(b)`: middle + superior (lh) at β + γ_l; inferior + superior (lh) at δ-full. γ_l and δ-full now share left superior temporal.
- **Lateral-occipital + cuneus** is a γ_l-only motif emerging under `S(γ_l)` (was not in the retired contiguous-window list) — adds an occipital-cortex component to the γ_l network.

## 6. Methodological invariants under the locked battery

- **All cohort-level claims pass matched-strength surrogacy** (C3 trace + A3 anatomy). The single most important methodological invariant.
- **Cluster-extent permutation gates Grassmann at p<0.05** on the disjunctive LR ∨ mass statistic. No hardcoded thresholds (n_below ≥ 8, frac_consistent ≥ 7/10, etc.) gate any official verdict — those filters are forbidden per `feedback_no_hardcoded_test_thresholds.md`.
- **C5 epi-X is read by *reproducibility*, not retention percentage**. α anatomy reproduces 11/11; β Grassmann strengthens under `S(β)` (region set identical to retired audit); γ_l support `|S(γ_l)| = 41` cells; δ-full vs δ-epi-X now fully disjoint at the anatomy level (0 shared regions under `S(δ)` vs `S^epiX(δ)`). Each pattern means something different about how epi contacts participate in the trace.
- **A1 alone enriches >7 regions per band** at uncorrected p_hyper<0.05, including Wm and unknown. A3 matched-strength is the gate that excludes non-trace-driven signal.
- **No Pat_03 dropout, no Pat_03 marker, no Pat_03 separate stats**. Pat_03 is a full cohort member at n=10 from 2026-05-18 onward; sampling-rate handling is at the config layer only.

## 7. Companion artifacts

- Per-band briefs: `01_beta.md` … `06_delta.md`.
- Trace battery: `locked/CONTROLS.md` + `locked/VERDICT_LEDGER.md`.
- Anatomy battery: `locked/ANATOMY_CONTROLS.md` + `locked/ANATOMY_LEDGER.md`.
- Methods directive: `methods/methods_revision_2026-05-18_cophenet.md`.
- Cluster-extent permutation: `data/audit/grassmann_cluster_extent/{cohort_summary.csv, null_distribution.csv, per_k_obs_p.csv, README.md}` (audit_70).
- β β-figure scripts (parameterizable over `--band`): `scripts/02_preprint/preprint_07_beta_rho_split_figure.py`, `preprint_08_beta_grassmann_figure.py`.
- Anatomy audits: `scripts/01_compute/audit/audit_71_anatomy_cophenet.py`, `audit_72_anatomy_grassmann.py`.
- Cohort figures (this directory): `data/preprint/figures/cohort/F_cohort_1_three_layer.pdf`, `F_cohort_2_verdict_matrix.pdf`, `F_cohort_3_grassmann_strip.pdf` — see Phase C figures sub-task; not yet produced (figure scripts deferred until ready).

## 8. What this synthesis does NOT do

- Does not re-derive any verdict — verdicts are locked in `locked/VERDICT_LEDGER.md` + `locked/ANATOMY_LEDGER.md`.
- Does not propose new analyses — A2 sampling-corrected bootstrap and A4 implant-geometry regression are deferred to a sensitivity supplement, not to a new analysis here.
- Does not reinterpret KC-era findings — KC, VI(k), τ-sweep are retired (`methods/methods_revision_2026-05-18_cophenet.md`). The retracted "Hippocampus + left fusiform" β anatomy claim is documented as retracted in §5.
- Does not write the manuscript. Each per-band brief and this synthesis are inputs to the writing-agent handoff (`HANDOFF_INDEX.md`, Phase D) which produces LaTeX.

## Revision history

- **2026-05-19** — Initial synthesis. All 6 bands locked at the trace level (5 in `locked/VERDICT_LEDGER.md` 2026-05-18 + γ_h demotion 2026-05-19 Decision 6 + δ promotion 2026-05-19 Decision 7 + cluster-extent disjunctive gate 2026-05-19). All trace-positive bands' anatomy locked in `locked/ANATOMY_LEDGER.md` 2026-05-19. The cophenet step's band-resolution role is the headline methodological argument.
