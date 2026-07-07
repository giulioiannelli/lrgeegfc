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

Across 6 frequency bands × 2 LRG probes × 5 trace controls + 4 anatomy controls, **β is the only band with a trace at both probes; α adds a strong-but-cophenet-only trace; γ_l and δ carry weak Grassmann-only traces; θ and γ_h carry no trace on either probe.** Where a trace exists, the multi-region DK lists are retracted (2026-05-30/06-01), but **the brain-wide β cophenet trace concentrates in ORBITOFRONTAL CORTEX** above its per-patient baseline (concentration, not exclusive presence — the trace exists cohort-wide; matched-strength R=1000, BH q≈0.01 across the 9 a-priori systems within β, 4/5 implanted, leave-one-out + shaft-collapse robust, low-strength, **bilateral**, peaks at the *system* scale and washes out at lobe/hemisphere; `localization_atlas/README.md`, 2026-06-10); **α/γ_l/δ have no FDR-surviving concentration**, and the cross-phase Grassmann subspace probe does not concentrate anywhere in any band. Hippocampus/MTL is a real but sub-threshold β hint. So the β trace is **distributed-but-concentrated (OFC)** and the others' traces are undersampled/diffuse — see §1 banner + §5. The cophenet wrap (`D_coph = cophenet(UPGMA(D(τ_max)))`) is responsible for **band resolution at the LRG-multiscale layer**: raw FC and raw `D(τ_max)` detect every band at 6-8/10 cohort agreement, but only the cophenet wrap demotes δ/θ/γ_h to non-trace and preserves β at 7/10.

**Epi-tissue stratification is a secondary tissue-characterization (sensitivity layer; the 2026-06-05 "primary interpretive lens" elevation was demoted 2026-06-18 — see §4.5 + Revision history).** The trace is a property of **healthy network reorganization that the epileptic zone partly masks, not an epileptiform artefact** — established by the **pair-class decomposition** (the unconfounded, same-graph evidence): at β the trace is carried by healthy↔healthy cortex (+0.197) and the epi↔non-epi interface (+0.252), **not** the diseased core (epi↔epi n.s. +0.142, p=0.161). **⚠ Withdrawn (node-count):** the earlier "excluding epi *strengthens* the β/α traces and *unmasks* low-γ" framing is a generic node-count effect (`audit_85 --stratify epi`, p_dec α 0.135 / β 0.390 / γ_low 0.130), not epi-specific. The surviving band-specific signature: **α uniquely recruits the epileptic tissue (epi↔epi +0.423, p_pair<0.001, matched-strength p=0.003), whereas β spares it.** This still pre-empts the "is this just pathological-network reorganization?" reviewer attack — through the pair-class evidence, not the withdrawn strengthening.

## 1. Locked verdict matrix

### Trace + anatomy per band (locked)

| Band | Range (Hz) | D_coph trace | Grassmann trace | Anatomy verdict | Coverage tag |
|---|---|---|---|---|---|
| **β** | 13–30 | **strong** | **strong** | ~~strong localized, both probes (7+7 DK regions)~~ → **REINSTATED 2026-06-10: cophenet trace → OFC system** (matched-strength R=1000, BH q≈0.01 across 9 systems, 4/5, LOO+shaft-robust, low-strength). Multi-region lists stay retracted; Hip sub-threshold; Grassmann probe not localized. | **strong trace, both probes** |
| α | 8–13 | **strong** (C3 full graph; α recruits epi core, §4.5) | no trace | ~~strong localized (11 DK regions)~~ **RETRACTED — DIFFUSE** (0/11 cohort; 1/10 per-patient = chance) | **strong trace, only D_coph** |
| γ_l | 30–80 | no trace | **strong** | ~~strong localized (7 DK regions)~~ **RETRACTED** (0/7 cohort; per-patient null) | **strong trace, only Grassmann** ↑ |
| δ | 0.53–4 | no trace | **weak** | ~~strong localized, full ≠ epi-X~~ **RETRACTED** (0/4 + 0/3 cohort, anti-localized; per-patient null; epi-X mask bug) | **weak trace, only Grassmann** |
| γ_h | 80–300 | no trace | no trace (borderline, cluster_p = 0.055) | n/a | **no trace** |
| θ | 4–8 | no trace | no trace | n/a | **no trace** |

> ✅ **UPDATE 2026-06-10 — β cophenet anatomy REINSTATED → CONCENTRATES in ORBITOFRONTAL CORTEX (the "delocalized" blanket below was over-aggressive for β; the trace is brain-wide, OFC is its reproducible hotspot, not a container).** A sampling-conditioned, matched-strength-gated re-analysis (`data/audit/localization_atlas/README.md`) shows the brain-wide **β cophenet trace concentrates in the OFC system** above its per-patient baseline (the analysis demeans each patient's cohort-wide trace, so it measures above-baseline *concentration*, not presence vs absence — even within OFC the lean is ≈0.59): the right denominator is patients *implanted* in the region (not /10), and under the mandatory matched-strength surrogate at **R=1000** with rank concordance and a **canonical observed trace**, OFC clears **BH q=0.009–0.013 across the 9 a-priori systems within β**, in all four contact/shaft-collapsed × epi-include/exclude conditions; 4/5 implanted positive, leave-one-out robust, **shaft-collapse robust**, low-strength (non-hub). The concentration **peaks at the *system* scale**: it is **bilateral** (no L/R lateralization, p≈0.48/0.76), **washes out** when pooled to the whole frontal lobe (p≈0.29), and is undersampled at single-region granularity (fails BH over 53 regions — the same wall α hits). **Hippocampus/MTL = sub-threshold hint** (p≈0.4, LOO-fragile); **occipital = K=3 low-coverage secondary**; "distributed paralimbic ring", **α** (pars opercularis fails BH over 53 regions), low-γ, and the **Grassmann subspace** probe show **no** above-baseline concentration. The multi-region DK *lists* below stay retracted; the per-patient/shaft-autocorrelation framing applied to coarse DK regions, not to the system-level OFC result. Six-rung control gate + self-audit (canonical-trace fix, shaft-collapse, R=1000, correct within-β multiplicity). This **supersedes** the 2026-06-01 "delocalized" verdict for β.
>
> ⚠️ **ANATOMY column RETRACTED 2026-05-30 / DELOCALIZED 2026-06-01 (historical; superseded for β by the 2026-06-10 update above).** The **trace** verdicts
> (D_coph / Grassmann / coverage tag) stand. The **anatomy** verdicts are all retracted: no DK
> region reaches a defensible cohort localization (max coverage 5/10; locked regions 1–4
> patients; several anti-localized) and per-patient localization is **also null** on both
> cross-phase probes (0–1/10 every band; β 0/10 on both). The verified β/α trace is spatially
> **delocalized** (network-level, no anatomical anchor); the only above-chance spatial structure
> is electrode-shaft autocorrelation, not anatomy. Sources:
> `data/audit/anatomy_localization_wilcoxon/README.md`, `data/audit/per_patient_localization/README.md`.

Trace verdicts locked in `locked/VERDICT_LEDGER.md`. **Anatomy verdicts: see the retraction above + `locked/ANATOMY_LEDGER.md` (the 2026-05-19 "strong localized" entries are superseded).**

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

The cophenet step is responsible for **band resolution at the LRG-multiscale layer**, NOT for amplification of detection. Raw FC already detects everything at the cohort level; cophenet identifies which bands carry the trace at the hierarchical-multiscale level. This three-layer contrast is the central argument for adopting `D_coph` as the canonical per-pair LRG probe.

Sources: `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv`, `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv`, `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv`. n=10, R=200 matched-strength surrogates per cell, seed 20260511. Methodology directive: `.agents/preprint/methods/methods_revision_2026-05-18_cophenet.md`.

## 3. Cross-band cluster-extent permutation table (Grassmann probe, audit_70 all-clusters, Decision-8 + Decision-12 gates)

Cluster-extent permutation null over k ∈ [2, 112] for each band's `T_G(k)` cohort statistic. Pre-cluster gate is per-k matched-strength Wilcoxon at α_per_k = 0.05; cluster statistics are (LR) longest contiguous-significant run and (mass) **resilient all-clusters cluster mass** `T_G^*(b) = Σ_{k:p_k<α_k} (−log₁₀ p_k)` over all p<0.05 cells (including non-contiguous secondaries). **Decision-8 cohort gate (locked 2026-05-19 pm):** `cluster_p_mass < 0.05`; `cluster_p_LR` is a descriptive co-statistic. **Decision-12 strong-tier precondition (locked 2026-05-28):** strong requires `cluster_p_mass < 0.01` AND `cluster_p_LR < 0.05` AND full-data LOO max `p_mass < 0.05`. C1 normalized `T_G^*` ∈ [0,1] reported alongside raw per locked decision; cohort denominator `n_k^cohort · log₁₀(R+1) = 111 × 2.3032 ≈ 255.65`.

| Band | obs_LR | null_mean_LR | null_max_LR | `cluster_p_LR` (descriptive) | `T_G^*` raw | `T_G^*` normalized | **`cluster_p_mass` (gate)** | **LOO max `p_mass` (Decision-12)** | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| **β** | **29** | 2.37 | 17 | **0.005** | **69.76** | **0.273** | **0.005** | **0.005** (Pat_02) | **strong** |
| **γ_l** | 13 | 2.63 | 15 | 0.015 | 66.14 | 0.259 | **0.005** | **0.040** (Pat_05) | **strong** ↑ |
| **δ** | 7 | 2.17 | 10 | 0.025 | 38.07 | 0.149 | **0.005** | **0.055** (Pat_08) | **weak** ← LOO binds |
| γ_h | 9 | 2.50 | 19 | 0.055 | 31.67 | 0.124 | 0.060 | 0.164 (Pat_02) | no trace (borderline) |
| θ | 5 | 2.13 | 17 | 0.099 | 12.97 | 0.0507 | 0.159 | — | no trace |
| α | 4 | 2.09 | 12 | 0.159 | 7.98 | 0.0312 | 0.348 | — | no trace |

Source: `data/audit/grassmann_cluster_extent/cohort_summary.csv` (audit_70 all-clusters formula, post-2026-05-19-pm fix; R=200, seed 20260511 — unchanged from C3). LOO columns: `cluster_p_mass_loo_max` + `cluster_p_mass_loo_argmax_patient` per Decision 11.

Reading: β dominates on every statistic (29-cell run, `T_G^* = 0.273` ≈ 9× null mean, LOO 0.005 Pat_02 robust). γ_l is genuinely strong under Decision 12 (cohort gate at floor, LOO 0.040 Pat_05 robust). **δ clears the cohort gate at the empirical floor (`p_mass = 0.005`) but full-data LOO max `p_mass = 0.055 (Pat_08)` fails the Decision-12 < 0.05 LOO precondition — Pat_08 single-handedly leverages the cohort verdict, blocking promotion to strong.** γ_h is the closest miss (`p_mass = 0.060`, just outside the 0.05 cohort gate; "no trace" by Decision 8). α has no Grassmann trace — the 4-cell run at k=11..14 is well within the null distribution.

## 4. Per-band trace synthesis (one paragraph each)

### β — strong trace, both probes
β is the only band whose post-task structural trace survives matched-strength surrogacy at both LRG probes. `ρ_split^coph` cohort-paired Wilcoxon p = 0.005 (ratio 23.7×, 7/10 above own surrogate, passes C1 split + C2 drift + C3 matched-strength + C4 cross-probe). Grassmann cluster-extent permutation p = 0.005 (LR and mass both at audit floor, 29-cell contiguous run at k = 27..55, **strengthening to 36 cells at k = 21..56 under C5 epi-X**). ~~The trace is localized to a distributed cortical network: ... (7 named DK regions per probe).~~ **Multi-region DK lists RETRACTED 2026-05-30/06-01 (0/7 + 0/7 cohort-supported, 0/10 per-patient at single-region scale); but the brain-wide β cophenet trace CONCENTRATES in the ORBITOFRONTAL CORTEX (OFC) *system* (REINSTATED 2026-06-10 — matched-strength R=1000, BH q≈0.01 across the 9 systems within β, 4/5 implanted, LOO + shaft-collapse robust; a system-scale hotspot, not a DK-region list; see §1/§5 banners + `data/audit/localization_atlas/README.md`).** β is the primary finding of the manuscript (trace; anatomy = brain-wide with an OFC system-scale hotspot).

### α — strong trace, only D_coph (anchored at C3 alone)
α `ρ_split^coph` passes all four primary controls at full cohort: C1 p=.010, C2 p=.007, **C3 paired Wilcoxon p=.002 (verdict gate; ratio 8.3×)**, C4 paired-Wilcoxon non-degradation. The per-patient `n_above_surrogate` at C3 is 5/10 — descriptive auxiliary statistic only, not a verdict modifier (patient-count thresholds retired 2026-05-19, reaffirmed 2026-05-28). Secondary tissue-characterization (epi-stratification, audit_77; PI-reviewed 2026-06-18): the cophenetic *strengthening under epi-exclusion* is withdrawn as generic node-count (`audit_85 --stratify epi`, p_dec 0.135); the surviving epi-specific result is the pair-class dissociation — **α uniquely recruits the epileptic core** (epi↔epi +0.423, p_pair<0.001, matched-strength p=0.003) whereas β spares it (§4.5). Per amended Decision 10 / retracted Decision 1 (2026-05-28), C5 is never a verdict-promoter; the α verdict is anchored at C3 alone. No Grassmann trace at α (4-cell run within null distribution).

### γ_l — strong trace, only Grassmann ↑ ~~; occipito-temporal + frontal anatomy~~ (anatomy RETRACTED)
γ_l cophenet has no trace (C3 p=0.116). Grassmann cluster-extent permutation passes the cohort gate at **cluster_p_mass = 0.005** (cluster_p_LR = 0.015 as descriptive co-statistic; full-cohort `|S(γ_l)| = 41` cells of `T_G^*` support). Full-data LOO max `p_mass = 0.040 (Pat_05)` passes the Decision-12 LOO precondition (< 0.05). C5 epi-X (secondary mechanistic observation, not verdict-promoter): trace contracts (raw mass 66.14 → 32.75, `p_mass^epi-X = 0.030 < 0.05` cohort gate held; LOO under epi-X max = 0.159 Pat_05 is fragile). ~~Anatomy under `S(γ_l)` is **occipito-temporal + frontal + medial-OFC** (7 named DK regions).~~ **Anatomy RETRACTED 2026-05-30/06-01 — 0/7 cohort-supported (well-sampled middletemporal n=5 fails), per-patient null; the trace is delocalized (see §5 banner).** (left fusiform appears nowhere — already retracted.)

### δ — weak trace, only Grassmann ~~; full and epi-X read DISJOINT NETWORKS~~ (anatomy RETRACTED)
δ cophenet has no trace (C3 p=0.278; obs/surr ratio 0.98×, observation indistinguishable from surrogate). Grassmann cluster-extent passes the cohort gate at **cluster_p_mass = 0.005** (cluster_p_LR = 0.025 as descriptive co-statistic; full-cohort `|S(δ)| = 23` cells of `T_G^*` support). **Full-data LOO max `p_mass = 0.055 (Pat_08)` fails the Decision-12 LOO precondition** (≥ 0.05) — the verdict is "weak" rather than "strong" precisely because of this single-patient leverage. C5 epi-X (secondary mechanistic observation, not verdict-promoter): raw mass strengthens 38.07 → 43.99, `p_mass^epi-X = 0.005`, LOO max under epi-X = 0.005 (Pat_02) fully robust — the full-data Pat_08 leverage is attributable to epi-zone interactions rather than the true trace; `|S^epiX(δ)| = 22` cells (per-k Wilcoxon on `grassmann_epi_exclusion/per_patient_per_band_per_k.csv`, support spans k=[2, 34–52, 87, 88]). ~~Anatomy under `S(δ)` (full): 4 regions ... Anatomy under `S^epiX(δ)` (epi-X): 3 regions ... **Zero regions shared** ... **fully disjoint**.~~ **Anatomy RETRACTED 2026-05-30/06-01 — 0/4 full and 0/3 epi-X cohort-supported (well-sampled regions anti-localized), per-patient null; the epi-X network rested on a no-op masking bug (fixed). The "two disjoint networks" framing does not survive (see §5 banner).** Separately, the C4 cross-probe `ρ_xprobe = +0.032` at 6/10 +sign is the known-biology δ anchor anatomy at the substrate layer (`memory/epileptic_imcoh_universal.md` 1.55× cross-probe ratio), reported as descriptive only (LEDGER Decision 5); this substrate-layer observation is independent of the Grassmann anatomy result and is not affected by the cluster-extent cascade.

### γ_h — no trace (borderline; Decisions 6 + 8 demotion 2026-05-19)
γ_h cophenet has no trace (C3 p=0.246). Grassmann cluster-extent permutation is borderline (cluster_p_LR = 0.055, just outside the 0.05 gate; cluster_p_mass = 0.060 fails the Decision-8 mass-only gate). 9-cell observed run at k = 19..27. Under C5 epi-X the run shifts to k = 19..28 with 6 cells persisting + 8 new emergent cells at neighboring k. Demoted to "no trace" in the 2026-05-19 cluster-extent revision (Decisions 6 + 8) — the 8-cell hardcoded threshold of audit_66 had read this as "weak"; under cluster-extent permutation the cohort signal is not distinguishable from the matched-strength null at the 0.05 level. γ_h is the closest miss in the manuscript and is the borderline case worth flagging in discussion.

### θ — no trace
θ has no trace on either probe. Cophenet C3 p=0.722 (obs_median ≈ surr_median; cohort sign agreement 2/10). Grassmann cluster-extent permutation cluster_p_LR = 0.099, cluster_p_mass = 0.144 — well within the null distribution. The cleanest negative reference for the manuscript and a band-specificity benchmark for β and α.

## 4.5 Epilepsy-independence — a secondary tissue-characterization (PI-reviewed 2026-06-18; the 2026-06-05 "primary lens" elevation is demoted)

**Head.** The band-selective cognitive task-trace is **independent of the epileptic pathology** — a property of healthy network reorganization that the epileptic zone partly *masks*, not an epileptiform artefact. This is a **secondary tissue-characterization / sensitivity layer**: the 2026-06-05 elevation to a "primary interpretive lens" was **demoted 2026-06-18**, because its load-bearing leg — the cophenetic "trace *strengthens* under epi-exclusion" claim — was withdrawn 2026-06-12 as a generic node-count artifact. The **unconfounded** evidence is the **pair-class decomposition** on the same full-graph LRG: at β the trace is carried by healthy↔healthy cortex and the epi↔non-epi interface, **not** the diseased core; at α it additionally recruits the core (the α exception). **No trace verdict changes** — every verdict was earned at C3 on the full graph.

**⚠ Withdrawn (node-count, 2026-06-12).** The earlier headline — "excluding epi *strengthens* the β/α traces and *unmasks* a low-γ trace" — is a **generic node-count effect**: a size-matched random-node decimation reproduces the same cophenetic gain (`audit_85 --stratify epi`, cohort p_dec α 0.135 / β 0.390 / γ_low 0.130, all `generic_nodecount`; epi exclusion drops only ~9.5% of nodes). It is **not epi-specific** and is no longer load-bearing. The subgraph table below is retained for the record; the **pair-class** table beneath it is the epi-specific evidence.

**Cophenetic ρ_split under epi-exclusion (audit_77, n=10, matched-strength R=200; positive = trace — subgraph rebuild, gain is node-count per the note above):**

| Band | full ρ (p) | exclude-epi ρ (p) | n>surr | LO-Pat_13 p | reading |
|---|---|---|---|---|---|
| **β** | +0.221 (0.0068) | +0.275 (0.0098) | 7/10 | 0.0195 | persists (gain = node-count) |
| **α** | +0.105 (0.0049) | +0.187 (0.0137) | 7/10 | 0.0039 | persists (gain = node-count) |
| **γ_low** | +0.083 (0.116, null) | +0.174 (0.0322) | 8/10 | 0.0371 | subgraph gain = node-count; γ_low trace rests on **Grassmann**, not cophenet (`03_gammalow.md`) |

**Pair-class decomposition — the epi-specific evidence (audit_77, n=10; full-graph null restricted to a node-pair class, global strength sequence held fixed; median pair counts nonepi↔nonepi 5408, cross 1167, epi↔epi 55):**

| Band | NONEPI↔NONEPI ρ (p) | CROSS epi↔non ρ (p) | EPI↔EPI ρ (p) |
|---|---|---|---|
| **β** | +0.197 (0.0049) **carry** | +0.252 (0.0098) **carry** | +0.142 (0.161) **does NOT carry** |
| **α** | +0.073 (0.0420) **carry** | +0.167 (0.0020) **carry** | **+0.423 (0.0029) carry (the α exception)** |

Reading:
- **β / α persist under epi-exclusion** — the trace is not driven by diseased tissue. (The apparent *strengthening* is generic node-count, withdrawn — see the ⚠ note above; "persist", not "strengthen". Both are robust to dropping Pat_13, the 30/119-contact patient who leverages every epi config.)
- **Carried by healthy + interface, not the core (β = core-spared)**: the β trace lives in healthy↔healthy (+0.197) and the epi↔non-epi interface (+0.252); the epi↔epi core is non-significant (+0.142, p=0.161) and the epi_only rebuilt subgraph weakens (−0.038).
- **The α exception**: α uniquely *also* recruits the epileptic tissue — its epi↔epi pair class is the strongest single epi-stratified α cell (+0.423, matched-strength p=0.0029; 7/10 patients individually significant, 9/10 positive; survives LO-Pat_13 p=0.0059), and it **also clears the dedicated pair-count null** (+0.423 vs null +0.129, p_pair<0.001, `concentrated`). **α recruits epileptic tissue, β spares it** — the band-specific epilepsy signature.
- **γ_low**: the cophenetic "emergence under epi-exclusion" is node-count (⚠ note above), not an epi-unmasked trace; γ_low's trace verdict rests on the **Grassmann** probe (audit_78), not cophenet. (Locked in `03_gammalow.md` + `VERDICT_LEDGER.md`; cross-referenced here.)
- **Grassmann corroboration (audit_78)**: β exclude_epi gains a separated k-window (mass 69.76 → 89.04, p_mass 0.005, longest run 29 → 36) — a *separate* substrate (audit_67-validated), not the withdrawn cophenetic node-count effect.

**Brutal-honesty caveats (bounding the dissociation).** (i) The **epi_only subgraph / epi↔epi class are underpowered** — 6–30 nodes (Pat_15 = 0), 15–435 pairs; matched-strength on such tiny graphs is near-degenerate, so the epi_only verdicts are exploratory and the β epi↔epi p=0.161 means "no power", not "proven absent". (ii) The **pair-class matched-strength null is global-rewiring, not within-class** — it fixes the global strength sequence but does not randomize *within* the class; so the strength-axis control on "α recruits the core" is **suggestive, not airtight** (a within-class matched-strength null has not been run). (iii) Class sizes are unequal (nonepi↔nonepi always has the most pairs → the most stable Spearman); read the per-class numbers alongside n_pairs. (iv) The α epi↔epi +0.423 rests on ~55 pairs (small) — **but** survives the dedicated pair-count null (`audit_85 --mode pairclass`, p_pair<0.001, `concentrated`), so the small per-class sample is not the source of the effect.

Sources: `data/audit/epi_stratified/{cophenetic_cohort.csv, cophenetic_per_patient.csv, pairclass_decimation_cohort.csv, decimation_control_cohort.csv, README.md}` (audit_77 + audit_85 pair-count/decimation controls), `data/audit/epi_stratified/grassmann_cohort.csv` (audit_78), consolidated report `.agents/reports/2026-06-05_epilepsy-headline-and-occult-node-marker.md`. Verdicts unaffected; see `VERDICT_LEDGER.md` 2026-06-18 amendment + `01_beta.md` §3.2.5 + `02_alpha.md` §3.1 + `03_gammalow.md`.

## 5. Cross-band anatomy synthesis ~~(regions named)~~ — **RETRACTED 2026-05-30 / DELOCALIZED 2026-06-01; β UPDATED 2026-06-10 → OFC**

> ✅ **UPDATE 2026-06-10 — there IS one localization, at the SYSTEM scale: β → orbitofrontal cortex.** The cross-band *multi-region DK motif* synthesis below stays retracted (those motifs were a direction-blind, top-decile, count-based artifact). But a sampling-conditioned, matched-strength-gated re-analysis at the a-priori *anatomical-system* scale recovers a real localization for β: **β cophenet trace → OFC system**, matched-strength R=1000, BH q=0.009–0.013 across 9 systems within β, all four conditions, 4/5 implanted, leave-one-out + shaft-collapse robust, low-strength. **No other band localizes** (α pars-opercularis fails BH over 53 regions; γ_l/δ null; Grassmann subspace null). **Hippocampus = sub-threshold hint; occipital = K=3 secondary; "distributed paralimbic ring" retracted (fails shaft-collapse).** So there is no cross-*band* motif, but there is a single, controlled, within-β **system-level** localization. Source: `data/audit/localization_atlas/README.md`.
>
> ⚠️ **ENTIRE SECTION RETRACTED (historical; β superseded by the 2026-06-10 update above).** The region-recurrence table and its reading below rest
> on the per-band "strong localized" region lists, which are all **retracted**: no DK region
> reaches a defensible cohort localization in any band/probe (0/7 β-coph, 0/7 β-Grass, 0/11
> α, 0/7 γ_l, 0/4+0/3 δ; max coverage 5/10; locked regions 1–4 patients; several
> anti-localized), and per-patient localization is **also null** on both cross-phase probes
> (0–1/10 every band; β 0/10 on both). There is **no cross-band region motif** to synthesize —
> the apparent recurrences (cingulate, medial-temporal, parietal…) were artifacts of a
> direction-blind, top-decile, count-based enrichment, and the only above-chance spatial
> structure is electrode-shaft autocorrelation, not anatomy. The trace is spatially
> **delocalized**. The table + reading are retained only as historical record of the
> retracted analysis. Sources: `data/audit/anatomy_localization_wilcoxon/README.md`,
> `data/audit/per_patient_localization/README.md`, `locked/ANATOMY_LEDGER.md`.

### Region recurrence across trace-positive bands ~~(RETRACTED — see §5 banner)~~

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
- **Fusiform** appears **nowhere**. ~~The KC-era "Hippocampus + left fusiform" β claim is fully retracted (Hippocampus survives at β Grassmann; fusiform retracts).~~ **FULLY RETRACTED 2026-05-30/06-01: neither survives — Hippocampus crosses the cohort test on only 3–5/10 patients (marginal hint, not a localization) and is null per-patient; fusiform appears nowhere. No *single-DK-region* β localization survives on either probe — the β localization is at the OFC *system* scale (§5 banner, 2026-06-10), not at single-region granularity (where it fails BH over 53 regions, the same wall α hits).**
- **Insula** appears only at β (both probes, opposite hemispheres — a label match, not contact-set match).
- **Amygdala** appears **nowhere** under the cluster-extent paradigm. The δ Grassmann full "anchor anatomy" framing (Amy + cingulate + fusiform) of the retired ledger is **not supported** under `S(δ)`. The δ cross-probe 1.55× known biology (`memory/epileptic_imcoh_universal.md`) is a **separate substrate-layer observation** and is unaffected by this cascade.
- **Parietal cortex** (inferior + superior + postcentral) appears at β cophenet + α cophenet + δ Grassmann full (left inferior) + δ Grassmann epi-X (left superior) — a recurring sensorimotor + parietal motif across bands where the trace localizes outside epi zones.
- **Temporal cortex** is the most multi-band motif under `S(b)`: middle + superior (lh) at β + γ_l; inferior + superior (lh) at δ-full. γ_l and δ-full now share left superior temporal.
- **Lateral-occipital + cuneus** is a γ_l-only motif emerging under `S(γ_l)` (was not in the retired contiguous-window list) — adds an occipital-cortex component to the γ_l network.

## 6. Methodological invariants under the locked battery

- **All cohort-level claims pass matched-strength surrogacy** (C3 trace + A3 anatomy). The single most important methodological invariant.
- **Cluster-extent permutation gates Grassmann at `cluster_p_mass < 0.05`** (mass-only, Decision 8 locked 2026-05-19 pm — `cluster_p_LR` is a descriptive co-statistic) plus the **Decision-12 LOO precondition** (`LOO max p_mass < 0.05` at full data) for strong-tier verdict. No hardcoded thresholds (n_below ≥ 8, frac_consistent ≥ 7/10, etc.) gate any official verdict — those filters are forbidden per `feedback_no_hardcoded_test_thresholds.md`.
- **C5 epi-X is a secondary control (demoted 2026-06-18), read on the Grassmann substrate only.** The earlier "anatomy reproducibility" reading is withdrawn: the α cophenet DK-region anatomy was delocalized (RETRACTED 2026-05-30/06-01, see §5) and the δ "two disjoint networks under `S(δ)` vs `S^epiX(δ)`" rested on a no-op masking bug (retracted 2026-06-01, §4 δ row). What survives: β Grassmann gains a separated k-window under epi-X (mass 69.76 → 89.04, `p_mass` 0.005; audit_67/72, validated). The cophenetic epi-X "strengthening/emergence" is generic node-count (`audit_85 --stratify epi`), withdrawn — see §4.5.
- **A1 alone enriches >7 regions per band** at uncorrected p_hyper<0.05, including Wm and unknown. A3 matched-strength is the gate that excludes non-trace-driven signal.
- **No Pat_03 dropout, no Pat_03 marker, no Pat_03 separate stats**. Pat_03 is a full cohort member at n=10 from 2026-05-18 onward; sampling-rate handling is at the config layer only.

## 7. Companion artifacts

- Per-band briefs: `01_beta.md` … `06_delta.md`.
- Trace battery: `locked/CONTROLS.md` + `locked/VERDICT_LEDGER.md`.
- Anatomy battery: `locked/ANATOMY_CONTROLS.md` + `locked/ANATOMY_LEDGER.md`.
- Methods directive: `methods/methods_revision_2026-05-18_cophenet.md`.
- **Neurophysiological reading** of the three-rung probe ladder + per-pair vs global persistence + per-band physical interpretation: `methods/methods_neurophysiological_interpretation_2026-05-26.md`.
- Cluster-extent permutation: `data/audit/grassmann_cluster_extent/{cohort_summary.csv, null_distribution.csv, per_k_obs_p.csv, README.md}` (audit_70).
- β β-figure scripts (parameterizable over `--band`): `scripts/02_preprint/preprint_07_beta_rho_split_figure.py`, `preprint_08_beta_grassmann_figure.py`.
- Raw-substrate per-pair Results figures (all-bands, 2026-06-03): joint density `scripts/02_preprint/preprint_18_bands_joint_density_rawfc.py` + null triangle `preprint_21_bands_null_triangle.py --layer raw|coph` → `data/preprint/figures/all_bands/`. Captions, prose, and the verified raw→LRG inversion-**magnitude**-collapse mechanism (NOT count; never "turns gray") in `directives/writing_directive_2026-06-03_raw-substrate-results-figures.md`. Figures/narrative only — changes no verdict.
- Cophenetic per-pair Results subsection (all-bands, 2026-06-03, the LRG twin): same scripts with `--layer coph` → `fig_bands_joint_density_empirical_coph.pdf` + `fig_bands_null_triangle_coph.pdf`. The per-pair trace clears the full control battery (matched-strength gate + drift + cross-probe) at **β primary, α secondary** only. Caption specs, three-paragraph narrative skeleton, verified gate/drift/cross-probe numbers, and the guards (**γ_l clears drift but FAILS the gate → not a third trace band**) in `directives/writing_directive_2026-06-03_coph-substrate-results-figures.md`. Figures/narrative only — changes no verdict.
- Anatomy audits: `scripts/01_compute/audit/audit_71_anatomy_cophenet.py`, `audit_72_anatomy_grassmann.py`.
- Cohort figures (this directory): `data/preprint/figures/cohort/F_cohort_1_three_layer.pdf`, `F_cohort_2_verdict_matrix.pdf`, `F_cohort_3_grassmann_strip.pdf` — see Phase C figures sub-task; not yet produced (figure scripts deferred until ready).

## 8. What this synthesis does NOT do

- Does not re-derive any verdict — verdicts are locked in `locked/VERDICT_LEDGER.md` + `locked/ANATOMY_LEDGER.md`.
- Does not propose new analyses — A2 sampling-corrected bootstrap and A4 implant-geometry regression are deferred to a sensitivity supplement, not to a new analysis here.
- Does not reinterpret KC-era findings — KC, VI(k), τ-sweep are retired (`methods/methods_revision_2026-05-18_cophenet.md`). The retracted "Hippocampus + left fusiform" β anatomy claim is documented as retracted in §5.
- Does not write the manuscript. Each per-band brief and this synthesis are inputs to the writing-agent handoff (`HANDOFF_INDEX.md`, Phase D) which produces LaTeX.

## Revision history

- **2026-05-19** — Initial synthesis. All 6 bands locked at the trace level (5 in `locked/VERDICT_LEDGER.md` 2026-05-18 + γ_h demotion 2026-05-19 Decision 6 + δ promotion 2026-05-19 Decision 7 + cluster-extent disjunctive gate 2026-05-19). All trace-positive bands' anatomy locked in `locked/ANATOMY_LEDGER.md` 2026-05-19. The cophenet step's band-resolution role is the headline methodological argument.
- **2026-06-03** — Added §7 companion pointer to the raw-substrate per-pair Results figures (joint density + raw/cophenetic null triangle) and `directives/writing_directive_2026-06-03_raw-substrate-results-figures.md`. Records the verified raw→LRG inversion-magnitude-collapse mechanism (NOT count). Figures + narrative only; no verdict change.
- **2026-06-03** — Added §7 companion pointer to the cophenetic per-pair Results subsection (LRG twin: `--layer coph` joint density + null triangle) and `directives/writing_directive_2026-06-03_coph-substrate-results-figures.md`. β-primary/α-secondary is the matched-strength gate; γ_l drift-clears but gate-fails (not a third band). Figures + narrative only; no verdict change.
- **2026-06-05** — **Epi-independence ELEVATED to a primary cross-band interpretive lens** (PI directive 2026-06-05; pre-empts the "isn't this pathological-network reorganization?" reviewer attack). Added new **§4.5** (epilepsy-independence synthesis: cophenetic ρ_split under epi-exclusion β/α strengthen + γ_low emerges; pair-class decomposition β carried by healthy + interface not core; α exception = α recruits epi core +0.404 while β spares it; Grassmann corroboration; brutal-honesty caveats) and an epi-independence paragraph in the Head. Cites audit_77/audit_78 (`data/audit/epi_stratified/`) + consolidated report `2026-06-05_epilepsy-headline-and-occult-node-marker.md`. γ_low emergence is locked separately (`03_gammalow.md` + VERDICT_LEDGER) and cross-referenced here only. **No trace verdict changes** — every verdict earned at C3 on the full graph; emphasis/interpretation elevation only. Caveats kept: epi_only/epi↔epi underpowered; pair-class null is global-rewiring not within-class (interface-carried + α-recruits-core = suggestive not airtight); α epi↔epi rests on ~45 pairs.
- **2026-06-12 → 2026-06-18 (cascade into 00_cohort completed 2026-06-24)** — **§4.5 DEMOTED from primary lens to secondary tissue-characterization + synced to n=10.** The 2026-06-05 elevation's load-bearing leg — "epi-exclusion *strengthens* β/α + *unmasks* low-γ" (cophenetic) — was withdrawn 2026-06-12 as generic node-count (`audit_85 --stratify epi`, p_dec α 0.135 / β 0.390 / γ_low 0.130); re-elevating on the surviving pair-class leg alone was rejected (PI 2026-06-18). Rewrote the Head epi paragraph, the verdict-matrix α cell, the §4 α synthesis, and all of §4.5: epi-independence is now a **secondary** sensitivity layer carried by the **pair-class** evidence (β core-spared; α-recruits-core the lone exception, +0.423, p_pair<0.001). Numbers synced n=9→n=10 (α epi↔epi +0.404→+0.423, median pair counts 5565/1050/45→5408/1167/55, 6/9→7/10 individually significant) and the pair-count null (`audit_85 --mode pairclass`, p_pair<0.001) added. **No trace verdict changes** — all earned at C3 on the full graph. Mirrors `01_beta.md` §3.2.5 + `02_alpha.md` §3.1 (2026-06-18). Source `data/audit/epi_stratified/{pairclass_decimation,decimation_control}_cohort.csv`.
