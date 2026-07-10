---
name: preprint-handoff-index
era: IMCOH_ABS_COHORT_N10
status: locked_2026-05-19
kind: writing-agent-handoff
cohort: [Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15]
n_patients: 10
fc_method: imcoh_abs
tau: 1/lambda_max
canonical_lrg_object: "D_coph = cophenet(UPGMA(D(tau_max)))"
canonical_subspace_object: "U_k = span{phi_2, ..., phi_{k+1}}"
---

# Handoff index — writing-agent entry point (locked 2026-05-19)

**Head.** You are reading the **single entry point for the writing agent** that will produce the LaTeX manuscript on LRG analysis of sEEG functional connectivity. Everything you need is in `.agents/preprint/`. Every numerical claim, every verdict, every region, every probe is locked in a CSV or a ledger inside this directory. **You do not run code, you do not propose new analyses, you do not re-derive verdicts.** You read the locked artifacts, you write LaTeX prose, you cite CSV row pointers. If a number you want is not in a brief, the answer is that the number does not exist for the preprint — not that you should compute it.

## Read order

Read in this exact order. Skipping any step means you will produce prose that gets rejected at the evaluation pass (`EVALUATION_PROTOCOL.md`).

1. **This file** (`HANDOFF_INDEX.md`) — full briefing.
2. **Methods directive** (`directives/archive/2026-05/methods_revision_2026-05-18_cophenet.md`) — binding methodological language. KC retired, VI(k) retired, τ-sweep retired, `D_coph` adopted as canonical per-pair object. Read this before writing any methods paragraph.
3. **Grassmann methods companion** (`directives/methods_grassmann_cluster_extent.md`) — complete locked methodology for the Grassmann probe (subspace `U_k`, chordal distance, per-`k` Wilcoxon, cluster-extent permutation null with co-primary `LR` and `cluster_mass` statistics, disjunctive verdict gate, and the **band-level k-independent scalar `T_G^*(band) = Σ_{k ∈ C*} (−log10 p_k)`**). Read this before writing any Grassmann methods paragraph or before quoting a `T_G^*` value.
4. **Trace control battery** (`locked/CONTROLS.md`) — locked 5 controls (C1 split, C2 drift, C3 matched-strength, C4 cross-probe, C5 epi-X sensitivity). Read the disjunctive Grassmann gate (LR ∨ mass) in §C3, locked 2026-05-19.
5. **Trace verdict ledger** (`locked/VERDICT_LEDGER.md`) — per-band locked verdicts + transcribed CSV evidence + decision log. Decisions 1–7 are the audit trail for every verdict.
6. **Anatomy control battery** (`locked/ANATOMY_CONTROLS.md`) — locked 4 controls (A1 hypergeometric, A2 sampling-corrected [deferred], A3 matched-strength surrogate **mandatory**, A4 implant-geometry regression [deferred]).
7. **Anatomy verdict ledger** (`locked/ANATOMY_LEDGER.md`) — per-(band, probe) locked anatomy verdicts + locked region lists.
8. **Cross-band synthesis** (`bands/00_cohort.md`) — read the headline three-layer cohort table (§2) and the cross-band cluster-extent permutation table (§3). This is the headline methodological argument of the manuscript.
9. **Per-band briefs in order**: `bands/01_beta.md` (primary) → `bands/02_alpha.md` (primary) → `bands/03_gammalow.md` → `bands/04_theta.md` (negative reference) → `bands/05_gammah.md` → `bands/06_delta.md`. Read all six even for bands you might not write a separate paragraph on; they constrain how you frame the positive results.

After reading all 9, draft Results + Methods + Discussion sections; cite CSV rows for every number.

## ⚠️ ANATOMY — multi-region DK lists RETRACTED (2026-05-30); β REINSTATED → OFC system (2026-06-10); α/γ_l/δ no FDR-surviving localization

**Every locked multi-region DK *region-list* is retracted** (2026-05-30 signed sampling-aware
audit): no DK region reaches a defensible **cohort** localization (max coverage 5/10; locked
regions rest on 1–4 patients; several single-patient or anti-localized), and per-patient
localization is null at the single-DK-region level. The only above-chance *fine-grained* spatial
structure is **electrode-shaft autocorrelation** (the anatomy-free shaft partition reproduces it;
NMI 0.65), not single-region anatomy. **BUT the 2026-06-01 "spatially delocalized / no anatomy
anywhere" blanket was over-aggressive for β and is SUPERSEDED** (2026-06-10, `locked/ANATOMY_LEDGER.md`):
a sampling-conditioned, matched-strength-gated re-analysis at the a-priori *anatomical-system* scale
shows the brain-wide **β cophenet trace CONCENTRATES in the ORBITOFRONTAL CORTEX (OFC) system** —
concentration above a per-patient-demeaned baseline (a reproducible *hotspot*, not a container; even
within OFC the lean is ≈0.59), matched-strength **R=1000**, **BH q=0.009–0.013** across the 9 systems
within β, all four contact/shaft × incl/excl conditions, 4/5 implanted, leave-one-out + shaft-collapse
robust, low-strength, **bilateral**, peaks at the *system* scale (washes out at lobe/hemisphere;
undersampled at single-region). **β is the only band that localizes**; **α/γ_l/δ have no FDR-surviving
localization** (α pars opercularis fails BH over 53 regions; γ_l/δ null; the Grassmann *subspace* probe
localizes nowhere). Hippocampus/MTL = sub-threshold β hint (p≈0.4); occipital = K=3 secondary;
"distributed paralimbic ring" retracted (fails shaft-collapse). The multi-region DK *lists* below stay
retracted (historical record). Sources: `locked/ANATOMY_LEDGER.md`,
`data/audit/localization_atlas/README.md`, `data/audit/anatomy_localization_wilcoxon/README.md`,
`data/audit/per_patient_localization/README.md`. **Writing guidance:** for **β** write "**brain-wide
but concentrates / over-expresses in the orbitofrontal cortex (OFC) system**" (NOT "delocalized", NOT
"localized to OFC" as a container, NOT a multi-region DK list); for **α/γ_l/δ** "**no anatomical
localization** (undersampled / diffuse)" is correct. The anti-pattern is the old multi-region "localized
cortical network" DK list, and — for β only — the bare "delocalized / no anatomical anchor" blanket.

## What this manuscript IS about

- **The cophenet step is responsible for band resolution at the LRG-multiscale layer.** Raw FC and raw `D(τ_max)` detect a structural signal at every band at 6-8/10 cohort agreement. The cophenetic ultrametric `D_coph = cophenet(UPGMA(D(τ_max)))` selectively demotes δ/θ/γ_h to non-trace (obs ≈ surrogate under matched-strength) and preserves β + α as the trace-positive bands.
- **Two LRG probes**: per-pair `ρ_split^coph` on `D_coph` (cophenet per-pair multiscale) + whole-network `d_G(k)` on `U_k = span{φ_2..φ_{k+1}}` (Grassmann chordal distance on leading-k eigenmodes).
- **5 trace controls + 4 anatomy controls + 1 sensitivity layer (epi-X)**. Matched-strength surrogacy is **mandatory** for both batteries.
- **β is the primary finding** (strong trace, both probes; ~~strong localized anatomy at 7+7 named DK regions~~ **anatomy RETRACTED — delocalized, no cohort or per-patient localization**). α is the principal second finding (strong trace, only D_coph; strengthens under C5 epi-X to ratio 27.7×; ~~11-region anatomy reproduces identically under C5~~ **anatomy RETRACTED**). γ_l carries a strong Grassmann-only trace; δ a weak one (full-data LOO Pat_08-fragile, `p_mass = 0.005` at floor). θ + γ_h are negative references. **The trace verdicts stand; all per-band anatomy/region-lists are retracted (see banner above).**

## What this manuscript is NOT about

- **NOT about KC tree distance**. KC is retired (2026-05-18). Do not cite `T_KC(λ)`, do not refer to "matching cluster" or "Robinson-Foulds" distances in the body. If KC appears in an old artifact, ignore it.
- **NOT about VI(k) partitioning** at any k-cell level. VI is retired (2026-05-18).
- **NOT about τ-sweep multiscale structure**. The LRG spectrum on our continuous-spectrum FC is gap-less; τ-sweep is retired (2026-05-18); only `τ_max = 1/λ_max` is the canonical scale; the cophenetic dendrogram replaces gap-based scale identification.
- **NOT about Pat_03 outlier status**. Pat_03 is a full cohort member at n=10 from 2026-05-18 onward. Do not write "Pat_03 was excluded", "Pat_03 sensitivity test", "Pat_03 acquired at 1024 Hz (treated separately)". The 1024-Hz sampling rate is handled at the config layer; Pat_03 is treated identically at the analysis layer.
- **NOT about a multi-region DK-list localization.** The KC-era "Hippocampus + left fusiform" (retracted 2026-05-19), the cluster-extent 7+7 / 11-region DK lists (retracted 2026-05-30), and "Hippocampus survives at β Grassmann" (withdrawn 2026-06-01 — Hip = 3–5/10 marginal hint, not a localization) are **all retracted**. **But the trace DOES have one audited localization: β concentrates in the ORBITOFRONTAL CORTEX (OFC) *system*** (2026-06-10 — system-scale, brain-wide-with-OFC-hotspot, NOT a DK-region list; cite as concentration above per-patient baseline). α/γ_l/δ have **no FDR-surviving localization**; the Grassmann subspace probe localizes nowhere. See `locked/ANATOMY_LEDGER.md` + `data/audit/localization_atlas/README.md`.

## Anti-pattern checklist (run before submitting any draft)

These are the patterns that will cause your draft to be rejected at evaluation.

### Object-naming anti-patterns

- ❌ Writing `D(τ)` when the analysis is on `D_coph`. `D(τ)` is the raw propagator distance; the per-pair LRG probe is the cophenetic wrap of UPGMA of `D(τ_max)`. Always write `D_coph` or, in long form, `D_coph = cophenet(UPGMA(D(τ_max)))`.
- ❌ Writing "Grassmann distance on the LRG" without specifying `U_k`. The Grassmann probe is on the **leading-k Laplacian eigenmode subspace** `U_k = span{φ_2, …, φ_{k+1}}`, NOT on the LRG propagator or the cophenetic distance. Always state `U_k`.
- ❌ Writing "diffusion time" or "LRG time" when you mean `τ_max = 1/λ_max`. Use the latter.
- ❌ Writing "ultrametric matrix" without specifying it is `D_coph` (or its equivalent in code, `lrg.ultrametric_matrix` — which is `D_coph`, NOT raw `D(τ)`). See `feedback_never_confuse_D_with_cophenet.md`.
- ❌ Writing "linkage matrix `Z`" when you mean `D_coph`. `Z` is the linkage (tree structure); `D_coph = cophenet(Z)` is the dendrogram-derived per-pair ultrametric distance matrix.

### Methodology anti-patterns

- ❌ Calling cophenet "denoising". The cophenet step does **band resolution**, not denoising. Raw FC already detects every band; cophenet selects which bands carry the trace at the multiscale level.
- ❌ Claiming τ-sweep yields multiscale information. The LRG spectrum on our FC is gap-less (continuous-spectrum outlier case); the τ-sweep is degenerate; multiscale structure is captured by `D_coph` (per-pair across dendrogram merge heights) and by `U_k` (across subspace dimensions), not by a τ-sweep.
- ❌ Reintroducing KC or VI(k) at any layer.
- ❌ Stacking hardcoded patient-count thresholds (`n_below ≥ 8`, `frac_consistent ≥ 7/10`, `|med_surr|<0.05·|med_obs|`) on top of a statistical test. The Wilcoxon p-value or the cluster-extent permutation p-value **is** the gate; no additional counts gate any official verdict. See `feedback_no_hardcoded_test_thresholds.md`.
- ❌ Citing a number from memory or from an older report. Every number in your LaTeX must trace back to a CSV row cited in a brief. The brief tells you which CSV.

### Framing anti-patterns

- ❌ Framing `d_P = 1 − Pearson(triu A_a, triu A_b)` as "volume + topology" or as "orthogonal" to `d_S`. `d_P` is a magnitude-weighted complement to `d_S`; they correlate at ρ = 0.85–0.95 per band. The defensible triad is rank-only / magnitude-only / magnitude-weighted-complement. See `feedback_dP_framing.md`.
- ❌ Using bare "persistence" for `T_d < 0` findings. Use the **trace / anchor / reset / emergent** taxonomy. "Trace" = task reorganized AND change persists into rsPost; "anchor" = unchanged across all phases; "reset" = task reorganized AND module reverts; "emergent" = did not exist in rsPre. See `.agents/guides/01_project/terminology.md`.
- ✅ For **α / γ_l / δ**: "diffuse / no anatomical localization" is correct (no FDR-surviving region or system; undersampled, ≤5/10 coverage, per-patient null at single-region scale). ❌ Do NOT name DK regions or claim a "band-specific cortical network" for these bands (the multi-region lists are retracted). **β is different — β localizes to the OFC *system* (see the next bullet); do NOT call β "delocalized".**
- ✅ For **β**: write "the brain-wide trace **concentrates / over-expresses in the orbitofrontal cortex (OFC) system**" (2026-06-10 audited verdict). ❌ Do NOT write β is "delocalized / has no anatomical anchor" (superseded), nor "localized *to* OFC" as a container (it is brain-wide with an OFC hotspot — within-OFC lean ≈0.59), nor a multi-region DK "cortical network" list (retracted). For **α/γ_l/δ**: "**no anatomical localization** (undersampled / diffuse)" is correct; do NOT name DK regions. The precise β statement is *brain-wide across the sampled connectivity graph, concentrating above per-patient baseline in the OFC system at the anatomical-system scale (washes out at lobe/hemisphere, undersampled at single-region)*.
- ❌ Citing any of the **retracted** localizations — the KC-era "Hippocampus + left fusiform", the cluster-extent 7+7 / 11-region DK lists, or "Hippocampus survives at β Grassmann" (all retracted; Hip = 3–5/10 marginal hint, not a localization, null per-patient). The **one** anatomical localization you may cite is **β → OFC system** (2026-06-10, audited): cite it as a system-scale *concentration* above per-patient baseline, never as a single DK region or region list.

### Sycophancy anti-patterns

- ❌ Confidence laundering on unverified measures. Any measure that has not been matched-strength-controlled is **unverified** and cannot be cited as a primary claim. See `feedback_matched_strength_mandatory.md` + `feedback_brutal_honesty_no_sycophancy.md`.
- ❌ Writing "p < 0.05" without specifying which control. Every cohort-level p in the manuscript is paired one-sided Wilcoxon against a specific null; state the null.
- ❌ Writing "the analysis decisively shows X" without citing the controls that support X. The 5-control battery is the basis; cite which controls support the claim.

## Per-band one-paragraph cheatsheet

Use these as starting points for the band-specific Results paragraphs. The numbers are locked; the prose is yours.

- **β** (`bands/01_beta.md`): strong trace, both probes. `ρ_split^coph` p=0.005, ratio 23.7×, 7/10 above own surrogate, passes C1+C2+C3+C4. Grassmann cluster-extent permutation p=0.005, 29-cell run k=27..55, strengthening to 36 cells k=21..56 under C5. ~~Anatomy: 7+7 named DK regions across the two probes.~~ **Anatomy RETRACTED — delocalized: no cohort localization (Hip only, 3–5/10, marginal) and no per-patient localization on either probe (0/10).** Primary finding (trace, not anatomy).
- **α** (`bands/02_alpha.md`): strong trace, only D_coph; **anchored at C3 alone**. `ρ_split^coph` C3 paired Wilcoxon p=0.002 (verdict gate), ratio 8.3×; per-patient `n_above_surrogate` = 5/10 is descriptive only, NOT a verdict modifier (patient-count thresholds retired 2026-05-19, reaffirmed 2026-05-28). Secondary mechanistic observation (C5 epi-X): ratio strengthens 8.3× → 27.7×, Wilcoxon-on-epi-X p = 0.0098, 9 of 11 cophenet DK regions reproduce under audit_71 epi-X rerun — supportive of non-epi-cortex contribution; per amended Decision 10 / retracted Decision 1 (2026-05-28) NOT a verdict-promoter. Grassmann no trace (cluster_p=0.159; 4-cell run within null). ~~Anatomy: 11 named DK regions cophenet.~~ **Anatomy RETRACTED — diffuse (0/11 cohort-supported; all locked regions n≤3; per-patient also null).**
- **γ_l** (`bands/03_gammalow.md`): **strong trace, only Grassmann ↑**. Cophenet no trace (C3 p=0.116). Grassmann cluster-extent permutation **cluster_p_mass = 0.005** (gate; cluster_p_LR = 0.015 as descriptive co-statistic), support `|S(γ_l)| = 41` cells under the all-clusters paradigm. Full-data LOO max `p_mass = 0.040 (Pat_05)` — passes Decision-12 LOO precondition. C5 epi-X: trace contracts (raw mass 66.14 → 32.75, `p_mass^epi-X = 0.030 < 0.05` passes cohort gate; LOO under epi-X max = 0.159 Pat_05 is fragile, reported as secondary mechanistic observation). ~~Anatomy under `S(γ_l)`: 7 named DK regions Grassmann.~~ **Anatomy RETRACTED 2026-05-30 — 0/7 cohort-supported (well-sampled middletemporal n=5 fails); per-patient also null. (left fusiform appears nowhere, retracted.)**
- **θ** (`bands/04_theta.md`): no trace either probe. Cophenet C3 p=0.722 (obs_median ≈ surr_median; sign agreement 2/10). Grassmann cluster_p_LR=0.099, cluster_p_mass=0.144. The cleanest negative reference and band-specificity benchmark.
- **γ_h** (`bands/05_gammah.md`): no trace (Decisions 6 + 8 demotion 2026-05-19; `cluster_p_mass = 0.060` fails Decision-8 mass-only gate). Cophenet C3 p=0.246. Grassmann cluster_p_LR=0.055, cluster_p_mass=0.060 — borderline-just-outside-gate. Discussion-worthy as the closest miss.
- **δ** (`bands/06_delta.md`): **weak trace, only Grassmann**. Cophenet no trace (C3 p=0.278, obs ≈ surr ratio 0.98×). Grassmann **cluster_p_mass = 0.005** (cohort gate passes at floor; cluster_p_LR = 0.025 as descriptive co-statistic), support `|S(δ)| = 23` cells under the all-clusters paradigm. **Full-data LOO max `p_mass = 0.055 (Pat_08)` fails Decision-12 LOO precondition** — verdict is "weak" rather than "strong". C5 epi-X (secondary mechanistic observation, not verdict-promoter): raw mass strengthens 38.07 → 43.99, `p_mass^epi-X = 0.005`, LOO max under epi-X = 0.005 (Pat_02) fully robust — the full-data Pat_08 leverage is attributable to epi-zone interactions, not the true trace; `|S^epiX(δ)| = 22` cells. ~~Anatomy under `S(δ)` (full): 4 named DK regions; `S^epiX(δ)`: 3 named DK regions; "fully disjoint networks".~~ **Anatomy RETRACTED 2026-05-30 — 0/4 full and 0/3 epi-X cohort-supported (well-sampled regions anti-localized); per-patient also null; the epi-X mask was additionally a no-op bug (fixed). The "fully disjoint two networks" framing does not survive.** The C4 cross-probe `ρ_xprobe=+0.032` 6/10 +sign is descriptive δ anchor anatomy at the substrate layer only (1.55× ratio, known biology), not a Grassmann trace claim.

## Manuscript section guidance

### Methods section

- Define `imcoh_abs` (Ewald 2012; Bastos & Schoffelen 2016; Nolte 2004). Cite the volume-conduction-immune assumptions.
- Define LRG framework (Villegas 2023; Villegas 2025): `L̂`, `K̂(τ) = exp(−τL̂)`, `ρ̂_ij(τ) = K̂_ij / Z` with `Z = Σ_j K̂_ij`; `D(τ_max) = 1/ρ̂(τ_max)` with `τ_max = 1/λ_max`. State that the spectrum on our FC is continuous (no gap), so the canonical scale is `τ_max`; the cophenetic step replaces τ-sweep for multiscale structure.
- Define `D_coph = cophenet(UPGMA(D(τ_max)))`. State explicitly: in our codebase `lrg.ultrametric_matrix` is `D_coph`, **not** raw `D(τ)`.
- Define `U_k = span{φ_2, …, φ_{k+1}}` (skip trivial mode φ_1) and chordal Grassmann distance `d_G²(U_k^a, U_k^b) = k − ‖Q_k^a · Q_k^b‖_F²` where Q is the QR-orthonormalized basis.
- Define `ρ_split^coph`: per-pair Spearman correlation of Δ_task = `D_coph^task − D_coph^rsPre_A` with Δ_rest = `D_coph^rsPost − D_coph^rsPre_B` across all `N(N−1)/2` pairs.
- Define the 5-control battery from `locked/CONTROLS.md`. Lock the matched-strength algorithm (4-cycle ±δ swap, R=200, swap_target=20, seed=20260511).
- Define cluster-extent permutation (audit_70): two statistics LR (longest contiguous-significant run) and cluster_mass (Σ_k −log₁₀ p_k over all p<0.05); disjunctive gate.
- Define the 4-control anatomy battery from `locked/ANATOMY_CONTROLS.md`. A1 hypergeometric + A3 matched-strength surrogate are the primary controls; A2 sampling-corrected and A4 implant-geometry are deferred to sensitivity supplement.

### Results section

Lead with the three-layer cohort table (`bands/00_cohort.md` §2). State the central finding: cophenet does band resolution at LRG layer. Then per-band paragraphs in coverage-tag order: β > α > γ_l > δ > θ > γ_h. End with the locked verdict matrix (`bands/00_cohort.md` §1). **Anatomy: do NOT write per-band localization sub-sections — anatomy is retracted (no cohort or per-patient localization). Instead, a single short "spatial extent" subsection should report the trace as anatomically DELOCALIZED (no region anchor; the only spatial structure is electrode-shaft autocorrelation), citing `data/audit/anatomy_localization_wilcoxon/README.md` + `data/audit/per_patient_localization/README.md`.**

**Per-pair substrate figures (raw vs LRG, added 2026-06-03).** For the Results subsection that opens at the raw `|ImCoh|` substrate — the per-band joint-density portrait and the raw/cophenetic "null-triangle" pair — read `directives/writing_directive_2026-06-03_raw-substrate-results-figures.md`. It holds the approved captions, the raw-only and cophenetic prose, and the **verified raw→LRG mechanism: an inversion-MAGNITUDE collapse at α/β/γ_l, NOT a change in inversion count. Never write "the inversions turn gray / vanish" or "the cohort becomes unanimous"** — the dissenting patients persist, their inversions merely shrink. Figures: `data/preprint/figures/all_bands/`. Changes no verdict.

**Cophenetic per-pair subsection (the LRG twin, added 2026-06-03).** The all-band cophenetic-distance counterpart — the same joint-density portrait and null triangle on `D`, showing the per-pair trace clears the full control battery (matched-strength gate + drift + cross-probe) at **β (primary) and α (secondary) only** — is specced in `directives/writing_directive_2026-06-03_coph-substrate-results-figures.md`: caption specs, a three-paragraph narrative skeleton, the verified gate/drift/cross-probe numbers, and the guards. **Two guards the writing agent must respect: (i) "mainly β, secondarily α" is the matched-strength GATE verdict (β p=0.005, α p=0.002); (ii) γ_l clears the drift control but FAILS the gate and sits below its strength floor → it is NOT a third trace band.** Scripts: `preprint_18_bands_joint_density_rawfc.py --layer coph` + `preprint_21_bands_null_triangle.py --layer coph`. Changes no verdict.

### Discussion section

- **β concentrates in the orbitofrontal cortex (OFC) *system*** (2026-06-10 audited; brain-wide with an OFC hotspot at the anatomical-system scale — not a DK-region list). Discuss what an OFC-concentrated-but-brain-wide β structural trace means for task-related reorganization (a graph-level phenomenon with an anatomical-system hotspot, not a focal DK-region one); **α/γ_l/δ have no anatomical localization**. The retracted multi-region DK-list framing (β "medial-temporal + cingulate + insula network") must NOT be used, **nor** the superseded "β is delocalized / no anchor" blanket.
- α non-epi-cortex strengthening under C5 — most informative single methodological point about how epi contacts can mask cortical signals (a trace-strength observation, not a localization).
- The anatomy retraction itself is a worth-discussing methodological result: a direction-blind, top-decile, count-based enrichment manufactured apparent localization that a signed, threshold-free, cohort-consistency + per-patient test (and an electrode-shaft control) dissolved. (γ_l left-fusiform / KC-era claims all retracted.)
- δ cross-probe ratio 1.55× as known biology at the substrate (anchor) layer — not a Grassmann trace or localization claim.
- γ_h borderline (closest miss); θ as the cleanest negative reference.
- Methodological limitations: A2 + A4 deferred; epi-X is sensitivity not gate; cohort is sEEG-only; **localization is bounded by sEEG sampling — "no localization" means none detectable at achievable power within implanted tissue**.

## CSV row reference index (every number cited in any per-band brief)

Every numerical citation in your LaTeX must trace back to one of these CSV rows. The brief lists the row; the brief is the source of truth. If you need a number not in any brief, you do not have authorization to write it.

Master CSVs (cite from briefs, never from memory):

| CSV | Contents |
|---|---|
| `data/audit/ctm_triangle/cohort_summary.csv` | C1 split, C2 drift, C4 cross-probe per band; rho_split_median, n_above_drift, rho_xprobe_median, n_trace_xprobe_int |
| `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` | C3 matched-strength on `ρ_split^coph` per band; paired_wilcoxon_p, obs_median_rho, surr_median_rho, ratio, n_above_surrogate |
| `data/audit/alpha_epi_exclusion/cohort_summary.csv` | C5 cophenet epi-X (α only) |
| `data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv` | C3 Grassmann per-k T_G + surr stats |
| `data/audit/grassmann_cluster_extent/cohort_summary.csv` | C3 Grassmann cluster-extent permutation (audit_70) |
| `data/audit/grassmann_epi_exclusion/cohort_summary.csv` + `sensitivity.csv` | C5 Grassmann epi-X per-k + summary |
| `data/audit/grassmann_regate_no_filter/contig_summary.csv` | Grassmann contiguous-significant windows full + epi-X (audit_66 + audit_67) |
| `data/audit/raw_fc_matched_strength/cohort_summary_all_bands.csv` | Substrate layer of three-layer cohort table |
| `data/preprint/rho_split_raw_D/all_bands_matched_strength_cohort.csv` | Raw D layer of three-layer cohort table |
| `data/audit/anatomy_<band>_<probe>[_epiX]/cohort_summary.csv` | A1 hypergeometric + A3 matched-strength enrichment per DK region |

## What happens after you submit the draft

The user pastes your LaTeX into a Claude Code session. Claude runs the protocol in `EVALUATION_PROTOCOL.md` — number-by-number CSV cross-check, anti-pattern scan, verdict consistency check vs `locked/VERDICT_LEDGER.md` + `locked/ANATOMY_LEDGER.md`, anti-pattern checklist scan, figure-path coverage check. Findings come back as a delta vs the briefs.

If a number in your LaTeX doesn't match a CSV row, you are asked to fix it.
If a verdict in your LaTeX contradicts a ledger entry, you are asked to revert.
If an anti-pattern is detected, you are asked to rephrase.

## Revision history

- **2026-05-19** — Locked. Trace verdicts (5+1 controls) locked 2026-05-18 + 2026-05-19 cluster-extent revision. Anatomy verdicts (A1+A3) locked 2026-05-19. Per-band briefs and cohort synthesis locked 2026-05-19.
- **2026-06-03** — Added Results-section pointer to `directives/writing_directive_2026-06-03_raw-substrate-results-figures.md` (raw `|ImCoh|` per-pair joint density + raw/cophenetic null-triangle figures, captions, prose, and the verified raw→LRG inversion-magnitude-collapse mechanism). Figures only; changes no verdict.
- **2026-06-03** — Added the sibling pointer to `directives/writing_directive_2026-06-03_coph-substrate-results-figures.md` (cophenetic per-pair subsection: all-band joint density + null triangle on `D`, generalized `preprint_18 --layer coph` + `preprint_21 --layer coph`, caption specs, narrative skeleton, verified gate/drift/cross-probe numbers, β-primary/α-secondary guards). Figures + narrative only; changes no verdict.
