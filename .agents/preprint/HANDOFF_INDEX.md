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
2. **Methods directive** (`methods/methods_revision_2026-05-18_cophenet.md`) — binding methodological language. KC retired, VI(k) retired, τ-sweep retired, `D_coph` adopted as canonical per-pair object. Read this before writing any methods paragraph.
3. **Grassmann methods companion** (`methods/methods_grassmann_cluster_extent.md`) — complete locked methodology for the Grassmann probe (subspace `U_k`, chordal distance, per-`k` Wilcoxon, cluster-extent permutation null with co-primary `LR` and `cluster_mass` statistics, disjunctive verdict gate, and the **band-level k-independent scalar `T_G^*(band) = Σ_{k ∈ C*} (−log10 p_k)`**). Read this before writing any Grassmann methods paragraph or before quoting a `T_G^*` value.
4. **Trace control battery** (`locked/CONTROLS.md`) — locked 5 controls (C1 split, C2 drift, C3 matched-strength, C4 cross-probe, C5 epi-X sensitivity). Read the disjunctive Grassmann gate (LR ∨ mass) in §C3, locked 2026-05-19.
5. **Trace verdict ledger** (`locked/VERDICT_LEDGER.md`) — per-band locked verdicts + transcribed CSV evidence + decision log. Decisions 1–7 are the audit trail for every verdict.
6. **Anatomy control battery** (`locked/ANATOMY_CONTROLS.md`) — locked 4 controls (A1 hypergeometric, A2 sampling-corrected [deferred], A3 matched-strength surrogate **mandatory**, A4 implant-geometry regression [deferred]).
7. **Anatomy verdict ledger** (`locked/ANATOMY_LEDGER.md`) — per-(band, probe) locked anatomy verdicts + locked region lists.
8. **Cross-band synthesis** (`bands/00_cohort.md`) — read the headline three-layer cohort table (§2) and the cross-band cluster-extent permutation table (§3). This is the headline methodological argument of the manuscript.
9. **Per-band briefs in order**: `bands/01_beta.md` (primary) → `bands/02_alpha.md` (primary) → `bands/03_gammalow.md` → `bands/04_theta.md` (negative reference) → `bands/05_gammah.md` → `bands/06_delta.md`. Read all six even for bands you might not write a separate paragraph on; they constrain how you frame the positive results.

After reading all 9, draft Results + Methods + Discussion sections; cite CSV rows for every number.

## What this manuscript IS about

- **The cophenet step is responsible for band resolution at the LRG-multiscale layer.** Raw FC and raw `D(τ_max)` detect a structural signal at every band at 6-8/10 cohort agreement. The cophenetic ultrametric `D_coph = cophenet(UPGMA(D(τ_max)))` selectively demotes δ/θ/γ_h to non-trace (obs ≈ surrogate under matched-strength) and preserves β + α as the trace-positive bands.
- **Two LRG probes**: per-pair `ρ_split^coph` on `D_coph` (cophenet per-pair multiscale) + whole-network `d_G(k)` on `U_k = span{φ_2..φ_{k+1}}` (Grassmann chordal distance on leading-k eigenmodes).
- **5 trace controls + 4 anatomy controls + 1 sensitivity layer (epi-X)**. Matched-strength surrogacy is **mandatory** for both batteries.
- **β is the primary finding** (strong trace, both probes; strong localized anatomy at 7+7 named DK regions). α is the principal second finding (strong trace, only D_coph; strengthens under C5 epi-X to ratio 27.7×; 11-region anatomy reproduces identically under C5). γ_l carries a strong Grassmann-only trace; δ a weak one (full-data LOO Pat_08-fragile, `p_mass = 0.005` at floor). θ + γ_h are negative references.

## What this manuscript is NOT about

- **NOT about KC tree distance**. KC is retired (2026-05-18). Do not cite `T_KC(λ)`, do not refer to "matching cluster" or "Robinson-Foulds" distances in the body. If KC appears in an old artifact, ignore it.
- **NOT about VI(k) partitioning** at any k-cell level. VI is retired (2026-05-18).
- **NOT about τ-sweep multiscale structure**. The LRG spectrum on our continuous-spectrum FC is gap-less; τ-sweep is retired (2026-05-18); only `τ_max = 1/λ_max` is the canonical scale; the cophenetic dendrogram replaces gap-based scale identification.
- **NOT about Pat_03 outlier status**. Pat_03 is a full cohort member at n=10 from 2026-05-18 onward. Do not write "Pat_03 was excluded", "Pat_03 sensitivity test", "Pat_03 acquired at 1024 Hz (treated separately)". The 1024-Hz sampling rate is handled at the config layer; Pat_03 is treated identically at the analysis layer.
- **NOT about the KC-era anatomy claim "Hippocampus + left fusiform"**. Retracted 2026-05-19; further refined under the cluster-extent paradigm 2026-05-19 pm. Hippocampus survives at β Grassmann; **left fusiform appears nowhere under the cluster-extent paradigm `S(b)`** (was tentatively re-attributed to γ_l + δ under the retired `K*(b)` windows; the cluster-extent rerun shows it appears at neither). The retraction is in `locked/ANATOMY_LEDGER.md` 2026-05-19 pm revision and `bands/00_cohort.md` §5.

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
- ❌ Framing the β anatomy as "diffuse, not localized". Retired 2026-05-19. The β trace is **localized to a band-specific distributed cortical network** of 7+7 named DK regions across both probes under matched-strength surrogacy. "Localized" here means anatomically structured under A1+A3 (cophenet) or A3 alone (Grassmann); NOT single-region, NOT diffuse-brain-wide.
- ❌ Writing "the trace is brain-wide" — the trace is localized to a band-specific cortical network at every trace-positive band. See `locked/ANATOMY_LEDGER.md` for the locked region lists.
- ❌ Citing the KC-era β "Hippocampus + left fusiform" claim. Hippocampus survives at β Grassmann; left fusiform does NOT — under the locked cluster-extent paradigm `S(b)`, left fusiform appears at **neither** β, γ_l, nor δ. Fully retracted.

### Sycophancy anti-patterns

- ❌ Confidence laundering on unverified measures. Any measure that has not been matched-strength-controlled is **unverified** and cannot be cited as a primary claim. See `feedback_matched_strength_mandatory.md` + `feedback_brutal_honesty_no_sycophancy.md`.
- ❌ Writing "p < 0.05" without specifying which control. Every cohort-level p in the manuscript is paired one-sided Wilcoxon against a specific null; state the null.
- ❌ Writing "the analysis decisively shows X" without citing the controls that support X. The 5-control battery is the basis; cite which controls support the claim.

## Per-band one-paragraph cheatsheet

Use these as starting points for the band-specific Results paragraphs. The numbers are locked; the prose is yours.

- **β** (`bands/01_beta.md`): strong trace, both probes. `ρ_split^coph` p=0.005, ratio 23.7×, 7/10 above own surrogate, passes C1+C2+C3+C4. Grassmann cluster-extent permutation p=0.005, 29-cell run k=27..55, strengthening to 36 cells k=21..56 under C5. Anatomy: 7+7 named DK regions across the two probes (isthmus cingulate + rostral anterior cingulate + parahippocampal + entorhinal + insula + postcentral + superior frontal on cophenet; Hippocampus + temporal + orbitofrontal + insula + rostral middle frontal on Grassmann). Primary finding.
- **α** (`bands/02_alpha.md`): strong trace, only D_coph; **anchored at C3 alone**. `ρ_split^coph` C3 paired Wilcoxon p=0.002 (verdict gate), ratio 8.3×; per-patient `n_above_surrogate` = 5/10 is descriptive only, NOT a verdict modifier (patient-count thresholds retired 2026-05-19, reaffirmed 2026-05-28). Secondary mechanistic observation (C5 epi-X): ratio strengthens 8.3× → 27.7×, Wilcoxon-on-epi-X p = 0.0098, 9 of 11 cophenet DK regions reproduce under audit_71 epi-X rerun — supportive of non-epi-cortex contribution; per amended Decision 10 / retracted Decision 1 (2026-05-28) NOT a verdict-promoter. Grassmann no trace (cluster_p=0.159; 4-cell run within null). Anatomy: 11 named DK regions cophenet (bilateral cingulate + parahippocampal + medial OFC + caudal middle frontal + postcentral + precuneus + superior parietal).
- **γ_l** (`bands/03_gammalow.md`): **strong trace, only Grassmann ↑**. Cophenet no trace (C3 p=0.116). Grassmann cluster-extent permutation **cluster_p_mass = 0.005** (gate; cluster_p_LR = 0.015 as descriptive co-statistic), support `|S(γ_l)| = 41` cells under the all-clusters paradigm. Full-data LOO max `p_mass = 0.040 (Pat_05)` — passes Decision-12 LOO precondition. C5 epi-X: trace contracts (raw mass 66.14 → 32.75, `p_mass^epi-X = 0.030 < 0.05` passes cohort gate; LOO under epi-X max = 0.159 Pat_05 is fragile, reported as secondary mechanistic observation). Anatomy under `S(γ_l)`: 7 named DK regions Grassmann (**occipito-temporal + frontal + medial-OFC**: lateral occipital + cuneus + middle/superior temporal + rostral middle frontal + medial OFC + pars triangularis). The retired `K*(γ_l)` window had picked up left fusiform; **left fusiform retracts entirely under `S(γ_l)`** (and under `S(b)` everywhere).
- **θ** (`bands/04_theta.md`): no trace either probe. Cophenet C3 p=0.722 (obs_median ≈ surr_median; sign agreement 2/10). Grassmann cluster_p_LR=0.099, cluster_p_mass=0.144. The cleanest negative reference and band-specificity benchmark.
- **γ_h** (`bands/05_gammah.md`): no trace (Decisions 6 + 8 demotion 2026-05-19; `cluster_p_mass = 0.060` fails Decision-8 mass-only gate). Cophenet C3 p=0.246. Grassmann cluster_p_LR=0.055, cluster_p_mass=0.060 — borderline-just-outside-gate. Discussion-worthy as the closest miss.
- **δ** (`bands/06_delta.md`): **weak trace, only Grassmann**. Cophenet no trace (C3 p=0.278, obs ≈ surr ratio 0.98×). Grassmann **cluster_p_mass = 0.005** (cohort gate passes at floor; cluster_p_LR = 0.025 as descriptive co-statistic), support `|S(δ)| = 23` cells under the all-clusters paradigm. **Full-data LOO max `p_mass = 0.055 (Pat_08)` fails Decision-12 LOO precondition** — verdict is "weak" rather than "strong". C5 epi-X (secondary mechanistic observation, not verdict-promoter): raw mass strengthens 38.07 → 43.99, `p_mass^epi-X = 0.005`, LOO max under epi-X = 0.005 (Pat_02) fully robust — the full-data Pat_08 leverage is attributable to epi-zone interactions, not the true trace; `|S^epiX(δ)| = 22` cells. Anatomy under `S(δ)` (full): 4 named DK regions — **temporal + parietal + frontal** (inferior + superior temporal + inferior parietal + pars triangularis). Anatomy under `S^epiX(δ)`: 3 named DK regions — superior parietal + rostral middle frontal + superior frontal. **Zero regions shared** between the two networks (was 2/6 under the retired windows; dissociation strengthens to fully disjoint). δ Grassmann is a **mixture of two distinct phenomena** with no anatomical overlap. The retired "anchor anatomy" framing (Amy + cingulate + OFC + bankssts) is **not supported** under `S(δ)`. The C4 cross-probe `ρ_xprobe=+0.032` 6/10 +sign is descriptive δ anchor anatomy at the substrate layer only (1.55× ratio, known biology), not a Grassmann trace claim.

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

Lead with the three-layer cohort table (`bands/00_cohort.md` §2). State the central finding: cophenet does band resolution at LRG layer. Then per-band paragraphs in coverage-tag order: β > α > γ_l > δ > θ > γ_h. End with the locked verdict matrix (`bands/00_cohort.md` §1). Per-band anatomy goes in dedicated Results sub-sub-sections (one per trace-positive band) citing `locked/ANATOMY_LEDGER.md`.

### Discussion section

- β medial-temporal + cingulate + insula network and its meaning for task-related reorganization.
- α non-epi-cortex strengthening under C5 — most informative single methodological point about how epi contacts can mask cortical signals.
- γ_l left-fusiform + temporal-cortex network (KC-era retraction worked example).
- δ Grassmann two-phenomena reading (anchor + physiological); cross-probe ratio 1.55× as known biology.
- γ_h borderline (closest miss); θ as the cleanest negative reference.
- Methodological limitations: A2 + A4 deferred; epi-X is sensitivity not gate; cohort is sEEG-only.

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
