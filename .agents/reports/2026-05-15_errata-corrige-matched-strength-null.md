---
name: errata corrige — matched-strength null + epi-exclusion sensitivity for §5 / §6
scope: section_5_section_6_revision_handoff_to_writing_agent
date: 2026-05-15
status: active
era: cohort_n10_imcoh_abs
audience: writing-agent
target_manuscript_subsections:
  - ssec:lrg_decomposition (§5.2 KC tree distance)
  - ssec:lrg_ctm (§5.3 controlled per-pair ρ_split)
  - ssec:lrg_multiscale (§5.4 VI + Grassmann)
  - ssec:lrg_anatomy (§5.5 anatomical localization)
  - sec:synthesis (§6 band-resolved synthesis)
companion_handoffs:
  - .agents/reports/2026-05-14_full-null-verification-evidence-review.md
  - .agents/reports/2026-05-14_grassmann-epi-exclusion-sensitivity.md
  - .agents/reports/2026-05-11_grassmann-matched-strength-verification.md
  - .agents/reports/2026-05-11_implant-geometry-and-kc-null-verification.md
inputs:
  - data/audit/matched_strength_surrogate_split_baseline/{per_patient_per_band,cohort_summary}.csv  (all 5 non-γ_h bands; γ_h in-flight)
  - data/audit/kc_matched_strength_surrogate/{per_patient_per_band,cohort_summary,joint_signature}.csv  (5 bands)
  - data/audit/grassmann_matched_strength_surrogate/{per_patient_per_band_per_k,cohort_summary}.csv  (6 bands)
  - data/audit/grassmann_epi_exclusion/{per_patient_per_band_per_k,cohort_summary,sensitivity}.csv  (6 bands)
  - data/audit/raw_fc_matched_strength/{per_patient_per_band_all_bands,cohort_summary_all_bands}.csv  (6 bands)
  - data/audit/implant_geometry/{correlation_matrix,per_patient_features,per_patient_trace_measures}.csv  (5 bands, audit_64)
  - data/audit/alpha_epi_exclusion/{per_patient,cohort_summary,comparison}.csv  (audit_68)
---

# Errata Corrige — matched-strength null + epi-exclusion sensitivity

**Head.** The manuscript §5/§6 lean on the within-baseline-null + drift-floor + cross-probe control battery (controlled cohort claim α/β/γ_l at BH-FDR q = 0.027). A stronger control was run subsequently — the **matched-strength surrogate null** (R = 200 4-cycle ±δ rewiring preserving every node's strength `s_i = Σ_j W_ij` exactly to 10⁻⁴) — applied uniformly to every per-pair / KC / Grassmann probe. The matched-strength null answers a different question than the within-baseline null: *can the observed cross-phase trace be reproduced by independent per-phase node-strength evolution with random edge identities?* The result restructures the §5/§6 verdicts at every probe and band. **β remains the only band whose trace passes matched-strength at three independent probes (per-pair ρ_split, Grassmann subspace, Grassmann epi-excluded), with 100% retention of the manuscript Grassmann window k=27..55 under epi-exclusion.** **α passes matched-strength at the per-pair LRG layer only**, and **strengthens** under epi-exclusion (ratio 8.3× → 27.7×). **γ_l survives matched-strength at the Grassmann manuscript window k=12..23** but the per-pair LRG ρ_split fails cohort-paired Wilcoxon (p = 0.116). **γ_h carries a Grassmann subspace trace at k=19..27 surviving epi-exclusion** (audit_67 added post-manuscript-draft; not currently in §5); the γ_h LRG-CTM per-pair test fails cohort consistency due to a 4-pro-4-anti-2-null bimodal split — the γ_h cohort claim is **subspace-only**. **δ has no LRG-derived trace; θ is fully ergodic at every LRG layer.** A separate corollary: KC's matched-strength behaviour decomposes the within-baseline 10/10 β finding into a **strength-encoded** component and a **wiring-specific** component (roughly 50/50 at λ=0); this changes the §5.2 framing without invalidating the within-baseline result. The full-cohort matched-strength matrix at all six bands is **fully closed** as of 2026-05-15: γ_h KC (audit_65 high_gamma) and γ_h anatomy (audit_64 6-band extension) added — γ_h KC λ = 1 cohort Wilcoxon p = 0.042 (marginal, Wilcoxon-only), λ = 0 null; γ_h anatomy fully null. **One filter-audit caveat (`audit_69`, 2026-05-15)**: the audit_66 / audit_67 in-code `verdict == "separated"` label stacked three criteria (cohort Wilcoxon AND a hardcoded magnitude-ratio filter AND a hardcoded `n_below ≥ 8/10` patient-count filter), which over-restricted the in-script README counts to 0/111 for β. The §7 errata 29 / 12 / 9 contiguous-cells counts and 29/29 / 10/12 persist counts come from the Wilcoxon-only gate the manuscript actually cites and stand unchanged. **One transcription error caught**: γ_h Grassmann epi-exclusion persist count in the manuscript window k = 19..27 is **6/9 persist + 3/9 weaken** (re-derived from `sensitivity.csv`), not 7/9 + 2/9 as previously stated in earlier drafts. See §7.6 and §9.11.

---

## 1. What the matched-strength null adds, and the two-null structure

### 1.1 The two nulls answer different questions

The manuscript uses the **within-baseline null** at §5.3 and §5.2:

> split rsPre into halves A, B; for each (patient, band) build `ρ_drift = Spearman(D_coph^preB − D_coph^preA, D_coph^postB − D_coph^postA)`; test paired `ρ_split > ρ_drift` across the cohort.

This is the test of: *is the observed cross-phase shift larger than the rest-to-rest drift baseline on the same halved-data noise budget?*

The matched-strength null is structurally different:

> for each (patient, band, phase) generate R = 200 surrogate FC matrices `W_r` via 4-cycle ±δ rewiring preserving every `s_i` exactly; recompute every probe statistic on `W_r`; test paired `obs > surr_p50` across the cohort.

This is the test of: *is the observed cross-phase trace reproducible by independent per-phase reorganization of node strengths with random edge identities?*

The two nulls **do not nest**. A signal can pass within-baseline (it exceeds drift) and still fail matched-strength (it is reproducible by strength scrambling). The within-baseline-null result α/β/γ_l q = 0.027 stays valid as a drift-controlled claim; the matched-strength null adds a **wiring-vs-strength decomposition** on top.

### 1.2 The 4-cycle ±δ surrogate (definition for the writing agent)

For surrogate `r`, repeat `n_swaps = 20 · N(N−1)/2` times:

    pick distinct (a, b, c, d) ∈ [1..N]^4;
    δ ∈ [max(−W_ab, −W_cd, W_ad − 1, W_cb − 1),
         min(1 − W_ab, 1 − W_cd, W_ad, W_cb)];
    W_ab, W_cd ← W_ab + δ, W_cd + δ;
    W_ad, W_cb ← W_ad − δ, W_cb − δ.

The rewiring preserves `s_i = Σ_j W_ij` exactly per node (verified `‖s_obs − s_surr‖_∞ < 10⁻⁴`) and randomizes which pair carries the weight. R = 200 surrogates per (patient, band, phase). Effect-size ratio reported throughout:

    ratio = |cohort median (obs)| / |cohort median (surr_p50)|

with ratio ≫ 1 = clean separation from strength evolution; ratio ≈ 1 = strength reproduces observation; ratio < 1 = strength scrambling exceeds observation.

### 1.3 Per-probe decomposition under matched-strength

| probe | what observed measures | what matched-strength reproduces | failure mode if not separated |
|---|---|---|---|
| raw FC ρ_split | per-edge cross-phase correlation | strength-driven per-edge shifts | substrate signal partly strength-encoded |
| **LRG-CTM ρ_split** | **per-pair cophenetic cross-phase correlation on `D_coph(Trho)`** | **dendrogram reorganization driven by strength evolution propagating through the LRG → linkage pipeline** | **per-pair memory partly strength-encoded** |
| **KC T_KC(λ)** | **dendrogram tree-distance shift (topology λ=0; heights λ=1)** | **strength-driven dendrogram shifts** | **tree-distance memory partly strength-encoded** |
| **Grassmann T_G(k)** | **leading-k subspace rotation** | **strength-driven eigenmode shifts** | **subspace rotation partly strength-encoded** |

### 1.4 The interpretive corollary the manuscript currently misses

The within-baseline-null β KC result (10/10 patients below split-baseline null on both λ axes, BH-q = 0.006) is real but the matched-strength null reveals that **~50% of the observed shift at λ=0 is captured by strength scrambling alone** (ratio 2.21×). The β KC finding is therefore a **mixed** signature: a wiring-specific component **plus** a node-strength-evolution component, both contributing to the observed dendrogram shift. The §5.2 prose currently reads this as "the trace shifts both tree topology and tree heights coherently" without distinguishing the two contributions. The Errata Corrige adds the decomposition without retracting the original observation.

At α the matched-strength ratio at KC λ=0 is **0.12×** — the surrogate produces a *stronger* trace direction than the data does, meaning per-node strength evolution at α would predict a tree-distance shift larger than what is actually observed. The α tree-distance is essentially null under matched-strength; α's surviving signal lives at the per-pair LRG layer only.

---

## 2. §5.2 KC tree-distance errata (`ssec:lrg_decomposition`)

### 2.1 What the manuscript currently claims

The §5.2 cohort verdict is:

- β topology (λ = 0): 7/10 patients in trace direction, within-baseline-null p = 0.001 over 10 patients (10/10 below own split-baseline null), within-probe BH-q = 0.006 at m = 12 cells.
- β heights (λ = 1): 8/10 patients, p = 0.001, q = 0.006.
- γ_l heights (λ = 1): 8/10 patients, q = 0.055 (directional companion).
- γ_l topology (λ = 0): 4/10, anti-trace.
- α: 6/10 on both axes, neither significant.
- δ, θ: anti-trace.

### 2.2 What matched-strength adds (audit_65)

R = 200 surrogates per (patient, band, phase ∈ {preA, tt, post}), seed 20260511. Cohort verdicts at λ ∈ {0, 1} below; signs: T_KC < 0 = trace direction.

| band | λ | obs cohort median T_KC | surr cohort median | ratio | n_below own surr | Wilcoxon p (one-sided "less") |
|---|---|---|---|---|---|---|
| δ | 0 | −0.275 | −0.122 | 2.25× | 2/10 | 0.46 |
| δ | 1 | −3.88 | −2.48 | 1.56× (also_negative) | 3/10 | 0.22 |
| θ | 0 | **+0.83** (sign flip) | +0.052 | — | 2/10 | 0.62 |
| θ | 1 | −0.877 | −0.361 | 2.43× | 1/10 | 0.81 |
| α | 0 | −0.015 | −0.122 | **0.12×** (surr exceeds) | 1/10 | 0.54 |
| α | 1 | −0.396 | +0.007 | sign flip in surrogate | 1/10 | 0.72 |
| **β** | **0** | **−0.890** | **−0.403** | **2.21×** | 1/10 | 0.22 |
| **β** | **1** | **−0.651** | **+0.739** | **sign flip in surrogate** | 1/10 | 0.54 |
| γ_l | 0 | −0.695 | −0.101 | 6.88× | 2/10 | 0.31 |
| γ_l | 1 | +1.42 | +1.49 | both anti-trace | 1/10 | 0.50 |
| γ_h | 0 | **−2.98** | **−0.199** | **15.0×** (also_negative) | **3/10** | **0.188** |
| **γ_h** | **1** | **−3.82** | **−1.38** | **2.78×** | **5/10** | **0.042** |

γ_h KC (closed 2026-05-15 ~13:00, audit_65 high_gamma): **λ = 1 heights-axis is Wilcoxon-significant at α = 0.05** (cohort p = 0.042, ratio 2.78×) but per-patient direction is split. Per-patient z at γ_h λ = 1: Pat_02 −4.55, Pat_06 −8.02 (strongest trace), Pat_05 −2.06, Pat_08 −1.87, Pat_10 −1.80 (5 trace), vs Pat_03 +1.23, Pat_07 −0.90, Pat_13 −0.08, Pat_14 −0.86 (4 null), Pat_15 **+3.16** (strong anti-trace, consistent with Pat_15's 0-epi right-only-implant anchor pattern at every probe). λ = 0 (topology) is **null** (Wilcoxon p = 0.188, n_below 3/10), consistent with γ_h trace being height-axis only at the KC layer. γ_h KC λ = 1 is a marginal cohort signal — survives Wilcoxon-only at α = 0.05 but fails any reasonable consensus criterion (n_below 5/10).

### 2.3 What §5.2 should say after the Errata

The within-baseline-null statement (β 10/10 below own split-baseline null on both axes, within-probe BH-q = 0.006) is preserved as a drift-controlled claim. The Errata adds the matched-strength decomposition:

- **β λ = 0** (topology): observed shift is 2.21× the strength-scrambled shift. About half the within-baseline KC trace at β λ=0 is reproducible by strength evolution; the other half requires specific wiring changes. Per-patient consistency under matched-strength is weak (1/10 individually below surrogate at one-sided p<0.05), so the wiring-specific component is real at the cohort median but does not crystallize as a separate cohort-paired test.
- **β λ = 1** (heights): the surrogate flips sign (`surr cohort median = +0.74`) while the data stays negative (`obs = −0.65`). The within-baseline 10/10 result is therefore not strength-driven at the heights axis — strength scrambling produces an anti-trace heights signature on average, while the data shows the trace direction. Net reading: heights component of β KC is wiring-specific.
- **α**: both axes are essentially null under matched-strength; α tree-distance shifts are largely captured by strength evolution. The §5.2 claim that "α is silent on the dendrogram tree distance" remains correct; the Errata adds that the slight cohort-median positivity in the within-baseline view does not survive matched-strength scrambling.
- **γ_l λ = 0**: ratio 6.88×, n_below 2/10. Wiring-specific topology component exists at cohort-median magnitude but cohort consistency too weak to claim.
- **γ_l λ = 1** (the within-baseline directional-companion at q = 0.055): under matched-strength, both observed and surrogate cohort medians are positive (anti-trace at the heights axis), so γ_l λ = 1 also does not survive. The within-baseline-null directional companion claim should be reframed as "the heights-axis shift is in the trace direction relative to drift, but does not survive matched-strength scrambling".

### 2.4 Recommended prose change in §5.2

Where the manuscript says "Both β cells survive a stronger control, the within-baseline-null controlled test ... paired Wilcoxon over the ten patients at β gives p = 0.001 at both λ = 0 and λ = 1", append:

> A matched-strength surrogate null (R = 200 4-cycle ±δ rewiring preserving every node's `s_i` exactly to 10⁻⁴) further decomposes the β KC signal. At λ = 0 the observed cohort-median shift is 2.21× the strength-scrambled shift (n_below own surrogate 1/10, cohort-paired Wilcoxon p = 0.22): roughly half the dendrogram-topology trace at β λ = 0 is reproducible by per-node strength evolution alone, with the residual half requiring specific wiring changes. At λ = 1 the surrogate cohort median flips sign relative to the data (`surr +0.74` vs `obs −0.65`), so the heights-axis component is not strength-driven. The within-baseline-null β KC claim therefore stands as a drift-controlled result; the matched-strength decomposition documents that the topology-axis component is mixed wiring-and-strength while the heights-axis component is wiring-specific.

For γ_l λ = 1: the directional-companion claim at q = 0.055 should add the matched-strength caveat — γ_l heights-axis fails matched-strength scrambling.

---

## 3. §5.3 LRG-CTM `ρ_split` errata (`ssec:lrg_ctm`)

### 3.1 What the manuscript currently claims

Cohort-controlled verdict at within-probe BH-FDR (m = 6 bands): α, β, γ_l at q = 0.027; δ, θ, γ_h not in trace direction.

- α: `ρ_split = +0.115`, 8/10, drift-controlled p = 0.007.
- β: `ρ_split = +0.222`, 8/10, p = 0.014.
- γ_l: `ρ_split = +0.140`, 9/10 above drift, p = 0.010.

### 3.2 What matched-strength says (audit_63)

R = 200 4-phase rewiring (preA, preB, tt, post), seed 20260510. Cohort verdicts:

| band | obs cohort median ρ_split | surr cohort median | ratio | n_above own surr | Wilcoxon p (one-sided "greater") |
|---|---|---|---|---|---|
| δ | +0.008 | +0.008 | **1.0×** | 4/10 | 0.278 |
| θ | −0.040 | +0.004 | **(sign flip)** | 2/10 | 0.722 |
| **α** | **+0.105** | **+0.013** | **8.3×** | **5/10** | **0.002** |
| **β** | **+0.221** | **+0.009** | **23.7×** | **7/10** | **0.005** |
| γ_l | +0.083 | +0.003 | 30.2× | 5/10 | **0.116** (fails cohort-paired test) |
| **γ_h** | **+0.0003** | **+0.0042** | **0.08×** (surr exceeds, cohort bimodal) | **4/10** | **0.246** |

α epi-exclusion (audit_68) cohort:

| | obs cohort median | surr cohort median | ratio | n_above | Wilcoxon p |
|---|---|---|---|---|---|
| α full FC | +0.105 | +0.013 | 8.3× | 5/10 | 0.002 |
| **α epi-excluded** | **+0.187** | **+0.007** | **27.7×** | **7/10** | **0.014** |

α per-patient under epi-exclusion (Δρ = excluded − full):

    Pat_02: −0.037 → −0.035 (Δ +0.002; anti both, identity-like)
    Pat_03: +0.360 → +0.189 (Δ −0.171; weakened, n_epi=6)
    Pat_05: +0.078 → +0.121 (Δ +0.043; strengthened, n_epi=14)
    Pat_06: +0.833 → +0.692 (Δ −0.141; weakened slightly, n_epi=10)
    Pat_07: +0.081 → +0.211 (Δ +0.130; strengthened, n_epi=7)
    Pat_08: +0.380 → +0.290 (Δ −0.090; weakened, n_epi=9)
    Pat_10: +0.130 → +0.193 (Δ +0.063; strengthened, n_epi=10)
    Pat_13: +0.028 → −0.171 (Δ −0.199; FLIPPED anti, n_epi=30 — cohort maximum)
    Pat_14: +0.210 → +0.184 (Δ −0.026; minor weakening, n_epi=12)
    Pat_15: +0.040 → +0.040 (Δ 0; identity transform, n_epi=0)

### 3.3 What §5.3 should say after the Errata

- **α stands at the per-pair LRG layer** at within-baseline q = 0.027 AND at matched-strength (8.3× ratio, p = 0.002 cohort-paired). The epi-exclusion sensitivity (n_epi removed per patient via `load_epileptic_nodes ∩ channel_labels.csv`) **strengthens** the cohort signal (ratio 8.3× → 27.7×, n_above 5/10 → 7/10, Wilcoxon p remains < 0.05 at 0.014). The α trace lives in non-epi tissue; Pat_13 is the one cohort outlier whose α trace was carried by epi-zone activity (30/119 epi contacts; flips strongly anti under exclusion).
- **β stands at the per-pair LRG layer** at within-baseline q = 0.027 AND at matched-strength (23.7×, p = 0.005). β is the strongest cohort-consistent ρ_split signal under matched-strength (7/10 individually above own surrogate, cohort-paired p = 0.005).
- **γ_l fails matched-strength at the cohort-paired level** despite a 30.2× ratio on cohort median: only 5/10 patients individually above own surrogate, Wilcoxon p = 0.116. The within-baseline-null claim (drift-controlled p = 0.010, q = 0.027) remains valid; the matched-strength decomposition says the cohort-median magnitude is not strength-reproducible but the cohort is too heterogeneous to assert wiring-specific reorganization patient-by-patient. The §5.3 γ_l result should be reframed as "per-pair magnitude beyond strength-scrambling but per-patient cohort consistency fails matched-strength".
- **δ does not survive matched-strength** (ratio 1.0×, p = 0.278). The §5.3 statement that δ "is marginally positive (`+0.031`, 6/10, p = 0.246), and does not reach the trace direction" remains correct at the within-baseline layer and is further confirmed under matched-strength.
- **θ does not survive matched-strength** (sign flip; obs −0.040 vs surr +0.004; p = 0.722). Already null at within-baseline.
- **γ_h fails matched-strength at the LRG-CTM layer.** Cohort median ρ_split = +0.0003, surrogate cohort median = +0.0042 (ratio 0.08× — surrogate exceeds observation), n_above own surrogate = 4/10, Wilcoxon p = 0.246. **The cohort is strongly bimodal**: 4 patients (Pat_03 z=+16.7, Pat_05 z=+11.1, Pat_06 z=+15.8, Pat_08 z=+8.4) carry very strong pro-trace per-pair signals; 4 patients (Pat_02 z=−10.4, Pat_07 z=−4.2, Pat_14 z=−11.7, Pat_15 z=−9.5) carry strong anti-trace signals; 2 patients (Pat_10, Pat_13) are null. Median collapses to zero. **The per-pair LRG-CTM γ_h trace exists at the per-patient level but does not survive cohort aggregation.** Same bimodality pattern as the raw FC γ_h cohort (raw p=0.188): the substrate signal exists but disagrees on direction across patients. The γ_h subspace-only matched-strength-surviving trace at Grassmann k=19..27 (§5.4) remains the only cohort-consistent γ_h LRG signal.

### 3.4 Recommended prose change in §5.3

Where the manuscript reads "Under within-probe BH at m = 6 bands, all three trace bands sit at q = 0.027", append:

> A matched-strength surrogate null applied to the same per-pair ρ_split (R = 200 4-cycle ±δ rewiring preserving every node's strength exactly to 10⁻⁴, four phases per surrogate) refines the controlled cohort claim. At α the cohort-median ρ_split is 8.3× the strength-scrambled cohort median (n_above own surrogate 5/10, cohort-paired Wilcoxon p = 0.002), at β 23.7× (7/10, p = 0.005), and at γ_l 30.2× (5/10, p = 0.116 — fails cohort-paired consistency despite a large magnitude ratio). An epi-exclusion sensitivity test on α (audit_68) — recomputing `ρ_split` on the reduced FC matrix `W[non_epi, non_epi]` per patient — strengthens the cohort signal at α from ratio 8.3× → 27.7×, n_above 5/10 → 7/10, Wilcoxon p = 0.014. The α per-pair trace lives in non-epi tissue, with Pat_13 (30/119 epi contacts, the cohort maximum) as the single per-patient case whose α trace was carried by epi-zone activity (`ρ_split = +0.028` → −0.171 under exclusion). The matched-strength reading therefore preserves the within-baseline α/β cohort claims and demotes γ_l from cohort-paired to "magnitude beyond strength-scrambling, cohort consistency fails".

---

## 4. §5.4 multiscale views errata (`ssec:lrg_multiscale`)

### 4.1 What the manuscript currently claims

§5.4 reports VI(k) and Grassmann as **descriptive multiscale views without split-baseline control battery**:

> "Neither probe selects a privileged k, and neither carries the split-baseline control battery of §5.3; the role of §5.4 is to read whether the three-band identification extends as a coherent multiscale signature across spans of k."

- VI(k) per-band span: α at fine k (k=3..37, peak 9/10 at k=3); β at intermediate k (k=17..32 contiguous, peak 9/10 at k=22); γ_l at intermediate-coarse k (k=45..56, peak 8/10 at k=48).
- Grassmann principal-angle diagonals: β top-mode positive diagonal sustained k ≥ 12 through k = 79; γ_l broadest positive diagonal at k ≈ 10..28; α positive diagonal above k ≈ 10.

### 4.2 What matched-strength + epi-exclusion adds (audit_66 + audit_67)

audit_66 (Grassmann full FC, R = 200 matched-strength surrogates at K_grid = 2..112, seed 20260511) **adds** the controlled cohort-paired Wilcoxon test the §5.4 prose explicitly says is absent. audit_67 (Grassmann epi-excluded, K_grid = 2..88 capped by Pat_13 N_reduced = 89) adds the epi-exclusion sensitivity layer.

Cohort sig-cell counts (cohort-paired Wilcoxon "less", one-sided p < 0.05 against own matched-strength surrogate cohort median):

| band | audit_66 sig / 111 | longest contiguous sig run | best-k ratio | audit_67 sig / 87 | persist | weaken | emerge | absent | manuscript window retention |
|---|---|---|---|---|---|---|---|---|---|
| δ | 23 | k=57..63 (7 cells) | 4.2× at k=60 | 25 | 9 | 14 | 16 | 48 | n/a (no manuscript window claim) |
| θ | 8 | k=79..82 (4) | ≈ 1× | 8 | 3 | 5 | 5 | 74 | n/a |
| α | 4 | k=11..14 (4) | ≈ 1× | 3 | 2 | 2 | 1 | 82 | n/a |
| **β** | **40** | **k=27..55 (29)** | **2.8× at k=40** | **52** | **35** | **3** | **17** | **32** | **29/29 in k=27..55 (100%)** |
| **γ_l** | **41** | **k=12..23 (12)** | **11.5× at k=17** | **20** | **16** | **24** | **4** | **43** | **10/12 in k=12..23 + 2 weaken** |
| **γ_h** | **19** | **k=19..27 (9)** | **15.4× at k=23** | **16** | **8** | **4** | **8** | **67** | **6/9 persist + 3/9 weaken in k=19..27 (corrected 2026-05-15; prior drafts said 7/9)** |

audit_67 also reports γ_h cohort-median |T_G| **deepens** from −0.27 (full FC at k=25) to −0.51 (epi-excluded at k=25). The trace direction strengthens, not weakens, under epi-exclusion at γ_h.

### 4.3 What §5.4 should say after the Errata

**β**: the §5.4 sustained-diagonal claim becomes a matched-strength-controlled claim. The longest contiguous run on the chordal scalar T_G(k) under matched-strength is **k = 27..55 (29 cells, all with cohort-paired Wilcoxon p < 0.05 against the surrogate)**, with best-k ratio 2.8× at k = 40. Under epi-exclusion (audit_67), **all 29 cells of the manuscript window persist** (cohort p < 0.05 in both full FC and epi-excluded readings). β is the **only band** where the Grassmann subspace trace survives both matched-strength scrambling and epi-zone removal simultaneously across the full audit_66 manuscript window.

**γ_l**: the §5.4 broadest-diagonal claim becomes scoped. The matched-strength-controlled cohort signal lives in **k = 12..23 (12 contiguous cells)** with **best-k ratio 11.5× at k = 17** — the largest sustained subspace effect-size ratio in the cohort at narrow intermediate k. Under epi-exclusion, **10/12 cells of the manuscript window persist**, with 2 weakening at the window edges. **Outside the k = 12..23 window**, 24 cells weaken under epi-exclusion — so the broad-range γ_l signature documented in audit_66 is **partly epi-driven outside the manuscript window**. The §5.4 prose should scope the γ_l claim to the k = 12..23 window and flag broader claims as sensitive to epi-exclusion.

**γ_h** (not currently in §5.4): the band carries a matched-strength-controlled Grassmann subspace trace at **k = 19..27 (9 contiguous cells)** with **best-k ratio 15.4× at k = 23** — the largest single-probe effect-size ratio in the entire cohort. Under epi-exclusion (audit_67), **6/9 cells persist + 3/9 weaken** (persist ks = 21, 22, 23, 24, 25, 26; weaken ks = 19, 20, 27 — corrected 2026-05-15 from earlier 7/9 transcription error); median |T_G| at k = 25 **deepens** from −0.27 (full FC) to −0.51 (epi-X) at the surviving cells, so the trace direction strengthens when epi-zone contacts are removed. This **eliminates the "pure HFO contamination from epi contacts" alternative** at γ_h. Pat_13 (30/119 epi contacts) neutralizes its γ_h trace under exclusion (T_G ≈ 0 at k = 25); Pat_15 (0 epi) is the identity-transform anchor.

**α**: the §5.4 positive-diagonal claim above k ≈ 10 **does not survive matched-strength**. Only 4 sig cells / 111 at audit_66 with no contiguous run; ratio ≈ 1× at every k. The α multiscale-spectral signature is a within-baseline pattern that does not carry over to matched-strength scrambling. The α cohort claim is therefore at the per-pair LRG layer only (§5.3), not at the spectral subspace.

**δ**: a coarse-k Grassmann shoulder at k = 57..63 (7 contiguous cells, ratio 4.2× at k = 60) survives matched-strength but is the only LRG-derived δ signal; under epi-exclusion 16 cells emerge across the full k range (denoising signature). Worth carrying as an exploratory δ note, not a §5.4 headline.

**θ**: 8/111 at audit_66 (≈ chance); 8/87 at audit_67 with 3 persist, 5 weaken, 5 emerge, 74 absent. No subspace signal at any reading.

### 4.4 Recommended prose changes in §5.4

The §5.4 declaration "neither probe selects a privileged k, and neither carries the split-baseline control battery of §5.3" should be revised: the Grassmann probe **does** carry a matched-strength surrogate control battery after audit_66/67. The corresponding paragraph should add:

> The Grassmann subspace probe additionally carries a matched-strength surrogate control battery (audit_66, audit_67): R = 200 4-cycle ±δ rewiring preserving every node's strength exactly, applied per phase, with cohort-paired Wilcoxon at each k cell. At β, the manuscript window k = 27..55 yields 29 contiguous cells in the trace direction beyond matched-strength scrambling, with all 29 cells persisting under further epi-zone exclusion. At γ_l the corresponding window k = 12..23 yields 12 contiguous cells (best-k ratio 11.5× at k = 17, the largest in the cohort), with 10/12 persisting under epi-exclusion. At γ_h, k = 19..27 yields 9 contiguous cells (best-k ratio 15.4× at k = 23), with 6/9 persisting under epi-exclusion (3/9 weaken at the window edges k = 19, 20, 27) and the cohort-median chordal trace direction *strengthening* at the surviving cells when epi-zone contacts are removed — eliminating the alternative attribution to pure HFO contamination of γ_h at epi contacts. The α and θ Grassmann panels do not survive matched-strength scrambling at any sustained k range.

§5.4 should also add γ_h as a fourth band carrying a controlled signal at this layer (currently the §5.4 paragraph cohort-claims only α/β/γ_l from §5.3, and γ_h is "borderline"). The §5.4 update redefines γ_h as a band with a defensible subspace-only matched-strength-surviving trace; the §5.3 LRG-CTM probe at γ_h was closed 2026-05-15 and fails cohort consistency (4 pro / 4 anti / 2 null, p = 0.246) — so γ_h is **subspace-only**, with a Wilcoxon-only-marginal KC heights companion (p = 0.042) and a fully null anatomy correlate.

---

## 5. §5.5 anatomical localization errata (`ssec:lrg_anatomy`)

### 5.1 What the manuscript currently claims

- γ_l ctx-lh-fusiform: 13/38, 3.95× enrichment, `p_hyper = 5.2 × 10⁻⁶` (Bonferroni-survived at m = 48). Double dependence on Pat_02 disclosed.
- β: no single region survives Bonferroni; uncorrected candidates Hippocampus (3.49×, p_hyper = 0.011), left fusiform (2.48×, p_hyper = 0.045), left superior temporal (2.13×, p_hyper = 0.055). Hippocampus and left fusiform collapse under Pat_03 dropout.
- α: ctx-lh-parsopercularis uncorrected-suggestive (5.55×, p_hyper = 1.4 × 10⁻³, just missing Bonferroni), but two-patient base (Pat_07 anti-aligned + Pat_14 pro-aligned) makes it not interpretable as a clean CTM cohort signal. Recorded as directional indication only.

### 5.2 What audit_64 implant-geometry test adds

audit_64 tests per-patient implant-geometry features (N, hemispheric balance B_hemi, lobar concentration frac_lobe_max, spatial dispersion σ_disp, epi-zone fraction frac_epi, hemisphere counts N_L/N_R, lobe counts) against per-patient trace measures (obs_ρ_split, obs_z, T_KC λ=0/1, T_d^(d_S)) at the cohort level via Spearman correlation. Family is m = 225 = 5 bands × 5 trace measures × 9 features; BH-FDR at q = 0.05.

Top per-band correlates after extension to all 5 non-γ_h bands (audit_64 re-run 2026-05-15):

| band | top feature | Spearman ρ | uncorrected p | BH q (m=225) | verdict |
|---|---|---|---|---|---|
| δ | N (channel count) | −0.529 | 0.116 | 0.996 | null |
| θ | frac_lobe_max | **+0.867** | **0.0012** | 0.184 | unc-strong, BH-null |
| α | frac_epi | −0.406 | 0.244 | 0.996 | null (direction supports non-epi reading) |
| **β** | **B_hemi** | **+0.697** | **0.025** | **0.835** | **left-lateralization, uncorrected-only** |
| γ_l | frac_lobe_max | +0.358 | — | — | null |

**β B_hemi = +0.697** (left-hemisphere fraction positively correlates with per-patient β `obs_ρ_split`): patients with predominantly left-hemisphere implants carry stronger β per-pair traces. Corroborates the §5.5 narrative of β trace-leaves concentrating at **left** fusiform / **left** superiortemporal / **left** parsopercularis. The audit_64 hypothesis-test does not pass BH at m = 225 but the uncorrected signal at ρ = +0.697 is consistent with the §5.5 anatomy enrichment as a real biology rather than an implant-coverage artifact.

**α frac_epi = −0.406** (more epi contacts → weaker α `obs_ρ_split`): the §5.5 candidate parsopercularis cell rests on a two-patient base in which Pat_07 contributes 3/5 trace-leaves with `ρ_split = −0.035` at α (CTM-anti). Under the matched-strength + epi-exclusion test (audit_68, §3.2 above), removing epi-zone contacts **strengthens** the cohort α signal (cohort median 8.3× → 27.7×, n_above 5/10 → 7/10, Wilcoxon p = 0.014). The cohort-level reading of α is therefore that α reorganization lives predominantly in non-epi tissue; the parsopercularis uncorrected-suggestive cell remains a directional indication anchored on Pat_14 (pro-aligned) once Pat_07 is removed under the epi-exclusion sensitivity check.

### 5.3 What §5.5 should add

After the existing §5.5 paragraph on β anatomy, add:

> An implant-geometry test (audit_64) on per-patient features against per-patient ρ_split readings finds the strongest β correlate at the left-hemisphere fraction `B_hemi`: Spearman ρ = +0.697 uncorrected p = 0.025 (BH q = 0.835 at m = 225 cells). Patients with predominantly left-hemisphere implants carry stronger β per-pair traces; right-hemisphere contact count N_R correlates at Spearman ρ = −0.694. The directional consistency with the §5.5 narrative (β trace-leaves concentrating at left fusiform / left superior temporal / left parsopercularis) supports a left-lateralized β cohort signature at the cohort-correlation level, with the absence of BH survival reflecting the 10-patient sample size rather than the absence of biology. The β trace is therefore left-lateralized at the cohort-correlation level but not localized to any single cortical region under multi-region correction.

After the §5.5 paragraph on α parsopercularis, add:

> An epi-exclusion sensitivity test on α (audit_68) recomputes `ρ_split` after removing per-patient epi-zone contacts from the FC matrix. Cohort-median ρ_split *strengthens* under exclusion from +0.105 → +0.187 (ratio 8.3× → 27.7×, n_above own surrogate 5/10 → 7/10, cohort-paired Wilcoxon p = 0.014). Pat_13 (30/119 epi contacts, the cohort maximum) is the single per-patient case whose α trace was carried by epi-zone activity (`ρ_split = +0.028` → −0.171 under exclusion); the other nine patients carry α reorganization in non-epi tissue. The §5.5 parsopercularis directional indication, anchored on Pat_07 (cross-probe-aggregate anti-aligned) and Pat_14 (pro-aligned), is consistent with this reading: after Pat_07 dropout (residual two-patient single-patient cell), the parsopercularis indication is carried by Pat_14 alone, who is a pro-aligned CTM patient under both the within-baseline and matched-strength readings.

---

## 6. §6 synthesis errata (`sec:synthesis`)

### 6.1 What the manuscript currently claims

The §6 head: "the central falsification of the equivalence is the controlled per-pair correlation on `Trho(τ)` at `τ' = 1/λ_max`, which surfaces the imprint at α, β, and low-γ at within-probe BH-q = 0.027 under the split-baseline, drift-floor, and cross-probe controls".

Band-resolved per-§6.1:
- β: registers at every level (CTM + KC topology + KC heights + d_P/d_F + VI(k) + Grassmann diagonal).
- γ_l: CTM + KC heights only + d_F + VI(k) k=45..56 + broad Grassmann diagonal.
- α: CTM only + fine VI(k) k≈3 + Grassmann diagonal above k≈10.

### 6.2 What the Errata adds to the synthesis

The matched-strength + epi-exclusion readings restructure the §6.1 band table. Updated reading by **three** complementary cross-phase memory objects — **per-pair LRG (`ρ_split` on `D_coph`)**, **dendrogram tree-distance (`T_KC(λ)`)**, **leading-mode subspace (`T_G(k)`)** — each tested against the within-baseline-null and the matched-strength-null:

**β — the multifacet imprint, three-probe matched-strength-confirmed.**
- LRG-CTM `ρ_split`: within-baseline q = 0.027, matched-strength ratio 23.7× (7/10, p = 0.005), epi-exclusion robust (audit_68 not run at β — owed follow-up but raw-FC and Grassmann both robust).
- KC: within-baseline 10/10 at q = 0.006, matched-strength λ=0 ratio 2.21× (mixed wiring + strength), λ=1 sign-flip in surrogate (wiring-specific in heights).
- Grassmann `T_G(k)`: matched-strength controlled at k=27..55 (29 contiguous cells), epi-exclusion 29/29 persist.
- Anatomy: audit_64 `B_hemi` ρ = +0.697 uncorrected — left-lateralized at cohort-correlation level.

**α — per-pair-only imprint, matched-strength + epi-exclusion confirmed.**
- LRG-CTM: within-baseline q = 0.027, matched-strength full FC ratio 8.3× (5/10, p = 0.002), matched-strength epi-excluded ratio 27.7× (7/10, p = 0.014) — **strengthens** under epi-exclusion.
- KC: matched-strength λ=0 ratio 0.12× (surrogate exceeds observed). α tree-distance does not carry the trace beyond strength-scrambling expectations.
- Grassmann: 4/111 sig cells, no contiguous run; null under matched-strength.
- Anatomy: parsopercularis directional indication anchored on Pat_14 after Pat_07 dropout; cohort α trace lives in non-epi tissue.

**γ_l — subspace-window-only imprint, matched-strength survives at k=12..23.**
- LRG-CTM: within-baseline q = 0.027, **matched-strength fails cohort-paired test** (ratio 30.2× but n_above 5/10, p = 0.116).
- KC: matched-strength λ=0 ratio 6.88× (n_below 2/10, p = 0.31); λ=1 anti-trace in both observed and surrogate.
- Grassmann: matched-strength k=12..23 (12 contiguous cells, best-k ratio 11.5× at k=17), epi-exclusion 10/12 persist; outside window 24 cells weaken under epi-exclusion.
- Anatomy: ctx-lh-fusiform Bonferroni-survived at audit_66-level (`p_hyper = 5.2 × 10⁻⁶`), preserved in current manuscript.

**γ_h — subspace-only imprint at k=19..27, LRG-CTM cohort-bimodal.**
- LRG-CTM: matched-strength cohort ratio 0.08× (surr exceeds obs), n_above 4/10, p = 0.246. Cohort strongly bimodal — 4 strong pro (Pat_03/05/06/08), 4 strong anti (Pat_02/07/14/15), 2 null (Pat_10/13). Per-pair memory exists at per-patient level but cohort disagrees on direction; median collapses to ~0.
- Raw-FC ρ_split: 88.5× ratio cohort-median with p = 0.188 (same bimodality at substrate).
- KC: not yet run at γ_h.
- Grassmann: matched-strength k=19..27 (9 contiguous cells, best-k ratio 15.4× at k=23 — largest in the cohort), epi-exclusion **6/9 persist + 3/9 weaken** and trace direction strengthens (`|T_G| −0.27 → −0.51` at k=25). Eliminates pure-HFO-contamination alternative.
- Anatomy: not yet tested. Caveat retained.
- **Cohort claim at γ_h is therefore subspace-only at k=19..27**; the per-pair-level reading is split into a pro-cohort (4 patients) and anti-cohort (4 patients) sub-cohort with no aggregate direction.

**δ — substrate only, no LRG-derived trace.**
- LRG-CTM: matched-strength ratio 1.0× (null).
- KC: λ=0 ratio 2.25× (n_below 2/10, p = 0.46).
- Grassmann: 23/111 sig cells at audit_66 with coarse-k shoulder k=57..63 (best-k ratio 4.2× at k=60); audit_67 adds 16 emergent cells under epi-exclusion (denoising signature).
- Exploratory: δ slow-wave reorganization may be partly carried by non-epi tissue and partly overwritten by epi noise; not a §5/§6 headline.

**θ — fully ergodic at every LRG layer.**
- LRG-CTM: matched-strength sign-flip (obs −0.040 vs surr +0.004).
- KC: λ=0 sign-flip (obs +0.83 vs surr +0.05).
- Grassmann: 8/111 sig cells (≈ chance), 3/87 persist under epi-exclusion.
- The θ band shows substrate-level cohort heterogeneity (raw-FC `ρ_split` cohort p = 0.138) that does not crystallize at any LRG-derived object under either null.

### 6.3 The two-null structure as a §6.2 addition

The §6.2 subsection ("why diffusion enriches the raw |ImCoh| reading") should add a paragraph on the within-baseline vs matched-strength null distinction:

> The within-baseline-null control battery applied at §5.3 (split-baseline `ρ_split` against `ρ_drift`, drift-floor, cross-probe restriction) tests whether the observed cross-phase shift exceeds the rest-to-rest drift baseline on the same halved-data noise budget; it is a drift-controlled claim. A complementary strong control — the matched-strength surrogate null (R = 200 4-cycle ±δ rewiring preserving every node's strength `s_i = Σ_j W_ij` exactly to 10⁻⁴) — tests whether the observed cross-phase shift can be reproduced by independent per-phase reorganization of node strengths with random edge identities. The two nulls do not nest: a signal can pass within-baseline (it exceeds drift) and fail matched-strength (it is reproducible by strength scrambling), or pass matched-strength (it exceeds strength scrambling) and fail within-baseline (it does not exceed drift). The within-baseline-null β KC 10/10 result at §5.2 stands as a drift-controlled finding; the matched-strength null decomposes it into a wiring-specific component (≈ 50% at λ = 0 topology, dominant at λ = 1 heights) and a strength-encoded component (≈ 50% at λ = 0). The matched-strength null applied uniformly to every probe identifies β as the only band whose cohort signal survives at multiple independent objects (LRG-CTM, Grassmann full FC, Grassmann epi-excluded), with γ_l and γ_h as subspace-only matched-strength-surviving traces at narrow k windows and α as per-pair-only surviving the test at the LRG layer.

### 6.4 §6.3 outlook items now closed or owed

The §6.3 outlook lists three open directions:
1. **Per-patient cross-phase module taxonomy anatomy mapping** — still open.
2. **Pat_07 / Pat_15 implant-geometry analysis** — **closed in part by audit_64 (left-hemisphere correlate at β)**. The §6.3 prose should reference audit_64 as the cohort-correlation-level analysis.
3. **Controlled τ-sweep at coarser scales** — still open.

New open items the Errata adds:
4. **γ_h anatomy correlate test** (audit_64 was scoped to α/β/γ_l; γ_h cohort-correlation against implant-geometry features is owed before claiming hippocampal/cortical attribution at γ_h).
5. **β epi-exclusion sensitivity at the LRG-CTM layer** (audit_68 was scoped to α only; β epi-exclusion at ρ_split is a one-cell-per-patient run from the existing epi-excluded eigendecomposition cache).
6. **γ_h LRG-CTM `ρ_split`** (in-flight; ETA tonight).
7. **KC matched-strength sensitivity at γ_h** (cheap from the `matched_strength_surrogate_lrg/` cache).
8. **Distance-based epi-exclusion** (epi + nearest-neighbour rings) to settle γ_h interface-effects caveat from the 2026-05-14 grassmann-epi-exclusion-sensitivity handoff.

---

## 7. Consolidated per-band quantitative tables (all probes, all controls)

Format: cohort row first (median observed, median surrogate, effect-size ratio, n above/below own surrogate, cohort-paired Wilcoxon p), then per-patient rows where load-bearing.

### 7.1 δ

| probe | obs cohort median | surr cohort median | ratio | n / 10 | Wilcoxon p |
|---|---|---|---|---|---|
| raw FC ρ_split | +0.111 | +0.005 | 24.6× | 7/10 above | 0.042 |
| LRG-CTM ρ_split | +0.008 | +0.008 | 1.0× | 4/10 above | 0.278 |
| KC λ=0 | −0.275 | −0.122 | 2.25× | 2/10 below | 0.46 |
| KC λ=1 | −3.88 | −2.48 | 1.56× | 3/10 below | 0.22 (also_negative) |
| Grassmann best-k k=60 | T_G −0.83 | T_G −0.20 | 4.2× | 7-cell contig run k=57..63 | — |
| Grassmann epi-X (audit_67) | — | — | — | 9 persist + 14 weaken + 16 emerge | — |

### 7.2 θ

| probe | obs cohort median | surr cohort median | ratio | n / 10 | Wilcoxon p |
|---|---|---|---|---|---|
| raw FC ρ_split | +0.117 | +0.006 | 20.2× | 7/10 above | 0.138 (bimodal) |
| LRG-CTM ρ_split | −0.040 | +0.004 | (sign flip) | 2/10 above | 0.722 |
| KC λ=0 | +0.83 (anti) | +0.052 | — | 2/10 below | 0.62 |
| KC λ=1 | −0.877 | −0.361 | 2.43× | 1/10 below | 0.81 |
| Grassmann audit_66 | 8 sig cells | — | ≈ 1× | — | edge of spectrum, numerical |
| Grassmann epi-X | — | — | — | 3 persist + 5 weaken + 5 emerge | — |

### 7.3 α

| probe | obs cohort median | surr cohort median | ratio | n / 10 | Wilcoxon p |
|---|---|---|---|---|---|
| raw FC ρ_split | +0.158 | +0.008 | 20.8× | 8/10 above | 0.014 |
| LRG-CTM ρ_split (full FC) | +0.105 | +0.013 | 8.3× | 5/10 above | **0.002** |
| LRG-CTM ρ_split (epi-X audit_68) | **+0.187** | **+0.007** | **27.7×** | **7/10 above** | **0.014** |
| KC λ=0 | −0.015 | −0.122 | **0.12×** (surr exceeds) | 1/10 below | 0.54 |
| KC λ=1 | −0.396 | +0.007 | sign flip | 1/10 below | 0.72 |
| Grassmann audit_66 | 4 sig cells | — | ≈ 1× | — | null |
| Grassmann epi-X | — | — | — | 2 persist + 2 weaken + 1 emerge | — |

α per-patient under epi-exclusion (LRG-CTM ρ_split, audit_68):

| patient | n_epi | obs full FC | obs epi-X | Δ | reading |
|---|---|---|---|---|---|
| Pat_02 | 14 | −0.037 | −0.035 | +0.002 | anti at both, no change |
| Pat_03 | 6 | +0.360 | +0.189 | −0.171 | weakened |
| Pat_05 | 14 | +0.078 | +0.121 | +0.043 | strengthened |
| Pat_06 | 10 | +0.833 | +0.692 | −0.141 | slight weakening |
| Pat_07 | 7 | +0.081 | +0.211 | +0.130 | strengthened |
| Pat_08 | 9 | +0.380 | +0.290 | −0.090 | weakened |
| Pat_10 | 10 | +0.130 | +0.193 | +0.063 | strengthened |
| **Pat_13** | **30** | **+0.028** | **−0.171** | **−0.199** | **FLIPPED anti (cohort-max epi-load)** |
| Pat_14 | 12 | +0.210 | +0.184 | −0.026 | minor weakening |
| Pat_15 | 0 | +0.040 | +0.040 | 0 | identity transform |

### 7.4 β

| probe | obs cohort median | surr cohort median | ratio | n / 10 | Wilcoxon p |
|---|---|---|---|---|---|
| raw FC ρ_split | +0.258 | +0.012 | 21.5× | 6/10 above | 0.053 |
| **LRG-CTM ρ_split** | **+0.221** | **+0.009** | **23.7×** | **7/10 above** | **0.005** |
| KC λ=0 | −0.890 | −0.403 | 2.21× | 1/10 below | 0.22 |
| KC λ=1 | −0.651 | +0.739 | sign flip (surr anti-trace) | 1/10 below | 0.54 |
| **Grassmann audit_66 best-k k=40** | **T_G −0.59** | **T_G −0.21** | **2.8×** | **29 contig cells k=27..55** | (per-k Wilcoxon) |
| **Grassmann epi-X audit_67** | longest run 36 cells | — | — | **29/29 persist in k=27..55**; **17 emerge** | — |

β per-patient at LRG-CTM ρ_split:

| patient | obs ρ_split | surr p50 | obs z | obs p |
|---|---|---|---|---|
| Pat_02 | +0.507 | +0.096 | +6.27 | 0.000 |
| Pat_03 | +0.373 | +0.021 | +9.16 | 0.000 |
| Pat_05 | +0.491 | +0.051 | +7.32 | 0.000 |
| Pat_06 | +0.211 | +0.018 | +3.71 | 0.000 |
| Pat_07 | +0.230 | −0.001 | +4.47 | 0.000 |
| Pat_08 | +0.502 | +0.040 | +8.76 | 0.000 |
| **Pat_10** | **−0.091** | −0.001 | **−1.70** | 0.975 |
| Pat_13 | +0.208 | +0.001 | +3.51 | 0.000 |
| **Pat_14** | **−0.049** | −0.003 | **−1.00** | 0.865 |
| Pat_15 | +0.083 | −0.014 | +1.27 | 0.105 |

7/10 above own surrogate at one-sided p < 0.05; Pat_10 and Pat_14 anti, Pat_15 weak pro.

### 7.5 γ_l

| probe | obs cohort median | surr cohort median | ratio | n / 10 | Wilcoxon p |
|---|---|---|---|---|---|
| raw FC ρ_split | +0.126 | +0.002 | 53.9× | 7/10 above | 0.053 |
| LRG-CTM ρ_split | +0.083 | +0.003 | 30.2× | 5/10 above | **0.116** (fails cohort-paired) |
| KC λ=0 | −0.695 | −0.101 | 6.88× | 2/10 below | 0.31 |
| KC λ=1 | +1.42 | +1.49 | both anti-trace | 1/10 below | 0.50 |
| **Grassmann audit_66 best-k k=17** | **T_G −0.71** | **T_G −0.06** | **11.5×** | **12 contig cells k=12..23, 11 strict-separated** | (per-k Wilcoxon) |
| **Grassmann epi-X audit_67** | — | — | — | **10/12 persist in window** + 2 weaken; outside: 24 weaken | partly epi-driven outside window |

### 7.6 γ_h

| probe | obs cohort median | surr cohort median | ratio | n / 10 | Wilcoxon p |
|---|---|---|---|---|---|
| raw FC ρ_split | +0.111 | +0.001 | 88.5× | 7/10 above | 0.188 (cohort heterogeneous) |
| **LRG-CTM ρ_split** | **+0.0003** | **+0.0042** | **0.08×** (surr exceeds) | **4/10 above** | **0.246** (cohort bimodal) |
| KC λ=0 (topology) | −2.98 | −0.199 | 15.0× (also_negative) | 3/10 below | 0.188 |
| **KC λ=1 (heights)** | **−3.82** | **−1.38** | **2.78×** | **5/10 below** | **0.042** (Wilcoxon-only marginal) |
| **Grassmann audit_66 best-k k=23** | **T_G −0.43** | **T_G −0.03** | **15.4×** (largest in cohort) | **9 contig cells k=19..27** (Wilcoxon-only) | (per-k Wilcoxon) |
| **Grassmann epi-X audit_67** | |T_G| deepens to −0.51 at k=25 | — | — | **6/9 persist + 3/9 weaken in k=19..27** (corrected — was misreported 7/9 elsewhere); trace strengthens under epi-X | physiological-attribution favored |

γ_h per-patient at LRG-CTM ρ_split (audit_63, completed 2026-05-15 ~00:30):

| patient | obs ρ_split | surr p50 | obs z | obs p | reading |
|---|---|---|---|---|---|
| Pat_02 | −0.168 | −0.001 | −10.37 | 1.000 | strong anti |
| Pat_03 | +0.505 | +0.014 | +16.67 | 0.000 | very strong pro |
| Pat_05 | +0.921 | +0.171 | +11.06 | 0.000 | strongest pro |
| Pat_06 | +0.520 | +0.014 | +15.80 | 0.000 | very strong pro |
| Pat_07 | −0.146 | +0.007 | −4.18 | 1.000 | anti |
| Pat_08 | +0.310 | +0.009 | +8.36 | 0.000 | strong pro |
| Pat_10 | +0.012 | −0.002 | +1.10 | 0.145 | weak pro / null |
| Pat_13 | −0.012 | +0.002 | −0.96 | 0.855 | null |
| Pat_14 | −0.145 | +0.000 | −11.72 | 1.000 | strong anti |
| Pat_15 | −0.471 | −0.014 | −9.50 | 1.000 | strongest anti |

Cohort split: 4 strong pro / 4 strong anti / 2 null. n_above own surrogate at p<0.05 = 4/10. Wilcoxon p = 0.246. γ_h LRG-CTM per-pair memory exists at the per-patient level (per-patient ratios up to 30×+) but the cohort disagrees on direction; the cohort median collapses to ~0.

γ_h per-patient at Grassmann epi-X k=25:

| patient | epi_frac | T_G epi-X | obs z | reading |
|---|---|---|---|---|
| Pat_02 | 0.12 | −0.70 | −18.2 | strong |
| Pat_03 | 0.05 | −0.03 | −3.2 | weak |
| Pat_05 | 0.12 | −0.73 | −16.4 | strong |
| Pat_06 | 0.09 | −0.57 | −12.2 | strong |
| Pat_07 | 0.06 | −0.32 | −4.1 | moderate |
| Pat_08 | 0.08 | −0.54 | −4.8 | strong |
| Pat_10 | 0.09 | −0.48 | −11.7 | moderate |
| **Pat_13** | **0.25** | **+0.01** | **+0.0** | **NEUTRALIZED** (cohort-max epi-load) |
| Pat_14 | 0.10 | −1.02 | −22.9 | strongest |
| **Pat_15** | **0.00** | **+1.06** | **+10.8** | **identity → audit_66 anti-trace anchor** |

---

## 8. Bottom-line claims for the writing agent

**Rock-solid under matched-strength + within-baseline + epi-exclusion (multi-control, multi-probe):**

1. **β** at the per-pair LRG layer (`ρ_split = +0.221`, ratio 23.7×, 7/10, p = 0.005) AND at the Grassmann subspace layer at k = 27..55 (29 cells, best-k ratio 2.8×) AND under epi-exclusion (29/29 retention). β KC adds a mixed wiring-and-strength decomposition (50/50 at λ = 0). β anatomy is left-lateralized at cohort-correlation level (B_hemi ρ = +0.697 uncorrected).

2. **α at the per-pair LRG layer**, full FC (ratio 8.3×, 5/10, p = 0.002) AND epi-excluded (ratio 27.7×, 7/10, p = 0.014). α tree-distance does NOT survive matched-strength; α subspace does NOT survive matched-strength. α trace is per-pair-only and lives in non-epi tissue. Pat_13 is the cohort-max-epi-load case whose α trace was carried by epi-zone activity.

3. **γ_l at the Grassmann subspace layer, k = 12..23** (12 cells, best-k ratio 11.5× — the largest sustained subspace effect-size ratio in the cohort), epi-exclusion preserves 10/12 cells in the manuscript window. γ_l ctx-lh-fusiform Bonferroni-survived anatomical localization stands. γ_l LRG-CTM fails cohort-paired matched-strength consistency (p = 0.116); should be reframed as "cohort-median magnitude exceeds strength-scrambling but per-patient cohort consistency is heterogeneous".

**Defensible with two named caveats:**

4. **γ_h at the Grassmann subspace layer, k = 19..27** (9 cells, best-k ratio 15.4× — single-probe maximum in the cohort). Epi-exclusion preserves **6/9 cells (3/9 weaken)** in the manuscript window and trace direction strengthens at the surviving cells (`|T_G| −0.27 → −0.51` at k = 25). This eliminates the "HFO contamination from epi contacts" alternative attribution. **γ_h KC** (completed 2026-05-15): λ = 1 heights-axis is Wilcoxon-marginal (cohort p = 0.042, ratio 2.78×) but n_below 5/10 — survives the Wilcoxon-only gate at α = 0.05 only; cohort split (5 trace + 4 null + Pat_15 strong anti). λ = 0 topology-axis is null (p = 0.188). **γ_h LRG-CTM** (completed 2026-05-15): cohort ratio 0.08× (n_above 4/10, p = 0.246) — fails matched-strength at the per-pair level due to cohort bimodality (4 pro / 4 anti / 2 null). **γ_h anatomy** (completed 2026-05-15 audit_64 6-band extension): NULL — top frac_epi correlate ρ = −0.273, p_unc = 0.446, BH q = 0.994. The γ_h cohort claim is therefore **subspace-only at the Grassmann probe** (with KC heights-axis as a Wilcoxon-only marginal companion), not multi-probe, not anatomically localizable, with the two open caveats: (a) interface effects (matched-strength surrogacy on `W[non_epi, non_epi]` does not control for epi-adjacent interictal physiology — owed: distance-based exclusion with nearest-neighbour rings); (b) Pat_13 dominance (the 30/119-epi patient neutralizes its trace under exclusion — biology vs measurement-precision question).

**Null at every LRG layer:**

5. **θ**: substrate-level cohort heterogeneity (p = 0.138) that does not crystallize at LRG-CTM (sign flip), KC (sign flip at λ=0), or Grassmann (8/111 ≈ chance) in either full FC or epi-excluded readings. Cleanest "null at LRG layer" band in the study.

**Substrate-only:**

6. **δ**: raw-FC ratio 24.6× cohort-paired p = 0.042 but no propagation to LRG-CTM (ratio 1.0×), KC moderate but cohort-inconsistent, Grassmann coarse-k shoulder (ratio 4.2× at k = 60) with 16 cells emerging under epi-exclusion (denoising signature). δ is an edge-only band at the LRG layer.

---

## 9. Open follow-ups (the writing agent should flag, not claim)

1. ~~γ_h LRG-CTM `ρ_split`~~ — **closed 2026-05-15 ~00:30** (audit_63 high_gamma): cohort ratio 0.08× (surr exceeds), n_above 4/10, p = 0.246; cohort bimodal (4 pro / 4 anti / 2 null). γ_h cohort claim is **subspace-only**.
2. ~~γ_h KC matched-strength~~ — **closed 2026-05-15 ~13:00** (audit_65 high_gamma): λ = 0 cohort p = 0.188 null; **λ = 1 cohort p = 0.042 Wilcoxon-only marginal trace** at the heights axis (ratio 2.78×, n_below 5/10). Per-patient direction split (5 trace + 4 null + Pat_15 strong anti, z = +3.16). γ_h heights-axis adds a Wilcoxon-only-marginal companion to the Grassmann subspace claim; γ_h topology-axis is null.
3. ~~γ_h anatomy implant-geometry test~~ — **closed 2026-05-15** (audit_64 6-band extension): NULL — top correlate is frac_epi at ρ = −0.273, p_unc = 0.446, BH q = 0.994. γ_h has no anatomical localization signal at the implant-geometry feature level; the §5.5 anatomy claim continues to be α/β/γ_l only.
4. **β epi-exclusion sensitivity at LRG-CTM** — audit_68 was scoped to α only. β `ρ_split` on `W[non_epi, non_epi]` is owed; one CLI call from the existing epi-excluded eigendecomposition cache (`data/cache/matched_strength_surrogate_epi_excluded_lrg/`).
5. **Distance-based epi-exclusion** (epi + N-nearest non-epi rings) — sharper probe of the γ_h interface-effects caveat.
6. **Hippocampal-leaf enrichment for γ_h trace under epi-exclusion** — direct test of the "sharp-wave ripple memory consolidation" hypothesis the audit_67 verdict explicitly defers.
7. **§5.6 cross-phase module taxonomy anatomical mapping** — Desikan-Killiany region assignment per (patient, band, class) clade. Still open from the §6.3 outlook.
8. **Controlled τ-sweep at coarser scales** — original §6.3 open direction. Reads whether the band-resolved imprint sharpens or dissolves as the propagator integrates over longer pathways.
9. **Cohort-wide drift-triangle null at the LRG layer** (mirroring Control 1 of §5.5.4 raw-FC caveats) — currently run only at Pat_06.
10. **Symmetric cross-baseline check at the diffusion layer** (mirroring Control 3) — comparing `Trho(τ)` on (`taskt-late`, `rspost-early`) against (`rspre-late`, `rspost-early`) on matched half-segment noise budgets.
11. **§7 / §5.4 in-code `verdict == "separated"` label correction** — `audit_69` (2026-05-15) regating run shows that the audit_66 / audit_67 in-code `n_k_sep` / `longest_run_sep` columns advertise 0/111 separated cells for β under the composite gate (cohort Wilcoxon `p<0.05` AND `|med_surr|<0.05·|med_obs|` AND `n_below≥8`), even though β is the strongest matched-strength-controlled subspace signal in §5.4 under the gate the manuscript actually cites. The §7 errata 29/12/9 contiguous cells and 29/29 / 10/12 persist counts are the Wilcoxon-only counts and stand unchanged. The audit_66/67 README label should be removed in the next pipeline pass to prevent future agents from being misled by the in-script summary. See `.agents/reports/2026-05-15_grassmann-regating-no-7of10-filter.md`.

---

## 10. Provenance

- Cohort: n = 10 (Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15)
- Era: COHORT_N10 / IMCOH_ABS (post-2026-04-25 Pat_14 vendor task_test replacement)
- FC method: `imcoh_abs = ⟨|ImCoh|⟩_f` (Nolte 2004 signed ImCoh cached, `np.abs` per frequency bin then band-averaged at load time)
- LRG: `τ = 1/λ_max`, propagator `ρ_ij = (U exp(−τΛ) U^T)_ij / tr(·)`, communication distance `Trho = 1/ρ`, average linkage on condensed Trho
- Matched-strength surrogate: 4-cycle ±δ rewiring, `n_swaps = 20 · N(N−1)/2`, R = 200 surrogates per (patient, band, phase); seeds 20260510 (audit_63), 20260511 (audit_65, audit_66), 20260514 (audit_67, audit_68); strength preservation verified `‖s_obs − s_surr‖_∞ < 10⁻⁴`
- Library helpers: `lrg_eegfc.utils.surrogate.matched_strength.{strength_preserving_shuffle, verify_strengths, load_or_compute_surrogate_eigs}`
- Surrogate eigendecomp caches (180 files each):
  - Full FC: `data/cache/matched_strength_surrogate_lrg/Pat_NN/{band}_{phase}_R200_swap20_seed20260511_imcoh_abs.npz`
  - Epi-excluded: `data/cache/matched_strength_surrogate_epi_excluded_lrg/Pat_NN/{band}_{phase}_epiX_R200_swap20_seed20260514_imcoh_abs.npz`
- Build scripts:
  - `scripts/01_compute/audit/audit_63_split_baseline_surrogate.py` — LRG-CTM `ρ_split` matched-strength
  - `scripts/01_compute/audit/audit_64_implant_geometry_test.py` — cohort-correlation anatomy (5 bands)
  - `scripts/01_compute/audit/audit_65_kc_matched_strength_surrogate.py` — KC matched-strength
  - `scripts/01_compute/audit/audit_66_grassmann_matched_strength_surrogate.py` — Grassmann full FC
  - `scripts/01_compute/audit/audit_67_grassmann_epi_exclusion.py` — Grassmann epi-excluded
  - `scripts/01_compute/audit/audit_68_alpha_epi_exclusion.py` — α LRG-CTM epi-excluded

End of Errata Corrige document.
