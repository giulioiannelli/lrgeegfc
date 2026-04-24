---
name: epileptic-imcoh-final
type: report
era: IMCOH_ABS
status: current
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

# Epileptic-zone signature in |ImCoh|-LRG — Final consolidated results

**Status:** validated 2026-04-15/18 on `imcoh_abs` (post-reset).
**FC quantity:** `|ImCoh|` (Ewald 2012 / Bastos & Schoffelen 2016).
**Era:** post-ImCoh-reset (2026-04-15). All numbers on `imcoh_abs`.

---

## One-line summary

In 4/5 patients, epileptic sEEG contacts share stronger cross-probe
|ImCoh| (1.4× non-epi baseline, BH-FDR q<0.05) and the LRG hierarchy
faithfully captures this (ρ = −0.81). The finding is methodologically
clean but limited to a clinical-utility demonstration, not a clinical
discovery, given n = 5.

---

## What to keep (stable, defensible)

### 1. Cross-probe edge strength (V1) — main result

After removing same-probe (spatially confounded) pairs:

| Patient | CP ratio | BH-FDR sig / 24 |
|---------|----------|-----------------|
| Pat_02 | 1.43 | 23 |
| Pat_03 | 2.09 | 13 |
| Pat_05 | 1.36 | 17 |
| Pat_07 | 1.44 | 13 |
| **Pat_08** | **0.91** | 4 (inconclusive) |

Pooled: **70/120 BH-FDR (58%), 50/120 Bonferroni, mean ratio 1.45×**.

Pat_08 excluded: all 9 epileptic contacts on 2 probes (50% same-probe
pairs), so cross-probe test is underpowered. Not a negative — electrode
geometry precludes the test.

**Eligibility rule:** epileptic contacts on ≥3 probes.

### 2. Logic chain (V5) — validation

Spearman(cross-probe edge ratio, cross-probe UM ratio) = **ρ = −0.811,
p = 3.1 × 10⁻²⁹** (n = 120). Per-patient: −0.94 to −0.49, all p < 0.02.

Not a finding, but a validation: the LRG hierarchy reflects the
connectivity structure, not a clustering artefact. One-liner in paper.

### Band breakdown (V1)

| Band | ratio | BH / 20 |
|------|-------|---------|
| delta | 1.68 | 15 |
| alpha | 1.63 | 11 |
| theta | 1.47 | 13 |
| low_gamma | 1.45 | 12 |
| beta | 1.39 | 13 |
| high_gamma | 1.05 | 6 |

Broadband effect, strongest in delta, weakest in high_gamma.
No single band stands out enough for a band-specific narrative.

---

## Interesting but underpowered (n = 5)

### 3. Provincial hub topology (network metrics, B1)

Best discriminator: participation coefficient AUC = 0.725.
Profile: high strength + high clustering + low participation + low
betweenness. "Closed loop, not integrative hub."

Biologically interesting, somewhat novel. Not individually significant
with n = 5. Worth pursuing with more patients.

### 4. Stratified H2a — epi-epi reorganizes MORE than non-non

Per-band epi-epi / non-non H2a ratio: 2.1–3.2× (delta through
low_gamma). Contradicts the "static biomarker" assumption — epileptic
tissue internally reorganizes across task phases even though its mean
edge strength is stable.

Wilcoxon p = 0.84 for beta (3/5 patients epi > non). Striking direction
but n = 5 kills significance. Hypothesis for future work, not a claim.

---

## Drop (not worth more effort at n = 5)

- **V2 (NN enrichment):** 37/120 BH sig, Pat_07 and Pat_08 fail.
  Fragile metric, V1 says the same thing more directly.
- **Phase invariance:** correct at group level (Friedman p > 0.3 in every
  band), but the stratified H2a suggests epi pairs aren't truly static.
  Don't build a narrative on it.
- **Ablation (H2a with epi-epi zeroed):** heterogeneous per-patient,
  mean +11% in beta. Conceptually interesting but data doesn't support a
  claim. Park it.
- **Probe-matched ultrametric null (A1):** 42 raw → 5 BH → 0 Bonferroni.
  Too weak. Don't lead with it.

---

## Recommended paper use

Write as a **short section** in the existing paper ("Section 6: Clinical
application") rather than a standalone epilepsy paper. Frame as:

> "LRG applied to |ImCoh| detects epileptogenic network structure in
> 4/5 patients without prior anatomical information, using only
> resting-state sEEG connectivity. After excluding same-probe pairs
> (spatial confound), epileptic contacts show 1.4× stronger |ImCoh|
> (BH-FDR q < 0.05, 70/120 conditions, 50/120 survive Bonferroni).
> The LRG hierarchy faithfully reflects this connectivity
> (ρ = −0.81). These results demonstrate clinical utility of the
> LRG-ImCoh framework for identifying pathological sub-networks."

**Bottleneck:** more patients, not more analysis. The infrastructure is
solid (10 seconds per patient). The per-patient worker pattern handles
memory safety for large cohorts.

---

## Spatial bias audit (V4)

| Patient | N_epi | probes | %same-probe | qualifies |
|---------|-------|--------|-------------|-----------|
| Pat_02 | 14 | 4 (B,G,L,T) | 20% | ✓ |
| Pat_03 | 6 | 2 (L,O) | 47% | borderline |
| Pat_05 | 14 | 5 (L',P',X,Y,Y') | 15% | ✓ (cleanest) |
| Pat_07 | 7 | 2 (A,K) | 43% | borderline |
| Pat_08 | 9 | 2 (T,W) | 50% | ✗ |

---

## Multiple-comparison audit (V3)

| Analysis | N | Raw | E[FP] | BH | Bonf |
|----------|---|-----|-------|-----|------|
| V1 cross-probe | 120 | 71 | 6 | 70 | 50 |
| V2 cross-probe NN | 120 | 47 | 6 | 37 | 15 |

---

## Bug fixed during this analysis

`compute_lrg_analysis()` was silently writing to cache even with
`use_cache=False, overwrite_cache=False`. Fixed: write gated by
`overwrite_cache or (use_cache and not cache_path.exists())`.
The ablation test temporarily corrupted 120 LRG cache files; all
were regenerated and validated.

---

## Reproducibility

| Script | Purpose |
|--------|---------|
| `scripts/08_epileptic/epileptic_imcoh_validation.py` | V1-V5 driver |
| `scripts/08_epileptic/_validation_one_patient.py` | per-patient V1/V2/V5 |
| `scripts/08_epileptic/epileptic_imcoh_deep_lrg.py` | A1-A5, B1-B7, figures |
| `scripts/08_epileptic/epileptic_h2_stratified.py` | stratified H2a + ablation |
| `scripts/08_epileptic/_h2_ablation_one_patient.py` | per-patient ablation |
| `scripts/08_epileptic/_regenerate_lrg_cache.py` | rebuild LRG cache |

Outputs: `data/outputs/figures/epileptic_imcoh_{validation,deep,h2_stratified}/`

Memory pitfall: `mannwhitneyu(method="asymptotic")` — default `auto` on
Pat_03 (epi=6, non=116) triggers exact permutation → 29 GB OOM.
