---
name: preprint-anatomy-control-battery
era: IMCOH_ABS_COHORT_N10
status: locked_2026-05-19
kind: anatomy-control-battery-lockdown
companion: ANATOMY_LEDGER.md
parent: CONTROLS.md
---

# Anatomy control battery for the LRG trace probes (locked 2026-05-19)

**Head.** The trace-detection battery (`CONTROLS.md`, 5 controls) is silent on
**where** the trace lives anatomically. This document locks the **anatomy
control battery**: four controls A1–A4 that constitute the entire anatomy-side
lockdown for the preprint. Verdicts are tagged `strong localized` / `weak
localized` / `not localized` per (band, probe), recorded in
`ANATOMY_LEDGER.md`. As with the trace battery, matched-strength is mandatory:
every anatomy claim is matched-strength-controlled or it does not earn a
verdict.

## What this battery does not change

The trace verdict (`VERDICT_LEDGER.md`) is independent of the anatomy verdict.
A band can be `strong trace, both probes` and `not localized` (the trace
exists but distributes uniformly across DK regions); a band can also be `weak
trace, only Grassmann` and `strong localized` (the trace is weak but
concentrates sharply in a small region set). The anatomy ledger only addresses
**where** an existing trace sits — never whether the trace exists.

## The four primary anatomy controls

### A1 — Hypergeometric per-region enrichment (Fisher exact + BH-FDR)
For each DK region `r`, count `k_r` = number of "trace-flagged" leaves/pairs
that intersect `r`, vs the cohort-wide marginals. Fisher exact one-sided
greater per region, BH-FDR over `R = N_DK_regions_with_coverage` per band.
- **Trace-flagged unit**:
  - **Cophenet probe**: top-decile per-pair `Δ(ρ_split^coph)` pairs (the pairs whose split-cophenet correlation contributes the most to the cohort signal).
  - **Grassmann probe**: top-decile per-node `participation` averaged across `S(b) = {k : p_k(b) < α_k}` — the **support of the cluster-mass statistic `T_G^*`**, i.e. the set of k-cells whose per-k cohort-paired Wilcoxon clears `α_k`. Under the locked all-clusters `T_G^*` paradigm there is no canonical "significant window"; the retired `K*(b)` (longest contiguous-significant window: β k=27..55, γ_l k=12..23, δ-full k=57..63, δ-epi-X k=33..39) is replaced by `S(b)` everywhere. The per-node statistic is the **unweighted average** `part(i; b) = (1/|S(b)|) · Σ_{k ∈ S(b)} ||U_k^T e_i||²_2`.
- **Gate**: `q_BH < 0.05` per region.
- **Output**: per-band per-probe ranked enrichment table; cohort-aggregated `n_pat_trace ≥ 6/10` for cohort-robust enrichment (not gate; descriptive).

### A2 — Sampling-corrected enrichment
sEEG implants are clinically motivated; DK regions have wildly different
contact counts (e.g., hippocampus heavily sampled in clinical practice; cortex
sparsely sampled). A2 corrects A1 by conditioning on the cohort-wide contact
distribution per DK region, computing the conditional enrichment.
- **Method**: bootstrap N=1000 over patient resampling within each DK region; compute the median enrichment per region; report 95% CI.
- **Gate**: A1-significant region's 95% CI lower bound > 1.0 (i.e., the enrichment is robust to sampling).
- **Purpose**: rejects enrichment claims that exist only because a region is sampled in many patients (e.g., hippocampus would trivially appear enriched without A2).

### A3 — Matched-strength surrogate (mandatory)
The A1/A2 enrichment pipeline is rerun on each of R=200 matched-strength
surrogate adjacencies (4-cycle ±δ rewiring, seed 20260511 — same as C3 for
the trace battery). For each region `r`:
- `z_r = (enrichment_obs_r − mean_surr_enrichment_r) / sd_surr_enrichment_r`
- Empirical p-value: `p_r = (1 + #(surr_enrichment_r ≥ obs_enrichment_r)) / (R + 1)`.
- **Gate**: `p_r < 0.05` AND `z_r > 2.0`.
- **Reason for mandate**: matched-strength preserves node strength but reshuffles topology; an enrichment claim that survives this is genuinely about topological organization, not about which contacts are strong on their own.

### A4 — Implant-geometry cohort regression
Each contact `i` has coordinates `(x_i, y_i, z_i)` in MNI / patient-native
space and a per-DK-region label. For each band×probe, regress per-contact
"trace flag" (boolean: was this contact in a top-decile trace pair / a
high-participation mode in `U_k`) against:
- DK region one-hot (already covered by A1/A2)
- Distance to nearest epileptogenic-zone contact (cohort-wide)
- Geometric centroid distance from the patient implant's center of mass
- Hemisphere (L/R)

Reports which **geometric covariates independently predict trace flag**,
independent of DK label. This protects against an enrichment that is
actually about being-close-to-the-epi-zone (which DK labels alone cannot
detect) or about implant-center proximity.
- **Output**: per-band per-probe regression table; flagged geometric covariates with `q_BH < 0.05`.
- **Not a gate**: A4 is a sensitivity layer flagging when an A1+A2+A3-passing enrichment is potentially attributable to geometry rather than to DK region per se.

## Verdict vocabulary (locked)

### Per-(band, probe)
- **strong localized** — A1+A2+A3 all pass for at least one DK region. A4 confirms the localization is not purely geometric (no geometric covariate alone explains the trace flag at q_BH<0.05).
- **weak localized** — A1+A2 pass but A3 fails for all regions (cohort enrichment exists but is not robust to matched-strength rewiring). Or A1+A2+A3 pass for one region but A4 attributes >50% of the trace-flag variance to a geometric covariate alone.
- **not localized** — A1 fails for all regions, or only trivial regions (Wm = white matter; Unk = uncertain) enrich. The trace exists in the network but does not concentrate anatomically.

### Per-band aggregation
- **strong localized, both probes** (β only candidate — has both trace probes)
- **strong localized, only D_coph** (α or β with cophenet localized but Grassmann not)
- **strong localized, only Grassmann** (γ_l or δ with Grassmann localized)
- **mixed** (one probe strong, other weak)
- **weak localized, ...** (per-probe combinations)
- **not localized** (no probe shows DK enrichment robust to A1+A2+A3)

## Manuscript ↔ lab label mapping (locked 2026-05-19 pm)

The manuscript methods §sssec:methods_anatomy uses a compressed label set
because the deferred controls (lab A2 + lab A4) are **not disclosed in the
main text**. The mapping is:

| Manuscript label | Lab label | Test |
|---|---|---|
| **A1** | A1 | Hypergeometric per-region enrichment (Fisher exact + BH-FDR) |
| **A2** | **A3** | Matched-strength surrogate enrichment (R=200, 4-cycle ±δ) |
| (not in manuscript) | A2 | Sampling-corrected bootstrap (deferred to sensitivity supplement) |
| (not in manuscript) | A4 | Implant-geometry cohort regression (deferred to sensitivity supplement) |

All lab-side documentation (this file, `ANATOMY_LEDGER.md`, audit CSV column
names like `passes_A3`, per-band briefs `01_beta.md`–`06_delta.md`) **keeps
the lab labels A1 / A3** to preserve traceability with the locked battery.
Only the manuscript LaTeX uses A1 / A2; the writing-agent directives
(`writing_directive_2026-05-19_methods_anatomy.md`,
`writing_directive_2026-05-19_anatomy_clusterext_rerun.md`) carry both
conventions explicitly. The deferred controls retain their original numbers
(lab A2, lab A4) in this battery for when the sensitivity supplement is
written.

## Anti-revisitation

A verdict in `ANATOMY_LEDGER.md` cannot be changed by reading a different CSV.
Re-evaluation requires:
1. A new dated audit producing fresh `data/audit/<anatomy_probe>/cohort_summary.csv`.
2. A dated revision entry in `ANATOMY_LEDGER.md`.
3. Cascading updates in per-band briefs.

## What this battery explicitly is NOT

- **NOT a mode-by-mode anatomy of every eigenmode** — that's intractable
  (`N ≈ 112–121` modes per patient). A3 controls aggregate participation
  across the band's significant `k`-window.
- **NOT a per-pair direction analysis** — direction of `D_coph` change is
  not part of the anatomy claim; only that the pair sits in the top-decile.
- **NOT an attempt to explain *why* a region enriches** (e.g., causal /
  mechanistic / connectivity-based reasoning). That is interpretive content
  for the writing agent, not part of the lockdown.

## Object-to-control mapping per band's trace probe(s)

| Band | Trace verdict | Anatomy audit(s) | A3 mandatory |
|---|---|---|---|
| **β** | strong trace, both probes | cophenet (top-decile pairs) + Grassmann (top participation over `S(β)`) | yes for both |
| α | strong trace, only D_coph | cophenet (top-decile pairs; full + C5 epi-X subset) | yes |
| γ_l | strong trace, only Grassmann ↑ | Grassmann (top participation over `S(γ_l)`) | yes |
| δ | weak trace, only Grassmann (LOO Pat_08 fails Decision-12 precondition) | Grassmann (top participation over `S(δ)` full + `S(δ)` epi-X) | yes for both windows |
| θ | no trace | — | n/a |
| γ_h | no trace | — | n/a |

Total: 5 distinct audit runs (β cophenet, β Grassmann, α cophenet, γ_l Grassmann, δ Grassmann × 2 ensembles). All Grassmann runs aggregate participation over `S(b) = {k : p_k(b) < α_k}`, the support of `T_G^*`, **not** over the retired longest-contiguous-significant window. **Pending re-run with the corrected aggregation** (audit_72 currently outputs verdicts under the retired window — verdicts in `ANATOMY_LEDGER.md` are anchored to the retired CSV and must be regenerated; see revision-history note below).

## Reuse of existing artifacts (descriptive only)

- `data/audit/lrg_localization_anatomy/region_enrichment.csv` — **KC-era**, retired by `VERDICT_LEDGER.md` 2026-05-18; **not citable** for preprint claims. Useful only as a sanity check that the new audits land in the same region neighborhood for β.
- `data/audit/anchor_anatomy_baseline/cohort_aggregate.csv` — older δ anchor anatomy; reproduces the known δ cross-probe 1.55× pattern. Citable for the **descriptive δ anchor anatomy note** only (LEDGER Decision 5 framing), not as a trace localization claim.
- `data/audit/epi_fraction_anatomy/region_epi_fractions.csv` — per-contact epi-zone fractions; **citable** as A4 covariate (distance-to-epi).
- `data/audit/implant_geometry/` — implant coordinates; **citable** as A4 covariate.

## What earns a citation in the preprint

Only verdicts in `ANATOMY_LEDGER.md` are citable as trace localization claims.
KC-era anatomy memory (`result_2_lrg_beta_trace.md` Hippocampus + left
fusiform) is **archaeology** — replaced by the new audits' verdicts once Phase
F lands.
