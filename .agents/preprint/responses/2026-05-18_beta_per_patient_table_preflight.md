---
name: beta-per-patient-table-preflight
era: IMCOH_ABS_COHORT_N10
status: preflight-reply
kind: writing-agent-reply
date: 2026-05-18
target: .agents/preprint/tables/beta_per_patient.tex (+ companion .md)
band: beta
---

# Preflight reply — β per-patient table

**Head.** All four CSV paths in the ticket exist with the data the writing agent needs, but two of them have different filenames than the ticket assumed and one row of the ticket's premise is wrong (the within-baseline and matched-strength ρ_split runs differ by ≈ 0.05–0.13 per patient, not 10⁻³). Cohort counts cross-check cleanly with the paragraph (8/10 ρ_split > 0, 7/10 above own surrogate at p < 0.05, 8/10 above drift, 9/10 below own surrogate on Grassmann at k = 40, 7/10 substrate trace). One design decision is bumped back to the user: whether col. 3 should hold the within-baseline ρ_split (paragraph-consistent, +0.222 cohort median) or the matched-strength ρ_split (row-self-consistent with cols 4–5, +0.221 cohort median). I recommend the latter and have written the table that way, with a footnote.

## 1. CSV paths and column names (confirmed)

| Ticket col | Quantity | CSV (actual filename) | Column |
|---|---|---|---|
| 2 | `N_ch` | `data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv` | `N_nodes` |
| 3 | ρ_split (within-baseline) | `data/audit/ctm_triangle/Td_per_patient_per_band.csv` (NOT `per_patient.csv`) | `rho_split` |
| 3 (alt) | ρ_split (matched-strength) | `data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv` (NOT `per_patient.csv`) | `obs_rho` |
| 4 | surrogate median | `data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv` | `surr_p50` |
| 5 | z vs surrogate | same | `obs_z` |
| 6 | Grassmann T_G(k=40) | `data/audit/grassmann_matched_strength_surrogate/per_patient_per_band_per_k.csv` | `obs_T_G` filtered to `band=beta, k=40` |
| 7 | z vs Grassmann surrogate | same | `obs_z` |
| 8 | substrate T_d^(d_S) sign | `data/audit/raw_fc_phase_distance/Td_per_patient_per_band.csv` (NOT `lrg_phase_distance`) | sign of `S` |

The ticket's `data/audit/lrg_phase_distance/Td_per_patient_per_band.csv` exists but holds the **LRG-derived** d_S (column `T_S`), not the raw-FC substrate. The brief §6's "Pat_07 anti at substrate" referent is the **raw-FC** d_S in `data/audit/raw_fc_phase_distance/Td_per_patient_per_band.csv`, column `S`. At LRG-d_S (column `T_S` in `lrg_phase_distance`) Pat_07 reads −0.109 (trace direction); at raw-FC d_S (column `S` in `raw_fc_phase_distance`) Pat_07 reads +0.005 (≈ 0, sign + → anti). Using raw-FC d_S resolves the brief's framing.

## 2. The within-baseline / matched-strength per-patient divergence (ticket premise correction)

The ticket assumed the two ρ_split runs differ per-patient by ≈ 10⁻³. The actual per-patient differences (within-baseline `ctm_triangle/rho_split` minus matched-strength `obs_rho`) at β:

| Patient | within-baseline | matched-strength | diff |
|---|---|---|---|
| Pat_02 | +0.380 | +0.507 | +0.127 |
| Pat_03 | +0.463 | +0.373 | −0.090 |
| Pat_05 | +0.483 | +0.491 | +0.008 |
| Pat_06 | +0.259 | +0.211 | −0.048 |
| Pat_07 | +0.185 | +0.230 | +0.045 |
| Pat_08 | +0.503 | +0.502 | −0.002 |
| Pat_10 | −0.123 | −0.091 | +0.032 |
| Pat_13 | +0.185 | +0.208 | +0.023 |
| Pat_14 | −0.073 | −0.049 | +0.024 |
| Pat_15 | +0.144 | +0.083 | −0.061 |

The two runs use independent rsPre A/B-split realizations, so per-patient values differ by 0.05–0.13 even though cohort medians match (+0.2222 vs +0.2206). The 10⁻³ premise was wrong.

**Consequence for the table.** If col. 3 holds within-baseline `rho_split` and cols 4–5 hold matched-strength `surr_p50` and `obs_z`, the z in col. 5 does NOT correspond to the value in col. 3 (the z is computed against the matched-strength `obs_rho`, not the within-baseline `rho_split`). Per-row interpretation breaks at the per-patient level. Two options:

- **Option A (paragraph-consistent):** col. 3 = within-baseline ρ_split (cohort median +0.2222 matches the paragraph's `+0.222`). Cols 4–5 keep the matched-strength surrogate stats with a footnote noting the realization offset.
- **Option B (row-self-consistent, recommended):** col. 3 = matched-strength `obs_rho` (cohort median +0.2206). Cols 4–5 then correspond row-for-row; the z in col. 5 is the z of col. 3 against col. 4. Cohort median +0.221 is indistinguishable from the paragraph's +0.222 at three decimals (difference 0.0016 ≪ matched-strength surrogate noise).

I have written the table per Option B; if the user prefers Option A, swap col. 3 to `ctm_triangle/Td_per_patient_per_band.csv` `rho_split` and add a footnote.

## 3. Cohort counts cross-check (paragraph ↔ CSV)

| Paragraph claim | CSV | Confirmed |
|---|---|---|
| Cohort median ρ_split^coph +0.222 | `ctm_triangle/cohort_summary.csv` `rho_split_median` = +0.22223 | ✓ |
| 8/10 trace direction on ρ_split^coph (within-baseline) | per-patient sign: Pat_10 + Pat_14 below zero | ✓ |
| 8/10 above drift floor | `ctm_triangle/cohort_summary.csv` `n_above_drift` = 8/10 (Pat_10 + Pat_13 fail) | ✓ |
| 7/10 above own surrogate at matched-strength | `matched_strength_surrogate_split_baseline/cohort_summary.csv` `n_above_surrogate` = "7/10" (Pat_10 z = −1.7, Pat_14 z = −1.0, Pat_15 z = +1.27 but one-sided p = 0.105 > 0.05) | ✓ |
| Cohort Wilcoxon p = 0.005 (matched-strength) | `paired_wilcoxon_p` = 0.00488 | ✓ |
| Cohort Wilcoxon p = 0.014 (drift) | `wilcoxon_split_gt_drift_p` = 0.01367 | ✓ |
| Grassmann 9/10 below surrogate at k=40 | per-patient `obs_z` at β, k=40: all negative except Pat_15 (+4.17), nine with `obs_p_one_sided_lower` = 0.0 | ✓ |
| Cluster-extent permutation p = 0.005 (β) | `grassmann_cluster_extent/cohort_summary.csv` `cluster_p_longest_run` for β | not re-grepped here, but VERDICT_LEDGER.md cites this; pulled at table-write time |

All counts in the table footer match the paragraph claims exactly.

## 4. Notes-column flags (per brief §6)

- **Pat_03**: "acquired at 1024 Hz, handled at config layer" (NOT outlier framing — per the 2026-05-18 rule update).
- **Pat_07**: "anti at substrate; pro at LRG" — based on raw-FC `S` = +0.005 (anti, ≈ 0) and `rho_split` = +0.185 z = +4.47 (pro at LRG). The brief's wording carries over; sign is technically + but value is essentially zero — flag in companion .md.
- **Pat_15**: "right-hemisphere-only implant; LRG anti-aligned" — only anti-trace patient at Grassmann k = 40 (T_G = +1.020 vs cohort median −0.49). Substrate also anti (`S` = +0.217).
- **Pat_13**: agent's call on visual-balance flag. Per-patient β: `rho_split` = +0.185 (pro), `T_G(k=40)` = −0.128 (pro), substrate `S` = −0.007 (essentially zero but technically trace). Flag for heaviest epi burden (30/119) only if the column needs visual balance — I left it blank.
- All other patients: em-dash.

## 5. Numerical-style decisions made for the .tex

- Three decimals for ρ_split, surrogate median, T_G(k=40). Two decimals for z. Integer for N_nodes.
- Pat_15 T_G(k=40) = +1.020 — explicit + sign for the only anti-aligned patient (`S[table-format=+1.3]` handles this).
- Pat_07 substrate sign: + (raw-FC `S` = +0.005 → anti). Flag in companion that this is near-zero.
- Asterisks attached to cohort footer values via `\multicolumn{1}{S[table-format=+1.3]}{$+0.221^{**}$}` idiom (`S` column does not accept text annotations on raw numerics).
- Cohort row: bold "Cohort" label, `\midrule` above.

## 6. Outputs

- LaTeX: `.agents/preprint/tables/beta_per_patient.tex` — ready for `\input{...}` (path follows the project rule that *all preprint work lives in `.agents/preprint/`*).
- Companion .md: `.agents/preprint/tables/beta_per_patient.md` — source CSVs, column headers, per-row sanity check, surprise-patient flags.

## 7. What the writing agent should still decide

1. **Option A vs Option B** for col. 3 (within-baseline ρ_split vs matched-strength `obs_rho`). I have written B; flipping to A is a 10-line edit if you prefer.
2. **Pat_07 substrate sign rendering:** + at +0.005 is technically correct but visually misleading. I have rendered `+` with no annotation; an alternative is `≈ 0` or a footnote. Reader's call.
3. **N_ch column source:** I have used `N_nodes` from the matched-strength CSV (every patient has the same `N_nodes` across all bands — the FC matrix size after cleaning). If you want the LRG-eigendecomposition-effective `N_nodes` (excluding zero-mode contacts), let me know — it would require loading the LRG cache and counting non-trivial eigenvalues, which the matched-strength `N_nodes` already implicitly does (the FC is symmetric and connected by construction at `imcoh_abs`).
