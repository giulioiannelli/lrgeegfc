---
name: beta-per-patient-table-companion
era: IMCOH_ABS_COHORT_N10
status: current
kind: table-companion
date: 2026-05-18
band: beta
target: tab:beta_per_patient
preflight: .agents/preprint/2026-05-18_beta_per_patient_table_preflight.md
location: .agents/preprint/tables/ (all preprint work lives under .agents/preprint/ per project rule)
---

# β per-patient table — companion notes

This file documents the data sources, columns, sanity checks, and
surprise-patient flags for `beta_per_patient.tex`. Every number in
the table traces to one of the CSVs listed below at the exact column
named.

## Source CSVs (by table column)

| Table col | Quantity | CSV path | Column name |
|---|---|---|---|
| 2 | `N_ch` | `data/audit/matched_strength_surrogate_split_baseline/per_patient_per_band.csv` | `N_nodes` |
| 3 | `ρ_split^coph` (matched-strength obs) | same | `obs_rho` |
| 4 | surrogate cohort median | same | `surr_p50` |
| 5 | z vs own surrogate | same | `obs_z` |
| 6 | `T_G^{*, p}` (normalized per-patient cluster mass) | `data/audit/grassmann_matched_strength_surrogate/per_patient_per_band_per_k.csv` | derived from `obs_p_one_sided_lower` over full `k` grid; see Methods Eq.~\eqref{eq:methods_TGstar_perpatient} |
| 7 | `n_sig, k` | same | count of `obs_p_one_sided_lower < α_k` cells |
| 8 | substrate `T_d^(d_S)` sign (sign of `−T_d^(d_S)`; `+` = trace, `−` = anti, `0` = `\|T_d^(d_S)\| < 0.01`) | `data/audit/raw_fc_phase_distance/Td_per_patient_per_band.csv` | sign of `S` (column `S` filtered `band=beta`) |

`N_nodes` is the count of contacts in the FC adjacency matrix after
patient-specific cleaning (Pat_10 task rows [53, 54, 55] dropped at
load; all other rows are vendor-canonical). The `|ImCoh|` adjacency
is dense and connected by construction, so `N_nodes` is also the
effective number of eigenmodes the LRG step operates on.

## Cohort-row sanity check (table footer ↔ paragraph claim)

| Footer cell | Table value | Source row | Paragraph claim | Match? |
|---|---|---|---|---|
| Col.\ 3 cohort median ρ_split | `+0.221` | `matched_strength_surrogate_split_baseline/cohort_summary.csv` `obs_median_rho` = 0.22057 | text cites `+0.222` (within-baseline ctm_triangle; 0.22223) — same cohort, different realization | ✓ (within 0.002) |
| Col.\ 3 asterisks | `**` | `paired_wilcoxon_p` = 0.00488 | `p = 0.005` | ✓ |
| Col.\ 4 surrogate cohort median | `+0.009` | `surr_median_rho_median` = 0.00931 | `+0.009` | ✓ |
| Col.\ 5 cohort count | `7/10` | `n_above_surrogate` = "7/10" | `7/10 above own surrogate` | ✓ |
| Col.\ 5 asterisks | `**` | `paired_wilcoxon_p` = 0.00488 | `p = 0.005` matched-strength | ✓ |
| Col.\ 6 cohort median `T_G^{*, p}` (normalized) | `0.653` | computed from `per_patient_per_band_per_k.csv` (cohort median of the per-patient normalized cluster mass) | implicit from paragraph "majority of patients carry the trace" | ✓ |
| Col.\ 6 asterisks | `**` | `grassmann_cluster_extent/cohort_summary.csv` β `cluster_p_cluster_mass` = 0.00498 | `p_mass = 0.005` cluster-mass permutation | ✓ |
| Col.\ 7 cohort count | `9/10` | patients with `T_G^{*, p} > floor` — 9 of 10 have ≥ 2 sig cells; Pat_15 has only 1 (floor) | implicit "majority on the Grassmann probe" | ✓ |
| Col.\ 7 asterisks | `**` | Anchored to band-level `p_mass = 0.005` (footnote c) | — | ✓ |
| Col.\ 8 substrate-comparison sign count | `7/10` | sign of `S` at β: Pat_02, 03, 05, 06, 08, 13, 14 in trace direction (`+`); Pat_10, Pat_15 anti (`−`); Pat_07 near zero (`0`, `\|T_d^(d_S)\| < 0.01`). **Not the LRG-trace gate** — trace classification is established only on the LRG measures (cols. 3–7). | brief §3.1 substrate `n_trace^(d_S) = 7/10` | ✓ |
| Bold-significance per-patient (LRG-trace gate) | 9/10 patients have at least one of cols. 3 or 6 bold | Pat_02, 03, 05, 06, 07, 08, 13: bold at both LRG probes; Pat_10, Pat_14: bold at T_G^{*, p} only (mixed); Pat_15: neither bold (LRG-anti-aligned) | derived; the cohort claim "majority of patients carry the trace at the LRG layer" is anchored here | ✓ |

The within-baseline ρ_split cohort median `+0.2222` (`ctm_triangle/cohort_summary.csv`)
and the matched-strength run cohort median `+0.2206` differ at the
fourth decimal but the paragraph cites the within-baseline value
(`+0.222`). The table footer reports the matched-strength value
(`+0.221`) so that the per-row z in col.\ 5 corresponds to the
ρ_split in col.\ 3 (same realization). Reader is alerted via the
caption clause "cohort median +0.221 ≈ +0.222 of the within-baseline
run cited in the text".

## Surprise patients — flag for the user

**Bolding convention (locked).** In the `.tex`, per-patient values in cols. 3 (`ρ_split^coph`) and 6 (`T_G^{*, p}`) are rendered **bold** when individually significant against the patient's own matched-strength surrogate (`z ≥ 1.96` on `ρ_split^coph`; `n_sig_k ≥ 20` on `T_G^{*, p}`, well above the null expectation of ≈ 5.5). A patient is LRG-trace iff at least one of the two cells is bold. Col. 8 (`d_S`) is the substrate analysis's sign-of-trace as comparison reference; it does **not** gate LRG-trace classification.

- **Pat_07** is **LRG-trace** — bold at both `ρ_split^coph` (+0.230,
  z = +4.47) and `T_G^{*, p}` (0.530, n_sig_k = 62). The substrate
  column shows `0` for Pat_07 (`T_d^(d_S) = +0.005`, below the
  `|T_d^(d_S)| < 0.01` near-zero threshold), reflecting the
  substrate analysis's view that there is no decisive sign at the
  edge-rank level. The substrate `0` is a comparison reference
  only; Pat_07's LRG-trace status is established by the bold values
  in cols. 3 + 6, not by col. 8.

- **Pat_15** is the **only patient with no LRG-significance** —
  neither col. 3 nor col. 6 is bold (`ρ_split^coph` = +0.083 with
  z = +1.27, not above own surrogate; `T_G^{*, p}` = 0.006, only 1
  marginal cell at the empirical-null floor). Pat_15 is therefore
  the lone LRG-anti-aligned patient. The right-hemisphere-only
  implant (0 epi contacts) is the most parsimonious mechanism —
  every other patient has at least partial left-hemisphere
  coverage. The **LRG-native `n = 9` robustness restriction**
  (drop Pat_15 only) preserves the trace at every probe:
  `ρ_split^coph` cohort median +0.230 with cohort-paired Wilcoxon
  p = 0.010; Grassmann `T_G(k=40)` strengthens to p = 0.002 (since
  Pat_15 was pulling the cohort median toward zero); drift-floor
  p = 0.027. The earlier legacy `n = 8` restriction also dropping
  Pat_07 is **retired** — Pat_07 is LRG-trace at both probes.

- **Pat_03** is fully in the cohort-typical direction at every LRG
  measure (`ρ_split^coph` +0.373 z = +9.16; `T_G^{*, p}` = 0.859,
  n_sig_k = 96 — bold at both probes). The Notes column flags the
  1024 Hz acquisition per current policy (handled at config layer;
  no dropout sensitivity test, no figure marker — see CLAUDE.md
  and `.agents/guides/04_rules/never-always-list.md` 2026-05-18
  updates).

- **Pat_10** and **Pat_14** are **mixed** at the LRG layer — bold
  at `T_G^{*, p}` (Pat_10: 0.612 with n_sig_k = 71; Pat_14: 0.695
  with n_sig_k = 78) but **not** bold at `ρ_split^coph` (Pat_10:
  −0.091, z = −1.70; Pat_14: −0.049, z = −1.00; both anti and
  below significance). The Notes column for these two rows reads
  "mixed". Both are LRG-trace under the at-least-one-LRG-bold rule
  (carried by the Grassmann probe). Per `bands/01_beta.md` §3.3,
  the two LRG probes read different facets of the same geometry,
  and a patient can carry one without the other.

- **Pat_13** has the highest epi-zone burden (30/119 contacts) but
  remains LRG-trace (`ρ_split^coph` +0.208 z = +3.51 → bold;
  `T_G^{*, p}` = 0.520 with n_sig_k = 61 → bold). The epi
  exclusion strengthens the cohort signal further (29/29
  manuscript cells persist + 7 emerge).

## Realization note (Option B chosen)

The ticket considered Option A (col.\ 3 within-baseline ρ_split,
+0.222 cohort median) and Option B (col.\ 3 matched-strength
obs_rho, +0.221 cohort median). The two runs use independent rsPre
A/B-split realizations; per-patient differences are 0.05–0.13 at β,
even though cohort medians are within 0.002. Option B was chosen
because cols.\ 4 and 5 (matched-strength surrogate median and z)
correspond row-for-row to col.\ 3 only under Option B. Switching to
Option A would require either (i) accepting that the z in col.\ 5 is
NOT against the ρ_split in col.\ 3, (ii) recomputing the surrogate
stats against the within-baseline ρ_split realization, which would
require a new R=200 surrogate run with the within-baseline seed
(audit_63 was seeded `20260511`; ctm_triangle uses a different
seed at the within-baseline A/B split).

## How to switch to Option A (if requested by user)

Replace col.\ 3 numbers from the matched-strength `obs_rho` with the
within-baseline `rho_split` from
`data/audit/ctm_triangle/Td_per_patient_per_band.csv`. The
per-patient mapping:

| Patient | within-baseline `rho_split` |
|---|---|
| Pat_02 | +0.380 |
| Pat_03 | +0.463 |
| Pat_05 | +0.483 |
| Pat_06 | +0.259 |
| Pat_07 | +0.185 |
| Pat_08 | +0.503 |
| Pat_10 | −0.123 |
| Pat_13 | +0.185 |
| Pat_14 | −0.073 |
| Pat_15 | +0.144 |

The cohort footer value changes from `+0.221` to `+0.222`. Cols.\ 4
and 5 stay as-is (they describe the matched-strength run); add a
caption footnote that the z in col.\ 5 is computed against
matched-strength `obs_rho`, not the within-baseline `rho_split` in
col.\ 3.

## Companion to figures

Patient ordering used: numerical by patient ID. If a future preprint
figure of the cohort (e.g. per-patient strip plot) uses a different
ordering (e.g. ascending ρ_split), update both this table and the
companion figure to share the same order.
