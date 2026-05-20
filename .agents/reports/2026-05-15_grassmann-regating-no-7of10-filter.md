---
name: 2026-05-15_grassmann-regating-no-7of10-filter
era: IMCOH_ABS_COHORT_N10
status: current
kind: report
audit: 69
upstream: [audit_66, audit_67]
---

# Grassmann regating — drop hardcoded patient-count + magnitude filters

The §5.4 audit_66 / audit_67 cohort verdict tagged a (band, k) cell **separated** only when **three** criteria hit at once: cohort-paired one-sided Wilcoxon `p < 0.05` (the statistical gate), `|cohort median surr T_G| < 0.05 × max(1, |cohort median obs T_G|)` (a hardcoded magnitude-ratio filter), and `n_below ≥ 8/10` patients individually below their own surrogate at one-sided `p < 0.05` (a hardcoded patient-count direction-agreement filter). The two extra filters were holdovers from the §5.4 small-multiples visualization convention and were never disclosed as gating criteria in the §7 errata or in the §5.4 manuscript text. This regating run drops both and gates every (band, k) cell on the cohort Wilcoxon `p < 0.05` alone.

## Per-band side-by-side: old gate vs Wilcoxon-only gate

| audit | band | n_k | n_sig old | n_sig new | Δ | old run (k_start–k_end, len) | new run (k_start–k_end, len) | best-k ratio |
|---|---|---|---|---|---|---|---|---|
| audit_66_full | delta | 111 | 0 | 23 | +23 **⚠** | — | k=57–63, L=7 | k=87 (+4.89) |
| audit_66_full | theta | 111 | 0 | 8 | +8 **⚠** | — | k=79–82, L=4 | k=82 (+3.72) |
| audit_66_full | alpha | 111 | 0 | 4 | +4 **⚠** | — | k=11–14, L=4 | k=11 (+1.77) |
| audit_66_full | low_gamma | 111 | 11 | 41 | +30 **⚠** | k=13–16, L=4 | k=12–23, L=12 | k=85 (-74.72) |
| audit_66_full | beta | 111 | 0 | 40 | +40 **⚠** | — | k=27–55, L=29 | k=21 (-207.21) |
| audit_66_full | high_gamma | 111 | 9 | 19 | +10 **⚠** | k=22–27, L=6 | k=19–27, L=9 | k=25 (+341.32) |
| audit_67_epiX | delta | 87 | 0 | 25 | +25 **⚠** | — | k=33–39, L=7 | k=72 (+3.51) |
| audit_67_epiX | theta | 87 | 0 | 8 | +8 **⚠** | — | k=63–64, L=2 | k=68 (+2.95) |
| audit_67_epiX | alpha | 87 | 0 | 3 | +3 | — | k=10–12, L=3 | k=10 (+2.26) |
| audit_67_epiX | low_gamma | 87 | 0 | 20 | +20 **⚠** | — | k=19–28, L=10 | k=7 (+9.72) |
| audit_67_epiX | beta | 87 | 4 | 52 | +48 **⚠** | k=21–23, L=3 | k=21–56, L=36 | k=17 (+41.56) |
| audit_67_epiX | high_gamma | 87 | 4 | 16 | +12 **⚠** | k=25–26, L=2 | k=21–26, L=6 | k=18 (-1267.53) |

**Δ flag (⚠)** marks cells where the new gate changes the count by more than ±3 cells. Sign of Δ is `new − old` — positive means the Wilcoxon-only gate accepts cells the old composite gate rejected.

## Manuscript-window load-bearing claims (§7 errata)

The §7 errata cites three Grassmann manuscript windows under the old composite gate. Recomputed under the Wilcoxon-only gate inside the same window:

| audit | band | window k=k0..k1 | window cells | n_sig old in window | n_sig new in window | longest contig run (new) within window |
|---|---|---|---|---|---|---|
| audit_66_full | beta | 27..55 | 29 | 0 | 29 | L=29 (k=27–55) |
| audit_67_epiX | beta | 27..55 | 29 | 0 | 29 | L=29 (k=27–55) |
| audit_66_full | low_gamma | 12..23 | 12 | 4 | 12 | L=12 (k=12–23) |
| audit_67_epiX | low_gamma | 12..23 | 12 | 0 | 10 | L=5 (k=13–17) |
| audit_66_full | high_gamma | 19..27 | 9 | 6 | 9 | L=9 (k=19–27) |
| audit_67_epiX | high_gamma | 19..27 | 9 | 2 | 6 | L=6 (k=21–26) |

## Verdict per window (new vs old)

- **beta k=27..55** — audit_66 window count 0 → 29 (expanded, Δ=+29); audit_67 window count 0 → 29 (expanded, Δ=+29). Manuscript headline of **29 contiguous + 29/29 persist** is reproduced exactly by the Wilcoxon-only gate.
- **low_gamma k=12..23** — audit_66 window count 4 → 12 (expanded, Δ=+8); audit_67 window count 0 → 10 (expanded, Δ=+10). Manuscript headline of **12 contiguous + 10/12 persist** is reproduced exactly.
- **high_gamma k=19..27** — audit_66 window count 6 → 9 (expanded, Δ=+3); audit_67 window count 2 → 6 (expanded, Δ=+4). Manuscript headline of **9 contiguous** is reproduced exactly. **Persist count cited in §7 errata is "7/9"; current `sensitivity.csv` reports 6 persist + 3 weaken** (persist ks = 21, 22, 23, 24, 25, 26; weaken ks = 19, 20, 27). This is a single-cell discrepancy in the errata text; the correct manuscript number is **6/9 persist + 3/9 weaken**. Likely a tally error in the original §7 draft — not a filter artifact.

## Headline implication for §7 errata

- The three Grassmann manuscript windows (β k=27..55, γ_l k=12..23, γ_h k=19..27) and their cited contiguous-cells counts (29, 12, 9) **were already on the Wilcoxon-only gate** in the underlying audit_67 `sensitivity.csv`. The §7 errata is not contaminated by the patient-count or magnitude-ratio filter — those filters only entered the `verdict == "separated"` label and the README `longest_run_sep` column inside `data/audit/grassmann_*/README.md`, which the manuscript never cites.
- Two consequences:
  1. **The 29 / 12 / 9 contiguous counts and 29/29 / 10/12 persist counts stand.** No manuscript re-derivation needed for those numbers.
  2. **The single-cell γ_h "7/9 → 6/9" discrepancy must be corrected in the §7 errata text and downstream summaries.** It is not a regating artifact; it is a transcription/tally error in the original errata draft.
- The in-code `verdict == "separated"` label and `longest_run_sep` table in `audit_66`/`audit_67` README outputs remain misleading and should be removed in the next pipeline pass, since they advertise 0 separated cells for β under the composite gate even though β is the strongest matched-strength-controlled signal in §5.4 under the gate the manuscript actually uses.

## What the new gate means

Under the Wilcoxon-only gate, **sig_cell(band, k)** is the cohort-paired one-sided Wilcoxon signed-rank test of per-patient `T_G_obs` vs per-patient `T_G_surr_median` rejecting `H_0: T_G_obs ≥ T_G_surr` at α=0.05. `n_below_own_surrogate` is reported as a descriptive column but no longer enters the gate. The magnitude-ratio criterion is also dropped — Wilcoxon ranks already encode the magnitude evidence the ratio filter was duplicating.

## Provenance

- Inputs:
  - `data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv`
  - `data/audit/grassmann_epi_exclusion/cohort_summary.csv`
- Outputs:
  - `data/audit/grassmann_regate_no_filter/per_cell_audit66.csv`
  - `data/audit/grassmann_regate_no_filter/per_cell_audit67.csv`
  - `data/audit/grassmann_regate_no_filter/contig_summary.csv`
  - `data/audit/grassmann_regate_no_filter/window_check.csv`
- Build script:
  - `scripts/01_compute/audit/audit_69_grassmann_regate_no_filter.py`
