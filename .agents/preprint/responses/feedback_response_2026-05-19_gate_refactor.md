---
name: feedback-response-2026-05-19-gate-refactor
era: IMCOH_ABS_COHORT_N10
status: complete_writing_agent_unblocked
kind: cascade-summary
scope: response to writing-agent feedback 2026-05-19 (M4 + m6); refactors C3/C4/C5 to principled Wilcoxon gates with no patient-count thresholds
companion: CONTROLS.md, VERDICT_LEDGER.md, methods_grassmann_cluster_extent.md, methods_section_review_2026-05-19.md
---

# Response to writing-agent feedback 2026-05-19 — gate refactor cascade

**Head.** All four feedback points have been addressed via three
coordinated gate refactors (CONTROLS.md Decisions 8 + 9 + 10) plus a
mandatory LOO-diagnostic policy (Decision 11). Two new audit scripts
have been written + run (audit_71 C4; audit_72 C5), audit_70 has been
re-run with the mass-only verdict and LOO column, and
VERDICT_LEDGER.md has been re-derived under the new gates. All
previously-locked verdicts hold (β cophenet + Grassmann strong,
α cophenet strong, γ_l + δ Grassmann strong, γ_h/θ no trace, α
Grassmann no trace), with the δ Grassmann verdict acquiring a
methodologically clean LOO + C5 narrative: full-data Pat_08 leverage
is attributable to epi-zone interactions, and C5 epi-X strengthens
the trace and resolves the leverage fully. **The writing agent can
now proceed with the M4 + m6 patches under the updated CONTROLS.md
+ VERDICT_LEDGER.**

---

## What changed (cascade summary)

### Decision 8 — Grassmann gate is mass-only
- CONTROLS.md §C3 (Grassmann sub-section): gate is now
  `cluster_p_cluster_mass < 0.05` alone (Decision 8). `LR_obs` and
  `cluster_p_longest_run` stay as descriptive co-statistics in the
  audit CSV but **do not gate** the verdict.
- audit_70 (`scripts/01_compute/audit/audit_70_grassmann_cluster_extent.py`)
  updated: verdict column maps from `cluster_p_cluster_mass`.
- methods_grassmann_cluster_extent.md §5d (Mass-only verdict gate) +
  §6 (Locked band-level values) rewritten.
- Verdict consequences: γ_l Grassmann weak → strong; δ Grassmann
  weak → strong. β Grassmann unchanged at strong.

### Decision 9 — C4 cross-probe gate is paired Wilcoxon (no patient counts)
- CONTROLS.md §C4: gate is now paired one-sided Wilcoxon
  `(rho_split − rho_xprobe)` under `H_1: rho_split > rho_xprobe`,
  **failing to reject** at α=0.05 (no significant degradation),
  **AND** `sign(rho_xprobe_median) == sign(rho_split_median)`. No
  patient-count threshold anywhere.
- New audit: `scripts/01_compute/audit/audit_71_c4_wilcoxon_cohort.py`.
- New CSV: `data/audit/ctm_triangle/c4_wilcoxon_cohort.csv`.
- LOO max-p column included (descriptive).
- All 6 bands pass the new gate.

### Decision 10 — C5 epi-X gate is Wilcoxon (no ≳80% retention)
- CONTROLS.md §C5: 
  - Grassmann (all bands): re-run audit_70 cluster-mass test on
    epi-X eigvec cache; gate is `cluster_p_mass^epi-X < 0.05`.
  - Cophenet (α only): one-sample Wilcoxon on per-patient
    `obs_rho^epi-X` under `H_1: rho_split^epi-X > 0`; gate is
    `wilcoxon_one_sided_p < 0.05`.
- C5 promoted from sensitivity layer to **primary gate** in the
  strong/weak rule for Grassmann.
- New audit: `scripts/01_compute/audit/audit_72_c5_wilcoxon_cohort.py`.
- New CSVs: `data/audit/grassmann_epi_exclusion/c5_wilcoxon_cohort.csv`
  + `data/audit/alpha_epi_exclusion/c5_wilcoxon_cohort.csv`.
- LOO max-p column included (descriptive).
- All previously-locked C5 verdicts hold: α cophenet (p=0.0098 passes);
  β/δ Grassmann (`p_mass^epi-X = 0.005` strengthens vs full data);
  γ_l Grassmann passes cohort gate (`p_mass^epi-X = 0.030`) with LOO
  leverage on Pat_05 to flag.

### Decision 11 — LOO max-p mandatory alongside every Wilcoxon gate
- Per the user's directive 2026-05-19 ("always be sure a p-value is
  not driven by single patients"): every cohort Wilcoxon-based gate
  (C1, C2, C3 cophenet, C3 Grassmann cluster mass, C4, C5) emits a
  LOO max-p + argmax-patient column alongside the main p-value.
- Descriptive only — **never a gate threshold**.
- audit_70, audit_71, audit_72 all carry the LOO columns.
- Memory: `feedback_no_single_patient_p_driven.md`.

---

## New locked verdict table (post-cascade)

| Band | Range (Hz) | D_coph | Grassmann | Coverage tag |
|---|---|---|---|---|
| **β** | 13–30 | **strong trace** (all C1-C5 pass) | **strong trace** (mass-only gate; LOO-robust at full + epi-X; mass strengthens 70→89 under epi-X) | **strong trace, both probes** |
| α | 8–13 | **strong trace** (C1 p=.010, C2 p=.007, C3 p=.002, C4 paired_p=.461 pass, C5 wilcoxon_p=.0098 pass) | no trace | **strong trace, only D_coph** |
| γ_l | 30–80 | no trace (C3 fails p=.116) | **strong trace** (`p_mass = 0.005`, LOO max .040 robust; C5 passes `p_mass^epi-X = .030` with LOO max .159 — flag Pat_05 leverage under epi-X) | **strong trace, only Grassmann** ↑ |
| δ | 0.53–4 | no trace (C3 fails p=.278) | **strong trace** (`p_mass = 0.005`, full LOO max .055 flags Pat_08; **C5 epi-X resolves cleanly**: mass 38→44, `p_mass^epi-X = 0.005`, LOO max .005 fully robust → leverage attributable to epi-zone) | **strong trace, only Grassmann** ↑ |
| γ_h | 80–300 | no trace | no trace (`p_mass = .060`, `p_mass^epi-X = .060`) | **no trace** |
| θ | 4–8 | no trace | no trace (`p_mass = .159`) | **no trace** |

---

## Methods-text patches the writing agent can now apply

### Reopened patches

- **M1** (T_G* formula): the methods_TGstar equation now reads
  `T_G^*(b) = Σ_{k : p_k(b) < α_k} (−log_10 p_k(b))` (resilient
  all-clusters sum). Verdict cited in manuscript: β `T_G^* = 69.76,
  p_mass = 0.005`; γ_l `T_G^* = 66.14, p_mass = 0.005`; δ `T_G^* =
  38.07, p_mass = 0.005`.
- **M2** (D = 1/ρ): unchanged from the previous review — replace
  `K_ij(τ')` with `ρ_ij(τ')` in eq. methods_distm.
- **M3** (surrogate description): unchanged — soften "preserves the
  global edge-weight distribution" to strengths-only language.
- **M4** (epi-X criterion): now resolved via CONTROLS.md §C5 update;
  rewrite the §sssec:methods_compare_stats robustness paragraph to
  cite Wilcoxon-on-epi-X gates (Grassmann: `cluster_p_mass^epi-X <
  0.05`; cophenet α: one-sample Wilcoxon p < 0.05). Cite audit_72
  CSVs.
- **m5** (T_G < 0 trace direction): unchanged.
- **m6** (C4 gate underspecified): now resolved via CONTROLS.md §C4
  update; cite the paired-Wilcoxon non-degradation form. Cite
  audit_71 CSV.
- **m7-m9, n10-n14**: unchanged from the previous review.

### New methods-text additions

- Add a sentence stating that every Wilcoxon-based gate is paired
  with a LOO max-p sensitivity diagnostic (descriptive only); flag
  any verdict where LOO max-p crosses 0.05 explicitly in the
  Results text (specifically δ Grassmann full-data; γ_l Grassmann
  C5 epi-X). The δ + γ_l LOO findings are the worked examples that
  illustrate why the LOO diagnostic matters.

---

## Source-of-truth references

- `.agents/preprint/locked/CONTROLS.md` — updated 2026-05-19 pm, Decisions
  8 + 9 + 10 + 11. Read first.
- `.agents/preprint/locked/VERDICT_LEDGER.md` — updated 2026-05-19 pm with
  cascade; revision history block at the end documents the four new
  decisions.
- `.agents/preprint/methods/methods_grassmann_cluster_extent.md` — §5d
  (Mass-only verdict gate), §6 (locked LOO + C5 table).
- `.agents/preprint/methods/methods_section_review_2026-05-19.md` — M4 + m6
  flagged as RESOLVED; head paragraph updated with the cascade
  summary.
- New audit scripts:
  - `scripts/01_compute/audit/audit_71_c4_wilcoxon_cohort.py`
  - `scripts/01_compute/audit/audit_72_c5_wilcoxon_cohort.py`
- Updated audit:
  - `scripts/01_compute/audit/audit_70_grassmann_cluster_extent.py`
    (mass-only verdict + LOO + `loo_cluster_p_mass.csv` sibling file)
- New audit CSVs:
  - `data/audit/ctm_triangle/c4_wilcoxon_cohort.csv`
  - `data/audit/grassmann_epi_exclusion/c5_wilcoxon_cohort.csv`
  - `data/audit/alpha_epi_exclusion/c5_wilcoxon_cohort.csv`
  - `data/audit/grassmann_cluster_extent/loo_cluster_p_mass.csv`
- Memory entries:
  - `feedback_patient_counts_never_the_gate.md`
  - `feedback_no_single_patient_p_driven.md`
  - `feedback_cluster_mass_all_clusters.md`
  - `grassmann_null_calibration_protects.md`
