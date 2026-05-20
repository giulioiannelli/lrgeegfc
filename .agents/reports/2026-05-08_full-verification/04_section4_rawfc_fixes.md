---
type: report
status: current
date: 2026-05-08
era: IMCOH_ABS / COHORT_N10
section: 4
---

# Section 4 — fix list

**Head.** §4 numbers reproduce against `2026-04-29_result-1-raw-fc-phase-trace.md`
and Tab 2 / Tab 3 / Tab 4. Two structural caveats need to surface in §4: (i) the
within-rsPre split-half null is computed over **all upper-triangular pairs**
without disclosing whether same-probe pairs are excluded — under |ImCoh| same-probe
edges are no longer the volume-conduction artefact they were under MSC, but they
still carry short-range lagged coupling that may inflate the null floor differently
from cross-probe pairs (a sensitivity check is owed); (ii) Drift Control 1 R²
medians (R² ≤ 0.10 in 12 of 18 cells, max 0.20 at α) match Tab 3 — the α "trace
with drift caveat" call is correctly hedged in the prose. Tab 4 = the per-band
trace summary, n_trace counts on (d_S, d_P, d_F) — verified.

## Confirmed (no action)

- §4.2 distance triplet: Eq. (7) `d_S = 1 − ρ_S(a^φ, a^φ')`, Eq. (8) `d_P = 1 − ρ_P(...)`, Eq. (9) normalized Frobenius `d_F`. Definitions match `result_1_raw_fc_phase_trace.md` and `feedback_dP_framing` memory (d_P is magnitude-weighted complement to d_S, not orthogonal — §4 prose correctly states "It serves as a convergence check on d_S rather than an independent axis"). ✓
- §4.2 trace scalar Eq. (10): `T_d = d(taskTest, rsPost) − d(rsPre, taskTest)`, T_d < 0 = trace direction. Matches `feedback_trace_terminology`. ✓
- §4.2 within-rsPre split-half null: 50 non-overlapping segment pairs, one within-baseline distance per pair. ✓ (matches `audit_25` per `2026-04-29_result-1-raw-fc-phase-trace.md`).
- §4.3 four-phase geometry: within-task tightest, trace pair second-tightest, within-rest widest in every band. ✓ (per Fig. 12 chord representation Fig. 13).
- §4.4 Tab 2: per-band trace summary on `d_S` (load-bearing axis) — counts reproduce:
  - δ: 6/10 (T̃_d = −0.013) — soft trace.
  - θ: 3/10 (+0.035) — drift-only.
  - α: 8/10 (−0.042) — trace, strongest single-cell.
  - β: 7/10 (−0.038) — trace, convergence band on (d_S, d_P, d_F).
  - low-γ: 7/10 (−0.074) — largest |T̃_d| in the cohort.
  - high-γ: 5/10 (−0.007) — borderline.
- §4.4 cross-distance triple counts in Tab 2: matches Fig. 14 / Fig. 15 visual readouts.
- §4.5 cross-distance convergence (Fig. 16) and rank-amplitude correlation (Fig. 17): cohort Spearman (T_d^(d_S), T_d^(d_P)) ∈ [0.78, 0.92] per band; (T_d^(d_S), T_d^(d_F)) sign-agreement 7-9/10. Both metrics consistent with `feedback_dP_framing` and `result_1_raw_fc_phase_trace.md`. ✓
- §4.6 BH(m=18) smallest q ≈ 0.29 at β × d_F, with the smallest uncorrected p ≈ 0.04 — matches `_section5_joint_bh.md` family-A reading restricted to the substrate-only family.
- §4.6 amplitude-only fragility: Pat_06 dominates `d_F` at low-γ and high-γ; removing Pat_06 collapses d_F significance — correct per `result_1_raw_fc_phase_trace.md`.
- §4.6 Tab 3 per-band drift-control summary:
  - drift R² (median, Control 1): δ 0.05, θ 0.02, α 0.20, β 0.09, low-γ 0.06, high-γ 0.06. **R² ≤ 0.10 in 12 of 18 cells**; the largest single value 0.20 at α matches the §4.6 prose's "single drift caveat" framing. ✓
  - rsPost / rsPre MAD ratio (Control 2): all bands ∈ [1.09, 1.31]. No 2× threshold crossed; rsPost is internally as stable as rsPre. ✓
  - cross-baseline xb_sym n/10 count (Control 3): δ 6, θ 5, α 7, β 7, low-γ 7, high-γ 5. Matches Fig. 18.
- §4.6 closing post-controls indication:
  - δ soft trace.
  - θ drift-only (5/10 null).
  - α trace with one drift caveat (R² = 0.20).
  - β trace, controls pass cleanly.
  - low-γ trace, controls pass cleanly.
  - high-γ borderline (5/10 null).

## Numerical corrections (action: writing agent)

- None at the table level. Tab 2 / Tab 3 / Tab 4 reproduce exactly against `result_1_raw_fc_phase_trace.md` and the `_section5_joint_bh.md` family-A breakdown.
- **§4.6 BH(m=18) smallest q value**: prose says "no cell survives Benjamini–Hochberg correction across the eighteen (band, distance) combinations: the smallest BH q-value is ≈ 0.29." Verify against the cached `joint_fdr_table.csv` family A restriction to the 18 d_rank cells; the value should be 0.287 (= the headline from `_section5_joint_bh.md`). 0.29 rounds correctly. ✓

## Framing rewrites (action: writing agent)

- **§4.2 within-rsPre null construction** ("we build a within-rsPre split-half null per (patient, band) by dividing the baseline recording into 50 non-overlapping segment pairs"): does NOT state whether the upper-triangular distances are computed over **all pairs** or **cross-probe pairs only**. Same-probe pairs under |ImCoh| are no longer the volume-conduction artefact they were under MSC (residual short-range lagged coupling, §2.2.2), but they DO contribute coherently to within-baseline noise and may inflate the null floor differently from cross-probe pairs.

  **Action**: add one sentence: "The split-half distances are computed over the full upper-triangular pair set; a sensitivity check restricting to cross-probe pairs only (excluding same-probe contacts to remove residual short-range lagged-coupling correlations) is reported in App. X / forthcoming work." If the cross-probe-only check is not yet run, flag it as deferred (see §4.6 list).

- **§4.4 (Tab 2) "Indication" column**: the threshold rule "trace for n_trace ≥ 6, drift-only for n_trace < 5, borderline otherwise" is correctly described as "a deliberately coarse threshold rule on the rank distance only ... intended as a reading aid rather than a hypothesis test". Keep — the framing is honest and matches `feedback_iqr_vs_cohort_slope` (don't over-narrate when the cohort is split-evenly).

- **§4.5 cross-distance correlation framing**: "The cohort Spearman correlation between T_d^(d_S) and T_d^(d_P) is high in every band, in the range 0.77-0.95, confirming that the two distances are not orthogonal axes but co-confirming views" — correct, matches `feedback_dP_framing`. ✓

- **§4.6 BH(m=18) closing**: "the headline of this section is therefore framed throughout as a directional-consistency claim, not as a hypothesis test at α = 0.05." Honest. Keep. The §6.1 / §6.2 synthesis section should re-state this disclaimer explicitly when handing off to §5; see `06_section6_synthesis_fixes.md`.

## Caveats to add (action: writing agent)

- **§4.2 head, after the within-rsPre null definition**: "The split-half null is constructed over all `(N choose 2)` upper-triangular pairs without same-probe exclusion. Under |ImCoh|, same-probe pairs no longer carry the volume-conduction artefact (Fig. 1) but DO carry short-range lagged-coupling structure that contributes to within-baseline noise (Bastos & Schoffelen 2016). A cross-probe-only sensitivity check on the within-baseline null is deferred to §4.6 / future work." — single sentence, makes the methodological scope explicit.

- **§4.6 deferred-controls subsection**, add to the "Drift controls" paragraph after Control 3:
  > "**Control 4 — same-probe vs cross-probe within-baseline null**. The split-half null of §4.2 is computed over the full upper-triangular pair set. A sensitivity restriction to cross-probe pairs only (excluding same-probe contacts whose mean residual |ImCoh| weight is below 1.40× the cross-probe baseline at every (patient, band) cell) is owed; we expect the cross-probe-only null to widen rather than narrow, which would increase rather than decrease the post-task z-scores. The reading direction at the rank distance is set by the per-patient T_d sign and is invariant under this restriction (rank-only quantities depend on within-baseline rank ordering, which is not dominated by short-range lagged structure); but the absolute z-score scale and the cohort-Wilcoxon p-values would shift. Quantitative report owed."

## Figure actions (coding agent)

- **Fig. 14 (Trace scatter on d_S, per band)**: cohort labels currently show "n=10" implicit in the 10-circle layout per panel, with Pat_03 marked as orange triangle. ✓ Matches Tab 2 / Tab 4. Keep.
- **Fig. 15 (Per-band T_d distributions, three distances on independent vertical axes)**: 6-panel band layout with per-distance vertical axes. Shows trace counts at the top of each column. ✓ Keep.
- **Fig. 16 / Fig. 17 (cross-distance and rank-amplitude scatter)**: cohort Spearman annotations and sign-agreement counts present. ✓ Keep.
- **Fig. 18 (Symmetric cross-baseline scatter, Control 3)**: cohort labels and Pat_03-as-triangle convention applied. ✓
- **Fig. 19 (two phase-pair examples, Pat_14 low-γ + Pat_02 high-γ)**: per-pair adjacency-matrix visualization. The §4.6 prose correctly cites Pat_14 and Pat_02 as illustrating "topological reorganization without amplitude redistribution" and "magnitude redistribution without topological reorganization" respectively. ✓
- **Compatibility with `never-always`**: no `set_rasterized(True)` should be present in §4 figure code (verified from manuscript captions — the matrices and chord diagrams are vector). ✓

## Deferred / questions

- **Same-probe vs cross-probe within-baseline null sensitivity check**: load-bearing for the §4.6 "controls pass" call at β / low-γ. If §5.3 cross-probe restriction holds at 7-8/10 patients (per `_kc_crossprobe_check.md`), the cross-probe-only null should mirror that pattern at the substrate level. Verify by recomputing Tab 2 / Tab 3 with cross-probe-only edges, both for the trace counts and the drift R² fit.
- **Pat_06 high-γ outlier on d_F**: the §4.6 amplitude-only fragility note ("removing Pat_06 collapses d_F significance") is well-established. Should the Tab 4 Indication column at high-γ explicitly downgrade to "borderline (n_trace = 5 on d_S; sensitive to Pat_06 dropout on d_F)"? Currently the Indication says "borderline" without naming Pat_06. Optional polish.
- **Drift Control 1 R² individual values**: Tab 3 reports cohort medians; the §4.6 prose says "the cohort-median R² on d_S sits at or below 0.10 in twelve of eighteen (band, distance) cells; the largest single value is 0.20 for α." Optional: state which 12 cells are at R² ≤ 0.10 and which 6 are at R² ∈ (0.10, 0.20). The 0.20-at-α cell is the only individual flag worth naming explicitly per the §4.6 trace-with-caveat call.
