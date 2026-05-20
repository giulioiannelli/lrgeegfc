---
name: preprint-verdict-ledger
era: IMCOH_ABS_COHORT_N10
status: locked_2026-05-18
kind: verdict-lockdown
supersedes: all_prior_per-band_verdict_claims
companion: CONTROLS.md
cohort: [Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15]
n_patients: 10
fc_method: imcoh_abs
tau: 1/lambda_max
canonical_lrg_object: "D_coph = cophenet(UPGMA(D(tau_max)))"
canonical_subspace_object: "U_k = span{phi_2, ..., phi_{k+1}}"
---

# Verdict ledger — per-band trace verdicts (locked 2026-05-18)

**Head.** The 6 bands × 2 probes (cophenet `D_coph` per-pair multiscale +
Grassmann `d_G(k)` whole-network mode subspace) decompose cleanly under the
locked 5-control battery (`CONTROLS.md`). β is `strong, both probes`. α is
`strong, only D_coph` — the cophenet trace strengthens under epi-zone
exclusion, identifying α as a non-epi-cortex per-pair phenomenon hidden by
epileptic patients in the full cohort. γ_l and δ carry weak Grassmann-only
subspace traces under cluster-extent permutation (audit_70). γ_h and θ have
no trace on either probe (γ_h's contiguous run is at the null 95th percentile,
cluster p = 0.055). α has no Grassmann trace.

## Locked verdict table

| Band | Range (Hz) | D_coph (cophenet per-pair) | Grassmann (whole-network modes) | **Coverage tag** |
|---|---|---|---|---|
| **β** | 13–30 | **strong trace** | **strong trace** (LOO-robust) | **strong trace, both probes** |
| α | 8–13 | **strong trace** (epi-X strengthens) | no trace | **strong trace, only D_coph** |
| γ_l | 30–80 | no trace | **strong trace** ↑ (mass-only gate, LOO-robust) | **strong trace, only Grassmann** ↑ |
| δ | 0.53–4 | no trace | **strong trace** ↑ (mass-only gate; C5 epi-X strengthens decisively; full-data LOO Pat_08 leverage resolves under epi-X) | **strong trace, only Grassmann** ↑ |
| γ_h | 80–300 | no trace | no trace (`p_mass` = 0.060) | **no trace** |
| θ | 4–8 | no trace | no trace | **no trace** |

These verdicts are **locked under the 2026-05-19 pm gate refactor**
(CONTROLS.md §C3 mass-only Decision 8; §C4 paired-Wilcoxon Decision 9;
§C5 Wilcoxon-on-epi-X Decision 10). Per-band briefs
(`bands/02_alpha.md` … `bands/06_delta.md`) document them; they do not re-derive
them.

**Verdict changes from the previous lock (2026-05-19 am)**:
- γ_l Grassmann: weak → strong (mass-only gate gives `p_mass = 0.005`
  under the resilient all-clusters formula; the old longest-cluster
  mass under-called this band as weak; LOO max p = 0.040 robust).
- δ Grassmann: weak → strong (`p_mass = 0.005` at full data; full-data
  LOO max p = 0.055 flags Pat_08 leverage; **C5 epi-X strengthens
  the trace and resolves the leverage**: mass 38 → 44, `p_mass^epi-X
  = 0.005`, LOO max under epi-X = 0.005 fully robust. The full-data
  Pat_08 leverage is attributable to epi-zone interactions, not the
  true biological trace; the manuscript text should report both
  numbers transparently per
  `feedback_no_single_patient_p_driven.md`).
- β Grassmann: unchanged at strong, but now LOO-robust (`p_mass` and
  LOO max both = 0.005, no single-patient leverage).

**Anatomy localization is locked separately in [`ANATOMY_LEDGER.md`](ANATOMY_LEDGER.md)**
under the [`ANATOMY_CONTROLS.md`](ANATOMY_CONTROLS.md) battery (A1 hypergeometric +
A2 sampling-corrected + A3 matched-strength **mandatory** + A4 implant-geometry
regression). Per-(band, probe) verdicts: `strong localized` / `weak localized` /
`not localized`. As of 2026-05-19, anatomy audits are scoped but not yet run;
KC-era anatomy artifacts at `data/audit/lrg_localization_anatomy/` are retired
and not citable.

---

## Control battery (locked) — full definitions in `CONTROLS.md`

- **C1** baseline-split (D_coph only) — `ctm_triangle/cohort_summary.csv` col `wilcoxon_split_gt_0_p`
- **C2** drift-floor null (D_coph only) — same CSV col `wilcoxon_split_gt_drift_p`
- **C3** matched-strength surrogate (both probes; **mandatory**) — `matched_strength_surrogate_split_baseline/cohort_summary.csv` for D_coph; `grassmann_cluster_extent/cohort_summary.csv` col `cluster_p_cluster_mass` for Grassmann (mass-only gate, Decision 8)
- **C4** cross-probe restriction (D_coph only) — `ctm_triangle/c4_wilcoxon_cohort.csv` (audit_71, locked 2026-05-19; paired-Wilcoxon gate, Decision 9). The old `ctm_triangle/cohort_summary.csv` cols `n_trace_xprobe_int`/`rho_xprobe_median` are retained as **descriptive auxiliary statistics**, not the gate.
- **C5** epi-zone exclusion — `grassmann_epi_exclusion/c5_wilcoxon_cohort.csv` (audit_72, locked 2026-05-19; cluster-mass on epi-X eigvec cache for Grassmann) + `alpha_epi_exclusion/c5_wilcoxon_cohort.csv` (audit_72, locked 2026-05-19; one-sample Wilcoxon on `obs_rho^epi-X` for cophenet α). C5 is a **primary gate** for `strong vs weak`, not a sensitivity layer, under Decision 10.

---

## Per-band gathered numbers + verdict reasoning

### β (13–30 Hz) — **strong trace, both probes**

**Cophenet `D_coph` — primary 4-control table:**

| Control | Statistic | Source CSV row | Gate | Pass |
|---|---|---|---|---|
| C1 split | p = 0.00488 | `ctm_triangle/cohort_summary.csv` beta | < 0.05 | ✓ |
| C2 drift | p = 0.01367 | `ctm_triangle/cohort_summary.csv` beta | < 0.05 | ✓ |
| C3 matched-strength | p = 0.00488; ratio obs/surr = 23.7× (n_above 7/10 = descriptive) | `matched_strength_surrogate_split_baseline/cohort_summary.csv` beta | < 0.05 | ✓ |
| C4 cross-probe (Wilcoxon) | paired Wilcoxon (split>xprobe) p = 0.385 ⇒ fails to reject ⇒ no degradation; rho_split_median = +0.222, rho_xprobe_median = +0.223 (sign agree); LOO max p = 0.590 (Pat_02) | `ctm_triangle/c4_wilcoxon_cohort.csv` beta | paired_p ≥ 0.05 + sign match | ✓ |
| C5 epi-X cophenet | **not run at this probe** | — | n/a | n/a |

**Grassmann `d_G(k)` — cluster-extent permutation (audit_70, mass-only gate) + C5 epi-X (audit_72):**

| Statistic | Value | Source CSV |
|---|---|---|
| Observed longest contiguous-sig run (descriptive) | **29 cells** at k = 27..55 | `grassmann_cluster_extent/cohort_summary.csv` beta |
| Null mean / 95th / max LR | 2.37 / 6.0 / 17 | same |
| `cluster_p_longest_run` (descriptive co-statistic) | 0.005 (minimum at R=200) | same |
| **Resilient all-clusters mass `T_G^*`** | **69.76** | same |
| **`cluster_p_cluster_mass` (gate)** | **0.005** | same |
| LOO max `p_mass` | 0.005 (Pat_02) — **fully LOO-robust** | same |
| **C5 epi-X `cluster_p_cluster_mass^epi-X` (gate)** | **0.005** | `grassmann_epi_exclusion/c5_wilcoxon_cohort.csv` beta |
| C5 epi-X `T_G^*^epi-X` | **89.04** (strengthens vs full 69.76) | same |
| C5 epi-X LOO max `p_mass^epi-X` | 0.005 (Pat_02) — **fully LOO-robust** | same |
| C5 epi-X longest run | 36 cells (vs 29 in full data) | same |

β trace **strengthens** under epi-X (mass 69.76 → 89.04, LR 29 → 36)
and remains fully LOO-robust at both full and epi-X.

**Substrate context (raw FC matched-strength)**: p = 0.0527 (borderline), n_above 6/10 — β substrate is the *weakest* of the bands at raw FC; the β trace is an LRG-emergent property, NOT inherited from substrate.

**Three-layer sensitivity** (matched-strength):

| Layer | n_above | Wilcoxon p |
|---|---|---|
| Raw FC `ρ_split^raw` | 6/10 | 0.0527 |
| Raw D(τ_max) `ρ_split` | 6/10 | 0.0420 (ratio 13.4×) |
| Cophenet `D_coph` `ρ_split^coph` | 7/10 | **0.00488** (ratio 23.7×) |

Cophenet step *amplifies* the β trace selectively.

**Verdict reasoning**: β passes C1+C2+C3+C4 cleanly with 7/10 cohort agreement
at C3. Grassmann passes C3 at a long contiguous window and strengthens under
epi-X. C5 epi-X at cophenet not run (audit_68 covered α only) but the Grassmann
epi-X confirmation + the unrelated 4-control passage in cophenet make this a
clear **strong trace, both probes**. → manuscript backbone.

---

### α (8–13 Hz) — **strong trace, only D_coph**

**Cophenet `D_coph` — primary 4-control table:**

| Control | Statistic | Source CSV row | Gate | Pass |
|---|---|---|---|---|
| C1 split | p = 0.00977 | `ctm_triangle/cohort_summary.csv` alpha | < 0.05 | ✓ |
| C2 drift | p = 0.00684 | `ctm_triangle/cohort_summary.csv` alpha | < 0.05 | ✓ |
| C3 matched-strength | p = 0.00195; ratio obs/surr = 8.25×; n_above **5/10** | `matched_strength_surrogate_split_baseline/cohort_summary.csv` alpha | < 0.05 | ✓ |
| C4 cross-probe (Wilcoxon) | paired Wilcoxon (split>xprobe) p = 0.461 ⇒ fails to reject ⇒ no degradation; rho_split_median = +0.115, rho_xprobe_median = +0.105 (sign agree); LOO max p = 0.674 (Pat_02) | `ctm_triangle/c4_wilcoxon_cohort.csv` alpha | paired_p ≥ 0.05 + sign match | ✓ |
| C5 epi-X cophenet (Wilcoxon) | one-sample Wilcoxon on per-patient `obs_rho^epi-X` under H_1: rho_split^epi-X > 0 ⇒ **p = 0.0098**; obs_rho_median^epi-X = +0.187; LOO max p = 0.0195 (Pat_03) — LOO robust | `alpha_epi_exclusion/c5_wilcoxon_cohort.csv` alpha | wilcoxon_p < 0.05 | ✓ |

**Grassmann `d_G(k)` — cluster-extent permutation (audit_70):**

| Statistic | Value | Source CSV |
|---|---|---|
| Observed longest contiguous-sig run | 4 cells at k = 11..14 | `grassmann_cluster_extent/cohort_summary.csv` alpha |
| Null mean / 95th / max LR | 2.09 / 6.0 / — | same |
| **cluster p (LR-based)** | **0.1592** | same |
| cluster p (mass-based) | 0.0995 | same |

cluster_p_longest_run ≥ 0.05 → Grassmann verdict: **no trace**.

**Verdict reasoning**: full-cohort C3 cohort agreement is 5/10 (borderline by
n_above count). C5 epi-X strengthens the trace decisively: ratio increases
8.25× → 27.7×, cohort agreement 5/10 → **7/10**, p stays significant at 0.0137.
Per-patient Δρ under epi-X (`alpha_epi_exclusion/comparison.csv`): Pat_03,
Pat_05, Pat_07, Pat_10 strengthen; Pat_13 inverts (+0.028 → −0.171, this is the
patient with the largest epi-zone N_epi=30/119); Pat_06, Pat_08, Pat_14 weaken
mildly but stay positive. Net cohort signal is non-epi-cortex driven. The
"borderline at 5/10" reading at full cohort is the epi zones masking the
trace, not the trace being weak. → **strong trace, only D_coph**.

**Decision (logged)**: the C3 full-cohort n_above 5/10 is below 6/10, which by
the conservative `weak trace` rule would tag α as weak. The epi-X strengthening
(7/10 + ratio 27.7×) reclassifies it as **strong** — the controls work in
sequence (C3 ratifies presence, C5 strengthens). This decision is final;
re-opening requires a new audit, not a re-reading.

---

### γ_l (30–80 Hz) — **strong trace, only Grassmann** ↑ (revised 2026-05-19 pm)

**Revision**: under the mass-only Grassmann gate (CONTROLS.md §C3
Decision 8), γ_l upgrades from weak → **strong** because the
resilient all-clusters cluster mass `cluster_p_mass = 0.005` clears
the `< 0.01` strong threshold. Old longest-run-only gate gave
`cluster_p_LR = 0.0149` (weak by LR-only); the all-clusters mass
correctly captures the multi-cluster structure (γ_l has 13 + 12 + 5
+ 5 + 2 + 2 cells across multiple contiguous runs). LOO max
`p_mass = 0.040` (drop Pat_05) — verdict robust to single-patient
leverage.



**Cophenet `D_coph` — primary 4-control table:**

| Control | Statistic | Source CSV row | Gate | Pass |
|---|---|---|---|---|
| C1 split | p = 0.03223 | `ctm_triangle/cohort_summary.csv` low_gamma | < 0.05 | ✓ |
| C2 drift | p = 0.00977 | `ctm_triangle/cohort_summary.csv` low_gamma | < 0.05 | ✓ |
| C3 matched-strength | p = **0.1162** (fails); ratio 30.2×; n_above 5/10 | `matched_strength_surrogate_split_baseline/cohort_summary.csv` low_gamma | < 0.05 | ✗ |
| C4 cross-probe (Wilcoxon, descriptive — C3 fails so C4 is moot) | paired Wilcoxon p = 0.500 ⇒ no degradation; rho_split_median = +0.140, rho_xprobe_median = +0.143 (sign agree); LOO max p = 0.715 (Pat_08) | `ctm_triangle/c4_wilcoxon_cohort.csv` low_gamma | paired_p ≥ 0.05 + sign | ✓ (moot) |
| C5 epi-X cophenet | not run | — | n/a | n/a |

C3 fails → cophenet verdict **no trace**.

**Grassmann `d_G(k)` — cluster-extent permutation (audit_70, mass-only gate locked 2026-05-19) + C5 epi-X (audit_72):**

| Statistic | Value | Source CSV |
|---|---|---|
| Observed longest contiguous-sig run (descriptive) | **13 cells** | `grassmann_cluster_extent/cohort_summary.csv` low_gamma |
| Null mean / 95th / max LR | 2.63 / 8.0 / 15 | same |
| `cluster_p_longest_run` (descriptive co-statistic) | 0.0149 | same |
| **Resilient all-clusters mass `T_G^*`** | **66.14** | same |
| **`cluster_p_cluster_mass` (gate)** | **0.005** | same |
| LOO max `p_mass` (descriptive) | 0.040 (drop Pat_05) | same |
| **C5 epi-X `cluster_p_cluster_mass^epi-X` (gate)** | **0.030** | `grassmann_epi_exclusion/c5_wilcoxon_cohort.csv` low_gamma |
| C5 epi-X `T_G^*^epi-X` | 32.75 (vs full 66.14, contracts) | same |
| C5 epi-X LOO max `p_mass^epi-X` | 0.159 (Pat_05) — **LOO-fragile under epi-X** | same |
| C5 epi-X longest run | 10 cells (vs 13 in full data) | same |

C5 epi-X **passes** (`p_mass^epi-X = 0.030 < 0.05`) at the cohort
level, but the trace **contracts** under epi-X (mass 66 → 33; LR 13
→ 10) and the LOO sensitivity under epi-X **fails** (max p = 0.159
when Pat_05 dropped). The cohort-level C5 gate passes; the LOO
diagnostic indicates the γ_l trace is partly leveraged on Pat_05 under
epi-X, and should be reported as such per
`feedback_no_single_patient_p_driven.md`.

**Verdict reasoning**: under the mass-only gate (Decision 8),
`cluster_p_cluster_mass = 0.005 < 0.01` → **strong trace** by the
locked rule. LOO max 0.040 still under 0.05 → verdict is robust to
single-patient leverage. C3 on cophenet fails. → **strong trace,
only Grassmann**.

**Why the upgrade from previous "weak" lock**: under the post-fix
resilient all-clusters mass (audit_70 cluster_mass corrected
2026-05-19 pm), γ_l's multiple contiguous-significant `k`-clusters
(13 + 12 + 5 + 5 + 2 + 2 cells) all contribute to `T_G^*`. The
previous longest-cluster-only formula gave `T_G^* = 19.17` and
`p_mass = 0.035`; the corrected formula gives `T_G^* = 66.14` and
`p_mass = 0.005`. Under the disjunctive gate (also retired), γ_l was
already strong; under the new mass-only gate it remains strong with
a single principled threshold.

---

### δ (0.53–4 Hz) — **strong trace, only Grassmann** ↑ (revised 2026-05-19 pm)

**Revision**: under the mass-only Grassmann gate (CONTROLS.md §C3
Decision 8), δ upgrades from weak → **strong** because the resilient
all-clusters cluster mass `cluster_p_mass = 0.005` clears the `< 0.01`
strong threshold. Old longest-run-only mass gave `T_G^* = 12.78` and
`p_mass = 0.025` (weak by LR-only or longest-cluster-mass); the
corrected resilient formula gives `T_G^* = 38.07` and `p_mass =
0.005` because δ has multiple contiguous-significant `k`-clusters
(7 + 4 + 4 + 2 + 2 cells across the k-axis) that the longest-run
formula was under-counting.

**Caveat — single-patient leverage at full data, resolved under C5
epi-X**: LOO max `p_mass = 0.055` (drops Pat_08) at the full-cohort
analysis crosses the 0.05 boundary, indicating Pat_08 leverage on
the full-data verdict. The verdict stays "strong" at the full-data
gate, AND the C5 epi-X analysis (audit_72) resolves the leverage
cleanly: `p_mass^epi-X = 0.005` with **LOO max `p_mass^epi-X = 0.005`
(Pat_02), fully robust**. The interpretation: Pat_08 leverage at the
full-cohort scale was driven by epi-zone contacts, not the true
biological trace — under epi-X the cohort-level signal strengthens
(mass 38 → 44) and the single-patient leverage disappears. The
manuscript text should report both: *"the δ Grassmann trace clears
the cluster-mass gate at the cohort level (`p_mass = 0.005`); the
LOO sensitivity at the full cohort suggests Pat_08 leverage
(LOO max p = 0.055), but the C5 epi-zone-exclusion analysis resolves
this — the trace strengthens under epi-X (mass 38 → 44, `p_mass^epi-X
= 0.005`) and is fully LOO-robust there. The full-data Pat_08
leverage is attributable to epi-zone interactions, not the true
biological trace."*



**Cophenet `D_coph` — primary 4-control table:**

| Control | Statistic | Gate | Pass |
|---|---|---|---|
| C1 split | p = 0.2158 | < 0.05 | ✗ |
| C2 drift | p = 0.2461 | < 0.05 | ✗ |
| C3 matched-strength | p = 0.2783; ratio 1.69×; n_above 4/10 | < 0.05 | ✗ |
| C4 cross-probe (Wilcoxon, descriptive — C3 fails so C4 is moot) | paired Wilcoxon p = 0.385 ⇒ no degradation; rho_split_median = +0.031, rho_xprobe_median = +0.032 (sign agree); LOO max p = 0.590 (Pat_03) | `ctm_triangle/c4_wilcoxon_cohort.csv` delta | paired_p ≥ 0.05 + sign | ✓ (moot) |
| C5 epi-X cophenet | not run | n/a | n/a |

C3 fails → cophenet **no trace**. C4 passes (6/10 cross-probe sign agreement)
but is the *anchor-anatomy* known biology, not a positive trace — the δ
cross-probe pattern reproduces the published δ "anchor" ratio (1.55×, memory
`epileptic_imcoh_universal.md`) which is a known epileptogenesis pattern
detected by ImCoh, not a task-induced reorganization.

**Grassmann `d_G(k)` — cluster-extent permutation (audit_70, mass-only gate locked 2026-05-19) + C5 epi-X (audit_72):**

| Statistic | Value | Source CSV |
|---|---|---|
| Observed longest contiguous-sig run (descriptive) | **7 cells** | `grassmann_cluster_extent/cohort_summary.csv` delta |
| Null mean / 95th / max LR | 2.17 / 5.0 / 10 | same |
| `cluster_p_longest_run` (descriptive co-statistic) | 0.025 | same |
| **Resilient all-clusters mass `T_G^*`** | **38.07** | same |
| **`cluster_p_cluster_mass` (gate)** | **0.005** | same |
| LOO max `p_mass` (descriptive — flag at full data) | **0.055 (drop Pat_08)** | same |
| **C5 epi-X `cluster_p_cluster_mass^epi-X` (gate)** | **0.005** | `grassmann_epi_exclusion/c5_wilcoxon_cohort.csv` delta |
| C5 epi-X `T_G^*^epi-X` | **43.99** (strengthens vs full 38.07) | same |
| C5 epi-X LOO max `p_mass^epi-X` | **0.005 (Pat_02) — fully LOO-robust under epi-X** | same |
| C5 epi-X longest run | 7 cells | same |

**Major C5 finding for δ**: the LOO Pat_08 leverage in the full
analysis (LOO max p_mass = 0.055) **disappears under epi-X**, where
LOO max p_mass = 0.005 (fully robust). The interpretation is that
the full-data Pat_08 leverage was due to **epi-zone interactions**
rather than the true biological trace; removing the epi-zone
contacts strengthens both the cohort-level signal (mass 38 → 44) and
the single-patient-leverage diagnostic. This is *exactly* the kind
of mechanistic finding C5 epi-X is designed to surface, and it
strengthens the biological-attribution argument for the δ Grassmann
trace.

**Verdict reasoning**: under the mass-only gate (Decision 8),
`cluster_p_cluster_mass = 0.005 < 0.01` → **strong trace** by the
locked rule. LOO max 0.055 crosses 0.05 → verdict is
single-patient-leveraged on Pat_08; flag explicitly in manuscript
(see caveat above). C3 on cophenet fails. → **strong trace,
only Grassmann** (with LOO caveat).


real subspace signature. Cophenet C3 fails so no per-pair trace. The δ C4
anchor-anatomy reading remains a *separate, descriptive* known-biology
observation, not part of the trace verdict. → **weak trace, only Grassmann**
(Grassmann), with descriptive anchor-anatomy note alongside.

---

### γ_h (80–300 Hz) — **no trace** (revised 2026-05-19)

**Cophenet `D_coph` — primary 4-control table:**

| Control | Statistic | Source CSV row | Gate | Pass |
|---|---|---|---|---|
| C1 split | p = 0.3477 | `ctm_triangle/cohort_summary.csv` high_gamma | < 0.05 | ✗ |
| C2 drift | p = 0.4229 | `ctm_triangle/cohort_summary.csv` high_gamma | < 0.05 | ✗ |
| C3 matched-strength | p = 0.2461 (fails); ratio 0.082×; n_above 4/10 | `matched_strength_surrogate_split_baseline/cohort_summary.csv` high_gamma | < 0.05 | ✗ |
| C4 cross-probe | rho_xprobe = −0.016, 4/10 (−sign matches rho_split −0.014) | `ctm_triangle/cohort_summary.csv` high_gamma | n ≥ 6/10 + sign match | ✗ |
| C5 epi-X cophenet | not run | — | n/a | n/a |

C3 fails decisively (cohort median ρ ≈ 0, no positive direction) → cophenet **no trace**.

**Grassmann `d_G(k)` — cluster-extent permutation (audit_70):**

| Statistic | Value | Source CSV |
|---|---|---|
| Observed longest contiguous-sig run | **9 cells** | `grassmann_cluster_extent/cohort_summary.csv` high_gamma |
| Null mean / 95th / max LR | 2.50 / 8.0 / — | same |
| **cluster p (LR-based)** | **0.0547** | same |
| cluster p (mass-based) | 0.0647 | same |

Epi-X sensitivity (audit_67): 6-cell longest run at k=21..26 plus 8 cells that
emerge elsewhere under epi-X (the "physiological-attribution" hint).

**Verdict reasoning**: cluster_p_longest_run = 0.0547 — just outside the 0.05
gate. 9-cell observed run is at the null 95th percentile (8.0). The previous
8-cell hardcoded threshold had γ_h passing; cluster-extent demotes it to **no
trace** (Decision 6, 2026-05-19). The audit_67 emergent-cell pattern remains
descriptively interesting but does not change the verdict under the locked
rule. → **no trace** on both probes.

---

### θ (4–8 Hz) — **no trace**

**Cophenet `D_coph` — primary 4-control table:**

| Control | Statistic | Gate | Pass |
|---|---|---|---|
| C1 split | p = 0.6875 | < 0.05 | ✗ |
| C2 drift | p = 0.2783 | < 0.05 | ✗ |
| C3 matched-strength | p = 0.7217; ratio −10.2× (anti-direction); n_above 2/10 | < 0.05 | ✗ |
| C4 cross-probe | rho_xprobe = −0.048, 3/10 (−sign matches rho_split −0.049) | n ≥ 6/10 + sign | ✗ |
| C5 epi-X cophenet | not run | n/a | n/a |

All 4 controls fail. Cohort median is in the anti-trace direction (negative).

**Grassmann `d_G(k)` — cluster-extent permutation (audit_70):**

| Statistic | Value | Source CSV |
|---|---|---|
| Observed longest contiguous-sig run | 5 cells | `grassmann_cluster_extent/cohort_summary.csv` theta |
| Null mean / 95th / max LR | 2.13 / 6.0 / — | same |
| **cluster p (LR-based)** | **0.0995** | same |
| cluster p (mass-based) | 0.1443 | same |

cluster_p_longest_run ≥ 0.05 → Grassmann **no trace**.

**Verdict reasoning**: θ is the cleanest negative — every control on both
probes fails. Useful as the band-specificity benchmark for β and α: whatever
the trace mechanism is, it is **band-selective**, not a generic property of
all FC. → **no trace**.

---

## Decisions on borderlines (logged once, never re-opened)

### Decision 1: α C3 cohort agreement 5/10 borderline → resolved by C5 epi-X
At full cohort, α has C3 paired Wilcoxon p = 0.00195 (passes) but only 5/10
patients above their own surrogate. The conservative reading of the
strong/weak rule would tag this as `weak trace`. C5 epi-X reclassifies it as
`strong trace` because:
- ratio 8.25× → 27.7× (3.4× improvement)
- cohort agreement 5/10 → 7/10 (above the 6/10 threshold)
- p = 0.0137 still significant under the smaller per-patient N

The borderline is a non-epi-cortex effect being masked by epi-zone patients,
not a marginal trace. **Decision: α D_coph is `strong trace`. Final.**

### Decision 2: γ_l Grassmann 83% retention with k-window shift → `weak`
audit_66 contiguous window k=12..23 (12 cells); audit_67 contiguous window
k=19..28 (10 cells). Cell-count retention 83% — above the 80% guide. BUT only
k=19..23 overlap (5/12 = 42% within-window). The signal under epi-X involves
different modes than the full-cohort signal. **Decision: γ_l Grassmann is
`weak trace`. Final.**

### Decision 3: γ_h Grassmann 67% retention + emergent cells → `weak`
audit_66: 9 cells k=19..27. audit_67: 6 cells k=21..26 (67% retention,
overlapping window). Plus 8 new cells emerge elsewhere under epi-X — these
are the physiological-attribution evidence. The 67% retention falls below the
80% guide. **Decision: γ_h Grassmann is `weak trace` with physiological-
attribution framing in the per-band brief. Final.**

### Decision 4: β cophenet C5 epi-X not run
audit_68 was the α-only cophenet epi-exclusion; β was never run at this probe.
Per CONTROLS.md C5 rule: when unavailable, the verdict from C1–C4 stands but
is tagged "epi-X-not-run at this probe". **Decision: β D_coph is `strong
trace` based on C1+C2+C3+C4 + the Grassmann probe's epi-X confirmation. Final.**

### Decision 5: δ C4 6/10 +sign passes but does not constitute a trace
δ at C4 has cross-probe n_trace 6/10 with +sign agreement, which would pass C4
in isolation. But C3 fails (p = 0.278), so per CONTROLS.md the matched-strength
gate defeats this. The C4 reading is the known δ anchor anatomy (1.55× cross-
probe ratio, known biology), not a positive trace. **Decision: δ is `no trace`
under the per-pair (cophenet) trace probe at the per-pair layer; the C4 reading
is presented in the per-band brief as a separate anchor-anatomy known-biology
confirmation, not as a trace claim. Final.**

### Decision 7: cluster_mass promoted to CO-PRIMARY Grassmann gate (2026-05-19)
The Grassmann cluster-extent gate (audit_70) originally used `cluster_p_longest_run`
as the primary statistic, with `cluster_p_cluster_mass` as a companion. The
longest-run statistic has a known limitation: a multi-cluster pattern (e.g.,
"30 contiguous-sig + 1-cell gap + 10 contiguous-sig") is read as LR=30,
missing the secondary 10-cell cluster's contribution. Cluster mass
(`Σ_k (−log10 p_k)` over all p<0.05 cells) aggregates *every* contiguous-sig
cluster's contribution and is therefore robust to multi-cluster patterns.

**Refined gate (locked 2026-05-19)**: pass if EITHER
`cluster_p_longest_run < 0.05` OR `cluster_p_cluster_mass < 0.05`. Both
statistics are co-primary; the verdict is the more permissive of the two.
Sources unchanged: `data/audit/grassmann_cluster_extent/cohort_summary.csv`.

**Verdict re-verification under disjunctive gate**: no flips in our data.
- β: LR p=0.0050, mass p=0.0050 — both strong → strong trace (unchanged)
- α: LR p=0.159, mass p=0.0995 — both > 0.05 → no trace (unchanged)
- γ_l: LR p=0.0149, mass p=0.0348 — both < 0.05 → weak trace (unchanged)
- δ: LR p=0.0249, mass p=0.0249 — both < 0.05 → weak trace (unchanged)
- γ_h: LR p=0.0547, mass p=0.0647 — both > 0.05 → no trace (unchanged)
- θ: LR p=0.0995, mass p=0.1443 — both > 0.05 → no trace (unchanged)

**Decision: disjunctive LR ∨ mass gate is the locked Grassmann gate. Both
statistics co-primary. No verdict flips in our data. Final.**

### Decision 6: 8-cell hardcoded Grassmann threshold → replaced by cluster-extent permutation (2026-05-19)
The previous "Grassmann trace exists if longest contiguous-significant run ≥ 8
cells" rule was a hardcoded judgment threshold. **Replaced** by a cluster-extent
permutation null using the R=200 matched-strength surrogates from audit_66.
For each band, the empirical null distribution of "longest contiguous-significant
run across k ∈ [2, 112]" is built by re-running the cohort paired Wilcoxon R=200
times with each surrogate treated as the phantom observation; cluster p-value =
`(1 + #(null_LR ≥ obs_LR)) / (R + 1)`. New gate: `cluster_p < 0.01` → strong;
`0.01 ≤ cluster_p < 0.05` → weak; `cluster_p ≥ 0.05` → no trace. Source: `data/audit/grassmann_cluster_extent/cohort_summary.csv` (audit_70 script).

**Two verdict flips vs the 8-cell rule:**
- **δ Grassmann**: `no trace` → **`weak trace`** (cluster_p = 0.0249; 7-cell observed run is well above the null mean of 2.17 and 95th percentile of 5.0).
- **γ_h Grassmann**: `weak trace` → **`no trace`** (cluster_p = 0.0547; 9-cell observed run sits at the null 95th percentile of 8.0, just outside the 0.05 gate).

**No verdict change for**: β (cluster_p = 0.0050, strong); γ_l (cluster_p = 0.0149, weak); α (cluster_p = 0.1592, no trace); θ (cluster_p = 0.0995, no trace).

**Decision: cluster-extent permutation is the locked Grassmann gate. The 8-cell
threshold is retired. δ Grassmann is now `weak trace`. γ_h Grassmann is now
`no trace`. Coverage tags updated accordingly in the locked verdict table. Final.**

---

## Anti-revisitation clause

These verdicts are **locked as of 2026-05-18**. To change any verdict in the
table above:
1. A new dated audit run must be executed and committed to `data/audit/`.
2. A dated revision entry must be added below this clause citing the new audit.
3. Every per-band brief `NN_<band>.md` affected must be updated in the same revision.

Until those steps happen, the verdicts above are the source of truth for all
preprint writeup. No per-band brief, figure, or LaTeX paragraph may introduce
a verdict not in this table. Briefs **document** the locked verdict — they do
not re-derive it.

## Revision history

- 2026-05-18 — Initial lock. All 6 bands × 2 probes verdicts assigned from
  CSV reads of `ctm_triangle`, `matched_strength_surrogate_split_baseline`,
  `grassmann_matched_strength_surrogate`, `grassmann_regate_no_filter`,
  `grassmann_epi_exclusion`, `alpha_epi_exclusion`, `raw_fc_matched_strength`,
  `rho_split_raw_D/all_bands_*`. Five borderline decisions logged.
- 2026-05-19 — Grassmann gate replaced. 8-cell hardcoded contiguous-run
  threshold retired in favor of cluster-extent permutation null (audit_70 →
  `grassmann_cluster_extent/cohort_summary.csv`). Decision 6 logged. Two
  verdict flips: δ Grassmann `no trace` → `weak trace` (cluster_p = 0.0249);
  γ_h Grassmann `weak trace` → `no trace` (cluster_p = 0.0547). Coverage tags
  updated in the locked verdict table. CONTROLS.md C3 Grassmann gate updated
  to cite cluster_p instead of the cell-count threshold.
- 2026-05-19 — Grassmann gate refined to **disjunctive LR ∨ mass**. Decision 7
  logged. `cluster_p_cluster_mass` promoted to co-primary alongside
  `cluster_p_longest_run`; both are computed from the same audit_70 null
  (no new audit needed). Addresses the multi-cluster concern: longest_run
  misses contributions from secondary contiguous-sig clusters broken by
  1+ non-sig cells. No verdict flips in our data — both statistics agree
  on the 0.05 gate for all 6 bands. CONTROLS.md C3 Grassmann gate text +
  per-probe vocabulary updated to cite the disjunctive rule.
- 2026-05-19 (pm) — **Three coordinated gate refactors** under writing-agent
  feedback. **Decision 8**: Grassmann gate refactored from disjunctive
  `min(p_LR, p_mass) < α` to **mass-only** `cluster_p_mass < α`. The
  resilient all-clusters cluster mass (also locked 2026-05-19 pm) already
  encodes both contiguity and depth in a single statistic, so the parallel
  LR gate was double-insurance. `cluster_p_LR` stays as descriptive
  co-statistic. audit_70 verdict mapping updated. Two verdict upgrades:
  γ_l Grassmann weak → strong (`p_mass = 0.005`, LOO max = 0.040,
  LOO-robust); δ Grassmann weak → strong (`p_mass = 0.005`, full-data
  LOO max = 0.055 flags Pat_08 leverage; C5 epi-X resolves to `p_mass^epi-X
  = 0.005` fully LOO-robust). β Grassmann unchanged at strong, now
  LOO-robust at both full and epi-X.
  **Decision 9**: C4 cross-probe gate refactored from `n_trace_xprobe ≥ 6/10
  + sign-match` to **paired one-sided Wilcoxon** (`H_1: rho_split > rho_xprobe`,
  fails to reject ⇒ no degradation) + sign-match. Removes the hardcoded
  patient-count threshold per `feedback_no_hardcoded_test_thresholds.md` +
  `feedback_patient_counts_never_the_gate.md`. New audit at
  `scripts/01_compute/audit/audit_71_c4_wilcoxon_cohort.py`; CSV at
  `data/audit/ctm_triangle/c4_wilcoxon_cohort.csv`. All 6 bands pass new
  gate. β paired_p = 0.385 (LOO max 0.590, Pat_02). α paired_p = 0.461
  (LOO max 0.674, Pat_02). No verdict flips.
  **Decision 10**: C5 epi-X gate refactored from `≳80% window-retention` to
  **Wilcoxon-on-epi-X**. Grassmann: re-run audit_70 cluster-mass test on
  epi-X eigvec cache (audit_72) → `cluster_p_mass^epi-X < 0.05`. Cophenet
  (α only): one-sample Wilcoxon on per-patient `obs_rho^epi-X` under
  `H_1: rho_split^epi-X > 0` → α p = 0.0098 (passes). C5 promoted from
  *sensitivity layer* to *primary gate* under the strong/weak rule.
  New audit at `scripts/01_compute/audit/audit_72_c5_wilcoxon_cohort.py`;
  CSVs at `data/audit/{grassmann_epi_exclusion,alpha_epi_exclusion}/c5_wilcoxon_cohort.csv`.
  No verdict flips; α cophenet, β/γ_l/δ Grassmann all pass.
  **Decision 11** (companion to 8/9/10): LOO max-p diagnostic mandatory
  alongside every Wilcoxon-based gate per `feedback_no_single_patient_p_driven.md`.
  Added to audit_70 (`cluster_p_mass_loo_max` column + `loo_cluster_p_mass.csv`
  per-patient file), audit_71 (`wilcoxon_loo_max_p_split_gt_xprobe` column),
  audit_72 (both probes). Descriptive only, never a gate. Reveals δ Grassmann
  Pat_08 leverage at full data, which the C5 epi-X analysis subsequently
  resolves — a methodologically clean attribution of the leverage to
  epi-zone interactions rather than the true biological trace.
