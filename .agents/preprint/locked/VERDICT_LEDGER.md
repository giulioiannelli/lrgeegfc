---
name: preprint-verdict-ledger
era: IMCOH_ABS_COHORT_N10
status: locked_2026-05-18 (revised 2026-05-19 cluster-extent gate; addendum 2026-05-26 β cophenet C5 — Decision 4 resolved, verdict tags unchanged)
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
| δ | 0.53–4 | no trace | **weak trace** (cohort gate clears at `p_mass = 0.005` floor; full-data LOO max p_mass = 0.055 Pat_08 fails Decision-12 < 0.05 LOO precondition — Decision-8 strong promotion retracted; C5 epi-X strengthens and LOO resolves, secondary mechanistic observation per Decision 10) | **weak trace, only Grassmann** |
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

**Anatomy localization is locked separately in [`ANATOMY_LEDGER.md`](ANATOMY_LEDGER.md).**

> ⚠️ **ANATOMY FULLY RETRACTED (2026-05-30) + DELOCALIZED (2026-06-01).** The trace
> verdicts in the table above are **unaffected** and stand — only the *where* is
> retracted. The signed, threshold-free localization audit found **no DK region reaches
> a defensible cohort localization** in any band/probe (max coverage 5/10; locked regions
> rest on 1–4 patients; several anti-localized), and **per-patient localization is also
> null** on both cross-phase probes (cophenet + Grassmann trace: 0–1/10 patients beat
> their own implant-shuffle null in every band; β 0/10 on both). The verified β/α trace
> is **spatially DELOCALIZED** — a distributed network reorganization with no anatomical
> anchor at cohort or single-patient level; the only above-chance spatial structure is
> electrode-shaft autocorrelation (NMI 0.65 region↔shaft), not anatomy. The
> `strong/weak/not localized` tiers and all region lists are **superseded**. Sources:
> `ANATOMY_LEDGER.md`, `data/audit/anatomy_localization_wilcoxon/README.md`,
> `data/audit/per_patient_localization/README.md`.

The original battery: [`ANATOMY_CONTROLS.md`](ANATOMY_CONTROLS.md) (A1 hypergeometric +
A2 sampling-corrected + A3 matched-strength **mandatory** + A4 implant-geometry
regression); KC-era anatomy artifacts at `data/audit/lrg_localization_anatomy/` are
retired and not citable.

---

## Control battery (locked) — full definitions in `CONTROLS.md`

- **C1** baseline-split (D_coph only) — `ctm_triangle/cohort_summary.csv` col `wilcoxon_split_gt_0_p`
- **C2** drift-floor null (D_coph only) — same CSV col `wilcoxon_split_gt_drift_p`
- **C3** matched-strength surrogate (both probes; **mandatory**) — `matched_strength_surrogate_split_baseline/cohort_summary.csv` for D_coph; `grassmann_cluster_extent/cohort_summary.csv` col `cluster_p_cluster_mass` for Grassmann (mass-only gate, Decision 8)
- **C4** cross-probe restriction (D_coph only) — `ctm_triangle/c4_wilcoxon_cohort.csv` (audit_71, locked 2026-05-19; paired-Wilcoxon gate, Decision 9). The old `ctm_triangle/cohort_summary.csv` cols `n_trace_xprobe_int`/`rho_xprobe_median` are retained as **descriptive auxiliary statistics**, not the gate.
- **C5** epi-zone exclusion — `grassmann_epi_exclusion/c5_wilcoxon_cohort.csv` (audit_72, locked 2026-05-19; cluster-mass on epi-X eigvec cache for Grassmann) + `alpha_epi_exclusion/c5_wilcoxon_cohort.csv` (audit_72, locked 2026-05-19; one-sample Wilcoxon on `obs_rho^epi-X` for cophenet α). C5 is a **secondary mechanistic observation** documenting whether epi-zone exclusion strengthens or weakens the trace; **never a verdict-promoter** for either probe (Decision 10 amended 2026-05-28 — the earlier "primary gate" framing was retracted as contradictory to Decision 12, which already locked C5 as secondary for the Grassmann probe; the same rule now applies symmetrically to the cophenet probe).

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
| C5 epi-X cophenet (added 2026-05-26 post-2026-05-20 audit_68_beta) | one-sample Wilcoxon on per-patient `obs_rho^epi-X` under H_1: `rho_split^epi-X > 0` ⇒ **p = 0.003**; obs_rho_median^epi-X = +0.275 (vs full +0.221, strengthens); ratio 26.0× (vs full 23.7×); LOO max p = 0.006 (Pat_02) — LOO robust | `beta_epi_exclusion/c5_wilcoxon_cohort.csv` beta + `cohort_summary.csv` beta | wilcoxon_p < 0.05 | ✓ |

**Grassmann `d_G(k)` — cluster-extent permutation (audit_70, mass-only gate) + C5 epi-X (audit_72):**

| Statistic | Value | Source CSV |
|---|---|---|
| Observed longest contiguous-sig run (descriptive) | **29 cells** at k = 27..55 | `grassmann_cluster_extent/cohort_summary.csv` beta |
| Null mean / 95th / max LR | 2.37 / 6.0 / 17 | same |
| `cluster_p_longest_run` (descriptive co-statistic) | 0.005 (minimum at R=200) | same |
| **Resilient all-clusters mass `T_G^*`** | **69.76** (normalized **0.273** per C1, denominator 255.65) | same |
| **`cluster_p_cluster_mass` (gate)** | **0.005** | same |
| LOO max `p_mass` | 0.005 (Pat_02) — **fully LOO-robust** | same |
| **C5 epi-X `cluster_p_cluster_mass^epi-X` (gate)** | **0.005** | `grassmann_epi_exclusion/c5_wilcoxon_cohort.csv` beta |
| C5 epi-X `T_G^*^epi-X` | **89.04** (normalized **0.444** per C1, denominator 200.38 with n_k^epi-X=87; strengthens vs full 69.76 / 0.273) | same |
| C5 epi-X LOO max `p_mass^epi-X` | 0.005 (Pat_02) — **fully LOO-robust** | same |
| C5 epi-X longest run | 36 cells (vs 29 in full data) | same |

β trace **strengthens** under epi-X (mass 69.76 (0.273) → 89.04 (0.444), LR 29 → 36)
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

**Verdict reasoning**: C3 paired Wilcoxon p = 0.00195 clears α decisively
at the locked Wilcoxon-is-the-gate rule (per
`feedback_patient_counts_never_the_gate.md`). The per-patient `n_above` =
5/10 is a descriptive cohort-agreement statistic, not a verdict modifier;
patient-count thresholds were retired at the 2026-05-19 lock and reaffirmed
2026-05-28. C1, C2, C3, C4 all pass at full cohort. The Grassmann probe
shows no contiguous-significant window at α. → **strong trace, only D_coph**
at C3 alone.

**Secondary mechanistic observation (C5 epi-X, audit_68)**: under epi-zone
exclusion the α effect-size ratio strengthens (obs_rho/surr_p50 8.25× →
27.7×) and the C5 one-sample Wilcoxon on per-patient `obs_rho^epi-X` clears
at p = 0.0098 (LOO max p = 0.0195 Pat_03). Per-patient Δρ under epi-X
(`alpha_epi_exclusion/comparison.csv`): Pat_03/05/07/10 strengthen;
Pat_13 — the patient with the largest epi-zone burden (N_epi = 30/119) —
inverts from +0.028 to −0.171; Pat_06/08/14 weaken mildly but stay positive.
This is a supportive mechanistic narrative consistent with α-relevant
non-epi cortex carrying the trace; it is **not a verdict-promoter**, per
Decision 10 (amended 2026-05-28). The α verdict was already strong at C3
alone.

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
| **Resilient all-clusters mass `T_G^*`** | **66.14** (normalized **0.259** per C1, denominator 255.65) | same |
| **`cluster_p_cluster_mass` (gate)** | **0.005** | same |
| LOO max `p_mass` (descriptive) | 0.040 (drop Pat_05) | same |
| **C5 epi-X `cluster_p_cluster_mass^epi-X` (gate)** | **0.030** | `grassmann_epi_exclusion/c5_wilcoxon_cohort.csv` low_gamma |
| C5 epi-X `T_G^*^epi-X` | 32.75 (normalized **0.163** per C1, denominator 200.38 with n_k^epi-X=87; vs full 66.14 / 0.259, contracts ~50%) | same |
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
`p_mass = 0.035`; the corrected formula gives `T_G^* = 66.14 (normalized 0.259)` and
`p_mass = 0.005`. Under the disjunctive gate (also retired), γ_l was
already strong; under the new mass-only gate it remains strong with
a single principled threshold.

---

### δ (0.53–4 Hz) — **weak trace, only Grassmann** (revised 2026-05-19 pm; Decision-12 cascade 2026-05-28)

**Revision history**: under Decision 6 (cluster-extent revision, 2026-05-19 am)
δ Grassmann was promoted from no trace → weak. Under Decision 8 (mass-only gate,
2026-05-19 pm) δ was tentatively further promoted weak → strong because
`cluster_p_mass = 0.005` cleared the < 0.01 strong threshold by the mass-only
mechanical rule. Under **Decision 12 (LOO + extent preconditions, 2026-05-28)**
the strong promotion is **retracted**: δ stays **weak** because full-data
LOO max `p_mass = 0.055` (drops Pat_08) fails the Decision-12 < 0.05 LOO
robustness precondition. The cohort gate is held — `cluster_p_mass = 0.005`
at the empirical floor — but a single patient (Pat_08) drags the cohort
verdict over the gate, which under brutal-honesty + no-single-patient-p-driven
rules cannot be called "strong". The verdict is therefore **weak trace,
only Grassmann**: cohort gate clears but LOO robustness fails at full data.

**Why the strong promotion was wrong (Decision 12 rationale)**: Decision 8
was a mechanical rule (`cluster_p_mass < 0.01`) without an LOO precondition.
Empirically `p_mass = 0.005` is the empirical-null floor `1/(R+1)` for `R=200`
surrogates — it can't go lower no matter how decisive the cohort signal is.
So the p-value alone cannot discriminate marginal-clear from decisive-clear
signal. The LOO sensitivity is what supplies that discrimination: β passes
LOO 0.005 (Pat_02) and γ_l passes LOO 0.040 (Pat_05), both genuinely
robust to single-patient removal. δ at full-data LOO 0.055 (Pat_08) does
not, so it is not in the same robustness class as β or γ_l. The "strong"
tag would falsely advertise that equivalence.

**C5 epi-X — secondary mechanistic observation (Decision 10, not verdict-promoter)**:
under epi-zone exclusion the trace strengthens decisively — `p_mass^epi-X = 0.005`,
`T_G^*^epi-X = 43.99` (vs full 38.07), and LOO max `p_mass^epi-X = 0.005`
(Pat_02), fully robust. The interpretation is that Pat_08 leverage at the
full-cohort scale was driven by epi-zone interactions, not the true
biological trace. This is **interesting biology to report in Discussion**
but per Decision 10 C5 epi-X is never a verdict-promoter — the full-data
verdict is what gates the strong/weak/no-trace tag. The full-data Pat_08
LOO failure is the binding constraint; C5 resolution is descriptive of
*why* (epi-zone interactions), not a recipe for promotion.

**Manuscript text should report**: *"The δ Grassmann probe clears the
cluster-mass gate at the cohort level (`cluster_p_mass = 0.005`,
`T_G^* = 38.07` raw / 0.149 normalized) but full-data LOO sensitivity
identifies Pat_08 as the leveraging patient (LOO max p = 0.055). Under
the locked verdict-tier rule (Decision 12), this single-patient leverage
prevents promotion to strong tier — the verdict is weak. The C5
epi-zone-exclusion analysis (secondary, mechanistic) shows the trace
strengthens and LOO becomes fully robust under epi-X (mass 38 → 44,
`p_mass^epi-X = 0.005`, LOO under epi-X = 0.005), indicating Pat_08's
full-data leverage is attributable to epi-zone interactions rather than
the true biological trace."*

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
| **Resilient all-clusters mass `T_G^*`** | **38.07** (normalized **0.149** per C1, denominator 255.65) | same |
| **`cluster_p_cluster_mass` (gate)** | **0.005** | same |
| LOO max `p_mass` (descriptive — flag at full data) | **0.055 (drop Pat_08)** | same |
| **C5 epi-X `cluster_p_cluster_mass^epi-X` (gate)** | **0.005** | `grassmann_epi_exclusion/c5_wilcoxon_cohort.csv` delta |
| C5 epi-X `T_G^*^epi-X` | **43.99** (normalized **0.220** per C1, denominator 200.38 with n_k^epi-X=87; strengthens vs full 38.07 / 0.149) | same |
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

**Verdict reasoning (Decision 12 cascade)**: under the Decision-8 mass-only
gate `cluster_p_mass = 0.005 < 0.01` clears the strong-tier p-value threshold,
and the cluster-extent co-statistic `cluster_p_LR = 0.025 < 0.05` clears the
extent threshold. **However, the Decision-12 LOO precondition (full-data
LOO max p_mass < 0.05) fails: LOO max = 0.055 (Pat_08).** Under Decision 12
the strong tier requires all three conditions; δ fails one, so the verdict
is **weak**. C3 on cophenet fails. δ C4 anchor-anatomy reading is a
*separate, descriptive* known-biology observation (LEDGER Decision 5),
not part of the trace verdict. C5 epi-X is a secondary mechanistic
observation per Decision 10 (strengthens trace, resolves LOO under epi-X
to 0.005 Pat_02, fully robust) — interpretable as Pat_08 full-data leverage
being epi-zone-driven, but never a verdict-promoter. → **weak trace,
only Grassmann** (Grassmann), with descriptive anchor-anatomy note
alongside.

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

> **Note on ordering**: decisions appear in file-order for narrative
> coherence (verdict-by-verdict at the top, methodology refinements at the
> bottom), not strict chronology. Each Decision header carries its own date
> stamp; consult those for chronological ordering. Decisions 8 (mass-only
> Grassmann gate), 9 (paired-Wilcoxon C4), and 10 (Wilcoxon-on-epi-X C5)
> live in `CONTROLS.md` rather than this file — they are control-battery
> refinements; this file logs verdict-side decisions. Decision 11 (LOO
> max-p diagnostic mandatory) is companion to 8/9/10 and is described in
> the 2026-05-19 pm revision-history entry below.

### Decision 1: α C3 cohort agreement 5/10 borderline → resolved by C5 epi-X — **RETRACTED 2026-05-28**

**Retraction note (2026-05-28)**: the original Decision 1 used C5 epi-X to
"reclassify" α from weak → strong, anchored on the 5/10 → 7/10 per-patient
count transition and the 8.25× → 27.7× ratio strengthening under epi-X.
This decision is **retracted** for two reasons:

1. **Patient-count thresholds were retired at the 2026-05-19 lock** (per
   `feedback_patient_counts_never_the_gate.md` and CONTROLS.md §C3) — the
   5/10 vs 7/10 distinction was never a principled verdict modifier and
   could not legitimately rescue a tag.
2. **C5 epi-X is a secondary mechanistic observation, not a verdict-
   promoter** — this rule was already locked for the Grassmann probe under
   Decision 12 (2026-05-28) and is now extended symmetrically to the
   cophenet probe under amended Decision 10. C5 cannot reclassify
   anything; it can only document whether epi exclusion strengthens or
   weakens the trace.

**The α verdict tag (`strong trace, only D_coph`) is preserved** because
C3 paired Wilcoxon p = 0.00195 was already passing at full cohort under
the locked Wilcoxon-as-gate rule. The verdict was earned at C3 alone; the
rescue narrative was an unsanctioned overlay. C5 epi-X strengthening is
reported in the α §verdict reasoning above as a supportive mechanistic
observation, not as a promoter.

Original Decision 1 text removed; this retraction note retained for
audit history.

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

### Decision 4: β cophenet C5 epi-X — RESOLVED 2026-05-26 (audit_68_beta_epi_exclusion run 2026-05-20)
On 2026-05-19 the ledger recorded β cophenet C5 epi-X as not run. The β analogue
of `audit_68_alpha_epi_exclusion.py` (`audit_68_beta_epi_exclusion.py`) was
created and run the next day, 2026-05-20 — one day post-lockdown. Outputs live
at `data/audit/beta_epi_exclusion/`. Outcome: cohort median `ρ_split^coph`
strengthens from +0.221 (full) to +0.275 (epi-X), effect-size ratio rises
23.7× → 26.0×, Methods-locked C5 one-sample Wilcoxon on
`ρ_split^epi-X > 0` returns p = 0.003 (LOO max p = 0.006, Pat_02-driven),
n_above_surrogate stays at 7/10, c5_pass = True. Per-patient detail in
`bands/01_beta.md` §3.2.5: 6 of 10 patients strengthen or flip positive under
epi-X (notably Pat_10 inverts from −0.091 to +0.029, the only LRG-anti
patient at full β); Pat_02 retains a strong but halved effect (+0.507 →
+0.266), consistent with Pat_02 being the LOO argmax-p driver — a substantial
part of Pat_02's full β trace lived in epi-zone-coupled topology, and the
cohort verdict still passes after that contribution is stripped. **Decision:
β D_coph is `strong trace` based on C1+C2+C3+C4+C5 cophenet (audit_68_beta)
+ the Grassmann probe's C5 confirmation. The β verdict tag is unchanged
(`strong trace, both probes`); the C5 cophenet line now reinforces rather
than is owed. Final.**

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

### Decision 12: LOO + extent preconditions for strong-tier verdict (2026-05-28)

Decision 8 (mass-only gate, 2026-05-19 pm) was a mechanical rule:
`cluster_p_mass < 0.01` → strong; `0.01 ≤ cluster_p_mass < 0.05` → weak;
`cluster_p_mass ≥ 0.05` → no trace. Empirically `cluster_p_mass = 0.005` is
the empirical-null floor `1/(R+1)` for `R=200` surrogates, so the p-value
alone cannot discriminate marginal-clear from decisive-clear signal once
the floor is reached. The LR null and LOO sensitivity supply the missing
discrimination.

**Refined rule (locked under Decision 12):** strong-tier verdict requires
**all three** conditions at full data:
1. `cluster_p_mass < 0.01` (Decision-8 cluster-mass-null clearance)
2. `cluster_p_longest_run < 0.05` (cluster-extent-null clearance; weak
   threshold suffices since mass is the gate)
3. **`LOO max p_mass < 0.05` (full data) — cohort verdict robust to
   single-patient removal**

C5 epi-X resolutions are **secondary mechanistic observations** (per
Decision 10) and **never** promote a band to strong if conditions (1)–(3)
fail at full data.

**Per-band check under Decision 12** (data from
`grassmann_cluster_extent/cohort_summary.csv`, 2026-05-26):

| Band | `mass_p < 0.01` | `LR_p < 0.05` | LOO < 0.05 | Verdict |
|---|---|---|---|---|
| β | ✓ 0.005 | ✓ 0.005 | ✓ 0.005 (Pat_02) | **strong** |
| γ_l | ✓ 0.005 | ✓ 0.015 | ✓ 0.040 (Pat_05) | **strong** ↑ |
| **δ** | ✓ 0.005 | ✓ 0.025 | **✗ 0.055 (Pat_08)** | **weak** ← LOO binds |
| γ_h | ✗ 0.060 | — | — | no trace |
| θ | ✗ 0.159 | — | — | no trace |
| α | ✗ 0.348 | — | — | no trace |

**One verdict flip vs Decision 8:**
- **δ Grassmann**: tentative strong (Decision 8) → **weak** (Decision 12).
  Cohort gate clears at `p_mass = 0.005` floor and `LR_p = 0.025 < 0.05`,
  but full-data LOO max `p_mass = 0.055 (Pat_08)` fails the < 0.05 LOO
  precondition. Pat_08 single-handedly leverages the cohort verdict over
  the gate at full data. C5 epi-X strengthens the trace and resolves the
  LOO (mass 38.07 → 43.99, LOO under epi-X = 0.005 Pat_02) — secondary
  mechanistic observation interpretable as Pat_08 full-data leverage
  being epi-zone-driven, but never a verdict-promoter per Decision 10.

**No verdict change for**: β (LOO 0.005 robust); γ_l (LOO 0.040 robust);
α (cophenet probe, unaffected by Grassmann gate); γ_h, θ (already no trace
under Decision 8).

**Decision: LOO + extent preconditions are locked. δ Grassmann is now
`weak trace`. Coverage tags updated accordingly in the locked verdict
table. Final.**

Rationale ties: `feedback_no_single_patient_p_driven.md` (Wilcoxon at n=10
vulnerable to direction-outliers; LOO mandatory) + `feedback_brutal_honesty_no_sycophancy.md`
(don't overclaim "strong" when robustness checks reveal fragility) +
`feedback_no_hardcoded_test_thresholds.md` (LOO max is principled robustness
on the same statistical test, not a hardcoded patient-count threshold —
the test IS the gate, LOO sensitivity is part of the test).

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
  **Decision 10**: C5 epi-X test refactored from `≳80% window-retention` to
  **Wilcoxon-on-epi-X**. Grassmann: re-run audit_70 cluster-mass test on
  epi-X eigvec cache (audit_72) → `cluster_p_mass^epi-X < 0.05`. Cophenet
  (α only): one-sample Wilcoxon on per-patient `obs_rho^epi-X` under
  `H_1: rho_split^epi-X > 0` → α p = 0.0098.
  **Amended 2026-05-28**: the original Decision 10 promoted C5 from
  *sensitivity layer* to *primary gate*. This promotion is **retracted**.
  C5 epi-X is a **secondary mechanistic observation** for both probes
  (cophenet and Grassmann) — it documents whether epi-zone exclusion
  strengthens or weakens the trace but **never promotes a band to strong**
  if the C1–C4 primary controls do not already pass at full data. This
  brings Decision 10 into alignment with Decision 12 (locked 2026-05-28
  for the Grassmann probe), and triggers retraction of Decision 1
  (α cophenet C5 rescue). Under the amended Decision 10, the α verdict
  is anchored by C3 paired Wilcoxon p = 0.00195 alone; Decision 12 governs
  the Grassmann strong/weak gate.
  New audit at `scripts/01_compute/audit/audit_72_c5_wilcoxon_cohort.py`;
  CSVs at `data/audit/{grassmann_epi_exclusion,alpha_epi_exclusion}/c5_wilcoxon_cohort.csv`.
  No verdict flips: α cophenet stays strong (anchored at C3); β/γ_l/δ
  Grassmann verdicts governed by Decision 12.
  **Decision 11** (companion to 8/9/10): LOO max-p diagnostic mandatory
  alongside every Wilcoxon-based gate per `feedback_no_single_patient_p_driven.md`.
  Added to audit_70 (`cluster_p_mass_loo_max` column + `loo_cluster_p_mass.csv`
  per-patient file), audit_71 (`wilcoxon_loo_max_p_split_gt_xprobe` column),
  audit_72 (both probes). Descriptive only, never a gate. Reveals δ Grassmann
  Pat_08 leverage at full data, which the C5 epi-X analysis subsequently
  resolves — a methodologically clean attribution of the leverage to
  epi-zone interactions rather than the true biological trace.
- **2026-05-30 → 2026-06-01 — ANATOMY localization RETRACTED + trace shown
  DELOCALIZED. No trace-verdict flips (the table above stands).** A signed,
  threshold-free localization audit (`diag_anatomy_localization_wilcoxon.py`)
  retracted all 6 locked anatomy cells: no DK region reaches a defensible cohort
  localization (max coverage 5/10; locked regions rest on 1–4 patients; several
  anti-localized; Hip crosses on only 3–5/10 = marginal hint). The per-patient
  test (`diag_per_patient_localization.py`) found localization **also null at the
  single-patient level** on both cross-phase probes — cophenet AND Grassmann trace
  (`s_i = (p^task−p^pre)(p^post−p^pre)`): 0–1/10 patients beat their own
  implant-shuffle null in every band; β 0/10 on both; holds at lobe granularity.
  An apparent between-region η² "concentration" is electrode-shaft spatial
  autocorrelation (anatomy-free shaft partition reproduces it; region↔shaft
  NMI 0.65), not anatomy. **Verdict: the verified β/α trace is spatially
  DELOCALIZED — distributed reorganization, no anatomical anchor at cohort or
  single-patient level.** Verified by workflow `wf_ddbfbe43-da7` (5 adversarial
  lenses + shaft control). Reports: `data/audit/anatomy_localization_wilcoxon/README.md`,
  `data/audit/per_patient_localization/README.md`. Anatomy `strong/weak/not
  localized` tiers + region lists superseded; see `ANATOMY_LEDGER.md`.
