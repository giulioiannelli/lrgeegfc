---
name: preprint-control-battery
era: IMCOH_ABS_COHORT_N10
status: locked_2026-05-18
kind: control-battery-lockdown
supersedes: any_prior_implicit_control_definitions
companion: VERDICT_LEDGER.md
---

# Control battery for the LRG trace probes (locked 2026-05-18)

**Head.** Four primary controls + two sensitivity layers (C5 epi-zone exclusion,
C6 white-matter exclusion) constitute the entire control battery for the preprint
trace claims. The two sensitivity layers are mechanistic/robustness observations,
not verdict gates. No claim is upgraded above what these controls support.
Verdicts are tagged `strong trace` / `weak trace` / `no trace` per probe
(cophenet `ρ_split^coph` on `D_coph`; Grassmann `d_G(k)` on `U_k`) and combined
into a per-band coverage tag (both probes / only `D_coph` / only Grassmann / none).

## The four primary controls

### C1 — Baseline split (`ρ_split^coph > 0`)
Constructs the per-pair trace by correlating Δ_task = `D_coph^task − D_coph^rsPre_A`
with Δ_rest = `D_coph^rsPost − D_coph^rsPre_B` across all `N(N−1)/2` contact pairs,
with **independent** rsPre halves. Addresses common-reference contamination in rsPre.
- **Gate**: paired one-sided Wilcoxon `wilcoxon_split_gt_0_p < 0.05` against zero.
- **Source**: `data/audit/ctm_triangle/cohort_summary.csv` col `wilcoxon_split_gt_0_p`.
- **Applies to**: `ρ_split^coph` (D_coph probe).
- **Not applicable to** Grassmann (no split-baseline construction).

### C2 — Drift-floor null (`ρ_split > ρ_drift`)
Rest-only counterpart of `ρ_split`. `ρ_drift = ρ_S(D_coph^rsPre_B − D_coph^rsPre_A,
D_coph^rsPost_B − D_coph^rsPost_A)`. Captures within-session drift on the same
halved-data noise budget as `ρ_split`.
- **Gate**: paired one-sided Wilcoxon `wilcoxon_split_gt_drift_p < 0.05`.
- **Source**: `data/audit/ctm_triangle/cohort_summary.csv` col `wilcoxon_split_gt_drift_p`.
- **Applies to**: `ρ_split^coph` (D_coph probe).
- **Not applicable to** Grassmann (no within-rest drift Grassmann analog audited).

### C3 — Matched-strength surrogate (`obs > surrogate`)
Strength-preserving 4-cycle ±δ rewiring of each adjacency matrix; preserves node
strength `s_i` to 10⁻⁶. R=200 surrogates, swap_target=20, seed=20260511. Mandatory
per `feedback_matched_strength_mandatory.md`.
- **Gate (cophenet)**: paired one-sided Wilcoxon `paired_wilcoxon_p < 0.05`. Source: `data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv` col `paired_wilcoxon_p`. LOO diagnostic: `wilcoxon_loo_max_p` column to be added per `feedback_no_single_patient_p_driven.md`.
- **Gate (Grassmann)**: **cluster-extent permutation, mass-only** (audit_70, refined 2026-05-19 Decision 8). For each band, build the empirical null distribution of the **resilient all-clusters cluster mass** `T_G^*(b) = Σ_{k : p_k < α_k} (−log_10 p_k)` across all contiguous-significant `k`-cells (`α_k = 0.05`). R=200 phantom-surrogate tests (each surrogate `r` treated as the observation against the remaining R−1). Empirical cluster p-value `cluster_p_mass(b) = (1 + #{null ≥ obs}) / (R + 1)`. **Gate**: `cluster_p_mass(b) < 0.05`. The resilient all-clusters mass encodes both contiguity (long contiguous runs contribute many `−log_10 p_k` terms with consistent depth → large mass) and depth (deep cells contribute exponentially more than marginal cells) in a single statistic, so the parallel `cluster_p_LR` gate is double-insurance against a failure mode the mass formula has already absorbed (see writing-agent feedback 2026-05-19). The longest-run length `L_obs` and its empirical `cluster_p_LR` stay in the CSV as a **descriptive co-statistic** — useful to read alongside `T_G^*` to characterise whether the trace concentrates in one contiguous window or distributes across multiple — but **do not gate** the verdict. Source: `data/audit/grassmann_cluster_extent/cohort_summary.csv` col `cluster_p_cluster_mass` (gate) + `cluster_p_longest_run` (descriptive). LOO diagnostic: `cluster_p_mass_loo_max` + `cluster_p_mass_loo_argmax_patient` columns added 2026-05-19 per `feedback_no_single_patient_p_driven.md`.
- **Applies to**: both probes. This is the *only* primary control for Grassmann.

### C4 — Cross-probe restriction (paired Wilcoxon non-degradation + sign agreement)
`ρ_split^coph` recomputed on cross-probe pairs only (drops same-probe pairs);
denote per-patient `ρ_xprobe`. The gate asks whether the cross-probe
restriction *does not significantly degrade* the per-pair correlation
**and** whether the cross-probe cohort median preserves the sign of the
full-pair cohort median.
- **Gate**: paired one-sided Wilcoxon (`rho_split − rho_xprobe`) under
  `H_1: rho_split > rho_xprobe`, **fails to reject** at `α = 0.05`
  (i.e., `paired_wilcoxon_p ≥ 0.05` ⇒ no significant degradation),
  **AND** `sign(rho_xprobe_median) == sign(rho_split_median)`.
- **Source**: `data/audit/ctm_triangle/c4_wilcoxon_cohort.csv`
  (computed 2026-05-19 from per-patient `rho_split` + `rho_xprobe` in
  `data/audit/ctm_triangle/Td_per_patient_per_band.csv`). Cols
  `paired_wilcoxon_p_split_gt_xprobe`, `rho_split_median`,
  `rho_xprobe_median`, `sign_agreement`, `c4_pass`,
  `wilcoxon_loo_max_p_split_gt_xprobe`, `wilcoxon_loo_argmax_patient`.
- **Applies to**: `ρ_split^coph` (D_coph probe).
- **Not applicable to** Grassmann (subspaces are on all N nodes, not pair-level).
- **Replaces** (2026-05-19 Decision 9, per writing-agent feedback): the
  earlier `n_trace_xprobe ≥ 6/10 AND sign-match` gate. Patient-count
  thresholds (`n_trace ≥ N/10`) are arbitrary and couple the verdict to
  cohort size in a non-portable way; the paired Wilcoxon answers
  "does the per-patient distribution sit on the trace side, accounting
  for sign and magnitude?" without any hardcoded count cutoff.
  See `feedback_no_hardcoded_test_thresholds.md`.

## The sensitivity layers

### C5 — Epi-zone exclusion (paired Wilcoxon on per-patient trace under epi-X)
Recompute the probe with epileptogenic-zone contacts removed per patient
(`load_epileptic_nodes` × `channel_labels.csv`); cached surrogate eigvecs at
`data/cache/matched_strength_surrogate_epi_excluded_lrg/Pat_NN/{band}_{phase}_epiX_R200_swap20_seed20260514_imcoh_abs.npz`.
The gate asks whether the trace direction survives epi-zone removal —
biological-attribution layer.
- **Gate (Grassmann)**: re-run the audit_70 cluster-extent permutation
  test on the epi-X surrogate ensemble, producing
  `cluster_p_mass^epi-X(b)`. C5 passes iff `cluster_p_mass^epi-X < 0.05`
  (i.e., the resilient all-clusters cluster mass under epi-X clears
  the same gate as the full-data analysis). Same statistic as C3 for
  Grassmann, applied to the epi-X-restricted graph.
- **Gate (cophenet, α only)**: one-sample one-sided Wilcoxon on
  per-patient `rho_split^epi-X` under `H_1: rho_split^epi-X > 0`
  (trace direction at cophenet); C5 passes iff
  `wilcoxon_one_sided_p < 0.05`. Where C5 is not available at the
  cophenet probe for other bands (no audit run), the verdict is
  tagged "epi-X-not-run at this probe" in VERDICT_LEDGER.md and the
  verdict from C1–C4 stands.
- **Source for Grassmann (all bands)**:
  `data/audit/grassmann_epi_exclusion/c5_wilcoxon_cohort.csv`
  (computed 2026-05-19; cluster_p_mass^epi-X + LOO max columns).
- **Source for cophenet (α only)**:
  `data/audit/alpha_epi_exclusion/c5_wilcoxon_cohort.csv`
  (computed 2026-05-19; one-sample Wilcoxon p + LOO max columns).
- **LOO diagnostic**: `wilcoxon_loo_max_p_epiX` + `loo_argmax_patient`
  columns in both CSVs, descriptive only (per
  `feedback_no_single_patient_p_driven.md`).
- **Replaces** (2026-05-19 Decision 10, per writing-agent feedback):
  the earlier "≳80% retention of the original significant window"
  rule. The 80% number is arbitrary and outside the statistical
  framework; replacing with a Wilcoxon-on-epi-X gate gives a
  principled answer to the same question ("does the trace direction
  survive without the epi zone?") with no hardcoded threshold.
  See `feedback_no_hardcoded_test_thresholds.md`.
- **Node-count control on the cophenet "strengthens under epi-X" claim
  (2026-06-12, `audit_85_wm_decimation_control.py --stratify epi` — NB distinct
  from the same-numbered `audit_85_epi_propagator_recovery.py`).** The 2026-06-05
  PI directive elevated epi-X to a
  PRIMARY interpretive lens partly on "the cophenet trace *strengthens* when the
  diseased zone is excluded". A size-matched random-node-decimation control shows
  **that subgraph strengthening is NOT epi-specific** — it is a generic
  node-count effect. epi-X drops only ~9.5 % of nodes (median), yet removing the
  same number of *random* nodes reproduces the cohort ρ_split: cophenet
  `cohort_p_dec` = α 0.135, β 0.390, γ_l 0.130 (all `generic_nodecount`, none
  < 0.05; `data/audit/epi_stratified/decimation_control_cohort.csv`). Raw δ/θ are
  `epi_carried` (removing epi *hurts* — consistent with δ as the epi channel).
  **What this means:** the *subgraph* `exclude_epi` strengthening/emergence
  (β/α strengthen, γ_l emerge) is demoted to "not epi-specific". The interpretation
  is **re-anchored on the same-graph PAIR-CLASS decomposition** (no node-count
  confound).
- **PI resolution + pair-count control (2026-06-18, `audit_85 --mode pairclass`).**
  PI decision: epi-exclusion is **demoted from the 2026-06-05 "primary interpretive
  lens" to a SECONDARY control** (mechanistic / tissue-characterization layer,
  parallel to C6) — cohort-wide (α/β/γ_l); the genuine pair-class findings are
  retained as reported secondary characterizations. The pair classes were re-run at
  the current **n=10** (Pat_15 now SOZ-labelled; `audit_77` regenerated — no
  matched-strength flag flips vs n=9) and put through a **pair-count control** (is a
  class more trace-bearing than a random subset of the same number of pairs?).
  Resolved β reading (stated on the cophenetic object `D_coph`; raw |ImCoh| is the
  comparison baseline only, never a result): **the β cophenetic trace is gray-
  cortical-dominant** — the genuine pair-count hotspot is gray↔gray coupling
  (`p_pair`<0.001; wm↔wm depleted; C6 gray/WM split). The **diseased seizure core
  does NOT carry the trace** — epi↔epi fails its matched-strength null (`weaken`,
  +0.142 p=0.16) and is not a pair-count hotspot; it *trends* depleted (cophenet
  `p_pair`=0.945) but does **not** clear the 0.95 gate, so "core spared" is
  directional only, **not** a significant result. **A second withdrawal:** the
  "epi↔non-epi CROSS interface is the strongest carrier" sub-claim does NOT survive
  the pair-count control (interface not a hotspot: cophenet `p_pair`=0.093). α
  `epi_epi` is genuinely `concentrated` (pair-count `p_pair`<0.001): **α RECRUITS the
  diseased core whereas β spares it** — a reported secondary tissue-characterization
  (the α-specific contrast), not part of the β headline. Cascaded into
  `bands/01_beta.md` (§3.2.5 + Head + frontmatter). Source
  `data/audit/epi_stratified/{cophenetic_cohort.csv, pairclass_decimation_cohort.csv}`;
  handoff `.agents/reports/2026-06-12_decimation-controls-handoff.md`. The IDENTICAL
  node-count confound was found and resolved for C6/WM the same day (2026-06-12).

### C6 — White-matter exclusion (paired Wilcoxon on per-patient trace under exclude_wm)
Recompute the probe with **dominant-white-matter contacts removed** per patient
(`load_channel_regions(pat)["region"] == "Wm"`, atlas argmax of Desikan-Killiany
tissue weights, PTD excluded); cached surrogate eigvecs at
`data/cache/matched_strength_surrogate_wm_excluded_lrg/Pat_NN/{band}_{phase}_wmX_R200_swap20_seed20260608_imcoh_abs.npz`.
WM is **30–57 % of every montage** (cohort ≈40 %; median 117→74 nodes under
`exclude_wm`), so this is a structural cut, not a fine perturbation. Because
ImCoh is **pairwise**, the gray-only FC equals the gray-only submatrix of the
cached FC to numerical precision — `exclude_wm` is therefore the exact "drop WM
channels before computing FC" test; only the node-coupled LRG step changes. The
gate asks whether the trace direction survives a gray-only montage —
recording-substrate-attribution layer (the WM analogue of C5's
biological-attribution layer).
- **Status — SECONDARY, mechanistic (not a verdict gate).** Like the amended C5
  (Decision 10, 2026-05-28), C6 is a mechanistic/robustness observation
  documenting whether the trace strengthens, persists, or weakens when
  white-matter contacts are removed. It is **NOT a primary verdict gate** and
  changes **no** locked verdict. It is **not** elevated to a primary interpretive
  lens — the node-count caveat below precludes that. (Epi-exclusion's 2026-06-05
  elevation to a primary lens was likewise demoted to secondary on 2026-06-18, so
  C5 and C6 now sit at the same secondary tier.)
- **Gate (cophenet)**: cohort one-sided paired Wilcoxon of per-patient
  `ρ_split^coph` (exclude_wm) vs own matched-strength surrogate median, regenerated
  on the gray-only submatrix; passes iff `paired_wilcoxon_p < 0.05`. LOO line:
  leave-Pat_15-out p (Pat_15 has the most WM, 57 %, and is the β-LRG anti-aligned
  patient).
- **Gate (Grassmann) — now on the locked C3 gate (audit_86, 2026-06-12).** The
  C6-Grassmann verdict is the audit_70 cluster-extent permutation **mass** gate
  (`cluster_p_mass^wmX < 0.01` strong / `< 0.05` weak) re-run on the `exclude_wm`
  per-`k` T_G curves — the SAME gate as the locked C3 Grassmann verdict.
  **β and low-γ re-pass the gate** (`cluster_p_mass^wmX = 0.005` both, strong,
  LOO-robust); **δ fails** (0.105, no_trace — the WM-dependent cross-probe
  epi-biology channel, LEDGER Decision 5); γ_h emerges (0.005) but LOO-fragile
  (0.050); α/θ absent both. The earlier per-`k` Wilcoxon count (β 40→40 etc.) is
  retained as a co-statistic in `grassmann_cohort.csv`, but the gate verdict is
  the cluster mass. Source: `data/audit/wm_stratified/grassmann_cluster_extent_wmX.csv`.
- **Source (cophenet + raw)**: `data/audit/wm_stratified/cophenetic_raw_cohort.csv`
  (cols `paired_wilcoxon_p`, `lopat15_wilcoxon_p`, `obs_median`, `sensitivity_flag`,
  per (substrate ∈ {cophenetic, raw}, band, config ∈ {full, exclude_wm, wm_only,
  gray_gray, cross, wm_wm})). **Source (Grassmann)**:
  `data/audit/wm_stratified/grassmann_cohort.csv`. Full report:
  `.agents/reports/2026-06-08_white-matter-exclusion.md`. Audits: `audit_83`
  (cophenetic + raw), `audit_84` (Grassmann per-`k`), `audit_85`
  (random-node-decimation control → `decimation_control_cohort.csv`, flag A),
  `audit_86` (Grassmann C3 cluster-extent gate → `grassmann_cluster_extent_wmX.csv`,
  flag B).
- **What it shows (2026-06-08, amended 2026-06-12).** WM removal does **not**
  weaken the trace — the trace **survives** a gray-only montage: cophenet α/β
  still clear matched-strength on the gray-only submatrix (α p 0.005→0.014,
  LO-P15 0.027; β obs +0.22→+0.32, full-cohort p=0.053 marginal but LO-P15 0.037,
  `gray_gray` p=0.032, raw-β p=0.024 — a small-n power crossing, preserved not
  lost), and **β/low-γ Grassmann re-pass the locked C3 cluster-extent gate**
  (audit_86, `cluster_p_mass^wmX`=0.005 both, strong, LOO-robust). The **only**
  WM-dependent trace component is **δ Grassmann** (fails the gate, 0.105; per-`k`
  23→7) — the cross-probe epi-biology channel (LEDGER Decision 5), not the task
  trace.
  **What we do NOT claim (amended 2026-06-12, flag A resolved by audit_85).** The
  apparent cophenetic *sharpening* (α/β rise, γ_l/θ "emerge") is **not
  WM-specific** — a size-matched random-node-decimation control reproduces it by
  removing *any* ~40 % of nodes (cohort `p_dec`: β 0.145, α 0.305, γ_l 0.185 —
  none < 0.05 ⇒ `generic_nodecount`). So WM is **neither carrying nor diluting**
  the cophenetic trace (it is not WM_carried either: `p_dec`≁1). The earlier
  "WM removal sharpens the trace / gray-resident-because-it-sharpens" framing is
  **withdrawn** for the cophenetic substrate; the load-bearing C6 claim is
  **survival/robustness**, not sharpening. The **one** genuinely WM-specific
  effect is the **raw |ImCoh| substrate**: raw β/γ_l clear matched-strength only
  after WM removal and ARE WM-specific (decimation `p_dec` 0.000/0.030,
  `WM_specific`) — the raw substrate's weakness genuinely was WM strength-
  structure. `wm_only` carries weaker traces of its own, so the trace is
  gray-matter-**dominant**, not gray-matter-**exclusive**.

## Verdict vocabulary (locked)

### Per-probe (D_coph and Grassmann independently)

For **D_coph** (gate via Wilcoxon p):
- **strong trace** — C3 paired Wilcoxon p < 0.05 AND every applicable PRIMARY control passes its own Wilcoxon-based gate (C1 split-vs-zero p<0.05; C2 paired-split-vs-drift p<0.05; C4 paired non-degradation under epi-X-style cross-probe restriction passes per §C4 above).
- **weak trace** — C3 passes BUT one applicable primary control's Wilcoxon-based gate fails. The verdict text records which control failed.
- **no trace** — C3 Wilcoxon p ≥ 0.05. All other controls become moot.
- **C5 epi-X is SECONDARY** (amended 2026-05-28): C5 is a mechanistic observation documenting whether epi-zone exclusion strengthens or weakens the trace; it is NOT a primary verdict gate for either probe. The earlier "C5 promoted to primary gate" framing (Decision 10, 2026-05-19) is retracted in favor of symmetric treatment with the Grassmann probe under Decision 12 — see VERDICT_LEDGER.md retraction notes for Decision 1 + amended Decision 10.
- **Cohort-agreement language retired** (2026-05-19): the earlier "borderline ≈5/10 patients" downgrade rule is removed because patient-count thresholds are arbitrary; the relevant question — *"does the per-patient distribution sit on the trace side?"* — is answered by the C3 paired Wilcoxon already. LOO max-p per `feedback_no_single_patient_p_driven.md` is the descriptive sensitivity layer for single-patient leverage.

For **Grassmann** (gate via cluster-extent permutation, audit_70; **mass-only** as of 2026-05-19 Decision 8):
- **strong trace** — `cluster_p_cluster_mass < 0.01`. Where C5 epi-X is available, `cluster_p_cluster_mass^epi-X < 0.05` (i.e., the resilient all-clusters cluster mass under epi-X also clears the gate).
- **weak trace** — `0.01 ≤ cluster_p_cluster_mass < 0.05`. The all-clusters mass distinguishes the observed signature from the empirical null but not at strong strength.
- **no trace** — `cluster_p_cluster_mass ≥ 0.05`. The all-clusters cluster mass is not distinguishable from the empirical null.

`cluster_p_longest_run` is reported alongside the verdict as a descriptive co-statistic (does the trace concentrate in one contiguous window or spread across multiple?) but does **not** gate the verdict. The mass-only gate replaces the earlier disjunctive `min(p_LR, p_mass) < α` rule because the resilient all-clusters cluster mass already encodes contiguity (via the sum-of-many-deep-`p_k` terms) and depth (via `−log_10 p_k` per cell) in a single statistic — adding `p_LR` as a second disjunctive arm enlarges the rejection region without principled basis (writing-agent feedback 2026-05-19).

### Per-band coverage tag

Combines the two probes:
- **strong trace, both probes** (D_coph + Grassmann both strong)
- **strong trace, only D_coph** (strong on D_coph, no trace on Grassmann)
- **strong trace, only Grassmann** (strong on Grassmann, no trace on D_coph)
- **weak trace, both probes** (weak on at least one, no strong anywhere)
- **weak trace, only D_coph** (weak on D_coph, no trace on Grassmann)
- **weak trace, only Grassmann** (weak on Grassmann, no trace on D_coph)
- **no trace** (no trace on both probes)

Each band gets exactly one verdict from this list. Verdicts are written into
`VERDICT_LEDGER.md` once and not changed silently.

## Anti-revisitation

A verdict in `VERDICT_LEDGER.md` cannot be changed by reading a different CSV or
running a new sensitivity test. Re-evaluation requires:
1. A new dated audit run reflected in `data/audit/<probe>/cohort_summary.csv`.
2. A dated revision entry in `VERDICT_LEDGER.md` citing the new audit.
3. Cascading updates in the per-band brief `NN_<band>.md`.

Sensitivity tests outside the locked battery (coverage-matched permutation,
symmetric cross-baseline at LRG layer, τ-sweep) are descriptive only and do
not change the verdict. (Pat_03 dropout was previously listed here as a
sensitivity test; **retired 2026-05-18** — Pat_03 is now a full cohort
member treated identically to every other patient at the analysis layer,
with sampling-rate handled at the config layer only.)
