---
name: cohort-coverage-matrix
type: scope
era: COHORT_N10
status: draft
created: 2026-04-29
updated: 2026-04-29
pointers:
  - .agents/guides/task-persistence-investigation/README.md
  - .agents/guides/task-persistence-investigation/2026-04-29_decision-rules.md
  - .agents/guides/task-persistence-investigation/2026-04-29_measure-correctness-audit.md
  - .agents/reports/2026-04-27_task-persistence-reconciliation.md
---

# Cohort-coverage matrix — schema, artefacts, and figure specification

**The data-flow / artefact side of the rebuild's central deliverable. Defines
the per-cell CSV (one row per `patient × band × rung × scale_axis × scale_value`),
the aggregated per-(band × rung) CSV that materialises `V(b, r)`, the per-band
triangulation CSV that materialises `T(b)`, and the band × rung heatmap PDF
that supersedes `task_trace_band_k_n10_imcoh_abs.pdf` as the headline figure.
Pairs strictly with `2026-04-29_decision-rules.md` (the predicate side):
this scope says *what artefacts to write and which inputs to read*, the
decision-rules scope says *which thresholds and triangulation logic to apply*.
Together they specify the producer (`audit_24_cohort_coverage_matrix.py`)
end-to-end. Both must be committed before the producer runs.**

---

## Notation

- `P = {Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08, Pat_10, Pat_13, Pat_14, Pat_15}` — cohort, `|P| = 10`. Source: `data/audit/cohort_metadata.csv`.
- `B = {delta, theta, alpha, beta, low_gamma, high_gamma}` — six bands as stored in input CSVs (Greek glyphs in figures and prose, ASCII in CSV columns).
- `R = {L1, L3, L4, L4_aux, L5_k, L5_hrel, L6, L7}` — eight rung-IDs (seven ladder rungs from `2026-04-29_decision-rules.md` §1 plus `L4_aux` for the Cohesion-CBR side-by-side row; L0 closed dead; **L2 (entropy curve / specific heat) permanently removed** — `S(τ)` and `C(τ)` are non-informative for our continuous-spectrum outlier case, see `lrg-framework-guide.md` §6 + memory `lrg_outlier_case_fully_connected.md`). The CSV `rung` column uses these exact strings.
- `s` — scale value along the rung's scale axis. For rungs without a scale axis (`L1, L3, L4, L4_aux, L6`), `scale_axis = "none"` and `scale_value = NaN`. For scale-axis rungs:
  - `L5_k`: `scale_axis = "k"`, `scale_value ∈ {2, 3, …, 119}` (integer, partition cardinality from h2_partition_multiscale).
  - `L5_hrel`: `scale_axis = "h_rel"`, `scale_value ∈ [0.05, 0.95]` (float; from `dvi_hrel_n10_imcoh_abs.csv` rows; `linspace` grid is primary, `logspace` is sensitivity).
  - `L7`: `scale_axis = "k"`, `scale_value ∈ {2, …, 119}` (integer, MSPC k axis).
  - `L3`: `scale_axis = "lambda"`, `scale_value ∈ {0.0, 0.1, …, 1.0}` (KC λ-mix, owed by §7.2 scope).
- For per-leaf rungs (`L4, L4_aux, L6, L7`): the per-cell value `Q_r(p, b [, s])` is a **count of T-leaves** `n_T(p, b [, k])` per the §2 schema, *not* a per-leaf row. Per-leaf evidence lives in the producer CSVs (CNP/MSPC/Cohesion); the matrix carries only the cohort-level count.
- `audit_verdict ∈ {current, stale-numeric, broken, duplicate, archive}` — from `data/audit/measure_audit/measure_audit_n10_imcoh_abs.csv` (Phase 0 deliverable). The matrix only ingests rows from measures with `audit_verdict == "current"`.

---

## Definitions

### 1. The per-cell CSV — `cohort_coverage_matrix_n10_imcoh_abs.csv`

One row per `(patient, band, rung, scale_axis, scale_value)`. Columns:

| Column | Type | Domain | Description |
|---|---|---|---|
| `patient` | str | `P` | cohort patient ID |
| `band` | str | `B` | one of six band names (ASCII) |
| `rung` | str | `R` | ladder rung ID per Notation |
| `scale_axis` | str | `{none, k, h_rel, lambda}` | scale axis name (or `none`) |
| `scale_value` | float | varies | scale coordinate (`NaN` if `scale_axis = none`) |
| `measure_value` | float | rung-specific | `Q_r(p, b [, s])` from the producer artefact |
| `null_value` | float | rung-specific | `Q_r^{null}(p, b [, s])` from the rung-matched null artefact (`NaN` if no null available) |
| `passes_null` | bool | `{0, 1, NaN}` | `1` iff `measure_value > null_value` in the trace direction; `NaN` if `null_value is NaN` |
| `passes_zero` | bool | `{0, 1, NaN}` | `1` iff `measure_value > 0` for continuous statistics; `NaN` for count-style rungs (L4, L4_aux, L6, L7) where the threshold is `≥ τ_leafcount` rather than `> 0` |
| `passes_count` | bool | `{0, 1, NaN}` | `1` iff `n_T-leaves(p, b [, k]) ≥ τ_leafcount(b [, k])` for count-style rungs; `NaN` for continuous-statistic rungs |
| `eligible` | bool | `{0, 1}` | `1` iff cell is geometrically evaluable (h_rel ≥ cohort `h_rel_min` for L5_hrel; `n_leaves` ≥ k for L5_k / L7; etc. — see §3 caveats) |
| `pi` | str | `{0, 1, ineligible, uncontrolled}` | per-cell positivity verdict per `2026-04-29_decision-rules.md` §2.1 |
| `dissenter_marker` | str | `{none, x, +}` | `x` for Pat_02 (anatomy-dissenter), `+` for Pat_03 (1024 Hz outlier), `none` for the rest |
| `producer_artefact` | str | path | source CSV or NPZ (relative to repo root) |
| `null_artefact` | str | path | source CSV/NPZ for the null (`""` if no null) |
| `audit_verdict` | str | `{current, …}` | from Phase 0 audit (must be `current` for ingestion) |
| `notes` | str | free text | any per-cell flag (e.g. `"k > n_leaves"`, `"h_rel below dmin"`, `"null missing"`) |

**Row count expectation (v1, before §7 controls land):**
- `L1`: `10 × 6 = 60` rows (no scale axis).
- `L4`: `10 × 6 = 60` rows (k-bin scope flattened to per-cell `n_T`).
- `L4_aux`: `10 × 6 = 60` rows.
- `L5_k`: `10 × 6 × |k_grid|` rows (k_grid 2 … min(119, n_leaves)). At max ≈ `10 × 6 × 118 = 7080` rows.
- `L5_hrel`: `10 × 6 × |h_rel_grid|` rows (`linspace` primary, ~60 points = ~3600 rows).
- `L6`: `10 × 6 = 60` rows (per-cell `n_T-leaves`).
- `L7`: `10 × 6 × |k_grid|` rows.
- `L3`: empty in v1 (no producer yet); rows added in matrix v2 once §7.2 lands.
- **Total v1 cells:** O(15k) rows. Compressible to a Parquet sibling if read-time becomes a bottleneck.

### 2. The aggregated per-(band × rung) CSV — `band_verdict_n10_imcoh_abs.csv`

One row per `(band, rung, scale_value)` for scale-axis rungs and one row per `(band, rung)` for non-scale rungs. Columns:

| Column | Type | Description |
|---|---|---|
| `band` | str | one of six band names |
| `rung` | str | ladder rung ID |
| `scale_axis` | str | `{none, k, h_rel, lambda}` |
| `scale_value` | float | `NaN` if no scale axis |
| `n_pos` | int | `|{p ∈ P : pi(p, b, r [, s]) = 1}|` |
| `n_eligible` | int | `|{p ∈ P : pi != "ineligible"}|` |
| `frac_pos` | float | `n_pos / n_eligible` (`NaN` if `n_eligible = 0`) |
| `frac_pos_robust` | float | `n_pos / n_eligible` after dropping Pat_02 + Pat_03 (sub-claim, never headline) |
| `wilcoxon_z` | float | one-sided `Q_full > Q_null` Z statistic (per band, paired) |
| `wilcoxon_p` | float | one-sided p-value |
| `wilcoxon_q` | float | BH-FDR adjusted within rung, m=6 bands |
| `ridge_length` | int | longest contiguous run of `frac_pos ≥ 0.8` along scale axis (`NaN` if no scale axis) |
| `V` | str | `{positive, negative, silent, ineligible, uncontrolled}` per `2026-04-29_decision-rules.md` §2.3 |
| `n_total_cells` | int | `|P|` (always 10) |
| `notes` | str | per-cell flags collapsed (e.g. `"3 cells flagged k > n_leaves"`) |

For scale-axis rungs, `V(b, r)` is the *band-level* verdict computed as: positive iff `cohort_threshold ∧ wilcoxon_q ≤ 0.05 ∧ ridge_length ≥ 3`. The per-`s` rows carry per-cell `frac_pos` for figure rendering; the band-level aggregate appears as the per-(b, r) row with `scale_value = NaN` and the ridge / wilcoxon collapsed across `s`.

### 3. The triangulation CSV — `triangulation_n10_imcoh_abs.csv`

One row per band. Columns:

| Column | Type | Description |
|---|---|---|
| `band` | str | one of six band names |
| `V_L1, V_L3, V_L4, V_L4_aux, V_L5_k, V_L5_hrel, V_L6, V_L7` | str | `{positive, negative, silent, ineligible, uncontrolled}` |
| `T` | str | `{headline-triangulated, partition-resolution-locked, continuous-only, anatomy-suspect, ergodic, uncontrolled-exploratory, inconsistent}` per `2026-04-29_decision-rules.md` §2.4 |
| `T_pending_controls` | str | comma-separated list of rungs missing controls in v1 (e.g. `"L3,L7"`) — empty in v2 once §7 controls land |
| `dissenter_sensitivity` | str | `{stable, dissenter-driven, dissenter-suppressed}` from comparing strict vs robust frac_pos across all rungs |
| `narrative_pointer` | str | path to per-band paragraph in `2026-04-29_task-trace-rebuild-verdict.md` (rebuild plan §10 step 19) |

### 4. The headline figure — `cohort_coverage_matrix_n10_imcoh_abs.pdf`

Band × rung heatmap. Specification in §"Visualisation spec".

---

## Consumer list — which input each rung pulls from

This table is the **read contract**. The producer (`audit_24_cohort_coverage_matrix.py`) reads only these paths, only after Phase 0 verdict is `current` for the producing measure.

| Rung | `measure_value` source | `null_value` source | Cell key in source | Notes |
|---|---|---|---|---|
| `L1` | `data/reports/imcoh_continuous_trace/per_cell_summary_split.csv` (`rho` column = `ρ_split` = Run A) | `data/reports/imcoh_continuous_trace/controls_summary.csv` (`rho_null_drift` column = Run C drift floor) | `(patient, band)` | both sources are 60-row CSVs; join on `(patient, band)` |
| `L3` | **NOT YET** — KC λ-sweep, owed by `2026-04-29_kc-lambda-sweep.md` (§7.2) | `d_KC(λ; T_rpre_A, T_rpre_B)` | `(patient, band, λ)` | v1: emit `pi = uncontrolled` for all L3 cells. (L2 entropy-curve rung permanently removed — `S(τ)`/`C(τ)` non-informative for our continuous-spectrum case.) |
| `L4` | `data/audit/trace_modules/trace_subtrees_n10_imcoh_abs.csv` (count `regime == "T"` rows per `(patient, band)`) | **NOT YET** — half-baseline T-count, owed by `2026-04-29_strict-identity-half-baseline.md` (§7.3) | `(patient, band)` | v1: `pi = 1` iff `n_T ≥ 1` (no null comparison); flagged `uncontrolled` until null lands |
| `L4_aux` | `data/audit/per_patient_hierarchy_cohesion/leaf_assignment.csv` (count `dominant == "T"` rows per `(patient, band)`) | **NOT YET** — half-baseline cohesion T-count | `(patient, band)` | v1: same caveat as L4 |
| `L5_k` | `data/reports/imcoh_vi/h2_partition_multiscale_raw.csv` (`d_VI` column) | `data/audit/dvi_split_baseline/dvi_split_baseline_n10_imcoh_abs.csv` (`drift_dVI` column) | `(patient, band, k)` | join on `(patient, band, k)`; filter `eligible` by checking `k ≤ n_leaves(patient, band, phase)` from `dmin_dmax_inventory.csv` |
| `L5_hrel` | `data/audit/dvi_split_baseline/dvi_hrel_n10_imcoh_abs.csv` (`d_VI` column, `grid == "linspace"` primary) | **NOT YET** — h_rel drift floor, owed by `2026-04-29_dvi-hrel-drift-floor.md` (§7.4) | `(patient, band, h_rel)` | v1: emit `pi = uncontrolled`; eligibility = `h_rel ≥ h_rel_min(patient, band)` from dmin inventory |
| `L6` | `data/audit/per_patient_hierarchy_cnp/leaf_assignment.csv` (count `dominant == "T"` rows per `(patient, band)`) | **NOT YET** — per-leaf half-baseline ρ, owed by `2026-04-29_cnp-half-baseline.md` (§7.5) | `(patient, band)` | v1: emit `pi = uncontrolled` for L6 cells; `n_T-leaves` reported as `measure_value` |
| `L7` | `data/audit/per_patient_hierarchy_mspc/multiscale_assignment.csv` (count `dominant_at_k == "T"` rows per `(patient, band, k)`) | **NOT YET** — per-(leaf, k) half-baseline J, owed by `2026-04-29_mspc-half-baseline.md` (§7.6) | `(patient, band, k)` | v1: emit `pi = uncontrolled`; eligibility = `k ≤ n_leaves(patient, band, phase)` |
| **dissenter flags** | `data/audit/cohort_metadata.csv` (`is_core_N5`, `notes` columns) | n/a | `patient` | Pat_02 → `x`, Pat_03 → `+` |
| **eligibility gate (k, h_rel)** | `data/audit/dvi_split_baseline/dmin_dmax_inventory.csv` (`n_leaves, dmin, dmax, h_rel_min` columns) | n/a | `(patient, band, phase)` | use `n_seconds_*` columns from cohort_metadata to confirm phase coverage |
| **audit verdicts** | `data/audit/measure_audit/measure_audit_n10_imcoh_abs.csv` (`audit_verdict` column) | n/a | `measure_id` | preflight check: refuse to ingest any input whose row is not `audit_verdict == "current"` |

**Out-of-scope inputs (explicitly NOT consumed by v1):**
- `data/audit/four_phase_taumin/`, `data/audit/four_phase_tau_star/` — closed by `2026-04-29_measure-correctness-audit.md` (audit verdict `archive` for the τ-star cohort claim).
- `data/audit/psi_tau_scan/` — closed by `2026-04-28_psi-tau-scan-verdict.md` (Ψ-as-selector dead).
- `data/audit/residual_subspace/` — closed by `2026-04-28_residual-subspace-diagnostic.md` (rank-1 trivial on dense FC).

---

## Properties

- **Per-cell evaluation is a pure function of (input CSV row, decision-rules scope, audit verdict).** Reproducible bit-identical given fixed inputs; no RNG without explicit seed.
- **Joinability is the load-bearing schema property.** Every per-cell CSV row carries the (patient, band, rung, scale_value) tuple; downstream filters (drop ineligibles, drop dissenters, restrict to controlled rungs) are one-line pandas operations.
- **Append-only schema across versions.** Matrix v2 (after §7 controls land) adds rows for L3 and replaces `null_value = NaN` with computed nulls for L4/L4_aux/L5_hrel/L6/L7. Existing v1 rows must remain interpretable when v2 lands; only the `pi`, `passes_null`, `null_value`, `null_artefact` fields change for these cells. No column is removed; new columns may be appended.
- **`uncontrolled` is a first-class verdict.** It is *not* the same as `silent`. A cell can be `pi = uncontrolled` (no null available, value reported), `pi = ineligible` (geometrically degenerate), `pi = 0` (controlled and below null), or `pi = 1` (controlled and above null). Triangulation predicates per `2026-04-29_decision-rules.md` §2.4 distinguish these explicitly.
- **Long-form CSV preferred over wide.** Wide layout (one column per rung) breaks scale-axis rungs and forces sentinel encoding for missing nulls. Long form is one row per cell; pivoting for figure rendering is trivial.
- **Schema validates against `2026-04-29_decision-rules.md` §2 contract.** Each `pi` value is computable from this CSV's `(measure_value, null_value, eligible)` triplet by the decision-rules positivity test for that rung. Violation is a producer bug.

---

## Caveats and failure modes

| Caveat | Mitigation |
|---|---|
| **Sparse vs dense rungs.** L4 / L4_aux / L6 producer CSVs are *sparse* (only T-rows recorded; absent `(patient, band)` ⇒ `n_T = 0`). L5_k / L5_hrel are dense (every cell present). | Producer must left-join the `cohort_metadata × bands × scales` Cartesian product against the producer CSV; absent rows fill with `n_T = 0` for L4 / L4_aux / L6. Failure ⇒ silent under-counting. |
| **Two h_rel grids in `dvi_hrel_n10_imcoh_abs.csv` (linspace, logspace).** Decision-rules pre-registers `linspace` as primary. | v1 matrix ingests only `grid == "linspace"`. `logspace` runs feed a sensitivity row in `2026-04-29_task-trace-rebuild-verdict.md` (rebuild plan §10 step 19), not the headline matrix. |
| **k-axis degeneracy near singletons / giant.** `Δ_VI(k)` collapses for `k ≥ n_leaves(patient, band, rest_post_phase)`. | Eligibility gate: cells with `k > n_leaves(rest_post)` ⇒ `eligible = 0, pi = "ineligible"`. The `n_leaves` per phase comes from `dmin_dmax_inventory.csv`. |
| **MSPC table is huge (~407k rows).** Materialising in memory may stress smaller machines. | Producer streams MSPC: groupby `(patient, band, k)` with chunked `pd.read_csv(chunksize=...)`; emits the per-(p, b, k) `n_T-leaves` count without materialising the full per-leaf assignment. |
| **`trace_subtrees_n10_imcoh_abs.csv` is currently 1 row** (Pat_02 alpha only at the moment). | This is correct behaviour (sparse + null cohort-wide); producer left-joins against full Cartesian product and emits `n_T = 0` for the other 59 cells. Audit reports `n_T = 0` for all but Pat_02 alpha. |
| **L1 `rho_split` and `rho_null_drift` come from different CSVs (per_cell_summary_split.csv, controls_summary.csv).** Schema mismatch ⇒ silent join failure. | Schema check: both CSVs have `(patient, band)` as a unique key with `|P| × |B| = 60` rows each. Producer asserts row counts match before join; aborts on mismatch. |
| **Pat_03 has 1024 Hz sampling and recomputes use `nperseg = 2048`** — its `n_leaves` and FC matrix may differ in shape from the n=10 default. | Eligibility computed per-(patient, band, phase) from `dmin_dmax_inventory.csv`, which is sampling-rate-correct. No additional handling needed. Pat_03 is flagged `+` always. |
| **Pat_02 is a strong contributor to L5_k δ ridge AND a known anatomy dissenter.** Robust-w/o-Pat_02 verdicts can flip cells from `positive` to `silent`. | `frac_pos_robust` reported in band_verdict CSV alongside `frac_pos`. Headline always uses `frac_pos`. The `dissenter_sensitivity` column in triangulation CSV records whether `T(b)` is stable under dissenter removal. |
| **Audit verdict `current` from Phase 0 is bound to the audit run timestamp** (`audit_run_metadata.json`). If a producer is re-run after Phase 0, the audit verdict is stale. | Producer reads `audit_run_metadata.json` and asserts that all consumed CSVs have `mtime ≤ audit_timestamp`. Stale ⇒ refuse to run; user must re-run audit_21 first. |
| **`triangulation_n10_imcoh_abs.csv` v1 will have `T_pending_controls` non-empty for most bands** (L3, L4_null, L5_hrel_null, L6_null, L7_null missing). | Expected. The `T_pending_controls` column is the *map of what controls v2 needs to add*. v1 verdicts are valid as "best-effort with current control coverage"; v2 verdicts (after §7) are the publication headline. |
| **The visualisation may suggest false certainty.** Heatmap colour ⇒ implies a cohort statistic; in v1 several rungs are uncontrolled. | Figure must render uncontrolled cells with a hatched overlay (not just colour) and explicitly annotate "uncontrolled — sensitivity analysis only". The figure caption names every uncontrolled rung. |

---

## Pseudocode

```
INPUT  : decision_rules_thresholds  # from 2026-04-29_decision-rules.md §2
       : measure_audit_csv           # data/audit/measure_audit/measure_audit_n10_imcoh_abs.csv
       : cohort_metadata             # data/audit/cohort_metadata.csv
       : dmin_dmax_inventory         # data/audit/dvi_split_baseline/dmin_dmax_inventory.csv
       : per-rung input CSVs (per Consumer list)
OUTPUT : cohort_coverage_matrix_n10_imcoh_abs.csv
       : band_verdict_n10_imcoh_abs.csv
       : triangulation_n10_imcoh_abs.csv
       : cohort_coverage_matrix_n10_imcoh_abs.pdf

# Step 0: Preflight
audit_df = read(measure_audit_csv)
assert all(audit_df.audit_verdict == "current") for measures in CONSUMER_LIST
audit_ts = read(audit_run_metadata.json).timestamp
for path in CONSUMER_LIST.values():
    assert mtime(path) <= audit_ts, f"stale: re-run audit_21"

# Step 1: Build the Cartesian skeleton
patients = read(cohort_metadata).patient_id   # 10
bands    = ["delta","theta","alpha","beta","low_gamma","high_gamma"]
rungs    = ["L1","L3","L4","L4_aux","L5_k","L5_hrel","L6","L7"]
scales   = build_scale_grid(rungs, dmin_dmax_inventory)
skeleton = cartesian(patients, bands, rungs, scales)   # ~15k rows

# Step 2: Populate measure_value / null_value per rung
matrix = []
for rung in rungs:
    if rung == "L1":
        m = read(per_cell_summary_split.csv).rename(rho -> measure_value)
        n = read(controls_summary.csv).rename(rho_null_drift -> null_value)
        matrix.extend(join(skeleton[L1], m, n, on=(patient, band)))
    elif rung == "L5_k":
        m = read(h2_partition_multiscale_raw.csv).rename(d_VI -> measure_value)
        n = read(dvi_split_baseline_n10_imcoh_abs.csv).rename(drift_dVI -> null_value)
        matrix.extend(join(skeleton[L5_k], m, n, on=(patient, band, k)))
    elif rung in {"L4","L4_aux","L6"}:
        m = read(producer_csv).filter(dominant=="T" or regime=="T").groupby(patient,band).size()
        matrix.extend(left_join(skeleton[rung], m, fillna=0))   # sparse
        # null_value = NaN until §7 controls land
    elif rung == "L7":
        for chunk in stream(mspc_assignment.csv, chunksize=100_000):
            n_T = chunk[chunk.dominant_at_k=="T"].groupby(patient,band,k).size()
            matrix.extend(left_join(skeleton[L7][chunk_keys], n_T, fillna=0))
    elif rung == "L5_hrel":
        m = read(dvi_hrel_n10_imcoh_abs.csv).filter(grid=="linspace").rename(d_VI -> measure_value)
        matrix.extend(join(skeleton[L5_hrel], m, on=(patient, band, h_rel)))
    elif rung == "L3":
        # v1: emit skeleton rows with measure_value=NaN, pi=uncontrolled
        matrix.extend(skeleton[rung].assign(measure_value=NaN, pi="uncontrolled"))

# Step 3: Eligibility per cell
for row in matrix:
    if row.scale_axis == "k":
        n_lv = dmin_dmax_inventory[(row.patient,row.band,"rest_post")].n_leaves
        row.eligible = (row.scale_value <= n_lv)
    elif row.scale_axis == "h_rel":
        h_min = dmin_dmax_inventory[(row.patient,row.band,"rest_post")].h_rel_min
        row.eligible = (row.scale_value >= h_min)
    else:
        row.eligible = True

# Step 4: Per-cell positivity per decision-rules §2.1
for row in matrix:
    if not row.eligible:
        row.pi = "ineligible"; continue
    if isnan(row.null_value):
        row.pi = "uncontrolled"; continue
    if row.rung in {"L1","L5_k","L5_hrel"}:  # continuous
        passes_null = row.measure_value > row.null_value
        passes_zero = row.measure_value > 0
        row.pi = 1 if (passes_null and passes_zero) else 0
    elif row.rung in {"L4","L4_aux","L6","L7"}:  # count
        tau_lc = leafcount_threshold(row.band [, row.scale_value])  # 95th pct of half-baseline
        row.pi = 1 if (row.measure_value >= tau_lc) else 0

# Step 5: Dissenter markers
for row in matrix:
    row.dissenter_marker = {"Pat_02":"x","Pat_03":"+"}.get(row.patient, "none")

# Step 6: Emit per-cell CSV
write(matrix, cohort_coverage_matrix_n10_imcoh_abs.csv)

# Step 7: Aggregate to per-(band, rung [, scale]) — band_verdict CSV
# (consumes decision-rules.md §2.3 directly)
band_verdict = aggregate(matrix, decision_rules_thresholds)
write(band_verdict, band_verdict_n10_imcoh_abs.csv)

# Step 8: Triangulation per band
triangulation = triangulate(band_verdict, decision_rules.md §2.4)
write(triangulation, triangulation_n10_imcoh_abs.csv)

# Step 9: Render figure
render_heatmap(band_verdict, triangulation, cohort_coverage_matrix_n10_imcoh_abs.pdf)
```

---

## Visualisation spec

The figure `data/outputs/figures/section6/cohort_coverage_matrix_n10_imcoh_abs.pdf` renders the matrix as a band × rung heatmap with a triangulation strip on the right.

**Layout**
- Single PDF page, landscape, 1 main heatmap axis + 1 right-side annotation strip.
- Rows: 6 bands top-to-bottom in canonical order (`delta, theta, alpha, beta, low_gamma, high_gamma`), labels in Greek glyphs (`δ, θ, α, β, γ_l, γ_h`).
- Columns: 8 rungs left-to-right in geometric order (`L1, L3, L4, L4_aux, L5_k, L5_hrel, L6, L7`). L2 is permanently absent.
- Each cell's main fill is the band-aggregated `frac_pos(b, r)` on a viridis colormap normalized `[0, 1]`. For scale-axis rungs (`L5_k, L5_hrel, L7`), the displayed `frac_pos` is the *maximum over the eligible scale range* of `frac_pos(b, r, s)`; the per-`s` ridge is shown in a smaller per-band sparkline strip below the row (optional v2 addition; v1 just shows the max).

**Annotations inside cells**
- Cell text: `"n_pos/n_eligible"` as a small fraction (e.g. `9/10`).
- Dissenter glyphs (`x` for Pat_02, `+` for Pat_03) appear in the corner if the cell is positive *only when* the dissenter is dropped (i.e. `frac_pos < 0.8` strict but `frac_pos_robust ≥ 0.8`). Glyphs are *informational*; they do not change the strict verdict.
- Border:
  - **black solid 2pt**: `V(b, r) = positive`
  - **red dashed 2pt**: `V(b, r) = negative`
  - **thin grey 0.5pt**: `V(b, r) = silent`
  - **none**: `V(b, r) = ineligible`
- Hatched overlay (forward-slash `\\\\`): `V(b, r) = uncontrolled` (v1: most cells outside L1, L5_k will be hatched). Hatch is *additive* on top of fill colour and border.

**Right-side triangulation strip**
- One vertical strip per band, color-coded by `T(b)`:
  - green = headline-triangulated
  - yellow = partition-resolution-locked
  - blue = continuous-only
  - orange = anatomy-suspect
  - grey = ergodic
  - purple = uncontrolled-exploratory
  - red = inconsistent
- Strip width 1 cell, height 1 row. Label inside the strip: `T(b)` short code (`HT, PRL, CO, AS, ER, UE, IN`).
- Below each strip in tiny font: `T_pending_controls` (e.g. `"pending: L3, L7"`) — only if non-empty.

**Caption**
- Mandatory caption (rendered as figure-text below the axes, NOT as `fig.suptitle`):
  > Cohort-coverage matrix at n=10 under `imcoh_abs`. Rows = bands; columns = LRG geometric ladder rungs (L1 = ultrametric pair shifts; L5_k = integer-k partition; L5_hrel = fractional-depth partition; L4 = strict subtree identity; L6 = per-leaf cophenetic vector; L7 = per-(leaf, k) cluster-mate sets; L4_aux = Cohesion-CBR; L3 v1 uncontrolled; L2 entropy/specific-heat permanently removed for the continuous-spectrum outlier case). Cell colour: cohort fraction `n_pos/n_eligible` (≥ 8/10 ⇒ positive). Border: solid black = positive, red dashed = negative, grey = silent. Hatch: uncontrolled (no rung-matched null). `x` = Pat_02 (anatomy-dissenter); `+` = Pat_03 (1024 Hz outlier) — both counted in denominator. Right strip: per-band triangulation verdict `T(b)`. v1 reflects current control coverage; v2 (after rebuild §7) replaces hatched cells with controlled verdicts.

**Style invariants (per `CLAUDE.md`)**
- PDF only; no PNG sibling.
- Heatmap full vector (no `set_rasterized`).
- No `fig.suptitle`.
- ≥ 3 bands × 3 rungs visible (always satisfied: 6 × 9 grid).
- Use `lrg_eegfc.config.const.BRAIN_BAND_TEX_DICT` for band glyph rendering.
- Colourbar: vertical, right of main axis (before the triangulation strip), labeled `n_pos / n_eligible`.

---

## Connection to prior tools

This scope is *strictly downstream* — it specifies how to materialise the join of all upstream artefacts into one canonical surface.

| Prior tool | Relationship |
|---|---|
| `2026-04-29_decision-rules.md` (predicate side) | **Paired-PR sibling.** Decision-rules specifies *how* to label a cell (`pi`, `V`, `T`); this scope specifies *what cells exist and where their values come from*. Both must be committed before audit_24 runs. |
| `2026-04-29_measure-correctness-audit.md` (Phase 0) | **Hard precondition.** This scope's preflight refuses to ingest any input whose Phase 0 verdict is not `current`. The `audit_verdict` column carries this through into per-cell rows; downstream consumers of the matrix can filter on it. |
| `2026-04-26_continuous-trace-matrix.md` | The L1 row of this matrix consumes the three-control-consolidation outputs from that scope (`per_cell_summary_split.csv` for `ρ_split`, `controls_summary.csv` for `ρ_null_drift`). The continuous-trace matrix figure remains valid for L1-only deep dive; this matrix supersedes it as the headline. |
| `2026-04-25_task-trace-canonical.md` (P/T/R/RA leaf-set regimes) | The L4 / L4_aux / L6 / L7 rungs of this matrix all aggregate the canonical T-regime predicate. The canonical defines the leaf-level regime; this matrix computes the cohort-level count of T-leaves and gates by `τ_leafcount`. |
| `2026-04-27_task-persistence-reconciliation.md` §8.8 plain-English verdicts | **Cross-check target.** Matrix v1 (rebuild plan §10 step 5) must reproduce §8.8 verdicts when restricted to currently-controlled rungs (L1, L5_k, L4). Mismatch ⇒ matrix bug, not new finding. After §7 controls land, v2 supersedes §8.8. |
| `2026-04-25_measure-ledger.csv` | Inputs to this matrix are exactly the rows where `ledger.status == "current"` *AND* `audit_verdict == "current"`. Ledger is the catalogue; audit is the validation; matrix is the consumer. |
| `task_trace_band_k_n10_imcoh_abs.pdf` (current headline figure) | **Will be retired by rebuild plan §10 step 18** in favour of this matrix figure. Until step 18, both coexist; era-map.md flags the band-k figure as superseded once the matrix figure has cohort approval. |
| MRL / Cohesion-CBR / CNP / MSPC scope reports | This matrix consumes their per-leaf CSVs but does NOT replace them. Per-patient deep-dive figures continue to render from the per-leaf CSVs; the matrix is the cohort-level census. |

---

## Implementation plan

**Producer:** `scripts/01_compute/audit/audit_24_cohort_coverage_matrix.py` (rebuild plan §10 step 3 — named `audit_24` to avoid collision with the existing `audit_18_taustar_four_phase.py`; the rebuild plan's "audit_18" label refers to the same script under its previous numbering). Single CLI:

```
lrg-eegfc audit cohort-coverage-matrix --cohort n10 --fc-method imcoh_abs -v
# OR direct:
python scripts/01_compute/audit/audit_24_cohort_coverage_matrix.py --cohort n10 --fc-method imcoh_abs -v
```

Flags:
- `--cohort {n10}` — single cohort string (defensively rejects `n9`).
- `--fc-method {imcoh_abs}` — single method string (defensively rejects others).
- `--h_rel-grid {linspace, logspace}` — default `linspace`; sensitivity runs use `logspace`.
- `--include-uncontrolled / --no-include-uncontrolled` — default include; useful for v1.
- `-v` verbose logging.

**Outputs** (paths absolute from repo root):
- `data/audit/cohort_coverage_matrix/cohort_coverage_matrix_n10_imcoh_abs.csv`
- `data/audit/cohort_coverage_matrix/band_verdict_n10_imcoh_abs.csv`
- `data/audit/cohort_coverage_matrix/triangulation_n10_imcoh_abs.csv`
- `data/audit/cohort_coverage_matrix/audit_24_run_metadata.json` — git SHA, audit_21 timestamp, decision-rules SHA, input mtimes, scale grids used, dissenter map.
- `data/outputs/figures/section6/cohort_coverage_matrix_n10_imcoh_abs.pdf` (rebuild plan §10 step 4 — separate PR if cleanly separable; otherwise bundle).

**Library entrypoints (must reuse — promotion-to-library justification):**
- `lrg_eegfc.utils.metrics.hypothesis.{wilcoxon_z, bh_fdr, rank_biserial, boot_ci_mean}` — for per-(b, r) Wilcoxon and BH-FDR within rung.
- `lrg_eegfc.utils.metrics.tree.{simpson_neff, dmax_from_Z, h_log_grid}` — for eligibility checks if `dmin_dmax_inventory.csv` lacks a needed field.
- `lrg_eegfc.workflow.fc.load_fc_matrix` — only if eligibility recomputation is needed (default: read pre-computed inventory).
- `lrg_eegfc.config.const.{BRAIN_BANDS, BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT}` — band ordering and figure labels.
- `lrg_eegfc.config.paths.{IMCOH_LRG_CACHE, FIGURES_ROOT}` — never hardcode paths.
- **No private re-implementations.** Per `coding-rules.md`. If a new helper is needed (e.g. `longest_contiguous_ridge`), promote on second use.

**Acceptance gate.** No matrix-cell evaluation runs until BOTH this scope (`2026-04-29_cohort-coverage-matrix.md`) AND `2026-04-29_decision-rules.md` are committed. The two are a paired-PR. Audit_18 v1 lands as a separate PR after both scopes.

**v1 vs v2 split.**
- v1 (rebuild plan §10 step 3): builds matrix from current-coverage inputs; L3, L4_null, L5_hrel_null, L6_null, L7_null cells emit `pi = uncontrolled`. v1 verdicts are valid for figure rendering but always carry the `T_pending_controls` annotation.
- v2 (rebuild plan §10 step 16): re-emits matrix after §7.1–§7.6 controls land. `T_pending_controls` becomes empty for the headline-triangulated bands; the figure replaces hatched cells with proper verdict borders.

---

## Open questions

1. **Should the per-cell CSV carry the per-leaf evidence path?** v1 schema records `producer_artefact` as the cohort-level CSV (e.g. `per_patient_hierarchy_cnp/leaf_assignment.csv`). For per-leaf rungs (L4, L6, L7) a reader who wants to drill down still needs the per-leaf row indices. Trade-off: pointer column adds noise to most rows. Provisional answer: defer; the producer artefact path is sufficient since per-leaf evidence is filterable by `(patient, band)` directly in the per-leaf CSV.
2. **Should `band_verdict_n10_imcoh_abs.csv` carry the per-`s` rows AND the band-level aggregate, or split into two CSVs?** v1 plan: single CSV with `scale_value = NaN` rows for band-level, non-NaN rows for per-`s`. Split would simplify pivoting but doubles the file count. Provisional: single CSV, document the convention in the column descriptions.
3. **MSPC stream-vs-load.** MSPC is ~407k rows; stream-with-chunksize is the safe default. If the matrix becomes a routine recompute, materialising MSPC once into a per-(p, b, k) `n_T-leaves` Parquet would be a one-line speedup. Provisional: stream in v1; revisit if recompute time exceeds 30s.
4. **What `T_pending_controls` value qualifies for `headline-triangulated`?** Decision-rules §2.4 defines `headline-triangulated` strictly as `V(L1) = positive AND V(L5_k) = positive AND (V(L3) = positive λ ≥ 0.5 OR V(L7) = positive)`. v1 cannot satisfy the L3/L7 condition (both uncontrolled). Triangulation CSV will emit `T = uncontrolled-exploratory` with `T_pending_controls = "L3,L7"` for any band where v1 has L1 + L5_k positive. v2 resolves to `headline-triangulated` once L7 control lands. This is correct behaviour; document explicitly so the user is not confused by v1 having no `headline-triangulated` bands.
5. **Schema versioning.** v1 → v2 schema is forward-compatible (append-only). If a v3 ever needs to *rename* a column, do it as a new dated scope report (`2026-MM-DD_cohort-coverage-matrix-v3.md`) and supersede this one. Provisional: v1 / v2 share schema; no version field in the CSV header (use the run metadata JSON for versioning).
6. **Rasterised heatmap dpi.** Default `dpi=200` per `CLAUDE.md`. For a 6 × 9 cell grid this is overkill; could drop to `dpi=100`. Provisional: stick to `dpi=200` for consistency with project convention.
7. **Should the figure also carry per-patient strips below the matrix?** Showing each patient's `pi` per (band, rung) as a tile-strip below the cohort heatmap would let a reader spot whether a positive verdict is broad cohort or 8/10 with two near-misses. Trade-off: figure complexity. Provisional: defer to a supplementary figure in `2026-04-29_task-trace-rebuild-verdict.md` (rebuild plan §10 step 19).
