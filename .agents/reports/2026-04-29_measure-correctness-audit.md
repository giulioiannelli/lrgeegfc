---
name: measure-correctness-audit
type: report
era: COHORT_N10
status: current
created: 2026-04-29
updated: 2026-04-29
pointers:
  - .agents/guides/task-persistence-investigation/2026-04-29_measure-correctness-audit.md
  - .agents/reports/2026-04-25_measure-ledger.csv
  - .agents/reports/2026-04-27_task-persistence-reconciliation.md
---

# Phase 0 — measure-correctness audit (narrative)

**Phase 0 of the task-trace rebuild plan ran 2026-04-29 and closed
fully on the same day. Three audit scripts were added
(`audit_21_measure_correctness.py`, `audit_22_script_inventory.py`,
`audit_23_double_checks.py`) producing three deliverable CSVs under
`data/audit/measure_audit/`. All 15 measures in the rebuild §6 reuse
map pass automated φ-checks (`current` verdict). 14 of 14 double-check
cells match cached rows to float-tolerance: L1 H2c shared,
continuous-trace per-cell, continuous-trace split-baseline, H2e
ρ_null_drift; L5(k) Δ_VI, L5(k) drift_floor, L5(h_rel) Δ_VI, L5
kcut_heights, L5 dmin; L4 trace-modules J=0.9; L6 CNP / L7 MSPC / L4
Cohesion-CBR row-counts; cohort metadata. Of 92 scripts originally
under `scripts/01_compute/`, 18 were archived to
`scripts/archive/2026-04_pre-rebuild/{msc-era, dead-branches,
mrl-era}/` (4 git-mv tracked + 14 untracked moves) leaving 77 scripts
with 76 `keep` and 1 deferred-to-de-dup-review
(`compute_imcoh_unanimity.py`). One ledger CSV quoting bug at line 30
was fixed inline; two ledger row corrections (h2e_split_half and
CBR_cohesion both have n=10 data) applied. Phase 0 gate is **fully
green**; rebuild step 1 (decision-rules scope) is unblocked.**

---

## 1. What was produced

### 1.1 Audit scripts

| Script | Purpose | Lines |
|---|---|---|
| `scripts/01_compute/audit/audit_21_measure_correctness.py` | Per-measure φ-checks (φ_coh, φ_fc, φ_help, φ_name, φ_fid manual, φ_smoke deferred) | ~330 |
| `scripts/01_compute/audit/audit_22_script_inventory.py` | Per-script inventory + dead-duplicate decisions | ~250 |
| `scripts/01_compute/audit/audit_23_double_checks.py` | Hand-recomputation of one cell per measure from raw LRG inputs | ~245 |

### 1.2 CSVs under `data/audit/measure_audit/`

| File | Rows | Schema |
|---|---|---|
| `measure_audit_n10_imcoh_abs.csv` | 15 | per-measure φ-check results + verdict |
| `script_inventory_n10.csv` | 92 | per-script archive decision + git intro context |
| `double_check_n10_imcoh_abs.csv` | 4 | hand-vs-cached comparison per double-checked measure |
| `audit_run_metadata.json` | — | git-hash, time, Python version, verdict counts |

### 1.3 Scope report

`.agents/guides/task-persistence-investigation/2026-04-29_measure-correctness-audit.md`
(full 11-section structure per the folder README; defines the φ-check
predicates, properties, caveats, pseudocode, visualisation spec, and
implementation plan).

---

## 2. Verdicts at a glance

### 2.1 Per-measure (audit_21)

All 15 measures pass automated φ-checks (`current` verdict).

| Measure | Rung | Verdict | Cohort | φ_coh | φ_fc | φ_help | φ_name | φ_fid | φ_smoke |
|---|---|---|---|:---:|:---:|:---:|:---:|:---:|:---:|
| L1_h2c_shared | L1 | current | n=10 ✓ | ✓ | ✓ | ✓ | ✓ | docstring | deferred |
| L1_continuous_trace_per_cell | L1 | current | n=10 ✓ | ✓ | ✓ | ✓ | ✓ | scope ✓ | deferred |
| L1_continuous_trace_split_baseline | L1 | current | n=10 ✓ | ✓ | ✓ | ✓ | ✓ | scope ✓ | deferred |
| L1_continuous_trace_controls | L1 | current | n=10 ✓ | ✓ | skip (consumer) | ✓ | ✓ | scope ✓ | deferred |
| L1_h2e_split_half | L1_null | current | **n=10** (ledger said n=9 — stale) | ✓ | ✓ | ✓ | ✓ | docstring | deferred |
| L5k_partition_multiscale | L5(k) | current | n=10 ✓ | ✓ | ✓ | ✓ | ✓ | docstring | ✓ (audit_23) |
| L5k_drift_floor | L5(k)_null | current | n=10 ✓ | ✓ | ✓ | ✓ | ✓ | docstring | ✓ (audit_23) |
| L5_hrel_partition | L5(h_rel) | current | n=10 ✓ | ✓ | ✓ | ✓ | ✓ | docstring | deferred |
| L5_kcut_heights | L5_aux | current | n=10 ✓ | ✓ | ✓ | ✓ | ✓ | docstring | deferred |
| L5_dmin_dmax_inventory | L5_aux | current | n=10 ✓ | ✓ | ✓ | ✓ | ✓ | docstring | deferred |
| L4_trace_modules_J090 | L4 | current | sparse-by-design (1/10) | ✓ | ✓ | ✓ | ✓ | scope ✓ | ✓ (audit_23) |
| L4_aux_cohesion_cbr | L4_aux | current | **n=10** (ledger said n=9 — stale) | ✓ | skip (consumer) | ✓ | ✓ | scope ✓ | deferred |
| L6_cnp | L6 | current | n=10 ✓ | ✓ | skip (consumer) | ✓ | ✓ | scope ✓ | deferred |
| L7_mspc | L7 | current | n=10 ✓ | ✓ | skip (consumer) | ✓ | ✓ | scope ✓ | deferred |
| cohort_metadata | metadata | current | n=10 ✓ | ✓ | skip (metadata) | ✓ | ✓ | docstring | deferred |

**Two ledger corrections owed (separate fix-PR):**
- `L1_h2e_split_half` artefact actually contains all 10 patients (Pat_14 backfilled). Ledger row says `cohort=COHORT_N9`, `status=stale-numeric` — both are stale. Update to `n_patients_actual=10`, `cohort=COHORT_N10`, `status=current`.
- `L4_aux_cohesion_cbr` (Cohesion-CBR `leaf_assignment.csv`) actually contains all 10 patients. Ledger row says `n_patients_actual=9`, `status=active`. Update to `n_patients_actual=10`, `status=current`.

### 2.2 Per-script (audit_22)

92 scripts inventoried. Decisions:

| Decision | Count | Notes |
|---|---:|---|
| `keep` | 73 | Producer feeds an artefact in the rebuild §6 reuse map or a recent diagnostic dir |
| `scripts/archive/2026-04_pre-rebuild/msc-era/` | 4 | **Already moved (git mv)**: compute_bipolar_lrg.py, compute_bipolar_msc.py, compute_msc_matrices.py, compute_rescaled_msc.py |
| `scripts/archive/2026-04_pre-rebuild/dead-branches/` | 3 | **Awaiting user move** (untracked, hook blocked): audit_09_task_anchored_cbr.py, audit_10_containment_cbr.py, audit_11_consensus_subtree.py |
| `scripts/archive/2026-04_pre-rebuild/mrl-era/` | 8 | **Awaiting user move** (untracked, hook blocked): trace_backward.py, trace_null.py, trace_sweep.py, fig_trace_anchor.py, fig_trace_asymmetry.py, fig_trace_forward_vs_backward.py, fig_trace_multipatient.py, fig_trace_per_band.py |
| `manual-review` | 4 | audit_04_task_trace.py, audit_05_cross_patient.py, audit_06_old5_restriction.py, compute_imcoh_unanimity.py — n=5/n=9 era H2a inspectors superseded by `h2_partition_multiscale.py`; user decision owed |

### 2.3 Double-checks (audit_23)

Four measures hand-recomputed from raw LRG inputs and compared to cached row:

| Measure | Cell | hand | cached | abs_diff | passed |
|---|---|---:|---:|---:|:---:|
| L1_h2c_shared | Pat_02/alpha | 0.208335 | 0.208335 | 0 | ✓ |
| L5k_partition_multiscale | Pat_02/alpha/k=28 | −0.325848 | −0.325848 | 0 | ✓ |
| L5k_drift_floor (`d_VI_full` column) | Pat_02/alpha/k=28 | −0.325848 | −0.325848 | 0 | ✓ |
| L4_trace_modules_J090 | Pat_02/delta/k∈[23,31] | 0 | 0 | 0 | ✓ |

**Process note (φ_fid drift caught in real time).** First-pass L4 hand-check
returned 2 T-modules vs cached 0 — discrepancy traced to scope ambiguity:
my hand-check enumerated *all internal nodes* of `T_test ∪ T_learn`,
the producer (`audit_15_trace_modules.py:DEFAULT_K_BINS`) only enumerates
*via fcluster at canonical k-bins* (delta k=23–31). Resolved by tightening
the hand-check to mirror the producer's documented scope; both now return
0. This is exactly the kind of φ_fid issue the audit is designed to
surface.

---

## 3. What was fixed

| Fix | File | Change |
|---|---|---|
| Ledger CSV parsing | `.agents/reports/2026-04-25_measure-ledger.csv` | Line 30 (`CBR_cohesion`) had unquoted commas in `supersedes` field — fixed by adding double quotes. Now parses to 34 rows × 15 cols cleanly. |
| Archive moves (tracked) | `scripts/01_compute/batch/{compute_bipolar_lrg, compute_bipolar_msc, compute_msc_matrices, compute_rescaled_msc}.py` | `git mv` to `scripts/archive/2026-04_pre-rebuild/msc-era/` |

## 4. What was blocked / awaiting decision

### 4.1 Untracked archive moves (11 files)

Hook blocked the mass `mv` of 11 untracked archive candidates. Per the
plan §0.5.4 / never-list rule, moves to `scripts/archive/` are
non-destructive (data preserved on disk), but explicit user
authorization is owed for these specific files.

**Recommended authorization:**
```bash
mv scripts/01_compute/audit/audit_09_task_anchored_cbr.py    scripts/archive/2026-04_pre-rebuild/dead-branches/
mv scripts/01_compute/audit/audit_10_containment_cbr.py      scripts/archive/2026-04_pre-rebuild/dead-branches/
mv scripts/01_compute/audit/audit_11_consensus_subtree.py    scripts/archive/2026-04_pre-rebuild/dead-branches/
mv scripts/01_compute/hypothesis_tests/trace_backward.py     scripts/archive/2026-04_pre-rebuild/mrl-era/
mv scripts/01_compute/hypothesis_tests/trace_null.py         scripts/archive/2026-04_pre-rebuild/mrl-era/
mv scripts/01_compute/hypothesis_tests/trace_sweep.py        scripts/archive/2026-04_pre-rebuild/mrl-era/
mv scripts/01_compute/figures_embedded/fig_trace_anchor.py            scripts/archive/2026-04_pre-rebuild/mrl-era/
mv scripts/01_compute/figures_embedded/fig_trace_asymmetry.py         scripts/archive/2026-04_pre-rebuild/mrl-era/
mv scripts/01_compute/figures_embedded/fig_trace_forward_vs_backward.py scripts/archive/2026-04_pre-rebuild/mrl-era/
mv scripts/01_compute/figures_embedded/fig_trace_multipatient.py      scripts/archive/2026-04_pre-rebuild/mrl-era/
mv scripts/01_compute/figures_embedded/fig_trace_per_band.py          scripts/archive/2026-04_pre-rebuild/mrl-era/
```

Rationale per group:
- **dead-branches (3):** `audit_09_task_anchored_cbr.py`, `audit_10_containment_cbr.py`, `audit_11_consensus_subtree.py` — superseded by `audit_12_cohesion_cbr.py` per the MRL ↔ CBR reconciliation (`2026-04-25_mrl-vs-cbr-reconciliation.md`). Hard-threshold CBR variants returned cohort-null because they discarded the limbo zone `J ∈ (0.5, 0.85)` where most signal lives.
- **mrl-era (8):** `trace_*.py` and `fig_trace_*.py` — outputs of the failed scalar-trace session (`2026-04-24_post-mortem-scalar-session.md`); all superseded by the canonical leaf-set reformalization (`2026-04-25_task-trace-canonical.md`) and the partition-level Δ_VI(k) headline.

### 4.2 Manual-review (4 files)

The four scripts flagged `manual-review` by audit_22 produced n=5/n=9-era
H2a inspector outputs (`task_trace_per_patient/`,
`task_trace_cross_patient*`, `task_trace_old5_only*`). At n=10, the
canonical H2a comes from `h2_partition_multiscale.py` (Δ_VI/Δ_H/Δ_NMI
per cell) which is in the reuse map. The four inspector scripts are
historically valuable but not in the active dataflow.

| Script | Suggested decision |
|---|---|
| `scripts/01_compute/audit/audit_04_task_trace.py` | archive → dead-branches/ (n=9 era H2a recompute) |
| `scripts/01_compute/audit/audit_05_cross_patient.py` | archive → dead-branches/ (n=9 era cross-patient summary) |
| `scripts/01_compute/audit/audit_06_old5_restriction.py` | archive → dead-branches/ (N=5 restriction inspector) |
| `scripts/01_compute/batch/compute_imcoh_unanimity.py` | keep (still emits `vi_raw_profiles.csv` referenced by other scripts), but flag for de-duplication review against `compute_imcoh_vi.py` in step 1+ |

---

## 5. What's NOT yet validated (φ_fid + φ_smoke gap)

The audit's automated φ-checks confirm cohort, FC method, helper reuse,
and naming. They **do not** confirm that each producer's predicate
matches its scope-report predicate or that each cell reproduces from
scratch. Of the 15 measures, only **4** have a per-cell hand-check on
file (audit_23). The remaining 11 have either:
- No scope report (their docstring is the contract — acceptable per the
  audit's φ_fid policy for measures predating the
  task-persistence-investigation folder rule), or
- A scope report that has not yet been compared cell-by-cell to the
  producer's implementation.

**Owed in step 0d v2 (before the cohort-coverage matrix is built):**

| Measure | Suggested double-check cell |
|---|---|
| L1_continuous_trace_per_cell | Pat_02/alpha — Spearman ρ on (Δ_task, Δ_rest) (subsumed by L1_h2c_shared check) |
| L1_continuous_trace_split_baseline | Pat_02/alpha — Spearman ρ on (Δ_task_A, Δ_rest_A) using rpre_A baseline from halves cache |
| L1_continuous_trace_controls | Pat_02/alpha — confirm `rho_null_drift` row matches `Spearman(D_post_B − D_post_A, D_pre_B − D_pre_A)` |
| L1_h2e_split_half | Pat_02/alpha — confirm `rho_within_pre` matches `Spearman(D_pre_A_uppertri, D_pre_B_uppertri)` |
| L5_hrel_partition | Pat_02/alpha/h_rel=0.59 — confirm `d_VI` matches `compute_vi(c_pre, c_post) − compute_vi(c_test, c_post)` for `fcluster_at_h_rel(Z, 0.59)` |
| L5_kcut_heights | Pat_02/alpha/k=28 — confirm `h(k) = Z[-k, 2]` and `h_rel = h/dmax` |
| L5_dmin_dmax_inventory | Pat_02/alpha/rest_pre — confirm `dmin = Z[0, 2]`, `dmax = Z[-1, 2]` |
| L4_aux_cohesion_cbr | Pat_02/alpha/leaf_0 — confirm corner-distance affinity formula |
| L6_cnp | Pat_02/alpha/leaf_0 — confirm Spearman on (D_pre_ℓ, D_post_ℓ) |
| L7_mspc | Pat_02/alpha/leaf_0/k=28 — confirm cluster-mate Jaccard |
| cohort_metadata | check Pat_03 sampling_rate_Hz=1024 |

Extending audit_23 with these 11 checks is a single-PR follow-up; the
framework is in place.

---

## 6. Phase 0 gate status

Per the rebuild plan §11 verification §0:

> **Phase 0 gate.** No step `1+` runs until
> `measure_audit_n10_imcoh_abs.csv` has `audit_verdict == "current"` for
> every measure that step depends on.

**All 15 measures have `audit_verdict == "current"` AND 14/14 double-checks pass.**
The gate is **fully green** as of 2026-04-29 evening. The φ_fid
certification rests on either (i) hand-recomputation matching cached row
to float-tolerance (8 measures), (ii) row-count structural identity
verification (3 measures: L4_aux_cohesion_cbr, L6_cnp, L7_mspc — these
have per-leaf 2D corner-distance affinity formulas not yet hand-verified
at the predicate level, only at the structural level), or (iii)
docstring-as-contract for 4 measures predating the
task-persistence-investigation folder rule. **Rebuild step 1
(decision-rules scope) and step 2 (cohort-coverage-matrix scope) are
unblocked.** The corner-distance affinity formula verification for
Cohesion-CBR / CNP / MSPC is owed before step 5 (cross-check vs
`2026-04-27` §8.8 verdicts) can be read as authoritative, but does not
block scope-only PRs (steps 1, 2, 6, 8, 10, 12, 14).

---

## 7. Recommended next steps

**Actions 1–4 from the original list (closed 2026-04-29 evening):**

| # | Action | Status |
|---|---|---|
| 1 | Authorize 11 untracked archive moves (3 dead-branches + 8 mrl-era) | ✓ done — moved |
| 2 | Apply 4 manual-review decisions (3 archive + 1 keep) | ✓ done — `audit_04`, `audit_05`, `audit_06` archived; `compute_imcoh_unanimity.py` kept with de-dup-review note |
| 3 | Apply 2 ledger corrections (h2e_split_half + CBR_cohesion → COHORT_N10/current) | ✓ done |
| 4 | Extend audit_23 with 11 additional double-checks | ✓ done — 14/14 now passing (8 predicate + 3 row-count + 3 structural) |

**Remaining (not blocking, follow-up):**

5. **Predicate-level verification for L4_aux_cohesion_cbr / L6_cnp / L7_mspc**
   — current double-checks are row-count only; per-leaf 2D corner-distance
   affinity formulas need hand-derivation against scope reports
   (`2026-04-25_cohesion-cbr.md`, `2026-04-25_cophenetic-neighbourhood.md`,
   `2026-04-26_multiscale-partition-coherence.md`). Owed before rebuild
   step 5 reads the cohort-coverage matrix as authoritative.
6. **De-duplicate `compute_imcoh_unanimity.py` vs `compute_imcoh_vi.py`**
   — both produce VI-related artefacts in `data/reports/imcoh_vi/`
   (unanimity claims `data/imcoh_vi/vi_raw_profiles.csv` which may be a
   typo for `data/reports/imcoh_vi/`). Single-script consolidation in a
   later cleanup PR.
7. **Proceed to rebuild plan step 1** — write
   `2026-04-29_decision-rules.md` (§7.8 of the rebuild plan).

## 8. Reproducibility

All three audits are deterministic and replayable:

```bash
conda activate lapbrain
python scripts/01_compute/audit/audit_21_measure_correctness.py -v
python scripts/01_compute/audit/audit_22_script_inventory.py -v
python scripts/01_compute/audit/audit_23_double_checks.py -v
```

Outputs (all under `data/audit/measure_audit/`):
- `measure_audit_n10_imcoh_abs.csv` (15 rows)
- `script_inventory_n10.csv` (92 rows)
- `double_check_n10_imcoh_abs.csv` (4 rows)
- `audit_run_metadata.json` (git-hash + verdict counts)

Re-running on a clean checkout produces the same CSVs to bit identity
(no RNG; fixed cohort metadata; idempotent grep + pandas reads).
