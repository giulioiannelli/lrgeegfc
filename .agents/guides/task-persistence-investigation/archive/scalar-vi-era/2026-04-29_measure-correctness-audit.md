---
name: measure-correctness-audit
type: scope
era: COHORT_N10
status: draft
created: 2026-04-29
updated: 2026-04-29
pointers:
  - .agents/guides/task-persistence-investigation/README.md
  - .agents/reports/2026-04-25_measure-ledger.csv
  - .agents/reports/2026-04-27_task-persistence-reconciliation.md
  - .agents/guides/04_rules/coding-rules.md
  - .agents/guides/04_rules/never-always-list.md
---

# Phase 0 — measure-correctness audit (BLOCKING)

**Phase 0 of the task-trace rebuild: before any new measure or matrix-cell
evaluation, audit every existing measure end-to-end (script↔scope-report
fidelity, cohort/era correctness, library-helper reuse, smoke
reproducibility) and inventory every script under `scripts/01_compute/`
for dead-duplicate archival. Two CSVs (`measure_audit_n10_imcoh_abs.csv`,
`script_inventory_n10.csv`) plus a narrative report gate every downstream
rebuild step. A cohort-coverage matrix built on a silently-broken or
duplicate measure is worse than no matrix.**

---

## Notation

- `M` — set of measures listed in the §6 reuse map of the rebuild plan
  (`/home/giulio/.claude/plans/i-think-there-are-binary-puppy.md`):
  H2c controls, H2e split-half, partition-multiscale, ΔVI fixed-h_rel,
  ΔVI split-baseline, dmin/dmax, trace-modules J=0.9, Cohesion-CBR, CNP,
  MSPC.
- `S` — set of scripts under
  `scripts/01_compute/{hypothesis_tests, diagnostics, batch, audit, figures_embedded}/`.
- For each `m ∈ M`: `producer(m) ∈ S` is the script that wrote the
  cached artefact; `scope(m)` is the scope report under
  `.agents/guides/task-persistence-investigation/` that defines the
  predicate `m` should compute; `artefact(m)` is the CSV / NPZ on disk.
- For each `s ∈ S`: `predicate(s)` is the docstring / first-40-line
  summary of what `s` computes; `output(s)` is the artefact path.
- `helpers_lib = {compute_vi, conditional_entropy, wilcoxon_z,
  rank_biserial, boot_ci_mean, bh_fdr, cluster_stats,
  jaccard_leafsets, fcluster_at_h_rel, tree_internal_nodes,
  dmax_from_Z, h_log_grid, simpson_neff, cluster_size_stats,
  partition_vi_on_subset, kc_distance, matching_cluster_distance,
  weighted_rf_distance, load_fc_matrix}` — the canonical library
  primitives every measure must reuse (see plan §6).
- Cohort: `P_n10 = {Pat_02, Pat_03, Pat_05, Pat_06, Pat_07, Pat_08,
  Pat_10, Pat_13, Pat_14, Pat_15}`.
- Bands: `B = {δ, θ, α, β, γ_l, γ_h}`. Phases:
  `Φ = {rest_pre, task_learn, task_test, rest_post}`.
- FC method: `imcoh_abs = ⟨|ImCoh|⟩_f` per `CLAUDE.md`.

## Definitions

For each measure `m ∈ M` the audit produces a tuple
`A(m) = (φ_fid, φ_coh, φ_fc, φ_help, φ_smoke, φ_name, V) ∈ {0, 1}^6 × Verdict`
with components defined below. A measure passes iff
`φ_fid = φ_coh = φ_fc = φ_help = φ_smoke = φ_name = 1`.

**φ_fid (script ↔ scope-report fidelity).** Equals 1 iff the predicate
implemented by `producer(m)` matches the predicate stated in `scope(m)`.
Concretely: aggregator (mean / median), correlation type
(Spearman / Pearson / Kendall), axis of reduction (over patients vs over
k vs over leaves), masking (band-mask, probe-mask, NaN handling), and
sign convention (e.g. `Δ_VI = VI(pre,post) − VI(test,post)` not the
reverse). Any mismatch ⇒ `φ_fid = 0`, regardless of magnitude.

**φ_coh (cohort actual vs claimed).** Equals 1 iff the unique patient set
in `artefact(m)` equals `P_n10` AND the era flag in the file (`fc_method`,
`cohort` columns or filename token) is `imcoh_abs` × `n10`. n=9 artefacts
with explicit `era=COHORT_N9` flag pass this check (they are not lying),
but consumers of `m` at n=10 must be re-flagged as `stale-numeric`.

**φ_fc (FC method actual vs claimed).** Equals 1 iff the cache path the
producer reads from is `data/cache/imcoh_lrg/` (or `imcoh_lrg_halves/`)
AND any band-averaging step honours Jensen's inequality
(per-frequency-bin `np.abs` / squaring before band-mean — see
`workflow/fc.py:205-258`). A producer that averages over signed `imcoh`
then takes `np.abs` violates the taxonomy
(`memory/imcoh_taxonomy.md`).

**φ_help (helpers used vs duplicated).** Equals 1 iff every primitive in
`helpers_lib` that the script needs is *imported* (not re-implemented).
A local `def _wilcoxon_z(...)` or `def _vi(...)` ⇒ `φ_help = 0` and
the duplicate is enumerated in the audit `notes` field with line range.

**φ_smoke (smoke reproducibility).** Equals 1 iff re-running
`producer(m)` on `Pat_02 / band ∈ {alpha} / phase ∈ {rest_pre, task_test, rest_post}`
(or the minimal cell the measure needs) produces a row that matches
the corresponding cached CSV row to float-tolerance
`abs(Δ) ≤ max(1e-9, 1e-6 · |x|)`. Stochastic helpers
(`boot_ci_mean`, sign-flip permutation) re-run with their original seed
or pass with looser tolerance noted.

**φ_name (naming-sensitivity check).** Equals 1 iff cache filenames
match the canonical patterns in `CLAUDE.md`:
- ImCoh freq-resolved: `{band}_{phase}_imcoh_freqresolved_nperseg-{N}.npy`
  with `N == nperseg_for_fs(fs)` (4096 at 2048 Hz, 2048 at Pat_03 1024 Hz).
- ImCoh LRG: `{band}_{phase}_lrg_imcoh-{abs|sq}.npz`.
- ImCoh LRG halves: `{band}_{phase}_{A|B}_lrg_imcoh-abs.npz`.
- Correlation: `filter-{none|abs}`, `zero_diag-{true|false}`.
- MSC: `sparsify-{none|soft}`, `nperseg-{value}`, `n_surrogates-{value}`.

**Verdict V ∈ {current, stale-numeric, broken, duplicate, archive}.**
- `current` ⇔ `φ_fid = φ_coh = φ_fc = φ_help = φ_smoke = φ_name = 1`.
- `stale-numeric` ⇔ `φ_fid = φ_fc = φ_help = φ_smoke = φ_name = 1` but
  `φ_coh = 0` solely because the artefact is at n=9 — fixable by Pat_14
  backfill alone.
- `broken` ⇔ `φ_fid = 0` OR `φ_smoke = 0` (predicate drift or
  irreproducibility) — blocks downstream consumers; a fix-PR is owed
  before the consumer can run.
- `duplicate` ⇔ `producer(m)` re-implements a `helpers_lib` primitive
  (`φ_help = 0`); the duplicate is removed and replaced by the library
  import in the same fix-PR.
- `archive` ⇔ `producer(m)` produces only superseded artefacts (per
  `.agents/reports/2026-04-25_measure-ledger.csv` `superseded` rows or
  this audit's reclassification) — `git mv` to
  `scripts/archive/2026-04_pre-rebuild/<subdir>/`.

For each script `s ∈ S`, the audit produces an inventory tuple
`I(s) = (path, mtime, lines, predicate_summary, produces_artefact,
supersedes, archive_destination)`. `archive_destination` ∈ {`keep`,
`scripts/archive/2026-04_pre-rebuild/<subdir>/`}. A script is archived
iff `produces_artefact ∉ {current artefacts in ledger}` AND no other
current measure depends on it.

## Properties

- **Two-CSV deliverable, one verdict per row.** `measure_audit_n10_imcoh_abs.csv`
  has `|M|` rows; `script_inventory_n10.csv` has `|S|` rows. Both are
  recoverable from a fresh checkout in O(min) without re-running the
  underlying measures (only smoke checks re-run a single cell).
- **Idempotent.** Re-running the audit on the same checkout produces
  bit-identical CSVs (mtime is read once at the start of the audit run
  and stamped in a sidecar JSON).
- **Conservative.** A measure is `current` only if **all six** φ-checks
  pass. Single-check failures default to the most-blocking verdict
  (`broken` over `duplicate` over `stale-numeric`).
- **Order-of-operations preserved.** φ_fc explicitly verifies the
  Jensen-correct ImCoh transform per
  `memory/imcoh_taxonomy.md`: the cache holds *signed*
  `ImC(f)`, transforms apply *per-frequency-bin first*, then
  band-average. Reverse-order producers (`np.abs(np.mean(signed))`) are
  flagged `broken`, not `current`.
- **Non-destructive.** No script is deleted. `archive` ⇒ `git mv` only.
  Local-helper duplicates are removed from the script body but the
  original helper definition remains in git history.

**What this audit does NOT do.**
- It does **not** validate the *scientific predicate* of any measure
  beyond its scope-report definition. If a scope report says "Spearman ρ"
  and the script computes "Spearman ρ", `φ_fid = 1` even if Spearman
  is the wrong correlation for the scientific question.
- It does **not** introduce new statistical claims. Verdicts are
  bookkeeping, not new evidence.
- It does **not** replace per-rung *control* gaps (those are §7.1–§7.6
  of the rebuild plan).

## Caveats and failure modes

| Caveat | Mitigation |
|---|---|
| **Smoke reproducibility for stochastic measures** (bootstrap CI, sign-flip permutation) — re-run with same RNG seed produces identical output, but if the seed is hard-coded the smoke check is trivially circular | For RNG measures, smoke checks re-run with `--seed` matching the cached run; tolerance is widened to `1e-3` and the verdict footnote names the seed |
| **Naming convention drift** — a script may use a non-canonical cache subdir (e.g. `imcoh_abs_lrg/` instead of `imcoh_lrg/`) yet still be functionally correct | `φ_name = 0` is a **soft** failure: the verdict is `current` with a `notes` flag; rename in step 0e |
| **Library helpers absent** — a script that needs a missing library function (e.g. fixed-h_rel KC λ-sweep) is not "broken", it is "scope-only" | Such scripts are rare; verdict is `current` if the script doesn't exist yet (no producer to audit) |
| **Cross-script duplicate predicates** — two scripts implement the same predicate slightly differently (e.g. two flavours of ΔVI(k)) | The newer of the two keeps verdict `current`; the older is reclassified `archive` in `script_inventory_n10.csv` regardless of its individual `audit_verdict` |
| **Single-patient smoke is too narrow** — Pat_02 is the persistent dissenter; some measures may pass on Pat_02 but fail elsewhere | Smoke check uses Pat_02 first; if `φ_smoke = 0` only on Pat_02, retry on Pat_05 (cohort-typical) before declaring `broken` |
| **Halves cache reads** — measures that consume `imcoh_lrg_halves/` need to assert the cache contains all 10 patients × 6 bands × 2 phases × 2 halves (`rest_pre/post` only — task halves do not exist) | A producer that tries to read `task_test_A_lrg_imcoh-abs.npz` is `broken`; only `rest_pre` / `rest_post` halves exist |
| **MSC era artefacts present in cache** — pre-2026-04-15 MSC caches still on disk; MSC-era scripts may inadvertently be flagged `current` | Filter `M` to `imcoh_abs`-bound producers only; MSC producers default-archive to `scripts/archive/2026-04_pre-rebuild/msc-era/` |
| **Pat_03 1024 Hz handling** — `nperseg_for_fs(fs)` must dispatch correctly per patient | φ_fc explicitly checks Pat_03's nperseg=2048 cache file exists with that token in the filename |

## Pseudocode

```
INPUT  : reuse-map M from rebuild plan §6
       : script-tree S = scripts/01_compute/{hypothesis_tests, diagnostics,
                                              batch, audit, figures_embedded}/
       : ledger L = .agents/reports/2026-04-25_measure-ledger.csv
       : library helpers_lib (lrg_eegfc.utils.metrics.* + workflow.fc)
       : cohort_metadata = data/audit/cohort_metadata.csv
OUTPUT : measure_audit_n10_imcoh_abs.csv  (rows: |M|)
         script_inventory_n10.csv         (rows: |S|)

# Step A: per-measure audit (§0.5.1)
for each measure m in M:
    locate producer(m) by walking S for the script that writes artefact(m)
    locate scope(m)    by reading the ledger pointer or grepping
                          .agents/guides/task-persistence-investigation/
    phi_fid   = compare_predicates(producer(m), scope(m))
    phi_coh   = (unique patients in artefact(m)) == cohort_metadata.patient_set
                AND artefact has imcoh_abs token
    phi_fc    = producer(m) reads from imcoh_lrg/ or imcoh_lrg_halves/
                AND respects Jensen-order in workflow/fc.py:205-258
    phi_help  = no local re-implementation of any helpers_lib symbol in producer(m)
    phi_smoke = rerun_one_cell(producer(m), Pat_02, alpha, ...)
                  matches cached row within float-tolerance
    phi_name  = filename(artefact(m)) matches canonical patterns in CLAUDE.md
    verdict   = classify(phi_fid, phi_coh, phi_fc, phi_help, phi_smoke, phi_name)
    write_row(measure_audit_csv, m, producer, scope, phi_*, verdict, notes)

# Step B: per-script inventory (§0.5.2)
for each script s in S:
    path             = relpath(s)
    mtime            = stat(s).st_mtime
    lines            = wc -l s
    predicate_summary = first_docstring_or_top_40_lines(s)
    produces_artefact = grep "to_csv\|np.savez\|savefig" s
    supersedes        = grep "supersedes:" first 40 lines (or git log)
    archive_destination = decide_archive(produces_artefact, ledger, supersedes)
    write_row(script_inventory_csv, path, mtime, lines, ..., archive_destination)

# Step C: cross-check
assert every measure m in M has at most one producer s in S
assert every "current" measure has all six phi-checks = 1
assert every "archive" script has produces_artefact in ledger.superseded
       OR is older-of-duplicate-pair
```

## Visualization spec

Single PDF `data/audit/measure_audit/measure_audit_n10_imcoh_abs.pdf`
with two stacked panels:

**Panel A — measure verdict heatmap.** Rows: measures `m ∈ M` ordered
by rung (L1 → L7) then by ledger date. Columns: the six φ-checks plus
the final Verdict. Colour: green=1, red=0; verdict column uses
`{current=green, stale-numeric=yellow, broken=red, duplicate=orange,
archive=grey}`. Each cell carries a small text glyph if the
corresponding `notes` field is non-empty.

**Panel B — script-inventory bar chart.** One row per
`scripts/01_compute/<subdir>/`; bars stacked by `archive_destination`
(`keep` vs `archive_<subdir>`). Total bar length = file count in subdir.
Annotation: count of `current`-feeding vs `archive`-only scripts per
subdir.

Reading rules: vertical green columns = an entire φ-check passed by all
measures (suggests a healthy invariant); horizontal red rows = the
worst-offending measure (top priority for fix-PR). Panel B answers
"how much script-archival cleanup is owed" at a glance.

No `fig.suptitle`. PDF only, with imshow rasterised at dpi=200 per
`feedback_no_png_duplicates.md`. ≥3 measures × 3 subdirs visible per
panel; if `|M| > 30`, panel A splits across two pages.

## Connection to prior tools

This audit is **prior** to every measure listed in
`.agents/reports/2026-04-25_measure-ledger.csv`, not a replacement.
Where the ledger says `status: current`, this audit re-verifies that
status by code-execution rather than self-report.

Concrete relationships:

- **Subsumes** the `notes` column of the measure ledger as the
  authoritative record of script ↔ scope drift. After Phase 0, the
  ledger gets a new column `audit_verdict_2026_04_29` matching
  `measure_audit_n10_imcoh_abs.csv`'s Verdict column.
- **Complements** the rebuild plan §0.5: this scope formalises the
  predicates and outputs; the rebuild plan owns the PR sequence and
  blocking rules.
- **Does not replace** the per-rung control gaps (§7.1–§7.6 of the
  rebuild plan). Those are *new measures*; this audit only checks
  *existing measures*.
- **Reuses** all primitives in `helpers_lib`. The audit script itself
  is forbidden from re-implementing any of them — a meta-application
  of the library-first rule.

## Implementation plan

**Producers (audit code).**
- `scripts/01_compute/audit/audit_21_measure_correctness.py` — per-measure
  audit (§0.5.1). Outputs:
  `data/audit/measure_audit/measure_audit_n10_imcoh_abs.csv`.
- `scripts/01_compute/audit/audit_22_script_inventory.py` — per-script
  inventory (§0.5.2). Outputs:
  `data/audit/measure_audit/script_inventory_n10.csv`.
- `scripts/01_compute/audit/audit_23_double_checks.py` — per-measure
  hand-computed double-checks (§0.5.3). Outputs append to the
  measure-audit CSV (`double_check_passed` boolean column).

**Library entry points (must reuse, do not re-implement).**
- `lrg_eegfc.utils.metrics.tree.{tree_internal_nodes, dmax_from_Z, jaccard_leafsets, h_log_grid, fcluster_at_h_rel, simpson_neff, cluster_size_stats, partition_vi_on_subset}`
- `lrg_eegfc.utils.metrics.vi.{compute_vi, conditional_entropy}`
- `lrg_eegfc.utils.metrics.hypothesis.{wilcoxon_z, rank_biserial, boot_ci_mean, bh_fdr, cluster_stats}`
- `lrg_eegfc.utils.metrics.tree_distance.{kc_distance, matching_cluster_distance, weighted_rf_distance}`
- `lrg_eegfc.workflow.fc.load_fc_matrix`
- `lrg_eegfc.workflow.lrg.LRGResult` (npz schema)

**Outputs.**
- `data/audit/measure_audit/measure_audit_n10_imcoh_abs.csv`
- `data/audit/measure_audit/script_inventory_n10.csv`
- `data/audit/measure_audit/measure_audit_n10_imcoh_abs.pdf`
- `data/audit/measure_audit/audit_run_metadata.json` (audit-time, git-hash,
  Python version, env hash, RNG seed)
- `.agents/reports/2026-04-29_measure-correctness-audit.md` (narrative
  report — Phase 0 step 0f).

**CLI invocation (smoke).**
```
conda activate lapbrain
python scripts/01_compute/audit/audit_21_measure_correctness.py \
    --cohort n10 --fc-method imcoh_abs -v
python scripts/01_compute/audit/audit_22_script_inventory.py -v
python scripts/01_compute/audit/audit_23_double_checks.py \
    --cohort n10 --fc-method imcoh_abs -v
```

**Acceptance gate.** Phase 1+ of the rebuild plan
(`/home/giulio/.claude/plans/i-think-there-are-binary-puppy.md` §10
steps 1+) does NOT begin until `measure_audit_n10_imcoh_abs.csv` has
either `Verdict = current` or `Verdict = stale-numeric` for every
measure in the §6 reuse map. `broken` / `duplicate` rows trigger fix-PRs
and re-runs; `archive` rows trigger the `git mv` PR (step 0e).

## Open questions

1. **Smoke tolerance for `compute_lrg_analysis`.** The LRG pipeline's
   ultrametric depends on a τ-grid; tiny float drift accumulates.
   Default `1e-6` may be too strict for L4 / L5 measures derived from
   linkage-matrix cuts. Decision: relax to `1e-4` for measures whose
   inputs include `LRGResult.ultrametric_matrix`; document the
   relaxation in the per-measure `notes`.
2. **Stale-numeric vs broken for n=9 artefacts.** H2a′, H2b, H2e,
   H1-topo, H2a-topo are all flagged `stale-numeric` in the ledger
   (need Pat_14 backfill). Question: do they require a Pat_14
   smoke-cell check, or does the smoke check only run on cells already
   present in the cached CSV? Decision: smoke check runs on whatever
   `(patient, band, phase)` tuple is present — Pat_02 by default. n=9
   artefacts with `Verdict = stale-numeric` are tagged with a separate
   `pat14_backfill_owed` boolean flag.
3. **Cohesion-CBR n=9 active.** Per the ledger, Cohesion-CBR runs at
   n=9 not n=10 (per-patient hierarchy figures). Audit
   classification: `stale-numeric` (consumer of CBR will need n=10
   re-run before §5 matrix can use it as L4-aux input). Confirmed in
   `2026-04-25_cohesion-cbr.md` frontmatter.
4. **Predicate-comparison heuristic for φ_fid.** Comparing a script's
   numerical implementation to a scope's mathematical predicate is
   semi-automatic at best. Decision: φ_fid is human-judged in v1
   (audit script pretty-prints script body + scope predicate side by
   side; auditor flags); v2 may add AST-level heuristics (e.g. detect
   `scipy.stats.spearmanr` vs `pearsonr` mismatches automatically).
5. **Where does the audit's *own* correctness audit live?** Meta-rule:
   `audit_21_measure_correctness.py` itself must not duplicate
   library helpers. It is its own first reader. The double-check
   protocol (§0.5.3) is its self-application; if the audit script
   silently miscounts, the double-check protocol on a known-current
   measure will fail and surface the bug.
