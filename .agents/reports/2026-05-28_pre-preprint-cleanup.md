---
name: 2026-05-28_pre-preprint-cleanup
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-05-28
updated: 2026-05-29
pointers:
  - /home/giulio/.claude/plans/goofy-crafting-treehouse.md
---

# Pre-preprint cleanup — rolling TLDR

One block per phase, updated as each phase completes. Plan is at
`/home/giulio/.claude/plans/goofy-crafting-treehouse.md` (user-approved
2026-05-28).

---

## Phase 1 — Preprint-unblock sign fixes — **DONE** (with 1 follow-up question)

**What changed:**

1. `scripts/01_compute/audit/audit_48_epileptic_n10_compute.py:29` — docstring
   flipped to `T_KC_epi = d_KC(rest_pre, task) − d_KC(task, rest_post)`,
   "Positive = trace". Compute code (line 414) was already correct.
2. `scripts/01_compute/audit/audit_48b_epileptic_n10_figures.py:32` — same
   docstring flip.
3. `scripts/02_preprint/preprint_08_beta_grassmann_figure.py:192–193` —
   LaTeX y-label flipped to
   `T_G(k) = d_chord^{rsPre,taskT} − d_chord^{taskT,rsPost}` (the cached
   `obs_T_G` column from `audit_70` already carries the new convention, so
   only the label was wrong).
4. `.agents/reports/2026-04-29_result-1-raw-fc-phase-trace.md` — archived
   to `.agents/reports/archive/2026-04/` (`git mv`), `status: superseded`,
   sign-convention banner added explaining how to read it in the current
   convention. Numerical content unchanged.
5. `.agents/START_HERE.md` — full rewrite. Now reflects:
   `IMCOH_ABS × COHORT_N10` (n=10), locked T_d sign convention, current
   Section-5 headline (`2026-05-05_result-2-lrg-beta-trace.md`), preprint
   routing via `.agents/preprint/WRITING_GUIDE.md`, and a flagged item:
   `move_to_root` doesn't exist in any module — only `move_to_rootf`.
6. `.agents/era-map.md` — `COHORT_N9` block annotated as superseded by
   `COHORT_N10` on 2026-04-25; new `COHORT_N10` block added with locked
   sign convention + landmark doc pointers.

**What was NOT changed (audit was wrong):**

- `.agents/preprint/directives/writing_directive_2026-05-20_methods_audit_application.md`
  — the audit agent flagged "line ~380" (file is 118 lines). The actual
  sign discussion at lines 35–36 is **correct**: line 35 quotes the
  WRONG manuscript text under "Current text:", and line 36 prescribes
  the right fix. The directive is a TODO instruction for a future LaTeX
  writing pass, not a claim itself. No edit needed.
- `scripts/01_compute/audit/audit_48_rf_k_clade_persistence.py` and
  `audit_48b_rf_k_soft_jaccard.py` — both already use the NEW sign
  convention (`T_RF > 0 = trace direction`). The audit's "minor reword
  for lexical order" suggestion is cosmetic, not a correctness issue.
- `scripts/02_preprint/preprint_02_beta_rho_split_raw_D.py:35` — the
  audit called `move_to_rootf` a typo, but `move_to_rootf` is the actual
  function in `lrgsglib.config.funcs:28`. `move_to_root` (no 'f') doesn't
  exist anywhere. Real fix is either to add an alias or to correct
  CLAUDE.md's canonical header — see follow-up below.

**Smoke test verdict (live tree, excluding archives):**

- 4 categories of residual `T_d < 0` references found:
  - **3 LIVE reports** (status: current) carrying OLD-sign text:
    `.agents/reports/2026-05-05_result-2-lrg-beta-trace.md`,
    `.agents/reports/2026-05-06_section-5-critical-review.md`,
    `.agents/reports/2026-05-08_trace-minus-epi-resection.md`.
    User decision needed — see Q1 below.
  - **1 preprint hub doc** already marked `status: mostly_superseded`
    (`.agents/preprint/directives/archive/2026-05/methods_section_review_2026-05-19.md`) —
    OLD-sign references inside a superseded doc are acceptable; no
    action.
  - **1 writing bundle**
    (`.agents/writing-bundles/raw-fc/README.md`) — substrate-specific
    T_d; the audit accepted this with an era-tag clarification
    pending. Tracked for Phase 2 or later.
  - **Educational/historical references** in
    `.agents/preprint/HANDOFF_INDEX.md`, `METHODS_AUDIT_ISSUES.md`,
    superseded `2026-04-28_raw-fc-phase-distance-verdict.md` —
    intentional; left alone.
- 1 stale active plan (`2026-04-29_eigenvector-direct-pivot-plan.md`)
  uses OLD-sign notation. Already scheduled for archive in Phase 6
  (Batch 6 decision).

**Follow-up answers applied (2026-05-28 afternoon):**

- **Q1 (user choice: archive all 3 with banner)** — applied:
  - `2026-05-05_result-2-lrg-beta-trace.md`,
    `2026-05-06_section-5-critical-review.md`, and
    `2026-05-08_trace-minus-epi-resection.md` each got
    `status: superseded`, `sign_convention: pre-2026-05-26`, and a
    top-of-file banner explaining how to read the report in current
    convention.
  - All three `git mv`'d to `.agents/reports/archive/2026-05/`.
  - Bulk `sed` updated 16 citers across `.agents/` and
    `data/reports/notes_verification/` to point at the new archive
    paths. Diary entries (`2026-05-07.md`, `2026-05-08.md`) kept OLD
    paths intentionally — they are dated snapshots.
- **Q2 (user choice: update CLAUDE.md + AGENTS.md to say
  `move_to_rootf`)** — deferred to Phase 2, which already touches
  CLAUDE.md and AGENTS.md.

Phase 1 is now fully complete. No outstanding sign-convention questions.

---

## Phase 2 — Lock meta-rules into CLAUDE.md + AGENTS.md + memory — **DONE**

**What changed:**

1. `CLAUDE.md` Library-first section — added new bullet: library
   module and function names under `src/lrg_eegfc/` NEVER reference
   manuscript-local scope (`section3`, `figureN`, `preprint`, etc.).
   Use general graph / network / statistics / I/O concepts
   (`network_layouts`, `network_drawing`, `tree_metrics`,
   `surrogate_helpers`, `patient_io`). Locked 2026-05-28.
2. `CLAUDE.md` Never list — added two new entries:
   - Never use `imshow_colorbar_caxdivider` for a colorbar shared
     across multiple columns in the same row (keep explicit
     `make_axes_locatable` / `fig.add_axes` pattern for that case).
   - Never name a library module / function after a manuscript-local
     token (cross-link to the Library-first bullet).
3. `CLAUDE.md` Plotting rule 2 (colorbars) — appended a "**Scope:**"
   paragraph describing the multi-column-shared-cbar exception.
4. `CLAUDE.md` Technical invariants + Standard imports — flipped the
   canonical notebook header from `move_to_root(pathname=...)` to
   `move_to_rootf(pathname=...)` (the actual function that exists in
   `lrgsglib.config.funcs:28`). Both occurrences updated via
   `replace_all`. Closes Q2 from Phase 1.
5. `AGENTS.md` ≡ `CLAUDE.md` mirrored via `cp` (project rule). `diff`
   confirmed zero-byte difference.
6. `.agents/guides/04_rules/never-always-list.md` — added the two
   never entries (library-local-names + colorbar-multi-column-scope),
   each citing the new feedback memories.
7. `.agents/guides/04_rules/coding-rules.md` — added hard-rule 7
   (library names are general, not local-scope).
8. `.agents/guides/05_plotting/README.md` rule 2 — appended the
   `imshow_colorbar_caxdivider` scope clarification.
9. Two new memory files saved:
   - `feedback_library_names_general.md` — library naming meta-rule
     with examples of forbidden tokens and good names.
   - `feedback_imshow_colorbar_caxdivider_scope.md` — colorbar helper
     scope rule, mechanism (`make_axes_locatable` attaches to one
     axis), and concrete row-shared-cbar pattern.
10. `MEMORY.md` index updated with the two new entries.
11. `START_HERE.md` — removed the temporary "move_to_root is incorrect"
    note now that CLAUDE.md says `move_to_rootf`. Notebook header block
    is canonical and self-consistent.
12. Also caught while doing Phase 2: 3 stale OLD-sign references inside
    `CLAUDE.md` itself (lines 47, 50, 96 in original). All flipped to
    NEW sign (`T_d > 0` = trace at the cross-phase taxonomy block;
    `T_d < 0` = reset at the reset definition). Locked T_d rule at
    lines 125–128 was already correct and unchanged.

**Smoke test verdict:**

- `grep -n "move_to_root\b\|move_to_rootf" CLAUDE.md AGENTS.md` returns
  only `move_to_rootf` (2 occurrences each, identical).
- `diff CLAUDE.md AGENTS.md` is empty.
- The 4 docs touched are internally consistent: never-list, coding-rules,
  plotting README, and CLAUDE.md all carry both new rules with the same
  wording and date stamp (2026-05-28).
- No outstanding questions before Phase 3.

---

## Phase 3 — Default flips (MSC → imcoh_abs, PNG → PDF) — **DONE**

**What changed:**

1. `src/lrg_eegfc/cli/_common.py`:
   - `fc_method_option(default=...)` flipped to `"imcoh_abs"` (line 72).
   - `output_options` `--format` default flipped to `"pdf"` (line 119).
     Choices reordered to `["pdf", "png", "svg", "eps"]` so `--format`
     help shows the canonical format first. Help text updated:
     "PDF is canonical (vector, no PNG sibling). PNG is opt-in."
2. `src/lrg_eegfc/cli/compute.py` — 3 sites (`lrg`, `time-windows`,
   `reorganization`) flipped via `replace_all` (`@fc_method_option(default="msc")` → `"imcoh_abs"`).
3. `src/lrg_eegfc/cli/plot.py` — 4 sites flipped via `replace_all`.
4. `src/lrg_eegfc/utils/scripting.py` — `load_all_lrg` and `load_all_fc`
   defaults flipped + docstrings updated to name all 5 FC methods and
   call out `imcoh_abs` as the current-era default.
5. `src/lrg_eegfc/visuals/spatial.py` — 5 function signatures + 5
   docstring lines flipped (all `fc_method = "msc"` → `"imcoh_abs"`, all
   `Default is "msc"` doc lines updated to name `imcoh_abs` as the
   current-era default and list the other 4 options).
6. `src/lrg_eegfc/cli/bundle.py` — `rglob("*.png")` replaced with the
   `--png` opt-in pattern: default collects `*.pdf`; pass `--png` to
   collect `*.png`. Docstring clarifies the canonical format is PDF.

**Smoke test verdict:**

- `grep -rnE 'fc_method *[:=] *"msc"' src/` returns **empty**.
- `grep -rn 'default="png"' src/lrg_eegfc/cli/` returns **empty**.
- `lrg-eegfc plot lrg --help` → `--fc-method` default `imcoh_abs`,
  `--format` default `pdf`.
- `lrg-eegfc compute lrg --help` → `--fc-method` default `imcoh_abs`.
- `lrg-eegfc bundle overleaf --help` → `--png` opt-in flag visible;
  docstring explains PDF is canonical.

**Notes:**

- `cli/_common.py` doc TODO at lines 79–81 (`migrate to two-level
  --coh-method / --coh-transform API`) was left in place — that's a
  larger interface redesign tracked at
  `/home/giulio/.claude/plans/enumerated-napping-codd.md`, out of
  scope for this cleanup.
- The MSC subpackage in `src/lrg_eegfc/utils/fc/msc/` was kept intact
  per Batch 3 decision E.2 (user wanted MSC accessible). Defaults
  now make `imcoh_abs` the path of least resistance, but `--fc-method
  msc` still works for diagnostics.

---

## Phase 4-A — Library promotions + colorbar canonicalization + PatientRecording unification — **DONE**

(The long-file splits, section F of the plan, are tracked as Phase 4-B
below and will run after Phase 5.)

**G.1 — Library helper promotions:**

- `lrg_eegfc.utils.metrics.hypothesis.surrogate_p_value(obs, surr_array,
  tail='upper'|'lower', *, add_one=True)` — canonical Phipson–Smyth
  upper/lower-tail empirical p-value. Replaces ad-hoc
  `(1 + (surr >= obs).sum()) / (R + 1)` in audit_62/65/66/67/70.
- `lrg_eegfc.utils.metrics.hypothesis.loo_sensitivity(cohort_values,
  test_fn, *, labels=None) -> dict` — leave-one-out wrapper enforcing
  the 2026-05-19 rule (see `feedback_no_single_patient_p_driven`).
  Returns `{full, loo, worst, worst_patient}`. Descriptive only.
- `lrg_eegfc.utils.metrics.tree.cophenet_matrix(Z, condensed=False)` —
  canonical `squareform(cophenet(Z))` pattern (4+ in-script copies in
  audit_48/54/71). Handles the scipy 2-tuple vs 1-tensor return.
- `lrg_eegfc.utils.io.patient.PatientMasks` dataclass +
  `build_epi_masks(patient, root_path=SEEG_DATAPATH) -> PatientMasks` —
  one canonical call replaces 4 in-script `_build_masks` copies. Uses
  string set membership (avoids the `np.isin(int, str)` silent-False
  bug fixed in audit_71 commit e97d249; see
  `audit_epi_mask_pattern` memory).

**Audit-mistake found and corrected:** the plan named six helpers; only
five exist as named. `compute_vi` is already in the library at
`lrg_eegfc.utils.metrics.vi:82` (imported by `partition_vi_on_subset`),
and `tree_distance_bootstrap` does not exist in any source script (the
bootstrap patterns under `scripts/01_compute/audit/` are bespoke
`bootstrap_null` / `bootstrap_wilcoxon_z_p` per-script — not worth
premature abstraction). Memoised in the audit follow-up section below.

**Tests:** 22 new tests across two new test files
(`tests/test_metrics_tree.py`, `tests/test_io_patient.py`) + extended
`tests/test_hypothesis.py`. Total suite 37 passing,
0 failing. New test count by helper: surrogate_p_value (6),
loo_sensitivity (4), cophenet_matrix (5), build_epi_masks (4),
existing helpers regression (18).

**D — Colorbar canonical:**

- `src/lrg_eegfc/visuals/correlation.py:25–33` — removed the
  in-module re-implementation of `imshow_colorbar_caxdivider`; replaced
  with `from lrgsglib.plotlib.colorbars import imshow_colorbar_caxdivider`.
  All 5 call sites in the same file now use the canonical helper.
- `src/lrg_eegfc/visuals/__init__.py` — added the canonical re-export
  so `from lrg_eegfc.visuals import imshow_colorbar_caxdivider` works,
  with a comment block citing the scope rule and the locked
  `feedback_imshow_colorbar_caxdivider_scope` memory.
- `visuals/msc.py`, `visuals/spatial.py`, `visuals/lrg.py` —
  **inspected, not changed**. Those modules use
  `make_axes_locatable` in shared-row-colorbar patterns (multi-column
  panels) which is the explicit exception the new rule documents.
  Replacing them would break the layout. Left alone.

**G.3 — PatientRecording unification:**

- Audit-finding correction: the audit claimed three identical class
  definitions in `patient.py`, `patient_robust.py`,
  `patient_bootstrap.py`. Reality: only two
  (`patient.py:63` + `patient_robust.py:29`) — `patient_bootstrap.py`
  defines `MigrationAction` / `MigrationPlan` (data-layout migration
  tooling), not a `PatientRecording` variant.
- The two `PatientRecording` classes had **identical fields**
  (`timeseries`, `parameters`, `channel_metadata`); strategy-kwarg
  merge was unnecessary. Replaced the `patient_robust.py` definition
  with `from .patient import PatientRecording` and a comment
  explaining the 2026-05-28 unification.
- Smoke check: `from .patient import PatientRecording as A;
  from .patient_robust import PatientRecording as B; assert A is B`
  passes. All existing tests still pass.

**Smoke test verdict:**

- `pytest tests/ -q` → 37 passing, 1 skipped (MSC surrogate scaling,
  env-gated by `LRG_EEGFC_SURROGATE_TEST`).
- `python -c "from lrg_eegfc.utils.metrics.hypothesis import
  surrogate_p_value, loo_sensitivity; from
  lrg_eegfc.utils.metrics.tree import cophenet_matrix; from
  lrg_eegfc.utils.io.patient import PatientMasks, build_epi_masks"` →
  no errors.
- `python -c "from lrg_eegfc.visuals import
  imshow_colorbar_caxdivider"` → no errors, points to
  `lrgsglib.plotlib.colorbars`.

## Phase 4-B — Long-file splits (F) — **DONE (2026-05-29)**

7 long files split into single-purpose library modules per the F
directive of the plan. Started after Phases 5-A + 6 + 7 landed so
F operates on the cleaned tree (fewer false-positive churn paths).
Ordering: safest first (single-script caller surface), risky last
(big library callers, click registration).

### Split 1/7 — `figures_for_notes/_shared.py` (1246 → 458 LOC) — DONE

- **Promoted to** `src/lrg_eegfc/visuals/network_layouts.py` (246 LOC)
  — `compute_network_layout` (the 7-strategy dispatch:
  `mds_ultrametric`, `mds_log_ultrametric`, `mds_lrg_continuous`,
  `kk_ultrametric`, `spring_lrg_distance`, `sbm_nested`, `spring`),
  `compute_percolation_threshold`, internal helpers
  `_seeded_initial_pos` + `_compute_sbm_data`.
- **Promoted to** `src/lrg_eegfc/visuals/network_drawing.py` (539 LOC)
  — `EDGE_GAMMA`, `scale_edge_weights`, `draw_network_edges` (the
  bucketed-LineCollection edge renderer with same-probe vs
  cross-probe highlighting), `_probe_color_map`, `render_sbm_panel`,
  `render_lrg_panel` (graph-tool driven panels).
- **Module names** follow the locked library-naming meta-rule
  (`network_layouts` / `network_drawing` not `section3_layout`):
  general graph/network concepts, importable by Section-7 figures
  tomorrow without renaming.
- **`_shared.py` keeps script-local helpers** (`apply_pub_style`,
  `set_memory_limit`, `load_channel_labels`, `load_fc_for_patient`,
  `save_fig`, network metrics, constants `METHOD_COLORS`,
  `SECTION2_ROOT`, etc.) and re-exports the promoted symbols at the
  top of the module so all 23 caller scripts continue to work
  unchanged (no shim, no deprecation warning — the library is the
  single canonical home, `_shared.py` re-imports it like any other
  script module).
- **Public surface** updated: `lrg_eegfc.visuals.__init__` adds the 7
  new exports.
- **Verification**: identity check confirms
  `_shared.compute_network_layout is
  lrg_eegfc.visuals.network_layouts.compute_network_layout` (one
  canonical object). `pytest tests/` 37 passing / 1 skipped. All 23
  caller files (19 figures_for_notes + 3 section3 + audit
  `q_bundle_figures.py`) compile clean.

### Split 2/7 — `audit_round3_section5_redo.py` (2236 LOC, 10 redos) — DONE

- **Family scripts** generated under `scripts/01_compute/audit/`:
  - `audit_round3_section5_psi.py` (redo1 + redo2, 306 LOC)
  - `audit_round3_section5_kc.py` (redo3 + redo6..redo10 + dendrogram helpers, 1577 LOC)
  - `audit_round3_section5_ctm.py` (redo4, 127 LOC)
  - `audit_round3_section5_cross_probe.py` (redo5, 100 LOC)
- **Shared header** at `_audit_round3_shared.py` (240 LOC) holds imports,
  output-dir constants, `asterisks` helper, and the 6 dendrogram helpers
  (`_reorder_linkage_to_match`, `_enumerate_clades`, `_max_jaccard_match`,
  `_clade_leaf_xrange`, `_draw_clade`, `_best_jaccard_partner`).
- **Orchestrator** `audit_round3_section5_redo.py` (40 LOC) imports each
  family's `main()` and runs them in sequence — preserves the original
  bundle behaviour for callers that want everything in one go.
- **Verification**: all 6 files compile clean; runtime import-check
  confirms `_audit_round3_shared` exposes all 21 expected symbols and
  all 10 redo functions exist across the family modules.

### Split 3/7 — `q_bundle_figures.py` (1317 LOC, 9 figures) — DONE

- **Per-figure scripts** generated under `scripts/01_compute/audit/queries/`:
  `q_fig1_trace_scatter.py`, `q_fig2_phase_geometry.py`,
  `q_fig3_dS_dP_convergence.py`, `q_fig4_drift.py`,
  `q_fig5_distance_class.py`, `q_fig6_chord_arc.py`,
  `q_figS1_z_inflation.py`, `q_figS2_swarm.py`,
  `q_figS3_stoplight.py` (each wraps its original procedural section
  in `def main(): ... ; if __name__ == "__main__": main()`).
- **Shared header** at `_q_bundle_shared.py` (141 LOC) holds CLI
  substrate dispatch (`resolve_substrate`), data loader
  (`load_bundle_data`), phase-pair palette (`PAIR_COLORS` /
  `PAIR_LABELS` / `PAIR_ORDER`), `lookup_pair` helper, and an
  inline `load_channel_labels`.
- **`sys.path.insert(scripts/archive/2026-02_imcoh-dev-notes)` removed.**
  `compute_network_layout`, `draw_network_edges`, `_probe_color_map`
  now imported from `lrg_eegfc.visuals.network_layouts` /
  `network_drawing` (canonical post split 1/7). `load_channel_labels`
  inlined into the shared header (25 LOC).
- **Orchestrator** `q_bundle_figures.py` (61 LOC) imports each fig's
  `main()` and runs them. Original CLI substrate arg (`imcoh_abs` /
  `imcoh_sq`) still supported.
- **Verification**: all 11 files compile clean; runtime import-check
  confirms `_q_bundle_shared` exposes all 20 expected symbols and
  each fig module's `main()` is reachable.

### Split 4/7 — `visuals/lrg.py` (1106 → 711 LOC) — DONE

- **New module** `src/lrg_eegfc/visuals/lrg_panels.py` (423 LOC) holds
  the composite `plot_lrg_full_panel`. The 5-panel composite is
  self-contained — it only uses three internal helpers
  (`_load_channel_labels`, `_find_psi_optimal_partition`,
  `compute_partition_stability_index`), which it imports from `.lrg`.
- **`lrg.py` keeps** the individual-panel helpers
  (`plot_lrg_entropy_curves`, `plot_lrg_dendrogram`,
  `plot_lrg_dendrogram_shaft_colored`, `plot_ultrametric_heatmap`) +
  the three private utilities + `compute_partition_stability_index`.
- **Public surface** updated in `visuals/__init__.py`:
  `plot_lrg_full_panel` now comes from `.lrg_panels`. Identity check
  confirms `lrg_eegfc.visuals.plot_lrg_full_panel is
  lrg_eegfc.visuals.lrg_panels.plot_lrg_full_panel`.
- **Caller rewire**: 1 site updated.
  `src/lrg_eegfc/cli/plot.py:693` flipped from
  `from lrg_eegfc.visuals.lrg import plot_lrg_full_panel` to
  `from lrg_eegfc.visuals.lrg_panels import plot_lrg_full_panel`.
  Other 5 callers (`notebook.py`, `cli/plot.py:264`,
  `scripts/02_visualize/visualize_lrg.py:35`, etc.) all use
  `from lrg_eegfc.visuals import plot_lrg_full_panel` and continue
  to work via the re-export.
- **Verification**: `pytest tests/` 37 passing / 1 skipped; identity
  check + symbol absence in `lrg.py` both confirmed.

### Split 5/7 — `cli/plot.py` (1025 → package) — DONE

- **Module → package conversion**: `cli/plot.py` deleted; new
  `cli/plot/` package with one submodule per command family:
  - `__init__.py` (33 LOC) — defines the `@click.group() plot`
    + side-effect imports of the 4 submodules to register commands.
  - `_fc.py` (505 LOC, 7 cmds) — `corr`, `msc`, `msc-grid`,
    `msc-all-patients`, `msc-validation`, `cleaning`, `comparison`.
  - `_lrg.py` (223 LOC, 3 cmds) — `lrg`, `lrg-phase-grid`,
    `lrg-video`.
  - `_reorg.py` (235 LOC, 5 cmds) — `reorganization`, `reorg-metrics`,
    `reorg-summary`, `metric-correlation`, `metastable`.
  - `_misc.py` (178 LOC, 2 cmds) — `cross`, `time-windows`.
- **Click registration**: Each submodule does `from . import plot`
  to fetch the group; `@plot.command(...)` decorators register
  on import (`__init__` triggers `from . import _fc, _lrg, _reorg,
  _misc`). Trial-swap test: round-trip enumerates exactly the
  17 expected commands, no extras, no missing.
- **Defaults preserved**: `lrg-eegfc plot lrg --help` still shows
  `--fc-method` default `imcoh_abs` and `--format` default `pdf`.
- **No external caller edits needed**: package vs single-module
  swap is transparent to all importers (`from lrg_eegfc.cli.plot
  import plot` still works).

### Split 6/7 — `visuals/spatial.py` (1529 → 678 LOC) — DONE

- **New module** `visuals/spatial_coords.py` (446 LOC) — coordinate
  loaders + MNI transforms: `_normalize_label`,
  `load_spatial_metadata`, `prepare_spatial_coordinates`,
  `_estimate_mni_transform`.
- **New module** `visuals/spatial_nilearn.py` (477 LOC) — glass-brain
  wrappers: `view_brain_connectome`, `plot_brain_connectome`,
  `plot_brain_connectome_at_n`. Imports the two coord helpers from
  `.spatial_coords`.
- **`spatial.py` keeps** the core Plotly + Matplotlib 3D network
  plots: `build_edge_traces`, `_get_cluster_colors`,
  `plot_spatial_network_3d`, `plot_spatial_network_3d_mpl`,
  `plot_spatial_clusters_comparison`. Re-exports the 5 moved-out
  symbols at the top of the module for backwards compat.
- **No external caller edits needed**: 6 preprint-cited scripts
  doing `from lrg_eegfc.visuals.spatial import load_spatial_metadata,
  prepare_spatial_coordinates` continue to work via the re-export.
  Public surface in `visuals/__init__.py` unchanged.

### Split 7/7 — `visuals/network_templates.py` (2453 → 1406 LOC) — DONE

- **New module** `visuals/network_chord.py` (757 LOC) — hierarchy
  chord plots (graph-tool curvy edges through LRG): the
  `HIERARCHY_*` + `DENDROGRAM_*` constants,
  `build_chord_depth2_layout`, `render_hierarchy_chord`,
  `_chord_px_transform`, `draw_circular_dendrogram_overlay`,
  `draw_radial_leaf_labels`, `_chord_default_edge_arrays`,
  `plot_chord_with_dendrogram`, `plot_chord_with_highlight`.
- **`visuals/network_layouts.py` extended** (246 → 609 LOC) — gains
  the 11 pure `(A, ...) → (N, 2)` layouts from `network_templates.py`
  (`layout_spring`, `layout_kk`, `layout_spectral`,
  `layout_laplacian_pca`, `layout_circular_by_shaft`,
  `layout_community_grouped`, `layout_backbone_guided`,
  `layout_lrg_kk`, `layout_sfdp`, `layout_arf`, `layout_lrg_sfdp`)
  + `LAYOUT_REGISTRY` + `DEFAULT_GALLERY` + `is_lrg_layout` +
  `compute_layout`. Now sits alongside the
  `compute_network_layout` (Section-2 dispatch) from split 1/7.
- **`network_templates.py` keeps** the non-layout / non-chord core:
  constants (`SAME_PROBE_RGB`, `CROSS_PROBE_RGB`, `DEFAULT_GAMMA`,
  …), loaders (`load_probe_labels`, `matrix_to_gt`, `nx_to_gt`),
  preprocessing (`rank_transform`, `power_transform`,
  `disparity_backbone`, `top_k_backbone`), edge renderer
  (`draw_gamma_edges`), node decoration (`shaft_colors`,
  `community_colors`, `draw_nodes`), top-level templates
  (`plot_fc_network`, `plot_fc_network_row`, `plot_fc_network_grid`,
  `plot_fc_network_lrg`, `plot_layout_gallery`), matrix+network
  (`spring_auto_k`, `plot_fc_matrix_and_network`,
  `plot_fc_matrix_and_network_rows`).
- **Re-exports** the moved-out symbols at the **top** of the module
  (not bottom — `DEFAULT_GALLERY` is used as a default arg value in
  `plot_layout_gallery`, evaluated at def-time, not call-time).
  All 23 external callers (`.agents/guides/05_plotting/`,
  `audit_47..52`, `diag_*`, `preprint_07_test2`) compile clean
  without per-caller edits.
- **`visuals/__init__.py` public surface unchanged** — the same
  symbols re-export, just from a now-properly-split internal
  structure.

### Phase 4-B totals

Original combined LOC: 11206 (7 files).
Post-split: still ~12K LOC across 28 files (some growth from
explicit imports + docstrings in new modules), but every file is
now single-purpose. The four library-facing splits (4, 5, 6, 7)
preserved the public surface via top-of-module or `__init__`
re-exports — zero external caller edits needed.

---

## Phase 5-A — Library style sweep + Pat_03 marker removal — **DONE**

**Library `src/lrg_eegfc/visuals/`:**

- 6 `fig.suptitle(...)` calls removed on publication-grade plot
  functions, replaced with one-line comments citing the rule:
  - `visuals/compare.py:180` (`plot_fc_comparison`)
  - `visuals/lrg.py:901` (`plot_lrg_full_panel`)
  - `visuals/reorganization.py:409` (`plot_phase_reorganization`)
  - `visuals/reorganization.py:639` (`plot_phase_distance_matrix`)
  - `visuals/msc.py:348` (`plot_msc_comparison_dense_vs_validated`)
  - `visuals/msc.py:575` (`plot_msc_summary`)
- Final library state: `grep -rnE "^[^#]*\.suptitle\(" visuals/`
  returns 0 actual calls (only doc-comments and the rule-explainer
  comments). `set_rasterized` already absent from code (3 mentions
  in doc strings only). No default PNG-only saves in library code.

**Pat_03 outlier-marker removal (C.sub):**

- `scripts/07_figures/gen_radar_final.py` — `OUTLIER = "Pat_03"`
  constant deleted; the polygon loop now uses uniform `ls="-"`,
  `marker="o"`, no asterisk in legend; the dashed-line footnote
  `"*Pat_03 (dashed): 1024 Hz acquisition — outlier at $\gamma_h$"`
  removed; markdown companion + docstring updated to flag the n=5
  cohort as an era artifact and to cite `feedback_pat03_no_dropout`.
- `scripts/07_figures/gen_h2b_figures.py` — same edits applied:
  `OUTLIER` constant removed, uniform line style in the radar loop,
  footnote text removed, narrative updated.

**Library import smoke:** `python -c "import lrg_eegfc.visuals"`
returns no errors.

## Phase 5-B — Preprint-cited script style sweep — **DONE (2026-05-29)** (see end of report)

Inventory of remaining script-level violations (post-archive):

- `fig.suptitle` in 15 scripts across
  `scripts/01_compute/figures_embedded/` (`fig_evidence_chain`,
  `fig_dendrogram_grid`, `fig_band_persistence_trajectories`,
  `fig_band_unanimity_curves`, `fig_trace_vs_ergodic`) +
  `scripts/07_figures/` (10 files including `consolidate_report_figures`,
  `hypotheses_final_figures`, `fig_S5_*`, etc. — most going to
  archive in Phase 6).
- `set_rasterized(True)` in 7 scripts: 5 in `figures_embedded/`
  (`fig_trace_step1`, `fig_task_trace_band_k_n10`,
  `fig_continuous_trace_per_cell`, `fig_continuous_trace_cohort`,
  `fig_mrl_landscape`) + 2 in `07_figures/`
  (`gen_corr_vs_msc_figures`, `gen_fc_figures_fast`).

Phase 5-B will run **after Phase 6** so it operates on the
already-pruned set (Phase 6 archives most of the offending 07_figures
files: `gen_corr_vs_msc_*`, `gen_mslcd_*`, `plot_crema_*`, plus the
v1/v2 duplicates). The residual sweep is then ~7 files, not 22.

---

## Phase 6 — Archive sweep + v1/v2 collapse + dead stubs — **DONE**

**MSC-era exploratory folders → `scripts/archive/2026-05_metric-exploration-era/`:**

- `scripts/03_analysis/` (18 files, ~10K LOC) — `git mv`'d.
- `scripts/04_reorganization/` (13 files) — `git mv`'d.
- `scripts/06_metric_sweep/` (10 files) — `git mv`'d.
- `scripts/wp0/` (7 files + `_common.py`) — `git mv`'d.
- `scripts/wp1/` (7 files) — `git mv`'d.
- Stale `scripts/05_multiscale/` files (4) reading the legacy
  `partition_multiscale/results.csv`: `multiscale_community_flow`,
  `multiscale_affinity_analysis`, `multiscale_h2_profile`,
  `explore_raw_vi` → `05_multiscale_legacy/` subfolder.

**Dead library stubs:** `git rm` of `src/lrg_eegfc/visuals/lrg_backup.py`
and `src/lrg_eegfc/visuals/lrg_revised.py` (3 LOC each, zero callers).

**Notebook archive:** `ipynb/02_fc_msc/` → `ipynb/90_archive/2026-04_fc_msc_era/02_fc_msc/`.

**Stale active plans:** `2026-04-25_surface-multiscale-trace.md` and
`2026-04-29_eigenvector-direct-pivot-plan.md` → `.agents/plans/archive/2026-04/`.

**V1/V2 duplicate collapse:**

- `audit_round2_section5_v2.py` + `_v2_fixes.py` → archived.
  `audit_round3_section5_redo.py` remains live and will be split in
  Phase 4-B.
- `gen_cross_phase_comparison_v2.py` → archived; v1 stays live.
  Parametrization (`--row3-metric` flag) deferred to Phase 4-B.
- `gen_brain_connectome_multiscale.py` (MSC variant) → archived;
  `gen_brain_connectome_imcoh.py` stays live. Parametrization
  deferred to Phase 4-B.
- `scripts/08_epileptic/epileptic_lrg_report_v2.py` renamed to
  `epileptic_lrg_report.py` (no v1 existed; `_v2` suffix dropped).

**Shell scripts (folded in from Phase 7):** broken `run_full_analysis.sh`
and `run_step.sh` (called non-existent `scripts/py/*.py`) → archived to
`scripts/archive/2026-04_legacy-shell/`.

**Archive READMEs added:**

- `scripts/archive/2026-05_metric-exploration-era/README.md` —
  describes the 6 sub-collections + the rescue list of helpers that
  were promoted to the library BEFORE the moves (Phase 4-A).
- `scripts/archive/2026-04_pre-rebuild/README.md` — first-time
  documentation of the 3 sub-eras (msc-era, mrl-era, dead-branches).

**Smoke test verdict:**

- `grep -rn "from scripts.03_analysis\|from scripts.04_reorganization\|
  from scripts.06_metric_sweep\|from scripts.wp0\|from scripts.wp1\|
  from lrg_eegfc.visuals.lrg_backup\|from lrg_eegfc.visuals.lrg_revised\|
  epileptic_lrg_report_v2"` over live tree → **zero hits**. Only one
  comment in `gen_fig_J_concordance_comparison.py:129` references
  `scripts/wp0/_common.py` historically (comment, not import) — left
  alone since it's documentation.
- `pytest tests/ -q` → 37 passing, 1 skipped. No regressions.
- Final live `scripts/` directory: 8 folders (00_prepare, 01_compute,
  02_preprint, 02_visualize, 05_multiscale, 07_figures, 08_epileptic,
  09_surrogate) + archive/. Down from 14 folders, ~38% reduction in
  live-script surface.

---

## Phase 7 — Docs, notebooks, root hygiene — **DONE**

**Root hygiene (J + LICENSE):**

- `pyproject.toml`:
  - `license = {file = "LICENSE"}` → `license = {text = "GPL-3.0-or-later"}`.
  - classifier flipped from `License :: OSI Approved :: MIT License` →
    `License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)`
    so PyPI metadata matches the LICENSE file.
  - `package-data` updated from `["*.txt", "*.md"]` to
    `["*.txt", "*.md", "visuals/styles/*.mplstyle"]` so the project
    mplstyle ships with installed wheels.
- `.gitignore` — appended `.pytest_cache/`, `.mypy_cache/`,
  `.ipynb_checkpoints/`. Verified no tracked `__pycache__` /
  `.pytest_cache` remain.
- `.gitmodules` — branch stays `main` (already correct); the
  README-vs-`.gitmodules` mismatch was resolved by the README
  rewrite (see below) which no longer instructs `git checkout
  lrg_eegfc`.

**`README.md` — full rewrite.** New structure:

1. One-paragraph what-and-why.
2. **Current era (locked)** block: FC carrier `imcoh_abs`, n=10
   cohort, T_d sign convention (positive = trace), preprint hub
   pointer, agent landing page pointer.
3. Quick start.
4. Installation (5 steps; `lrgsglib` submodule on `main` — the
   spurious `git checkout lrg_eegfc` instruction removed).
5. **Canonical CLI flow** with `imcoh_abs` examples — replaces the
   old MSC-only examples. Lists the 7 subcommand groups and a typical
   pipeline.
6. Dataset layout (4-bucket `data/`).
7. Python API overview with `imcoh_abs`-era imports + the 4 helpers
   promoted in Phase 4-A (`surrogate_p_value`, `loo_sensitivity`,
   `cophenet_matrix`, `build_epi_masks`).
8. Plotting utilities (canonical colorbar helper + 8-rule summary).
9. Agent + developer guide pointers (`START_HERE.md`, `era-map.md`,
   `CLAUDE.md`, `never-always-list.md`).
10. License notice — GPL-3.0-or-later.

**Guide catch-up (I.guides):**

- `cli-reference.md` — resynced to 40 subcommands (was 33–36). Added
  `compute.diagnostics`, `compute.threshold-analysis`, `plot.cross`,
  `data.normalize`. `updated:` bumped.
- `caching-guide.md` — added 6 missing cache subdir rows
  (`matched_strength_surrogate_lrg/`,
  `matched_strength_surrogate_epi_excluded_lrg/`,
  `matched_strength_per_pair/`, `imcoh_halves_fc/`,
  `imcoh_lrg_halves/`, `spectral_pctl/`) with a note that those are
  written by audit scripts, not via the CLI, so they don't have
  first-class `config.paths` constants yet.
- `function-map.md` — added a 2026-05-28 frontmatter note
  (`coverage: representative-not-exhaustive`) plus a banner
  explaining the ~30% coverage and pointing readers to the source
  code as canonical.

**Notebook header sweep (L.2):**

- Audit said 13 notebooks had the legacy `from lrgsglib import
  move_to_rootf` header with `pathname='lrg_eegfc'`. Reality: only
  **3 active notebooks** (`00_intake/01_sanity_imports`,
  `01_preprocessing/time_windows_analysis`,
  `04_reorganization/distance_measures_comparison`) carried the bad
  pathname; the other 33 active notebooks already use
  `move_to_rootf(pathname="lrgeegfc")`. Plus 1 loose `.py` file in
  `ipynb/dev/`. All 4 updated via a Python notebook-JSON script that
  preserves cell structure; the `from lrgsglib import move_to_rootf`
  legacy import was replaced with `from lrg_eegfc.notebook import *`
  and the pathname argument was flipped to `"lrgeegfc"`.
- The wider Phase 5-B `use_lrg_style()` activation across 18 figure
  notebooks + `fig.suptitle` strip in 6 cells is deferred with the
  other script-level style work (Phase 5-B runs after Phase 4-B and
  before final verification).

**`ipynb/INDEX.md` update:**

- 5 `02_fc_msc/` rows (NB-002, NB-003, NB-004, NB-050, NB-051)
  flipped to `archived` with updated paths under
  `ipynb/90_archive/2026-04_fc_msc_era/`.
- Added 9 missing entries: 2 `05_figures/` (NB-063 `10_spatial_network_3d`,
  NB-064 `11_spatial_cluster_comparison`) + 7 `06_presentation_figures/`
  (NB-065..NB-071: `01_msc_matrices_per_band` through
  `07_fc_section_figures`).

**Smoke test verdict:**

- `pytest tests/` → 37 passing, 1 skipped.
- `python -c "import lrg_eegfc, lrg_eegfc.visuals, lrg_eegfc.cli"` →
  no errors.
- `diff -q CLAUDE.md AGENTS.md` → no output (still identical).

## Phase 5-B — Script-level style sweep (preprint-cited subset) — **DONE (2026-05-29)**

After Phases 4-B's splits and 6's archive sweep landed, the residual
style-rule violations in preprint-cited scripts and `figures_embedded/`
were swept in a focused pass.

- **5 `fig.suptitle()` calls removed** in `figures_embedded/`:
  `fig_band_persistence_trajectories.py`, `fig_band_unanimity_curves.py`,
  `fig_dendrogram_grid.py`, `fig_evidence_chain.py`,
  `fig_trace_vs_ergodic.py`. Each replaced with a one-line comment
  pointing to the publication-figure no-suptitle rule.
- **2 `set_rasterized(True)` calls removed** in `figures_embedded/`:
  `fig_mrl_landscape.py:159`, `fig_trace_step1.py:252`. PDFs are now
  fully vector per the locked rule.
- **`use_lrg_style()` activation inserted in 7 preprint-cited scripts**
  that lacked it (`preprint_02..06` + `gen_h2b_figures.py` +
  `gen_radar_final.py`). Insertion site: just after the top import
  block, before any module-level code.
- **Verification**: all 14 touched files compile clean; `pytest tests/`
  37 passing / 1 skipped throughout.

## Final 10-point verification — **DONE (2026-05-29)**

| # | Check | Result |
|---|---|---|
| 1 | Library + CLI imports clean | [OK] |
| 2 | `lrg-eegfc plot` surface (17 subcommands round-trip) | [OK] |
| 3 | CLI defaults (`--fc-method=imcoh_abs`, `--format=pdf`) | [OK] |
| 4 | `pytest tests/` → 37 passed / 1 skipped | [OK] |
| 5 | No live OLD-sign `T_d` formula in compute sites | [OK] |
| 6 | `CLAUDE.md == AGENTS.md` | [OK] |
| 7 | Public `visuals` surface intact (12+ representative symbols) | [OK] |
| 8 | Phase 4-B split files all exist (11 new modules) | [OK] |
| 9 | No archive `sys.path.insert` in LIVE scripts | [OK] |
| 10 | No `fig.suptitle` / `set_rasterized` in preprint-scope | [OK] |

Repository state at completion:

- `src/lrg_eegfc/visuals/` — 19 library modules (was ~10).
- `src/lrg_eegfc/cli/plot/` — 5 subcommand files (was a single 1025-LOC `plot.py`).
- 37 tests passing throughout the entire cleanup.
- Branch `audit/cohort-n10-diagnostic`, HEAD `6d7374a`.

## Out-of-scope post-preprint (no work here)

- Cache cleanup (~2.4 GB stale).
- Broader test coverage push for `workflow.fc/lrg`, `utils.metrics.tree*`,
  `utils.surrogate.matched_strength`, `visuals/*`, `cli/*`.
- `_legacy.py` CLI removal.
- `data/audit/` reorganization into the canonical 4-bucket layout.
