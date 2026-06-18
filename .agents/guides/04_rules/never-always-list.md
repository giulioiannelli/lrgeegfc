---
name: never-always-list
type: guide
era: CROSS_ERA
status: current
created: 2026-04-24
updated: 2026-04-24
pointers:
  - CLAUDE.md
  - .agents/guides/04_rules/coding-rules.md
---

# Never / always list

**The enforced list of user preferences. Seeded from feedback memories.
When the user says "never X" or "always Y", it lands here on first
mention and gets a matching `feedback_<short>.md` memory saved.**

## Never

- Never use `lrg.optimal_threshold` as a diffusion time τ — they live
  in different spaces.
- Never call `fig.suptitle` on publication figures — context goes in
  the `.md` sidecar, not the figure.
- Never mock FC data in hypothesis-level tests — mocks masked a prod
  migration failure once already.
- Never run scripts outside the `lapbrain` conda env.
- Never start a new scalar hypothesis test when the signal is visible
  in existing VI(k) / partition-multiscale / H2d-θ artifacts — surface
  the existing signal first.
- Never build opaque matrix-level surrogate nulls (e.g. Haar-rotation
  coherency surrogates). Default null is **matched-strength**; a null must
  be explainable in one sentence and not swingable by an unconstrained
  modelling knob. **Validate any NEW null on synthetic ground-truth
  (a must-be-positive case + a must-be-null case) BEFORE running it on
  real data and reporting verdicts.** Phase-randomized physical surrogates
  do not help an `|ImCoh|` pipeline (multivariate PR preserves the
  cross-spectrum exactly → no-op; univariate destroys all connectivity) and
  are not a substitute for matched-strength. The "task-specific trace vs
  stable trait fingerprint" separation has no clean null without a task-free
  control session — state it as a one-sentence limitation, don't surrogate
  around it. (CRC + coordinated SB-CRC archived 2026-06-12;
  `feedback_no_opaque_matrix_nulls`,
  `scripts/archive/2026-06_opaque-matrix-nulls/POSTMORTEM.md`.)
- Never delete files that document research history — `git mv` to
  `<parent>/archive/YYYY-MM/` instead.
- Never pool metrics into a consensus scalar (explicitly forbidden by
  the user).
- Never skip frontmatter on a new `.agents/` .md file.
- Never invent metric names — cite literature or existing code.
- Never apply BH-FDR / Bonferroni / Holm (or any multiple-comparison
  correction) unless all three checks pass: (1) the corrected `p` /
  `q` is actually used as a verdict gate in the methods or control
  battery, (2) the family of tests is a coordinated unit of
  inference (not just "tested across N independent claims"), and (3)
  the per-test gate does not already address the multiple-testing
  concern (matched-strength surrogate at R=200 is itself a
  calibrated empirical p — layering BH on top is scaffolding without
  function). Cross-band BH-FDR on per-band `ρ_split^coph` (m=6) was
  retired 2026-05-20 by exactly this check; anatomy A1 hypergeometric
  across DK regions per (band, probe) is kept because the
  multi-region family IS coordinated. See
  `feedback_no_unmotivated_bh_fdr.md`.
- Never plot Δ_ARI(k) in partition-multiscale band×k publication
  figures — canonical 3 are `Δ_VI, Δ_H, Δ_NMI`; CSV may retain `d_ARI`.
- **Never use colormaps whose interior has a near-white band**
  (`Spectral`, `Spectral_r`, `RdYlBu`, `RdYlGn`, `RdBu`, `twilight`,
  `twilight_shifted`, `BrBG`, `PiYG`, `PuOr`) on white-background
  figures, and especially never use them when colours are
  interpolated or weighted-averaged downstream (kernel smoothing,
  size-weighted subtree-RGB mean up a dendrogram, etc.) — far-apart
  colormap positions collapse to the near-white mid-band and become
  invisible against the paper.  Prefer `turbo` (saturated rainbow,
  no near-white), the viridis family (`viridis`, `plasma`, `magma`,
  `inferno`), `cividis`, or a custom `LinearSegmentedColormap`.  If
  a diverging cmap with a neutral centre is unavoidable, use a
  saturated light-gray centre (e.g. `#dddddd`) instead of pure
  white.  See `feedback_no_near_white_cmaps.md`.
- **Never use ARI / NMI / Fowlkes-Mallows or any sklearn
  `*_rand_score` / partition-cut metric for cross-phase LRG
  dendrogram similarity, and never colour leaves / branches by
  `fcluster(Z, t=K, ...)` at any single K.** The canonical metric
  is the per-pair cophenet distance correlation
  `ρ^coph(A, B) = Spearman(D_coph_A[triu], D_coph_B[triu])` — same
  family as the §5.3 `ρ_split^coph` headline.  Reason: partition-cut
  metrics collapse the continuous tree onto one K and are dominated
  by the one large cluster that any single K produces on FC-derived
  hierarchies; they are not in the project vocabulary anywhere.  For
  cross-phase persistence visualisations, use a CONTINUOUS leaf
  colouring (e.g. task's `leaves_list` slot index → rainbow).  See
  `feedback_no_partition_metrics_use_rho_coph.md`.
- Never use channel-label letter prefixes (A/B/.../Q) for cohort-level
  implant analysis — they are arbitrary clinical labels with no
  cross-patient anatomical meaning. Use `(x, y, z)` coordinates and
  Desikan-Killany regions from `implant_pat_NN.csv` instead.
- Never title a figure with `fig.suptitle` (duplicate of above for
  emphasis — it keeps slipping through).
- Never put descriptive text in figures.  Figures are for visual
  content; only axis labels and single-letter panel tags (``"(a)"``,
  ``"β"``) are acceptable.  Forbidden: ``ax.set_title`` summary
  sentences, boxed annotations with ρ/p/N, italic floating "noise
  floor" / "null" callouts, multi-line headers with model+source+
  Wilcoxon p baked in, in-axes stat legends.  Push every number to
  stdout and let the manuscript caption / companion .md carry the
  context.  See [[feedback-no-text-in-figures]].
- Always save figures with ``transparent=True`` in ``fig.savefig(...)``.
  Default white background is wrong for our workflow (slide overlays,
  journal compositing, dark/light theme reuse).  No
  ``ax.set_facecolor(...)`` unless the user explicitly asks for an
  opaque axes background; if added, comment why the opt-in overrides
  the default.  See [[feedback-default-transparent-figures]].
- Never save figures as both PDF and PNG. **PDF only** is the default
  and only format. PNG is opt-in with explicit user request.
- Never call `im.set_rasterized(True)` (or any `set_rasterized`).
  PDFs are **fully vector** for every artist — including FC matrices,
  audit heatmaps, scatter plots, `pcolormesh`, dendrograms. The
  earlier "rasterise heavy artists" rule is **withdrawn**; vector
  is sharper and the file size for typical FC matrices (`N ≤ 130`)
  is small. If a future figure ever produces a genuinely
  unmanageable PDF, escalate to the user before changing the rule.
- Never add a grey provenance watermark / footer by default. The
  file name carries patient / band / phase / fc_method, and the
  watermark just duplicates that. **Watermark is opt-in only**:
  every plotting helper exposes `watermark=False` (or call
  `lrg_eegfc.visuals.layout.add_provenance_footer(fig, label)`).
  Use `watermark=True` only when the figure will be detached from
  its file name (slide deck, screenshot).
- Never label adjacency-matrix axes with the word "contact" /
  "channel". Use math indices `$i$` (x) and `$j$` (y).
- Never use entropy `S(τ)` or specific heat `C(τ)` as a task-trace
  comparison metric. They are non-informative for our fully-connected
  weight-heterogeneous outlier case (continuous spectrum, no discrete
  `C(τ)` peaks). The L2 entropy-curve rung is permanently removed from
  the geometric ladder; any new "entropy-curve trace" proposal is
  theatre. See `lrg-framework-guide.md` §6 + memory entry
  `lrg_outlier_case_fully_connected.md`.
- Don't confuse the four cross-phase phenomena. The taxonomy at
  `.agents/guides/01_project/terminology.md` defines **trace** (task
  reorganized AND change stuck — our `T_d < 0` finding), **anchor**
  (module unchanged across all phases), **reset** (task changed AND
  came back), **emergent** (module that did not exist in `RPre`).
  Use the explicit taxonomy in **cross-phase taxonomy tables**,
  **mixed-band figures**, **captions of cross-phase summary figures**,
  and any context where the reader could otherwise read "anchor" or
  "reset" when we mean "trace". (Softened 2026-05-18: in unambiguous
  task-trace sections — e.g. the β results subsection where the
  Methods has already defined "persistence into rsPost" as the
  target phenomenon — bare "persistence" is acceptable and reads
  cleanly. Save the explicit taxonomy for the cases where ambiguity
  matters.)
- Never frame `d_P = 1 − Pearson(triu A_a, triu A_b)` as
  "volume + topology" or as the "orthogonal pair" with `d_S`. Pearson
  on `triu(A)` is a magnitude-weighted linear correlation that mixes
  rank and amplitude; it does NOT decompose into "topology + volume",
  and cohort Spearman ρ between `T_d^(d_S)` and `T_d^(d_P)` is
  0.85–0.95 per band — they are strongly correlated, not orthogonal.
  Defensible triad: rank-only (`d_S`) / magnitude-only (`d_F`) /
  magnitude-weighted complement (`d_P`); β is the **convergence**
  cell, not "orthogonal agreement".
- Never pre-register an acceptance gate when scoping a new measure
  (IQR floor, direction-count threshold, BH-FDR cutoff, regime-width
  rule, predicted-direction one-sided test). The scope says what is
  computed and what is plotted; judging signal-vs-noise happens
  *post-hoc* with the user once the figures land. Acceptable to
  include in the scope: per-patient observation, per-patient null
  reference, per-patient z-score, cohort descriptive stats (median,
  IQR, min/max), Wilcoxon p / BH-FDR q reported as descriptive
  footnote columns. The "Acceptance" / "Cohort regularity" /
  "Verdict" / "Decision rule" sub-section does not exist. See
  `feedback_no_pre_registered_acceptance.md` (2026-05-10 sharpening
  of `feedback_regularity_over_bh_null.md` — even the IQR > 20%
  floor was an invented gate the user later objected to).
- Never lead a result writeup with "BH-FDR null / no significant
  cells" when the cohort has structured per-band / per-λ / per-phase
  shape. Describe the regularities first (medians + IQR per cell,
  per-patient direction counts, the per-band signature), apply the
  IQR > 20% × |median| floor to identify reportable shape, attach
  theoretical reading. Demote BH-FDR + Wilcoxon to a methodological
  footnote. p-values gate which cells reach the manuscript;
  cohort-shape phenomenology is what's actually in the data. See
  `feedback_regularity_over_bh_null.md` and the sister rules
  `feedback_iqr_vs_cohort_slope.md` (over-claim guardrail) and
  `feedback_dont_rerun_scalar_tests.md` (don't reach for new scalar
  tests when existing data has the picture).
- Never restrict a new investigation (e.g. epilepsy phenomenology)
  to the project's existing analysis frame by default. The trace /
  anchor / reset / rearrange (TARR) taxonomy was built for Section
  5's β-trace question; for any other target it is *one*
  connectivity-pattern family among spectral (Grassmann angles,
  eigenmode IPR), path-integral (heat kernel ρ̂(τ),
  communicability centrality, heat flux), and
  distance-distribution (cophenetic, effective resistance / commute
  time, ρ̂-leakage). Lead a new investigation plan with a
  *primitive map* (which graph object does each direction use?)
  rather than a question hierarchy lifted from Section 5. See
  `feedback_epilepsy_not_trace_locked.md`.
- **Never present any FC-derived cohort claim — raw FC, LRG
  ultrametric, KC tree-distance, Grassmann, eigenmode, partition
  metrics, anything built on the connectivity matrix — without first
  running a strength-preserving matched-strength surrogate null and
  reporting its verdict explicitly.** Within-baseline split-half
  nulls / drift triangles / cohort sampling jitter nulls are
  **exploratory diagnostics**, not verification. They reject only
  "task-rest moves more than rest-rest under cohort heterogeneity"
  and are silent on the much more dangerous alternative "matched-
  strength edge randomization reproduces this signal by construction".
  The KC β 10/10 / q=0.006 within-baseline headline collapsed on
  2026-05-11 because matched-strength produced ≈50% of the observed
  tree-distance shift on its own. See `audit_65_kc_matched_strength_verdict.md`,
  `feedback_matched_strength_mandatory.md`.
- **Never attribute an "exclude subset X → the trace strengthens /
  emerges" result to tissue X without first running a size-matched
  random-node-decimation control.** Rebuilding the LRG on a node
  submatrix changes the renormalization for *any* subset removed — the
  cophenetic ρ_split rises as the graph shrinks regardless of which
  nodes go. So the full-vs-exclude *change* is confounded with node
  count. The control: drop K = #X nodes **at random** R≥200 times,
  recompute the statistic, and compare the observed `exclude_X` value
  to that distribution (cohort `p_dec`); `< 0.05` = tissue-specific,
  `≈0.5` = generic node-count (demote the claim), `≈1` = tissue-carried.
  This is the COMPLEMENT of matched-strength (which holds per-node
  strength fixed, not node count) — both are required for an
  exclusion-contrast claim. Survival (does the trace clear
  matched-strength on the reduced graph) is a SEPARATE within-subset
  question and is not threatened by this. Same-graph **pair-class /
  pair-restriction** analyses have NO node-count confound and are
  exempt. Caught 2026-06-12: the C6 WM-X "sharpening" AND the C5 epi-X
  "strengthens/emerges" (a PRIMARY interpretive lens) were BOTH generic
  node-count, not tissue-specific (audit_85: WM/epi cophenet `p_dec`
  0.13–0.39, none `< 0.05`); only the raw-substrate WM β/γ_l (`p_dec`
  0.000/0.030) was tissue-specific. See
  `feedback_decimation_control_for_subset_exclusion.md`.
- **Never stack hardcoded patient-count / fractional-agreement /
  magnitude-ratio filters on top of a statistical test for any
  cohort verdict that ships to a writeup, manuscript, errata, or
  memory entry.** The statistical test (Wilcoxon `p<0.05`, BH-FDR
  `q<0.10`, etc.) IS the gate. Descriptive counts like
  `n_below_own_surrogate`, `n_trace`, `frac_pro_trace` are reported
  as *columns* in the per-cell output, not AND-ed into the verdict
  label. The audit_66/67 `verdict == "separated"` label stacked
  three criteria (`p<0.05` AND `|med_surr|<0.05·|med_obs|` AND
  `n_below>=8`) — the README it generated reported 0/111 separated
  cells for β even though β is §5.4's strongest matched-strength-
  controlled signal under the gate the manuscript actually cites.
  Before citing any count from an audit pipeline, grep the audit
  script for any `>= N`, `>= 0.X`, `n_below`, `frac_`, `n_sig`
  thresholds and verify the gate matches what the manuscript
  describes. See `feedback_no_hardcoded_test_thresholds.md` +
  regating audit `audit_69` (2026-05-15).
- **Never produce sycophantic answers, soft hedges, or confidence
  laundering.** Default posture is brutal scientific honesty:
  scientific questioning of every methodology in play, every null,
  every result. If a result depends on a control that hasn't been
  run, say so on the first line. If the null is weaker than the
  alternative explanations require, say so. If three iterations of
  a figure are sitting on a measure that hasn't been matched-
  strength tested, refuse to draw a fourth until it is. The user
  explicitly prefers an honest "this is currently unverified" to a
  confident headline that gets retracted. See
  `feedback_brutal_honesty_no_sycophancy.md`.
- **Never editorialize without evidence. Qualitative terms must be
  weighted with the number that justifies them.** There is no banned-
  word list — terms like *borderline, fails, clears, separated,
  promote, strengthens, marginal* are acceptable when they sit next
  to the principled-test number and per-patient structure that
  justifies them, and unacceptable as standalone verdicts. The sin
  is laundering uncertainty through an adjective that points one
  direction while the underlying numbers point another (or both
  ways). Default reporting form is still **per-patient table + exact
  principled-test `p` + LOO max-p + structural counts at operational
  thresholds (`n_obs_ρ > 0`, `n_obs_z > +2`, `n_obs_z < -2`,
  `min(obs_z)`, `range(obs_ρ)`)**; qualitative summaries on top are
  allowed when each one is tied to a specific cited number and the
  same adjective applied symmetrically to another band/result would
  not invert the verdict. With `n = 10` the principled-test `p` alone
  cannot discriminate between qualitatively distinct per-patient
  signatures (`5 strong-pro + 5 null` vs `5 strong-pro + 1 strong-
  anti + 4 null` can yield the same `n_above_surrogate`); the per-
  patient table is the load-bearing object. Caught 2026-05-26 in the
  α vs γ_low cophenet comparison: "α borderline → resolved by C5"
  vs "γ_low fails decisively" — same 5/10 count, opposite
  adjectives, no quantitative support for the asymmetry. User
  correction (2026-05-26 evening): not a banned-word problem, an
  evidence-weighting problem. See
  `feedback_no_qualitative_editorializing.md`.
- **Never write LaTeX prose or suggest manuscript edits unless the
  user explicitly asks.** The manuscript is the user's domain.
  Quantitative reports stop at the per-patient table + test `p` + LOO.
  Do not write replacement sentences, do not suggest framings, do not
  rewrite cited paragraphs. If asked for LaTeX edits, do them; never
  volunteer them.
- **Never compute a triangle scalar with the old T_d sign convention.**
  Locked 2026-05-26: every triangle scalar `T_d` MUST be
  `T_d = d(rest_pre, task) − d(task, rest_post)` so that **T_d > 0 =
  TRACE**, **T_d < 0 = ANTI-TRACE**, **T_d = 0 = NO TRACE**. Applies at
  every layer (raw FC, LRG D_coph, KC, Grassmann, eigenmode E1).
  Wilcoxon one-sided trace tests use `alternative='greater'`. Per-
  patient trace counts are `(T > 0).sum()`. Surrogate upper-tail p is
  `mean(s_finite >= obs_T)`. Plot ylabels / titles say
  `positive = trace` (never `negative = trace`). Do not reintroduce
  `d(task, rsPost) − d(rsPre, task)` in any new or refactored compute
  site. See `feedback_td_sign_convention.md`.
- **Never surface T_d triangle scalars or `d_P`/`d_F` variants in
  preprint context.** Locked 2026-05-26: for raw FC and LRG D_coph
  layers, the preprint locked probe is `ρ_split` (the Spearman
  correlation `Spearman(Δ_task, Δ_rest)` from audit_33 / audit_63 /
  Methods CTM block) under the C3 matched-strength gate. `ρ_split` is
  intrinsically rank-only (it IS a Spearman correlation by
  construction) — there is no Pearson or Frobenius variant of
  `ρ_split` to drop. The audit_25 / audit_35 *triangle scalars*
  `T_d^(d_S)`, `T_d^(d_P)`, `T_d^(d_F)` are a separate, three-
  distance internal diagnostic family at a different layer of the
  ladder; they do NOT appear in preprint tables, prose, per-band
  summaries, figures, captions, verdict ledgers, or any chat/report
  Claude writes about the preprint cross-band picture. The
  /tmp/td_loo.py three-distance LOO table is internal diagnostic.
  Manuscript-facing β raw FC = `ρ_split^raw` median +0.258 (p=0.053);
  β D_coph = `ρ_split^coph` median +0.221 (p=0.005). When user says
  "focus on Spearman" they mean `ρ_split` (already Spearman by
  construction), not "pick d_S within the T_d triangle". See
  `feedback_preprint_d_s_only.md`.
- **Never name a library module / function after a manuscript-local
  token** (`section3`, `figureN`, `preprint`, `chapter`, `H2c`, etc.).
  Library names under `src/lrg_eegfc/` reflect general graph / network /
  statistics / I/O concepts (`network_layouts`, `network_drawing`,
  `tree_metrics`, `surrogate_helpers`, `patient_io`). A helper used by
  Section-3 figures today must be importable from Section-7 figures
  tomorrow without renaming. Locked 2026-05-28. See
  `feedback_library_names_general.md`.
- **Never use `imshow_colorbar_caxdivider` for a colorbar shared
  across multiple columns in the same row.** The helper attaches the
  colorbar to a single axis via `make_axes_locatable`; it has no
  notion of a multi-column shared cbar and will misalign or resize
  the wrong axis. For row-shared cbars keep the explicit
  `make_axes_locatable` / `fig.add_axes([...])` pattern. The helper is
  the canonical choice for any *single-imshow* axis. Locked 2026-05-28.
  See `feedback_imshow_colorbar_caxdivider_scope.md`.
- **Never use the bare letter `β` for a regression slope, OLS
  coefficient, or any non-band-name numerical quantity.** The letter
  is reserved for the 13–30 Hz EEG band throughout the codebase,
  manuscript, CSV files, figure axes, prose, and memory. The
  asymmetric pair-trace recovery slope `⟨Δ_task, Δ_rest⟩ / ‖Δ_task‖²`
  is named `s_TR` (math, `s_{\mathrm{TR}}`) / `slope` (code, CSV
  columns, helper args) / "asymmetric pair-trace slope" or "recovery
  slope" (prose). The helper at `src/lrg_eegfc/utils/metrics/hypothesis.py`
  is `regression_slope_through_origin`, NEVER `beta_slope`. Same rule
  for any future asymmetric measure — pick a letter that does not
  collide with the seven band symbols. Locked 2026-05-30. See
  `feedback_no_beta_for_regression_slope.md`.

## Always

- Always show ≥3 patients / bands / phases in published figures (no
  Pat_02-only / beta-only / rest_pre-only plots).
- Always route data loading through `workflow.fc.load_fc_matrix`.
- Always import statistical helpers from
  `lrg_eegfc.utils.metrics.hypothesis` (wilcoxon_z, bh_fdr,
  rank_biserial, boot_ci_mean).
- Always set dendrogram y-limits as
  `tmin = merge_heights[0]*0.8, tmax = merge_heights[-1]*1.05`.
  Never 0.5× / 2.0×.
- Always zoom nilearn glass-brain panels to electrode bbox, not full
  default brain.
- Pat_03 is acquired at 1024 Hz (others at 2048 Hz). Sampling-rate
  handling is **config-level only** — `nperseg_for_fs(fs)` and
  `FS_OVERRIDES` at `config/const.py` adapt the spectral estimator to
  the actual sampling rate. Pat_03 is a full cohort member at n = 10
  and is treated identically to every other patient at the analysis
  level. **Do not** run Pat_03-dropout sensitivity tests, **do not**
  mark Pat_03 distinctly in figures, **do not** report Pat_03 values
  separately in tables, **do not** describe Pat_03 as an "outlier" or
  "negative control". (Softened 2026-05-18: the previous "always flag
  Pat_03 as 1024 Hz outlier (negative control)" framing was retired
  on user direction — the sampling-rate difference is fully absorbed
  by the config layer.)
- **Patient-dropout policy: don't drop patients in analysis.** The
  default cohort for every test is the full `n = 10`. Two narrow
  exceptions are allowed:
  1. **A single patient anti-aligned at the cohort level** for the
     specific probe under test, biology-driven (not magnitude-driven).
     The canonical example is Pat_15 at the LRG layer for β (right-
     hemisphere-only implant, no epileptic contacts, the only patient
     with `T_G(k=40) > 0` and the only pro-cohort-direction substrate
     value above noise floor). Reported as an `n=9` sensitivity row,
     not the default cohort.
  2. **Genuinely problematic data** (vendor corruption, sampling-rate
     difference handled at config layer, etc.). None currently active
     in the n=10 cohort.
  **Retired dropouts (do not reintroduce):**
  - **Pat_03 dropout** (1024 Hz "outlier") — retired 2026-05-18 (am).
    Sampling-rate handled at config layer only.
  - **Pat_07 dropout** (substrate "anti" at `S = +0.005`, near-zero) —
    retired 2026-05-18 (pm). Pat_07 is solidly pro at both LRG probes
    (`ρ_split^coph` z = +4.47; Grassmann `T_G(k=40)` z = −3.90).
    Dropping Pat_07 was inherited from the §4 substrate-layer analysis
    and has no LRG-layer rationale.
  - **Legacy `n=8` "pro-cohort restriction"** (drop Pat_07 + Pat_15) —
    retired 2026-05-18 (pm), replaced by the LRG-native `n=9` drop of
    Pat_15 only.
  General rule: "anti-direction" must mean the *probe* under test
  registers the patient as substantially anti-aligned (not "marginally
  anti-aligned at a different probe"). Apply this when justifying any
  future `n < 10` restriction.
- Always write a renormalization-style head before any body (see
  `renormalization-style.md`).
- Always file new task-trace / task-persistence investigation tooling
  under `.agents/guides/task-persistence-investigation/<YYYY-MM-DD>_<measure>.md`
  as a mathematically rigorous scope report **before writing any code**
  (notation → predicates → properties → caveats → pseudocode →
  visualization → prior-tool connection → open questions). See that
  folder's `README.md` for the required structure.
- Always lead substantive answers (more than ~5 sentences) with a
  table of contents of 3-7 plain-English declarative statements; then
  expand each below. Reader skims TOC, picks what to read in detail.
- Always expand abbreviations on first use in any chat reply (VI =
  variation of information, BH-FDR = Benjamini-Hochberg false discovery
  rate, PRL = partition-resolution-locked, KC = Kendall-Colijn, FC =
  functional connectivity, etc.). Don't drop jargon the user hasn't
  seen written out.
- Always use the **trace / anchor / reset / emergent** taxonomy when
  describing how a module behaves across phases *in any context where
  the reader could otherwise be confused between the four phenomena*
  (cross-phase taxonomy tables, mixed-band paragraphs, cross-phase
  summary captions, variable names in code: `n_trace` instead of
  `n_persist`). See `.agents/guides/01_project/terminology.md`. In
  unambiguous task-trace sections (e.g. the β results subsection)
  bare "persistence" is acceptable per 2026-05-18 softening — the
  taxonomy is the disambiguation tool, not a global ban on the word.
- Always activate the project matplotlib style at the top of every
  figure script:
  `from lrg_eegfc.visuals.styles import use_lrg_style; use_lrg_style()`.
  Single source: `src/lrg_eegfc/visuals/styles/lrg_eegfc.mplstyle`
  (font sizes, tick widths, `pdf.fonttype=42` TrueType embed). Do
  not redefine these defaults inline; for one-off overrides use
  `with rc_context({...}):`, never edit the mplstyle for a single
  figure. On any LogNorm colorbar, end with
  `_apply_factored_sci_format(clb, axis_orientation=...)` so inline
  `2 × 10ⁿ` minor-tick labels never appear (full rule:
  `.agents/guides/05_plotting/colormaps-and-styles.md`).
- **Always think deeply about the test before running it. Open every
  new test / null / methodology with a 5-point critical preamble
  (script docstring or scope doc, before any code):**
  1. *Claim:* the exact thing to be tested at the cohort and/or
     per-patient scale, in one sentence.
  2. *Null:* the specific null hypothesis being run.
  3. *Strongest alternative:* the most plausible non-trace explanation
     the null *should* control for (strength heterogeneity, drift,
     volume conduction, cohort sampling jitter, edge-identity
     coincidence, etc.).
  4. *Does the null actually control for it?* — mechanism, not vibes.
     If matched-strength preserves the strength distribution, does it
     also preserve the artifact you're worried about? If split-half
     reshuffles time, does it reach the relevant alternative? Spell
     out what the null *cannot* reject.
  5. *Falsification + limitations:* what outcome would falsify the
     claim, and what limitations remain even if the test passes
     (e.g. independent per-phase rewiring decorrelates cross-phase
     identity by construction → ρ_split surviving does not isolate
     shape from strength; it only rejects the matched-strength null).
  Non-optional. The point is to catch the KC-style mistake at the
  design stage, before two weeks of figures sit on a null that doesn't
  reach the relevant alternative. See
  `feedback_critical_null_preamble.md`.
- **Always state limitations of unverified methods up front, in the
  first paragraph of the result writeup, not buried in a "caveats"
  appendix.** Any measure not yet tested against matched-strength
  surrogacy is marked **"unverified"** in writeups, memory, and chat,
  and cannot be cited as the load-bearing claim. The phrase
  "this currently has no matched-strength control" should appear in
  every preliminary writeup of any FC-derived cohort claim.
- **Always pair every cohort-level p-value with a leave-one-out (LOO)
  sensitivity report.** Wilcoxon at `n = 10` is robust to magnitude
  outliers but vulnerable to direction outliers — a single very-
  pro-cohort patient can shift `p` from `0.5` to `0.005`. For every
  Wilcoxon-based gate (C1, C2, C3, C4, C5; anatomy A1, A3) and every
  cluster-permutation cluster_p_mass, also report the LOO max-p (worst-
  case after dropping each patient once) and the identity of the
  most-influential patient. The LOO is a **trust diagnostic** for
  the reader, not a gate — never threshold it (no "LOO max-p < 0.05"
  cutoff). When LOO max-p crosses 0.05, the verdict stays as the
  original test result, but the manuscript text must flag it
  explicitly (e.g., *"p = 0.005 overall, p = 0.07 after dropping
  Pat_XX"*). See `feedback_no_single_patient_p_driven.md`.

## Meta-rule

User says "never X" or "always Y" → append to this list on first
mention AND save a `feedback_<short>.md` memory. This list is the
single source of truth; CLAUDE.md mirrors it inline so every agent
session sees the rules in context.
