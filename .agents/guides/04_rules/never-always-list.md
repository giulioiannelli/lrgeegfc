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
- Never place repository information or code / implementation
  technicalities in the paper — no numba/JIT, run-times or speedups
  ("98× faster"), bit-identical / reproducibility notes, code-QA
  parentheticals ("verified to 1e-6 per surrogate"), script / audit
  names, or cache paths. The paper reports the scientific method + result;
  that detail belongs in the repo, not the manuscript. Keep the
  statistical parameters (e.g. `R = 200` surrogates, strength preserved to
  machine precision); drop how-it-ran-fast and how-it-was-verified-in-code.
  (Added 2026-07-07 — user deleted the numba "Reproducibility" paragraph
  from the methods; `feedback_no_repo_or_code_technicalities_in_paper`.)
- Never start a new scalar hypothesis test when the signal is visible
  in existing VI(k) / partition-multiscale / H2d-θ artifacts — surface
  the existing signal first.
- Never collapse a τ/scale-resolved result to a scale-max (or any
  single scalar) and report that scalar as the verdict. ALWAYS report
  and PLOT the per-scale τ-dependence first; a scale-max is a
  conservative supplement only, and valid only when the null is a
  proper distribution — for a single-realization null (drift arc) the
  fair test is per-scale paired `obs(s) vs null(s)` at the scale where
  the signal emerges, never max-over-scale. And OPEN the rendered
  figure before stating what it shows — never infer the shape from a
  CSV or from memory. (Added 2026-07-12 — I reported the backbone drift
  null as "β = drift, all bands fail" from an unfair scale-max gate;
  per-scale, β is above drift at every scale and high_γ beats it at
  coarse scales. `feedback_report_tau_dependence_no_scalar_collapse`.)
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
- Never run or present a spatial-/proximity-matched null on `|ImCoh|`
  (`imcoh_abs`) findings, and never frame physical proximity as a
  connectivity confound. `|ImCoh|` is zero-lag-immune (Nolte 2004) so
  proximity **cannot** enter connectivity; a spatial null wrongly
  removes **real** local connectivity — it destroys signal, not a
  confound, and is the wrong control (not "the hardest null"). The only
  mandatory null is **matched-strength**. For a SOZ marker the
  CONNECTIVITY claim = all-contacts AUC vs matched-strength (audit_132,
  multiband δ/β/low-γ); off-shaft leave-one-shaft-out answers a separate
  CLINICAL question ("finds SOZ a distance ruler can't") and must NOT be
  the default headline — it over-handicaps the connectivity claim.
  (Locked 2026-06-22, refined 2026-06-23; `feedback_no_spatial_proximity_null`.)
- Never delete files that document research history — `git mv` to
  `<parent>/archive/YYYY-MM/` instead.
- Never pool metrics into a consensus scalar (explicitly forbidden by
  the user).
- Never present the ρ_sym cophenet **cohort gate** as cleanly separating
  trace bands from non-trace bands. At n=10 the signed-rank p is discrete
  and γ_l/γ_h/δ sit at p=.080 (~1 patient from the p<.05 line); β itself
  is only 7/10 positive and clears on **magnitude**, not count. ρ_sym also
  COMPRESSES the separation vs ρ_split (α/β less extreme, γ/δ pulled up to
  marginal). Present the **graded tiers** (clear α/β · marginal γ_l/γ_h/δ
  · absent θ) + the per-patient spread — never a clean "α/β trace vs rest
  no-trace" binary. β's flagship status rests on magnitude + Grassmann +
  OFC convergence, NOT this single gate. (Locked 2026-07-06;
  `feedback_rho_sym_gate_marginal_not_clean_bands`.)
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
- **Never hardcode per-band colours in a figure.** Every plot that
  distinguishes frequency bands by colour MUST call `band_color(band)`
  (or index `BAND_COLORS`) from `lrg_eegfc.visuals.styles` (re-exported
  from `lrg_eegfc.visuals`), so the whole codebase restyles from one
  constant. The canonical palette is a spectrum — slow bands red, fast
  bands blue, rainbow between (`turbo` reversed; the single swap point is
  `BAND_CMAP_NAME` / `BAND_CMAP_SPAN`). Use `band_color(band, shade<1)`
  to darken the bright midtones (alpha, beta) for text / thin lines.
  Locked 2026-07-09. See `feedback_band_color_palette.md`.
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
- Never truncate / "length-match" a recording (e.g. head-truncate
  `task_test` to `task_learn` length). No analysis series is ever
  shortened. The truncation-based duration null (audit_113b/161) is
  **retired-invalid** — a contiguous window injects a common-mode shift
  into `D_TT` and the δ negative-control false-positives. The clean,
  locked duration control is **duration-ratio regression** (regress the
  per-patient effect on the `task_test/task_learn` ratio; non-destructive,
  uses all data), plus the fact that encoding (`e = c_learn − c_pre`) is
  length-symmetric by construction. See `feedback_no_series_truncation.md`,
  `arc_inference_consolidation_2026_06_18.md`. Locked 2026-07-13.
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
- **Never state a result on the raw functional-connectivity substrate.**
  Results are given ONLY in the introduced LRG/Laplacian framework —
  the cophenetic communication distance `D_coph` (`ρ_split^coph`) and
  the Grassmann leading-mode subspace (`d_G(k)`/`T_G`). Raw FC (raw
  |ImCoh|, `d_S`) is the **comparison baseline only** ("does the LRG
  layer add over raw?"), never a result and never load-bearing evidence
  for a claim. If a finding is significant on raw but only a *trend* on
  the cophenetic/Grassmann objects, the honest verdict is the
  cophenetic/Grassmann trend, NOT the raw significance — do not write
  "decisive on raw, a trend on cophenet" as if raw rescues the claim.
  Locked 2026-06-18 (user directive). Caught the same day: the epi
  "core spared" claim had been leaned on raw `p_pair`=1.0; on the
  cophenetic framework it is only a trend (`p_pair`=0.945) + matched-
  strength `weaken`, so the real result is gray-dominance (gray↔gray
  cophenetic pair-count hotspot), not "core spared". See
  `feedback_results_only_in_laplacian_framework.md`.
- **Never switch the core Laplacian operator without explicit user
  sign-off.** The analysis operator is the combinatorial Villegas
  "fluid" Laplacian `L̂ = D̂ − W` and its propagator `e^{−τL̂}` — the LRG
  framework. Do NOT substitute a symmetric-normalized
  (`I − D^{-1/2}WD^{-1/2}`), random-walk (`I − D^{-1}W`), signed, or
  magnetic Laplacian to "debias strength" or for any other reason on
  your own initiative — a result on a different operator is not an LRG
  result and leaves the framework. When the combinatorial `D_τ` is
  strength-confounded, address it with the matched-strength null
  *inside* the framework, never by changing the operator. Offering an
  operator change as an option ≠ authorization to run it. Locked
  2026-07-07 (user washed out an unauthorized symmetric-normalized
  diffusion detour). See `feedback_combinatorial_laplacian_only.md`.
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
- **Never source a load-bearing (especially unwelcome) number from a
  previous agent's report, a cached summary table, or already-written
  prose — recompute it fresh from the raw substrate with code you can
  inspect, and use the cache only as a cross-check, never as the
  source.** Prior agents and prior writeups hallucinate and mis-derive.
  Caught 2026-07-10: two subagents reconstructing audit_166 disagreed on
  cophenetic's per-patient β-inference breadth (6/10 vs 10/10) — at least
  one fabricated a number about to enter a verdict. Any result that gates
  a claim, contradicts the user's expectation, or will touch the preprint
  must be re-derived by independent computation (reuse validated library
  primitives + the FC matrices, not the derived CSV), print intermediates,
  sanity-check one cell by hand, and diff against the cache; if they
  disagree the cache was wrong and the fresh number stands. User directive
  (2026-07-10): "never rely on text already written... don't rely on
  previous agents, remember that we could have made hallucinations." See
  `feedback_fresh_recompute_not_prior_text.md`.
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
- **Never compress a band/cohort verdict to a single number (cohort
  median, gate `p`, pass/fail tag, or a one-glyph figure marker) when
  the per-patient fluctuations are large — surface the spread, it may
  carry the result.** This extends the editorializing rule above from a
  reporting-hygiene point to a *scientific* one: large cross-patient
  heterogeneity is itself a candidate finding (the trace may concentrate
  in a subset of patients — possibly implant-coverage / anatomy-driven —
  not a uniform cohort effect), never noise to average away. Do not say
  "trace / no trace" with ease. Evidence (2026-06-25 figure session):
  cophenetic α *passes* (median Д=0.091, 9/10) while γ_l *fails* (median
  Д=0.081) though γ_l carries **3 patients above β's max** effect —
  magnitude-blind sign-consistency; on Grassmann γ_h has the **largest**
  median per-patient effect (z≈9, 5 above β's max) yet is "no trace", and
  a naive per-patient collapse flips the band ordering. In FIGURES: plot
  the 10 per-patient effects so the spread is the dominant visual; the
  verdict is subordinate annotation, never a hiding glyph (the
  forest-corner-glyph antipattern that made δ/γ_l/γ_h look like β). The
  locked verdicts still stand as the correct cohort tests; this governs
  how they are *reported and drawn*. Builds on
  `feedback_no_qualitative_editorializing` +
  `feedback_no_single_patient_p_driven`. See
  `feedback_fluctuations_are_signal.md`.
- **Never narrate the study or use essayistic / literary register in
  paper prose — present results, not a story.** When writing manuscript
  prose the user HAS asked for, three registers are banned (PI 2026-07-06,
  "never ever from now on"): (1) stating an absence — "there is no behavioral
  record", "task performance was not collected / cannot be obtained", "X was
  not available"; say what IS available, never what is not, and NOT even a
  terse caveat in Methods (PI sharpened 2026-07-06). *Enact* a ceiling by not
  overclaiming (make no behavioral claim); convey an absence by silence. (2) essayistic meta-framing of our own reasoning —
  "Two cautions keep that reading within its evidence", "The second caution
  is more fundamental and bounds the whole section"; state the finding, then
  mark interpretation with a plain "we interpret X as Y". (3) literary /
  poetic register — "mundane alternative", "licensed by", "a representation
  we decode", metaphors like "echo". Register = Nature Neuroscience results
  section. Caught 2026-07-06 in R2.1 para 2. See
  `feedback_results_not_story.md`.
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

- **Always frame anatomical localization as an OVEREXPRESSION /
  ACCUMULATION hotspot on top of a trace distributed across the whole
  network — never as an exclusive container or as "absent elsewhere".**
  "β trace → OFC" means OFC accumulates *more* trace than matched-strength
  predicts (within-OFC lean ≈0.59 — most trace mass is OUTSIDE OFC), not
  that the trace lives only in OFC. Report concentration *above the
  distributed baseline* (per-patient demeaned, matched-strength R=1000),
  never exclusive region membership. Language: "concentrates / accumulates
  / overexpressed in S", never "localized to S exclusively" / "absent
  outside S". Figures show the distributed field across all nodes with the
  hotspot glowing on top, never a binary region mask implying emptiness
  elsewhere. Carries into the β inference-mark localization and every
  future localization. Locked 2026-06-18 (user directive). See
  `feedback_localization_is_overexpression_not_container.md`.
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
- **Always optimize heavy compute, time-estimate it on a small batch
  first, and surface live progress + timing.** Three-part discipline for
  any pipeline that could run more than ~1 minute (surrogates, sweeps,
  per-patient×per-band grids, LOO, bootstraps):
  1. **Optimize before launching, don't after.** Hot numeric loops
     (surrogate shuffles, cophenetic reconstructions, rank stats) get a
     numba `@njit` / vectorization / cached-eigendecomposition pass FIRST.
     A version that would run for hours unoptimized is not acceptable to
     launch — the numba shuffle is 98× over pure Python and bit-identical
     (`feedback_numba_for_surrogates`). Extend that reflex to every heavy
     kernel, not just the surrogate shuffle.
  2. **Time a small batch, then extrapolate — never launch a long run
     blind.** Run 1–2 units (one patient, one band, R=5) wall-clock-timed,
     multiply out to the full job, and state the estimate BEFORE launching
     the full pipeline. If the extrapolation is hours, stop and optimize
     (rule 1) or shrink the job.
  3. **Surface progress + timing in the code.** Long scripts print a live
     `[i/N] label  elapsed=…s  ETA=…s` line per unit with **`flush=True`**
     (block-buffered stdout under redirection hides un-flushed prints —
     always flush), plus a total wall-clock at the end. The user must be
     able to see where the run is at any moment, not guess. See
     `feedback_optimize_time_and_surface_progress.md`.

## Drift is orthogonal to matched-strength (check both)

Cross-phase / arc trace claims must clear a **within-recording drift null**, not
only matched-strength — MS preserves node strength but ignores temporal structure,
so a slow-drift-driven trace passes MS unscathed ("raw FC always shows trace"). A
drift control already exists: **C2** in `locked/CONTROLS.md` (`ρ_split > ρ_drift`),
covering the **whole-task** trace (α/β/γ_low pass). It does **not** cover the
four-phase arc functionals (`T_infspec·e` inference-specific, `T_learn` encoding);
those need their own **duration-matched** drift gate. `audit_167` (2026-07-10):
whole-task β drift-robust under a 2nd construction; α/γ_low construction-dependent;
the **inference-specific decomposition is drift-unverified** (β p=0.22); raw-FC
traces are **entirely drift**. **Grep `locked/CONTROLS.md` before declaring a
control missing.** See `feedback_drift_null_mandatory.md`.

## ρ^coph ALWAYS means ρ_sym — never the non-symmetric raw form

**Locked 2026-07-13.** The notation `ρ^coph` — with or without a `sym`
subscript, in code, figures, captions, prose, or memory — ALWAYS denotes the
symmetric, split-baseline-subtracted estimator

```
ρ_sym = ½[ Spearman(D_task − D_preA, D_post − D_preB)
         + Spearman(D_task − D_preB, D_post − D_preA) ]
```

Never report, plot, or bar the **non-symmetric raw** cophenetic correlation
`Spearman(D_X[triu], D_test[triu])` (a bare phase→phase similarity such as
"phase → task_test") as a result or as a `ρ^coph` value again. The raw form
(a) is **not baseline-subtracted**, so it is dominated by the shared anchor
structure any two phases trivially have in common — its matched-strength null
is HIGH (~0.2), not ~0.01 — and (b) is **asymmetric** in the arbitrary A/B
split-half assignment. ρ_sym removes both.

Consequences for every ρ^coph axis/bar: its reference is the matched-strength
**null value** (mark it), and its **reachable maximum is the split-half
reliability** `ρ(A,B) ≈ 0.5` (β), NEVER 1 — anchor the bar to that reachable
ceiling, cut the unreachable region, and mark the null; never draw a full bar
to 1.0. Supersedes the raw definition previously given in
`feedback_no_partition_metrics_use_rho_coph.md`. See
`feedback_rho_coph_always_rho_sym.md` and
`feedback_rho_sym_canonical_estimator.md`.

## Meta-rule

User says "never X" or "always Y" → append to this list on first
mention AND save a `feedback_<short>.md` memory. This list is the
single source of truth; CLAUDE.md mirrors it inline so every agent
session sees the rules in context.
