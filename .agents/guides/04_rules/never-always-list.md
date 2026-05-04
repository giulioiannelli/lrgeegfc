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
- Never delete files that document research history — `git mv` to
  `<parent>/archive/YYYY-MM/` instead.
- Never pool metrics into a consensus scalar (explicitly forbidden by
  the user).
- Never skip frontmatter on a new `.agents/` .md file.
- Never invent metric names — cite literature or existing code.
- Never plot Δ_ARI(k) in partition-multiscale band×k publication
  figures — canonical 3 are `Δ_VI, Δ_H, Δ_NMI`; CSV may retain `d_ARI`.
- Never use channel-label letter prefixes (A/B/.../Q) for cohort-level
  implant analysis — they are arbitrary clinical labels with no
  cross-patient anatomical meaning. Use `(x, y, z)` coordinates and
  Desikan-Killany regions from `implant_pat_NN.csv` instead.
- Never title a figure with `fig.suptitle` (duplicate of above for
  emphasis — it keeps slipping through).
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
- Never call our `T_d < 0` / `d(TT, RPost) < d(RPre, TT)` finding the
  bare "persistence". It is a **trace** (task reorganized AND change
  stuck) — distinct from **anchor** (module unchanged across all
  phases), **reset** (task changed AND came back), **emergent**
  (module that did not exist in `RPre`). Use the taxonomy at
  `.agents/guides/01_project/terminology.md`. Bare "persistence"
  silently flips the reading from trace to anchor.
- Never frame `d_P = 1 − Pearson(triu A_a, triu A_b)` as
  "volume + topology" or as the "orthogonal pair" with `d_S`. Pearson
  on `triu(A)` is a magnitude-weighted linear correlation that mixes
  rank and amplitude; it does NOT decompose into "topology + volume",
  and cohort Spearman ρ between `T_d^(d_S)` and `T_d^(d_P)` is
  0.85–0.95 per band — they are strongly correlated, not orthogonal.
  Defensible triad: rank-only (`d_S`) / magnitude-only (`d_F`) /
  magnitude-weighted complement (`d_P`); β is the **convergence**
  cell, not "orthogonal agreement".

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
- Always flag Pat_03 as 1024 Hz outlier (negative control) — include
  in analyses, mark distinctly in figures, report separately in tables.
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
  describing how a module behaves across phases. See
  `.agents/guides/01_project/terminology.md`. Replace `n_persist` with
  `n_trace` in new code and figure annotations.

## Meta-rule

User says "never X" or "always Y" → append to this list on first
mention AND save a `feedback_<short>.md` memory. This list is the
single source of truth; CLAUDE.md mirrors it inline so every agent
session sees the rules in context.
