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
- Never title a figure with `fig.suptitle` (duplicate of above for
  emphasis — it keeps slipping through).

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
- Always exclude Pat_14 task_test (file corrupt).
- Always drop rows [53, 54, 55] of Pat_10 task data at load.
- Always write a renormalization-style head before any body (see
  `renormalization-style.md`).

## Meta-rule

User says "never X" or "always Y" → append to this list on first
mention AND save a `feedback_<short>.md` memory. This list is the
single source of truth; CLAUDE.md mirrors it inline so every agent
session sees the rules in context.
