---
type: report
status: current
date: 2026-05-08
era: IMCOH_ABS / COHORT_N10
scope: index for the full-verification package on notes_imcoh.tex (PDF .agents/references/notes_imcoh_260508.pdf)
---

# Full verification package — index

**Head.** Per-section verification + cross-cutting fix files for `notes_imcoh.tex` (PDF at `.agents/references/notes_imcoh_260508.pdf`). Per-section files (`01`–`06`) hold numerical / framing / caveat actions for the writing agent. Cross-cutting files (`10`–`16`) hold notation drift, terminology disambiguation, the m=48 overload, figure regeneration / redesign queues for the coding agent, new analyses owed, and the epilepsy-track pointer. Cohort is **n=10** (Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15) since 2026-04-25; Pat_03 is the 1024 Hz outlier; Pat_14 was vendor-replaced.

## Contents of this directory

Per-section verification (consumed first by the writing agent):

- `01_section1_intro_fixes.md` — §1 Introduction
- `02_section2_fc_fixes.md` — §2 FC estimation (still n=6 in figures — coding-agent action)
- `03_section3_lrg_fixes.md` — §3 LRG pipeline + comparison
- `04_section4_raw_fc_fixes.md` — §4 raw |ImCoh| phase-distance audit
- `05_section5_lrg_fixes.md` — §5 multiscale Laplacian phase-distance investigation
- `06_section6_synthesis_fixes.md` — §6 synthesis + outlook

Cross-cutting fix files (this batch):

- `10_notation_inconsistencies.md` — small notation drift across sections
- `11_terminology_disambiguation.md` — **CRITICAL** trace-direction / trace-module / trace-leaf overload
- `12_bonferroni_m48_overload.md` — two distinct `m=48` families
- `13_figures_to_regenerate.md` — coding-agent regeneration queue (cohort fixes; §2 still n=6)
- `14_figures_to_redesign.md` — coding-agent redesign queue (visual quality, NOT cohort)
- `15_new_analyses_owed.md` — ranked list of analyses owed
- `16_epilepsy_pointer.md` — parallel epilepsy track (NOT folded into manuscript)

## Connected scope reports (read for context, do not modify)

`.agents/guides/task-persistence-investigation/`:

- `2026-04-25_task-trace-canonical.md` — TARR taxonomy + P/T/R/RA reformalization
- `2026-05-07_kc-leaf-topological-memory.md` — KC trace-leaf scope
- `2026-05-08_epi-cross-phase-rigidity.md` — Direction A scope (parallel track)
- `2026-05-08_epi-eigenmode-localization.md` — Direction C scope (parallel track)

## Master writeups (current, read for context, do not modify)

- `.agents/reports/archive/2026-05/2026-05-05_result-2-lrg-beta-trace.md` — Result 2 source
- `.agents/reports/2026-05-07_kc-cross-phase-taxonomy.md` — §5.6 taxonomy source
- `.agents/reports/2026-05-08_section5_joint_bh.md` — joint-BH family arithmetic
- `.agents/reports/2026-05-08_section5_5_verify.md` — anatomy m=48 source
- `.agents/reports/2026-05-08_section5_6_verify.md` — taxonomy class counts
- `.agents/reports/2026-05-08_section5_4_verify.md` — VI(k) + Grassmann verification
- `.agents/reports/2026-05-08_psi_argmax_pinning.md` — Ψ-irrelevance numbers
- `.agents/reports/2026-05-08_kc_beta_pvalue_check.md` — KC β p-value cohort agreement
- `.agents/reports/2026-05-08_kc_crossprobe_check.md` — KC cross-probe restriction
- `.agents/reports/2026-05-08_ctm_lowgamma_counts.md` — CTM low-γ patient counts
- `.agents/reports/2026-05-07_epileptic-n10-revisit.md` — n=10 epi verdict
- `.agents/plans/active/2026-05-08_lrg-epilepsy-research-directions.md` — 4-direction epi plan

## Consumption order

1. **Writing agent** reads `01`–`06` (per-section) + `10`, `11`, `12` (notation, terminology, m=48). Apply rewrites in order.
2. **Coding agent** reads `13`, `14` (regenerate + redesign queues) and the relevant per-section files for figure-cohort actions. After regeneration the writing agent re-checks captions for any cohort-relabel.
3. **Both agents** read `15` (new analyses owed) for the deferred-controls + multiscale-taxonomy roadmap, and `16` (epilepsy pointer) to know **not** to fold the parallel track into the manuscript.
4. **Scope reports** (`task-persistence-investigation/`) are reference documents for new pipelines that come out of `15`; they are not direct inputs to the verification cycle.

## Single-line invariants for any agent touching the manuscript

- **Cohort = n=10** since 2026-04-25 (Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15). §2 figures are the only place where n=6 still leaks in — see `13_figures_to_regenerate.md`.
- **Trace** is overloaded; default to the precise term (direction / module / leaf) on first use in each subsection — see `11_terminology_disambiguation.md`.
- **m=48** appears in two unrelated places (§5.3 LRG joint-BH; §5.5 anatomy Bonferroni) — see `12_bonferroni_m48_overload.md`.
- **Pat_03** is the 1024 Hz outlier; included in cohort summaries, plotted distinctly, dropout robustness reported.
- **The epilepsy track is parallel and not folded into this manuscript** — see `16_epilepsy_pointer.md`.
