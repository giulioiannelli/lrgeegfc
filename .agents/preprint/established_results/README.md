---
name: established-results-folder-readme
era: IMCOH_ABS_COHORT_N10
status: historical-stub
status_updated: 2026-05-28
kind: index
scope: aspirational frozen-claims system; superseded in practice by bands/ briefs + locked/ ledgers
---

# `.agents/preprint/established_results/` — historical methodology Q&A

> **Note (2026-05-28).** This folder was designed as a frozen-claim
> provenance system (one .md per number in the preprint, with full
> reproducible computation chain). In practice the **bands/** briefs and
> **locked/** ledgers absorbed that role — every number in the manuscript
> is now traceable via `band brief → CSV row → audit script` without
> needing a separate established-results file per claim.
>
> The folder is kept as a **historical methodology Q&A archive** for
> open-question files (e.g., `00_open_methodology_question_lrg_D_convention.md`
> which is `status: withdrawn` and documents the Pipeline-1-vs-Pipeline-2
> cophenet normalization episode). The aspirational example list that
> previously appeared here (`beta_rho_split_within_baseline.md`,
> `beta_rho_split_matched_strength.md`, `beta_grassmann_window.md`,
> `beta_grassmann_epi_X.md`) **was never populated** and is dropped
> because those numbers live in `bands/01_beta.md` + the locked CSVs.
>
> Live source of truth for verdict-bearing numbers: `bands/<band>.md`
> brief + `locked/VERDICT_LEDGER.md` + `locked/ANATOMY_LEDGER.md` +
> the `data/audit/*/cohort_summary.csv` row each cites.

## Current contents

- `00_open_methodology_question_lrg_D_convention.md` — `status: withdrawn`.
  Documents the resolved Pipeline-1 (normalized cophenet) vs Pipeline-2
  (unnormalized cophenet) discrepancy at β ρ_split^coph. The cohort
  verdict is invariant to the within-patient monotone D rescaling;
  Pipeline-1 (normalized) is the locked convention.

## Cross-references

- `.agents/preprint/README.md` — top-level preprint folder index
- `.agents/preprint/locked/VERDICT_LEDGER.md` — locked verdicts + decisions log
- `.agents/preprint/locked/ANATOMY_LEDGER.md` — locked anatomy verdicts
- `.agents/preprint/bands/` — per-band briefs (the live verdict source)
- `.agents/preprint/METHODS_AUDIT_ISSUES.md` — live methods-section audit items
