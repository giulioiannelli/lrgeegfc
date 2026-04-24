---
name: start-here
type: guide
era: COHORT_N9
status: current
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

# Start Here

This file is the entry point for agents. It points to the current state of
the scientific pipeline and the key reports that carry the load-bearing
numerical results.

## Current scientific state (2026-04-24)

- **Cohort**: 10 sEEG patients (`Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15`).
  Cross-phase tests run on n=9 (Pat_14 excluded — vendor `task_test.mat`
  corrupt). H3 and H4 run on n=10.
- **Canonical FC metric**: `imcoh_abs` = `⟨|ImCoh(f)|⟩_f` (Ewald 2012
  convention of the Nolte 2004 imaginary coherency). Volume-conduction
  immune. MSC is retained as a method baseline only; don't use it for the
  current hypotheses.
- **Canonical report for the writing agent**: `.agents/reports/2026-04-24_multiscale-task-trace.md`.
  Supersedes `2026-04-24_imcoh-results-for-writing.md` (n=5 era, flagged stale).

## Where to begin

1. **If you are writing a paper section**: read
   `.agents/reports/2026-04-24_multiscale-task-trace.md` first. Cross-reference
   with `.agents/reports/2026-04-24_h1-h4-vi-results.md` for the canonical
   H1/H2a/H2b/H3/H4 table at n=9/10.
2. **If you are analyzing new data**: read
   `.agents/guides/01_project/agent-playbook.md` for the session workflow,
   then `.agents/guides/03_implementation/data-layout.md` for per-patient
   quirks (Pat_03 outlier, Pat_10 channel drop, Pat_13/14 phase-gap notes).
3. **If you are adding a method**: read
   `.agents/guides/02_methods/imcoh-guide.md` for the FC metric,
   `.agents/guides/02_methods/probe-bias-guide.md` for the volume-conduction
   bias that motivated the MSC → ImCoh switch, and
   `.agents/reports/2026-04-24_pipeline-status.md` for the era index flagging which
   artefacts are current vs superseded.

## Quick references

- `.agents/guides/INDEX.md` — full guide index
- `.agents/guides/03_implementation/cli-reference.md` — `lrg-eegfc` CLI
- `.agents/guides/03_implementation/function-map.md` — function lookup
- `.agents/guides/03_implementation/caching-guide.md` — cache layout
- `.agents/guides/03_implementation/figure-patterns.md` — plot templates

## Notebook header (required)

```python
from lrg_eegfc.notebook import *
move_to_root(pathname="lrgeegfc")
```

## Key invariants

- Cache-first: never recompute in visualization notebooks/scripts.
- Cache naming is parameter-sensitive (MSC uses `sparsify-*` and `nperseg-*`;
  ImCoh uses `nperseg-*` only — no surrogates).
- `fs` defaults to 2048 Hz; Pat_03 is the only override (1024 Hz) via
  `FS_OVERRIDES` in `src/lrg_eegfc/config/const.py`.
- Cross-phase analyses currently exclude **Pat_14** (missing `task_test`).
  Pat_06 is fully 4-phase as of 2026-04-22 — the old "exclude Pat_06" rule
  no longer applies. Pat_13 `rest_pre` was vendor-replaced 2026-04-23 and is
  now usable.
- Pat_10 is canonically 113-channel in every phase via load-time drop
  (`PATIENT_CHANNEL_DROP` in `config/const.py`). Raw files unmodified. See
  `memory/pat10_channel_mask.md` for the identification procedure.
