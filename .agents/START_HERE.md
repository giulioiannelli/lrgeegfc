# Start Here

This file is the entry point for agents. It points to the active plans and the
current state of the refactor.

## Current focus
- Use `.agents/plans/INDEX.md` for the canonical list of plans.
- Active refactor + analysis plans live in `.agents/plans/active/`.
- Completed refactor decisions are in `.agents/plans/developed/`.

## Where to begin
1) Read `.agents/plans/developed/2026-01-10_stage-06_refactor-packaging.md` (refactor status + layout).
2) Read `.agents/plans/active/2026-01-10_stage-01_data-qc.md` (data inventory + QC).
3) Read `.agents/plans/active/2026-01-10_stage-02_fc-msc-corr.md` (MSC/corr plan).
4) Read `.agents/plans/active/2026-01-10_stage-00_notebook-consolidation.md` (notebook cleanup).
5) Read `.agents/plans/active/2026-01-10_stage-03_lrg.md` and `2026-01-10_stage-04_reorganization-metrics.md`.
6) Read `.agents/plans/active/2026-01-10_stage-05_figures-overleaf.md` and `2026-01-10_stage-05u_spatial-embedding.md`.
7) Read `.agents/guides/MSC_METHOD_GUIDE.md`, `.agents/guides/CACHING_GUIDE.md`,
   and `.agents/guides/TIME_WINDOW_GUIDE.md`.

## Quick references (for fast lookup)
- `.agents/guides/CLI_REFERENCE.md` - CLI command reference (32 subcommands, examples, workflows)
- `.agents/guides/FUNCTION_MAP.md` - Complete function reference (130+ functions)
- `.agents/guides/FIGURE_PATTERNS.md` - Figure templates and patterns
- `.agents/guides/AGENT_TASKS.md` - Common autonomous task procedures

## Notebook header (required)
```python
%matplotlib inline
from lrgsglib.config.funcs import move_to_rootf
move_to_rootf(pathname="lrgeegfc")
from lrg_eegfc.notebook import *
```

## Key invariants
- Cache-first: never recompute in visualization notebooks/scripts.
- Use `lrg_eegfc.notebook` to avoid long imports and to keep CWD correct.
- Cache naming is parameter-sensitive (MSC uses `sparsify-*` and `nperseg-*`).
- Default missing `fs` to 2048 Hz; auto-transpose phases if channel count differs.
- Exclude patients missing phases (currently Pat_06, Pat_07) from cross-phase runs.
