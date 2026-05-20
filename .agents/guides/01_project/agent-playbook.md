---
name: agent-playbook
type: guide
era: CROSS_ERA
status: current
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

# Agent Playbook: How to Work Here

## Session Kick-off
- Read `.agents/START_HERE.md` and the latest `ToDo.md` in `.agents/plans/active/`.
- Skim `.agents/plans/active/*ANALYSIS_COMMANDS*.md` to see the canonical pipeline and expected cache layout.
- Check `.agents/plans/active/*DATA_STATUS*.md` before running analyses; missing fs/phase files are documented there.

## Do / Don't
- **Do** reuse caches; **don't** recompute inside visualization/CLI helpers. If data is missing, fail with a clear message and point to the cache path.
- **Do** use `config/const.py` for bands/phases; **don't** introduce new constants modules.
- **Do** route data loading through `utils/io/patient_robust.py` when you need fs extraction; **don't** access `.mat` files ad hoc.
- **Do** keep outputs under `data/{corr_cache,msc_cache,lrg_cache,figures}/` and `results/`; **don't** write into `ipynb/` or repo root.
- **Do** add docstrings/type hints and keep functions pure; **don't** add hidden state or global configuration.
- **Do** update exports (`__init__.py` or `visuals/__init__.py`) when adding modules; **don't** leave orphan utilities.

## Implementation Patterns to Mirror
- Workflows return dataclasses capturing metadata (`CorrResult`, `MSCResult`, `LRGResult`, `CleanedCorrResult`)—extend the pattern when adding new workflows.
- CLI scripts follow the argparse + cache-first pattern (see `scripts/py/compute_corr_matrices.py`, `scripts/py/compute_msc_matrices.py`); keep user-facing defaults and verbose logging consistent.
- Visual modules accept cached matrices/graphs and plain `Path` destinations; keep plotting code notebook-faithful (see `visuals/correlation.py` for style).
- Notebook ergonomics go through `lrg_eegfc.notebook` + `setup_notebook()`; when adding helpers, make them available there.
- Standard notebook header (avoid empty patient list): `move_to_rootf(pathname="lrgeegfc")` then `from lrg_eegfc.notebook import *`.

## Patient Inclusion / Exclusion

**Cohort (locked 2026-04-25, n = 10):** Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15.

| Patient | fs (Hz) | Phases available | Status |
|---------|---------|-----------------|--------|
| Pat_02 | 2048 | rest_pre, task_learn, task_test, rest_post | included |
| Pat_03 | **1024** | rest_pre, task_learn, task_test, rest_post | included |
| Pat_05 | 2048 | rest_pre, task_learn, task_test, rest_post | included |
| Pat_06 | 2048 | rest_pre, task_learn, task_test, rest_post | included (data-layout §6: phases backfilled 2026-04-22) |
| Pat_07 | 2048 | rest_pre, task_learn, task_test, rest_post | included |
| Pat_08 | 2048 | rest_pre, task_learn, task_test, rest_post | included |
| Pat_10 | 2048 | rest_pre, task_learn, task_test, rest_post | included (task rows [53, 54, 55] dropped at load) |
| Pat_13 | 2048 | rest_pre, task_learn, task_test, rest_post | included (rest_pre vendor-replaced 2026-04-23) |
| Pat_14 | 2048 | rest_pre, task_learn, task_test, rest_post | included (task_test vendor-replaced 2026-04-25) |
| Pat_15 | 2048 | rest_pre, task_learn, task_test, rest_post | included |

### Pat_03 sampling-rate handling (NOT exclusion)

Pat_03 is acquired at 1024 Hz (every other patient at 2048 Hz). The
sampling-rate difference is absorbed **at the config layer only** —
`nperseg_for_fs(fs)` returns `nperseg = 2048` for Pat_03 versus
`nperseg = 4096` for the 2048 Hz patients, so the Welch-segment
duration (≈ 2 seconds) is the same. `FS_OVERRIDES` in
`config/const.py` pins Pat_03's `fs = 1024`.

**At the analysis level, Pat_03 is identical to every other cohort
member.** There is no separate patient list, no dropout sensitivity
test, no figure marker distinguishing Pat_03, no "report values
separately" rule. The previous "EXCLUDED / negative control / 3×
denser FC" framing was retired on 2026-05-18; the sampling-rate
difference does not propagate beyond the config layer.

### nperseg and sampling rate

Always use `nperseg_for_fs(fs)` from `config.const` to compute nperseg:
```python
from lrg_eegfc.config.const import nperseg_for_fs
nperseg = nperseg_for_fs(fs)  # 4096 at 2048 Hz, 2048 at 1024 Hz
```
This ensures 2-second Welch segments (Δf = 0.5 Hz) for all patients.

## Known Rough Edges (flag for refactors)
- Keep the workflow and utils boundaries clean; avoid pushing visualization logic into `workflow/`.

## When Adding Features
- Decide the layer: `utils/*` (low-level), `workflow/*` (cached orchestration), `visuals/*` (plots), or `scripts/py/*.py` (CLI wrapper). Keep concerns separated.
- Add small, well-named helper functions instead of inlining notebook code; if extracted from a notebook, note the source cell in a comment if it affects behaviour.
- Update relevant guides/plan files with new commands or expectations.
- Run spot checks where possible (`pytest` if present, or a small CLI invocation) and record what was run in your notes/commit message.
