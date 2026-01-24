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
- CLI scripts follow the argparse + cache-first pattern (see `src/compute_corr_matrices.py`, `src/compute_msc_matrices.py`); keep user-facing defaults and verbose logging consistent.
- Visual modules accept cached matrices/graphs and plain `Path` destinations; keep plotting code notebook-faithful (see `visuals/correlation.py` for style).
- Notebook ergonomics go through `lrg_eegfc.notebook` + `setup_notebook()`; when adding helpers, make them available there.
- Standard notebook header (avoid empty patient list): `move_to_rootf(pathname="lrgeegfc")` then `from lrg_eegfc.notebook import *`.

## Known Rough Edges (flag for refactors)
- Keep the workflow and utils boundaries clean; avoid pushing visualization logic into `workflow/`.

## When Adding Features
- Decide the layer: `utils/*` (low-level), `workflow/*` (cached orchestration), `visuals/*` (plots), or `src/*.py` (CLI wrapper). Keep concerns separated.
- Add small, well-named helper functions instead of inlining notebook code; if extracted from a notebook, note the source cell in a comment if it affects behaviour.
- Update relevant guides/plan files with new commands or expectations.
- Run spot checks where possible (`pytest` if present, or a small CLI invocation) and record what was run in your notes/commit message.
