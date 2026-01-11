# Plan (2026-01-10): Stage 06 - Refactor and packaging

## Goals
- Clarify module boundaries without breaking existing imports.
- Ensure scripts are thin wrappers over reusable modules.

## Current status
- Refactor pass completed: workflows split into `workflow/`, utils into `utils/{io,fc,lrg,metrics,pipelines}`, visuals kept top-level.
- Compatibility stubs added for `workflow_*`, `plotting.py`, `compare.py`, `batch_compute.py`, and `io.py`.
- Constants unified under `config/const.py` (lazy patient listing).
- Docs and guides updated; sanity notebook validated cache reads.

## Tasks
1) Inventory + duplication audit
- Identify duplicate or overlapping modules and decide the canonical version:
  - `constants.py` vs `config/const.py` (keep `config/const.py`, remove `constants.py`).
  - `visuals/lrg.py`, `visuals/lrg_revised.py`, `visuals/lrg_backup.py` (merge or archive).
  - `compare.py` vs `utils/distances/comparison.py` (clarify compute vs metrics).
  - `shared.py` usage (remove; update any import sites).
- Record decisions in `2026-01-10_refactor_decisions.md`.
  - Include a compatibility matrix: old import -> new import.

2) Create subpackages (compatibility-first)
- Keep `config/` at top level.
- Keep `visuals/` at top level.
- Move domain utilities under `utils/`:
  - `utils/io/`, `utils/fc/`, `utils/lrg/`, `utils/metrics/`, `utils/pipelines/`
- Initial step: add `__init__.py` re-exports that point to existing code.
- Do not move files until all re-exports are in place.
- Create `workflow/` package (top-level) and split workflows into:
  - `workflow/corr.py`, `workflow/msc.py`, `workflow/lrg.py`, `workflow/cleaning.py`, `workflow/core.py`

3) Compatibility layer
- Keep existing imports working via re-exports in `__init__.py`.
- Preserve file names in `src/*.py` scripts; update their imports only.

4) Documentation updates
- Align `docs/overview.md` and `.agents/README.md` with the new layout.
- Update notebook references in `ipynb/INDEX.md` if module paths change.

## Compute vs visualize
- Keep compute logic in `utils/pipelines/` and workflows; keep visuals in
  top-level `visuals/`.
- Enforce load-first usage in visualization modules to avoid recompute.

## File-level migration map (first pass)
- `src/lrg_eegfc/utils/io/`
  - `__init__.py` re-exports from `utils/datamanag/{loaders,patient,patient_robust}.py`
  - Later move those files into `utils/io/` and update imports.
- `src/lrg_eegfc/utils/fc/`
  - `corr.py` wraps `workflow_corr` and `utils/corrmat/*`
  - `msc.py` wraps `workflow_msc` and `utils/coherence/*`
  - `cleaning.py` wraps `workflow_cleaning`
- `src/lrg_eegfc/utils/lrg/`
  - `analysis.py` wraps `workflow_lrg` (LRGResult, compute_lrg_analysis)
- `src/lrg_eegfc/utils/metrics/`
  - `reorganization.py` wraps `visuals/reorganization.py`
  - `compare.py` wraps `compare.py` and `utils/distances/*`
- `src/lrg_eegfc/utils/pipelines/`
  - New modules that `src/*.py` scripts call (thin orchestration)
- `src/lrg_eegfc/workflow/`
  - Split existing `workflow_*.py` files into submodules
  - Keep thin wrapper modules at top-level for compatibility

## Detailed moves (second pass)
- Move `utils/datamanag/*` -> `utils/io/` and update imports in:
  - `workflow/corr.py`, `workflow/msc.py`, `workflow/cleaning.py`, `cli.py`
- Move `utils/corrmat/*` -> `utils/fc/corr/` (or keep in `utils/fc/corr.py` as wrappers)
- Move `utils/coherence/*` -> `utils/fc/msc/` (or keep in `utils/fc/msc.py` as wrappers)
- Keep `visuals/*` at top level (no move); update imports to use `visuals/`.
- Move `compare.py` -> `utils/metrics/compare.py` and re-export from root
- Create `utils/metrics/reorganization.py` for non-plot computations currently in
  `visuals/reorganization.py` (split compute vs plot)
- Deprecate `lrg_revised.py` and `lrg_backup.py` after merging missing pieces
  into `visuals/lrg.py` (keep files as thin wrappers for compatibility)
- Move `batch_compute.py` into `utils/pipelines/` and keep a top-level stub.
- Keep `notebook.py` at top level; move `plotting.py` under `visuals/`.
- Remove `shared.py` and clean any `from lrg_eegfc.shared import *` usage.
- Normalize constants:
  - Remove `constants.py`.
  - Keep `config/const.py` and ensure it does not scan the filesystem on import;
    move `PATIENTS_LIST` to a helper function if needed.

## Compatibility checklist
- Update `lrg_eegfc/cli.py` imports (currently references missing `io` module).
- Ensure `lrg_eegfc/__init__.py` still exposes current public API.
- Maintain cache paths and filenames unchanged.
- Ensure top-level `workflow_*.py` modules remain importable (stubs).

## Deliverables
- New subpackage layout with compatibility shims.
- Updated documentation and scripts.
- Refactor decision log with duplicate resolution.
 - Compatibility matrix document under `.agents/plans/active/`.

## Exit criteria
- All existing pipelines run without changing external commands.
- Documentation reflects the new structure.

## Status update (2026-01-10)
### Completed
- Moved workflows to `src/lrg_eegfc/workflow/` and added stubs at top level.
- Moved utils into `src/lrg_eegfc/utils/{io,fc,lrg,metrics,pipelines}`.
- Moved plotting helpers to `src/lrg_eegfc/visuals/plotting.py`.
- Removed `constants.py` and `shared.py`; added `utils/common.py` for shared imports.
- Updated imports across `src/` scripts and tests.
- Updated `docs/overview.md`, `.agents/README.md`, and agent guides.
- Added sanity notebook `ipynb/00_intake/01_sanity_imports.ipynb`.

### Key learnings to carry forward
- Notebooks must call `move_to_rootf(pathname="lrg_eegfc")` or list_patients will
  return empty results due to relative paths.
- Cache filenames encode parameters (e.g., `msc` uses `sparsify-*` and `nperseg-*`);
  notebook loaders should match those defaults or explicitly fall back.
- MSC and correlation caches are distinct; avoid recomputation in notebooks.

### Open follow-ups (optional)
- Decide whether to keep absolute imports (`lrg_eegfc.*`) or switch to relative
  imports inside the package for style consistency.
## Sequential workflow (do in order)
1) Produce `2026-01-10_refactor_decisions.md` with duplicate resolutions.
2) Add re-export shims in the new `utils/*` subpackages (no moves yet).
3) Update imports in `src/*.py` scripts to use the shims.
4) Move the lowest-risk modules first (`utils/datamanag` -> `utils/io`).
5) Update workflow imports (`workflow_corr`, `workflow_msc`, `workflow_cleaning`, `cli`).
6) Migrate FC utilities (`utils/corrmat`, `utils/coherence`) to `utils/fc/`.
7) Migrate metrics (`compare.py`, parts of `visuals/reorganization.py`) to `utils/metrics/`.
8) Merge/retire `visuals/lrg_revised.py` + `visuals/lrg_backup.py` into `visuals/lrg.py`.
9) Update docs + `ipynb/INDEX.md` and run a small cache read to validate imports.

## Stop conditions
- If any refactor step breaks cached loads or CLI scripts, pause and fix before
  moving to the next step.
