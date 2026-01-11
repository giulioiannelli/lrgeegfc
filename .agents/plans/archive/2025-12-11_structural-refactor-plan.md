# Structural Refactor Plan — 2025-12-11
Goal: polish structure and make constraints obvious to future agents (documentation + code hygiene), without breaking the existing cached workflows or visualization flow.

## Pain Points Observed
- Dual constant sources (`constants.py` vs `config/const.py`) and wildcard exports in `lrg_eegfc/__init__.py`/`utils/__init__.py` make it unclear which symbols are canonical.
- `lrg_eegfc/cli.py` still references a non-existent `io` module; some docs (`docs/overview.md`) reference modules that no longer exist.
- Visualization and workflow patterns are well-defined but not summarized centrally; agents rely on scattered plan files for invariants.
- Minimal safety nets: no smoke tests to catch missing imports/caches, and CLI defaults aren’t exercised automatically.

## Plan (ordered)
1) **Choose a single constants source**
   - Designate `config/const.py` as canonical; make `constants.py` a thin re-export (or alias) with deprecation note.
   - Update imports in workflows/visuals/CLIs to avoid mixing the two.
   - Document the chosen source in README + docs/overview + notebook helpers.

2) **Tighten package exports**
   - Replace `from .utils import *` / `from .config.const import *` in `lrg_eegfc/__init__.py` with explicit exports for the public API.
   - Mirror the explicit list in `utils/__init__.py` and `visuals/__init__.py`; add a short “public API” section in README.
   - Verify notebooks that rely on `from lrg_eegfc.notebook import *` still import the same symbols (adjust notebook.py re-exports if needed).

3) **Fix CLI wiring + align with robust loaders**
   - Update `lrg_eegfc/cli.py` to import loaders from `utils.datamanag.patient`/`patient_robust`; ensure entrypoint works without missing modules.
   - Normalize CLI defaults (dataset roots, cache roots, band lists) to match `compute_*` scripts; share helper functions where possible.
   - Add error messages that point users to cache-generation commands rather than recomputing inside visualization.

4) **Refresh architecture docs**
   - Rewrite `docs/overview.md` to match the current module layout (workflows, utils subpackages, visuals) and dataflow.
   - Add a short “structural invariants” section (caching, dataset layout, visualization rules) and link to `.agents/guides`.
   - Update `docs/developer-guide.md` / README developer section with the new constants decision and export policy.

5) **Add lightweight safety nets**
   - Introduce smoke tests (or a `make check-imports` script) that import CLIs, run a tiny synthetic workflow (corr/msc/lrg) without disk I/O, and assert cache path creation.
   - Consider a `python -m lrg_eegfc.workflows_smoke` script for quick validation in agent sessions.

6) **Codify visualization/cache conventions**
   - Centralize figure/cache naming rules in a single module or doc snippet used by `src/visualize_*.py` scripts.
   - Ensure all visualization CLIs share the same argument set for cache roots, cleaned vs raw selection, and output directories.

## Outputs
- Updated constants path + exports, fixed CLI import, refreshed docs, and a small smoke-test harness.
- `.agents/guides` (added) serve as the quick-start structural guide; keep them in sync when structure changes.

