# Plan (2026-01-10): Stage 00 - Notebook consolidation

## Goals
- Convert the notebook set into a clean, ordered catalog that maps to the
  pipeline modules.
- Preserve legacy work while extracting reusable logic.
- Prevent duplicated analyses from diverging.

## Inputs
- `ipynb/INDEX.md` (current catalog)
- `ipynb/90_archive/.old/` (legacy notebooks to review)

## Current status
- Notebook folders created (`ipynb/00_intake` ... `ipynb/05_figures`, `ipynb/90_archive`).
- Active notebooks moved into ordered folders.
- Legacy notebooks are still under `ipynb/90_archive/.old/` until extraction is complete.
- `ipynb/INDEX.md` exists and lists all notebooks with keep/merge/drop tags.
- Stage 04 + Stage 05 legacy notebooks were merged into new notebooks and moved to `.old_reviewed/`.
- Stage 02 legacy notebooks were merged into updated single-case notebooks and moved to `.old_reviewed/`.
- Stage S time-window notebooks were merged into updated notebooks and moved to `.old_reviewed/`.

## Tasks
1) Inventory and tagging
- Completed; see updated `ipynb/INDEX.md`.
2) Legacy cleanup
- Do not move notebooks to `.old_reviewed/` until their code is extracted.

2) Extraction to modules or new notebooks
- For each "keep" or "merge" notebook, extract reusable code into
  `src/lrg_eegfc/` modules or a new ordered notebook.
- Add a short note in the new notebook indicating which legacy notebook it
  replaced.
- Standard notebook header (to avoid empty patient list):
  - `move_to_rootf(pathname="lrgeegfc")`
  - `from lrg_eegfc.notebook import *`
- Target destination folders:
  - `ipynb/01_preprocessing/` (data checks, time windows)
  - `ipynb/02_fc_msc/` (single-patient FC construction)
  - `ipynb/03_lrg/` (LRG derivations, PSI)
  - `ipynb/04_reorganization/` (phase distance metrics)
  - `ipynb/05_figures/` (single-case figure builders)

3) Archive lifecycle
- Move fully extracted notebooks to `ipynb/90_archive/.old_reviewed/`.
- Only delete notebooks after all their functionality is migrated.

## Compute vs visualize
- Notebooks are for single-case, cache-driven visualization or quick checks.
- Heavy computation should move into `src/` scripts and write caches.
- Notebook cells should use the load-first pattern from `CACHING_GUIDE.md`.

## Deliverables
- Updated `ipynb/INDEX.md` with tags and migration notes.
- New or updated notebooks with ordered IDs and descriptive names.
- Reviewed legacy notebooks moved under `.old_reviewed`.

## Exit criteria
- Each active notebook maps to a specific pipeline stage or figure.
- Legacy notebooks are either in `.old_reviewed` or flagged for extraction.
