---
name: lrg-params
type: plan
era: IMCOH_SQ
status: superseded
created: 2026-01-10
updated: 2026-04-24
pointers: []
---

# LRG Parameters (2026-01-10)

This file records the current LRG defaults and where they are defined.
These values should be confirmed against the FIGMNTGN notebooks.

## Current defaults (code)
- `compute_lrg_analysis`:
  - `entropy_steps=400`
  - `entropy_t1=-3.0`
  - `entropy_t2=5.0`
  - Dendrogram cut: `compute_optimal_threshold(linkage_matrix)`
- Giant component extraction is always applied before LRG entropy.
- Cached outputs stored under `data/lrg_cache/Pat_XX/`.

## Dev vs full run (proposed)
- Dev: keep defaults but run a single patient/band/phase.
- Full: keep defaults unless notebooks indicate a different tau range or step count.

## Notebook references (to confirm)
- `ipynb/90_archive/.old/2025-12-11/FIGMNTGN01.ipynb`
- `ipynb/90_archive/.old/2025-12-11/FIGMNTGN02.ipynb`
- `ipynb/90_archive/.old/2025-12-11/FIGMNTGN03.ipynb`
- `ipynb/90_archive/.old/2025-12-11/FIGMNTGN04.ipynb`

## Open checks
- Confirm if any notebooks override `entropy_steps` or tau range.
- Confirm if any notebooks use a non-default dendrogram cut.
