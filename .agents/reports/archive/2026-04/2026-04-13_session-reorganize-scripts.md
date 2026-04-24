---
name: session-zesty-roaming-lecun
type: report
era: IMCOH_SQ
status: dead
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

# Session: reorganize-data-scripts-config

**Date:** 2026-03-31 to 2026-04-13
**Resume name:** `reorganize-data-scripts-config`

## Scope

Full reorganization of the lrgeegfc repository: data folder, scripts, agent guides, and code paths.

## What was done

Reorganized `data/` from 31 flat directories down to 4 clean categories (`raw/`, `cache/`, `reports/`, `outputs/`). Created `config/paths.py` as the single source of truth for all data paths, replacing ~120 hardcoded `Path("data/...")` literals across the library. Created `utils/scripting.py` with shared helpers (`setup_script_env`, `iter_patient_band_phase`, `save_figure`, `load_all_lrg`). Added new constants to `const.py` (`PATIENTS_4PHASE`, `ALL_PHASE_PAIRS`, `classify_pair`). Moved 132 scripts from flat `scripts/py/` into 9 numbered directories (`01_compute` through `09_surrogate`) and refactored all of them to use centralized paths and helpers. Reorganized `.agents/guides/` into numbered subdirectories with an INDEX.md. Added ImCoh cache paths throughout.
