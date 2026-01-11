# Data Inventory Report (2026-01-10)

This report summarizes the output of `src/inspect_patient_data.py` for the
current `data/stereoeeg_patients` tree.

## Summary
- Patients inspected: 6
- Patients with issues: 4 (Pat_05, Pat_06, Pat_07, Pat_08)
- Total issues flagged: 19

## Sampling rates
- `2048 Hz`: Pat_02
- `1024 Hz`: Pat_03
- Missing `fs`: Pat_05, Pat_06, Pat_07, Pat_08

## Phase availability
- Missing phases:
  - Pat_06: `taskLearn`, `taskTest`
- Corrupted/missing data variable:
  - Pat_07: `taskTest` (no `Data` variable)

## Data shape anomalies
- Pat_03:
  - `taskLearn` shape `(856595, 122)` (likely transposed)
  - `taskTest` shape `(1235200, 122)` (likely transposed)
  - Other phases are `(122, T)`, so verify axis ordering in loaders.

## Metadata availability
- `channel_labels.csv` present: Pat_02, Pat_03, Pat_05
- `channel_labels.csv` missing: Pat_06, Pat_07, Pat_08
- `ChannelNames.mat` present: Pat_02 only
- Implant coordinate files: none detected (all patients)

## Follow-up actions
- Decide how to handle missing `fs` for Pat_05/06/07/08 (default value or metadata fix).
- Decide whether to transpose Pat_03 `taskLearn`/`taskTest` in the loader.
- Confirm how to handle missing phases (skip vs exclude patient in cross-phase analyses).
- Collect or map channel labels / implant coordinates for spatial plots.
