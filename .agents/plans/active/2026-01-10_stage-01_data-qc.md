# Plan (2026-01-10): Stage 01 - Data inventory and QC

## Goals
- Verify dataset completeness and metadata for every patient/phase.
- Capture sampling rate availability and channel counts.
- Produce a machine-readable report that downstream steps can trust.

## Inputs
- `data/stereoeeg_patients/Pat_XX/{phase}.mat`
- Optional metadata CSVs (`Implant_pat_XX.csv`, `channel_labels.csv`)
- Existing helpers: `src/inspect_patient_data.py`,
  `utils/datamanag/patient_robust.py` (post-refactor: `lrg_eegfc.utils.io`)

## Current status
- `src/inspect_patient_data.py` exists but does not emit a dated report in `.agents/`.
- Robust loader is available in `utils/datamanag/patient_robust.py`.
- No consolidated report under `.agents/plans/active/` yet.

## Tasks
1) Run data inventory
- Create or update a script to list patients, phases, channel counts, and
  sampling rate extraction status.
- Log missing phases and missing metadata explicitly.
- Track which patients have implant coordinate files (Implant_pat_XX.csv).
- Implementation target: update `src/inspect_patient_data.py` to emit both
  Markdown and CSV reports.

2) Sampling rate checks
- Confirm whether `fs` is present in the `Parameters` struct or alternate keys.
- Record whether the loader needs to fall back to default `sample_rate`.

3) Output report
- Store a summary report under `.agents/plans/active/` with date.
- Capture warnings and per-patient anomalies.
- Expected outputs:
  - `.agents/plans/active/2026-01-10_data_inventory.md`
  - `.agents/plans/active/2026-01-10_data_inventory.csv`

## Compute vs visualize
- This stage is compute-only (report generation). No figures are produced here.
- Downstream visualization should read the report instead of re-scanning data.

## Deliverables
- Data report file (Markdown or CSV) in `.agents/plans/active/`.
- Clear list of missing files or fields that block pipelines.

## Exit criteria
- For each patient/phase: data path, channel count, fs status, and metadata
  availability are known.
- Patients with electrode coordinate files are identified for spatial plots.
