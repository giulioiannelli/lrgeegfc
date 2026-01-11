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
- Data inventory run completed and summarized in:
  - `.agents/plans/active/2026-01-10_data_inventory.md`
- Issues detected in Pat_05/06/07/08 (missing `fs`, missing phases, missing labels).
- Pat_03 has likely transposed `taskLearn`/`taskTest` shapes.
- Decisions: default missing `fs` to 2048 Hz; auto-transpose when channel
  dimension mismatches other phases; exclude patients missing phases.

## Tasks
1) Run data inventory
- Completed; see report file above.

2) Sampling rate checks
- Pending decision: which default `fs` to use for Pat_05/06/07/08.
- Decision: use 2048 Hz default for missing `fs` (no per-patient overrides).

3) Output report
- Markdown summary is in place.
- CSV report not found in repo; regenerate if needed for tooling.

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
- Handling plan for missing phases / missing fs / transposed data is agreed.
  - Missing phases: exclude Pat_06 and Pat_07 from cross-phase analysis runs.
