# Refactor Decisions (2026-01-10)

## Canonical modules
- Constants: keep `src/lrg_eegfc/config/const.py`; remove `src/lrg_eegfc/constants.py`.
- Shared imports: replace `src/lrg_eegfc/shared.py` with `src/lrg_eegfc/utils/common.py`.
- Workflows: move to `src/lrg_eegfc/workflow/` and keep top-level stubs.
- Utilities: move under `src/lrg_eegfc/utils/{io,fc,lrg,metrics,pipelines}`.
- Visuals: keep under `src/lrg_eegfc/visuals/`; move plotting helpers to `visuals/plotting.py`.
- LRG visuals: keep `visuals/lrg.py` as canonical; `visuals/lrg_revised.py` and `visuals/lrg_backup.py` are wrappers.

## Compatibility matrix (old -> new)
- `lrg_eegfc.constants` -> `lrg_eegfc.config.const`
- `lrg_eegfc.shared` -> `lrg_eegfc.utils.common`
- `lrg_eegfc.utils.datamanag.*` -> `lrg_eegfc.utils.io.*`
- `lrg_eegfc.utils.corrmat.*` -> `lrg_eegfc.utils.fc.corr.*`
- `lrg_eegfc.utils.coherence.*` -> `lrg_eegfc.utils.fc.msc.*`
- `lrg_eegfc.utils.clustering.*` -> `lrg_eegfc.utils.lrg.*`
- `lrg_eegfc.utils.distances.*` -> `lrg_eegfc.utils.metrics.*`
- `lrg_eegfc.compare` -> `lrg_eegfc.utils.metrics.compare` (stub kept)
- `lrg_eegfc.plotting` -> `lrg_eegfc.visuals.plotting` (stub kept)
- `lrg_eegfc.batch_compute` -> `lrg_eegfc.utils.pipelines.batch_compute` (stub kept)
- `lrg_eegfc.workflow_*` -> `lrg_eegfc.workflow.*` (stubs kept)
- `lrg_eegfc.io` -> `lrg_eegfc.utils.io` (stub kept)

## Notes
- `PATIENTS_LIST` is now lazy-loaded to avoid filesystem scans on import.
- `BRAIN_BAND_LABELS` is retained as an alias for `BRAIN_BAND_TEX_DICT`.
