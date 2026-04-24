---
name: time-window-guide
type: guide
era: IMCOH_ABS
status: current
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

# Time-Window Guide (Stage S)

This guide defines the time-window analysis defaults and how to cache and
visualize windowed FC/LRG results.

## Why time windows
- Track network reorganization dynamics within each phase.
- Keep windows long enough that spectral estimates are stable (T >> N).
- Use the same protocol for correlation and MSC to compare behavior.

## Default windowing
- Window length (auto): `max(10 s, 10 cycles at lowest band frequency)`.
- Overlap: 25% (balance stability + independence).
- If you need strict alignment across bands, pass an explicit `--window-sec`.

## Stability notes (T >> N)
- Let `T` = window samples, `N` = channels.
- If `T/N < 10`, expect noisy FC estimates and MP overlap.
- Use longer windows or fewer channels when possible.

## MSC vs correlation in windows
- Correlation: bandpass per window + Pearson correlation.
- MSC: compute Welch coherence per window, then band-average.
- Final run uses **per-window surrogates** for MSC validation.
- Dev run uses `sparsify=none` and `filter_time` to keep runtime short.

## Caching layout
- Window FC caches live under:
  - `data/corr_cache_windows/Pat_XX/<phase>/<band>/<run_id>/`
  - `data/msc_cache_windows/Pat_XX/<phase>/<band>/<run_id>/`
- Window LRG caches live under:
  - `data/lrg_cache_windows/Pat_XX/<phase>/<band>/<run_id>/`
- Each run directory contains:
  - `run_meta.json` (params)
  - `windows.csv` (start/stop indices)
  - `win-0000.npy` matrices
  - `win-0000_lrg.npz` if LRG is computed
- Dev runs with `filter_time` go to `*_cache_windows_dev`.

## Compute vs visualize
- Compute: `scripts/py/compute_time_windows.py` (writes caches).
- Visualize: `scripts/py/visualize_time_windows.py` (reads caches only).
- Notebooks should load cached windows and avoid recomputation.

## Recommended dev run
```bash
python scripts/py/compute_time_windows.py \
  --patients Pat_02 --phases rest_post --bands alpha \
  --fc-method corr --window-sec 10 --overlap 0.25 \
  --filter-time 5000 --max-windows 5 --verbose
```

## Recommended final run (MSC validation)
```bash
python scripts/py/compute_time_windows.py \
  --patients Pat_02 --phases rest_post --bands alpha \
  --fc-method msc --window-sec 20 --overlap 0.25 \
  --sparsify soft --n-surrogates 200 --nperseg 1024 --lrg
```
