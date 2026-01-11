# MSC Method Guide

This guide captures the current MSC (magnitude-squared coherence) methodology,
with defaults and rationale aligned to the refactored codebase.

## Core idea
MSC uses Welch's method, so variance depends on the number of segments, not on
Marchenko–Pastur noise. Correlation matrices may require MP cleaning; MSC does not.

## Welch settings
- `nperseg`: window length (samples) for spectral estimation.
- `noverlap`: overlap between segments; default `nperseg // 2`.
- Frequency resolution: `df = fs / nperseg`.
  - With `fs=2048` and `nperseg=1024`, `df ≈ 2 Hz`.

## Default policy (current)
- Dev runs (fast):
  - `sparsify=none`, `n_surrogates=0`, `nperseg=512`.
- Full runs (final):
  - `sparsify=soft`, `n_surrogates=200` (target), `nperseg=1024`.
  - Increase `n_surrogates` to 500 only if scaling test shows saturation.

## Surrogates (validation)
- Use circular-shift surrogates to build a null distribution per band.
- Final time-window analysis uses **per-window surrogates** (no reuse across windows).
- Dev runs can use 0–50 surrogates to keep runtime under 60 seconds.

## Time-window analysis
- Window length rule: at least 10 cycles of the lowest band frequency, and a
  minimum of 10 seconds.
- Window overlap: 25% between time windows (separate from Welch `noverlap`).
- Rationale: ensure enough samples (T >> N) to reduce MSC noise.

## Cache naming (important)
- Dense MSC: `{band}_{phase}_msc_sparsify-none_nperseg-{nperseg}.npy`
- Validated MSC: `{band}_{phase}_msc_sparsify-soft_nsurr-{n_surrogates}_nperseg-{nperseg}.npy`

## Where defaults live
- `src/lrg_eegfc/workflow/msc.py` (current defaults `nperseg=1024`).
- `src/compute_msc_matrices.py` (CLI defaults).
- `src/lrg_eegfc/config/const.py` (`DEFAULT_N_SURROGATES=200`).

## TODO
- Run surrogate scaling test and update `n_surrogates` target in
  `.agents/plans/active/2026-01-10_fc_params.md`.
