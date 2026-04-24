---
name: msc-method-guide
type: guide
era: MSC
status: superseded
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

# MSC Method Guide

> **WARNING: MSC is contaminated by volume conduction.**
>
> Magnitude-squared coherence captures zero-lag volume conduction and
> common-reference artifacts, producing same-probe coherence 2-8x higher
> than cross-probe coherence. This bias dominates network structure and
> confounds all downstream analyses (LRG hierarchies, community detection,
> metastability).
>
> **For any analysis where same-probe bias matters, use ImCoh instead.**
> ImCoh (imaginary part of coherency) is immune to volume conduction by
> physics and is available via the same `compute_msc_welch()` function
> with `metric="imcoh"`.
>
> See `.agents/guides/IMCOH_GUIDE.md` for full details.
> See `.agents/guides/PROBE_BIAS_GUIDE.md` for the bias quantification.

---

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
- `scripts/py/compute_msc_matrices.py` (CLI defaults).
- `src/lrg_eegfc/config/const.py` (`DEFAULT_N_SURROGATES=200`).

## TODO
- Run surrogate scaling test and update `n_surrogates` target in
  `.agents/plans/active/2026-01-10_fc_params.md`.
