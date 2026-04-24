#!/usr/bin/env python3
"""Pre-compute expensive data for FC section figures.

Caches:
  1. Full coherence tensor (freqs, Coh) for fig:msc_example
  2. Surrogate null distribution for fig:surrogates

Run:
    python scripts/precompute_fc_figures.py
"""
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
from lrg_eegfc.config.const import (
    BRAIN_BANDS, DEFAULT_SAMPLE_RATE,
)
from lrg_eegfc.config.paths import FC_FIG_CACHE, SEEG_DATAPATH
from lrg_eegfc.utils.io import load_timeseries
from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch, band_average_msc

# ── Config ────────────────────────────────────────────────────────────
PATIENT = "Pat_02"
PHASE = "rest_pre"
NPERSEG = 4096
N_SURROGATES = 200
BAND_FIG = "beta"
CACHE_DIR = FC_FIG_CACHE
CACHE_DIR.mkdir(parents=True, exist_ok=True)

fs = DEFAULT_SAMPLE_RATE

# ── 1. Full coherence tensor ─────────────────────────────────────────
coh_path = CACHE_DIR / f"{PATIENT}_{PHASE}_coh_nperseg-{NPERSEG}.npz"
if coh_path.exists():
    print(f"[1/2] Coherence tensor already cached: {coh_path}")
    data = np.load(coh_path)
    freqs, Coh = data["freqs"], data["Coh"]
else:
    print(f"[1/2] Computing full coherence tensor for {PATIENT} {PHASE}...")
    X = load_timeseries(PATIENT, PHASE, SEEG_DATAPATH)
    print(f"  Timeseries shape: {X.shape}")
    freqs, Coh = compute_msc_welch(X, fs, nperseg=NPERSEG)
    print(f"  Coherence shape: {Coh.shape}")
    np.savez_compressed(coh_path, freqs=freqs, Coh=Coh)
    print(f"  Saved: {coh_path} ({coh_path.stat().st_size / 1e6:.1f} MB)")

# Pick representative channel pair (95th percentile of beta coherence)
W_bands = band_average_msc(Coh, freqs, BRAIN_BANDS)
W_beta = W_bands["beta"].copy()
np.fill_diagonal(W_beta, 0)
triu_idx = np.triu_indices_from(W_beta, k=1)
beta_vals = W_beta[triu_idx]
target = np.percentile(beta_vals, 95)
best_idx = np.argmin(np.abs(beta_vals - target))
ch_i, ch_j = int(triu_idx[0][best_idx]), int(triu_idx[1][best_idx])
observed = W_bands[BAND_FIG][ch_i, ch_j]
print(f"  Representative pair: ({ch_i}, {ch_j}), observed {BAND_FIG} MSC = {observed:.4f}")

# ── 2. Surrogate null distribution ───────────────────────────────────
null_path = CACHE_DIR / f"{PATIENT}_{PHASE}_{BAND_FIG}_null_ch{ch_i}-{ch_j}_nsur-{N_SURROGATES}.npz"
if null_path.exists():
    print(f"[2/2] Null distribution already cached: {null_path}")
else:
    print(f"[2/2] Computing {N_SURROGATES} surrogates for pair ({ch_i}, {ch_j})...")
    X = load_timeseries(PATIENT, PHASE, SEEG_DATAPATH)
    fmin, fmax = BRAIN_BANDS[BAND_FIG]
    rng = np.random.default_rng(42)
    null_values = np.empty(N_SURROGATES)

    for s in range(N_SURROGATES):
        X_sur = np.empty_like(X)
        for ch in range(X.shape[0]):
            shift = rng.integers(0, X.shape[1])
            X_sur[ch] = np.roll(X[ch], shift)

        X_pair = X_sur[[ch_i, ch_j], :]
        f_sur, Coh_sur = compute_msc_welch(X_pair, fs, nperseg=NPERSEG)
        band_mask = (f_sur >= fmin) & (f_sur <= fmax)
        null_values[s] = Coh_sur[0, 1, band_mask].mean()

        if (s + 1) % 25 == 0:
            print(f"  {s+1}/{N_SURROGATES}")

    p_val = (1 + np.sum(null_values >= observed)) / (N_SURROGATES + 1)
    np.savez_compressed(
        null_path,
        null_values=null_values, observed=observed,
        ch_i=ch_i, ch_j=ch_j, p_val=p_val,
        band=BAND_FIG,
    )
    print(f"  p-value: {p_val:.4f}")
    print(f"  Saved: {null_path}")

print("\nDone! All precomputed data cached in:", CACHE_DIR)
