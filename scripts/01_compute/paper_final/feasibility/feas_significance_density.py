#!/usr/bin/env python3
"""FEASIBILITY PROBE (2026-09-09) -- NOT A TEST.
Question for the definition page: if edges are kept only where |ImCoh|_band exceeds the no-coupling floor from
per-channel phase-randomised surrogates (autospectra preserved, cross-channel phase destroyed), what edge density
does each band keep on a whole resting phase? If ~everything is significant the protocol does not sparsify and a
different criterion is needed; if the density is band-dependent that becomes a data-driven variable of the paper.
One patient (arg), rest_pre, NSURR univariate FT-phase-randomised surrogates, same Welch recipe as the cache.
Reports per band: observed |ImCoh| quantiles, pooled null 95/99/99.9th percentiles, fraction of pairs above them,
and the fraction above the per-pair null max (NSURR values per pair, coarse), all with same-shaft pairs excluded.
Run: PYTHONPATH=src <lapbrain python> -u <this file> Pat_02
"""
import sys, time, numpy as np, pandas as pd
from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BANDS_NAMES, FS_OVERRIDES, DEFAULT_SAMPLE_RATE, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.fc.coherence.windowed import segment_ffts, imcoh_abs_cube, band_abs_average
from lrg_eegfc.utils.io.patient import load_timeseries, load_channel_labels
from lrg_eegfc.utils.probe import build_probe_mask
from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env(); OUT = ROOT / "data" / "paper_final" / "feasibility"
pat = sys.argv[1] if len(sys.argv) > 1 else "Pat_02"; NSURR = 10; rng = np.random.default_rng(1)
edges = [(BRAIN_BANDS[b][0], BRAIN_BANDS[b][1]) if not isinstance(BRAIN_BANDS[b], dict) else (BRAIN_BANDS[b]["low"], BRAIN_BANDS[b]["high"]) for b in BRAIN_BANDS_NAMES]
x = np.asarray(load_timeseries(pat, "rest_pre", SEEG_DATAPATH), float); N, T = x.shape
fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE); nperseg = nperseg_for_fs(fs); fmax = max(e[1] for e in edges) + 5.0
iu = np.triu_indices(N, 1); keep = ~build_probe_mask(list(load_channel_labels(pat)))[iu]

def bands_of(sig):
    freqs, ff = segment_ffts(sig, fs, nperseg, fmax_keep=fmax); cube = imcoh_abs_cube(ff)
    return np.stack([band_abs_average(cube, freqs, lo, hi)[iu] for lo, hi in edges])      # (n_band, n_pairs)

t0 = time.time(); obs = bands_of(x); print(f"{pat} N={N} T={T/fs:.0f}s obs pass {time.time()-t0:.0f}s", flush=True)
X = np.fft.rfft(x, axis=1); amp = np.abs(X); null = np.zeros((NSURR, len(edges), iu[0].size), np.float32)
for r in range(NSURR):
    ph = rng.uniform(0, 2 * np.pi, X.shape); ph[:, 0] = 0
    if T % 2 == 0: ph[:, -1] = 0
    xs = np.fft.irfft(amp * np.exp(1j * ph), n=T, axis=1); null[r] = bands_of(xs)
    print(f"  surrogate {r+1}/{NSURR} {time.time()-t0:.0f}s", flush=True)
rows = []
for bi, b in enumerate(BRAIN_BANDS_NAMES):
    o = obs[bi][keep]; nu = null[:, bi][:, keep]; pooled = nu.ravel(); pmax = nu.max(0)
    q = np.percentile(pooled, [50, 95, 99, 99.9])
    rows.append(dict(patient=pat, band=b, n_pairs=int(keep.sum()), obs_q10=np.percentile(o, 10), obs_med=np.median(o), obs_q90=np.percentile(o, 90),
                     null_med=q[0], null_p95=q[1], null_p99=q[2], null_p999=q[3],
                     frac_above_p95=float((o > q[1]).mean()), frac_above_p99=float((o > q[2]).mean()), frac_above_p999=float((o > q[3]).mean()),
                     frac_above_pairmax=float((o > pmax).mean()), z_med=float(np.median((o - nu.mean(0)) / (nu.std(0) + 1e-9)))))
df = pd.DataFrame(rows); df.to_csv(OUT / f"significance_density_{pat}.csv", index=False)
pd.set_option("display.width", 220); print(df.to_string(index=False, float_format=lambda v: f"{v:.3f}")); print(f"wall {time.time()-t0:.0f}s")
