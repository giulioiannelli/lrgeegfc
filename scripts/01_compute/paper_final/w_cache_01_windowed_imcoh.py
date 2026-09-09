#!/usr/bin/env python3
"""Windowed |ImCoh| cache: non-overlapping WIN-second windows, all six bands, four phases, n = 10.
Not a test. It is the substrate for (a) the per-level ordered sham of Lane M1 (reliability-gradient
control), (b) the objective-3 windowed constructions (drift kernel, rest_post lag profile, early/late DiD).
Recipe = library windowed backend (segment FFTs -> per-bin |Im C| -> band average), which reproduces the
cached full-phase imcoh_abs in the full-window limit (audit_124). One FFT pass per window, six band averages.
Output: CACHE_ROOT/windowed_imcoh/{patient}_{phase}_w{WIN}s.npz with W (n_win, n_band, N, N) float32,
bands, t_start (s), fs, nperseg. Run: PYTHONPATH=src <lapbrain python> -u <this file> [--workers 8 --win 30]
"""
import argparse, time, numpy as np
from multiprocessing import Pool
from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS, BRAIN_BANDS_NAMES, FS_OVERRIDES, DEFAULT_SAMPLE_RATE, nperseg_for_fs
from lrg_eegfc.config.paths import CACHE_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.fc.coherence.windowed import segment_ffts, imcoh_abs_cube, band_abs_average
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.scripting import setup_script_env
setup_script_env()
PHASES = ("rest_pre", "task_learn", "task_test", "rest_post"); OUT = CACHE_ROOT / "windowed_imcoh"
WIN = 30

def band_edges(name):
    b = BRAIN_BANDS[name]; return (b[0], b[1]) if not isinstance(b, dict) else (b["low"], b["high"])

def per_job(job):
    pat, ph, win = job; t0 = time.time()
    out = OUT / f"{pat}_{ph}_w{win}s.npz"
    if out.exists(): return f"{pat}/{ph} cached"
    x = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), float); N, T = x.shape
    fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE); nperseg = nperseg_for_fs(fs); L = int(win * fs)
    n_win = T // L; edges = [band_edges(b) for b in BRAIN_BANDS_NAMES]; fmax = max(e[1] for e in edges) + 5.0
    W = np.zeros((n_win, len(edges), N, N), np.float32)
    for w in range(n_win):
        freqs, ff = segment_ffts(x[:, w * L:(w + 1) * L], fs, nperseg, fmax_keep=fmax)
        cube = imcoh_abs_cube(ff)
        for bi, (lo, hi) in enumerate(edges): W[w, bi] = band_abs_average(cube, freqs, lo, hi)
    np.savez(out, W=W, bands=np.array(BRAIN_BANDS_NAMES), t_start=np.arange(n_win) * win, fs=fs, nperseg=nperseg, N=N, T=T)
    return f"{pat}/{ph} N={N} T={T/fs:.0f}s n_win={n_win} {time.time()-t0:.0f}s"

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=8); ap.add_argument("--win", type=int, default=WIN); a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(p, ph, a.win) for p in PATIENTS_4PHASE for ph in PHASES]
    print(f"[wcache] {len(jobs)} jobs | win={a.win}s | {a.workers} workers -> {OUT}", flush=True); t0 = time.time()
    with Pool(a.workers) as pool:
        for i, msg in enumerate(pool.imap_unordered(per_job, jobs), 1):
            el = time.time() - t0; print(f"[{i}/{len(jobs)}] {msg} | elapsed {el/60:.1f} min ETA {el/i*(len(jobs)-i)/60:.1f} min", flush=True)
    print(f"wall {(time.time()-t0)/60:.1f} min", flush=True)
