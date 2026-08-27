#!/usr/bin/env python3
"""W0-B pre-flight: validate (a) the observed rho_sym reproduction and (b) the
circular-shift <-> frequency-domain phase-rotation identity.

Nothing here is a null. This script exists so that every rung of the ladder is
attached to the SAME observed statistic that script 13
(``sparsified_arc/13_matched_strength_mst020.py``) gates, and so that the
cheap frequency-domain implementation of N1 is known to be a faithful stand-in
for a genuine time-domain circular shift rather than an assumed one.

Checks
------
A. Observed rho_sym(s) recomputed from the SAME cached FC matrices script 13
   loads, through the SAME backbone + LRG path, compared to the values script 13
   wrote. Must agree to float tolerance. (Not a load-bearing science number --
   this is a harness-equivalence check, which is exactly what a prior CSV is
   legitimately used for.)

B. FC reconstruction: ``mean_f |Im C(f)|`` computed from the raw timeseries via
   ``complex_coherency_band`` compared to the cached ``imcoh_abs`` adjacency.
   If these disagree, a timeseries-level null is nulling a different pipeline
   than the one being gated.

C. The phase-rotation identity. For a per-channel circular shift by Delta_c
   samples, the frequency-domain claim is

       S_ij(f) -> S_ij(f) * exp(-2*pi*i*f*(Delta_i - Delta_j)/fs)

   i.e. |C_ij(f)| exactly preserved, phase rotated by the shift DIFFERENCE. We
   compare, on real data:
     C1. an exact-DFT single-segment case (nperseg == L, no Welch averaging):
         the identity should hold to float64 round-off. This is the regime
         where the identity is a theorem.
     C2. the actual Welch regime used by the pipeline (many overlapping Hann
         segments): here the identity is NOT exact, because a whole-recording
         roll moves samples across segment boundaries and the shifted signal is
         a different set of segments. We quantify the discrepancy rather than
         assume it away.

   The honest statement we need for the report is C1 (identity is exact in the
   single-segment DFT sense) plus a measured bound from C2 (how far the Welch
   estimator departs from it). N1 is then implemented in whichever form the
   evidence supports.
"""
from __future__ import annotations
import sys, time
import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))

from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, DEFAULT_SAMPLE_RATE, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, rho_sym_over_scales
from lrg_eegfc.utils.surrogate.coherency_surrogate import complex_coherency_band
from lrg_eegfc.utils.fc.coherence._common import welch_csd

from audit_150_rho_sym_gate import load_phase

PHASES = ("A", "B", "task_test", "rest_post")
SGRID = np.logspace(0.0, np.log10(180.0), 16)
FRAC = 0.20
BACKBONE = "mst020"


def eig_mst(W):
    return laplacian_eig(select_backbone(W, BACKBONE, frac=FRAC))


def check_A(cells):
    print("\n=== A. observed rho_sym reproduction vs script 13 ===", flush=True)
    ref = pd.read_csv(ROOT / "data" / "sparsified_arc" / "ms_mst020" / "per_patient_scale.csv")
    worst = 0.0
    for pat, band in cells:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES}
        obs = rho_sym_over_scales({ph: eig_mst(Ws[ph]) for ph in PHASES}, SGRID)
        r = ref[(ref.patient == pat) & (ref.band == band)].sort_values("s")
        if len(r) != len(SGRID):
            print(f"  {pat} {band}: no reference rows"); continue
        d = np.nanmax(np.abs(obs - r.obs_rho.values))
        worst = max(worst, float(d))
        print(f"  {pat:7s} {band:11s} max|drho| = {d:.3e}   "
              f"obs[s=1]={obs[0]:+.4f} ref={r.obs_rho.values[0]:+.4f}", flush=True)
    print(f"  -> worst discrepancy across cells: {worst:.3e}", flush=True)
    return worst


def check_B(pat, bands):
    print(f"\n=== B. FC from timeseries vs cached imcoh_abs ({pat}) ===", flush=True)
    fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE)
    nps = nperseg_for_fs(fs)
    X = np.asarray(load_timeseries(pat, "rest_post", SEEG_DATAPATH), float)
    if X.shape[0] > X.shape[1]:
        X = X.T
    print(f"  X = {X.shape}, fs = {fs}, nperseg = {nps}", flush=True)
    for band in bands:
        C = complex_coherency_band(X, fs, BRAIN_BANDS[band], nps)   # (Fb,N,N)
        W_ts = np.abs(np.imag(C)).mean(axis=0)
        np.fill_diagonal(W_ts, 0.0)
        W_ts = 0.5 * (W_ts + W_ts.T)
        W_cache = load_phase(pat, "rest_post", band)
        d = float(np.max(np.abs(W_ts - W_cache)))
        rel = d / float(np.max(np.abs(W_cache)))
        print(f"  {band:11s} max|dW| = {d:.3e}  (rel {rel:.2e})  "
              f"mean_cache={W_cache.mean():.4f} mean_ts={W_ts.mean():.4f}", flush=True)
    return X, fs, nps


def check_C(X, fs, nps, band):
    print(f"\n=== C. circular shift <-> phase rotation identity ({band}) ===", flush=True)
    rng = np.random.default_rng(0)
    lo, hi = BRAIN_BANDS[band]

    # ---- C1: exact single-segment DFT regime (nperseg == L, one segment) ----
    # welch_csd with nperseg = L gives exactly one Hann-windowed segment. A
    # circular roll of the *windowed* segment is what the identity describes;
    # rolling the raw signal then windowing is NOT the same operation, so for
    # the theorem-regime check we roll a signal whose length equals nperseg and
    # compare against the analytic ramp applied to the un-windowed DFT.
    Nch, L = X.shape
    Lseg = 4096
    Xs = X[:12, :Lseg]                                   # 12 channels is enough
    shifts = rng.integers(0, Lseg, size=Xs.shape[0])
    Xr = np.stack([np.roll(Xs[c], int(shifts[c])) for c in range(Xs.shape[0])])

    freqs = np.fft.rfftfreq(Lseg, 1.0 / fs)
    F0 = np.fft.rfft(Xs, axis=1)
    Fr = np.fft.rfft(Xr, axis=1)
    ramp = np.exp(-2j * np.pi * np.outer(shifts, np.arange(freqs.size)) / Lseg)
    d_dft = float(np.max(np.abs(Fr - F0 * ramp)) / np.max(np.abs(F0)))
    print(f"  C1 raw-DFT shift theorem : max rel err = {d_dft:.3e}   (expect ~1e-15)", flush=True)

    # coherency built from that single un-windowed DFT, observed vs rotated
    def coh_from_fft(Fx):
        S = np.einsum("nf,mf->nmf", Fx, np.conj(Fx))
        P = np.real(np.einsum("nnf->nf", S))
        den = np.sqrt(P[:, None, :] * P[None, :, :])
        return np.divide(S, den, out=np.zeros_like(S), where=den > 0)

    C0 = coh_from_fft(F0)
    Cr = coh_from_fft(Fr)
    m = (freqs >= lo) & (freqs <= hi)
    d_mag = float(np.max(np.abs(np.abs(Cr[:, :, m]) - np.abs(C0[:, :, m]))))
    print(f"  C1 |coherency| preserved : max|d|C||   = {d_mag:.3e}   (expect ~1e-15)", flush=True)

    # and the analytic per-pair ramp exp(-2pi i f (Di - Dj)/fs) reproduces Cr
    k = np.arange(freqs.size)
    pair_ramp = np.exp(-2j * np.pi * (shifts[:, None, None] - shifts[None, :, None])
                       * k[None, None, :] / Lseg)
    d_pair = float(np.max(np.abs(Cr - C0 * pair_ramp)))
    print(f"  C1 pair-ramp identity    : max|dC|     = {d_pair:.3e}   (expect ~1e-15)", flush=True)

    # ---- C2: the actual Welch regime the pipeline uses ----
    Xw = X[:24, :]
    sh = rng.integers(0, Xw.shape[1], size=Xw.shape[0])
    Xw_roll = np.stack([np.roll(Xw[c], int(sh[c])) for c in range(Xw.shape[0])])
    C_obs = complex_coherency_band(Xw, fs, (lo, hi), nps)
    C_roll = complex_coherency_band(Xw_roll, fs, (lo, hi), nps)
    a_obs = np.abs(C_obs); a_roll = np.abs(C_roll)
    iu = np.triu_indices(Xw.shape[0], 1)
    mo = a_obs[:, iu[0], iu[1]].mean()
    mr = a_roll[:, iu[0], iu[1]].mean()
    print(f"  C2 Welch, true time-domain roll:", flush=True)
    print(f"     mean|C| observed = {mo:.4f}  rolled = {mr:.4f}  ratio = {mr/mo:.3f}", flush=True)
    print(f"     -> magnitude is NOT preserved by a whole-recording roll under Welch", flush=True)

    # what a whole-recording roll does to <|Im C|>: the quantity the pipeline uses
    w_obs = np.abs(np.imag(C_obs)).mean(axis=0)
    w_roll = np.abs(np.imag(C_roll)).mean(axis=0)
    print(f"     <|ImC|> observed = {w_obs[iu].mean():.4f}  rolled = {w_roll[iu].mean():.4f}"
          f"  ratio = {w_roll[iu].mean()/w_obs[iu].mean():.3f}", flush=True)

    # ---- C3: the frequency-domain implementation actually used by N1 ----
    # Apply the per-pair phase ramp directly to the cached complex coherency.
    # This is exact by construction for |C| and is what makes N1 cheap.
    Cc = complex_coherency_band(Xw, fs, (lo, hi), nps)
    fb = np.fft.rfftfreq(nps, 1.0 / fs)
    fb = fb[(fb >= lo) & (fb <= hi)]
    dl = rng.integers(0, nps, size=Xw.shape[0])
    ramp2 = np.exp(-2j * np.pi * np.outer(fb, dl) / fs)              # (Fb, N)
    C_rot = Cc * ramp2[:, :, None] * np.conj(ramp2)[:, None, :]
    d_absC = float(np.max(np.abs(np.abs(C_rot) - np.abs(Cc))))
    print(f"  C3 freq-domain ramp on cached C: max|d|C|| = {d_absC:.3e}  (exact by construction)",
          flush=True)
    w_rot = np.abs(np.imag(C_rot)).mean(axis=0)
    print(f"     <|ImC|> observed = {w_obs[iu].mean():.4f}  ramped = {w_rot[iu].mean():.4f}"
          f"  ratio = {w_rot[iu].mean()/w_obs[iu].mean():.3f}", flush=True)
    print(f"     (2/pi)*<|C|> prediction = {(2/np.pi)*a_obs[:, iu[0], iu[1]].mean():.4f}", flush=True)


def main():
    t0 = time.time()
    check_A([("Pat_05", "beta"), ("Pat_06", "delta"), ("Pat_13", "alpha")])
    X, fs, nps = check_B("Pat_05", ["alpha", "beta"])
    check_C(X, fs, nps, "beta")
    print(f"\n[preflight] {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
