#!/usr/bin/env python3
"""W0-B: validate every primitive in ``utils/surrogate/timeseries_nulls.py``
against the production FC path, and calibrate the surrogate magnitudes.

Checks
------
V1. ``segment_fft`` + ``csd_from_segment_subset`` over ALL segments reproduces
    ``welch_csd`` (and hence the cached ``imcoh_abs`` adjacency) to float32
    precision. If this fails, every rung is nulling the wrong pipeline.
V2. ``lag_randomized_coherency`` preserves ``|C_ij(f)|`` exactly and leaves the
    coherency diagonal real-unit.
V3. ``phase_randomized_coherency`` same invariants.
V4. Magnitude calibration: what each rung does to ``<|ImC|>`` vs the observed,
    and how close N1/N2 come to the ``(2/pi)<|C|>`` prediction.
V5. Global-rescale invariance of the READOUT (backbone -> LRG -> rho_sym over
    scales). N1/N2 inflate W by a several-fold constant; this must be provably
    absorbed or those rungs are not comparable to the observed statistic.
V6. ``segment_shifted_csd`` preserves each channel's Welch PSD exactly.
V7. ``block_partition_indices`` size / contiguity / coverage invariants.
V8. Timing for each rung -> cohort extrapolation.
"""
from __future__ import annotations
import sys, time
import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))

from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, DEFAULT_SAMPLE_RATE, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, rho_sym_over_scales
from lrg_eegfc.utils.fc.coherence._common import welch_csd
from lrg_eegfc.utils.surrogate.timeseries_nulls import (
    band_bin_frequencies, block_partition_indices, coherency_from_csd,
    csd_from_segment_subset, global_rescale_invariance_residual,
    imcoh_abs_from_coherency, lag_randomized_coherency,
    phase_randomized_coherency, segment_fft, segment_shifted_csd,
)
from audit_150_rho_sym_gate import load_phase

PAT, BAND = "Pat_05", "beta"
SGRID = np.logspace(0.0, np.log10(180.0), 16)
PHASES = ("A", "B", "task_test", "rest_post")


def main():
    fs = FS_OVERRIDES.get(PAT, DEFAULT_SAMPLE_RATE)
    nps = nperseg_for_fs(fs)
    lo, hi = BRAIN_BANDS[BAND]
    rng = np.random.default_rng(20260825)

    X = np.asarray(load_timeseries(PAT, "rest_post", SEEG_DATAPATH), float)
    if X.shape[0] > X.shape[1]:
        X = X.T
    print(f"[{PAT}/{BAND}] X={X.shape} fs={fs} nperseg={nps}", flush=True)

    # ---------------- V1 ----------------
    t = time.time()
    freqs_b, F, scale = segment_fft(X, fs, nps, band=(lo, hi))
    t_fft = time.time() - t
    print(f"\nV1 segment_fft: F={F.shape} {F.nbytes/1e6:.0f} MB in {t_fft:.1f}s", flush=True)
    t = time.time()
    S = csd_from_segment_subset(F, None, scale)
    t_csd = time.time() - t
    C = coherency_from_csd(S)
    W_lib = imcoh_abs_from_coherency(C)
    W_cache = load_phase(PAT, "rest_post", BAND)
    d = float(np.max(np.abs(W_lib - W_cache)))
    print(f"   full-subset CSD in {t_csd:.2f}s ; max|W_lib - W_cache| = {d:.3e} "
          f"(rel {d/W_cache.max():.2e})", flush=True)

    # ---------------- V2 / V3 ----------------
    freqs = band_bin_frequencies(fs, nps, (lo, hi))
    assert np.allclose(freqs, freqs_b), "band bin mismatch"
    C1 = lag_randomized_coherency(C, freqs, fs, rng, nperseg=nps)
    C2 = phase_randomized_coherency(C, rng)
    for tag, Cs in (("V2 N1", C1), ("V3 N2", C2)):
        dm = float(np.max(np.abs(np.abs(Cs) - np.abs(C))))
        dd = float(np.max(np.abs(np.diagonal(Cs, axis1=1, axis2=2) - 1.0)))
        print(f"{tag}: max|d|C||={dm:.3e}  max|diag-1|={dd:.3e}", flush=True)

    # ---------------- V4 ----------------
    iu = np.triu_indices(C.shape[1], 1)
    absC = np.abs(C).mean(axis=0)
    print(f"\nV4 magnitude calibration (mean over off-diagonal pairs):", flush=True)
    print(f"   observed <|ImC|>      = {W_cache[iu].mean():.5f}", flush=True)
    print(f"   <|C|>                 = {absC[iu].mean():.5f}", flush=True)
    print(f"   (2/pi)<|C|> predicted = {(2/np.pi)*absC[iu].mean():.5f}", flush=True)
    W1 = imcoh_abs_from_coherency(C1); W2 = imcoh_abs_from_coherency(C2)
    print(f"   N1  <|ImC|> = {W1[iu].mean():.5f}  (x{W1[iu].mean()/W_cache[iu].mean():.2f})",
          flush=True)
    print(f"   N2  <|ImC|> = {W2[iu].mean():.5f}  (x{W2[iu].mean()/W_cache[iu].mean():.2f})",
          flush=True)
    Sb = segment_shifted_csd(F, rng, scale, min_shift=1)
    Wb = imcoh_abs_from_coherency(coherency_from_csd(Sb))
    absCb = np.abs(coherency_from_csd(Sb)).mean(axis=0)
    print(f"   N1b <|ImC|> = {Wb[iu].mean():.5f}  (x{Wb[iu].mean()/W_cache[iu].mean():.2f})"
          f"   and <|C|> = {absCb[iu].mean():.5f} (vs {absC[iu].mean():.5f} observed)", flush=True)
    print("   -> N1/N2 INFLATE (lag randomized, |C| kept); N1b COLLAPSES |C| "
          "(independence null). Different H0s.", flush=True)

    # realization spread of N1 vs N2 (the variance claim in the module docstring)
    m1 = [imcoh_abs_from_coherency(lag_randomized_coherency(C, freqs, fs, rng, nps))[iu].mean()
          for _ in range(20)]
    m2 = [imcoh_abs_from_coherency(phase_randomized_coherency(C, rng))[iu].mean()
          for _ in range(20)]
    print(f"   realization sd over 20 draws: N1 {np.std(m1):.2e}   N2 {np.std(m2):.2e}"
          f"   (N2 is the low-variance limit of N1)", flush=True)

    # ---------------- V5 ----------------
    print(f"\nV5 global-rescale invariance of the readout:", flush=True)
    Wobs = {ph: load_phase(PAT, ph, BAND) for ph in PHASES}

    def readout(Wd):
        return rho_sym_over_scales(
            {ph: laplacian_eig(select_backbone(Wd[ph], "mst020", frac=0.20))
             for ph in PHASES}, SGRID)

    r0 = readout(Wobs)
    res = global_rescale_invariance_residual(readout, Wobs, c=6.5)
    print(f"   max|rho_sym(W) - rho_sym(6.5W)| over 16 scales = {res:.3e}", flush=True)
    print(f"   rho_sym(s=1)={r0[0]:+.4f}  rho_sym(s=180)={r0[-1]:+.4f}", flush=True)

    # ---------------- V6 ----------------
    print(f"\nV6 segment-shift preserves per-channel PSD:", flush=True)
    P_obs = np.real(np.diagonal(S, axis1=1, axis2=2))
    P_sh = np.real(np.diagonal(Sb, axis1=1, axis2=2))
    print(f"   max rel |dPSD| = {np.max(np.abs(P_sh-P_obs)/np.maximum(P_obs,1e-30)):.3e}",
          flush=True)

    # ---------------- V7 ----------------
    print(f"\nV7 block_partition_indices invariants:", flush=True)
    sizes = {"A": 10, "B": 10, "task_test": 37, "rest_post": 21}
    nb = sum(sizes.values())
    for mode in ("identity", "free", "order_preserving"):
        bi = block_partition_indices(nb, sizes, rng, mode=mode)
        allb = np.sort(np.concatenate(list(bi.values())))
        ok_sz = all(len(bi[k]) == sizes[k] for k in sizes)
        ok_cov = np.array_equal(allb, np.arange(nb))
        contig = all(np.all(np.diff(np.sort(v)) == 1) or mode == "free"
                     for v in bi.values())
        print(f"   {mode:17s} sizes_ok={ok_sz} covers_all={ok_cov} "
              f"contiguous_or_free={contig}  A={bi['A'][:4]}...", flush=True)
    rots = {tuple(block_partition_indices(nb, sizes, rng, "order_preserving", rotation=r)["A"])
            for r in range(nb)}
    print(f"   order_preserving distinct rotations = {len(rots)} (= n_blocks, exact test)",
          flush=True)

    # ---------------- V8 ----------------
    print(f"\nV8 timing per surrogate ({BAND}, N={X.shape[0]}, n_seg={F.shape[1]}):", flush=True)
    def t_of(fn, n=10):
        t = time.time()
        for _ in range(n):
            fn()
        return (time.time() - t) / n
    t1 = t_of(lambda: imcoh_abs_from_coherency(lag_randomized_coherency(C, freqs, fs, rng, nps)))
    t2 = t_of(lambda: imcoh_abs_from_coherency(phase_randomized_coherency(C, rng)))
    tb = t_of(lambda: imcoh_abs_from_coherency(coherency_from_csd(
        segment_shifted_csd(F, rng, scale))), n=3)
    t3 = t_of(lambda: imcoh_abs_from_coherency(coherency_from_csd(
        csd_from_segment_subset(F, np.sort(rng.choice(F.shape[1], F.shape[1]//3, False)), scale))),
        n=3)
    tro = t_of(lambda: readout(Wobs), n=3)
    print(f"   N1  coherency ramp   : {t1*1000:7.1f} ms  (x4 phases = {4*t1*1000:.0f} ms)",
          flush=True)
    print(f"   N2  phase random     : {t2*1000:7.1f} ms", flush=True)
    print(f"   N1b segment shift    : {tb*1000:7.1f} ms  (x4 phases = {4*tb*1000:.0f} ms)",
          flush=True)
    print(f"   N3  subset CSD (1/3) : {t3*1000:7.1f} ms  (x4 phases = {4*t3*1000:.0f} ms)",
          flush=True)
    print(f"   READOUT backbone+LRG : {tro*1000:7.1f} ms  <- dominates", flush=True)
    for tag, tt in (("N1", 4*t1), ("N2", 4*t2), ("N1b", 4*tb), ("N3", 4*t3)):
        per_surr = tt + tro
        print(f"   -> {tag:4s} R=200: {per_surr*200:6.0f} s/cell ; "
              f"60 cells /12 workers = {per_surr*200*60/12/60:5.1f} min", flush=True)


if __name__ == "__main__":
    main()
