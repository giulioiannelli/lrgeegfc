#!/usr/bin/env python3
"""Timing probe for W0-A: cost of one Welch complex-coherency pass and one trace cell.

Not a deliverable. Times (a) load_timeseries + complex_coherency_bands on one
rest_pre half, and (b) one (patient, band) matched-strength trace cell at R=10,
so the full A1 / A2 runtimes can be extrapolated BEFORE launch
(feedback_optimize_time_and_surface_progress).
"""
from __future__ import annotations
import time, sys
import numpy as np

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))

from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, DEFAULT_SAMPLE_RATE, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.surrogate.coherency_surrogate import complex_coherency_bands

t0 = time.time()
X = load_timeseries("Pat_02", "rest_pre", SEEG_DATAPATH)
print(f"load_timeseries: {time.time()-t0:.1f}s  shape={X.shape}", flush=True)
if X.shape[0] > X.shape[1]:
    X = X.T
fs = FS_OVERRIDES.get("Pat_02", DEFAULT_SAMPLE_RATE)
nper = max(256, nperseg_for_fs(fs) // 2)
T = X.shape[1]
t0 = time.time()
C = complex_coherency_bands(np.ascontiguousarray(X[:, : T // 2]), fs, BRAIN_BANDS, nper)
dt = time.time() - t0
tot_mb = sum(v.nbytes for v in C.values()) / 1e6
print(f"complex_coherency_bands (half A): {dt:.1f}s  bins={{{', '.join(f'{k}:{v.shape[0]}' for k,v in C.items())}}}  mem={tot_mb:.0f}MB", flush=True)
print(f"  -> per patient (2 halves) ~{2*dt:.0f}s; cohort n=10 ~{20*dt/60:.1f} min", flush=True)

# --- (b) one trace cell timing ---
from audit_150_rho_sym_gate import load_phase
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, rho_sym_over_scales
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle

SGRID = np.logspace(0.0, np.log10(180.0), 16)
PHASES = ("A", "B", "task_test", "rest_post")
Ws = {ph: load_phase("Pat_02", ph, "beta") for ph in PHASES}
N = Ws["A"].shape[0]
_ = matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))
rng = np.random.default_rng(0)
n_swaps = 20 * N * (N - 1) // 2

t0 = time.time()
Wsh = {ph: matched_strength_shuffle(Ws[ph], n_swaps, rng, 1.0) for ph in PHASES}
t_shuf = time.time() - t0
t0 = time.time()
eig = {ph: laplacian_eig(select_backbone(Wsh[ph], "mst020", frac=0.20)) for ph in PHASES}
t_eig = time.time() - t0
t0 = time.time()
_ = rho_sym_over_scales(eig, SGRID)
t_rho = time.time() - t0
print(f"\nper surrogate: shuffle(4ph)={t_shuf*1000:.0f}ms  backbone+eig(4ph)={t_eig*1000:.0f}ms  "
      f"rho_sym_over_16scales={t_rho*1000:.0f}ms  total={(t_shuf+t_eig+t_rho)*1000:.0f}ms", flush=True)
per_cell_R = t_shuf + t_eig + t_rho
print(f"  cell at R=200: {per_cell_R*200:.0f}s ; 60 cells / 12 workers: {per_cell_R*200*60/12/60:.1f} min per backbone config", flush=True)
print(f"  if shuffles are PRE-GENERATED and shared across configs, marginal per config: "
      f"{(t_eig+t_rho)*200*60/12/60:.1f} min", flush=True)
