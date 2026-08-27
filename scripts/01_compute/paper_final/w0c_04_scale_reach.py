#!/usr/bin/env python3
r"""W0-C: the physical companion to the scale axis -- how far the diffusion has reached.

Descriptive geometry, no null. ``N_eff(s)`` says how many components the process
still resolves at scale ``s``; this says how many millimetres it spans. Both are
needed before a scale may be described in words.

  ell(s) = max{ ||x_i - x_j|| : K_ij(s) >= 1/N }

the largest physical separation between two contacts whose heat exchange has
reached at least the uniform-equilibrium share ``1/N`` -- an actual length, not
a heat-weighted average. Reported against the electrode pitch (median
nearest-contact spacing) and the implant span (largest pairwise separation).

Substrate injected through the same environment variables as the master grid.
Output: data/paper_final/w0c_gate_tau/scale_reach.csv
"""
from __future__ import annotations

import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig
from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.phase_graphs import load_phase_fc

ROOT = setup_script_env()
BASE = Path(os.environ.get("W0C_BASE", ROOT / "data" / "paper_final" / "w0c_gate_tau"))
FC_METHOD = os.environ.get("W0C_FC_METHOD", "imcoh_abs")
BACKBONE = os.environ.get("W0C_BACKBONE", "mst020")
FRAC = float(os.environ.get("W0C_FRAC", "0.20"))
PHASE = "rest_post"
SGRID = np.logspace(np.log10(0.05), np.log10(180.0), 28)


def main():
    BASE.mkdir(parents=True, exist_ok=True)
    rows = []
    t0 = time.time()
    for ip, pat in enumerate(PATIENTS_4PHASE, 1):
        rdf = load_channel_regions(pat)
        xyz = rdf[["x", "y", "z"]].to_numpy(float) / 1000.0          # um -> mm
        mask = np.isfinite(xyz).all(axis=1)
        P = xyz[mask]
        Dm = cdist(P, P)
        span = float(Dm.max())
        Dnn = Dm.copy()
        np.fill_diagonal(Dnn, np.inf)
        pitch = float(np.median(Dnn.min(1)))
        for band in BRAIN_BANDS_NAMES:
            W = load_phase_fc(pat, PHASE, band, fc_method=FC_METHOD)
            if W.shape[0] != xyz.shape[0]:
                continue
            B = W if BACKBONE == "dense" else select_backbone(W, BACKBONE, frac=FRAC)
            ev, V = laplacian_eig(B)
            thr = 1.0 / ev.size
            for js, s in enumerate(SGRID):
                K = (V * np.exp(-(s / ev[-1]) * ev)) @ V.T
                comm = K[np.ix_(mask, mask)] >= thr
                np.fill_diagonal(comm, False)
                rows.append(dict(patient=pat, band=band, s=float(s), i_s=js,
                                 ell_mm=float(Dm[comm].max()) if comm.any() else np.nan,
                                 pitch_mm=pitch, span_mm=span, N=int(ev.size)))
        el = time.time() - t0
        print(f"[{ip}/{len(PATIENTS_4PHASE)}] {pat} {el:5.1f}s "
              f"ETA {el/ip*(len(PATIENTS_4PHASE)-ip):5.1f}s", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(BASE / "scale_reach.csv", index=False)
    pitch = df.pitch_mm.median()
    span = df.span_mm.median()
    print(f"\n  pitch ~ {pitch:.1f} mm   implant span ~ {span:.0f} mm", flush=True)
    print(f"  {'s':>8}  {'ell (mm, cohort median)':>24}  {'% of span':>10}", flush=True)
    for js, s in enumerate(SGRID):
        x = df[df.i_s == js].ell_mm
        v = np.nanmedian(x)
        if js % 3 == 0 or js == SGRID.size - 1:
            print(f"  {s:8.2f}  {v:24.1f}  {100*v/span:9.0f}%  "
                  f"(connected in {int(np.isfinite(x).sum())}/{len(x)} cells)", flush=True)
    print(f"\n[w0c-reach] -> {BASE/'scale_reach.csv'}", flush=True)


if __name__ == "__main__":
    main()
