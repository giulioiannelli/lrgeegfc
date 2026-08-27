#!/usr/bin/env python3
r"""W0-C: how much of the hierarchy is left after residualizing it on the raw edges.

Interpretation control for the DETECTION test. A null ``coph|raw`` means "the
hierarchy adds no persistent structure beyond raw" only if the residual still
*has* structure to test. If residualizing on raw destroys 97% of the cophenetic
variance, a null result is a statement about power, not about the hierarchy.

Per (patient, band, scale), observed graphs only, no surrogates:

  var_keep_lin  variance of rank(C) surviving the rank-linear residualization
  var_keep_np   the same for the binned conditional-mean residualization
  rho_CA        Spearman between the cophenetic distances and the raw edges

Output: data/paper_final/w0c_gate_tau/residual_diagnostics.csv
"""
from __future__ import annotations

import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import cophenetic_at_scale, laplacian_eig
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.phase_graphs import CROSS_PHASE_PHASES, phase_graphs

ROOT = setup_script_env()
BASE = Path(os.environ.get("W0C_BASE", ROOT / "data" / "paper_final" / "w0c_gate_tau"))

import importlib.util as _ilu                                          # noqa: E402
_spec = _ilu.spec_from_file_location(
    "w0c_grid", Path(__file__).with_name("w0c_01_gate_and_tau_grid.py"))
_g = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_g)

FC_METHOD = _g.FC_METHOD
BACKBONE, FRAC, SGRID = _g.BACKBONE, _g.FRAC, _g.SGRID


def main():
    BASE.mkdir(parents=True, exist_ok=True)
    rows = []
    t0 = time.time()
    n = len(PATIENTS_4PHASE)
    for ip, pat in enumerate(PATIENTS_4PHASE, 1):
        for band in BRAIN_BANDS_NAMES:
            try:
                Ws = phase_graphs(pat, band, fc_method=FC_METHOD)
            except Exception:                                          # noqa: BLE001
                continue
            for ph in CROSS_PHASE_PHASES:
                A = _g._triu(Ws[ph])
                B = Ws[ph] if BACKBONE == "dense" else select_backbone(
                    Ws[ph], BACKBONE, frac=FRAC)
                ev, V = laplacian_eig(B)
                bins = _g.bin_index(A)
                rA = _g._cr(A)
                for s in SGRID:
                    try:
                        C = cophenetic_at_scale(ev, V, s)
                    except Exception:                                  # noqa: BLE001
                        continue
                    rC = _g._cr(C)
                    v0 = float(rC @ rC)
                    if v0 <= 0:
                        continue
                    lin = _g.resid_rank_linear(C, A)
                    npr = _g.resid_binned(C, bins)
                    rows.append(dict(
                        patient=pat, band=band, phase=ph, s=float(s),
                        rho_CA=float(rA @ rC / np.sqrt((rA @ rA) * v0)),
                        var_keep_lin=float(lin @ lin) / v0,
                        var_keep_np=float(npr @ npr) / v0))
        el = time.time() - t0
        print(f"[{ip}/{n}] {pat} {el:5.1f}s ETA {el/ip*(n-ip):5.1f}s", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(BASE / "residual_diagnostics.csv", index=False)
    print("\n  cophenetic variance surviving residualization on the raw edges "
          "(cohort median over patients and phases):", flush=True)
    print(f"  {'band':11s}" + "".join(f"  s={v:<7.2f}" for v in
                                      (0.05, 1.0, 5.6, 30.0, 180.0)), flush=True)
    for band in BRAIN_BANDS_NAMES:
        x = df[df.band == band]
        cells = []
        for v in (0.05, 1.0, 5.6, 30.0, 180.0):
            sv = SGRID[int(np.argmin(np.abs(SGRID - v)))]
            y = x[np.isclose(x.s, sv)]
            cells.append(f"  {y.var_keep_lin.median():.2f}/{y.var_keep_np.median():.2f}")
        print(f"  {band:11s}" + "".join(f"{c:<11s}" for c in cells), flush=True)
    print("    (rank-linear / binned; 1.00 = nothing removed, 0.00 = fully absorbed)",
          flush=True)
    print(f"\n  |Spearman(cophenetic, raw edges)| cohort median: "
          f"{df.rho_CA.abs().median():.3f}", flush=True)
    print(f"\n[w0c-resid] -> {BASE/'residual_diagnostics.csv'}", flush=True)


if __name__ == "__main__":
    main()
