#!/usr/bin/env python3
"""D1 -- entropy S(tau) and specific heat C(tau) on the percolation backbone.

For every (patient, band, phase): percolation backbone -> combinatorial Laplacian
-> normalised entropy Shat(tau) in [0,1] (auto-ranged so Shat spans the full
0->1 sigmoid) and specific heat C(tau) = -dShat/dlog10 tau. The C(tau) peak
count is the multiscale diagnostic (Villegas 2025: 1 peak = collapsed single
scale; >1 = multiscale ladder) -- i.e. whether the backbone even HAS the
multiscale structure the tau-sweep needs.

Outputs:
    data/sparsified_arc/entropy/{band}/{patient}_{phase}.npz   (s, tau, S, C, peaks)
    data/sparsified_arc/entropy/summary.csv
    data/sparsified_arc/entropy/config.json
"""
from __future__ import annotations
import json
import sys
import time

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import load_phase, COHORT, BANDS
from lrg_eegfc.utils.fc.backbone import percolation_backbone
from lrg_eegfc.utils.fc.heat_multiscale import (
    laplacian_eig, entropy_specific_heat, specific_heat_peaks,
)

PHASES = ("A", "B", "task_test", "rest_post")
OUT = ROOT / "data" / "sparsified_arc" / "entropy"


def run_cell(pat, phase, band):
    W = load_phase(pat, phase, band)
    B, _ = percolation_backbone(W)
    ev, _ = laplacian_eig(B)
    res = entropy_specific_heat(ev, n=300)
    pk = specific_heat_peaks(res)
    cell = OUT / band
    cell.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cell / f"{pat}_{phase}.npz",
        s=res["s"], tau=res["tau"], S=res["S"], C=res["C"],
        tau_min=res["tau_min"], tau_max=res["tau_max"],
        lambda_max=res["lambda_max"], lambda_2=res["lambda_2"],
        s_peaks=pk["s_peaks"], n_peaks=pk["n_peaks"], dominant_s=pk["dominant_s"],
    )
    return dict(patient=pat, phase=phase, band=band, N_nodes=int(res["N"]),
                lambda_max=res["lambda_max"], lambda_2=res["lambda_2"],
                s_window=float(res["lambda_max"] / max(res["lambda_2"], 1e-30)),
                n_peaks=int(pk["n_peaks"]), dominant_s=float(pk["dominant_s"]),
                S_at_tau_min=float(res["S"][np.argmin(np.abs(res["s"] - 1.0))]))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(p, ph, b) for b in BANDS for p in COHORT for ph in PHASES]
    print(f"[D1-entropy] {len(jobs)} cells", flush=True)
    t0 = time.time(); rows = []
    for i, (p, ph, b) in enumerate(jobs, 1):
        try:
            rows.append(run_cell(p, ph, b))
        except Exception as e:
            print(f"[{i}/{len(jobs)}] SKIP {p}/{ph}/{b}: {e}", flush=True)
        if i % 40 == 0 or i == len(jobs):
            print(f"[{i}/{len(jobs)}] {time.time()-t0:.0f}s", flush=True)
    df = pd.DataFrame(rows).sort_values(["band", "patient", "phase"])
    df.to_csv(OUT / "summary.csv", index=False)
    print("\n=== C(tau) peak count per band (multiscale ladder diagnostic) ===", flush=True)
    print(f"{'band':<11} {'median n_peaks':>15} {'frac >1 peak':>14} {'median dom. s':>15}", flush=True)
    for b in BANDS:
        x = df[df.band == b]
        if x.empty:
            continue
        print(f"{b:<11} {x.n_peaks.median():>15.1f} {(x.n_peaks > 1).mean()*100:>13.0f}% "
              f"{x.dominant_s.median():>15.1f}", flush=True)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="D1_entropy_specific_heat", backbone="percolation (theta*)",
        cohort=COHORT, bands=BANDS, phases=list(PHASES), n_tau=300,
        entropy="normalised von Neumann Shat=S/lnN in [0,1]; grid auto-ranged to span Shat~[0.001,0.999]",
        peak_rule="find_peaks on C(tau), rel_prominence=0.05*Cmax",
    ), indent=2))
    print(f"\n[D1-entropy] {len(df)} cells in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
