#!/usr/bin/env python3
"""FEASIBILITY PROBE (2026-09-08) -- NOT A TEST. No null, no knob integration, no LOO.
Establishes nothing; it exists so the numbers quoted in
.agents/plans/active/2026-09-08_new-paths-to-objectives.md can be regenerated.
The methodology it probes is scoped (5-point preamble) in
.agents/guides/task-persistence-investigation/2026-09-08_hierarchy-level-decomposition.md.
Run from the worktree with PYTHONPATH=src and the lapbrain python. Output lands under
data/paper_final/feasibility/ (copy from the ROOT the script resolves).
"""
"""Spectral-range audit: how many decades of non-trivial diffusion scale each graph has.

Per (patient, band, phase in {A, rest_post}, backbone) row:
  N, lambda_2, lambda_max, ratio, decades=log10(ratio), spectral PR over nonzero ev,
  n_zero_eigs (connected components), n_peaks of C(tau), s_peak1..4 (ascending s)
  + prom_peak1..4, dominant_s, frac ev < lambda_max/10, < lambda_max/100,
  n_eff_heat = participation ratio of the 12x12 Spearman matrix of triu rho(s)
  across s in logspace(0, log10(200), 12).
"""
from __future__ import annotations

import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.stats import spearmanr

from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.workflow.phase_graphs import load_phase_fc
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import (
    laplacian_eig, entropy_specific_heat, specific_heat_peaks,
)

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("spectral_range_audit.csv")
PATIENTS = list(PATIENTS_4PHASE)
if len(sys.argv) > 2:                       # optional patient subset for timing
    PATIENTS = sys.argv[2].split(",")
BANDS = list(BRAIN_BANDS_NAMES)
PHASES = ("A", "rest_post")
BACKBONES = [
    ("dense", dict(kind="dense")),
    ("mst0.07", dict(kind="mst", frac=0.07)),
    ("mst0.10", dict(kind="mst", frac=0.10)),
    ("mst0.14", dict(kind="mst", frac=0.14)),
    ("mst0.20", dict(kind="mst", frac=0.20)),
    ("tmfg", dict(kind="tmfg")),
    ("disparity0.20", dict(kind="disparity", disparity_alpha=0.20)),
]
S_GRID = np.logspace(0.0, np.log10(200.0), 12)
MAX_PEAKS = 4
REL_PROM = 0.05                              # same as specific_heat_peaks default


def n_eff_heat(ev: np.ndarray, V: np.ndarray, lam_max: float) -> float:
    """Participation ratio of the 12x12 Spearman matrix of triu(rho(s)) across s."""
    N = ev.size
    r, c = np.triu_indices(N, k=1)
    X = np.empty((S_GRID.size, r.size))
    for i, s in enumerate(S_GRID):
        tau = s / lam_max
        rho = (V * np.exp(-tau * ev)) @ V.T
        rho /= np.trace(rho)
        X[i] = rho[r, c]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        R = spearmanr(X, axis=1).correlation
    R = np.asarray(R, float)
    R = np.nan_to_num(R, nan=0.0)
    np.fill_diagonal(R, 1.0)
    eig = np.linalg.eigvalsh(0.5 * (R + R.T))
    return float(eig.sum() ** 2 / np.sum(eig ** 2))


def audit_one(B: np.ndarray) -> dict:
    ev, V = laplacian_eig(B)
    N = ev.size
    lam_max = float(ev[-1])
    lam2 = float(ev[1])
    tol = 1e-10 * lam_max
    nz = ev[ev > tol]
    n_zero = int(N - nz.size)
    ratio = lam_max / lam2 if lam2 > tol else np.inf
    row = dict(
        N=N, lambda_2=lam2, lambda_max=lam_max, ratio=ratio,
        decades=float(np.log10(ratio)) if np.isfinite(ratio) else np.inf,
        n_zero_eigs=n_zero,
        spectral_pr=float(nz.sum() ** 2 / np.sum(nz ** 2)),
        frac_ev_below_lmax_10=float(np.mean(ev < lam_max / 10.0)),
        frac_ev_below_lmax_100=float(np.mean(ev < lam_max / 100.0)),
        density=float((B[np.triu_indices(N, 1)] > 0).mean()),
    )
    # specific heat + peaks (library), prominences re-derived with the same rule
    res = entropy_specific_heat(ev)
    pk = specific_heat_peaks(res, rel_prominence=REL_PROM)
    C = res["C"]
    prom_thr = max(REL_PROM * float(C.max()), 1e-9)
    idx, props = find_peaks(C, prominence=prom_thr)
    if idx.size == 0:                        # library fallback: argmax, no prominence
        idx = np.array([int(np.argmax(C))])
        proms = np.array([np.nan])
    else:
        proms = props["prominences"]
    order = np.argsort(res["s"][idx])        # ascending s
    s_pk = res["s"][idx][order]
    p_pk = proms[order]
    row["n_peaks"] = int(pk["n_peaks"])
    row["dominant_s"] = float(pk["dominant_s"])
    row["fallback_argmax"] = bool(props is None or len(props.get("prominences", [])) == 0)
    for k in range(MAX_PEAKS):
        row[f"s_peak{k+1}"] = float(s_pk[k]) if k < s_pk.size else np.nan
        row[f"prom_peak{k+1}"] = float(p_pk[k]) if k < p_pk.size else np.nan
    row["s_range_lo"] = float(res["s"][0])
    row["s_range_hi"] = float(res["s"][-1])
    row["n_eff_heat"] = n_eff_heat(ev, V, lam_max)
    return row


def main() -> None:
    cells = [(p, b, ph) for p in PATIENTS for b in BANDS for ph in PHASES]
    n_tot = len(cells) * len(BACKBONES)
    rows, i, t0 = [], 0, time.time()
    for p, b, ph in cells:
        W = load_phase_fc(p, ph, b)
        for name, kw in BACKBONES:
            i += 1
            row = dict(patient=p, band=b, phase=ph, backbone=name, error="")
            try:
                B = select_backbone(W, **kw)
                row.update(audit_one(B))
            except Exception as e:           # record NaN and continue
                row["error"] = f"{type(e).__name__}: {e}"[:200]
            rows.append(row)
            el = time.time() - t0
            print(f"[{i}/{n_tot}] {p} {b} {ph} {name} elapsed {el:6.1f}s "
                  f"ETA {el / i * (n_tot - i):6.1f}s", flush=True)
    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"wrote {OUT} ({len(df)} rows) wall {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
