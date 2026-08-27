#!/usr/bin/env python3
"""Equivalence check for the new library cross-phase functionals.

Asserts, on real FC (not mocked -- hypothesis-level tests never mock FC):
  (a) cross_phase_functionals_over_scales(...)["T_probe"] == rho_sym_over_scales(...)
      on the same four phases (the incumbent estimator is unchanged);
  (b) the four five-phase functionals reproduce the reference implementation in
      scripts/01_compute/sparsified_arc/05_enc_inf_arc.py (repeated pairwise
      spearmanr + partial), to floating-point tolerance;
  (c) the four-phase call path works with task_learn absent (phase set as data).
"""
from __future__ import annotations
import sys
import numpy as np
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import (
    cophenetic_at_scale, cross_phase_functionals_over_scales, laplacian_eig,
    rho_sym_over_scales,
)
from lrg_eegfc.workflow.fc import load_fc_matrix

HALVES = ROOT / "data" / "paper_final" / "w0a_substrate" / "halves"
SGRID = np.logspace(0.0, np.log10(180.0), 16)
PH5 = ("A", "B", "task_learn", "task_test", "rest_post")


def load_phase(pat, phase, band, tk="abs"):
    if phase in ("A", "B"):
        W = np.load(HALVES / pat / f"{band}_rest_pre_{phase}_imcoh_{tk}.npy")
    else:
        W = load_fc_matrix(pat, phase, band, fc_method=f"imcoh_{tk}")
    W = np.asarray(W, float); np.fill_diagonal(W, 0.0); W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


def _rho(a, b):
    return spearmanr(a, b).statistic


def _partial(a, b, c):
    rab, rac, rbc = _rho(a, b), _rho(a, c), _rho(b, c)
    den = np.sqrt(max(0.0, (1 - rac ** 2) * (1 - rbc ** 2)))
    return (rab - rac * rbc) / den if den > 0 else np.nan


def reference_functionals(eig5, s_grid):
    """Verbatim reference: 05_enc_inf_arc.functionals_over_scales."""
    keys = ("T_test", "T_learn", "T_infspec", "T_infspec_pe")
    out = {k: np.full(len(s_grid), np.nan) for k in keys}
    for i, s in enumerate(s_grid):
        D = {ph: cophenetic_at_scale(*eig5[ph], s) for ph in PH5}
        e = D["task_learn"] - D["A"]; e2 = D["task_learn"] - D["B"]
        f = D["task_test"] - D["task_learn"]
        g = D["task_test"] - D["A"]; g2 = D["task_test"] - D["B"]
        p = D["rest_post"] - D["B"]; p2 = D["rest_post"] - D["A"]
        out["T_test"][i] = 0.5 * (_rho(g, p) + _rho(g2, p2))
        out["T_learn"][i] = 0.5 * (_rho(e, p) + _rho(e2, p2))
        out["T_infspec"][i] = 0.5 * (_rho(f, p) + _rho(f, p2))
        out["T_infspec_pe"][i] = 0.5 * (_partial(f, p, e) + _partial(f, p2, e2))
    return out


def main():
    bad = 0
    for pat, band in (("Pat_02", "beta"), ("Pat_13", "alpha"), ("Pat_15", "low_gamma")):
        Ws = {ph: load_phase(pat, ph, band) for ph in PH5}
        eig5 = {ph: laplacian_eig(select_backbone(Ws[ph], "mst", frac=0.20)) for ph in PH5}

        new = cross_phase_functionals_over_scales(eig5, SGRID)
        ref = reference_functionals(eig5, SGRID)
        inc = rho_sym_over_scales({k: eig5[k] for k in ("A", "B", "task_test", "rest_post")},
                                  SGRID)
        eig4 = {k: eig5[k] for k in ("A", "B", "task_test", "rest_post")}
        new4 = cross_phase_functionals_over_scales(eig4, SGRID)

        d_inc = np.nanmax(np.abs(new["T_probe"] - inc))
        d_4 = np.nanmax(np.abs(new4["T_probe"] - inc))
        pairs = (("T_probe", "T_test"), ("T_encode", "T_learn"),
                 ("T_probespec", "T_infspec"), ("T_probespec_pe", "T_infspec_pe"))
        ds = {a: float(np.nanmax(np.abs(new[a] - ref[b]))) for a, b in pairs}
        keys4 = sorted(new4)
        ok = (d_inc < 1e-9 and d_4 < 1e-9 and max(ds.values()) < 1e-9
              and keys4 == ["T_probe"])
        bad += (not ok)
        print(f"[{pat}/{band}] max|T_probe - rho_sym| = {d_inc:.2e} (5-phase call), "
              f"{d_4:.2e} (4-phase call); vs 05_enc_inf_arc reference: "
              + ", ".join(f"{k}={v:.2e}" for k, v in ds.items())
              + f"; 4-phase keys={keys4}  [{'OK' if ok else 'MISMATCH'}]", flush=True)
    print("\nALL EQUIVALENCE CHECKS PASSED" if bad == 0 else f"\n{bad} CELL(S) FAILED", flush=True)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
