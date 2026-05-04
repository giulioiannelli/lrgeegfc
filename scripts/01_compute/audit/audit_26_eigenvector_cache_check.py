#!/usr/bin/env python3
"""Audit 26 — eigenvector cache smoke check (Phase 0 gate for E-series).

For every cached LRG NPZ in the n=10 cohort (full-phase + halves),
verify the persisted ``eigenvalues`` / ``eigenvectors`` satisfy the
spectral identity::

    L_recon = V @ diag(λ) @ V.T  ≈  L̂ (to float tolerance)

with ``L̂`` recomputed on the fly from the original FC matrix via the
same lrgsglib pipeline that ``compute_lrg_analysis`` uses.

Three pass criteria per (patient, band, phase, source):
    * orthonormality:  ``‖V^T V - I‖_F < 1e-10``
    * trivial mode:    ``|λ_1| < 1e-8``
    * reconstruction:  ``‖L_recon - L̂‖_F / ‖L̂‖_F < 1e-10``

Usage
-----
    python scripts/01_compute/audit/audit_26_eigenvector_cache_check.py [-v]

Outputs
-------
    data/audit/eigenvector_cache_check/audit_26_eigenvector_cache_n10_imcoh_abs.csv
"""
from __future__ import annotations

import argparse
import sys

import networkx as nx
import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, FS_OVERRIDES, PATIENTS_LIST, nperseg_for_fs
from lrg_eegfc.config.paths import CACHE_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrgsglib.core import compute_laplacian_properties, get_giant_component

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "hypothesis_tests"))
from _fc_split_half import compute_imcoh_abs_halves  # noqa: E402

HALVES_CACHE = CACHE_ROOT / "imcoh_lrg_halves"
OUT_DIR = ROOT / "data" / "audit" / "eigenvector_cache_check"
OUT_CSV = OUT_DIR / "audit_26_eigenvector_cache_n10_imcoh_abs.csv"

ORTHO_TOL = 1e-10
TRIVIAL_TOL = 1e-8
RECON_TOL = 1e-10

FULL_PHASES = ("rest_pre", "task_learn", "task_test", "rest_post")
HALF_REST_PHASES = ("rest_pre", "rest_post")


def _check_one(L, eigvals, eigvecs, n_nodes):
    L = np.asarray(L)
    if eigvals is None or eigvecs is None:
        return dict(
            ortho_err=np.nan,
            trivial_eigval=np.nan,
            recon_err=np.nan,
            pass_ortho=False,
            pass_trivial=False,
            pass_recon=False,
            pass_all=False,
            note="no_eigvecs_in_npz",
        )
    ortho_err = float(np.linalg.norm(eigvecs.T @ eigvecs - np.eye(n_nodes)))
    trivial = float(abs(eigvals[0]))
    L_recon = eigvecs @ np.diag(eigvals) @ eigvecs.T
    recon_err = float(np.linalg.norm(L_recon - L) / max(np.linalg.norm(L), 1e-30))
    p_ortho = ortho_err < ORTHO_TOL
    p_trivial = trivial < TRIVIAL_TOL
    p_recon = recon_err < RECON_TOL
    return dict(
        ortho_err=ortho_err,
        trivial_eigval=trivial,
        recon_err=recon_err,
        pass_ortho=p_ortho,
        pass_trivial=p_trivial,
        pass_recon=p_recon,
        pass_all=p_ortho and p_trivial and p_recon,
        note="",
    )


def _laplacian_from_A(A: np.ndarray):
    g = nx.from_numpy_array(np.asarray(A))
    giant = get_giant_component(g)
    _, L, *_ = compute_laplacian_properties(giant)
    return np.asarray(L), giant.number_of_nodes()


def check_full_phase(pat: str, band: str, phase: str) -> dict:
    res = load_lrg_result(pat, phase, band, "imcoh_abs")
    if res is None:
        return dict(patient=pat, band=band, phase=phase, source="full",
                    pass_all=False, note="no_lrg_npz")
    A = load_fc_matrix(pat, phase, band, "imcoh_abs")
    if A is None:
        return dict(patient=pat, band=band, phase=phase, source="full",
                    pass_all=False, note="no_fc_matrix")
    L, n = _laplacian_from_A(A)
    out = _check_one(L, res.eigenvalues, res.eigenvectors, n)
    return dict(patient=pat, band=band, phase=phase, source="full",
                n_nodes=n, **out)


def check_halves(pat: str, verbose: bool = False) -> list[dict]:
    rows: list[dict] = []
    for phase in HALF_REST_PHASES:
        try:
            X = load_timeseries(pat, phase, SEEG_DATAPATH)
        except (FileNotFoundError, OSError):
            for band in BRAIN_BANDS_NAMES:
                for tag in ("A", "B"):
                    rows.append(dict(patient=pat, band=band,
                                     phase=f"{phase}_{tag}", source="half",
                                     pass_all=False, note="timeseries_missing"))
            continue
        if X is None:
            continue
        X = np.asarray(X, dtype=np.float64)
        if X.shape[0] > X.shape[1]:
            X = X.T
        fs = FS_OVERRIDES.get(pat, 2048.0)
        nperseg_half = max(256, nperseg_for_fs(fs) // 2)
        from lrg_eegfc.config.const import BRAIN_BANDS
        halves_fc = compute_imcoh_abs_halves(X, fs, nperseg_half, BRAIN_BANDS)
        for (band, tag), A in halves_fc.items():
            synth_phase = f"{phase}_{tag}"
            res = load_lrg_result(pat, synth_phase, band, "imcoh_abs",
                                  cache_root=HALVES_CACHE)
            if res is None:
                rows.append(dict(patient=pat, band=band, phase=synth_phase,
                                 source="half", pass_all=False,
                                 note="no_lrg_npz"))
                continue
            L, n = _laplacian_from_A(A)
            out = _check_one(L, res.eigenvalues, res.eigenvectors, n)
            rows.append(dict(patient=pat, band=band, phase=synth_phase,
                             source="half", n_nodes=n, **out))
        if verbose:
            print(f"  {pat}/{phase}: {len(halves_fc)} halves checked")
    return rows


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("--no-halves", action="store_true",
                   help="Skip halves check (full-phase only).")
    args = p.parse_args(argv)

    rows: list[dict] = []
    for pat in PATIENTS_LIST:
        if args.verbose:
            print(f"=== {pat} ===")
        for band in BRAIN_BANDS_NAMES:
            for phase in FULL_PHASES:
                rows.append(check_full_phase(pat, band, phase))
        if not args.no_halves:
            rows.extend(check_halves(pat, verbose=args.verbose))

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)

    n_total = len(df)
    n_pass = int(df["pass_all"].sum())
    n_fail = n_total - n_pass
    print(f"\nWrote {OUT_CSV}")
    print(f"  {n_pass}/{n_total} pass spectral identity")
    if n_fail:
        print(f"  {n_fail} FAIL — see CSV note column for details")
        bad = df[~df["pass_all"]]
        print(bad[["patient", "band", "phase", "source", "note"]].head(20).to_string(index=False))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
