#!/usr/bin/env python3
"""Audit Step 19 — same outputs as audit_17 but the per-phase dendrogram
is rebuilt at τ = 1/λ_idx where λ_idx is the *idx*-th eigenvalue of the
giant-component Laplacian counted from the largest (idx = 0 → λ_max,
idx = 1 → second-largest, …).

The zero eigenvalue (trivial mode of a connected graph) is excluded; the
last available index is therefore ``n_nonzero_eigvals − 1`` and is the
algebraic-connectivity inverse (longest LRG timescale).

Reuses every plot helper from ``audit_17_taumin_four_phase``.

Usage
-----
    python audit_19_lambda_idx_four_phase.py --lambda-idx N

Outputs
-------
data/audit/four_phase_lambda_idx{N}/<Patient>/<band>_dendrograms.pdf
data/audit/four_phase_lambda_idx{N}/<Patient>/<band>_psi.pdf
data/audit/four_phase_lambda_idx{N}/<Patient>/<band>_dendrogram_restpostref_k={K}.pdf
data/audit/four_phase_lambda_idx{N}/ledger.csv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List

import networkx as nx
import numpy as np
from scipy.spatial.distance import squareform

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from lrg_eegfc.config.const import (
    BRAIN_BANDS_NAMES,
    PATIENTS_4PHASE,
    PHASE_LABELS,
)
from lrg_eegfc.config.paths import DATA_ROOT, SEEG_DATAPATH
from lrg_eegfc.visuals.lrg import (
    _load_channel_labels,
    compute_partition_stability_index,
)
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrgsglib.nx_patches.funcs import get_giant_component
from lrgsglib.utils.lrg.clustering import compute_normalized_linkage
from lrgsglib.utils.lrg.spectral import compute_laplacian_properties

import audit_17_taumin_four_phase as A17
from audit_18_taustar_four_phase import _LRGAtTau

FC_METHOD = "imcoh_abs"


def _build_at_lambda_idx(pat: str, phase: str, band: str, lambda_idx: int):
    """Return ``(_LRGAtTau, tau, lambda_value)`` for the given cell, or
    ``(None, nan, nan)`` if the FC or cached LRG (used only for
    ``optimal_threshold``) is missing, or if ``lambda_idx`` is out of range
    for this graph's spectrum."""
    cached = load_lrg_result(pat, phase, band, FC_METHOD)
    fc = load_fc_matrix(pat, phase, band, FC_METHOD)
    if fc is None:
        return None, float("nan"), float("nan")

    A = np.abs(fc).astype(float)
    np.fill_diagonal(A, 0.0)
    G = get_giant_component(nx.from_numpy_array(A))

    spectrum = nx.laplacian_spectrum(G)
    spectrum = np.asarray(spectrum, dtype=float)
    nonzero = spectrum[spectrum > spectrum.max() * 1e-12]
    nonzero_desc = np.sort(nonzero)[::-1]
    if lambda_idx >= nonzero_desc.size or lambda_idx < 0:
        return None, float("nan"), float("nan")
    lam = float(nonzero_desc[lambda_idx])
    tau = 1.0 / lam

    _, _, _, Trho, _ = compute_laplacian_properties(G, tau=tau)
    Trho = np.asarray(Trho)
    dists = squareform(Trho, checks=False)
    Z, _, _ = compute_normalized_linkage(dists, G)
    n_nodes = G.number_of_nodes()
    fallback = float(cached.optimal_threshold) if cached is not None else float(np.median(Z[:, 2]))
    shim = _LRGAtTau(
        linkage_matrix=Z,
        n_nodes=n_nodes,
        optimal_threshold=fallback,
        ultrametric_matrix=Trho,
        patient=pat, phase=phase, band=band, fc_method=FC_METHOD,
    )
    return shim, tau, lam


def main(lambda_idx: int) -> None:
    out_dir = DATA_ROOT / "audit" / f"four_phase_lambda_idx{lambda_idx}"
    out_dir.mkdir(parents=True, exist_ok=True)
    ledger_path = out_dir / "ledger.csv"

    rows: List[dict] = []
    n_dend = n_psi = n_skip = 0

    for pat in PATIENTS_4PHASE:
        labels = _load_channel_labels(pat, SEEG_DATAPATH) or []
        print(f"{pat} (channels={len(labels)})  λ_idx={lambda_idx}", flush=True)

        for band in BRAIN_BANDS_NAMES:
            results: Dict[str, object] = {}
            tau_per_phase: Dict[str, float] = {}
            lam_per_phase: Dict[str, float] = {}
            for phase in PHASE_LABELS:
                try:
                    r, tau, lam = _build_at_lambda_idx(pat, phase, band, lambda_idx)
                except Exception as exc:  # noqa: BLE001
                    print(f"  [error] {band}/{phase}: {exc}")
                    r, tau, lam = None, float("nan"), float("nan")
                results[phase] = r
                tau_per_phase[phase] = tau
                lam_per_phase[phase] = lam
                if r is None:
                    n_skip += 1
                    continue
                psi_v, psi_nc = compute_partition_stability_index(r.linkage_matrix)
                psi_max = float(psi_v.max()) if psi_v.size else float("nan")
                argmax_n = int(psi_nc[np.argmax(psi_v)]) if psi_v.size else -1
                rows.append(dict(
                    patient=pat,
                    band=band,
                    phase=phase,
                    n_nodes=int(r.n_nodes),
                    lambda_idx=lambda_idx,
                    lambda_value=lam,
                    tau=tau,
                    psi_argmax_n=argmax_n,
                    psi_max=psi_max,
                    optimal_threshold=float(r.optimal_threshold),
                ))

            if any(v is not None for v in results.values()):
                dp = A17._plot_dendrograms(pat, band, results, labels, out_dir=out_dir)
                pp = A17._plot_psi(pat, band, results, out_dir=out_dir)
                ref_paths: List[Path] = []
                for k in A17.K_REF_LIST:
                    rp = A17._plot_dendrograms_phaseref(
                        pat, band, results, k, A17.REF_PHASE, out_dir=out_dir
                    )
                    if rp is not None:
                        ref_paths.append(rp)
                n_dend += 1
                n_psi += 1
                msg = f"  {band}: wrote {dp} + {pp}"
                for rp in ref_paths:
                    msg += f" + {rp}"
                print(msg)
            else:
                print(f"  {band}: skipped (no LRG cache or λ_idx out of range)")

    with ledger_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "patient", "band", "phase", "n_nodes",
            "lambda_idx", "lambda_value", "tau",
            "psi_argmax_n", "psi_max", "optimal_threshold",
        ])
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {ledger_path}  "
          f"({len(rows)} rows, {n_dend} dend pdfs, {n_psi} psi pdfs, "
          f"{n_skip} cells skipped)")


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Per-(patient, band) τ=1/λ_idx four-phase dendrograms."
    )
    p.add_argument(
        "--lambda-idx", type=int, default=0,
        help="0 = largest non-zero λ (= τ_min, default). 1 = second largest. "
             "Last index is the algebraic-connectivity inverse (= τ_max).",
    )
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    main(lambda_idx=args.lambda_idx)
