#!/usr/bin/env python3
"""Audit Step 18 — same outputs as audit_17 but the per-phase dendrogram
is rebuilt at the LRG characteristic timescale ``τ* = argmax_τ C(τ)``
(peak of spectral complexity), instead of the shortest scale 1/λ_max.

Reuses every plot helper from ``audit_17_taumin_four_phase`` (passing a
phase-specific ``out_dir``); the only new piece is constructing an
LRGResult-shaped object whose linkage was computed at τ*.

Outputs
-------
data/audit/four_phase_tau_star/<Patient>/<band>_dendrograms.pdf
data/audit/four_phase_tau_star/<Patient>/<band>_psi.pdf
data/audit/four_phase_tau_star/<Patient>/<band>_dendrogram_restpostref_k={K}.pdf
data/audit/four_phase_tau_star/ledger.csv
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional

import networkx as nx
import numpy as np
from scipy.spatial.distance import squareform

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

# Ensure sibling script importable regardless of cwd
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

FC_METHOD = "imcoh_abs"
OUT_DIR = DATA_ROOT / "audit" / "four_phase_tau_star"
LEDGER_PATH = OUT_DIR / "ledger.csv"


class _LRGAtTau:
    """Minimal LRGResult-like shim consumed by audit_17 plot helpers."""
    __slots__ = ("linkage_matrix", "n_nodes", "optimal_threshold",
                 "ultrametric_matrix", "patient", "phase", "band",
                 "fc_method")

    def __init__(self, linkage_matrix, n_nodes, optimal_threshold,
                 ultrametric_matrix, patient, phase, band, fc_method):
        self.linkage_matrix = linkage_matrix
        self.n_nodes = int(n_nodes)
        self.optimal_threshold = float(optimal_threshold)
        self.ultrametric_matrix = ultrametric_matrix
        self.patient = patient
        self.phase = phase
        self.band = band
        self.fc_method = fc_method


def _build_at_tau_star(pat: str, phase: str, band: str):
    """Return ``(_LRGAtTau, tau_star)`` for the given cell, or ``(None, nan)``
    if either FC or cached LRG (needed for τ*) is missing."""
    cached = load_lrg_result(pat, phase, band, FC_METHOD)
    if cached is None:
        return None, float("nan")
    fc = load_fc_matrix(pat, phase, band, FC_METHOD)
    if fc is None:
        return None, float("nan")

    A = np.abs(fc).astype(float)
    np.fill_diagonal(A, 0.0)
    G = get_giant_component(nx.from_numpy_array(A))

    tau_star = float(cached.entropy_tau[int(np.argmax(cached.entropy_C))])
    _, _, _, Trho, _ = compute_laplacian_properties(G, tau=tau_star)
    Trho = np.asarray(Trho)
    dists = squareform(Trho, checks=False)
    Z, _, _ = compute_normalized_linkage(dists, G)
    n_nodes = G.number_of_nodes()
    shim = _LRGAtTau(
        linkage_matrix=Z,
        n_nodes=n_nodes,
        optimal_threshold=cached.optimal_threshold,
        ultrametric_matrix=Trho,
        patient=pat, phase=phase, band=band, fc_method=FC_METHOD,
    )
    return shim, tau_star


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows: List[dict] = []
    n_dend = n_psi = n_skip = 0

    for pat in PATIENTS_4PHASE:
        labels = _load_channel_labels(pat, SEEG_DATAPATH) or []
        print(f"{pat} (channels={len(labels)})", flush=True)

        for band in BRAIN_BANDS_NAMES:
            results: Dict[str, object] = {}
            tau_star_per_phase: Dict[str, float] = {}
            for phase in PHASE_LABELS:
                try:
                    r, tau_star = _build_at_tau_star(pat, phase, band)
                except Exception as exc:  # noqa: BLE001
                    print(f"  [error] {band}/{phase}: {exc}")
                    r, tau_star = None, float("nan")
                results[phase] = r
                tau_star_per_phase[phase] = tau_star
                if r is None:
                    n_skip += 1
                    continue
                psi_values, n_comms = compute_partition_stability_index(
                    r.linkage_matrix
                )
                psi_max = float(psi_values.max()) if psi_values.size else float("nan")
                argmax_n = (
                    int(n_comms[np.argmax(psi_values)]) if psi_values.size else -1
                )
                rows.append(dict(
                    patient=pat,
                    band=band,
                    phase=phase,
                    n_nodes=int(r.n_nodes),
                    tau_star=tau_star,
                    psi_argmax_n=argmax_n,
                    psi_max=psi_max,
                    optimal_threshold=float(r.optimal_threshold),
                ))

            if any(v is not None for v in results.values()):
                dp = A17._plot_dendrograms(pat, band, results, labels, out_dir=OUT_DIR)
                pp = A17._plot_psi(pat, band, results, out_dir=OUT_DIR)
                ref_paths: List[Path] = []
                for k in A17.K_REF_LIST:
                    rp = A17._plot_dendrograms_phaseref(
                        pat, band, results, k, A17.REF_PHASE, out_dir=OUT_DIR
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
                print(f"  {band}: skipped (no LRG cache for any phase)")

    with LEDGER_PATH.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "patient", "band", "phase", "n_nodes", "tau_star",
            "psi_argmax_n", "psi_max", "optimal_threshold",
        ])
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {LEDGER_PATH}  "
          f"({len(rows)} rows, {n_dend} dend pdfs, {n_psi} psi pdfs, "
          f"{n_skip} cells skipped)")


if __name__ == "__main__":
    main()
