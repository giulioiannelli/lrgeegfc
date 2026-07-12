#!/usr/bin/env python3
"""D0 -- percolation thresholding protocol validation (P_inf / E_inf vs theta).

FOUNDATION of the sparsified arc. Before any diffusion/trace result we must show
the parameter-free percolation backbone behaves as designed for EVERY
(patient, band, phase):
  - the giant component spans all N nodes (P_inf = 1) right up to the bottleneck
    theta*, then erodes -- so the backbone at theta* is CONNECTED;
  - the surviving graph at theta* is CYCLE-RICH (cycle rank >> 0, non-zero
    transitivity), NOT a tree;
  - theta* / surviving-edge-fraction are STABLE across the cohort (no outlier
    cell that collapses to a near-tree or stays near-fully-connected).

theta* = weakest maximum-spanning-tree edge = largest theta with P_inf = 1
(percolation_backbone). Density is SET BY THE DATA, never chosen.

Outputs (cache; plotting is a separate step -- fig_percolation_protocol.py):
    data/sparsified_arc/percolation/{band}/{patient}_{phase}.npz
        theta, p_inf, e_inf, edge_frac, theta_star, edge_frac_star, w_max, ...
    data/sparsified_arc/percolation/summary.csv      (one row per cell, scalars)
    data/sparsified_arc/percolation/config.json      (provenance)
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
from audit_150_rho_sym_gate import load_phase, COHORT, BANDS  # dataset-specific loader
from lrg_eegfc.utils.fc.backbone import (
    percolation_sweep, percolation_backbone, backbone_density,
)
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig

PHASES = ("A", "B", "task_test", "rest_post")
OUT = ROOT / "data" / "sparsified_arc" / "percolation"


def _transitivity(B: np.ndarray) -> float:
    A = (B > 0).astype(float)
    A2 = A @ A
    tri = float(np.trace(A2 @ A)) / 6.0
    k = A.sum(1)
    triples = float(np.sum(k * (k - 1))) / 2.0
    return 3.0 * tri / triples if triples > 0 else 0.0


def run_cell(pat, phase, band):
    W = load_phase(pat, phase, band)
    N = W.shape[0]
    sweep = percolation_sweep(W)
    B, theta = percolation_backbone(W)
    ev, _ = laplacian_eig(B)
    n_edges_bb = int((np.triu(B, 1) > 0).sum())
    cycle_rank = n_edges_bb - N + 1                       # E - N + C, C=1 (connected)
    deg = (B > 0).sum(1)
    rec = dict(
        patient=pat, phase=phase, band=band, N_nodes=int(N),
        w_max=float(sweep["w_max"]),
        theta_star=float(sweep["theta_star"]),
        theta_star_frac=float(sweep["theta_star"] / sweep["w_max"]) if sweep["w_max"] > 0 else np.nan,
        edge_frac_star=float(sweep["edge_frac_star"]),
        p_inf_star=float(sweep["p_inf_star"]),
        density=float(backbone_density(B)),
        n_edges_backbone=n_edges_bb,
        cycle_rank=int(cycle_rank),
        transitivity=float(_transitivity(B)),
        mean_degree=float(deg.mean()), min_degree=int(deg.min()), max_degree=int(deg.max()),
        lambda_2=float(ev[1]) if N > 1 else np.nan,        # >0 confirms connected
        connected=bool(ev[1] > 1e-9) if N > 1 else False,
    )
    # cache the full percolation curve for the D0 figure
    cell_dir = OUT / band
    cell_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cell_dir / f"{pat}_{phase}.npz",
        theta=sweep["theta"], p_inf=sweep["p_inf"], e_inf=sweep["e_inf"],
        edge_frac=sweep["edge_frac"], theta_star=sweep["theta_star"],
        edge_frac_star=sweep["edge_frac_star"], w_max=sweep["w_max"],
        N_nodes=N, cycle_rank=cycle_rank, density=rec["density"],
    )
    return rec


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(p, ph, b) for b in BANDS for p in COHORT for ph in PHASES]
    print(f"[D0-percolation] {len(jobs)} cells (patient x phase x band)", flush=True)
    t0 = time.time()
    rows = []
    for i, (p, ph, b) in enumerate(jobs, 1):
        try:
            rec = run_cell(p, ph, b)
        except Exception as e:
            print(f"[{i}/{len(jobs)}] SKIP {p}/{ph}/{b}: {e}", flush=True)
            continue
        rows.append(rec)
        if i % 20 == 0 or i == len(jobs):
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {p}/{ph}/{b} "
                  f"theta*={rec['theta_star']:.3f} ({100*rec['theta_star_frac']:.0f}% of wmax) "
                  f"edge_frac*={100*rec['edge_frac_star']:.1f}% dens={rec['density']:.3f} "
                  f"cyc={rec['cycle_rank']} trans={rec['transitivity']:.2f} lam2={rec['lambda_2']:.1e} "
                  f"[{el:.0f}s]", flush=True)
    df = pd.DataFrame(rows).sort_values(["band", "patient", "phase"])
    df.to_csv(OUT / "summary.csv", index=False)

    # stability / sanity roll-up per band
    print("\n=== PERCOLATION PROTOCOL STABILITY (per band, over patients x phases) ===", flush=True)
    print(f"{'band':<11} {'edge_frac* [%]':>22} {'density':>16} {'cycle_rank':>16} "
          f"{'trans':>12} {'connected':>10}", flush=True)
    for b in BANDS:
        x = df[df.band == b]
        if x.empty:
            continue
        ef = 100 * x.edge_frac_star
        de = x.density
        cy = x.cycle_rank
        tr = x.transitivity
        print(f"{b:<11} {ef.median():6.1f} [{ef.min():5.1f},{ef.max():5.1f}] "
              f"{de.median():7.3f} [{de.min():.3f},{de.max():.3f}] "
              f"{cy.median():7.0f} [{cy.min():.0f},{cy.max():.0f}] "
              f"{tr.median():6.2f} [{tr.min():.2f},{tr.max():.2f}] "
              f"{int(x.connected.all())}/{'all' if x.connected.all() else 'SOME-FAIL'}",
              flush=True)

    n_tree = int((df.cycle_rank <= 0).sum())
    n_disc = int((~df.connected).sum())
    print(f"\n[flags] tree-like cells (cycle_rank<=0): {n_tree}/{len(df)} ; "
          f"disconnected cells: {n_disc}/{len(df)}", flush=True)

    config = dict(
        deliverable="D0_percolation_protocol",
        rule="keep edges >= theta* ; theta* = weakest max-spanning-tree edge = "
             "largest theta with P_inf==1 (parameter-free, data-driven density)",
        cohort=COHORT, bands=BANDS, phases=list(PHASES),
        fc_method="imcoh_abs", n_cells=len(df),
        tree_like_cells=n_tree, disconnected_cells=n_disc,
        note="A,B = split-half of rest_pre. theta* computed per (patient,phase,band).",
    )
    (OUT / "config.json").write_text(json.dumps(config, indent=2))
    print(f"\n[D0-percolation] {len(df)} cells in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
