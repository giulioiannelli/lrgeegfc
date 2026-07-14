#!/usr/bin/env python3
r"""26_diffusion_spatial_reach -- link the diffusion (communication) scale to the
system's PHYSICAL scale, from the contact pitch to the implant span.

Descriptive geometry, NOT a hypothesis test -- no null, no p-value. As we run the LRG
diffusion on the mst@0.20 backbone from fine to coarse, what PHYSICAL length does it
connect? We measure the physical extent (mm) the process has spanned at diffusion time tau
-- an actual physical length, NOT a heat-weighted average displacement.

Measure:
  B    = mst_union_top_fraction(W_restpost, 0.20)                   # connected backbone
  K    = e^{-tau L(B)} = V exp(-tau*ev) V^T, tau = s/lambda_max      # row-stochastic (L1=0)
  ell(s) = max{ ||x_i - x_j|| : K_ij(s) >= 1/N }   [mm]
           = the largest physical distance between two contacts that COMMUNICATE (the heat
             from one has reached at least its uniform-equilibrium share 1/N at the other).
  ell(s) grows monotonically from the finest connected functional cluster to the span, which
  it reaches when the two farthest contacts equilibrate. Reported against two physical refs:
    pitch = median_i min_{j!=i} ||x_i-x_j||   (~3.5 mm)  # electrode spacing (a single edge)
    span  = max_{i,j} ||x_i-x_j||   (~9-11 cm)           # the two farthest contacts = whole space

Two facts the curve makes concrete (both HONEST, one corrects an earlier draft):
  1. The finest connected functional length is ~15 mm, NOT the 3.5 mm pitch. The |ImCoh|
     backbone is FUNCTIONAL, not spatial: only ~2-4% of its edges are same-shaft (<=5 mm),
     median edge ~35 mm. Diffusion has almost no short edges to ride, so it never resolves
     the electrode pitch -- the finest functional cluster is already ~15 mm (delocalised).
  2. s=1 <=> tau=1/lambda_max is the e-folding time of the fastest Laplacian mode -- the
     numerically-principled finest diffusion time, and the raw<->communication boundary:
     for s<1, e^{-tau L} ~ I - tau L (off-diagonal ~ tau*W_ij) reads single raw edges only;
     multi-hop communication opens for s>1. So the anchor hides no finer *communication*
     scale (below it is raw FC = the baseline). Marked on the figure.

Output: data/sparsified_arc/spatial_reach_mst020/{per_cell,cohort}.csv
"""
from __future__ import annotations
import sys
import time

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, pdist

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, BANDS, load_phase          # noqa: E402
from lrg_eegfc.utils.fc.backbone import mst_union_top_fraction        # noqa: E402
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig          # noqa: E402
from lrg_eegfc.utils.io.regions import load_channel_regions           # noqa: E402

FRAC = 0.20
PHASE = "rest_post"                                   # phase Sec-1's trace is read on
# extended grid spanning the raw/single-edge regime (s<1) through the communication window
# up to the coarse end; s=1 = tau_min = 1/lambda_max is the raw<->communication boundary.
SGRID = np.logspace(np.log10(0.05), np.log10(180.0), 24)   # s = tau*lambda_max
OUT = ROOT / "data" / "sparsified_arc" / "spatial_reach_mst020"


def phys_coords_mm(pat: str):
    """(N,3) electrode coords in mm, FC-node order, + finite-coord mask."""
    rdf = load_channel_regions(pat)
    xyz = rdf[["x", "y", "z"]].to_numpy(float) / 1000.0      # micrometres -> mm
    mask = np.isfinite(xyz).all(axis=1)
    return xyz, mask


def connect_curve(ev, V, D, mask):
    """ell(s) [mm]: the physical EXTENT the diffusion has connected at scale s.

    ell(s) = max{ ||x_i - x_j|| : K_ij(s) >= 1/N } -- the largest physical distance
    between two contacts that COMMUNICATE, i.e. the heat released at one has reached at
    least its uniform-equilibrium share (1/N) at the other. This is a physical length in
    mm (NOT a heat-weighted average): it rises from the finest functional cluster (the
    shortest strong |ImCoh| edges, ~15 mm -- NOT the electrode pitch, since the backbone
    has almost no same-shaft edges) to the implant span (reached when the two farthest
    contacts equilibrate, s -> inf). ``D`` is the (Nc,Nc) distance matrix among coord-known
    nodes; the kernel runs on the FULL backbone, read out on placeable contacts.
    """
    lam_max = ev[-1]
    thr = 1.0 / V.shape[0]                                    # uniform-equilibrium share
    ell = np.full(len(SGRID), np.nan)                         # NaN = nothing connected yet (NOT 0 mm!)
    for k, s in enumerate(SGRID):
        K = (V * np.exp(-(s / lam_max) * ev)) @ V.T          # row-stochastic heat kernel
        comm = K[np.ix_(mask, mask)] >= thr                  # pairs that have equilibrated
        np.fill_diagonal(comm, False)                        # exclude self (d=0)
        if comm.any():
            ell[k] = float(D[comm].max())                    # physical diameter of connected set
    return ell


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    t0 = time.time()
    n_pat = len(COHORT)
    for ip, pat in enumerate(COHORT, 1):
        tp = time.time()
        xyz, mask = phys_coords_mm(pat)
        p = xyz[mask]
        D = cdist(p, p)                                   # (Nc,Nc) physical distances [mm]
        pdv = pdist(p)                                    # condensed pairwise distances
        span_mm = float(pdv.max())                        # true implant extent (max pairwise)
        Dnn = D.copy()
        np.fill_diagonal(Dnn, np.inf)
        pitch_mm = float(np.median(Dnn.min(1)))           # median nearest-contact spacing
        for band in BANDS:
            W = load_phase(pat, PHASE, band)
            if W.shape[0] != xyz.shape[0]:
                print(f"  [skip] {pat} {band}: FC N={W.shape[0]} != coords "
                      f"{xyz.shape[0]}", flush=True)
                continue
            ev, V = laplacian_eig(mst_union_top_fraction(W, FRAC))
            ell = connect_curve(ev, V, D, mask)
            for k, s in enumerate(SGRID):
                rows.append(dict(patient=pat, band=band, s=s, i_s=k,
                                 ell_mm=ell[k], pitch_mm=pitch_mm, span_mm=span_mm,
                                 n_coord=int(mask.sum()), N=int(xyz.shape[0])))
        el = time.time() - t0
        eta = el / ip * (n_pat - ip)
        print(f"[{ip}/{n_pat}] {pat}  n_coord={int(mask.sum())}/{xyz.shape[0]}  "
              f"({time.time()-tp:.1f}s)  elapsed={el:.1f}s  ETA={eta:.1f}s", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "per_cell.csv", index=False)

    # cohort curves: median + IQR over patients that have CONNECTED (ell finite); NaN/0
    # dead-zone patients are excluded, not averaged in as "0 mm". n_conn = #patients connected.
    g = df.groupby(["band", "i_s", "s"])
    coh = g.agg(ell_med=("ell_mm", lambda a: np.nanmedian(a) if np.isfinite(a).any() else np.nan),
                ell_lo=("ell_mm", lambda a: np.nanpercentile(a, 25) if np.isfinite(a).any() else np.nan),
                ell_hi=("ell_mm", lambda a: np.nanpercentile(a, 75) if np.isfinite(a).any() else np.nan),
                n_conn=("ell_mm", lambda a: int(np.isfinite(a).sum())),
                n=("ell_mm", "size")).reset_index()
    coh.to_csv(OUT / "cohort.csv", index=False)

    # headline: the connected physical extent ell(s) runs from the finest functional
    # cluster (shortest strong |ImCoh| edges, ~15 mm -- NOT the ~3.5 mm pitch) to the span.
    print(f"\n[done] {time.time()-t0:.1f}s  ->  {OUT}", flush=True)
    pitch = df["pitch_mm"].median()
    span = df["span_mm"].median()
    i1 = int(np.argmin(np.abs(SGRID - 1.0)))                    # nearest-grid index to s=1
    print(f"  physical endpoints (cohort median): pitch~{pitch:.1f}mm (electrode spacing, a "
          f"single edge)  |  implant span~{span:.0f}mm (two farthest contacts)", flush=True)
    for band in BANDS:
        sub = df[df.band == band]
        if sub.empty:
            continue
        e1 = sub[sub.i_s == i1]["ell_mm"].median()             # connected extent at s=1 (skips NaN)
        emax = sub[sub.i_s == len(SGRID) - 1]["ell_mm"].median()
        # finest cohort-level connected length (>= 5 patients connected)
        cb = coh[(coh.band == band) & (coh.n_conn >= 5)].sort_values("s")
        efine = cb["ell_med"].iloc[0] if len(cb) else np.nan
        print(f"  {band:>10s}: finest cohort connected~{efine:.0f}mm ({efine/pitch:.0f}x pitch)  "
              f"ell(s=1)~{e1:.0f}mm  ell(coarse)~{emax:.0f}mm (~{emax/span*100:.0f}% of span)",
              flush=True)


if __name__ == "__main__":
    main()
