#!/usr/bin/env python3
r"""PMFG-observed cophenetic-trace readout -- is alpha-loss a PLANAR-family property?

The TMFG (greedy chordal) matched-strength trace gate LOST alpha (0/16 scales,
was 12/16 on mst@0.20) and collapsed beta (16->3/16), while strength-based
sparsifiers (mst-family, disparity) KEEP alpha. That points the finger at
PLANARITY, not sparsity. But TMFG is only an *approximation* of the exact planar
filter (PMFG), and TMFG does not even contain the MST. So a reviewer can object:
"you disqualified the planar family on the greedy shortcut; the exact PMFG might
keep alpha."

This script settles that. It computes the OBSERVED cophenetic trace rho_sym per
(band, scale) on the EXACT PMFG backbone (networkx planarity, ~9 s/graph -> too
slow for a per-surrogate null, so observed-only), apples-to-apples with the same
mst@0.20 / TMFG exploratory readout (script 12). No matched-strength null here:
the exploratory rho_sym MAGNITUDE already discriminates cleanly (script 12:
alpha rho_sym mst@0.05=0.237 vs TMFG=0.116 at equal density). If PMFG's alpha
rho_sym sits with TMFG (~0.12) not with the strength-ranked mst-family (~0.24),
the whole PLANAR family -- exact included -- degrades the trace, and the
disqualification is not a greedy artifact.

Reuses script-12 primitives (load_phase, cophenetic_at_scale, rho_sym,
laplacian_eig) + lrg_eegfc.utils.fc.backbone. Parallel over cells.

Outputs (data/sparsified_arc/planar_trace_readout/):
  pmfg_rho_sym.csv    : patient, band, backbone, s, rho_sym  (backbone in {pmfg,tmfg,mst0.05,mst0.20})
  cohort_best.csv     : band x backbone -> best-scale cohort-median rho_sym + s
"""
from __future__ import annotations
import json, os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, BANDS, load_phase
from lrg_eegfc.utils.fc.backbone import (
    pmfg_backbone, tmfg_backbone, mst_union_top_fraction, backbone_density,
)
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, cophenetic_at_scale, rho_sym

OUT = ROOT / "data" / "sparsified_arc" / "planar_trace_readout"
PHASES = ("A", "B", "task_test", "rest_post")
SGRID = np.logspace(0.0, np.log10(180.0), 16)      # identical to script 12/13
BACKBONES = {
    "pmfg":    lambda W: pmfg_backbone(W),
    "tmfg":    lambda W: tmfg_backbone(W),          # re-confirm apples-to-apples
    "mst0.05": lambda W: mst_union_top_fraction(W, 0.05),   # strength-ranked, PMFG density
    "mst0.20": lambda W: mst_union_top_fraction(W, 0.20),   # current trace backbone
}


def per_cell(job):
    pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES}
    except Exception:
        return []
    rows = []
    for name, fn in BACKBONES.items():
        try:
            eig = {ph: laplacian_eig(fn(Ws[ph])) for ph in PHASES}
            dens = float(np.mean([backbone_density(fn(Ws[ph])) for ph in PHASES]))
        except Exception:
            continue
        for s in SGRID:
            try:
                D = {ph: cophenetic_at_scale(*eig[ph], float(s)) for ph in PHASES}
                if len({d.size for d in D.values()}) != 1:
                    raise ValueError("size")
                r, _ = rho_sym(D["A"], D["B"], D["task_test"], D["rest_post"])
            except Exception:
                r = np.nan
            rows.append(dict(patient=pat, band=band, backbone=name, density=dens,
                             s=float(s), rho_sym=float(r)))
    return rows


def cohort_best(df):
    out = []
    for band in BANDS:
        for bb in df.backbone.unique():
            x = df[(df.band == band) & (df.backbone == bb)]
            if x.empty:
                continue
            g = x.groupby("s").rho_sym.median()
            g = g[np.isfinite(g.values)]
            if g.empty:
                continue
            out.append(dict(band=band, backbone=bb, density=float(x.density.median()),
                            rho_best=float(g.max()), s_best=float(g.idxmax()),
                            n_pat=int(x.patient.nunique())))
    return pd.DataFrame(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(p, b) for b in BANDS for p in COHORT]
    ncpu = int(os.environ.get("SA_WORKERS", 12))
    print(f"[planar-readout] {len(jobs)} cells x {len(BACKBONES)} backbones "
          f"(pmfg ~9s/graph) | {ncpu} workers", flush=True)
    t0 = time.time(); rows = []
    with Pool(ncpu) as pool:
        for i, rl in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            rows.extend(rl)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {el:.0f}s ETA {el/i*(len(jobs)-i):.0f}s", flush=True)
    df = pd.DataFrame(rows)
    if df.empty:
        print("[planar-readout] no cells", flush=True); return
    df.to_csv(OUT / "pmfg_rho_sym.csv", index=False)
    cb = cohort_best(df); cb.to_csv(OUT / "cohort_best.csv", index=False)
    piv = cb.pivot(index="band", columns="backbone", values="rho_best").reindex(BANDS)
    cols = [c for c in ["mst0.20", "mst0.05", "pmfg", "tmfg"] if c in piv.columns]
    print(f"\n[planar-readout] {len(df)} rows in {time.time()-t0:.0f}s -> {OUT}\n", flush=True)
    print("=== best-scale cohort-median rho_sym (band x backbone) ===", flush=True)
    print("    density:", {c: round(float(cb[cb.backbone==c].density.median()),3) for c in cols}, flush=True)
    print(piv[cols].round(3).to_string(), flush=True)
    print("\n  KEY: if PMFG's alpha ~ TMFG (planar, low) not ~ mst0.05 (strength, high),"
          "\n       the EXACT planar filter also loses alpha -> planarity, not greedy shortcut.", flush=True)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="pmfg_observed_trace_readout", phases=list(PHASES),
        s_grid=[float(x) for x in SGRID], cohort=COHORT, bands=BANDS,
        backbones=list(BACKBONES), note="observed rho_sym magnitude, no MS null (PMFG too slow for null)"),
        indent=2))
    print(f"[planar-readout] done -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
