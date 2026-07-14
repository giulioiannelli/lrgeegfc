#!/usr/bin/env python3
"""Two settled controls on mst@0.20, from cached arcs (no FC/LRG recompute).

(1) LATERALITY vs SCALE (pillar-6, checkpoint item 7). Script 15 correlated the
    trace with implant left-fraction only at s=5.6. The per-scale trace is already
    cached (ms_mst020/per_patient_scale.csv, 16 scales); this re-reads it and asks
    at EVERY scale whether the inter-patient trace spread tracks left-hemisphere
    contact fraction — "is the laterality explanation scale-robust, or a single-scale
    coincidence?" Reports rho(s), p(s) for alpha and beta.

(2) DURATION control on the inference-specific trace (R2.9, checkpoint via audit_113).
    The persistent inference-specific beta component (T_infspec_pe, from the cached
    enc_inf_arc_mst020 arc) must not be a recording-length artifact of the
    test-vs-learn duration difference (f = D_test - D_learn). Correlates per-patient
    T_infspec_pe(beta) against the test/learn sample-count ratio (length_control.csv)
    at s1, the cohort mesoscale, and each patient's arc peak. Expect |rho| small,
    p large -> inference is not length-driven.

Output: data/sparsified_arc/laterality_duration/{laterality_vs_scale,duration_infspec}.csv
"""
from __future__ import annotations
import glob, os
import numpy as np, pandas as pd
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

MS = ROOT / "data/sparsified_arc/ms_mst020/per_patient_scale.csv"
GEO = ROOT / "data/audit/implant_geometry/per_patient_features.csv"
ENC = ROOT / "data/sparsified_arc/enc_inf_arc_mst020"
LEN = ROOT / "data/audit/inference_localization/length_control.csv"
OUT = ROOT / "data/sparsified_arc/laterality_duration"


def laterality_vs_scale():
    ms = pd.read_csv(MS)
    geo = pd.read_csv(GEO).set_index("patient")
    left_frac = (geo["N_L"] / (geo["N_L"] + geo["N_R"]).replace(0, np.nan))
    rows = []
    for band in ["alpha", "beta"]:
        for s in np.sort(ms.s.unique()):
            x = ms[(ms.band == band) & np.isclose(ms.s, s)].set_index("patient")
            lf = left_frac.reindex(x.index)
            ok = lf.notna() & x["obs_rho"].notna()
            if ok.sum() < 5:
                continue
            r, p = spearmanr(lf[ok], x["obs_rho"][ok])
            rows.append(dict(band=band, s=float(s), rho=float(r), p=float(p), n=int(ok.sum())))
    return pd.DataFrame(rows)


def duration_infspec():
    # per-patient duration ratio (test/learn sample counts)
    L = pd.read_csv(LEN)
    ratio = L.groupby("pat")["ratio"].first()
    # per-patient T_infspec_pe(beta) over scales
    recs = {}
    for f in sorted(glob.glob(str(ENC / "beta" / "*.npz"))):
        pat = os.path.basename(f).replace(".npz", "")
        d = np.load(f)
        recs[pat] = (d["s"], d["T_infspec_pe__obs"])
    pats = [p for p in recs if p in ratio.index]
    s_grid = recs[pats[0]][0]
    j1 = int(np.argmin(np.abs(s_grid - 1.0)))
    jmeso = int(np.argmin(np.abs(s_grid - 5.6)))
    rat = np.array([ratio[p] for p in pats])
    rows = []
    for label, getter in [
        ("s1", lambda v: v[j1]),
        ("mesoscale", lambda v: v[jmeso]),
        ("arc_peak", lambda v: v[np.nanargmax(v)]),
    ]:
        y = np.array([getter(recs[p][1]) for p in pats])
        ok = np.isfinite(y) & np.isfinite(rat)
        r, p = spearmanr(rat[ok], y[ok])
        rows.append(dict(scale_ref=label, rho=float(r), p=float(p), n=int(ok.sum()),
                         infspec_med=float(np.median(y[ok]))))
    return pd.DataFrame(rows), pats, rat


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lat = laterality_vs_scale(); lat.to_csv(OUT / "laterality_vs_scale.csv", index=False)
    dur, pats, rat = duration_infspec(); dur.to_csv(OUT / "duration_infspec.csv", index=False)

    print("=== (1) LATERALITY vs SCALE — Spearman(left_frac, trace rho) across patients ===")
    for band in ["alpha", "beta"]:
        b = lat[lat.band == band]
        sig = b[b.p < 0.05]
        print(f"  {band}: clears p<.05 at {len(sig)}/{len(b)} scales; "
              f"rho range [{b.rho.min():+.2f}, {b.rho.max():+.2f}]")
        for _, r in b.iterrows():
            flag = " *" if r.p < 0.05 else ""
            print(f"      s={r.s:6.2f}  rho={r.rho:+.3f}  p={r.p:.3f}{flag}")
    print("\n=== (2) DURATION control — Spearman(test/learn ratio, T_infspec_pe beta) ===")
    print(f"  patients (n={len(pats)}): ratio range [{rat.min():.2f}, {rat.max():.2f}]")
    for _, r in dur.iterrows():
        print(f"  {r.scale_ref:10s}: rho={r.rho:+.3f}  p={r.p:.3f}  "
              f"(infspec_med={r.infspec_med:+.3f}, n={int(r.n)})")
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
