#!/usr/bin/env python3
"""Is DENSE raw-FC non-selective across bands? (why the confound-fixed ladder shows fewer)

The 2026-07-13 ladder reads raw_fc on the mst@0.20 backbone (apples-to-apples read-out
control). This tests the OTHER comparison the reviewer/PI cares about: the STANDARD-PRACTICE
pairwise measure = raw FC on the FULL DENSE graph. Prediction (PI memory): dense raw FC
traces in essentially EVERY band incl theta -> it cannot discriminate cognitive from
non-cognitive bands, whereas the cophenetic hierarchy fires only in alpha/beta.

Per (patient, band in all 6), dense |ImCoh| upper-tri, trace rho_sym (T_test) vs the
IDENTICAL matched-strength null. Cohort Wilcoxon(obs - surr_p50, greater), n=10.

Output: data/sparsified_arc/raw_dense_selectivity/cohort.csv
"""
from __future__ import annotations
import os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import rankdata, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import COHORT, load_phase            # noqa: E402
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle  # noqa: E402

PHASES = ("A", "B", "task_test", "rest_post")
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
R = 200
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260713


def _sp(a, b):
    ra, rb = rankdata(a), rankdata(b)
    ra = ra - ra.mean(); rb = rb - rb.mean()
    da, db = np.sqrt(np.dot(ra, ra)), np.sqrt(np.dot(rb, rb))
    return float(np.dot(ra, rb) / (da * db)) if da > 0 and db > 0 else np.nan


def trace_dense(Ws):
    N = Ws["A"].shape[0]; iu = np.triu_indices(N, 1)
    A = {ph: Ws[ph][iu] for ph in PHASES}
    return 0.5 * (_sp(A["task_test"] - A["A"], A["rest_post"] - A["B"])
                  + _sp(A["task_test"] - A["B"], A["rest_post"] - A["A"]))


def per_cell(job):
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES}
    except Exception:
        return None
    N = Ws["A"].shape[0]
    o = trace_dense(Ws)
    rng = np.random.default_rng(BASE_SEED + idx)
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    s = np.full(R, np.nan)
    for r in range(R):
        Wsh = {ph: matched_strength_shuffle(Ws[ph], n_swaps, rng, W_MAX) for ph in PHASES}
        try:
            s[r] = trace_dense(Wsh)
        except Exception:
            continue
    col = s[np.isfinite(s)]
    p = float((np.sum(col >= o) + 1) / (col.size + 1)) if col.size else np.nan
    return dict(patient=pat, band=band, obs=float(o),
                surr_p50=float(np.nanmedian(col)) if col.size else np.nan, p=p)


def main():
    OUT = ROOT / "data" / "sparsified_arc" / "raw_dense_selectivity"
    OUT.mkdir(parents=True, exist_ok=True)
    _ = matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))
    jobs = [(i, p, b) for i, (p, b) in enumerate((p, b) for b in BANDS for p in COHORT)]
    print(f"[raw-dense selectivity] R={R} {len(jobs)} cells", flush=True)
    t0 = time.time(); recs = []
    nproc = min(12, max(1, (os.cpu_count() or 2) - 2))
    with Pool(nproc) as pool:
        for k, rl in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if rl:
                recs.append(rl)
    df = pd.DataFrame(recs)
    out = []
    for band in BANDS:
        x = df[df.band == band].dropna(subset=["obs", "surr_p50"])
        if len(x) < 5:
            continue
        diff = x.obs.values - x.surr_p50.values
        try:
            _, p = wilcoxon(diff, alternative="greater")
        except ValueError:
            p = 1.0
        out.append(dict(band=band, n=len(x), obs_med=float(x.obs.median()),
                        n_pos=int((diff > 0).sum()), gate_p=float(p)))
    c = pd.DataFrame(out); c.to_csv(OUT / "cohort.csv", index=False)
    nfire = int((c.gate_p < 0.05).sum())
    print("\n=== DENSE raw-FC TRACE (T_test) per band, matched-strength ===", flush=True)
    for _, r in c.iterrows():
        print(f"  {r.band:11s} gate_p={r.gate_p:.3f}{'*' if r.gate_p<0.05 else ' '}  "
              f"({int(r.n_pos)}/{int(r.n)}, obs {r.obs_med:+.3f})", flush=True)
    print(f"\n  DENSE raw FC fires {nfire}/6 bands  (cophenetic fires 2/6: only alpha, beta)", flush=True)
    print(f"[done] {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
