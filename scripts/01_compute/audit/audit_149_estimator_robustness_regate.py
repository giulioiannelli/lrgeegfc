#!/usr/bin/env python3
"""audit_149 — is the band taxonomy invariant to the ρ_split ESTIMATOR?

Settles the estimator-conditioning question (audit_148 follow-up): the split-
baseline ρ_split is heteroscedastic / A-B-split-dependent for low-reorganizability
patients. Does that instability change the COHORT MATCHED-STRENGTH GATE — the real
trace verdict — for any band? Re-runs the identical audit_63 gate under three
estimators, changing ONLY the per-cell ρ function; surrogate rewiring and the
cohort Wilcoxon are untouched.

Estimators (all share the same strength-preserving surrogate FCs):
  1. split : Spearman(D_task-D_preA, D_post-D_preB)            [published, audit_63]
  2. sym   : ½[ρ(A→task,B→rest) + ρ(B→task,A→rest)]           [removes A/B asymmetry]
  3. full  : Spearman(D_task-D_preFULL, D_post-D_preFULL)      [full reliable baseline;
             shared-baseline inflation is COMMON-MODE to obs & surrogate -> cancels in
             obs_p = mean(surr>=obs) and in obs_z = (obs-surr_mean)/surr_std]

5-point preamble
1. Claim: the cohort gate taxonomy (α,β CLEAR; δ,θ,low_γ,high_γ FAIL) is invariant
   to the estimator; β/α flagship status does not rest on the ill-conditioned split.
2. Null/baseline: the gate verdicts flip under a better estimator (the taxonomy was
   an artifact of the split choice).
3. Strongest alternative to beat: the full-baseline estimator is inflated and would
   false-positive every band. Controlled by giving it its OWN matched-strength
   surrogate (same shared full baseline) so the inflation is common-mode; the gate
   is obs-vs-surrogate, not raw ρ.
4. Cannot: R=200 gate p has Monte-Carlo error ~±1 rank; per-cell seeds differ from
   the published sequential RNG, so estimator-1 p reproduces the published VERDICT,
   not the bit value. Localization (β→OFC) is downstream and not re-run here.
5. Falsify: if α or β fails, or any of δ/θ/γ clears, under sym or full, the taxonomy
   is estimator-dependent and the flagship claim needs the better estimator to hold.

Outputs:
  data/audit/estimator_regate/per_cell.csv
  data/audit/estimator_regate/cohort_gate.csv
"""
from __future__ import annotations
import numpy as np, pandas as pd
from pathlib import Path
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import cophenet, linkage
from scipy.stats import spearmanr, wilcoxon
from multiprocessing import Pool
import os, sys

sys.path.insert(0, "/home/giulio/Documents/research/neural_networks/lrgeegfc/src")
from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.workflow.fc import load_fc_matrix

ROOT = Path("/home/giulio/Documents/research/neural_networks/lrgeegfc")
OUT = ROOT / "data/audit/estimator_regate"; OUT.mkdir(parents=True, exist_ok=True)
HALVES = CACHE_ROOT / "imcoh_halves_fc"
COH = ["Pat_02","Pat_03","Pat_05","Pat_06","Pat_07","Pat_08","Pat_10","Pat_13","Pat_14","Pat_15"]
BANDS = ["delta","theta","alpha","beta","low_gamma","high_gamma"]
R = 200; SWAP_FACTOR = 20; W_MAX = 1.0; BASE_SEED = 20260511


def load_phase(pat, phase, band):
    if phase in ("A", "B"):
        W = np.load(HALVES / pat / f"{band}_rest_pre_{phase}_imcoh_abs.npy")
    else:
        W = load_fc_matrix(pat, phase, band, fc_method="imcoh_abs")
    W = np.asarray(W, float); np.fill_diagonal(W, 0.0); W = np.clip(W, 0, 1)
    return 0.5 * (W + W.T)


def ultra(W):
    deg = W.sum(1); ev, V = np.linalg.eigh(np.diag(deg) - W); tau = 1.0 / ev[-1]
    rho = (V * np.exp(-tau * ev)) @ V.T; rho /= np.trace(rho)
    with np.errstate(divide="ignore"): T = 1.0 / rho
    T = np.maximum(T, T.T); np.fill_diagonal(T, 0.0); fin = np.isfinite(T)
    if not fin.all(): T = np.where(fin, T, np.nanmax(T[fin]))
    return cophenet(linkage(squareform(T, checks=False), method="average"))


def shuffle(W, n_swaps, rng, w_max=W_MAX):
    W = W.copy(); N = W.shape[0]
    s = rng.integers(0, N, size=(n_swaps, 4)); fr = rng.uniform(0, 1, size=n_swaps)
    for i in range(n_swaps):
        a, b, c, d = s[i]
        if len({a, b, c, d}) < 4: continue
        w1, w2, w3, w4 = W[a, b], W[c, d], W[a, d], W[c, b]
        lo = max(-w1, -w2, w3 - w_max, w4 - w_max); hi = min(w_max - w1, w_max - w2, w3, w4)
        if lo >= hi: continue
        dl = lo + fr[i] * (hi - lo)
        W[a, b] = W[b, a] = w1 + dl; W[c, d] = W[d, c] = w2 + dl
        W[a, d] = W[d, a] = w3 - dl; W[c, b] = W[b, c] = w4 - dl
    return W


def three_rhos(DA, DB, DF, Dt, Dp):
    r_ab = spearmanr(Dt - DA, Dp - DB)[0]
    r_ba = spearmanr(Dt - DB, Dp - DA)[0]
    r_full = spearmanr(Dt - DF, Dp - DF)[0]
    return float(r_ab), float(0.5 * (r_ab + r_ba)), float(r_full)


def worker(job):
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in ("A", "B", "rest_pre", "task_test", "rest_post")}
    except Exception as e:
        return None
    D = {k: ultra(v) for k, v in Ws.items()}
    obs = three_rhos(D["A"], D["B"], D["rest_pre"], D["task_test"], D["rest_post"])
    N = Ws["A"].shape[0]; n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    rng = np.random.default_rng(BASE_SEED + idx)
    surr = {"split": [], "sym": [], "full": []}
    for _ in range(R):
        Ds = {k: ultra(shuffle(v, n_swaps, rng)) for k, v in Ws.items()}
        rs = three_rhos(Ds["A"], Ds["B"], Ds["rest_pre"], Ds["task_test"], Ds["rest_post"])
        for est, val in zip(("split", "sym", "full"), rs): surr[est].append(val)
    row = {"patient": pat, "band": band}
    for est, o in zip(("split", "sym", "full"), obs):
        s = np.array([x for x in surr[est] if np.isfinite(x)])
        row[f"obs_{est}"] = o
        row[f"surrp50_{est}"] = float(np.median(s))
        row[f"obsp_{est}"] = float(np.mean(s >= o))
    return row


def main():
    jobs = [(i, p, b) for i, (b, p) in enumerate((b, p) for b in BANDS for p in COH)]
    ncpu = min(14, (os.cpu_count() or 4) - 2)
    print(f"[audit_149] {len(jobs)} cells, R={R}, {ncpu} workers")
    with Pool(ncpu) as pool:
        rows = [r for r in pool.map(worker, jobs) if r is not None]
    df = pd.DataFrame(rows); df.to_csv(OUT / "per_cell.csv", index=False)

    out = []
    for b in BANDS:
        x = df[df.band == b]
        rec = {"band": b}
        for est in ("split", "sym", "full"):
            diff = (x[f"obs_{est}"] - x[f"surrp50_{est}"]).values
            try: _, p = wilcoxon(diff, alternative="greater")
            except Exception: p = float("nan")
            rec[f"gate_p_{est}"] = round(float(p), 4)
            rec[f"nclear_{est}"] = int((x[f"obsp_{est}"] < 0.05).sum())
        out.append(rec)
    g = pd.DataFrame(out); g.to_csv(OUT / "cohort_gate.csv", index=False)
    print("\n=== COHORT MATCHED-STRENGTH GATE p (Wilcoxon, one-sided) per estimator ===")
    print(g.to_string(index=False))
    print("\nVerdict rule: CLEARS if gate_p < 0.05")
    for _, r in g.iterrows():
        v = {e: ("CLEAR" if r[f"gate_p_{e}"] < 0.05 else "fail") for e in ("split", "sym", "full")}
        flag = "" if len(set(v.values())) == 1 else "   <-- ESTIMATOR-DEPENDENT"
        print(f"  {r['band']:11s} split={v['split']:5s} sym={v['sym']:5s} full={v['full']:5s}{flag}")


if __name__ == "__main__":
    main()
