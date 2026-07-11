#!/usr/bin/env python3
"""audit_176 — is the raw-FC (no-diffusion) trace a property of the CONNECTIVITY,
or an artifact of the 1/A reciprocal distance form?

audit_174 showed that on the fully-connected graph at tau_min the DIFFUSION cophenetic
tree (D=1/K(tau_min)) and the RAW-FC tree (D=1/A) give the SAME beta gate (p=.032) but
diverge on alpha (diffusion clears, raw fails) and low_gamma (raw clears, diffusion
suppresses). This audit closes audit_174's open follow-up #1: the raw arm used ONE
distance functional (reciprocal 1/A). UPGMA average-linkage is NOT invariant to monotone
transforms of the distance, so the raw verdict could hinge on the reciprocal form (which
is exactly the tau_min heat-kernel linearization K_ij ~ tau*A_ij => 1/K ~ (1/tau)(1/A)).

We therefore rerun the canonical rho_sym cohort matched-strength gate for FOUR arms fed
the IDENTICAL surrogate draws (audit_150 seed/order -> diffusion arm reproduces audit_150
bit-for-bit as anchor):
    diff     : D = 1/rho(tau_min)          (the pipeline; == audit_150 ultra())
    raw_recip: D = 1/A                      (== audit_174 ultra_raw)
    raw_1mA  : D = 1 - A                    (bounded monotone-decreasing)
    raw_nlogA: D = -log(clip(A))            (log-reciprocal; another monotone form)

5-point critical preamble
1. Claim: the beta raw-FC clearing (and the low_gamma raw clearing) is a property of the
   FC rank structure, robust across monotone distance forms -- NOT specific to 1/A.
2. Null: the raw verdict is 1/A-specific (beta or low_gamma clear ONLY under the reciprocal
   form), i.e. the raw arm is a functional-form coincidence, not a connectivity property.
3. Strongest alternative it must control for: a shared-draw / estimator artifact. Controlled
   by feeding ALL four arms the SAME surrogate FC (one shuffle -> four trees), same rho_sym,
   same Wilcoxon; diff arm reproduces audit_150 exactly (anchor).
4. Cannot: this is at tau_min on the FULLY-CONNECTED graph (the degenerate LRG regime;
   Villegas 2025 App. F Wigner-semicircle). It does NOT test sparsification (audit_175) nor
   coarser tau (audit_172, where alpha's genuine mesoscale trace lives). It answers only:
   "given no diffusion, does the beta/low_gamma raw trace depend on the distance functional?"
5. Falsify: if beta clears under all three raw forms and diffusion, beta is a connectivity
   property (diffusion not load-bearing for beta). If low_gamma clears ONLY under 1/A, its
   raw 'trace' is a reciprocal-form artifact (further evidence the diffusion correctly
   filters it). If beta clears ONLY under 1/A, the reciprocal form is load-bearing.

Outputs:
    data/audit/rawform_robustness/per_patient_per_band.csv
    data/audit/rawform_robustness/cohort_summary.csv
    data/audit/rawform_robustness/README.md
"""
from __future__ import annotations
import argparse, os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import cophenet, linkage
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import (            # reuse the canonical gate machinery (DRY)
    load_phase, ultra, shuffle, rho_sym_split,
    COHORT, BANDS, N_SURROGATES, SWAP_FACTOR, BASE_SEED,
)
from audit_174_diffusion_vs_rawfc_gate import ultra_raw   # D = 1/A (reciprocal)

OUT = ROOT / "data" / "audit" / "rawform_robustness"
R = N_SURROGATES
ARMS = ("diff", "raw_recip", "raw_1mA", "raw_nlogA")


def _coph_from_D(T):
    """Shared tail: symmetrize, zero diagonal, replace non-finite with row-safe max, UPGMA."""
    T = np.maximum(T, T.T); np.fill_diagonal(T, 0.0); fin = np.isfinite(T)
    if not fin.all():
        T = np.where(fin, T, np.nanmax(T[fin]))
    return cophenet(linkage(squareform(T, checks=False), method="average"))


def raw_1mA(W):
    """Bounded monotone-decreasing distance D = 1 - A (A in [0,1])."""
    return _coph_from_D(1.0 - W)


def raw_nlogA(W):
    """Log-reciprocal distance D = -log(A), A clipped away from 0 (fully-connected => A>0)."""
    with np.errstate(divide="ignore"):
        T = -np.log(np.clip(W, 1e-12, 1.0))
    return _coph_from_D(T)


TREE = {"diff": ultra, "raw_recip": ultra_raw, "raw_1mA": raw_1mA, "raw_nlogA": raw_nlogA}


def _rho_arm(fn, Ws):
    D = {k: fn(v) for k, v in Ws.items()}
    if len({d.size for d in D.values()}) != 1:
        return None
    r, _ = rho_sym_split(D["A"], D["B"], D["task_test"], D["rest_post"])
    return r


def _pack(obs, surr, tag):
    s = surr[np.isfinite(surr)]
    if s.size == 0 or obs is None:
        return {f"obs_rho_{tag}": obs if obs is not None else float("nan"),
                f"surr_p50_{tag}": float("nan"), f"obs_p_{tag}": float("nan")}
    sm, sd = float(np.mean(s)), float(np.std(s, ddof=1))
    return {f"obs_rho_{tag}": obs, f"surr_p50_{tag}": float(np.quantile(s, 0.50)),
            f"obs_z_{tag}": (obs - sm) / sd if sd > 0 else float("nan"),
            f"obs_p_{tag}": float(np.mean(s >= obs))}


def per_cell(job):
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in ("A", "B", "task_test", "rest_post")}
    except Exception:
        return None
    obs = {a: _rho_arm(TREE[a], Ws) for a in ARMS}
    if obs["diff"] is None:
        return None
    N = Ws["A"].shape[0]; n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    rng = np.random.default_rng(BASE_SEED + idx)      # SAME seed/order as audit_150/174
    surr = {a: np.empty(R, float) for a in ARMS}
    for r in range(R):
        Wsurr = {k: shuffle(v, n_swaps, rng) for k, v in Ws.items()}   # one draw -> all arms
        for a in ARMS:
            v = _rho_arm(TREE[a], Wsurr)
            surr[a][r] = v if v is not None else np.nan
    rec = {"patient": pat, "band": band, "N_nodes": int(N)}
    for a in ARMS:
        rec.update(_pack(obs[a], surr[a], a))
    return rec


def cohort_gate(df):
    out = []
    for band in BANDS:
        x = df[df.band == band]
        if x.empty:
            continue
        rec = {"band": band, "n_patients": len(x)}
        for a in ARMS:
            try:
                _, p = wilcoxon(x[f"obs_rho_{a}"].values - x[f"surr_p50_{a}"].values,
                                alternative="greater")
            except Exception:
                p = float("nan")
            rec[f"gate_p_{a}"] = float(p)
            rec[f"median_{a}"] = float(x[f"obs_rho_{a}"].median())
            rec[f"n_above_{a}"] = int((x[f"obs_p_{a}"] < 0.05).sum())
            rec[f"verdict_{a}"] = "CLEAR" if p < 0.05 else "fail"
        out.append(rec)
    return pd.DataFrame(out)


def write_readme(gate, rt):
    L = ["---", "name: rawform_robustness",
         "scope: is_the_raw_FC_trace_a_reciprocal_form_artifact_or_a_connectivity_property",
         "date: 2026-07-11", "status: current", "---", "",
         "# Raw-FC trace: connectivity property or 1/A form artifact? (audit_176)", "",
         "**Head.** Canonical rho_sym matched-strength cohort gate on FOUR cophenetic arms fed",
         f"identical surrogate draws; R={R}, {rt:.0f}s. diff reproduces audit_150 (anchor).", "",
         "| band | diff | 1/A | 1-A | -logA |",
         "|---|---|---|---|---|"]
    for _, r in gate.iterrows():
        L.append(f"| {r['band']} | {r['gate_p_diff']:.4f} {r['verdict_diff']} "
                 f"| {r['gate_p_raw_recip']:.4f} {r['verdict_raw_recip']} "
                 f"| {r['gate_p_raw_1mA']:.4f} {r['verdict_raw_1mA']} "
                 f"| {r['gate_p_raw_nlogA']:.4f} {r['verdict_raw_nlogA']} |")
    L += ["", "cell = cohort Wilcoxon gate p (obs - surr_p50, one-sided greater); CLEAR = p<0.05.",
          "", "## Provenance",
          "- ablation of audit_174; adds 1-A and -logA raw distance forms.",
          f"- R={R} 4-cycle strength-preserving surrogate, SWAP_FACTOR={SWAP_FACTOR}, audit_150 seed/order.",
          "- one surrogate draw feeds all four arms (pure function-of-distance comparison).",
          "- Build: scripts/01_compute/audit/audit_176_rawform_robustness.py"]
    (OUT / "README.md").write_text("\n".join(L))


def main():
    global R
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=int, default=N_SURROGATES)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    R = a.R
    OUT.mkdir(parents=True, exist_ok=True)
    _ = shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))     # warm numba before fork
    jobs = [(i, p, b) for i, (b, p) in enumerate((b, p) for b in BANDS for p in COHORT)]
    if a.limit:
        jobs = jobs[:a.limit]
    ncpu = min(14, (os.cpu_count() or 4) - 2)
    print(f"[audit_176] {len(jobs)} cells, R={R}, 4 arms, {ncpu} workers", flush=True)
    t0 = time.time(); rows = []
    with Pool(ncpu) as pool:
        for i, r in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if r:
                rows.append(r)
                el = time.time() - t0; eta = el / i * (len(jobs) - i)
                print(f"[{i}/{len(jobs)}] {r['patient']}/{r['band']} "
                      f"diff={r['obs_rho_diff']:+.3f}(p{r['obs_p_diff']:.3f}) "
                      f"1/A={r['obs_rho_raw_recip']:+.3f}(p{r['obs_p_raw_recip']:.3f}) "
                      f"1-A={r['obs_rho_raw_1mA']:+.3f}(p{r['obs_p_raw_1mA']:.3f}) "
                      f"-logA={r['obs_rho_raw_nlogA']:+.3f}(p{r['obs_p_raw_nlogA']:.3f}) "
                      f"[{el:.0f}s ETA {eta:.0f}s]", flush=True)
    df = pd.DataFrame(rows).sort_values(["band", "patient"])
    df.to_csv(OUT / "per_patient_per_band.csv", index=False)
    gate = cohort_gate(df); gate.to_csv(OUT / "cohort_summary.csv", index=False)
    rt = time.time() - t0
    print(f"\n[audit_176] {len(rows)} cells in {rt:.0f}s")
    print("\n=== rho_sym matched-strength gate: 4 distance arms (gate p / verdict) ===")
    cols = ["band"] + [f"{k}_{a}" for a in ARMS for k in ("gate_p", "verdict")]
    print(gate[[c for c in cols if c in gate.columns]].to_string(index=False))
    if a.R == N_SURROGATES and not a.limit:
        write_readme(gate, rt)
    print(f"[audit_176] outputs -> {OUT}")


if __name__ == "__main__":
    main()
