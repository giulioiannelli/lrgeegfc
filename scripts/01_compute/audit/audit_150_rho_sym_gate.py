#!/usr/bin/env python3
"""audit_150 — canonical cohort trace gate under the rho_sym estimator (numba+MP).

Supersedes the split-baseline gate (audit_63 output) as the ESTIMATOR of record:
rho_sym removes the arbitrary split-half arm assignment that made bare rho_split
sign-unstable for near-zero patients (audit_149; feedback_rho_sym_canonical_estimator).
Does NOT edit audit_63 — new pipeline, new output dir. rho_split is kept as a
secondary column purely for the robustness supplement.

    rho_sym = 1/2 [ Spearman(D_task-D_preA, D_post-D_preB)
                  + Spearman(D_task-D_preB, D_post-D_preA) ]

Speed: the strength-preserving shuffle loop is JIT-compiled (numba, bit-identical,
98x faster: 882->9 ms) and cells run over a multiprocessing Pool. R=200 gate ~2 min
(was hours). See feedback_numba_for_surrogates.

5-point critical preamble
1. Claim: the cohort matched-strength gate (alpha, beta CLEAR; delta, theta,
   low_gamma, high_gamma FAIL) holds under rho_sym exactly as under rho_split.
2. Null: switching estimator flips a band verdict (the taxonomy depended on the
   arbitrary half assignment).
3. Strongest alternative: rho_sym is just rho_split relabelled and reproduces it
   trivially. Controlled by ALSO recording rho_split per cell (obs_rho_split) so the
   two are compared head-to-head, and by an independent per-cell RNG (not the
   audit_63 sequential rng) so agreement is not a shared-draw artifact.
4. Cannot: n=10 cohort Wilcoxon has coarse p-grid; per-cell seeds differ from
   audit_63 so this reproduces the VERDICT, not the bit value. Localization/arc are
   downstream, regenerated separately.
5. Falsify: if alpha or beta fails, or any of delta/theta/gamma clears, under
   rho_sym, the taxonomy is estimator-dependent.

Outputs:
    data/audit/rho_sym_gate/per_patient_per_band.csv   (obs_rho = rho_sym; + obs_rho_split)
    data/audit/rho_sym_gate/cohort_summary.csv         (per-band Wilcoxon gate)
    data/audit/rho_sym_gate/README.md
"""
from __future__ import annotations
import os, time
import numpy as np, pandas as pd, numba
from multiprocessing import Pool
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import cophenet, linkage
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.workflow.fc import load_fc_matrix

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
N_SURROGATES = 200
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260706
HALVES = CACHE_ROOT / "imcoh_halves_fc"
OUT = ROOT / "data" / "audit" / "rho_sym_gate"


# --- FC loading (halves cached by audit_63 pre-flight; full phases via loader) ---
def load_phase(pat, phase, band):
    if phase in ("A", "B"):
        W = np.load(HALVES / pat / f"{band}_rest_pre_{phase}_imcoh_abs.npy")
    else:
        W = load_fc_matrix(pat, phase, band, fc_method="imcoh_abs")
    W = np.asarray(W, float); np.fill_diagonal(W, 0.0); W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


def ultra(W):
    """LRG cophenetic condensed distances (tau=1/lam_max, T=1/rho, avg linkage)."""
    deg = W.sum(1); ev, V = np.linalg.eigh(np.diag(deg) - W); tau = 1.0 / ev[-1]
    rho = (V * np.exp(-tau * ev)) @ V.T; rho /= np.trace(rho)
    with np.errstate(divide="ignore"): T = 1.0 / rho
    T = np.maximum(T, T.T); np.fill_diagonal(T, 0.0); fin = np.isfinite(T)
    if not fin.all(): T = np.where(fin, T, np.nanmax(T[fin]))
    return cophenet(linkage(squareform(T, checks=False), method="average"))


@numba.njit
def _swap_loop(W, s, fr, w_max):
    for i in range(s.shape[0]):
        a = s[i, 0]; b = s[i, 1]; c = s[i, 2]; d = s[i, 3]
        if a == b or a == c or a == d or b == c or b == d or c == d:
            continue
        w1 = W[a, b]; w2 = W[c, d]; w3 = W[a, d]; w4 = W[c, b]
        lo = -w1
        if -w2 > lo: lo = -w2
        if w3 - w_max > lo: lo = w3 - w_max
        if w4 - w_max > lo: lo = w4 - w_max
        hi = w_max - w1
        if w_max - w2 < hi: hi = w_max - w2
        if w3 < hi: hi = w3
        if w4 < hi: hi = w4
        if lo >= hi:
            continue
        dl = lo + fr[i] * (hi - lo)
        W[a, b] = w1 + dl; W[b, a] = w1 + dl
        W[c, d] = w2 + dl; W[d, c] = w2 + dl
        W[a, d] = w3 - dl; W[d, a] = w3 - dl
        W[c, b] = w4 - dl; W[b, c] = w4 - dl
    return W


def shuffle(W, n_swaps, rng, w_max=W_MAX):
    """Bit-identical to the audit_63 pure-Python 4-cycle shuffle; rng draws the
    (samples, fracs) arrays outside the jitted loop -> reproducible + fast."""
    W = W.copy(); N = W.shape[0]
    s = rng.integers(0, N, size=(n_swaps, 4)).astype(np.int64)
    fr = rng.uniform(0.0, 1.0, size=n_swaps)
    return _swap_loop(W, s, fr, float(w_max))


def rho_sym_split(DA, DB, Dt, Dp):
    """(rho_sym, rho_split). rho_sym averages both arm assignments."""
    r_ab, _ = spearmanr(Dt - DA, Dp - DB)
    r_ba, _ = spearmanr(Dt - DB, Dp - DA)
    return float(0.5 * (r_ab + r_ba)), float(r_ab)


def per_cell(job):
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in ("A", "B", "task_test", "rest_post")}
    except Exception:
        return None
    D = {k: ultra(v) for k, v in Ws.items()}
    if len({d.size for d in D.values()}) != 1:
        return None
    rho_obs, rho_obs_split = rho_sym_split(D["A"], D["B"], D["task_test"], D["rest_post"])
    N = Ws["A"].shape[0]; n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    rng = np.random.default_rng(BASE_SEED + idx)
    surr = np.empty(N_SURROGATES, float)
    for r in range(N_SURROGATES):
        Ds = {k: ultra(shuffle(v, n_swaps, rng)) for k, v in Ws.items()}
        surr[r], _ = rho_sym_split(Ds["A"], Ds["B"], Ds["task_test"], Ds["rest_post"])
    s = surr[np.isfinite(surr)]
    if s.size == 0:
        return None
    sm, sd = float(np.mean(s)), float(np.std(s, ddof=1))
    return {
        "patient": pat, "band": band, "N_nodes": int(N), "n_pairs": int(D["A"].size),
        "n_surrogates": int(s.size), "obs_rho": rho_obs, "obs_rho_split": rho_obs_split,
        "surr_mean_rho": sm, "surr_std_rho": sd,
        "surr_p5": float(np.quantile(s, 0.05)), "surr_p50": float(np.quantile(s, 0.50)),
        "surr_p95": float(np.quantile(s, 0.95)),
        "obs_z": (rho_obs - sm) / sd if sd > 0 else float("nan"),
        "obs_p_one_sided": float(np.mean(s >= rho_obs)),
        "obs_p_split": float(np.mean(s >= rho_obs_split)),
    }


def cohort_gate(df):
    out = []
    for band in BANDS:
        x = df[df.band == band]
        if x.empty:
            continue
        rec = {"band": band, "n_patients": len(x),
               "obs_median_rho_sym": float(x.obs_rho.median()),
               "surr_median_rho": float(x.surr_p50.median()),
               "n_above": int((x.obs_p_one_sided < 0.05).sum())}
        for col, tag in [("obs_rho", "sym"), ("obs_rho_split", "split")]:
            try:
                _, p = wilcoxon(x[col].values - x.surr_p50.values, alternative="greater")
            except Exception:
                p = float("nan")
            rec[f"gate_p_{tag}"] = float(p)
        rec["verdict_sym"] = "CLEAR" if rec["gate_p_sym"] < 0.05 else "fail"
        out.append(rec)
    return pd.DataFrame(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # warm-compile numba in parent before fork
    _ = shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))
    jobs = [(i, p, b) for i, (b, p) in enumerate((b, p) for b in BANDS for p in COHORT)]
    ncpu = min(14, (os.cpu_count() or 4) - 2)
    print(f"[audit_150] {len(jobs)} cells, R={N_SURROGATES}, rho_sym, {ncpu} workers")
    t0 = time.time(); rows = []
    with Pool(ncpu) as pool:
        for i, r in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if r:
                rows.append(r)
                print(f"[{i}/{len(jobs)}] {r['patient']}/{r['band']} "
                      f"obs_sym={r['obs_rho']:+.3f} z={r['obs_z']:+.2f} p={r['obs_p_one_sided']:.3f}",
                      flush=True)
    df = pd.DataFrame(rows).sort_values(["band", "patient"])
    df.to_csv(OUT / "per_patient_per_band.csv", index=False)
    gate = cohort_gate(df); gate.to_csv(OUT / "cohort_summary.csv", index=False)
    rt = time.time() - t0
    print(f"\n[audit_150] {len(rows)} cells in {rt:.0f}s")
    print("\n=== COHORT MATCHED-STRENGTH GATE under rho_sym (own surrogate) ===")
    print(gate[["band", "obs_median_rho_sym", "gate_p_sym", "gate_p_split", "n_above", "verdict_sym"]].to_string(index=False))
    write_readme(gate, rt)
    print(f"[audit_150] outputs -> {OUT}")


def write_readme(gate, rt):
    L = ["---", "name: rho_sym_gate",
         "scope: canonical_cohort_trace_gate_under_rho_sym_estimator",
         "date: 2026-07-06", "status: current", "---", "",
         "# Cohort trace gate under rho_sym (canonical estimator)", "",
         "**Head.** rho_sym = mean of the two split-half arm assignments -- removes the",
         "arbitrary-half artifact of bare rho_split (audit_149). Matched-strength null",
         f"R={N_SURROGATES}, numba-accelerated ({rt:.0f}s total). Verdict per band below;",
         "rho_split kept as a supplement column (gate_p_split) for head-to-head.", "",
         "| band | obs median rho_sym | gate p (sym) | gate p (split) | n>surr | verdict |",
         "|---|---|---|---|---|---|"]
    for _, r in gate.iterrows():
        L.append(f"| {r['band']} | {r['obs_median_rho_sym']:+.3f} | {r['gate_p_sym']:.4f} "
                 f"| {r['gate_p_split']:.4f} | {r['n_above']}/{r['n_patients']} | **{r['verdict_sym']}** |")
    L += ["", "## Provenance",
          f"- estimator: rho_sym (primary); rho_split (supplement).",
          f"- R={N_SURROGATES} 4-cycle +/-delta strength-preserving surrogate, SWAP_FACTOR={SWAP_FACTOR}.",
          "- numba-JIT shuffle, bit-identical to audit_63 pure-Python (feedback_numba_for_surrogates).",
          "- independent per-cell rng (BASE_SEED + cell index).",
          "- cohort gate: Wilcoxon(obs - surr_p50) one-sided greater over 10 patients.",
          "- Build: scripts/01_compute/audit/audit_150_rho_sym_gate.py",
          "- Does NOT modify audit_63; supersedes it as estimator of record."]
    (OUT / "README.md").write_text("\n".join(L))


if __name__ == "__main__":
    main()
