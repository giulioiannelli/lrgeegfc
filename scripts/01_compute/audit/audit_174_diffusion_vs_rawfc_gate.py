#!/usr/bin/env python3
"""audit_174 — does the band-selective trace NEED the diffusion propagator?

Ablation of the canonical rho_sym cohort gate (audit_150): the DIFFUSION cophenetic
tree (D = 1/K(tau_min), the pipeline) vs a RAW-FC cophenetic tree (D = 1/A, diffusion
skipped). SAME matched-strength surrogate draws, SAME rho_sym estimator, SAME cohort
Wilcoxon — the ONLY thing that changes is the tree construction (1/rho(tau_min) vs 1/A).

Motivation. At tau = 1/lambda_max the propagator K(tau) is ~a monotone map of the FC
(Spearman 0.98, Pat_05 beta) and the diffusion cophenetic tree is ~0.93 rank-correlated
with the raw-FC tree — i.e. at tau_min the propagator does MODEST work. The question this
audit settles: does the published band gate (alpha p.024, beta p.032 CLEAR; delta/theta/
low_gamma/high_gamma FAIL, audit_150) SURVIVE when we cluster the FC directly, skipping
diffusion?

5-point critical preamble
1. Claim: the diffusion propagator is load-bearing — the band-selective gate requires the
   diffusion cophenetic tree and does NOT reproduce on raw-FC (1/A) trees.
2. Null: raw-FC reproduces the gate (same bands CLEAR/FAIL, matching gate p) => at tau_min
   the diffusion is ~cosmetic for the RESULT and the pipeline could be simplified to
   UPGMA-on-FC. That is the honest outcome to report if it happens.
3. Strongest alternative it must control for: a shared-draw / estimator artifact that makes
   the two arms agree (or disagree) trivially. Controlled by feeding BOTH arms the IDENTICAL
   surrogate FC matrices (same rng, same draw order as audit_150), the same rho_sym estimator,
   and the same Wilcoxon gate. The diffusion arm therefore reproduces audit_150 bit-for-bit
   (cross-check anchor); only 1/K vs 1/A differs.
4. Cannot: this is at tau_min (our cross-phase operating point; audit_121). It does NOT test
   whether a DIFFERENT tau gives a different trace — that is audit_172/173 (the tau-sweep:
   beta scale-broad, alpha mesoscale, coarse tau = collapse). A tau-divergence diagnostic is
   printed (tree(s) vs raw-FC as s = tau*lambda_max grows) to show the propagator IS a
   nontrivial tau-family, so "raw-FC vs diffusion" here means "raw-FC vs diffusion-at-tau_min".
5. Falsify: if raw-FC clears alpha AND beta and fails the rest with matching gate p's, the
   propagator is inessential at tau_min. If raw-FC loses alpha/beta, or lights other bands,
   the diffusion earns its place.

Outputs:
    data/audit/diffusion_vs_rawfc_gate/per_patient_per_band.csv
    data/audit/diffusion_vs_rawfc_gate/cohort_summary.csv   (per-band gate: diff vs raw)
    data/audit/diffusion_vs_rawfc_gate/README.md
"""
from __future__ import annotations
import argparse, os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import cophenet, linkage
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import (            # reuse the canonical gate machinery (DRY)
    load_phase, ultra, shuffle, rho_sym_split,
    COHORT, BANDS, N_SURROGATES, SWAP_FACTOR, BASE_SEED,
)

OUT = ROOT / "data" / "audit" / "diffusion_vs_rawfc_gate"
S_GRID = (1.0, 2.0, 4.0, 8.0)          # tau-divergence: s = tau * lambda_max (s=1 == tau_min)
R = N_SURROGATES                       # overridable via --R (set in main before Pool fork)


def ultra_raw(W):
    """Raw-FC cophenetic (diffusion SKIPPED): distance D = 1/W, average linkage.
    Identical UPGMA + cophenet to ultra(); ONLY the distance differs (1/W vs 1/rho(tau_min))."""
    with np.errstate(divide="ignore"):
        T = 1.0 / W
    T = np.maximum(T, T.T); np.fill_diagonal(T, 0.0); fin = np.isfinite(T)
    if not fin.all():
        T = np.where(fin, T, np.nanmax(T[fin]))
    return cophenet(linkage(squareform(T, checks=False), method="average"))


def ultra_s(W, s):
    """Diffusion cophenetic at scale s = tau * lambda_max (s=1 -> tau_min == ultra())."""
    deg = W.sum(1); ev, V = np.linalg.eigh(np.diag(deg) - W); tau = s / ev[-1]
    rho = (V * np.exp(-tau * ev)) @ V.T; rho /= np.trace(rho)
    with np.errstate(divide="ignore"):
        T = 1.0 / rho
    T = np.maximum(T, T.T); np.fill_diagonal(T, 0.0); fin = np.isfinite(T)
    if not fin.all():
        T = np.where(fin, T, np.nanmax(T[fin]))
    return cophenet(linkage(squareform(T, checks=False), method="average"))


def _pack(obs, obs_split, surr, tag):
    s = surr[np.isfinite(surr)]
    if s.size == 0:
        return {f"obs_rho_{tag}": obs, f"surr_p50_{tag}": float("nan"),
                f"obs_z_{tag}": float("nan"), f"obs_p_{tag}": float("nan")}
    sm, sd = float(np.mean(s)), float(np.std(s, ddof=1))
    return {
        f"obs_rho_{tag}": obs, f"obs_split_{tag}": obs_split,
        f"surr_mean_{tag}": sm, f"surr_p50_{tag}": float(np.quantile(s, 0.50)),
        f"obs_z_{tag}": (obs - sm) / sd if sd > 0 else float("nan"),
        f"obs_p_{tag}": float(np.mean(s >= obs)),
    }


def per_cell(job):
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in ("A", "B", "task_test", "rest_post")}
    except Exception:
        return None
    Dd = {k: ultra(v) for k, v in Ws.items()}         # diffusion arm (== audit_150)
    Dr = {k: ultra_raw(v) for k, v in Ws.items()}     # raw-FC arm (skip diffusion)
    if len({d.size for d in Dd.values()} | {d.size for d in Dr.values()}) != 1:
        return None
    obs_d, obs_d_s = rho_sym_split(Dd["A"], Dd["B"], Dd["task_test"], Dd["rest_post"])
    obs_r, obs_r_s = rho_sym_split(Dr["A"], Dr["B"], Dr["task_test"], Dr["rest_post"])

    N = Ws["A"].shape[0]; n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    rng = np.random.default_rng(BASE_SEED + idx)      # SAME seed/order as audit_150
    surr_d = np.empty(R, float); surr_r = np.empty(R, float)
    for r in range(R):
        Wsurr = {k: shuffle(v, n_swaps, rng) for k, v in Ws.items()}   # one draw -> both arms
        Sd = {k: ultra(m) for k, m in Wsurr.items()}
        Sr = {k: ultra_raw(m) for k, m in Wsurr.items()}
        surr_d[r], _ = rho_sym_split(Sd["A"], Sd["B"], Sd["task_test"], Sd["rest_post"])
        surr_r[r], _ = rho_sym_split(Sr["A"], Sr["B"], Sr["task_test"], Sr["rest_post"])

    # tau-divergence (obs, rest_post tree at scale s vs raw-FC tree): propagator IS a tau-family
    Wp = Ws["rest_post"]; raw_pref = ultra_raw(Wp)
    div = {}
    for s in S_GRID:
        try:
            div[f"treesim_s{s:g}"] = float(spearmanr(ultra_s(Wp, s), raw_pref).statistic)
        except Exception:
            div[f"treesim_s{s:g}"] = float("nan")

    rec = {"patient": pat, "band": band, "N_nodes": int(N), "n_pairs": int(Dd["A"].size),
           "n_surrogates": int(np.isfinite(surr_d).sum())}
    rec.update(_pack(obs_d, obs_d_s, surr_d, "diff"))
    rec.update(_pack(obs_r, obs_r_s, surr_r, "raw"))
    rec.update(div)
    return rec


def cohort_gate(df):
    out = []
    for band in BANDS:
        x = df[df.band == band]
        if x.empty:
            continue
        rec = {"band": band, "n_patients": len(x)}
        for tag in ("diff", "raw"):
            try:
                _, p = wilcoxon(x[f"obs_rho_{tag}"].values - x[f"surr_p50_{tag}"].values,
                                alternative="greater")
            except Exception:
                p = float("nan")
            rec[f"gate_p_{tag}"] = float(p)
            rec[f"obs_median_{tag}"] = float(x[f"obs_rho_{tag}"].median())
            rec[f"n_above_{tag}"] = int((x[f"obs_p_{tag}"] < 0.05).sum())
            rec[f"verdict_{tag}"] = "CLEAR" if p < 0.05 else "fail"
        out.append(rec)
    return pd.DataFrame(out)


def write_readme(gate, df, rt):
    sim_cols = [c for c in df.columns if c.startswith("treesim_s")]
    sim = df[sim_cols].mean()
    L = ["---", "name: diffusion_vs_rawfc_gate",
         "scope: ablation_does_the_trace_gate_need_the_diffusion_propagator",
         "date: 2026-07-11", "status: current", "---", "",
         "# Does the band-selective trace need the diffusion propagator? (audit_174)", "",
         "**Head.** rho_sym cohort matched-strength gate on DIFFUSION cophenetic trees",
         "(D=1/K(tau_min), the pipeline) vs RAW-FC trees (D=1/A, diffusion skipped). Same",
         f"surrogates, same estimator, same Wilcoxon; R={R}, {rt:.0f}s. The diffusion arm",
         "reproduces audit_150 (cross-check). Verdict per band below.", "",
         "| band | diff median | gate p (diff) | verdict | raw median | gate p (raw) | verdict |",
         "|---|---|---|---|---|---|---|"]
    for _, r in gate.iterrows():
        L.append(f"| {r['band']} | {r['obs_median_diff']:+.3f} | {r['gate_p_diff']:.4f} | "
                 f"**{r['verdict_diff']}** | {r['obs_median_raw']:+.3f} | {r['gate_p_raw']:.4f} | "
                 f"**{r['verdict_raw']}** |")
    L += ["", "## tau-divergence (cohort-mean Spearman of diffusion tree at s vs raw-FC tree)",
          "s = tau * lambda_max; s=1 is tau_min. Falling toward 0 as s grows = the propagator",
          "departs from raw-FC (it is a tau-family; raw-FC is a single fixed tree).", ""]
    for c in sim_cols:
        L.append(f"- {c.replace('treesim_s','s=')}: {sim[c]:.3f}")
    L += ["", "## Provenance",
          f"- ablation of audit_150 (canonical rho_sym gate); only ultra() -> ultra_raw() differs.",
          f"- R={R} 4-cycle strength-preserving surrogate, SWAP_FACTOR={SWAP_FACTOR}, same rng/order as audit_150.",
          "- cohort gate: Wilcoxon(obs - surr_p50) one-sided greater over the cohort.",
          "- Build: scripts/01_compute/audit/audit_174_diffusion_vs_rawfc_gate.py",
          "- Does NOT modify audit_150."]
    (OUT / "README.md").write_text("\n".join(L))


def main():
    global R
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=int, default=N_SURROGATES, help="surrogates per cell")
    ap.add_argument("--limit", type=int, default=0, help="limit #cells (0=all) for timing")
    a = ap.parse_args()
    R = a.R
    OUT.mkdir(parents=True, exist_ok=True)
    _ = shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))     # warm numba before fork
    jobs = [(i, p, b) for i, (b, p) in enumerate((b, p) for b in BANDS for p in COHORT)]
    if a.limit:
        jobs = jobs[:a.limit]
    ncpu = min(14, (os.cpu_count() or 4) - 2)
    print(f"[audit_174] {len(jobs)} cells, R={R}, diffusion vs raw-FC, {ncpu} workers", flush=True)
    t0 = time.time(); rows = []
    with Pool(ncpu) as pool:
        for i, r in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if r:
                rows.append(r)
                el = time.time() - t0; eta = el / i * (len(jobs) - i)
                print(f"[{i}/{len(jobs)}] {r['patient']}/{r['band']} "
                      f"diff={r['obs_rho_diff']:+.3f}(p{r['obs_p_diff']:.3f}) "
                      f"raw={r['obs_rho_raw']:+.3f}(p{r['obs_p_raw']:.3f}) "
                      f"sim@s1={r.get('treesim_s1', float('nan')):.2f}  "
                      f"[{el:.0f}s ETA {eta:.0f}s]", flush=True)
    df = pd.DataFrame(rows).sort_values(["band", "patient"])
    df.to_csv(OUT / "per_patient_per_band.csv", index=False)
    gate = cohort_gate(df); gate.to_csv(OUT / "cohort_summary.csv", index=False)
    rt = time.time() - t0
    print(f"\n[audit_174] {len(rows)} cells in {rt:.0f}s")
    print("\n=== GATE: diffusion (pipeline) vs raw-FC (1/A) — matched-strength, rho_sym ===")
    cols = ["band", "obs_median_diff", "gate_p_diff", "verdict_diff",
            "obs_median_raw", "gate_p_raw", "verdict_raw"]
    print(gate[cols].to_string(index=False))
    sim_cols = [c for c in df.columns if c.startswith("treesim_s")]
    print("\n=== tau-divergence: cohort-mean Spearman(diffusion tree at s, raw-FC tree) ===")
    print("  s=1 is tau_min; smaller = propagator departs from raw-FC")
    print(df[sim_cols].mean().to_string())
    if a.R == N_SURROGATES and not a.limit:
        write_readme(gate, df, rt)
    print(f"[audit_174] outputs -> {OUT}")


if __name__ == "__main__":
    main()
