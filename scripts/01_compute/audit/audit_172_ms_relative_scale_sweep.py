#!/usr/bin/env python3
"""audit_172 — MS-relative scale sweep of the cophenetic trace (rho_sym x s x matched-strength).

Fulfils the audit_159-scoped follow-up (2026-07-07_r1-tau-multiscale-reexamination.md
S6): the referee audit_121 lacked. audit_121 mapped the bare trace across the
diffusion scale s = tau*lam_max with only a trace-free PLACEBO; here we run the
canonical rho_sym gate (audit_150) SWEPT over s with a matched-strength null
REBUILT at every s, so the honest quantity is z(s) = (obs - surr_mean)/surr_std,
not raw rho_sym(s).

Question: does a band that FAILS the cohort gate at s=1 (alpha first, then delta,
low_gamma) CLEAR matched-strength at its own coarser scale -- an intrinsic
band<->scale match -- or does the null bump identically and kill the apparent gain?

Scope + 5-point preamble:
    .agents/guides/task-persistence-investigation/2026-07-11_ms-relative-scale-sweep.md

ANCHOR: seeds/swap/recipe are audit_150's exactly (BASE_SEED=20260706 + cell_idx,
R=200, phase order A->B->task->post, 4-cycle +/-delta numba swap, avg-linkage
cophenetic at tau=1/lam_max). Therefore the s=1 column reproduces audit_150's
per-cell obs_rho and surr mean BIT-FOR-BIT; asserted < 1e-9 before any reading.
Any s != 1 movement is real scale physics, not a code difference (audit_121 style).

Optimisation (feedback_optimize_time_and_surface_progress): each matrix is
eigendecomposed ONCE and D_coph(s) = V e^{-(s/lam_max) lam} V^T reused across all
s; numba swap (RNG drawn outside the jitted loop -> bit-identical); Pool over 60
cells; live [i/N] progress. Run `... audit_172_..py pilot` first to time one cell.

Outputs:
    data/audit/tau_sweep_ms_gate/per_patient_per_band_s.csv
    data/audit/tau_sweep_ms_gate/cohort_gate_s.csv
    data/audit/tau_sweep_ms_gate/README.md
"""
from __future__ import annotations
import os, sys, time
import numpy as np, pandas as pd, numba
from multiprocessing import Pool
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import cophenet, linkage
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.workflow.fc import load_fc_matrix

# --- config mirrors audit_150 exactly (anchor); adds the s-grid ------------------
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
N_SURROGATES = 200
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260706                      # == audit_150 (s=1 anchor)
PHASE_ORDER = ("A", "B", "task_test", "rest_post")   # RNG-consumption order (anchor)
# dense near the alpha lead (s~1.5-2.5), out to ~s* (=C-peak ~8-15); stay < collapse
SGRID = np.unique(np.concatenate([[0.5, 0.7, 0.85], np.geomspace(1.0, 12.0, 13)]))
S1_IDX = int(np.argmin(np.abs(SGRID - 1.0)))
VAR_FLOOR = 1e-12

HALVES = CACHE_ROOT / "imcoh_halves_fc"
OUT = ROOT / "data" / "audit" / "tau_sweep_ms_gate"
GATE150 = ROOT / "data" / "audit" / "rho_sym_gate" / "per_patient_per_band.csv"


def load_phase(pat, phase, band):
    """Identical to audit_150.load_phase."""
    if phase in ("A", "B"):
        W = np.load(HALVES / pat / f"{band}_rest_pre_{phase}_imcoh_abs.npy")
    else:
        W = load_fc_matrix(pat, phase, band, fc_method="imcoh_abs")
    W = np.asarray(W, float); np.fill_diagonal(W, 0.0); W = np.clip(W, 0.0, 1.0)
    return 0.5 * (W + W.T)


def eig(W):
    """Combinatorial-Laplacian eigenpairs, once per matrix; reused across all s."""
    deg = W.sum(1)
    ev, V = np.linalg.eigh(np.diag(deg) - W)
    return ev, V


def coph_s(ev, V, s):
    """Cophenetic condensed distances at scale s = tau*lam_max (tau = s/lam_max).
    Bit-identical to audit_150.ultra(W) at s=1 (tau = 1/lam_max)."""
    tau = s / ev[-1]
    rho = (V * np.exp(-tau * ev)) @ V.T
    rho /= np.trace(rho)
    with np.errstate(divide="ignore"):
        T = 1.0 / rho
    T = np.maximum(T, T.T); np.fill_diagonal(T, 0.0)
    fin = np.isfinite(T)
    if not fin.all():
        T = np.where(fin, T, np.nanmax(T[fin]))
    return cophenet(linkage(squareform(T, checks=False), method="average"))


@numba.njit
def _swap_loop(W, s, fr, w_max):
    """Verbatim from audit_150 / audit_63 (4-cycle +/-delta strength-preserving)."""
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
    W = W.copy(); N = W.shape[0]
    s = rng.integers(0, N, size=(n_swaps, 4)).astype(np.int64)
    fr = rng.uniform(0.0, 1.0, size=n_swaps)
    return _swap_loop(W, s, fr, float(w_max))


def rho_sym(DA, DB, Dt, Dp):
    """audit_150.rho_sym_split rho_sym leg (symmetric two-arm split-baseline)."""
    with np.errstate(all="ignore"):
        if (np.nanvar(Dt - DA) < VAR_FLOOR or np.nanvar(Dp - DB) < VAR_FLOOR or
                np.nanvar(Dt - DB) < VAR_FLOOR or np.nanvar(Dp - DA) < VAR_FLOOR):
            return np.nan
        r_ab = spearmanr(Dt - DA, Dp - DB).statistic
        r_ba = spearmanr(Dt - DB, Dp - DA).statistic
    return float(0.5 * (r_ab + r_ba))


def per_cell(job):
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASE_ORDER}
    except Exception:
        return None
    N = Ws["A"].shape[0]
    if any(w.shape != (N, N) for w in Ws.values()):
        return None

    # observed: eigendecompose each phase once, cophenetic at every s
    Eo = {ph: eig(Ws[ph]) for ph in PHASE_ORDER}
    Do = {ph: {si: coph_s(*Eo[ph], s) for si, s in enumerate(SGRID)} for ph in PHASE_ORDER}
    if len({Do["A"][0].size for _ in [0]}) and Do["A"][0].size != N * (N - 1) // 2:
        return None
    obs = np.array([rho_sym(Do["A"][si], Do["B"][si], Do["task_test"][si],
                            Do["rest_post"][si]) for si in range(len(SGRID))])

    # matched-strength surrogates: RNG order == audit_150 -> s=1 anchor bit-exact
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    rng = np.random.default_rng(BASE_SEED + idx)
    surr = np.full((len(SGRID), N_SURROGATES), np.nan)
    for r in range(N_SURROGATES):
        Es = {}
        for ph in PHASE_ORDER:                       # shuffle in the anchor order
            Es[ph] = eig(shuffle(Ws[ph], n_swaps, rng))
        for si, s in enumerate(SGRID):
            Ds = {ph: coph_s(*Es[ph], s) for ph in PHASE_ORDER}
            surr[si, r] = rho_sym(Ds["A"], Ds["B"], Ds["task_test"], Ds["rest_post"])

    rows = []
    for si, s in enumerate(SGRID):
        sr = surr[si][np.isfinite(surr[si])]
        if sr.size == 0:
            continue
        sm, sd = float(np.mean(sr)), float(np.std(sr, ddof=1))
        o = float(obs[si])
        rows.append(dict(
            patient=pat, band=band, s_idx=si, s=float(s), N_nodes=int(N),
            n_surr=int(sr.size), obs_rho=o,
            surr_mean=sm, surr_std=sd, surr_p50=float(np.median(sr)),
            z=(o - sm) / sd if (sd > 0 and np.isfinite(o)) else np.nan,
            p_one_sided=float(np.mean(sr >= o)) if np.isfinite(o) else np.nan,
            obs_finite=int(np.isfinite(o)),
        ))
    return rows


def anchor_check(df):
    """s=1 column must reproduce audit_150 per-cell obs_rho + surr mean bit-exactly."""
    if not GATE150.exists():
        print(f"[anchor] WARN {GATE150} missing -> cannot verify s=1 anchor", flush=True)
        return
    ref = pd.read_csv(GATE150).set_index(["patient", "band"])
    s1 = df[df.s_idx == S1_IDX].set_index(["patient", "band"])
    worst_o = worst_s = 0.0
    for k, row in s1.iterrows():
        if k not in ref.index:
            continue
        do = abs(row.obs_rho - ref.loc[k, "obs_rho"])
        ds = abs(row.surr_mean - ref.loc[k, "surr_mean_rho"])
        worst_o = max(worst_o, do); worst_s = max(worst_s, ds)
    print(f"[anchor] s=1 vs audit_150: max|Dobs|={worst_o:.2e}  max|Dsurr_mean|={worst_s:.2e}",
          flush=True)
    if worst_o > 1e-9 or worst_s > 1e-9:
        raise AssertionError("s=1 anchor mismatch vs audit_150 -> recipe drift, ABORT")
    print("[anchor] PASS — s=1 reproduces the certified gate; s!=1 is real scale physics.",
          flush=True)


def cohort_gate(df):
    out = []
    for band in BANDS:
        for si, s in enumerate(SGRID):
            x = df[(df.band == band) & (df.s_idx == si)]
            if len(x) < 3:
                continue
            d = (x.obs_rho.values - x.surr_p50.values)
            d = d[np.isfinite(d)]
            try:
                p = float(wilcoxon(d, alternative="greater").pvalue) if d.size >= 3 else np.nan
            except Exception:
                p = np.nan
            out.append(dict(
                band=band, s_idx=int(si), s=float(s), n_patients=int(len(x)),
                obs_median=float(np.nanmedian(x.obs_rho)),
                surr_median=float(np.nanmedian(x.surr_p50)),
                z_median=float(np.nanmedian(x.z)),
                n_above=int((x.p_one_sided < 0.05).sum()),
                gate_p_sym=p, verdict="CLEAR" if (np.isfinite(p) and p < 0.05) else "fail",
            ))
    return pd.DataFrame(out)


def run_pilot():
    _ = shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))  # warm numba
    job = (COHORT.index("Pat_08") + BANDS.index("alpha") * len(COHORT), "Pat_08", "alpha")
    print(f"[pilot] one cell {job[1]}/{job[2]}, R={N_SURROGATES}, |s|={len(SGRID)} ...", flush=True)
    t0 = time.time(); rows = per_cell(job); dt = time.time() - t0
    if rows:
        r = pd.DataFrame(rows)
        print(r[["s", "obs_rho", "surr_mean", "z", "p_one_sided"]].to_string(index=False))
    print(f"\n[pilot] {dt:.1f}s / cell  ->  60 cells serial ~{dt*60/60:.1f} min; "
          f"Pool(14) ~{dt*60/14/60:.1f} min", flush=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    _ = shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))  # warm numba before fork
    jobs = [(i, p, b) for i, (b, p) in
            enumerate((b, p) for b in BANDS for p in COHORT)]   # idx order == audit_150
    ncpu = min(14, (os.cpu_count() or 4) - 2)
    print(f"[audit_172] {len(jobs)} cells, R={N_SURROGATES}, |s|={len(SGRID)} "
          f"(s in [{SGRID[0]:.2f},{SGRID[-1]:.1f}]), {ncpu} workers", flush=True)
    t0 = time.time(); rows = []
    with Pool(ncpu) as pool:
        for i, r in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if r:
                rows.extend(r)
                rd = pd.DataFrame(r)
                s1 = rd[rd.s_idx == S1_IDX].iloc[0]
                best = rd.loc[rd.z.idxmax()] if rd.z.notna().any() else s1
                print(f"[{i}/{len(jobs)}] {r[0]['patient']}/{r[0]['band']} "
                      f"s=1 obs={s1.obs_rho:+.3f} z={s1.z:+.2f} | "
                      f"peak z={best.z:+.2f}@s={best.s:.2f}", flush=True)
    df = pd.DataFrame(rows).sort_values(["band", "patient", "s_idx"])
    df.to_csv(OUT / "per_patient_per_band_s.csv", index=False)
    anchor_check(df)
    gate = cohort_gate(df); gate.to_csv(OUT / "cohort_gate_s.csv", index=False)
    rt = time.time() - t0
    print(f"\n[audit_172] {len(df)} rows in {rt:.0f}s -> {OUT}", flush=True)

    # headline: per band, s=1 vs the best collapse-free s (gate_p minimum)
    print("\n=== COHORT GATE gate_p_sym(s): s=1 vs best s (min p) ===")
    print(" band         gate_p@s=1  verdict@1 | best_p   @s     verdict")
    for band in BANDS:
        g = gate[gate.band == band]
        if g.empty:
            continue
        g1 = g[g.s_idx == S1_IDX].iloc[0]
        gb = g.loc[g.gate_p_sym.idxmin()] if g.gate_p_sym.notna().any() else g1
        print(f" {band:<11}  {g1.gate_p_sym:.4f}    {g1.verdict:<5}   | "
              f"{gb.gate_p_sym:.4f}  {gb.s:5.2f}   {gb.verdict}", flush=True)
    write_readme(gate, rt)


def write_readme(gate, rt):
    L = ["---", "name: tau_sweep_ms_gate",
         "scope: ms_relative_scale_sweep_of_cophenetic_trace_rhosym",
         "date: 2026-07-11", "status: current",
         "implements: audit_172 (audit_159-scoped follow-up)", "---", "",
         "# MS-relative scale sweep of the cophenetic trace (rho_sym x s x matched-strength)",
         "",
         "**Head.** audit_150's rho_sym gate swept over the diffusion scale",
         "s = tau*lam_max with a matched-strength null rebuilt at each s. s=1 anchors",
         f"the certified gate bit-exactly. R={N_SURROGATES}. Runtime {rt:.0f}s.",
         "Question: does a band failing at s=1 clear MS at its own coarser scale?", "",
         "| band | gate_p @ s=1 | best gate_p | @ s | verdict |",
         "|---|---|---|---|---|"]
    for band in BANDS:
        g = gate[gate.band == band]
        if g.empty:
            continue
        g1 = g[g.s_idx == S1_IDX].iloc[0]
        gb = g.loc[g.gate_p_sym.idxmin()] if g.gate_p_sym.notna().any() else g1
        L.append(f"| {band} | {g1.gate_p_sym:.4f} | {gb.gate_p_sym:.4f} | "
                 f"{gb.s:.2f} | **{gb.verdict}** |")
    L += ["", "## Provenance",
          f"- s-grid: {np.array2string(SGRID, precision=2, separator=', ')}",
          f"- estimator rho_sym; MS 4-cycle +/-delta, R={N_SURROGATES}, SWAP_FACTOR={SWAP_FACTOR}.",
          "- s=1 anchored bit-exact to audit_150 (BASE_SEED+cell_idx, same swap order).",
          "- cophenetic reformed at each s from cached eigenpairs (V e^{-(s/lmax)lam} V^T).",
          "- cohort gate: Wilcoxon(obs - surr_p50) one-sided greater over 10 patients.",
          "- Build: scripts/01_compute/audit/audit_172_ms_relative_scale_sweep.py",
          "- Scope: .agents/guides/task-persistence-investigation/2026-07-11_ms-relative-scale-sweep.md",
          "- Does NOT modify audit_121/150; new pipeline, new output dir."]
    (OUT / "README.md").write_text("\n".join(L))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "pilot":
        run_pilot()
    else:
        main()
