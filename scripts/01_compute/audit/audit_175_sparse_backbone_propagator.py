#!/usr/bin/env python3
"""audit_175 — does a parameter-free sparse backbone restore the LRG propagator?

Sequel to audit_174 (which found that on the FULLY-CONNECTED |ImCoh| graph at
tau=1/lambda_max the diffusion cophenetic tree D=1/K is ~identical to the raw-FC
tree D=1/A: alpha p.024, beta p.032 CLEAR both). Villegas 2025 (PRR 7 013065)
Appendix F: high mean-connectivity networks have a Wigner-semicircle spectrum ->
a SINGLE C(tau) collapse peak, no multiscale ladder. Our fully-connected FC is
exactly that degenerate regime (guide section 6). LRG is designed for sparse /
topological graphs, where e^{-tau L} integrates ALL powers of L (all paths) that
a single adjacency scale cannot. This audit sparsifies to a PARAMETER-FREE
connected backbone and asks whether the propagator becomes load-bearing there.

Backbones (all spanning => keep all N nodes, connected, no giant-component
gymnastics; rho_sym per-pair alignment preserved for free):
  - maximum spanning tree (MST)     : minimal parameter-free connected backbone
  - TMFG (Massara 2016)             : 3N-6 edges, chordal, the canonical
                                      parameter-free multiscale backbone
  - MST union top-fraction grid     : robustness sweep over density; frac=1.0
                                      recovers the fully-connected ANCHOR.

Arms (per phase, on each backbone), UPGMA + cophenetic, rho_sym cross-phase:
  - diff : LRG diffusion, D=1/rho(tau_min), rho=e^{-tau L}/Tr (== audit_150 ultra
           at frac=1.0; underflow floored at +1e-30, K>=0 exactly -- Metzler).
  - geo  : graph geodesic (shortest-path on edge length 1/A) -- a NON-diffusion
           multi-step baseline. The SHARP control: if diff~=geo across the grid,
           topology is needed but shortest-path suffices (diffusion not uniquely
           necessary); if diff separates trace-from-null better, the all-paths
           integration is specifically valuable.
  - raw  : D=1/A single scale. At frac=1.0 == audit_174 raw arm (ANCHOR). On a
           backbone most pairs have no edge -> 1/A=inf -> degenerate: the
           single-scale control that BREAKS on a topological graph (reported,
           not a bug).

Null: matched-strength 4-cycle surrogate (audit_150 numba shuffle), SPARSIFIED
THE SAME WAY at the same backbone -- because sparsification imposes topology,
the null is MORE load-bearing here, not less (it separates a reproducible
data-multiscale trace from a sparsification artifact). Same rng/order as
audit_150/174 so the frac=1.0 diff arm reproduces them bit-for-bit.

5-point critical preamble
1. Claim: on a parameter-free backbone the diffusion propagator does work the
   raw single scale cannot (diff-tree departs from raw-tree AND from geodesic),
   a multiscale C(tau) ladder emerges (>1 peak vs 1 fully-connected), and the
   band-selective trace (alpha/beta) SURVIVES the matched-strength gate.
2. Null: sparsification does not change the picture (diff~=geo~=raw on the
   backbone too), OR sparsification destroys the trace (beta fails at every
   parameter-free backbone). Either => the diffusion framing is not rescued.
3. Strongest alternatives controlled: (a) p-hacking the density -> report every
   quantity vs density on a pre-committed grid AND at parameter-free MST/TMFG,
   effect must be robust not a magic point; (b) raw-FC ill-defined on sparse
   graphs is a weak win -> the sharp baseline is GEODESIC (well-defined,
   multi-step, non-diffusion); (c) sparsification removes trace-carrying edges
   -> report beta obs rho_sym vs backbone, must not collapse; (d) shared-draw
   artifact -> identical surrogate draws feed every arm; frac=1.0 diff
   reproduces audit_174 bit-for-bit (anchor).
4. Cannot: a rescued backbone does NOT validate the fully-connected pipeline
   (still degenerate) -- it motivates MIGRATING downstream results
   (beta->OFC, enc/inf, epi) to the backbone. Filter choice is a modeling
   choice -> robustness across MST/TMFG/grid required.
5. Falsify: if across the whole grid diff~=geo AND beta clears identically in
   both, the propagator is inessential even on the backbone. If beta dies at
   all reasonable backbones, the trace lives in the weak-edge bulk.

Outputs:
    data/audit/sparse_backbone_propagator/per_cell.csv
    data/audit/sparse_backbone_propagator/cohort_gate.csv   (band x backbone x arm)
    data/audit/sparse_backbone_propagator/README.md
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
from audit_150_rho_sym_gate import (            # reuse canonical gate machinery (DRY)
    load_phase, ultra, shuffle, rho_sym_split,
    COHORT, BANDS, N_SURROGATES, SWAP_FACTOR, BASE_SEED,
)
from lrg_eegfc.utils.fc.backbone import (
    maximum_spanning_tree, tmfg_backbone, mst_union_top_fraction, geodesic_distance,
)
from lrg_eegfc.workflow.diagnostics import compute_susceptibility_diagnostics

OUT = ROOT / "data" / "audit" / "sparse_backbone_propagator"
RHO_FLOOR = 1e-30           # heat-kernel underflow floor (K>=0 exactly; far pairs)
PHASES = ("A", "B", "task_test", "rest_post")

# Parameter-free backbones + a density grid for robustness. frac=1.0 == anchor.
BACKBONES: dict = {
    "full":   lambda A: A,                                 # fully-connected anchor
    "g0.35":  lambda A: mst_union_top_fraction(A, 0.35),
    "g0.20":  lambda A: mst_union_top_fraction(A, 0.20),
    "g0.10":  lambda A: mst_union_top_fraction(A, 0.10),
    "g0.05":  lambda A: mst_union_top_fraction(A, 0.05),
    "tmfg":   tmfg_backbone,
    "mst":    maximum_spanning_tree,
}
R = N_SURROGATES            # overridable via --R (set in main before Pool fork)


def diff_coph(W):
    """LRG diffusion cophenetic on a (possibly sparse) backbone.

    Identical to audit_150.ultra() on the full graph (RHO_FLOOR never triggers
    there); on a backbone it floors underflowed/roundoff rho at +RHO_FLOOR so
    far (K->0) pairs map to a large FINITE positive distance. K=e^{-tau L}>=0
    exactly (Metzler); rho_sym is Spearman-rank so the floor value is immaterial.
    """
    deg = W.sum(1); ev, V = np.linalg.eigh(np.diag(deg) - W); tau = 1.0 / ev[-1]
    rho = (V * np.exp(-tau * ev)) @ V.T; rho /= np.trace(rho)
    rho = np.where(rho > RHO_FLOOR, rho, RHO_FLOOR)
    T = 1.0 / rho
    np.fill_diagonal(T, 0.0); T = np.maximum(T, T.T)
    return cophenet(linkage(squareform(T, checks=False), method="average"))


def geo_coph(W):
    """Geodesic (shortest-path) cophenetic -- non-diffusion multi-step baseline."""
    D = geodesic_distance(W, mode="inverse")
    return cophenet(linkage(squareform(D, checks=False), method="average"))


def raw_coph(W):
    """Single-scale D=1/A cophenetic (== audit_174 ultra_raw on the full graph).

    On a sparse backbone most pairs have no edge -> inf -> clamped to the max
    finite distance: the control that degenerates on a topological graph.
    """
    with np.errstate(divide="ignore"):
        T = 1.0 / W
    np.fill_diagonal(T, 0.0); T = np.maximum(T, T.T); fin = np.isfinite(T)
    if fin.sum() < 3:
        return None
    if not fin.all():
        T = np.where(fin, T, np.nanmax(T[fin]))
    return cophenet(linkage(squareform(T, checks=False), method="average"))


def _npeaks(W):
    """C(tau) interior-peak count on the backbone Laplacian (regime restoration).

    Villegas 2025 Appendix A/F: 1 peak = collapsed single scale (semicircle);
    >1 = multiscale ladder. NOT a cross-phase trace metric (L2 rung retired)."""
    deg = W.sum(1); ev = np.linalg.eigvalsh(np.diag(deg) - W)
    ev = np.maximum(ev, 0.0)
    try:
        return int(compute_susceptibility_diagnostics(ev)["N_peaks"])
    except Exception:
        return -1


def _gate_pack(obs, surr, tag):
    s = surr[np.isfinite(surr)]
    if s.size == 0:
        return {f"obs_{tag}": obs, f"p_{tag}": float("nan"), f"surr_p50_{tag}": float("nan")}
    return {f"obs_{tag}": obs, f"surr_p50_{tag}": float(np.quantile(s, 0.50)),
            f"p_{tag}": float(np.mean(s >= obs))}


def per_cell(job):
    idx, pat, band = job
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES}
    except Exception:
        return None
    N = Ws["A"].shape[0]
    if any(v.shape[0] != N for v in Ws.values()):
        return None

    # observed per backbone
    rec = {"patient": pat, "band": band, "N_nodes": int(N)}
    obs = {}   # (bb) -> (rho_diff, rho_geo)
    for bb, filt in BACKBONES.items():
        Wb = {ph: filt(Ws[ph]) for ph in PHASES}
        Dd = {ph: diff_coph(Wb[ph]) for ph in PHASES}
        Dg = {ph: geo_coph(Wb[ph]) for ph in PHASES}
        rho_d, _ = rho_sym_split(Dd["A"], Dd["B"], Dd["task_test"], Dd["rest_post"])
        rho_g, _ = rho_sym_split(Dg["A"], Dg["B"], Dg["task_test"], Dg["rest_post"])
        obs[bb] = (rho_d, rho_g)
        # tree divergence + regime restoration (observed rest_post / rest_pre)
        rec[f"treesim_dg_{bb}"] = float(spearmanr(Dd["rest_post"], Dg["rest_post"]).statistic)
        rec[f"npeaks_{bb}"] = _npeaks(Wb["A"])
        rec[f"density_{bb}"] = float((np.triu(Wb["A"], 1) > 0).mean())
        if bb == "full":   # anchor cross-checks vs audit_174
            cr = raw_coph(Wb["rest_post"])
            rec["treesim_dr_full"] = (float(spearmanr(Dd["rest_post"], cr).statistic)
                                      if cr is not None else float("nan"))
            rho_r, _ = rho_sym_split(*[raw_coph(Wb[ph]) for ph in PHASES]) \
                if all(raw_coph(Wb[ph]) is not None for ph in PHASES) else (float("nan"), None)
            rec["obs_raw_full"] = rho_r

    # matched-strength null: one surrogate draw per r feeds ALL backbones+arms.
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    rng = np.random.default_rng(BASE_SEED + idx)     # SAME seed/order as audit_150/174
    surr = {bb: (np.empty(R), np.empty(R)) for bb in BACKBONES}
    for r in range(R):
        Wsurr = {ph: shuffle(Ws[ph], n_swaps, rng) for ph in PHASES}   # order preserved
        for bb, filt in BACKBONES.items():
            Wb = {ph: filt(Wsurr[ph]) for ph in PHASES}
            Dd = {ph: diff_coph(Wb[ph]) for ph in PHASES}
            Dg = {ph: geo_coph(Wb[ph]) for ph in PHASES}
            surr[bb][0][r], _ = rho_sym_split(Dd["A"], Dd["B"], Dd["task_test"], Dd["rest_post"])
            surr[bb][1][r], _ = rho_sym_split(Dg["A"], Dg["B"], Dg["task_test"], Dg["rest_post"])

    for bb in BACKBONES:
        rec.update(_gate_pack(obs[bb][0], surr[bb][0], f"diff_{bb}"))
        rec.update(_gate_pack(obs[bb][1], surr[bb][1], f"geo_{bb}"))
    return rec


def cohort_gate(df):
    out = []
    for band in BANDS:
        x = df[df.band == band]
        if x.empty:
            continue
        for bb in BACKBONES:
            rec = {"band": band, "backbone": bb, "n_patients": len(x),
                   "density": float(x[f"density_{bb}"].median()),
                   "npeaks_median": float(x[f"npeaks_{bb}"].median()),
                   "treesim_diff_geo": float(x[f"treesim_dg_{bb}"].median())}
            for arm in ("diff", "geo"):
                col, p50 = f"obs_{arm}_{bb}", f"surr_p50_{arm}_{bb}"
                try:
                    _, p = wilcoxon(x[col].values - x[p50].values, alternative="greater")
                except Exception:
                    p = float("nan")
                rec[f"gate_p_{arm}"] = float(p)
                rec[f"obs_median_{arm}"] = float(x[col].median())
                rec[f"n_above_{arm}"] = int((x[f"p_{arm}_{bb}"] < 0.05).sum())
                rec[f"verdict_{arm}"] = "CLEAR" if p < 0.05 else "fail"
            out.append(rec)
    return pd.DataFrame(out)


def write_readme(gate, df, rt):
    L = ["---", "name: sparse_backbone_propagator",
         "scope: does_a_parameter_free_sparse_backbone_make_the_LRG_propagator_load_bearing",
         "date: 2026-07-11", "status: current", "---", "",
         "# Sparse-backbone propagator recovery (audit_175)", "",
         "**Head.** rho_sym cohort matched-strength gate on the DIFFUSION vs GEODESIC",
         "cophenetic tree across parameter-free backbones (MST, TMFG) and an MST-union",
         f"density grid. Sequel to audit_174; R={R}, {rt:.0f}s. frac=1.0 diff arm",
         "reproduces audit_150/174 (anchor). Regime restoration = C(tau) peak count.", "",
         "| band | backbone | density | C(tau) peaks | diff median | gate p (diff) | geo median | gate p (geo) | sp(diff,geo) |",
         "|---|---|---|---|---|---|---|---|---|"]
    for _, r in gate.iterrows():
        L.append(f"| {r['band']} | {r['backbone']} | {r['density']:.3f} | {r['npeaks_median']:.0f} "
                 f"| {r['obs_median_diff']:+.3f} | {r['gate_p_diff']:.4f} "
                 f"| {r['obs_median_geo']:+.3f} | {r['gate_p_geo']:.4f} | {r['treesim_diff_geo']:.3f} |")
    L += ["", "## Anchors (frac=1.0 == audit_174)",
          f"- sp(diff-tree, raw 1/A tree) full = {df['treesim_dr_full'].mean():.3f} (audit_174 ~0.83-0.93)",
          f"- C(tau) peaks full (cohort median) = {df[[c for c in df if c=='npeaks_full']].median().values}",
          "", "## Provenance",
          f"- sequel to audit_174; reuses audit_150 load_phase/ultra/shuffle/rho_sym + numba surrogate.",
          f"- backbones from lrg_eegfc.utils.fc.backbone (parameter-free, spanning, connected).",
          f"- R={R} matched-strength surrogate SPARSIFIED the same way; same rng/order as audit_150.",
          "- Build: scripts/01_compute/audit/audit_175_sparse_backbone_propagator.py"]
    (OUT / "README.md").write_text("\n".join(L))


def main():
    global R
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=int, default=N_SURROGATES)
    ap.add_argument("--limit", type=int, default=0, help="limit #cells (timing)")
    ap.add_argument("--workers", type=int, default=0)
    a = ap.parse_args()
    R = a.R
    OUT.mkdir(parents=True, exist_ok=True)
    _ = shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))   # warm numba shuffle
    _ = tmfg_backbone(np.ones((6, 6)) - np.eye(6))               # warm numba tmfg
    jobs = [(i, p, b) for i, (b, p) in enumerate((b, p) for b in BANDS for p in COHORT)]
    if a.limit:
        jobs = jobs[:a.limit]
    ncpu = a.workers or min(14, (os.cpu_count() or 4) - 2)
    print(f"[audit_175] {len(jobs)} cells, R={R}, {len(BACKBONES)} backbones x (diff,geo), "
          f"{ncpu} workers", flush=True)
    t0 = time.time(); rows = []
    with Pool(ncpu) as pool:
        for i, r in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            if r:
                rows.append(r)
                el = time.time() - t0; eta = el / i * (len(jobs) - i)
                print(f"[{i}/{len(jobs)}] {r['patient']}/{r['band']} "
                      f"tmfg diff={r.get('obs_diff_tmfg', float('nan')):+.3f}"
                      f"(p{r.get('p_diff_tmfg', float('nan')):.3f}) "
                      f"geo={r.get('obs_geo_tmfg', float('nan')):+.3f}"
                      f"(p{r.get('p_geo_tmfg', float('nan')):.3f}) "
                      f"peaks[full/tmfg/mst]={r.get('npeaks_full')}/{r.get('npeaks_tmfg')}/{r.get('npeaks_mst')} "
                      f"[{el:.0f}s ETA {eta:.0f}s]", flush=True)
    df = pd.DataFrame(rows).sort_values(["band", "patient"])
    df.to_csv(OUT / "per_cell.csv", index=False)
    gate = cohort_gate(df); gate.to_csv(OUT / "cohort_gate.csv", index=False)
    rt = time.time() - t0
    print(f"\n[audit_175] {len(rows)} cells in {rt:.0f}s")
    print("\n=== GATE: diffusion vs geodesic, per band x backbone (matched-strength) ===")
    show = gate[["band", "backbone", "density", "npeaks_median", "obs_median_diff",
                 "gate_p_diff", "verdict_diff", "gate_p_geo", "verdict_geo", "treesim_diff_geo"]]
    print(show.to_string(index=False))
    if a.R == N_SURROGATES and not a.limit:
        write_readme(gate, df, rt)
    print(f"[audit_175] outputs -> {OUT}")


if __name__ == "__main__":
    main()
