#!/usr/bin/env python3
"""C1 — localization cohort-STATISTIC bake-off on the mst@0.20 backbone.

The "beta is delocalized on mst@0.20" verdict (scripts 16/17) and the whole-graph
"beta -> OFC" verdict (audit_83/151/158/160) rest on OPPOSITE broken cohort
statistics: scripts 16/17 use a per-patient Wilcoxon signed-rank across the 5-8
patients that sample each system (floor 1/2^K -> DOA at K=5), while the whole-graph
era used a pooled cohort-median-vs-R-surrogate-medians (floor 1/(R+1), liberal). This
script holds the substrate FIXED (mst@0.20, IDENTICAL machinery + matched-strength
null as script 17) and varies ONLY the cohort statistic, computing ALL FOUR side by
side per (target, band, scale, granularity, unit):

    (a) pooled_median   (b) wilcoxon   (c) stouffer   (d) sign_consistency

plus whole-grid AND per-cell BH, leave-one-patient-out worst-p, frac_pos, and K
coverage. No single statistic is the gate (user decision 2026-07-13): the reader
adjudicates by convergence. Three granularities so the coverage-fair supersystem
'limbic' (K=10, no Wilcoxon floor) is tested alongside 'system' and 'soz'.

NO series is ever truncated (locked 2026-07-13). Duration confound -> ratio
regression (script C3), not length-matching.

5-point preamble
1. Claim: on mst@0.20 the trace/encoding/inference cophenetic concordance concentrates
   on anatomical systems above matched-strength, at some scale.
2. Null: matched-strength R=200, IDENTICAL to scripts 13/16/17 (4-cycle +/-delta over
   every phase -> sparsify@0.20 -> LRG). Preserves node strength.
3. Strongest alternative: (a) node strength (controlled by MS); (b) THE COHORT
   STATISTIC ITSELF (the real adversary — not controlled by any null; this script is
   the control).
4. Cannot: the MS null says nothing about cohort aggregation. Whole-grid BH is the only
   honest multiplicity gate over 3 targets x 4 bands x 16 scales x systems. Coverage
   floors are structural, not effects.
5. Falsify: if under EVERY statistic (incl. the powered, consistency-favouring (d)) at
   supersystem granularity (limbic K=10) with whole-grid BH + LOO nothing concentrates,
   the trace is genuinely delocalized (a result), established by a valid statistic.

Two stages (resumable):
  Stage A  per (patient, band): dump obs + R-surrogate per-unit-mean arrays ->
           data/sparsified_arc/localization_bakeoff_mst020/arrays/{pat}_{band}.npz
  Stage B  cohort-combine across patients -> {system,supersystem,soz}_bakeoff.csv
Run:  [SA_FRAC=0.20] [SA_LIMIT=n] python 21_localization_statistic_bakeoff_mst020.py [--stage a|b|ab]
"""
from __future__ import annotations
import os, sys, time, argparse
import numpy as np, pandas as pd
from multiprocessing import Pool

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "sparsified_arc"))
from audit_150_rho_sym_gate import COHORT, load_phase                 # noqa: E402
from lrg_eegfc.utils.fc.backbone import mst_union_top_fraction         # noqa: E402
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, cophenetic_at_scale  # noqa: E402
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle  # noqa: E402
from lrg_eegfc.utils.metrics.node_localization import epi_keep_mask     # noqa: E402
from lrg_eegfc.utils.io.regions import load_channel_regions            # noqa: E402
from lrg_eegfc.utils.metrics import cohort_localization as CL          # noqa: E402
# concordance primitives (trace / encoding / inference) — reuse script 17 verbatim
from importlib import import_module                                    # noqa: E402
_m17 = import_module("17_localization_arc_mst020")
concordances_at_scale = _m17.concordances_at_scale
unit_means = _m17.unit_means

PHASES = ("A", "B", "task_learn", "task_test", "rest_post")
BANDS = ["delta", "alpha", "beta", "low_gamma"]
TARGETS = ("trace", "encoding", "inference")
# SA_GRANLADDER=1 -> scale-matched parcellation ladder (coarse->fine) + soz. The
# diffusion at scale s coarse-grains space; a coarse-scale trace must be tested at a
# coarse parcellation (hemisphere) and a fine one at a fine atlas (region). Testing a
# fixed 9-system atlas at every scale is a granularity mismatch (the "delocalized"
# read may just be the wrong resolution). Order coarse->fine: hemisphere < lobe <
# supersystem < system < region ; soz is an orthogonal binary tissue split.
GRAN_LADDER = os.environ.get("SA_GRANLADDER", "0") == "1"
# ladder rungs coarse->fine; `region` (DK) dropped — per-region coverage is 1–3 patients
# (below the min-3 cohort floor) so it adds cost without a testable cell.
GROUPINGS = (("hemisphere", "lobe", "supersystem", "system", "soz")
             if GRAN_LADDER else ("system", "supersystem", "soz"))
FRAC = float(os.environ.get("SA_FRAC", "0.20"))
SGRID = np.logspace(0.0, np.log10(180.0), 16)
R = 200
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260713
DROP = {"system": {"non_anatomical"},
        "supersystem": {"non_anatomical", "other"},
        "soz": set(),
        "hemisphere": {"unknown", "?", "nan", "none"},
        "lobe": {"unknown", "?", "nan", "none", "white_matter", "other"},
        "region": {"unknown", "Wm", "Unk", "nan", "Cerebellum-Cortex"}}
LIMIT = int(os.environ.get("SA_LIMIT", "0"))
_SUFFIX = ("_glad" if GRAN_LADDER else "") + ("" if FRAC == 0.20 else f"_frac{FRAC:0.2f}")
OUT = ROOT / "data" / "sparsified_arc" / f"localization_bakeoff_mst020{_SUFFIX}"
ARR = OUT / "arrays"


def eig_backbone(W):
    return laplacian_eig(mst_union_top_fraction(W, FRAC))


def label_vectors(rdf, pat, N):
    """node->unit label vectors (all ladder granularities + soz), aligned to FC order."""
    soz_vec = np.where(epi_keep_mask(rdf, pat), "non_soz", "soz")
    out = {"soz": soz_vec}
    for g in ("hemisphere", "lobe", "supersystem", "system", "region"):
        if g in rdf.columns:
            out[g] = rdf[g].to_numpy().astype(str)
    return {g: out[g] for g in GROUPINGS if g in out}


def per_cell(job):
    """One (patient, band): dump obs + surrogate per-unit-means for all
    (target, grouping, scale, unit) to an npz. Returns (pat, band, N, path|None)."""
    idx, pat, band = job
    tag = f"{pat}_{band}"
    outp = ARR / f"{tag}.npz"
    if outp.exists():
        return (pat, band, -1, str(outp))          # already done (resume)
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES}
    except Exception:
        return (pat, band, 0, None)
    N = Ws["A"].shape[0]
    if any(v.shape[0] != N for v in Ws.values()):
        return (pat, band, 0, None)
    rdf = load_channel_regions(pat)
    if len(rdf) != N:
        print(f"  [skip] {pat} {band}: regions {len(rdf)} != FC {N}", flush=True)
        return (pat, band, 0, None)
    labs = label_vectors(rdf, pat, N)
    iu = np.triu_indices(N, 1); ii, jj = iu[0].astype(int), iu[1].astype(int)
    NS = len(SGRID)

    def unit_means_at(eig5):
        """{target: {grouping: {scale_idx: {unit: mean}}}} for one FC set."""
        res = {t: {g: [None] * NS for g in GROUPINGS} for t in TARGETS}
        for si, s in enumerate(SGRID):
            C = {ph: cophenetic_at_scale(*eig5[ph], s) for ph in PHASES}
            conc = concordances_at_scale(C)
            for t in TARGETS:
                for g in GROUPINGS:
                    res[t][g][si] = unit_means(conc[t], ii, jj, labs[g], drop=DROP[g])
        return res

    eig_obs = {ph: eig_backbone(Ws[ph]) for ph in PHASES}
    obs = unit_means_at(eig_obs)

    # unit universe per grouping (from observed, sorted)
    units = {g: sorted({u for t in TARGETS for si in range(NS)
                        for u in obs[t][g][si]}) for g in GROUPINGS}
    uidx = {g: {u: k for k, u in enumerate(units[g])} for g in GROUPINGS}

    # obs tensors (target,g): (NS, U) ; surr tensors: (NS, U, R) NaN-filled
    obs_T = {(t, g): np.full((NS, len(units[g])), np.nan) for t in TARGETS for g in GROUPINGS}
    sur_T = {(t, g): np.full((NS, len(units[g]), R), np.nan) for t in TARGETS for g in GROUPINGS}
    for t in TARGETS:
        for g in GROUPINGS:
            for si in range(NS):
                for u, v in obs[t][g][si].items():
                    obs_T[(t, g)][si, uidx[g][u]] = v

    rng = np.random.default_rng(BASE_SEED + idx)
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    for r in range(R):
        eig_sh = {ph: eig_backbone(matched_strength_shuffle(Ws[ph], n_swaps, rng, W_MAX))
                  for ph in PHASES}
        try:
            sm = unit_means_at(eig_sh)
        except Exception:
            continue
        for t in TARGETS:
            for g in GROUPINGS:
                for si in range(NS):
                    for u, v in sm[t][g][si].items():
                        sur_T[(t, g)][si, uidx[g][u], r] = v

    save = {"scales": SGRID.astype(np.float32), "N": np.int32(N)}
    for g in GROUPINGS:
        save[f"{g}__units"] = np.array(units[g], dtype=object)
    for t in TARGETS:
        for g in GROUPINGS:
            save[f"{t}__{g}__obs"] = obs_T[(t, g)].astype(np.float32)
            save[f"{t}__{g}__surr"] = sur_T[(t, g)].astype(np.float32)
    ARR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(outp, **save)
    return (pat, band, N, str(outp))


def stage_a():
    ARR.mkdir(parents=True, exist_ok=True)
    _ = matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))  # warm numba
    jobs = [(i, pat, band) for i, (pat, band) in
            enumerate((p, b) for b in BANDS for p in COHORT)]
    if LIMIT:
        jobs = jobs[:LIMIT]
    print(f"[C1 stage A] frac={FRAC} R={R} {len(SGRID)} scales {len(jobs)} (pat,band) cells "
          f"-> {ARR}", flush=True)
    t0 = time.time()
    nproc = min(12, max(1, (os.cpu_count() or 2) - 2))
    with Pool(nproc) as pool:
        for k, (pat, band, N, path) in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            el = time.time() - t0
            stat = "cached" if N == -1 else (f"N={N}" if path else "SKIP")
            print(f"  [{k}/{len(jobs)}] {pat} {band} {stat}  elapsed {el:6.0f}s "
                  f"ETA {el/k*(len(jobs)-k):6.0f}s", flush=True)
    print(f"[C1 stage A done] {time.time()-t0:.0f}s", flush=True)


def stage_b():
    """Cohort-combine: for each (target,band,scale,grouping,unit) gather (obs, surr[R])
    across patients -> all 4 statistics + LOO + BH (whole-grid per band + per-cell)."""
    npzs = {}
    for band in BANDS:
        for pat in COHORT:
            p = ARR / f"{pat}_{band}.npz"
            if p.exists():
                npzs[(pat, band)] = np.load(p, allow_pickle=True)
    print(f"[C1 stage B] loaded {len(npzs)} (pat,band) arrays", flush=True)
    rows = []
    for g in GROUPINGS:
        for band in BANDS:
            pats = [pat for pat in COHORT if (pat, band) in npzs]
            for t in TARGETS:
                for si in range(len(SGRID)):
                    # gather per-unit across patients
                    per_unit = {}   # unit -> list of (obs, surr[R])
                    for pat in pats:
                        d = npzs[(pat, band)]
                        units = list(d[f"{g}__units"])
                        obsM = d[f"{t}__{g}__obs"]; surM = d[f"{t}__{g}__surr"]
                        for k, u in enumerate(units):
                            o = float(obsM[si, k])
                            if not np.isfinite(o):
                                continue
                            per_unit.setdefault(u, []).append((o, surM[si, k, :]))
                    for u, lst in per_unit.items():
                        if u in DROP[g] or len(lst) < 3:
                            continue
                        obs = np.array([o for o, _ in lst], float)
                        surr = np.vstack([s for _, s in lst]).astype(float)  # (K,R)
                        st = CL.all_statistics(obs, surr, upper=True)
                        st["loo_worst_signconsist"] = CL.loo_worst(obs, surr, "signconsist_p")
                        st["loo_worst_pooled"] = CL.loo_worst(obs, surr, "pooled_p")
                        st.update(dict(grouping=g, target=t, band=band,
                                       scale_idx=si, s=float(SGRID[si]), unit=u))
                        rows.append(st)
    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "all_cells.csv", index=False)
    pcols = ["pooled_p", "wilcoxon_p", "stouffer_p", "signconsist_p"]
    # per-cell BH (over units within target,band,scale,grouping) and whole-grid (per band,grouping)
    df = CL.bh_columns(df, pcols, group_cols=["grouping", "target", "band", "scale_idx"], suffix="_qcell")
    df = CL.bh_columns(df, pcols, group_cols=["grouping", "band"], suffix="_qgrid")
    for g in GROUPINGS:
        sub = df[df.grouping == g].copy()
        sub.to_csv(OUT / f"{g}_bakeoff.csv", index=False)
    _summarize(df)
    print(f"\n[C1 stage B done] -> {OUT}", flush=True)


def _summarize(df):
    key = {"trace": ["OFC", "cingulate", "PFC", "MTL", "insula"],
           "encoding": ["OFC", "cingulate", "PFC", "insula"],
           "inference": ["OFC", "cingulate", "PFC", "insula"]}
    print("\n" + "=" * 110)
    print("FOUR-STATISTIC BAKE-OFF — best-scale per pre-registered target/unit "
          "(qgrid = whole-grid BH within band x grouping)")
    print("=" * 110)
    for g in ("system", "supersystem"):
        print(f"\n########## GRANULARITY = {g} ##########")
        for t in TARGETS:
            print(f"  --- target={t} ---")
            for u in (key[t] if g == "system" else ["limbic", "PFC", "sensorimotor"]):
                sub = df[(df.grouping == g) & (df.target == t) & (df.unit == u)]
                if sub.empty:
                    continue
                b = sub.loc[sub.signconsist_p.idxmin()]
                def clears(col): return int((sub[col] < 0.05).sum())
                print(f"    {u:11s} K={int(b.n_patients)} best s{b.s:5.1f} fpos={b.frac_pos:.2f} | "
                      f"POOL p={b.pooled_p:.3f} q={b.pooled_qgrid:.3f} | WILC p={b.wilcoxon_p:.3f} "
                      f"q={b.wilcoxon_qgrid:.3f} | STOU p={b.stouffer_p:.4f} q={b.stouffer_qgrid:.3f} | "
                      f"SIGN p={b.signconsist_p:.4f} q={b.signconsist_qgrid:.3f} LOO={b.loo_worst_signconsist:.3f} "
                      f"|| grid-clears sign {clears('signconsist_qgrid')} pool {clears('pooled_qgrid')}")
    # SOZ line
    print(f"\n########## SOZ (soz vs non_soz) ##########")
    for t in TARGETS:
        sub = df[(df.grouping == "soz") & (df.target == t) & (df.unit == "soz")]
        if sub.empty:
            continue
        b = sub.loc[sub.signconsist_p.idxmin()]
        for band in BANDS:
            bb = sub[sub.band == band]
            if bb.empty:
                continue
            bx = bb.loc[bb.signconsist_p.idxmin()]
            print(f"    {t:10s} {band:10s} soz best s{bx.s:5.1f} K={int(bx.n_patients)} fpos={bx.frac_pos:.2f} "
                  f"SIGN p={bx.signconsist_p:.3f} q={bx.signconsist_qgrid:.3f} POOL p={bx.pooled_p:.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["a", "b", "ab"], default="ab")
    args = ap.parse_args()
    if args.stage in ("a", "ab"):
        stage_a()
    if args.stage in ("b", "ab"):
        stage_b()


if __name__ == "__main__":
    main()
