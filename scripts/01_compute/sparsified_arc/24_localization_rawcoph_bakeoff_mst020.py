#!/usr/bin/env python3
"""C4 — raw-vs-cophenetic localization value-add, four-statistic bake-off.

The load-bearing "localization is the cophenetic contribution" claim (audit_171:
coph inference->cingulate clears, raw edges are placeless) was pooled-median-based on
the DENSE degenerate graph. This ports it onto the mst@0.20 regime and, like C1
(script 23), reports ALL FOUR cohort statistics — so the value-add is judged by the
SAME statistic on BOTH representations, not by a statistic that flatters one.

Two representations of each phase's per-pair geometry (both dense per-pair vectors,
comparable):
  coph : the mst@0.20 cophenetic condensed distance at scale s   (from script 23's dump)
  raw  : the DENSE |ImCoh|_abs upper-triangle edges              (audit_171 'raw', no transform)
Identical per-pair concordance (trace/encoding/inference, rho_sym-symmetrised) + the
identical endpoint-incidence demeaned per-unit means + the identical matched-strength
null (R=200). The raw arm needs NO eigendecomposition (cheap). Value-add = a system
where coph clears (sign-consistency whole-grid q<0.05, LOO-robust) AND raw does not.

5-point preamble
1. Claim: the cophenetic transform CONCENTRATES the per-pair concordance on an
   anatomical system more than the raw edges do (raw comparatively placeless).
2. Null: matched-strength R=200 (identical to script 23), applied THROUGH each
   representation (shuffle dense FC -> raw triu for the raw arm; -> backbone -> coph for
   the coph arm).
3. Strongest alternative: co-localization — raw is EQUALLY concentrated (coph =
   deterministic f(backbone(W)) carries <= the info of W's edges), so the DEFAULT prior
   is that raw sees whatever coph sees. A clean 'coph sharp, raw flat' is the burden.
4. Cannot: this is a CONCENTRATION contrast, not detection (they tie on detection by
   construction). Whole-grid BH + LOO gate both arms identically.
5. Falsify: if raw's leading system clears at parity-or-better with coph under the same
   statistic, there is NO localization value-add — the hierarchy adds no anatomy.

Reuses script 23's coph arrays (data/.../localization_bakeoff_mst020/arrays). Run C1
(script 23) first. Output: data/sparsified_arc/localization_rawcoph_mst020/
Run: [SA_LIMIT=n] python 24_localization_rawcoph_bakeoff_mst020.py [--stage a|b|ab]
"""
from __future__ import annotations
import os, sys, time, argparse
import numpy as np, pandas as pd
from multiprocessing import Pool
from importlib import import_module

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "sparsified_arc"))
from audit_150_rho_sym_gate import COHORT, load_phase                 # noqa: E402
from lrg_eegfc.utils.metrics.surrogate import matched_strength_shuffle  # noqa: E402
from lrg_eegfc.utils.metrics.node_localization import epi_keep_mask     # noqa: E402
from lrg_eegfc.utils.io.regions import load_channel_regions            # noqa: E402
from lrg_eegfc.utils.metrics import cohort_localization as CL          # noqa: E402
_m17 = import_module("17_localization_arc_mst020")
sym_trace, sym_encoding, sym_inference = _m17.sym_trace, _m17.sym_encoding, _m17.sym_inference
unit_means = _m17.unit_means

PHASES = ("A", "B", "task_learn", "task_test", "rest_post")
BANDS = ["delta", "alpha", "beta", "low_gamma"]
TARGETS = ("trace", "encoding", "inference")
GROUPINGS = ("system", "supersystem", "soz")
R = 200
SWAP_FACTOR = 20
W_MAX = 1.0
BASE_SEED = 20260713          # SAME seed base as script 23 -> paired null draws
DROP = {"system": {"non_anatomical"},
        "supersystem": {"non_anatomical", "other"}, "soz": set()}
LIMIT = int(os.environ.get("SA_LIMIT", "0"))
COPH_DIR = ROOT / "data" / "sparsified_arc" / "localization_bakeoff_mst020"
OUT = ROOT / "data" / "sparsified_arc" / "localization_rawcoph_mst020"
ARR = OUT / "arrays"


def raw_triu(W):
    N = W.shape[0]; iu = np.triu_indices(N, 1)
    return np.asarray(W)[iu]


def concordances_raw(V):
    cA, cB, cL, cT, cP = (V["A"], V["B"], V["task_learn"], V["task_test"], V["rest_post"])
    return {"trace": sym_trace(cA, cB, cT, cP),
            "encoding": sym_encoding(cA, cB, cL, cP),
            "inference": sym_inference(cA, cB, cL, cT, cP)}


def label_vectors(rdf, pat):
    return {"system": rdf["system"].to_numpy().astype(str),
            "supersystem": rdf["supersystem"].to_numpy().astype(str),
            "soz": np.where(epi_keep_mask(rdf, pat), "non_soz", "soz")}


def per_cell(job):
    """RAW-arm per (patient, band): obs + R surrogate per-unit-means (scale-free)."""
    idx, pat, band = job
    outp = ARR / f"{pat}_{band}.npz"
    if outp.exists():
        return (pat, band, -1)
    try:
        Ws = {ph: load_phase(pat, ph, band) for ph in PHASES}
    except Exception:
        return (pat, band, 0)
    N = Ws["A"].shape[0]
    rdf = load_channel_regions(pat)
    if len(rdf) != N:
        return (pat, band, 0)
    labs = label_vectors(rdf, pat)
    iu = np.triu_indices(N, 1); ii, jj = iu[0].astype(int), iu[1].astype(int)

    def um(Vdict):
        conc = concordances_raw({ph: raw_triu(Vdict[ph]) for ph in PHASES})
        return {t: {g: unit_means(conc[t], ii, jj, labs[g], drop=DROP[g])
                    for g in GROUPINGS} for t in TARGETS}

    obs = um(Ws)
    units = {g: sorted({u for t in TARGETS for u in obs[t][g]}) for g in GROUPINGS}
    uidx = {g: {u: k for k, u in enumerate(units[g])} for g in GROUPINGS}
    obs_T = {(t, g): np.full(len(units[g]), np.nan) for t in TARGETS for g in GROUPINGS}
    sur_T = {(t, g): np.full((len(units[g]), R), np.nan) for t in TARGETS for g in GROUPINGS}
    for t in TARGETS:
        for g in GROUPINGS:
            for u, v in obs[t][g].items():
                obs_T[(t, g)][uidx[g][u]] = v

    rng = np.random.default_rng(BASE_SEED + idx)
    n_swaps = SWAP_FACTOR * N * (N - 1) // 2
    for r in range(R):
        sh = {ph: matched_strength_shuffle(Ws[ph], n_swaps, rng, W_MAX) for ph in PHASES}
        try:
            sm = um(sh)
        except Exception:
            continue
        for t in TARGETS:
            for g in GROUPINGS:
                for u, v in sm[t][g].items():
                    sur_T[(t, g)][uidx[g][u], r] = v

    save = {}
    for g in GROUPINGS:
        save[f"{g}__units"] = np.array(units[g], dtype=object)
    for t in TARGETS:
        for g in GROUPINGS:
            save[f"{t}__{g}__obs"] = obs_T[(t, g)].astype(np.float32)
            save[f"{t}__{g}__surr"] = sur_T[(t, g)].astype(np.float32)
    ARR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(outp, **save)
    return (pat, band, N)


def stage_a():
    ARR.mkdir(parents=True, exist_ok=True)
    _ = matched_strength_shuffle(np.zeros((5, 5)), 4, np.random.default_rng(0))
    jobs = [(i, pat, band) for i, (pat, band) in
            enumerate((p, b) for b in BANDS for p in COHORT)]
    if LIMIT:
        jobs = jobs[:LIMIT]
    print(f"[C4 raw-arm stage A] R={R} {len(jobs)} (pat,band) cells (no eig) -> {ARR}", flush=True)
    t0 = time.time()
    nproc = min(12, max(1, (os.cpu_count() or 2) - 2))
    with Pool(nproc) as pool:
        for k, (pat, band, N) in enumerate(pool.imap_unordered(per_cell, jobs), 1):
            el = time.time() - t0
            print(f"  [{k}/{len(jobs)}] {pat} {band} {'cached' if N==-1 else N}  "
                  f"elapsed {el:5.0f}s ETA {el/k*(len(jobs)-k):5.0f}s", flush=True)
    print(f"[C4 raw-arm stage A done] {time.time()-t0:.0f}s", flush=True)


def _raw_cohort():
    """Cohort-combine the RAW arm (scale-free) -> all 4 statistics per (target,band,g,unit)."""
    rows = []
    npzs = {(pat, band): np.load(ARR / f"{pat}_{band}.npz", allow_pickle=True)
            for band in BANDS for pat in COHORT if (ARR / f"{pat}_{band}.npz").exists()}
    for g in GROUPINGS:
        for band in BANDS:
            for t in TARGETS:
                per_unit = {}
                for pat in COHORT:
                    if (pat, band) not in npzs:
                        continue
                    d = npzs[(pat, band)]
                    units = list(d[f"{g}__units"]); obsV = d[f"{t}__{g}__obs"]; surV = d[f"{t}__{g}__surr"]
                    for k, u in enumerate(units):
                        o = float(obsV[k])
                        if np.isfinite(o):
                            per_unit.setdefault(u, []).append((o, surV[k, :]))
                for u, lst in per_unit.items():
                    if u in DROP[g] or len(lst) < 3:
                        continue
                    obs = np.array([o for o, _ in lst], float)
                    surr = np.vstack([s for _, s in lst]).astype(float)
                    st = CL.all_statistics(obs, surr, upper=True)
                    st["loo_worst_signconsist"] = CL.loo_worst(obs, surr, "signconsist_p")
                    st.update(dict(grouping=g, target=t, band=band, unit=u, repr="raw"))
                    rows.append(st)
    return pd.DataFrame(rows)


def stage_b():
    raw = _raw_cohort()
    OUT.mkdir(parents=True, exist_ok=True)
    raw = CL.bh_columns(raw, ["pooled_p", "wilcoxon_p", "stouffer_p", "signconsist_p"],
                        group_cols=["grouping", "band"], suffix="_qgrid")
    raw.to_csv(OUT / "raw_cohort.csv", index=False)

    # coph best-scale per (grouping,target,band,unit) from C1
    coph_all = COPH_DIR / "all_cells.csv"
    if not coph_all.exists():
        print(f"[C4] coph {coph_all} missing — run C1 (script 23) first; raw arm saved.", flush=True)
        return
    coph = pd.read_csv(coph_all)
    # recompute coph whole-grid q if not present
    if "signconsist_qgrid" not in coph.columns:
        coph = CL.bh_columns(coph, ["pooled_p", "signconsist_p"],
                             group_cols=["grouping", "band"], suffix="_qgrid")
    coph_best = (coph.sort_values("signconsist_p")
                 .groupby(["grouping", "target", "band", "unit"], as_index=False).first())
    merged = coph_best.merge(raw, on=["grouping", "target", "band", "unit"],
                             suffixes=("_coph", "_raw"))
    merged.to_csv(OUT / "rawcoph_merged.csv", index=False)
    _summarize(merged)


def _summarize(m):
    print("\n" + "=" * 108)
    print("RAW-vs-COPH LOCALIZATION VALUE-ADD (sign-consistency, whole-grid BH per band x grouping)")
    print("coph = best scale | raw = scale-free | VALUE-ADD = coph q<.05 & LOO<.05 & raw q>=.05")
    print("=" * 108)
    for g in ("system", "supersystem"):
        print(f"\n########## {g} ##########")
        sub = m[m.grouping == g]
        va = sub[(sub.signconsist_qgrid_coph < 0.05) & (sub.loo_worst_signconsist_coph < 0.05)
                 & (sub.signconsist_qgrid_raw >= 0.05)]
        both = sub[(sub.signconsist_qgrid_coph < 0.05) & (sub.signconsist_qgrid_raw < 0.05)]
        rawonly = sub[(sub.signconsist_qgrid_coph >= 0.05) & (sub.signconsist_qgrid_raw < 0.05)]
        print(f"  coph-value-add cells: {len(va)} | co-localized (both): {len(both)} | raw-only: {len(rawonly)}")
        for _, r in va.sort_values("signconsist_p_coph").head(12).iterrows():
            print(f"    +VALUE {r.target:9s} {r.unit:11s} {r.band:9s} s{r.s:5.1f}: "
                  f"coph q={r.signconsist_qgrid_coph:.3f} (LOO {r.loo_worst_signconsist_coph:.3f}) "
                  f"vs raw q={r.signconsist_qgrid_raw:.3f}")
        for _, r in both.sort_values("signconsist_p_coph").head(8).iterrows():
            print(f"    =BOTH  {r.target:9s} {r.unit:11s} {r.band:9s}: coph q={r.signconsist_qgrid_coph:.3f} "
                  f"raw q={r.signconsist_qgrid_raw:.3f}")


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
