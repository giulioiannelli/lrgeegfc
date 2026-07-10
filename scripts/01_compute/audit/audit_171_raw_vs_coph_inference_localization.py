#!/usr/bin/env python3
"""Raw-vs-cophenetic localization of the beta INFERENCE-SPECIFIC trace (rho_sym).

Tests the one drift-free axis that could separate cophenetic from raw edges on the
inference-specific cell -- where they TIE under matched-strength (0.0098) and drift is
invalid (conditional functional). Question: does the hierarchical transform
CONCENTRATE the inference-specific per-pair trace onto an anatomical system more
sharply than raw edges, or is raw equally (co-)localized / equally placeless?

Same machinery as audit_160/audit_110, run for TWO representations, everything else
identical (concordance_partial(f,p|e), rho_sym arm-symmetrized, demeaned per-system
unit_means, R=1000 matched-strength cached ensemble seed 20260511):
  coph : D[ph] = canonical cophenetic of W[ph]   (surrogate: cophenetic_from_eigs)
  raw  : D[ph] = raw upper-tri edges of W[ph]     (surrogate: raw_edges(adjacency_from_eigs))
Only phi(W) changes across the two runs; the null draws (eigs) are the SAME.

5-point preamble
1. Claim: the cophenetic inference-specific concordance concentrates on an anatomical
   system (cingulate/OFC) ABOVE matched-strength MORE than the raw-edge concordance
   does; raw is comparatively placeless.
2. Null: matched-strength (R=1000, strength-preserving, seed 20260511) reproduces any
   per-system concentration that is a function of node strength.
3. Strongest alternative: co-localization -- raw is ALSO concentrated on the same
   system, so there is NO separation (localization inherited, not a hierarchy value-add).
   coph = deterministic f(W) => raw carries >= coph's info => co-localization is the
   DEFAULT expectation.
4. Null's reach: MS controls strength; it does NOT make coph localize where raw doesn't
   -- only the DATA can. CANNOT reject: this is a CONCENTRATION/sharpness contrast, not
   detection (they tie on detection by construction). CEILING: cophenetic's own
   inference_pe localization is a DEMOTED HINT (cingulate BH-borderline, downgraded by
   the length control). Likely outcomes: both weak/placeless (honest null) or raw
   co-localizes; a clean 'coph sharp, raw flat' is possible but not the prior.
5. Falsify: if raw's leading inference system BH-clears at parity-or-better with coph,
   or neither concentrates, there is NO localization separator -- report plainly.
"""
from __future__ import annotations
import sys, time, argparse
from multiprocessing import Pool
import numpy as np
import pandas as pd
from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.node_localization import epi_keep_mask
from lrg_eegfc.utils.metrics.graph_descriptors import raw_edges
from lrg_eegfc.utils.surrogate import cophenetic_condensed_from_eigs, surrogate_cache_path
from lrg_eegfc.utils.surrogate.matched_strength import adjacency_from_laplacian_eigs
from audit_63_split_baseline_surrogate import ensure_half_fcs, load_phase_fc  # type: ignore
from audit_83_localization_matched_strength import (  # type: ignore
    unit_means_from_s, _canon_cophenet, DROP, PARACORE,
)
from audit_110_inference_mark_localization import _targets_from_D  # type: ignore

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
BAND = "beta"
PHASES5 = ("rest_pre_A", "rest_pre_B", "task_learn", "task_test", "rest_post")
SWAP, SEED, R = 20, 20260511, 1000


def repr_condensed(name, W):
    return _canon_cophenet(W) if name == "coph" else raw_edges(W)


def _swap(D):
    Ds = dict(D)
    Ds["rest_pre_A"], Ds["rest_pre_B"] = D["rest_pre_B"], D["rest_pre_A"]
    return Ds


def sym_infpe(D):
    """rho_sym per-pair inference_pe concordance (arm1+arm2)/2 (== audit_160)."""
    s1 = _targets_from_D(D)["inference_pe"]
    s2 = _targets_from_D(_swap(D))["inference_pe"]
    return 0.5 * (s1 + s2)


def surrogate_repr5(pat, name):
    out = {}
    for ph in PHASES5:
        path = surrogate_cache_path(pat, BAND, ph, R, SWAP, SEED, "imcoh_abs")
        if not path.exists():
            return None, ph
        with np.load(path) as d:
            evals, evecs = d["eigvals"], d["eigvecs"]
        vv = []
        for r in range(evals.shape[0]):
            ev = evals[r]
            if not np.all(np.isfinite(ev)):
                vv.append(None); continue
            if name == "coph":
                vv.append(cophenetic_condensed_from_eigs(ev, evecs[r]))
            else:
                vv.append(raw_edges(adjacency_from_laplacian_eigs(ev, evecs[r])))
        out[ph] = vv
    return out, None


def patient_systems(args):
    """{system: (M_obs_pat, surr_array)} for inference_pe rho_sym, this patient+repr."""
    pat, name, epi_excl = args
    surco, miss = surrogate_repr5(pat, name)
    if surco is None:
        return pat, name, epi_excl, None
    try:
        D = {ph: repr_condensed(name, load_phase_fc(pat, ph, BAND)) for ph in PHASES5}
    except FileNotFoundError:
        return pat, name, epi_excl, None
    s_obs = sym_infpe(D)
    N = int((1 + np.sqrt(1 + 8 * D["task_test"].size)) / 2)
    iu = np.triu_indices(N, k=1)
    iu_i, iu_j = iu[0].astype(int), iu[1].astype(int)
    rdf = load_channel_regions(pat)
    rdf["paracore"] = rdf["system"].map(
        lambda s: "paralimbic_core" if s in PARACORE
        else ("non_anatomical" if s == "non_anatomical" else "other_np"))
    keep = epi_keep_mask(rdf, pat) if epi_excl else None
    uv = rdf["system"].to_numpy().astype(str)
    om = unit_means_from_s(s_obs, iu_i, iu_j, uv, keep, DROP["system"])
    per = {u: [] for u in om}
    for r in range(R):
        Dr = {ph: surco[ph][r] for ph in PHASES5}
        if any(Dr[ph] is None for ph in PHASES5):
            continue
        sr = sym_infpe(Dr)
        sm = unit_means_from_s(sr, iu_i, iu_j, uv, keep, DROP["system"])
        for u in per:
            if u in sm:
                per[u].append(sm[u])
    return pat, name, epi_excl, {u: (om[u], np.array(per[u])) for u in om}


def _ms_p(obs_list, mat):
    M_obs = float(np.median(obs_list))
    M_surr = np.median(mat, axis=0)
    return (1 + int(np.sum(M_surr >= M_obs))) / (M_surr.size + 1)


def cohort_table(cells, name, epi_excl):
    """cells: list of (system-dict). Cohort M_obs + upper-tail MS_p + BH + LOO per system."""
    by_u = {}
    for c in cells:
        if c is None:
            continue
        for u, (mo, sa) in c.items():
            by_u.setdefault(u, []).append((mo, sa))
    rows = []
    for u, pairs in by_u.items():
        paired = [(mo, a) for mo, a in pairs if a is not None and len(a) > 0]
        if not paired:
            continue
        obsv = [mo for mo, _ in paired]
        Lmin = min(len(a) for _, a in paired)
        Rmat = np.vstack([a[:Lmin] for _, a in paired])
        p_full = _ms_p(obsv, Rmat)
        K = len(obsv)
        loo = [_ms_p([obsv[j] for j in range(K) if j != i],
                     np.vstack([Rmat[j] for j in range(K) if j != i]))
               for i in range(K)] if K > 1 else [p_full]
        rows.append({"repr": name, "epi": "exclude" if epi_excl else "include",
                     "system": u, "K": K, "M_obs": float(np.median(obsv)),
                     "surr_med": float(np.median(np.median(Rmat, axis=0))),
                     "MS_p": p_full, "MS_p_loo_max": max(loo), "MS_p_loo_min": min(loo)})
    df = pd.DataFrame(rows)
    if not df.empty:
        df["BH_q"] = bh_fdr(df["MS_p"].values)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="first N patients (timing)")
    args = ap.parse_args()
    cohort = COHORT[:args.limit] if args.limit else COHORT

    for pat in cohort:
        ensure_half_fcs(pat, [BAND])

    jobs = [(pat, name, epi) for name in ("coph", "raw")
            for epi in (False, True) for pat in cohort]
    t0 = time.time()
    results = {}
    ncpu = min(10, len(jobs))
    with Pool(ncpu) as pool:
        for i, (pat, name, epi, cell) in enumerate(pool.imap_unordered(patient_systems, jobs), 1):
            results.setdefault((name, epi), []).append(cell)
            el = time.time() - t0
            print(f"[{i}/{len(jobs)}] {pat} {name} epi={epi}  {el:.0f}s "
                  f"ETA {el/i*(len(jobs)-i):.0f}s", flush=True)

    alldf = []
    for (name, epi), cells in results.items():
        alldf.append(cohort_table(cells, name, epi))
    out = pd.concat(alldf, ignore_index=True)
    out.to_csv(ROOT / "data" / "audit" / "raw_vs_coph_inference_localization.csv", index=False)

    print("\n" + "=" * 78)
    print("BETA inference_pe (rho_sym) localization: COPH vs RAW  (R=1000 matched-strength)")
    print("per-system demeaned concentration; MS_p = upper tail; BH over systems")
    print("=" * 78)
    for epi in (False, True):
        etag = "exclude" if epi else "include"
        print(f"\n---- epi-{etag} ----")
        for name in ("coph", "raw"):
            sub = out[(out["repr"] == name) & (out.epi == etag)].sort_values("MS_p")
            print(f"  [{name}]  leading systems (by MS_p):")
            for _, r in sub.head(5).iterrows():
                star = "*BH" if r.BH_q < 0.05 else ("." if r.MS_p < 0.05 else " ")
                print(f"     {r.system:16s} M_obs={r.M_obs:+.3e} MS_p={r.MS_p:.3f} "
                      f"[LOO {r.MS_p_loo_min:.3f}..{r.MS_p_loo_max:.3f}] "
                      f"BH_q={r.BH_q:.3f} {star}")
    print(f"\n[done] {time.time()-t0:.0f}s -> data/audit/raw_vs_coph_inference_localization.csv")


if __name__ == "__main__":
    main()
