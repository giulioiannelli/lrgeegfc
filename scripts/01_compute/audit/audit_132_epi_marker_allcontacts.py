#!/usr/bin/env python3
"""Audit 132 — ALL-CONTACTS connectivity marker for epileptic nodes, judged ONLY against
the matched-strength null (no shaft holdout, no proximity control).

Why this design (the correction, 2026-06-23). |ImCoh| is zero-lag-immune (Nolte 2004), so
physical proximity CANNOT enter the connectivity. Holding shafts out (audit_101 LOSO) was
answering a *clinical* question ("does it find SOZ a distance ruler can't") by handicapping
the connectivity claim. For the CONNECTIVITY claim — "epileptic nodes are read out by the
|ImCoh| propagator beyond node strength" — the only mandatory null is matched-strength, on
ALL contacts. This script recomputes the marker that way: leave-one-SOZ-contact-out, seed
with every other SOZ (same-shaft included), score the held-out SOZ vs ALL healthy contacts,
strength-residualised, and compare to a strength-matched fake-SOZ label null.

Critical preamble
=================
(1) Claim: the |ImCoh| propagator ranks SOZ above healthy contacts beyond node strength,
    cohort-wide, with NO proximity control needed (|ImCoh| removes volume conduction).
(2) Null: matched-strength label shuffle — draw a fake-SOZ set of the SAME size and SAME
    strength profile as the true SOZ, run the identical leave-one-out marker, AUC. If the
    true SOZ is not more self-recoverable than strength-matched random sets, there is no
    connectivity-specific signal.
(3) Strongest alternative: SOZ are just hubs -> killed two ways here (strength-RESIDUAL AUC
    AND a strength-MATCHED null). A pure spatial-clustering tautology is NOT a concern: it
    would inflate a distance baseline, not the |ImCoh| graph, and is not what is tested.
(4) Reach: same operators/folds for every marker; markers not tuned per patient; the headline
    (heat tau5) must clear the matched-strength null per patient and cohort, not just beat 0.5.
(5) Falsification: if heat tau5 all-contacts strength-residual AUC does not exceed the
    matched-strength null, the connectivity does not read out SOZ. Reported flat.

Outputs (data/audit/epi_marker_allcontacts/):
  allcontacts_per_patient.csv ; allcontacts_cohort.csv ; allcontacts_nulls.csv ;
  allcontacts_stratnull_per_patient.csv ; README.md
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.utils.io.patient import build_epi_masks

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI  # type: ignore
from audit_101_epi_marker_library import (  # type: ignore
    build_operators, _resid, _conc, MIN_PAIRS,
)

OUT = ROOT / "data" / "audit" / "epi_marker_allcontacts"
HEAD = "heat_t5"
N_PERM = 300
N_BINS = 4


def loco_auc(marker_of_seeds, epi_idx, healthy, strength):
    """Leave-one-SOZ-contact-out, ALL-CONTACTS AUC (raw + strength-residual).

    seeds = every OTHER SOZ (same shaft included); held-out SOZ scored vs ALL healthy.
    """
    if epi_idx.size < 2 or healthy.size < 3:
        return np.nan, np.nan, 0
    raw = [0.0, 0]; res = [0.0, 0]
    for c in epi_idx:
        seeds = epi_idx[epi_idx != c]
        if seeds.size < 1:
            continue
        m = marker_of_seeds(seeds)
        mr = _resid(m, strength)
        cc, tt = _conc(m[[c]], m[healthy]); raw[0] += cc; raw[1] += tt
        cc, tt = _conc(mr[[c]], mr[healthy]); res[0] += cc; res[1] += tt
    if raw[1] < MIN_PAIRS:
        return np.nan, np.nan, raw[1]
    return raw[0] / raw[1], (res[0] / res[1] if res[1] >= MIN_PAIRS else np.nan), raw[1]


def _strength_matched_sets(strength, epi_idx, n_perm, rng, n_bins=N_BINS):
    """Random node sets of the same size + same strength-quantile profile as the true SOZ."""
    N = strength.size
    order = np.argsort(strength)
    ranks = np.empty(N, int); ranks[order] = np.arange(N)
    bins = np.minimum(ranks * n_bins // N, n_bins - 1)
    counts = [int((bins[epi_idx] == b).sum()) for b in range(n_bins)]
    pools = [np.where(bins == b)[0] for b in range(n_bins)]
    out = []
    for _ in range(n_perm):
        pick = []
        ok = True
        for b in range(n_bins):
            if counts[b] == 0:
                continue
            if pools[b].size < counts[b]:
                ok = False; break
            pick.append(rng.choice(pools[b], size=counts[b], replace=False))
        if ok and pick:
            out.append(np.concatenate(pick))
    return out


def per_patient_band(pat, band, rng):
    try:
        W = load_phase_fc(pat, "rest_post", band)
    except Exception:
        return None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None
    epi = np.asarray(pm.epi_mask, bool)
    epi_idx = np.where(epi)[0]
    healthy = np.where(~epi)[0]
    if epi_idx.size < max(MIN_EPI, 4) or healthy.size < 5:
        return None
    ops, _intrinsic, strength = build_operators(W)

    rows = []
    for name, O in ops.items():
        def marker(seeds, O=O):
            return O[:, seeds].mean(axis=1)
        a_raw, a_res, npairs = loco_auc(marker, epi_idx, healthy, strength)
        rows.append({"patient": pat, "band": band, "marker": name, "kind": "seedbased",
                     "auc_raw": a_raw, "auc_resid": a_res, "n_pairs": npairs})

    def str_marker(seeds, s=strength):
        return s
    a_raw, a_res, npairs = loco_auc(str_marker, epi_idx, healthy, strength)
    rows.append({"patient": pat, "band": band, "marker": "strength_baseline",
                 "kind": "baseline", "auc_raw": a_raw, "auc_resid": np.nan,
                 "n_pairs": npairs})

    # matched-strength null for the headline marker (heat tau5)
    O = ops[HEAD]
    obs = next(r for r in rows if r["marker"] == HEAD)["auc_resid"]
    null_aucs = []
    for fake in _strength_matched_sets(strength, epi_idx, N_PERM, rng):
        fake_heal = np.setdiff1d(np.arange(N), fake, assume_unique=False)

        def fmarker(seeds, O=O):
            return O[:, seeds].mean(axis=1)
        _, a_res_n, _ = loco_auc(fmarker, fake, fake_heal, strength)
        if np.isfinite(a_res_n):
            null_aucs.append(a_res_n)
    null_aucs = np.array(null_aucs, float)
    null_info = {"band": band, "patient": pat, "obs_resid": obs,
                 "null_p50": float(np.median(null_aucs)) if null_aucs.size else np.nan,
                 "null_p95": float(np.percentile(null_aucs, 95)) if null_aucs.size else np.nan,
                 "p_strat": (float(np.mean(null_aucs >= obs)) if (null_aucs.size and np.isfinite(obs)) else np.nan),
                 "beats_p95": (bool(obs >= np.percentile(null_aucs, 95)) if (null_aucs.size and np.isfinite(obs)) else False)}
    return pd.DataFrame(rows), null_info, null_aucs


def cohort(df):
    out = []
    for band in [b for b in es.ALL_BANDS if b in set(df.band)]:
        for marker in sorted(df.marker.unique()):
            d = df[(df.band == band) & (df.marker == marker)]
            for metric in ("auc_raw", "auc_resid"):
                v = d[metric].to_numpy(float); v = v[np.isfinite(v)]
                if v.size < 4:
                    continue
                p = np.nan
                try:
                    p = float(wilcoxon(v - 0.5, alternative="greater")[1])
                except ValueError:
                    pass
                out.append({"band": band, "marker": marker, "metric": metric,
                            "n_patients": int(v.size), "median_auc": float(np.median(v)),
                            "mean_auc": float(np.mean(v)),
                            "n_above_half": int((v > 0.5).sum()), "p_sign": p})
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=["delta", "beta", "low_gamma", "alpha"])
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rng = np.random.default_rng(20260623)
    frames, null_rows = [], []
    null_arrays = {}  # (band) -> list of per-patient null arrays
    for band in args.bands:
        null_arrays[band] = []
        for pat in args.patients:
            r = per_patient_band(pat, band, rng)
            if r is None:
                continue
            rows, ninfo, narr = r
            frames.append(rows)
            null_rows.append(ninfo)
            if narr.size:
                null_arrays[band].append(narr)
        print(f"[audit_132] {band}: done")
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(OUT / "allcontacts_per_patient.csv", index=False)
    coh = cohort(df)
    coh.to_csv(OUT / "allcontacts_cohort.csv", index=False)
    snp = pd.DataFrame(null_rows)
    snp.to_csv(OUT / "allcontacts_stratnull_per_patient.csv", index=False)

    # cohort-level matched-strength null for heat tau5 (median over patients per perm)
    nulls = []
    for band in args.bands:
        arrs = null_arrays.get(band, [])
        if not arrs or coh.empty:
            continue
        sub_obs = coh[(coh.band == band) & (coh.marker == HEAD)
                      & (coh.metric == "auc_resid")]
        if not len(sub_obs):
            continue
        m = min(a.size for a in arrs)
        M = np.vstack([a[:m] for a in arrs])           # patients x perms
        null_cohort_median = np.median(M, axis=0)       # per-perm cohort median
        obs_med = float(sub_obs.median_auc.iloc[0])
        nulls.append({"band": band, "marker": HEAD, "real_median_auc": obs_med,
                      "null_median_mean": float(np.mean(null_cohort_median)),
                      "null_median_p95": float(np.percentile(null_cohort_median, 95)),
                      "p_empirical": float(np.mean(null_cohort_median >= obs_med)),
                      "n_perm": int(null_cohort_median.size)})
    pd.DataFrame(nulls).to_csv(OUT / "allcontacts_nulls.csv", index=False)

    rt = time.time() - t0
    print("\n[audit_132] ALL-CONTACTS heat tau5 vs matched-strength null:")
    for n in nulls:
        nb = snp[(snp.band == n["band"])]
        nbeat = int(nb.beats_p95.sum())
        print(f"  {n['band']:10s} median AUC_resid={n['real_median_auc']:.3f}  "
              f"null p95={n['null_median_p95']:.3f}  p={n['p_empirical']:.3f}  "
              f"per-pt beat-own-p95={nbeat}/{len(nb)}")
    print("\n[audit_132] heat tau5 RAW (incl strength) all-contacts median AUC:")
    for band in args.bands:
        r = coh[(coh.band == band) & (coh.marker == HEAD) & (coh.metric == "auc_raw")] if not coh.empty else coh
        if not coh.empty and len(r):
            print(f"  {band:10s} median AUC_raw={float(r.median_auc.iloc[0]):.3f} "
                  f"({int(r.n_above_half.iloc[0])}/{int(r.n_patients.iloc[0])}>0.5)")
    print(f"[audit_132] done in {rt:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
