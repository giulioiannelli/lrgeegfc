#!/usr/bin/env python3
"""Audit 115 — recall GROWTH curve: as you include more of the top-ranked contacts (rank
depth N grows), what fraction of the remaining distant SOZ is recovered, recall(N)? And how
deep (as a fraction of the contact pool) must you read to reach 50/80/90/100% recovery —
versus a random list (where depth-to-X% = X%) and a node-strength list? The cumulative /
CMC view of the off-shaft delta heat-kernel marker; the honest baseline that any richer or
iterative use of the propagator must beat.

Marker = delta heat kernel e^{-tau L} at slow tau (heat_t5), strength-residualised
(audit_101/102 winner). Off-shaft (proximity removed): never score a contact sharing an
electrode with any seed (exactly audit_101/113/114).

Critical preamble
=================
(1) Claim: recall(N) climbs FASTER than a random list -> the marker reaches a given recovery
    fraction at a SMALLER read depth than chance (depth-to-X% << X% of the pool), with the
    early part of the curve (top contacts) strongly enriched.
(2) Null: recall(N) = N/|pool| (random ordering); equivalently depth-to-X% recall = X% of the
    pool. NOTE: recall(N=|pool|) = 1 for ANY ordering -> "full recall at full depth" is
    trivial and is NOT evidence; only the climb rate is.
(3) Strongest alternative: hubness (killed by the strength residual AND an explicit
    strength-RANKED growth curve at the same depths); proximity (killed off-shaft).
(4) Reach: k seeds drawn at random per patient (NDRAW draws), off pool rebuilt per draw;
    recall(N) at depths N=1..50 + depth-to-{50,80,90,100}% recovery per draw (absolute and
    as fraction of pool); per-patient + cohort (community responders); hub patients
    (Pat_10/15) reported, not hidden.
(5) Falsification: if recall(N) tracks the random diagonal (depth-to-X% ~ X%), the marker
    does not concentrate SOZ and the growth is reported as no-better-than-chance. No
    pre-registered acceptance gate (post-hoc evaluation with PI).

Outputs (data/audit/epi_marker_recall_growth/)
    recall_curve_long.csv ; recall_depth_summary.csv ; README.md
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.utils.io.patient import build_epi_masks

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI  # type: ignore
from audit_101_epi_marker_library import build_operators, _resid, _conc  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_marker_recall_growth"
ML = ROOT / "data" / "audit" / "epi_marker_library"
BAND = "delta"
MARKER = "heat_t5"
KSEEDS = [3, 5]
DEPTHS = list(range(1, 51))
XLEVELS = [0.5, 0.8, 0.9, 1.0]
MIN_OFF_CTRL = 3
RESP_THR = 0.568
SEED = 20260619


def _depth_to(cum, x, m, npool):
    """smallest read depth N (1-based) at which cumulative hits reach ceil(x*m)."""
    need = int(np.ceil(x * m))
    idx = int(np.searchsorted(cum, need, side="left"))
    return (idx + 1) if idx < len(cum) else npool


def per_patient(pat, ndraw, rng):
    try:
        W = load_phase_fc(pat, "rest_post", BAND)
    except Exception:
        return None, None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None, None
    epi = np.asarray(pm.epi_mask, bool)
    probes = np.asarray(pm.probes, object)
    epi_idx = np.where(epi)[0]
    if epi_idx.size < max(MIN_EPI, 4) or len(set(probes[epi_idx])) < 2:
        return None, None
    ops, _intr, strength = build_operators(W)
    O = ops[MARKER]

    curve_rows, sum_rows = [], []
    for k in KSEEDS:
        if epi_idx.size <= k:
            continue
        cp = {d: [] for d in DEPTHS}
        cs = {d: [] for d in DEPTHS}
        cr = {d: [] for d in DEPTHS}
        d2x_abs = {x: [] for x in XLEVELS}
        d2x_frac = {x: [] for x in XLEVELS}
        pools, ms, aucs = [], [], []
        for _ in range(ndraw):
            seeds = rng.choice(epi_idx, size=k, replace=False)
            seed_pr = set(probes[seeds])
            seed_set = set(seeds.tolist())
            off = np.array([j for j in range(N) if probes[j] not in seed_pr])
            tgt = np.array([j for j in off if epi[j] and j not in seed_set])
            ctl = np.array([j for j in off if not epi[j]])
            if tgt.size < 1 or ctl.size < MIN_OFF_CTRL:
                continue
            aff = _resid(O[:, seeds].mean(axis=1), strength)
            c, t = _conc(aff[tgt], aff[ctl])
            aucs.append(c / t if t else np.nan)
            m = tgt.size
            npool = off.size
            pools.append(npool)
            ms.append(m)
            order = off[np.argsort(-aff[off])]
            order_s = off[np.argsort(-strength[off])]
            tgt_set = set(tgt.tolist())
            cum = np.cumsum([1 if j in tgt_set else 0 for j in order])
            cum_s = np.cumsum([1 if j in tgt_set else 0 for j in order_s])
            for d in DEPTHS:
                cp[d].append((cum[d - 1] if d <= len(cum) else cum[-1]) / m)
                cs[d].append((cum_s[d - 1] if d <= len(cum_s) else cum_s[-1]) / m)
                cr[d].append(min(d, npool) / npool)
            for x in XLEVELS:
                Nx = _depth_to(cum, x, m, npool)
                d2x_abs[x].append(Nx)
                d2x_frac[x].append(Nx / npool)
        if not aucs:
            continue
        for d in DEPTHS:
            curve_rows.append(dict(patient=pat, k=k, depth=d,
                                   recall_prop=float(np.nanmean(cp[d])),
                                   recall_str=float(np.nanmean(cs[d])),
                                   recall_rand=float(np.nanmean(cr[d]))))
        s = dict(patient=pat, k=k, pool=float(np.mean(pools)),
                 n_off_soz=float(np.mean(ms)), auc=float(np.nanmean(aucs)))
        for x in XLEVELS:
            s[f"d{int(x*100)}_abs"] = float(np.mean(d2x_abs[x]))
            s[f"d{int(x*100)}_frac"] = float(np.mean(d2x_frac[x]))
        sum_rows.append(s)
    return (pd.DataFrame(curve_rows) if curve_rows else None,
            pd.DataFrame(sum_rows) if sum_rows else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ndraw", type=int, default=300)
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    t0 = time.time()

    try:
        ml = pd.read_csv(ML / "marker_library_per_patient.csv")
        mld = ml[(ml.band == BAND) & (ml.marker == MARKER)]
        auc_val = dict(zip(mld.patient, mld.auc_resid))
    except Exception:
        auc_val = {}
    resp = [p for p in args.patients if auc_val.get(p, 0) >= RESP_THR]

    cframes, sframes = [], []
    for p in args.patients:
        cf, sf = per_patient(p, args.ndraw, rng)
        if cf is not None:
            cframes.append(cf)
        if sf is not None:
            sframes.append(sf)
    curve = pd.concat(cframes, ignore_index=True)
    summ = pd.concat(sframes, ignore_index=True)
    curve.to_csv(OUT / "recall_curve_long.csv", index=False)
    summ.to_csv(OUT / "recall_depth_summary.csv", index=False)

    print(f"[audit_115] {BAND} {MARKER} off-shaft recall GROWTH curve "
          f"({args.ndraw} draws, {time.time()-t0:.1f}s)")
    print(f"  responders: {resp}\n")
    for k in KSEEDS:
        cr = curve[(curve.k == k) & (curve.patient.isin(resp))]
        if cr.empty:
            continue
        g = cr.groupby("depth")[["recall_prop", "recall_str", "recall_rand"]].mean()
        print(f"  ── k={k} seeds — COMMUNITY RESPONDERS, recall(N) mean over patients:")
        print("    depth N:   5     10    15    20    30    40    50")
        row = lambda col: "  ".join(f"{g.loc[d, col]:.2f}" for d in [5, 10, 15, 20, 30, 40, 50])
        print(f"    propag :  {row('recall_prop')}")
        print(f"    strength: {row('recall_str')}")
        print(f"    random :  {row('recall_rand')}")
        sd = summ[(summ.k == k) & (summ.patient.isin(resp))]
        print("    depth to reach recovery (median over responders, % of pool read):")
        print(f"      50%: {sd.d50_frac.median()*100:4.0f}%   80%: {sd.d80_frac.median()*100:4.0f}%   "
              f"90%: {sd.d90_frac.median()*100:4.0f}%   100%: {sd.d100_frac.median()*100:4.0f}%   "
              f"(random = 50/80/90/100%)")
        print(f"      median pool = {sd.pool.median():.0f} contacts, "
              f"median #off-SOZ = {sd.n_off_soz.median():.0f}\n")

    print("  per-patient depth-to-90% recall (fraction of pool, k=3):")
    s3 = summ[summ.k == 3].sort_values("d90_frac")
    for _, r in s3.iterrows():
        print(f"    {r.patient}: {r.d90_frac*100:4.0f}%  (auc {r.auc:.2f}, pool {r.pool:.0f}, "
              f"#SOZ {r.n_off_soz:.0f})")
    print(f"\n[audit_115] -> {OUT}")


if __name__ == "__main__":
    main()
