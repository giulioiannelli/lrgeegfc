#!/usr/bin/env python3
"""Audit 114 — recovery/RECALL curve vs #seeds: from k marked SOZ, what FRACTION of the
remaining DISTANT (off-shaft) SOZ do we actually pull back inside a top-B shortlist? The
deployment complement of audit_113 (which reported AUC / break-even lift, not recovered-%).
Answers the PI question "can we recover the whole epi cohort, and which % as a function of
the seeds?" — honestly, off-shaft, delta heat kernel.

Marker = delta heat kernel e^{-tau L} at slow tau (heat_t5), strength-residualised (the
audit_101/102 winner). Proximity removed exactly as audit_101/113: never score a contact
that shares an electrode with any seed.

Critical preamble
=================
(1) Claim: ranking off-shaft contacts by strength-residual delta-heat affinity from k
    seeds recovers a measurable FRACTION (recall) of the remaining distant SOZ inside a
    top-B shortlist; recall rises with k and with B, above a random shortlist of size B.
(2) Null: recall@B = B/|off pool| (a random top-B list captures only the prevalence
    fraction) -> lift = 1. "Whole cohort" null: pooled recall = prevalence.
(3) Strongest alternative: hubness (killed by the strength residual AND an explicit
    strength-RANKED recall baseline at the same budget); proximity (killed by scoring
    ONLY contacts on electrodes carrying no seed).
(4) Reach: k seeds drawn at random per patient (NDRAW draws), off pool rebuilt per draw;
    recall at fixed budgets B in {5,10,20} and at break-even B=m and B=2m; reported
    per-patient + cohort, split community-responders vs all-10; the two hub patients
    (Pat_10/15) are reported, not hidden; cohort coverage given as BOTH macro (per-patient
    mean) and target-weighted pooled.
(5) Falsification: if recall ~ random (lift ~ 1) cohort-wide, there is no recovery and it
    is reported flat. No acceptance gate is pre-registered (post-hoc evaluation with PI).

Outputs (data/audit/epi_marker_recall/)
    recall_per_patient.csv ; recall_cohort.csv ; README.md
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

OUT = ROOT / "data" / "audit" / "epi_marker_recall"
ML = ROOT / "data" / "audit" / "epi_marker_library"
BAND = "delta"
MARKER = "heat_t5"
KSEEDS = [1, 2, 3, 5, 8]
BUDGETS = [5, 10, 20]    # fixed shortlist sizes (contacts a clinician would read)
MIN_OFF_CTRL = 3
RESP_THR = 0.568         # heat_t5 label-shuffle null p95 (audit_102) -> responder threshold
SEED = 20260619


def per_patient(pat, ndraw, rng):
    try:
        W = load_phase_fc(pat, "rest_post", BAND)
    except Exception:
        return None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None
    epi = np.asarray(pm.epi_mask, bool)
    probes = np.asarray(pm.probes, object)
    epi_idx = np.where(epi)[0]
    if epi_idx.size < max(MIN_EPI, 4) or len(set(probes[epi_idx])) < 2:
        return None
    ops, _intr, strength = build_operators(W)
    O = ops[MARKER]

    rows = []
    for k in KSEEDS:
        if epi_idx.size <= k:
            continue
        rec = {b: [] for b in BUDGETS}
        rec_str = {b: [] for b in BUDGETS}
        rec_rand = {b: [] for b in BUDGETS}
        rec_m, rec_2m, a_prop, m_list = [], [], [], []
        for _ in range(ndraw):
            seeds = rng.choice(epi_idx, size=k, replace=False)
            seed_pr = set(probes[seeds])
            seed_set = set(seeds.tolist())
            off = np.array([j for j in range(N) if probes[j] not in seed_pr])
            tgt = np.array([j for j in off if epi[j] and j not in seed_set])
            ctl = np.array([j for j in off if not epi[j]])
            if tgt.size < 1 or ctl.size < MIN_OFF_CTRL:
                continue
            aff = _resid(O[:, seeds].mean(axis=1), strength)         # strength-orthogonal
            c, t = _conc(aff[tgt], aff[ctl])
            a_prop.append(c / t if t else np.nan)
            tgt_set = set(tgt.tolist())
            m = tgt.size
            m_list.append(m)
            npool = off.size
            order = off[np.argsort(-aff[off])]                       # propagator ranking
            order_str = off[np.argsort(-strength[off])]              # hubness baseline ranking
            for b in BUDGETS:
                rec[b].append(sum(1 for j in order[:b] if j in tgt_set) / m)
                rec_str[b].append(sum(1 for j in order_str[:b] if j in tgt_set) / m)
                rec_rand[b].append(min(b, npool) / npool)            # E[random recall@B] = B/|pool|
            rec_m.append(sum(1 for j in order[:m] if j in tgt_set) / m)          # break-even
            rec_2m.append(sum(1 for j in order[:2 * m] if j in tgt_set) / m)
        if not a_prop:
            continue
        row = dict(patient=pat, k=k, n_draws=len(a_prop),
                   auc_prop=float(np.nanmean(a_prop)),
                   n_off_soz=float(np.mean(m_list)),
                   recall_bm=float(np.nanmean(rec_m)),
                   recall_b2m=float(np.nanmean(rec_2m)))
        for b in BUDGETS:
            row[f"recall_b{b}"] = float(np.nanmean(rec[b]))
            row[f"recall_str_b{b}"] = float(np.nanmean(rec_str[b]))
            row[f"recall_rand_b{b}"] = float(np.nanmean(rec_rand[b]))
        rows.append(row)
    return pd.DataFrame(rows) if rows else None


def cohort(df, resp):
    out = []
    base = [f"recall_b{b}" for b in BUDGETS] + ["recall_bm", "recall_b2m"]
    aux = [f"recall_rand_b{b}" for b in BUDGETS] + [f"recall_str_b{b}" for b in BUDGETS]
    for k in KSEEDS:
        d = df[df.k == k]
        if d.empty:
            continue
        dr = d[d.patient.isin(resp)]
        row = dict(k=k, n_all=int(d.patient.nunique()), n_resp=int(dr.patient.nunique()),
                   n_off_soz=float(d.n_off_soz.median()))
        for c in base + aux:
            row[f"{c}_all"] = float(d[c].mean())
            row[f"{c}_resp"] = float(dr[c].mean()) if len(dr) else np.nan
        if len(dr):
            w = dr.n_off_soz.values
            row["recall_b10_resp_pooled"] = float(np.average(dr.recall_b10.values, weights=w))
            row["n_resp_above_rand10"] = int((dr.recall_b10 > dr.recall_rand_b10).sum())
        out.append(row)
    return pd.DataFrame(out)


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

    frames = [per_patient(p, args.ndraw, rng) for p in args.patients]
    df = pd.concat([f for f in frames if f is not None], ignore_index=True)
    df.to_csv(OUT / "recall_per_patient.csv", index=False)
    coh = cohort(df, resp)
    coh.to_csv(OUT / "recall_cohort.csv", index=False)

    print(f"[audit_114] {BAND} {MARKER} off-shaft RECALL curve ({args.ndraw} draws, "
          f"{df.patient.nunique()} patients, {time.time()-t0:.1f}s)")
    print(f"  responders (heat_t5 LOSO AUC>={RESP_THR}): {resp}\n")
    print("  COMMUNITY RESPONDERS — recall (% of distant SOZ recovered), mean over patients:")
    print("   k   top5   top10  top20  break-even(=m)  | rand@10  strength@10  median #off-SOZ")
    for _, r in coh.iterrows():
        print(f"   {int(r.k)}   {r.recall_b5_resp:.2f}   {r.recall_b10_resp:.2f}   "
              f"{r.recall_b20_resp:.2f}   {r.recall_bm_resp:.2f}            | "
              f"{r.recall_rand_b10_resp:.2f}    {r.recall_str_b10_resp:.2f}        {r.n_off_soz:.0f}")
    print("\n  ALL 10 — recall@top-10, mean over patients:")
    for _, r in coh.iterrows():
        print(f"   k={int(r.k)}: {r.recall_b10_all:.2f}   (responders {r.recall_b10_resp:.2f}, "
              f"random {r.recall_rand_b10_all:.2f})")
    print("\n  per-patient recall@top-10:")
    piv = df.pivot(index="patient", columns="k", values="recall_b10")
    print(piv.round(2).to_string())
    print(f"\n[audit_114] -> {OUT}")


if __name__ == "__main__":
    main()
