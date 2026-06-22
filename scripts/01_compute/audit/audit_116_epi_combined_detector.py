#!/usr/bin/env python3
"""Audit 116 — the diffusion-ONLY SOZ detector, evaluated on ALL held-out SOZ (same-shaft
INCLUDED), with NO spatial proximity anywhere in the method. A real detector must find the
nearby SOZ too — so we keep them as targets — but it must find them from the LRG propagator
alone (connectivity), never from implant coordinates. Proximity is a trivial label-bias; it
is not a feature here.

This corrects BOTH earlier errors:
  - off-shaft (audit_101/113/115) wrongly DELETED same-shaft SOZ from the targets -> only ~4
    targets, no power, not a detector.
  - the proximity-feature draft wrongly ADDED coordinates as a predictor -> trivial bias.
The honest design: targets = every remaining SOZ; method = diffusion affinity to the seeds
(+ a multiscale combination of all 6 heat-kernel scales); baselines to beat = node strength
(hubness), chance, and a label-shuffle null (count-matched FAKE-SOZ set -> must collapse to
chance, proving the detector tracks the TRUE SOZ community, not an arbitrary blob).

Rankers (every one is coordinate-free; scores all non-seed contacts):
    diff_slow       affinity to seeds at the slow heat scale  e^{-tau L}[:,S].mean      (marker)
    diff_resid      diff_slow residualised on node strength   (beyond-hubness validity)
    diff_multiscale z-sum of affinity across all 6 heat scales (uses MORE of the propagator)
    strength        node strength D_ii                         (hubness baseline to beat)
    labelshuffle    diff_multiscale from a random count-matched FAKE-SOZ set  (null)
    chance          prevalence floor

Critical preamble
=================
(1) Claim: the diffusion marker, coordinate-free, ranks ALL held-out SOZ (nearby + distant)
    above healthy with usable precision, beating strength and its own label-shuffle; a
    multiscale combination beats the single slow scale (more propagator information helps).
(2) Null: detector AUC = 0.5 / precision = prevalence; label-shuffle = detector (no SOZ-
    specific community); multiscale = slow (no extra scale information).
(3) Strongest alternative: hubness -> strength is an explicit ablation rung AND diff_resid
    removes it; "any spatial blob lights up" -> the label-shuffle null is exactly that test;
    proximity-via-coordinates -> not used at all.
(4) Reach: ALL held-out SOZ targeted (same-shaft included); k random seeds, N_DRAW draws;
    multiscale uses equal-weight z-sum (no fit -> no overfit); hub patients (Pat_10/15) kept;
    improvement test paired across patients.
(5) Falsification: if diff AUC ~ 0.5 / ~ label-shuffle, the propagator is not a detector and
    it is reported flat; if multiscale <= slow, extra scales add nothing. No pre-registered
    acceptance gate.

Outputs (data/audit/epi_combined_detector/)
    detector_per_patient.csv ; detector_cohort.csv ; improvement_test.csv ; README.md
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
from audit_85_epi_propagator_recovery import MIN_EPI, N_TAU  # type: ignore
from audit_99_epi_discovery_yield import _yield_curve  # type: ignore
from audit_101_epi_marker_library import build_operators, _resid, _conc  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_combined_detector"
KSEEDS = [3, 5]
BUDGETS = np.array([5, 10, 20])
N_DRAW = 120
SEED = 20260619
RANKERS = ["diff_slow", "diff_resid", "diff_multiscale", "strength", "labelshuffle"]


def _zc(v):
    v = np.asarray(v, float)
    ok = np.isfinite(v)
    if ok.sum() < 2 or v[ok].std() < 1e-12:
        return np.zeros(v.size)
    z = np.zeros(v.size)
    z[ok] = (v[ok] - v[ok].mean()) / v[ok].std()
    z[~ok] = z[ok].min()
    return z


def _auc(score, is_t):
    c, t = _conc(np.asarray(score, float)[is_t], np.asarray(score, float)[~is_t])
    return c / t if t else np.nan


def per_patient(pat, band, ndraw, rng):
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
    if epi_idx.size < max(MIN_EPI, max(KSEEDS) + 1) or healthy.size < 5:
        return None
    ops, _intr, strength = build_operators(W)
    heats = [ops[f"heat_t{i}"] for i in range(N_TAU)]
    all_idx = np.arange(N)

    rows = []
    for k in KSEEDS:
        if k >= epi_idx.size:
            continue
        acc = {r: {m: [] for m in ("auc", "p5", "p10", "r10", "r20")} for r in RANKERS}
        for _ in range(ndraw):
            seeds = rng.choice(epi_idx, size=k, replace=False)
            targets = np.array([e for e in epi_idx if e not in set(seeds.tolist())])
            pool = np.r_[targets, healthy]
            is_t = np.r_[np.ones(targets.size, bool), np.zeros(healthy.size, bool)]
            affs = [H[:, seeds].mean(axis=1)[pool] for H in heats]
            scores = {
                "diff_slow": affs[-1],
                "diff_resid": _resid(heats[-1][:, seeds].mean(axis=1), strength)[pool],
                "diff_multiscale": np.sum([_zc(a) for a in affs], axis=0),
                "strength": strength[pool],
            }
            # label-shuffle null: a random count-matched FAKE-SOZ set, same multiscale detector
            fake = rng.choice(all_idx, size=epi_idx.size, replace=False)
            fseeds, ftargets = fake[:k], fake[k:]
            fhealth = np.array([j for j in all_idx if j not in set(fake.tolist())])
            fpool = np.r_[ftargets, fhealth]
            fis_t = np.r_[np.ones(ftargets.size, bool), np.zeros(fhealth.size, bool)]
            faffs = [H[:, fseeds].mean(axis=1)[fpool] for H in heats]
            for r in RANKERS:
                if r == "labelshuffle":
                    sc, it = np.sum([_zc(a) for a in faffs], axis=0), fis_t
                else:
                    sc, it = scores[r], is_t
                acc[r]["auc"].append(_auc(sc, it))
                rec, prec = _yield_curve(sc, it, BUDGETS)
                acc[r]["p5"].append(prec[0]); acc[r]["p10"].append(prec[1])
                acc[r]["r10"].append(rec[1]); acc[r]["r20"].append(rec[2])
        for r in RANKERS:
            rows.append(dict(patient=pat, band=band, k=k, ranker=r,
                             n_soz=int(epi_idx.size), n_draws=ndraw,
                             auc=float(np.nanmean(acc[r]["auc"])),
                             prec5=float(np.nanmean(acc[r]["p5"])),
                             prec10=float(np.nanmean(acc[r]["p10"])),
                             rec10=float(np.nanmean(acc[r]["r10"])),
                             rec20=float(np.nanmean(acc[r]["r20"]))))
        prev = targets.size / pool.size  # using last draw's sizes (stable across draws)
        rows.append(dict(patient=pat, band=band, k=k, ranker="chance", n_soz=int(epi_idx.size),
                         n_draws=ndraw, auc=0.5, prec5=prev, prec10=prev,
                         rec10=float(BUDGETS[1] / pool.size), rec20=float(BUDGETS[2] / pool.size)))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=["delta"])
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    ap.add_argument("--ndraw", type=int, default=N_DRAW)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    t0 = time.time()

    frames = [per_patient(p, b, args.ndraw, rng) for b in args.bands for p in args.patients]
    df = pd.concat([f for f in frames if f is not None], ignore_index=True)
    df.to_csv(OUT / "detector_per_patient.csv", index=False)
    coh = (df.groupby(["band", "k", "ranker"])
             [["auc", "prec5", "prec10", "rec10", "rec20"]].mean().reset_index())
    coh.to_csv(OUT / "detector_cohort.csv", index=False)

    inc = []
    for band in args.bands:
        for k in KSEEDS:
            piv = df[(df.band == band) & (df.k == k)].pivot(
                index="patient", columns="ranker", values="auc")
            for base, comb, lbl in [("diff_slow", "diff_multiscale", "multiscale − slow"),
                                    ("strength", "diff_slow", "diff_slow − strength"),
                                    ("labelshuffle", "diff_multiscale", "multiscale − shuffle")]:
                if not {base, comb}.issubset(piv.columns):
                    continue
                d = (piv[comb] - piv[base]).dropna()
                try:
                    p = float(wilcoxon(d, alternative="greater")[1]) if (d != 0).any() else 1.0
                except ValueError:
                    p = np.nan
                inc.append(dict(band=band, k=k, test=lbl, median_delta_auc=float(d.median()),
                                n_pos=int((d > 0).sum()), n=int(d.size), p_greater=p))
    inc = pd.DataFrame(inc)
    inc.to_csv(OUT / "improvement_test.csv", index=False)

    print(f"[audit_116] diffusion-ONLY detector, ALL held-out SOZ, NO proximity "
          f"({time.time()-t0:.1f}s)\n")
    for band in args.bands:
        for k in KSEEDS:
            c = coh[(coh.band == band) & (coh.k == k)].set_index("ranker")
            if c.empty:
                continue
            print(f"  ── {band} k={k} — cohort mean (n=10), ALL contacts (no coords):")
            print("    ranker            AUC    prec@5  prec@10  rec@10  rec@20")
            for r in RANKERS + ["chance"]:
                if r in c.index:
                    x = c.loc[r]
                    print(f"    {r:15s}  {x.auc:.3f}  {x.prec5*100:5.1f}%  {x.prec10*100:5.1f}%  "
                          f"{x.rec10*100:5.1f}%  {x.rec20*100:5.1f}%")
            for _, r in inc[(inc.band == band) & (inc.k == k)].iterrows():
                print(f"    Δ {r.test}: median {r.median_delta_auc:+.3f} AUC, "
                      f"{r.n_pos}/{r.n}>0, p={r.p_greater:.3f}")
            print()
    print(f"[audit_116] -> {OUT}")


if __name__ == "__main__":
    main()
