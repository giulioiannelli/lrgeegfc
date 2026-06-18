#!/usr/bin/env python3
"""Audit 99 — DISCOVERY YIELD: mark a few SOZ, how many of the rest do you find per
contact reviewed? The clinician-facing performance test for the seed-based marker.

This is the deployment question stated as a yield curve: a clinician marks k seed SOZ,
the marker ranks every other contact, and we ask — if they review the top-N suggestions,
what fraction of the *remaining* SOZ do they recover (recall@N), and how pure is the list
(precision@N) — against the three baselines a skeptic will demand, plus a label-shuffle
null. Reuses the validated affinity math (audit_91/92), forks nothing.

Critical preamble
=================
(1) Claim: from k labelled SOZ, the strength-residual propagator-affinity ranking finds
    the remaining SOZ with higher recall-per-contact-reviewed than chance, node strength,
    OR spatial proximity to the seeds.
(2) Null: random review (recall@N = N/|pool|, precision = prevalence).
(3) Strongest alternatives the curve must beat, each a distinct objection:
    - node STRENGTH  → "you're re-reading hubness";
    - DISTANCE-to-nearest-seed → "epilepsy is focal, you just flag nearby contacts";
    - LABEL-SHUFFLE (relabel a random count-matched set as fake-SOZ, rerun) → "any marked
      set is recoverable; it's a graph artifact, not the real SOZ".
(4) Reach: marker = τ-mean affinity to seeds, residualised on strength (so the strength
    line is the fair hubness control); distance uses implant (x,y,z); random seeds are
    averaged over many draws; the clinically-seeded scheme (top-k strength SOZ) is a
    realistic-deployment companion, not the primary.
(5) Falsification: if the marker curve sits on the random/strength/distance curves, or if
    the label-shuffle curve matches the real one, the marker does NOT discover SOZ — it is
    chance, hubness, proximity, or a structural artifact.

Outputs (data/audit/epi_discovery_yield/)
    discovery_yield_per_patient.csv   (patient,band,k,scheme,method,budget,recall,precision)
    discovery_yield_cohort.csv        cohort mean/median per (band,k,scheme,method,budget)
    README_discovery_yield.md
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
from lrg_eegfc.utils.io.regions import load_channel_regions

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI, N_TAU  # type: ignore
from audit_92_epi_masked_recovery_relational import _aff_to_set  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_discovery_yield"
K_SEEDS = (2, 3, 5)
N_MAX = 20                  # review-budget axis: top-1 .. top-20 contacts
N_DRAW = 80                 # random seed draws (and label-shuffle draws)
BUDGETS = np.arange(1, N_MAX + 1)


def _yield_curve(score, is_target, budgets):
    """recall@N and precision@N over the candidate pool for each budget N."""
    s = np.asarray(score, float)
    t = np.asarray(is_target, bool)
    ok = np.isfinite(s)
    s, t = s[ok], t[ok]
    T = int(t.sum())
    if T == 0 or s.size == 0:
        return np.full(budgets.size, np.nan), np.full(budgets.size, np.nan)
    order = np.argsort(-s)
    hits = np.cumsum(t[order].astype(float))
    rec, prec = [], []
    for N in budgets:
        n = min(N, hits.size)
        rec.append(hits[n - 1] / T)
        prec.append(hits[n - 1] / N)
    return np.array(rec), np.array(prec)


def _dist_score(coords, pool, seeds):
    """Spatial baseline: higher = closer to nearest seed (−min distance)."""
    sc = np.full(coords.shape[0], -np.inf)
    cs = coords[seeds]
    for i in pool:
        if np.isfinite(coords[i]).all():
            sc[i] = -float(np.nanmin(np.linalg.norm(cs - coords[i], axis=1)))
    return sc


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
    probe = np.asarray(pm.probes, object)
    epi_idx = np.where(epi)[0]
    healthy = np.where(~epi)[0]
    if epi_idx.size < max(MIN_EPI, max(K_SEEDS) + 1) or healthy.size < 3:
        return None
    lam, U = np.linalg.eigh(np.diag(W.sum(1)) - W)
    if float(lam[-1]) <= 0:
        return None
    taus = np.geomspace(1.0 / lam[-1], 10.0 / lam[-1], N_TAU)
    strength = W.sum(1)
    reg = load_channel_regions(pat)
    coords = reg[["x", "y", "z"]].to_numpy(float)
    all_idx = np.arange(N)

    rows = []

    def _emit(k, scheme, method, rec, prec):
        for bi, Nb in enumerate(BUDGETS):
            rows.append({"patient": pat, "band": band, "k": k, "scheme": scheme,
                         "method": method, "budget": int(Nb),
                         "recall": float(rec[bi]), "precision": float(prec[bi])})

    for k in K_SEEDS:
        if k >= epi_idx.size:
            continue
        # ---- RANDOM seeds (primary): average curves over N_DRAW draws ----
        acc = {m: {"rec": [], "prec": []} for m in
               ("marker", "strength", "distance", "labelshuffle")}
        # OFF-SHAFT stratum: only contacts on a DIFFERENT probe than every seed
        # (where spatial proximity is useless — the marker's one shot at unique value)
        acc_off = {m: {"rec": [], "prec": []} for m in ("marker", "strength", "distance")}
        n_off_draws = 0
        for _ in range(N_DRAW):
            seeds = rng.choice(epi_idx, size=k, replace=False)
            targets = np.array([e for e in epi_idx if e not in set(seeds)])
            pool = np.r_[targets, healthy]
            is_t = np.r_[np.ones(targets.size, bool), np.zeros(healthy.size, bool)]
            _fr, f_res = _aff_to_set(lam, U, taus, seeds, strength)
            for m, sc in (("marker", f_res), ("strength", strength),
                          ("distance", _dist_score(coords, pool, seeds))):
                r, p = _yield_curve(sc[pool], is_t, BUDGETS)
                acc[m]["rec"].append(r); acc[m]["prec"].append(p)
            # off-shaft stratum
            seed_pr = set(probe[seeds])
            off = pool[np.array([probe[j] not in seed_pr for j in pool], bool)]
            off_t = np.array([j in set(targets.tolist()) for j in off], bool)
            if off.size and off_t.sum() >= 1:
                n_off_draws += 1
                for m, sc in (("marker", f_res), ("strength", strength),
                              ("distance", _dist_score(coords, off, seeds))):
                    r, p = _yield_curve(sc[off], off_t, BUDGETS)
                    acc_off[m]["rec"].append(r); acc_off[m]["prec"].append(p)
            # label-shuffle null: a random count-matched FAKE-SOZ set
            fake = rng.choice(all_idx, size=epi_idx.size, replace=False)
            fseeds = fake[:k]
            ftargets = fake[k:]
            fhealth = np.array([j for j in all_idx if j not in set(fake)])
            fpool = np.r_[ftargets, fhealth]
            fis_t = np.r_[np.ones(ftargets.size, bool), np.zeros(fhealth.size, bool)]
            _fr2, ff = _aff_to_set(lam, U, taus, fseeds, strength)
            r, p = _yield_curve(ff[fpool], fis_t, BUDGETS)
            acc["labelshuffle"]["rec"].append(r); acc["labelshuffle"]["prec"].append(p)
        for m in acc:
            _emit(k, "random", m,
                  np.nanmean(acc[m]["rec"], axis=0), np.nanmean(acc[m]["prec"], axis=0))
        if n_off_draws >= 0.2 * N_DRAW:        # only if off-shaft SOZ exist for this patient
            for m in acc_off:
                if acc_off[m]["rec"]:
                    _emit(k, "random_offshaft", m,
                          np.nanmean(acc_off[m]["rec"], axis=0),
                          np.nanmean(acc_off[m]["prec"], axis=0))
        # random-review reference (analytic)
        prev = (epi_idx.size - k) / ((epi_idx.size - k) + healthy.size)
        _emit(k, "random", "chance",
              BUDGETS / ((epi_idx.size - k) + healthy.size), np.full(BUDGETS.size, prev))

        # ---- CLINICALLY-seeded (companion): top-k strength SOZ, deterministic ----
        seeds = epi_idx[np.argsort(-strength[epi_idx])][:k]
        targets = np.array([e for e in epi_idx if e not in set(seeds)])
        pool = np.r_[targets, healthy]
        is_t = np.r_[np.ones(targets.size, bool), np.zeros(healthy.size, bool)]
        _fr, f_res = _aff_to_set(lam, U, taus, seeds, strength)
        for m, sc in (("marker", f_res), ("strength", strength),
                      ("distance", _dist_score(coords, pool, seeds))):
            r, p = _yield_curve(sc[pool], is_t, BUDGETS)
            _emit(k, "clinical", m, r, p)
    return pd.DataFrame(rows)


def cohort(df):
    g = (df.groupby(["band", "k", "scheme", "method", "budget"])
           .agg(n_patients=("patient", "nunique"),
                mean_recall=("recall", "mean"), med_recall=("recall", "median"),
                mean_precision=("precision", "mean"), med_precision=("precision", "median"))
           .reset_index())
    return g


def write_readme(coh, runtime):
    L = ["---", "name: epi_discovery_yield",
         "scope: seed_based_soz_discovery_yield_curve_with_baselines_and_labelshuffle_null",
         "era: COHORT_N10 / IMCOH_ABS",
         f"date: {time.strftime('%Y-%m-%d')}",
         "build_script: scripts/01_compute/audit/audit_99_epi_discovery_yield.py",
         "---", "",
         "# Discovery yield — mark k SOZ, recover the rest per contact reviewed", "",
         "**Head.** recall@N / precision@N for the propagator marker vs node-strength, "
         "distance-to-seed, a random-review floor, and a label-shuffle null (random "
         "count-matched fake-SOZ). Random seeds (primary, averaged over draws) + "
         "clinically-seeded (top-k strength SOZ) companion. The marker is good iff its "
         "curve beats strength AND distance AND the label-shuffle collapses to chance.", "",
         "## Recall @ review-budget (random seeds, cohort mean)", ""]
    for band in [b for b in es.ALL_BANDS if b in set(coh.band)]:
        for k in K_SEEDS:
            sub = coh[(coh.band == band) & (coh.k == k) & (coh.scheme == "random")]
            if sub.empty:
                continue
            L.append(f"### {band}, k={k}")
            L.append("| budget N | marker | strength | distance | label-shuffle | chance |")
            L.append("|---|---|---|---|---|---|")
            for Nb in (5, 10, 20):
                def gv(m):
                    r = sub[(sub.method == m) & (sub.budget == Nb)]
                    return f"{r.mean_recall.iloc[0]:.2f}" if not r.empty else "—"
                L.append(f"| {Nb} | {gv('marker')} | {gv('strength')} | {gv('distance')} "
                         f"| {gv('labelshuffle')} | {gv('chance')} |")
            L.append("")
    L += [f"- random seed draws per k: {N_DRAW}; budgets 1..{N_MAX}; "
          f"wall-clock {runtime:.1f}s"]
    (OUT / "README_discovery_yield.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=["delta", "beta", "low_gamma", "alpha"])
    ap.add_argument("--patients", nargs="+", default=list(es.COHORT))
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    frames = []
    for band in args.bands:
        for pat in args.patients:
            rng = np.random.default_rng(20260612 + sum(ord(c) for c in pat + band))
            r = per_patient_band(pat, band, rng)
            if r is not None:
                frames.append(r)
                print(f"[audit_99] {pat}/{band}: ok")
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(OUT / "discovery_yield_per_patient.csv", index=False)
    coh = cohort(df)
    coh.to_csv(OUT / "discovery_yield_cohort.csv", index=False)
    runtime = time.time() - t0
    write_readme(coh, runtime)
    print(f"\n[audit_99] recall@10 (random seeds, cohort mean): marker | strength | "
          f"distance | label-shuffle")
    for band in [b for b in es.ALL_BANDS if b in set(coh.band)]:
        for k in K_SEEDS:
            sub = coh[(coh.band == band) & (coh.k == k) & (coh.scheme == "random")
                      & (coh.budget == 10)]
            if sub.empty:
                continue
            def gv(m):
                r = sub[sub.method == m]
                return f"{r.mean_recall.iloc[0]:.2f}" if not r.empty else "—"
            print(f"  {band:10s} k={k}: {gv('marker')} | {gv('strength')} | "
                  f"{gv('distance')} | {gv('labelshuffle')}")
    print(f"[audit_99] done in {runtime:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
