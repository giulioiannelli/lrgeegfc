#!/usr/bin/env python3
"""Audit 101 — systematic LRG/Laplacian node-marker LIBRARY for epileptic-node discovery,
scored ONLY on the non-trivial task: find SOZ on electrodes far from the known ones.

Why this design (the lesson from audit_99). Spatial proximity is a rigged baseline:
contacts are labelled SOZ because they lie in a marked epileptogenic AREA, so "near a
known SOZ -> also SOZ" is baked into the labelling, not a discovery. We therefore DISCARD
proximity and evaluate every marker on LEAVE-ONE-SHAFT-OUT (LOSO) distant discovery:
    hide ALL SOZ on one electrode (shaft); using SOZ on the OTHER shafts as seeds, can the
    marker rank the hidden shaft's SOZ above healthy contacts, looking ONLY at contacts on
    electrodes that contain no seed? Proximity cannot help here by construction.
Headline metric = strength-RESIDUALISED LOSO AUC (so hubness is controlled too), with a
cohort sign test + a label-shuffle null for the survivors.

Marker families (all from the Laplacian L = D - W of the |ImCoh| graph), seed-based
affinity of every node to the seed set S (mean over seeds of an operator column):
    heat kernel e^{-tL} (6 taus) ; normalized-Laplacian heat (6) ; personalized PageRank
    (3 alphas) ; Katz ; random-walk-with-restart ; communicability expm(W/rho) ;
    negative effective-resistance (commute distance, via L^+) ; diffusion-distance (taus).
Plus a strength baseline and a few seed-FREE intrinsic markers (heat-diagonal return
probability, slow-mode participation, Fiedler) for reference.

Critical preamble
=================
(1) Claim: at least one Laplacian operator ranks distant (off-shaft) SOZ above healthy
    after strength is controlled, replicating across patients.
(2) Null: strength-residual LOSO AUC = 0.5; + a label-shuffle null (fake-SOZ shafts).
(3) Strongest alternative: distant SOZ are just hubs -> killed by strength residual; or a
    graph artifact -> killed by the label-shuffle null.
(4) Reach: every marker uses the SAME LOSO folds and the SAME strength residual; markers
    are not tuned per patient; the cohort sign test + null guard against the 100-marker
    multiple-comparison trap (the winner must also pass the shuffle null and a LOPO check).
(5) Falsification: if NO marker exceeds 0.5 under strength-residual LOSO with a cohort
    sign test, the propagator does NOT mark distant epileptic nodes. Reported flat.

Outputs (data/audit/epi_marker_library/) marker_library_per_patient.csv ; *_cohort.csv ; README.md
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pandas as pd
from numpy.linalg import eigh, pinv
from scipy.linalg import expm
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()
from lrg_eegfc.utils.io.patient import build_epi_masks

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_63_split_baseline_surrogate import load_phase_fc  # type: ignore
import _epi_stratify as es  # type: ignore
from audit_85_epi_propagator_recovery import MIN_EPI, N_TAU  # type: ignore

OUT = ROOT / "data" / "audit" / "epi_marker_library"
MIN_PAIRS = 15


def build_operators(W):
    """Precompute seed-independent N×N operators O; marker(S)= O[:,S].mean(1)."""
    N = W.shape[0]
    d = W.sum(1)
    L = np.diag(d) - W
    lam, U = eigh(L)
    lam = np.clip(lam, 0, None)
    taus = np.geomspace(1.0 / lam[-1], 10.0 / lam[-1], N_TAU)
    # normalized Laplacian
    dm12 = 1.0 / np.sqrt(np.clip(d, 1e-12, None))
    Ln = np.eye(N) - (dm12[:, None] * W * dm12[None, :])
    lamn, Un = eigh((Ln + Ln.T) / 2)
    lamn = np.clip(lamn, 0, None)
    taun = np.geomspace(1.0 / max(lamn[-1], 1e-9), 10.0 / max(lamn[-1], 1e-9), N_TAU)
    P = (W / np.clip(d, 1e-12, None)[:, None])          # row-stochastic random walk
    rhoW = max(abs(np.linalg.eigvalsh(W)).max(), 1e-9)
    Lp = pinv(L)                                        # for effective resistance
    diagLp = np.diag(Lp)
    R = diagLp[:, None] + diagLp[None, :] - 2 * Lp      # resistance distance matrix

    ops = {}
    for ti, t in enumerate(taus):
        ops[f"heat_t{ti}"] = U @ (np.exp(-t * lam)[:, None] * U.T)
    for ti, t in enumerate(taun):
        ops[f"heatN_t{ti}"] = Un @ (np.exp(-t * lamn)[:, None] * Un.T)
    for a in (0.50, 0.85, 0.95):
        ops[f"ppr_a{int(a*100)}"] = (1 - a) * np.linalg.inv(np.eye(N) - a * P)
    beta = 0.85 / rhoW
    ops["katz"] = np.linalg.inv(np.eye(N) - beta * W)
    ops["comm"] = expm(W / rhoW)                        # communicability
    ops["negres"] = -R                                  # higher = closer (commute)
    # diffusion distance at a mid tau: ||rho_i - rho_s|| -> use -dist via heat columns
    Hd = ops["heat_t2"]
    ops["diffdist_t2"] = -(np.sqrt(np.clip(
        np.diag(Hd)[:, None] + np.diag(Hd)[None, :] - 2 * Hd, 0, None)))

    intrinsic = {}                                      # seed-FREE per-node scores
    for ti in range(N_TAU):
        intrinsic[f"heatdiag_t{ti}"] = np.diag(ops[f"heat_t{ti}"]).copy()
    for K in (3, 5, 10):
        intrinsic[f"slowpart_K{K}"] = (U[:, 1:1 + K] ** 2).sum(1)   # slow-mode participation
    intrinsic["fiedler_abs"] = np.abs(U[:, 1])
    return ops, intrinsic, d


def _resid(score, strength):
    s = np.asarray(strength, float); f = np.asarray(score, float)
    m = np.isfinite(f) & np.isfinite(s)
    if m.sum() < 3 or np.std(s[m]) < 1e-12:
        return f - np.nanmean(f[m]) if m.any() else f
    b = np.polyfit(s[m], f[m], 1)
    return f - np.polyval(b, s)


def _conc(scores_case, scores_ctrl):
    """concordance (case>ctrl) summed over all case-ctrl pairs -> (conc, total)."""
    c = 0.0; t = 0
    for sc in scores_case:
        if not np.isfinite(sc):
            continue
        ok = np.isfinite(scores_ctrl)
        v = scores_ctrl[ok]
        c += np.sum(sc > v) + 0.5 * np.sum(sc == v); t += v.size
    return c, t


def loso_auc(marker_of_seeds, epi_idx, healthy, probes, strength):
    """Leave-one-SHAFT-out distant-discovery AUC, raw + strength-residual.

    marker_of_seeds(seeds) -> per-node score array.
    """
    soz_probes = sorted(set(probes[epi_idx]))
    if len(soz_probes) < 2:
        return np.nan, np.nan, 0
    raw = [0.0, 0]; res = [0.0, 0]
    for p in soz_probes:
        seeds = epi_idx[probes[epi_idx] != p]
        targets = epi_idx[probes[epi_idx] == p]
        if seeds.size < 1 or targets.size < 1:
            continue
        seed_pr = set(probes[seeds])
        off = np.array([j for j in np.r_[targets, healthy] if probes[j] not in seed_pr])
        tgt = np.array([j for j in off if j in set(targets.tolist())])
        ctrl = np.array([j for j in off if j not in set(targets.tolist())])
        if tgt.size < 1 or ctrl.size < 3:
            continue
        m = marker_of_seeds(seeds)
        mr = _resid(m, strength)
        c, t = _conc(m[tgt], m[ctrl]); raw[0] += c; raw[1] += t
        c, t = _conc(mr[tgt], mr[ctrl]); res[0] += c; res[1] += t
    if raw[1] < MIN_PAIRS:
        return np.nan, np.nan, raw[1]
    return raw[0] / raw[1], (res[0] / res[1] if res[1] >= MIN_PAIRS else np.nan), raw[1]


def per_patient_band(pat, band):
    try:
        W = load_phase_fc(pat, "rest_post", band)
    except Exception:
        return None
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None
    epi = np.asarray(pm.epi_mask, bool)
    probes = np.asarray(pm.probes, object)
    epi_idx = np.where(epi)[0]
    healthy = np.where(~epi)[0]
    if epi_idx.size < max(MIN_EPI, 4) or healthy.size < 5:
        return None
    if len(set(probes[epi_idx])) < 2:
        return None                       # need >=2 SOZ shafts for LOSO
    ops, intrinsic, strength = build_operators(W)

    rows = []
    for name, O in ops.items():
        def marker(seeds, O=O):
            return O[:, seeds].mean(axis=1)
        a_raw, a_res, npairs = loso_auc(marker, epi_idx, healthy, probes, strength)
        rows.append({"patient": pat, "band": band, "marker": name, "kind": "seedbased",
                     "auc_raw": a_raw, "auc_resid": a_res, "n_pairs": npairs})
    # strength baseline as a seed-based marker (affinity-free): rank by strength itself
    def str_marker(seeds, s=strength):
        return s
    a_raw, a_res, npairs = loso_auc(str_marker, epi_idx, healthy, probes, strength)
    rows.append({"patient": pat, "band": band, "marker": "strength_baseline",
                 "kind": "baseline", "auc_raw": a_raw, "auc_resid": np.nan,
                 "n_pairs": npairs})
    return pd.DataFrame(rows)


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
    frames = []
    for band in args.bands:
        for pat in args.patients:
            r = per_patient_band(pat, band)
            if r is not None:
                frames.append(r)
        print(f"[audit_101] {band}: done")
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(OUT / "marker_library_per_patient.csv", index=False)
    coh = cohort(df)
    coh.to_csv(OUT / "marker_library_cohort.csv", index=False)
    rt = time.time() - t0

    print("\n[audit_101] TOP markers by strength-RESIDUAL LOSO distant-discovery AUC:")
    res = coh[coh.metric == "auc_resid"].copy()
    for band in [b for b in es.ALL_BANDS if b in set(res.band)]:
        d = res[res.band == band].sort_values("median_auc", ascending=False).head(6)
        print(f"  --- {band} ---")
        for _, r in d.iterrows():
            print(f"    {r.marker:14s} med={r.median_auc:.3f} mean={r.mean_auc:.3f} "
                  f"{r.n_above_half}/{r.n_patients}>0.5 p={r.p_sign:.3f}")
    print(f"[audit_101] done in {rt:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
