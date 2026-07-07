#!/usr/bin/env python3
"""Audit 126 — N4 replay states RESCUE: unsupervised state decomposition.

audit_125 tested replay by template-matching each window to the full-phase task
config and found a cohort NEGATIVE (no post-specific burstiness; rest_pre as
task-like as rest_post). This script tries the genuinely different "states"
operationalization used in the dynamic-FC / replay literature: discover RECURRING
connectivity states from the data, then ask whether a TASK-LIKE state recurs MORE
in rest_post than rest_pre (occupancy + recurrence). Reuses the per-window
cophenetic cache from audit_125 (recipe == N1, anchored audit_124).

Two methods, per patient (within-subject; cohort aggregates the post-pre contrast):

METHOD 1 — nearest-centroid (in-framework, primary). Per window discriminant
    g_w = rho^coph(c_w, c_task_full) - rho^coph(c_w, c_pre_full)
  (how much more task-like than rest-like, both via Spearman = the rho^coph
  convention). Multiband = mean over {alpha, beta} of per-band z-scored g.
  TASK-STATE = g_w above the pooled(pre+post) median (so neither phase privileged).
  occupancy_pre/post = fraction of windows in the task-state; recurrence = number
  of separate task-state visits (runs).

METHOD 2 — k-means state discovery (secondary). Per-band column-standardized
  cophenetic vectors concatenated across {alpha, beta} -> PCA(10) -> KMeans(K).
  Fit on ALL windows (pre+task+post). TASK-STATE = cluster where task_test windows
  most concentrate (lift over chance reported). occupancy_pre/post as above.

REPLAY (scope N4.4): occupancy_post > occupancy_pre. Cohort one-sided Wilcoxon on
(occ_post - occ_pre), LOO-max (feedback_no_single_patient_p_driven). A positive,
LOO-robust occupancy excess (with task-state concentration lift > 1) would be the
first evidence for discrete task-state reinstatement; otherwise the scoped
negative stands.

Outputs
-------
    data/audit/replay_states/state_decomposition_per_patient.csv
    data/audit/replay_states/state_decomposition_cohort.csv
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.utils.surrogate.matched_strength import cophenetic_condensed_from_adjacency
from lrg_eegfc.workflow.fc import load_fc_matrix

OUT = ROOT / "data" / "audit" / "replay_states"
CACHE = OUT / "cache"

COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]


def coph(W):
    return cophenetic_condensed_from_adjacency(np.clip(np.asarray(W, float), 0.0, 1.0))


def load_cached(pat, band, phase):
    f = CACHE / f"{pat}_{band}_{phase}_coph.npz"
    if not f.exists():
        return None
    return np.load(f)["coph"].astype(np.float64)


def count_runs(mask):
    runs = 0
    prev = False
    for v in mask:
        if v and not prev:
            runs += 1
        prev = v
    return runs


def zscore_pooled(x, ref):
    mu, sd = ref.mean(), ref.std()
    return (x - mu) / sd if sd > 0 else x * 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patients", nargs="+", default=COHORT)
    ap.add_argument("--bands", nargs="+", default=["alpha", "beta"])
    ap.add_argument("--Ks", nargs="+", type=int, default=[3, 4, 5])
    ap.add_argument("--n-pc", type=int, default=10)
    ap.add_argument("--seed", type=int, default=20260622)
    args = ap.parse_args()

    t0 = time.time()
    rows = []
    for pat in args.patients:
        # references (full-phase cophenetic centroids)
        refs = {}
        ok = True
        for b in args.bands:
            Wpre = load_fc_matrix(pat, "rest_pre", b, fc_method="imcoh_abs")
            WT = load_fc_matrix(pat, "task_test", b, fc_method="imcoh_abs")
            if Wpre is None or WT is None:
                ok = False; break
            refs[b] = (coph(Wpre), coph(WT))
        # cached per-window cophenetics
        A = {b: {ph: load_cached(pat, b, ph) for ph in ("rest_pre", "task_test", "rest_post")}
             for b in args.bands}
        if not ok or any(A[b][ph] is None for b in args.bands for ph in A[b]):
            print(f"[skip] {pat}: missing refs/cache"); continue

        # ---------- METHOD 1: nearest-centroid discriminant ----------
        g = {ph: [] for ph in ("rest_pre", "rest_post")}
        for ph in ("rest_pre", "rest_post"):
            per_band = []
            for b in args.bands:
                cpre, ctask = refs[b]
                gb = np.array([spearmanr(c, ctask)[0] - spearmanr(c, cpre)[0] for c in A[b][ph]])
                per_band.append(gb)
            g[ph] = per_band  # list over bands
        # z-score each band on pooled pre+post, then average bands
        gm = {}
        for ph in ("rest_pre", "rest_post"):
            zb = []
            for bi in range(len(args.bands)):
                pooled = np.concatenate([g["rest_pre"][bi], g["rest_post"][bi]])
                zb.append(zscore_pooled(g[ph][bi], pooled))
            gm[ph] = np.mean(zb, axis=0)
        pooled_all = np.concatenate([gm["rest_pre"], gm["rest_post"]])
        theta = np.median(pooled_all)
        occ_pre_nc = float(np.mean(gm["rest_pre"] > theta))
        occ_post_nc = float(np.mean(gm["rest_post"] > theta))
        runs_pre = count_runs(gm["rest_pre"] > theta)
        runs_post = count_runs(gm["rest_post"] > theta)

        # ---------- METHOD 2: k-means state discovery ----------
        # multiband feature: per-band column-standardize on pooled-all-phase, concat
        feats = {}
        for ph in ("rest_pre", "task_test", "rest_post"):
            cols = []
            for b in args.bands:
                allp = np.vstack([A[b]["rest_pre"], A[b]["task_test"], A[b]["rest_post"]])
                mu = allp.mean(0); sd = allp.std(0); sd[sd == 0] = 1.0
                cols.append((A[b][ph] - mu) / sd)
            feats[ph] = np.hstack(cols)
        Xall = np.vstack([feats["rest_pre"], feats["task_test"], feats["rest_post"]])
        lab = np.array(["rest_pre"] * len(feats["rest_pre"]) +
                       ["task_test"] * len(feats["task_test"]) +
                       ["rest_post"] * len(feats["rest_post"]))
        npc = min(args.n_pc, Xall.shape[0] - 1, Xall.shape[1])
        Z = PCA(n_components=npc, random_state=args.seed).fit_transform(Xall)

        km_res = {}
        for K in args.Ks:
            km = KMeans(n_clusters=K, n_init=10, random_state=args.seed).fit(Z)
            assign = km.labels_
            is_task = lab == "task_test"
            # task-state = cluster where task windows most concentrate (P(cluster|task))
            frac_task = [np.mean(assign[is_task] == k) for k in range(K)]
            ts = int(np.argmax(frac_task))
            # lift = P(task in ts) / (size_ts/total)
            size_ts = np.mean(assign == ts)
            lift = (frac_task[ts] / size_ts) if size_ts > 0 else np.nan
            occ_pre = float(np.mean(assign[lab == "rest_pre"] == ts))
            occ_post = float(np.mean(assign[lab == "rest_post"] == ts))
            km_res[K] = (occ_pre, occ_post, float(lift))

        row = dict(patient=pat,
                   nc_occ_pre=occ_pre_nc, nc_occ_post=occ_post_nc,
                   nc_docc=occ_post_nc - occ_pre_nc,
                   nc_runs_pre=runs_pre, nc_runs_post=runs_post)
        for K in args.Ks:
            op, opo, lf = km_res[K]
            row[f"km{K}_occ_pre"] = op
            row[f"km{K}_occ_post"] = opo
            row[f"km{K}_docc"] = opo - op
            row[f"km{K}_lift"] = lf
        rows.append(row)
        print(f"  {pat}: NC occ pre={occ_pre_nc:.2f} post={occ_post_nc:.2f} (Δ{occ_post_nc-occ_pre_nc:+.2f}) "
              f"runs {runs_pre}->{runs_post} | "
              + " ".join(f"K{K}:Δ{km_res[K][1]-km_res[K][0]:+.2f}(lift{km_res[K][2]:.1f})" for K in args.Ks))

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "state_decomposition_per_patient.csv", index=False)

    # ---------- cohort verdict ----------
    def cohort_test(col):
        d = df[col].to_numpy()
        d = d[np.isfinite(d)]
        if len(d) < 5:
            return dict(n=len(d), median=np.nan, n_pos=0, wilcoxon_p=np.nan, loo_max_p=np.nan)
        try:
            p = wilcoxon(d, alternative="greater")[1]
        except Exception:
            p = np.nan
        loo = []
        for i in range(len(d)):
            dd = np.delete(d, i)
            try:
                loo.append(wilcoxon(dd, alternative="greater")[1])
            except Exception:
                loo.append(np.nan)
        return dict(n=len(d), median=float(np.median(d)),
                    n_pos=int((d > 0).sum()), wilcoxon_p=float(p),
                    loo_max_p=float(np.nanmax(loo)) if loo else np.nan)

    coh = []
    coh.append(dict(method="nearest_centroid", **cohort_test("nc_docc")))
    for K in args.Ks:
        coh.append(dict(method=f"kmeans_K{K}", **cohort_test(f"km{K}_docc")))
    cohdf = pd.DataFrame(coh)
    cohdf.to_csv(OUT / "state_decomposition_cohort.csv", index=False)

    print("\n=== STATE-DECOMPOSITION COHORT VERDICT (replay = occ_post > occ_pre) ===")
    print(cohdf.to_string(index=False))
    # mean k-means lift (clustering quality: does task_test actually concentrate?)
    for K in args.Ks:
        lifts = df[f"km{K}_lift"].to_numpy()
        print(f"  K={K} task-state concentration lift: median={np.nanmedian(lifts):.2f} "
              f"(>1 means task windows concentrate; ~1 means clustering is noise)")
    print(f"\n[audit_126] done in {time.time()-t0:.1f}s -> {OUT}")


if __name__ == "__main__":
    main()
