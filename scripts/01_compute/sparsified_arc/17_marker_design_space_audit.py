#!/usr/bin/env python3
r"""Design-space audit of the SOZ propagator marker — are we at the best setting? (coord descent)

The marker's pipeline inherited several choices from the TRACE side (backbone=mst@0.20,
phase=rest_post) that were never justified for SOZ DETECTION. This audits each knob by
coordinate descent from the current setting, one axis at a time, reporting the WHOLE profile
(not the max), so we can see whether the current setting is near-optimal or leaving
performance on the table. Axes:
  1. backbone   : dense | mst@{0.05,0.10,0.20,0.35,0.50} | percolation | tmfg
  2. phase      : rest_pre | task_learn | task_test | rest_post
  3. marker     : mean | max | sum affinity to seed (leave-self-out), strength-residual
  4. band-subset: all6 | delta+beta+low_gamma | +high_gamma

Metric per config (fused 6-band, nested LOPO logistic over the scale grid):
  auc_s10  = fused AUC at the fixed s=10 operating point   (fully honest, no scale pick)
  auc_pk   = fused AUC at the in-sample peak scale          (ceiling; labelled ORACLE)
  rprec_pk = R-precision (precision at k=n_soz; 'identify-all-SOZ') at that peak
  d_auc    = single-band delta multiscale-mean AUC          (reference)
Forking-path guard: cross-config ranking uses auc_s10 (fixed); auc_pk/rprec_pk are the
in-sample ceiling per config, shown for shape only. The winning config is re-checked with
the nested best-global-scale selection (module 16) before any claim.

Reuses 06/16 primitives + lrg_eegfc.utils.fc.backbone. Writes ms_marker_exploration/design_space_*.csv.
"""
from __future__ import annotations
import importlib.util
import sys
import time

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

_spec = importlib.util.spec_from_file_location(
    "epi_arc06", ROOT / "scripts/01_compute/sparsified_arc/06_epi_arc.py")
_epi = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_epi)
_resid, _seeded_marker, _conc = _epi._resid, _epi._seeded_marker, _epi._conc

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import load_phase, COHORT, BANDS   # noqa: E402
from lrg_eegfc.utils.io.patient import build_epi_masks          # noqa: E402
from lrg_eegfc.utils.fc import backbone as bb                   # noqa: E402
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, scale_grid  # noqa: E402

N_S = 16
MIN_EPI = 3
IDX_S10 = 6
OUT = ROOT / "data/sparsified_arc/ms_marker_exploration"
BAND_SUBSETS = {
    "all6": list(BANDS),
    "d_b_gl": ["delta", "beta", "low_gamma"],
    "d_b_gl_gh": ["delta", "beta", "low_gamma", "high_gamma"],
}


def make_backbone(W, kind):
    if kind == "dense":
        A = bb._clean(W); np.fill_diagonal(A, 0.0); return A
    if kind.startswith("mst"):
        return bb.mst_union_top_fraction(W, float(kind[3:]))
    if kind == "perc":
        return bb.percolation_backbone(W)[0]
    if kind == "tmfg":
        return bb.tmfg_backbone(W)
    raise ValueError(kind)


def _agg(K, seed, kind):
    Ks = K[:, seed]
    self_in = np.isin(np.arange(K.shape[0]), seed)
    if kind == "mean":
        return _seeded_marker(K, seed)
    if kind == "max":
        Km = Ks.copy()
        for a, j in enumerate(seed):        # blank the self column so max is leave-self-out
            Km[j, a] = -np.inf
        m = Km.max(1); m[np.isneginf(m)] = np.nanmin(m[np.isfinite(m)]); return m
    if kind == "sum":
        m = Ks.sum(1) - np.where(self_in, K[np.arange(K.shape[0]), np.arange(K.shape[0])], 0.0)
        return m
    raise ValueError(kind)


def build_stacks(backbone="mst0.20", phase="rest_post", marker="mean", bands=None):
    """{pat:{band:(16,N)}}, {pat:y}, s-grid — for one config."""
    bands = bands or list(BANDS)
    stacks, ys = {}, {}
    for p in COHORT:
        pm = build_epi_masks(p)
        for b in bands:
            try:
                W = load_phase(p, phase, b)
            except Exception:
                continue
            N = W.shape[0]
            if len(pm.channels) != N:
                continue
            y = np.asarray(pm.epi_mask, bool)
            if int(y.sum()) < MIN_EPI or (~y).sum() < MIN_EPI:
                continue
            B = make_backbone(W, backbone)
            ev, V = laplacian_eig(B)
            strength = B.sum(1); s_grid = scale_grid(ev, n=N_S); seed = np.where(y)[0]
            st = np.empty((N_S, N))
            for i, s in enumerate(s_grid):
                K = (V * np.exp(-(s / ev[-1]) * ev)) @ V.T
                m = _agg(K, seed, marker)
                st[i] = _resid(np.nan_to_num(m, nan=np.nanmin(m)), strength)
            stacks.setdefault(p, {})[b] = st; ys[p] = y
    return stacks, ys


def _zscore(X):
    mu = X.mean(0); sd = X.std(0); sd[sd == 0] = 1.0
    return (X - mu) / sd


def fused_profile(stacks, ys, bands):
    """LOPO fused logistic at each scale -> auc[16], rprec[16] (cohort mean)."""
    pats = [p for p in COHORT if p in stacks and all(b in stacks[p] for b in bands)]
    auc = np.full(N_S, np.nan); rpr = np.full(N_S, np.nan)
    for si in range(N_S):
        X = {p: _zscore(np.column_stack([stacks[p][b][si] for b in bands])) for p in pats}
        a = np.empty(len(pats)); r = np.empty(len(pats))
        for i, held in enumerate(pats):
            tr = [p for p in pats if p != held]
            Xtr = np.vstack([X[p] for p in tr]); ytr = np.concatenate([ys[p] for p in tr])
            clf = LogisticRegression(max_iter=400, class_weight="balanced").fit(Xtr, ytr)
            sc = clf.decision_function(X[held]); yh = ys[held]
            a[i] = _conc(sc[yh], sc[~yh])
            order = np.argsort(-sc); ns = int(yh.sum())
            r[i] = yh[order][:ns].mean()
        auc[si] = a.mean(); rpr[si] = r.mean()
    return auc, rpr, len(pats)


def single_band_delta(stacks, ys):
    pats = [p for p in COHORT if p in stacks and "delta" in stacks[p]]
    vals = []
    for p in pats:
        st = stacks[p]["delta"]; y = ys[p]
        vals.append(np.nanmean([_conc(st[i][y], st[i][~y]) for i in range(N_S)]))
    return float(np.nanmean(vals))


def eval_config(label, backbone, phase, marker, bands, rows):
    t0 = time.time()
    stacks, ys = build_stacks(backbone, phase, marker, bands)
    auc, rpr, npat = fused_profile(stacks, ys, bands)
    d_auc = single_band_delta(stacks, ys)
    pk = int(np.nanargmax(auc))
    rec = dict(config=label, backbone=backbone, phase=phase, marker=marker,
               bands="+".join(b[0] for b in bands), npat=npat,
               auc_s10=round(float(auc[IDX_S10]), 3), auc_pk=round(float(auc[pk]), 3),
               rprec_pk=round(float(rpr[pk]), 3), s_pk_idx=pk, d_auc=round(d_auc, 3))
    rows.append(rec)
    print(f"  {label:22s} auc_s10={rec['auc_s10']:.3f}  auc_pk={rec['auc_pk']:.3f}"
          f"(#{pk})  Rprec_pk={rec['rprec_pk']:.3f}  δ={rec['d_auc']:.3f}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    return auc, rpr


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    cur = dict(backbone="mst0.20", phase="rest_post", marker="mean", bands=list(BANDS))

    print("=== AXIS 1: backbone (phase=rest_post, marker=mean, all6) ===", flush=True)
    for bkb in ["dense", "mst0.05", "mst0.10", "mst0.20", "mst0.35", "mst0.50", "perc", "tmfg"]:
        eval_config(f"bkb={bkb}", bkb, cur["phase"], cur["marker"], cur["bands"], rows)
    best_bkb = max([r for r in rows if r["config"].startswith("bkb=")],
                   key=lambda r: r["auc_s10"])["backbone"]
    cur["backbone"] = best_bkb
    print(f"  -> best backbone by auc_s10: {best_bkb}", flush=True)

    print(f"\n=== AXIS 2: phase (backbone={best_bkb}, marker=mean, all6) ===", flush=True)
    for ph in ["rest_pre", "task_learn", "task_test", "rest_post"]:
        eval_config(f"phase={ph}", best_bkb, ph, cur["marker"], cur["bands"], rows)
    best_ph = max([r for r in rows if r["config"].startswith("phase=")],
                  key=lambda r: r["auc_s10"])["phase"]
    cur["phase"] = best_ph
    print(f"  -> best phase by auc_s10: {best_ph}", flush=True)

    print(f"\n=== AXIS 3: marker aggregation (backbone={best_bkb}, phase={best_ph}, all6) ===",
          flush=True)
    for mk in ["mean", "max", "sum"]:
        eval_config(f"marker={mk}", best_bkb, best_ph, mk, cur["bands"], rows)
    best_mk = max([r for r in rows if r["config"].startswith("marker=")],
                  key=lambda r: r["auc_s10"])["marker"]
    cur["marker"] = best_mk
    print(f"  -> best marker by auc_s10: {best_mk}", flush=True)

    print(f"\n=== AXIS 4: band subset (backbone={best_bkb}, phase={best_ph}, marker={best_mk}) ===",
          flush=True)
    for sname, subset in BAND_SUBSETS.items():
        eval_config(f"bands={sname}", best_bkb, best_ph, best_mk, subset, rows)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "design_space_audit.csv", index=False)
    print("\n=== FULL AUDIT TABLE ===", flush=True)
    print(df.to_string(index=False), flush=True)
    base = next(r for r in rows if r["config"] == "bkb=mst0.20")
    best = max(rows, key=lambda r: r["auc_s10"])
    print(f"\n  current setting (mst0.20/rest_post/mean/all6): auc_s10={base['auc_s10']:.3f}", flush=True)
    print(f"  best honest config: {best['config']} "
          f"[{best['backbone']}/{best['phase']}/{best['marker']}/{best['bands']}] "
          f"auc_s10={best['auc_s10']:.3f}  Rprec_pk={best['rprec_pk']:.3f}", flush=True)


if __name__ == "__main__":
    main()
