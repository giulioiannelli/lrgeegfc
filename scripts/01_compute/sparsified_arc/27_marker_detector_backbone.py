#!/usr/bin/env python3
r"""Fused leave-one-patient-out SOZ detector on a chosen backbone -- headline §3 numbers.

Reproduces the detector paragraph / fig:epi1c / fig:epi2 numbers on the SETTLED marker
backbone (TMFG), so §3 can move off the old dense audit_117 detector. Backbone-selectable
(SA_BACKBONE, default tmfg) so mst@0.20 can be run for a side-by-side.

Per (patient, band), phase=rest_post: backbone -> eig -> seeded heat-kernel affinity at the
single operating scale s=10 (tau=10/lambda_max), strength-residualised (06 primitives). The
6-band residual markers are z-scored per patient and fused by a leave-one-patient-out logistic
(class-balanced) -- identical fusion to the design-space audit (auc_s10). Reports, per patient
and cohort:
  fused AUC (concordance), precision@5, calibrated P(SOZ) at seizure vs healthy contacts,
  label-shuffle-LOPO null AUC, delta-only deployable AUC, and the two-population split at
  AUC 0.60 (community vs right-hemisphere hub).

Reuses 06 primitives (_resid,_seeded_marker,_conc) + select_backbone. Writes
data/sparsified_arc/marker_detector_<backbone>/{detector_lopo_per_patient.csv,cohort.json}.
"""
from __future__ import annotations
import importlib.util, json, os, sys, time
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
_spec = importlib.util.spec_from_file_location(
    "epi_arc06", ROOT / "scripts/01_compute/sparsified_arc/06_epi_arc.py")
_epi = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_epi)
_resid, _seeded_marker, _conc = _epi._resid, _epi._seeded_marker, _epi._conc

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import load_phase, COHORT, BANDS, BASE_SEED   # noqa: E402
from lrg_eegfc.utils.io.patient import build_epi_masks                     # noqa: E402
from lrg_eegfc.utils.fc.backbone import select_backbone                    # noqa: E402
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, scale_grid   # noqa: E402

BACKBONE = os.environ.get("SA_BACKBONE", "tmfg")
FRAC = float(os.environ.get("SA_FRAC", "0.20"))
DISP_ALPHA = float(os.environ.get("SA_DISP_ALPHA", "0.20"))
OUT = ROOT / "data" / "sparsified_arc" / f"marker_detector_{BACKBONE}"
N_S = 16
S_OP = 10.0            # single operating scale (design-audit s=10)
MIN_EPI = 3
HUB_AUC = 0.60        # fused AUC below -> hub-type implant (two populations)


def build_marker(pat, band):
    """strength-residual seeded marker at s=10 + labels, on the chosen backbone."""
    W = load_phase(pat, "rest_post", band)
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None
    y = np.asarray(pm.epi_mask, bool)
    if int(y.sum()) < MIN_EPI or (~y).sum() < MIN_EPI:
        return None
    B = select_backbone(W, BACKBONE, frac=FRAC, disparity_alpha=DISP_ALPHA)
    ev, V = laplacian_eig(B)
    strength = B.sum(1)
    s_grid = scale_grid(ev, n=N_S)
    s = s_grid[int(np.argmin(np.abs(s_grid - S_OP)))]
    K = (V * np.exp(-(s / ev[-1]) * ev)) @ V.T
    m = _seeded_marker(K, np.where(y)[0])
    mr = _resid(np.nan_to_num(m, nan=np.nanmin(m)), strength)
    return mr, y


def zc(x):
    mu, sd = x.mean(0), x.std(0); sd[sd == 0] = 1.0
    return (x - mu) / sd


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    # gather per-patient 6-band feature matrices + labels
    feats, ys = {}, {}
    for p in COHORT:
        cols = []
        ok = True
        for b in BANDS:
            r = build_marker(p, b)
            if r is None:
                ok = False; break
            cols.append(r[0]); y = r[1]
        if ok:
            feats[p] = zc(np.column_stack(cols)); ys[p] = y
    pats = [p for p in COHORT if p in feats]
    print(f"[detector-{BACKBONE}] {len(pats)} patients x {len(BANDS)} bands, s={S_OP} "
          f"[{time.time()-t0:.0f}s]", flush=True)

    rng = np.random.default_rng(BASE_SEED)
    di = BANDS.index("delta")
    rows = []
    prob_pool, y_pool = [], []        # for the pooled calibration curve
    node_rows = []                    # per-node LOPO p_soz (for the ROC / precision@k panels)
    d1_prec5, gbm_auc, gbm_prec5 = [], [], []
    for held in pats:
        tr = [p for p in pats if p != held]
        Xtr = np.vstack([feats[p] for p in tr]); ytr = np.concatenate([ys[p] for p in tr])
        yh = ys[held]; ns = int(yh.sum())
        clf = LogisticRegression(max_iter=500, class_weight="balanced").fit(Xtr, ytr)
        sc = clf.decision_function(feats[held])
        pr = clf.predict_proba(feats[held])[:, 1]
        auc = _conc(sc[yh], sc[~yh])
        order = np.argsort(-sc)
        p5 = yh[order][:5].mean(); rprec = yh[order][:ns].mean()
        prob_pool.append(pr); y_pool.append(yh)
        for ni, (isz, p) in enumerate(zip(yh.astype(int), pr)):   # per-node p_soz for ROC/prec@k
            node_rows.append(dict(patient=held, node=int(ni), is_soz=int(isz), p_soz=float(p)))
        # delta-only deployable AUC + precision@5
        Xtr_d = np.vstack([feats[p][:, [di]] for p in tr])
        clf_d = LogisticRegression(max_iter=500, class_weight="balanced").fit(Xtr_d, ytr)
        sc_d = clf_d.decision_function(feats[held][:, [di]])
        auc_d = _conc(sc_d[yh], sc_d[~yh]); d1_prec5.append(yh[np.argsort(-sc_d)][:5].mean())
        # nonlinear (GBM) 6-band: higher AUC / no precision gain -> rejected
        gbm = GradientBoostingClassifier(random_state=0).fit(Xtr, ytr)
        sg = gbm.predict_proba(feats[held])[:, 1]
        gbm_auc.append(_conc(sg[yh], sg[~yh])); gbm_prec5.append(yh[np.argsort(-sg)][:5].mean())
        # label-shuffle LOPO null
        ysh = rng.permutation(yh)
        auc_sh = _conc(sc[ysh], sc[~ysh])
        rows.append(dict(patient=held, n_soz=ns, N=len(yh), auc=float(auc),
                         prec5=float(p5), rprec=float(rprec), auc_delta_only=float(auc_d),
                         p_mean_soz=float(pr[yh].mean()), p_mean_healthy=float(pr[~yh].mean()),
                         auc_shuffle=float(auc_sh),
                         population="community" if auc >= HUB_AUC else "hub"))
    det = pd.DataFrame(rows).sort_values("auc", ascending=False)
    det.to_csv(OUT / "detector_lopo_per_patient.csv", index=False)
    pd.DataFrame(node_rows).to_csv(OUT / "detector_node_predictions.csv", index=False)

    # pooled reliability curve (equal-count bins over predicted P)
    P = np.concatenate(prob_pool); Y = np.concatenate(y_pool).astype(float)
    nb = 8
    q = np.quantile(P, np.linspace(0, 1, nb + 1)); q[-1] += 1e-9
    binid = np.clip(np.digitize(P, q[1:-1]), 0, nb - 1)
    cal = pd.DataFrame([dict(p_pred_mean=float(P[binid == b].mean()),
                             soz_observed=float(Y[binid == b].mean()),
                             n=int((binid == b).sum()))
                        for b in range(nb) if (binid == b).any()])
    cal.to_csv(OUT / "calibration_curve.csv", index=False)

    comm = det[det.population == "community"]; hub = det[det.population == "hub"]
    ablation = dict(delta_only=dict(auc=round(float(det.auc_delta_only.median()), 2),
                                    prec5=round(float(np.median(d1_prec5)), 2)),
                    six_band=dict(auc=round(float(det.auc.median()), 2),
                                  prec5=round(float(det.prec5.median()), 2)),
                    gbm=dict(auc=round(float(np.median(gbm_auc)), 2),
                             prec5=round(float(np.median(gbm_prec5)), 2)))
    summ = dict(
        backbone=BACKBONE, n_pat=len(det), s_op=S_OP,
        fused_auc_med=float(det.auc.median()), prec5_med=float(det.prec5.median()),
        rprec_med=float(det.rprec.median()), ablation=ablation,
        p_mean_soz=float(det.p_mean_soz.median()), p_mean_healthy=float(det.p_mean_healthy.median()),
        auc_shuffle_med=float(det.auc_shuffle.median()),
        n_above_chance=int((det.auc > 0.5).sum()),
        n_community=len(comm), n_hub=len(hub),
        auc_max=float(det.auc.max()), auc_ge_090=int((det.auc >= 0.90).sum()),
        rescued=det[det.patient.isin(["Pat_10", "Pat_15"])].set_index("patient").auc.round(3).to_dict(),
        hub_patients=hub.patient.tolist(),
    )
    (OUT / "cohort.json").write_text(json.dumps(summ, indent=2))
    print(f"\n=== FUSED DETECTOR ({BACKBONE}, LOPO, s={S_OP}) ===", flush=True)
    print(det.round(3).to_string(index=False), flush=True)
    print("\n--- cohort ---", flush=True)
    for k, v in summ.items():
        print(f"  {k}: {v}", flush=True)
    print(f"\n[detector-{BACKBONE}] done in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
