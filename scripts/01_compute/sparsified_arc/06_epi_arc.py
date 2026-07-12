#!/usr/bin/env python3
"""D4 -- tau-resolved epileptic-node (SOZ) marker on the percolation backbone.

Ports the seeded propagator SOZ marker (audit_101/132) onto the sparsified
percolation backbone with a full diffusion-scale sweep, and asks: does the
MULTISCALE marker beat the single-scale (heat_t5 = 10/lambda_max) baseline, and
does precision improve? PRECISION is foregrounded (precision@k, lift, PR-AUC),
not just ROC-AUC.

Per (patient, band), phase=rest_post:
  - percolation backbone B, eigendecompose;
  - at each diffusion scale s: K(tau)=V exp(-tau Lambda) V^T, tau=s/lambda_max;
    seeded marker  m_i = mean_{j in SOZ, j!=i} K_ij  (leave-self-out),
    strength-residualised (removes the trivial hub effect);
  - metrics vs SOZ labels y: ROC-AUC (concordance), precision@{5,10}, recall@10,
    lift@5, PR-AUC (sklearn average_precision_score).
  - single-scale baseline = s=10 (heat_t5 analog);
    multiscale marker = mean over the scale grid of the residualised markers.
  - strength-matched fake-SOZ null (R sets, quantile-matched) on the multiscale
    marker -> per-cell p. Occult candidates = top-ranked UNMARKED nodes (seeded).

Beat targets: audit_132 heat tau5 AUC_resid delta.80/gamma_low.74/beta.69;
audit_117 detector AUC.81 precision@5 60%; literature SOZ AUC 0.70-0.86.

Outputs:
    data/sparsified_arc/epi_arc/{band}/{patient}.npz   (s-curves of AUC/prec/prauc)
    data/sparsified_arc/epi_arc/per_cell.csv
    data/sparsified_arc/epi_arc/cohort_summary.csv
    data/sparsified_arc/epi_arc/occult_candidates.csv
    data/sparsified_arc/epi_arc/config.json
"""
from __future__ import annotations
import json, os, sys, time
import numpy as np, pandas as pd
from sklearn.metrics import average_precision_score

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_150_rho_sym_gate import load_phase, COHORT, BANDS, BASE_SEED
from lrg_eegfc.utils.io.patient import build_epi_masks
from lrg_eegfc.utils.fc.backbone import percolation_backbone, mst_union_top_fraction
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, scale_grid

# Backbone via env: BACKBONE=mst020 (recovered scheme) | perc (legacy percolation).
BACKBONE = os.environ.get("SA_BACKBONE", "mst020")
FRAC = 0.20
OUT = ROOT / "data" / "sparsified_arc" / ("epi_arc_mst020" if BACKBONE == "mst020" else "epi_arc")
N_S = 16
R_NULL = 200
MIN_EPI = 3
S_SINGLE = 10.0            # heat_t5 analog: tau = 10/lambda_max


def _resid(score, strength):
    """Linear strength-residual (degree-1 polyfit), audit_101 convention."""
    b = np.polyfit(strength, score, 1)
    return score - np.polyval(b, strength)


def _conc(case, ctrl):
    """Concordance = ROC-AUC numerator (P[case>ctrl] + 0.5 ties)."""
    if case.size == 0 or ctrl.size == 0:
        return np.nan
    gt = (case[:, None] > ctrl[None, :]).sum()
    eq = (case[:, None] == ctrl[None, :]).sum()
    return (gt + 0.5 * eq) / (case.size * ctrl.size)


def _seeded_marker(K, seed_idx):
    """m_i = mean_{j in seed, j!=i} K_ij (leave-self-out)."""
    m = K[:, seed_idx].sum(1)
    self_in = np.isin(np.arange(K.shape[0]), seed_idx)
    denom = np.full(K.shape[0], len(seed_idx), float)
    m = m - np.where(self_in, K[np.arange(K.shape[0]), np.arange(K.shape[0])], 0.0)  # K_ii if self seeded
    denom = denom - self_in.astype(float)
    denom[denom <= 0] = np.nan
    return m / denom


def _metrics(marker_resid, y):
    soz, hea = marker_resid[y], marker_resid[~y]
    auc = _conc(soz, hea)
    order = np.argsort(-marker_resid)
    y_sorted = y[order]
    n_soz = int(y.sum())
    prev = y.mean()
    p5 = y_sorted[:5].mean(); p10 = y_sorted[:10].mean()
    r10 = y_sorted[:10].sum() / n_soz
    prauc = average_precision_score(y.astype(int), marker_resid) if 0 < n_soz < y.size else np.nan
    return dict(auc=auc, prec5=float(p5), prec10=float(p10), recall10=float(r10),
                lift5=float(p5 / prev) if prev > 0 else np.nan, prauc=float(prauc))


def _strength_sets(strength, y, R, rng, nbins=4):
    order = np.argsort(strength)
    bins = np.array_split(order, nbins)
    binof = np.empty(len(strength), int)
    for bi, b in enumerate(bins):
        binof[b] = bi
    counts = np.bincount(binof[np.where(y)[0]], minlength=nbins)
    sets = []
    for _ in range(R):
        fake = []
        for bi in range(nbins):
            k = counts[bi]
            if k > 0:
                fake.extend(rng.choice(bins[bi], size=min(k, len(bins[bi])), replace=False))
        sets.append(np.array(fake, int))
    return sets


def run_cell(pat, band):
    W = load_phase(pat, "rest_post", band)
    N = W.shape[0]
    pm = build_epi_masks(pat)
    if len(pm.channels) != N:
        return None, None
    y = np.asarray(pm.epi_mask, bool)
    if int(y.sum()) < MIN_EPI or (~y).sum() < MIN_EPI:
        return None, None
    B = mst_union_top_fraction(W, FRAC) if BACKBONE == "mst020" else percolation_backbone(W)[0]
    ev, V = laplacian_eig(B)
    strength = B.sum(1)
    s_grid = scale_grid(ev, n=N_S)
    seed = np.where(y)[0]

    aucs = np.empty(N_S); p5 = np.empty(N_S); p10 = np.empty(N_S); prauc = np.empty(N_S)
    marker_stack = np.empty((N_S, N))
    for i, s in enumerate(s_grid):
        K = (V * np.exp(-(s / ev[-1]) * ev)) @ V.T
        m = _seeded_marker(K, seed)
        mr = _resid(np.nan_to_num(m, nan=np.nanmin(m)), strength)
        marker_stack[i] = mr
        mm = _metrics(mr, y)
        aucs[i] = mm["auc"]; p5[i] = mm["prec5"]; p10[i] = mm["prec10"]; prauc[i] = mm["prauc"]

    # single-scale baseline (s~10) and multiscale-mean marker
    j_single = int(np.argmin(np.abs(s_grid - S_SINGLE)))
    single = _metrics(marker_stack[j_single], y)
    multi_marker = marker_stack.mean(0)
    multi = _metrics(_resid(multi_marker, strength), y)
    j_best = int(np.nanargmax(aucs))

    # strength-matched fake-SOZ null on the multiscale marker
    rng = np.random.default_rng(BASE_SEED + hash((pat, band)) % (2**31))
    Kmean = np.zeros((N, N))
    for s in s_grid:
        Kmean += (V * np.exp(-(s / ev[-1]) * ev)) @ V.T
    Kmean /= N_S
    real_auc = multi["auc"]
    null_auc = []
    for fs in _strength_sets(strength, y, R_NULL, rng):
        if fs.size < MIN_EPI:
            continue
        yf = np.zeros(N, bool); yf[fs] = True
        mf = _resid(_seeded_marker(Kmean, fs), strength)
        null_auc.append(_conc(mf[yf], mf[~yf]))
    null_auc = np.array([a for a in null_auc if np.isfinite(a)])
    p_null = float(np.mean(null_auc >= real_auc)) if null_auc.size else np.nan

    cell = OUT / band
    cell.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cell / f"{pat}.npz", s=s_grid, auc=aucs, prec5=p5, prec10=p10,
                        prauc=prauc, n_soz=int(y.sum()), N=N)

    rec = dict(patient=pat, band=band, N_nodes=int(N), n_soz=int(y.sum()),
               prevalence=float(y.mean()),
               auc_single=single["auc"], prec5_single=single["prec5"],
               prec10_single=single["prec10"], prauc_single=single["prauc"],
               auc_multi=multi["auc"], prec5_multi=multi["prec5"],
               prec10_multi=multi["prec10"], recall10_multi=multi["recall10"],
               lift5_multi=multi["lift5"], prauc_multi=multi["prauc"],
               auc_best_s=float(aucs[j_best]), best_s=float(s_grid[j_best]),
               p_null_multi=p_null,
               auc_gain_multi_vs_single=float(multi["auc"] - single["auc"]))

    # occult candidates: top UNMARKED nodes by the multiscale marker
    mr = _resid(multi_marker, strength)
    hea_idx = np.where(~y)[0]
    top = hea_idx[np.argsort(-mr[hea_idx])[:5]]
    occ = [dict(patient=pat, band=band, rank=k + 1, node=int(n),
                channel=str(pm.channels[n]), probe=str(pm.probes[n]),
                marker_z=float((mr[n] - mr.mean()) / (mr.std() + 1e-12)))
           for k, n in enumerate(top)]
    return rec, occ


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    jobs = [(p, b) for b in BANDS for p in COHORT]
    print(f"[D4-epi] {len(jobs)} cells, R_null={R_NULL}", flush=True)
    t0 = time.time(); rows = []; occs = []
    for i, (p, b) in enumerate(jobs, 1):
        try:
            rec, occ = run_cell(p, b)
        except Exception as e:
            print(f"[{i}/{len(jobs)}] SKIP {p}/{b}: {e}", flush=True); continue
        if rec is None:
            continue
        rows.append(rec); occs.extend(occ)
        if i % 10 == 0 or i == len(jobs):
            print(f"[{i}/{len(jobs)}] {p}/{b} AUC single={rec['auc_single']:.2f} "
                  f"multi={rec['auc_multi']:.2f}(best {rec['auc_best_s']:.2f}@s{rec['best_s']:.0f}) "
                  f"p{rec['p_null_multi']:.3f} prec@5={rec['prec5_multi']:.2f} "
                  f"[{time.time()-t0:.0f}s]", flush=True)
    df = pd.DataFrame(rows).sort_values(["band", "patient"])
    df.to_csv(OUT / "per_cell.csv", index=False)
    pd.DataFrame(occs).to_csv(OUT / "occult_candidates.csv", index=False)

    summ = []
    for b in BANDS:
        x = df[df.band == b]
        if x.empty:
            continue
        summ.append(dict(band=b, n_pat=len(x),
                         auc_single_med=float(x.auc_single.median()),
                         auc_multi_med=float(x.auc_multi.median()),
                         auc_best_med=float(x.auc_best_s.median()),
                         prec5_single_med=float(x.prec5_single.median()),
                         prec5_multi_med=float(x.prec5_multi.median()),
                         prauc_multi_med=float(x.prauc_multi.median()),
                         lift5_multi_med=float(x.lift5_multi.median()),
                         n_beat_null=int((x.p_null_multi < 0.05).sum()),
                         multi_gt_single=int((x.auc_gain_multi_vs_single > 0).sum())))
    sdf = pd.DataFrame(summ); sdf.to_csv(OUT / "cohort_summary.csv", index=False)
    print("\n=== EPI MARKER: multiscale vs single-scale (percolation backbone, rest_post) ===", flush=True)
    print(sdf.to_string(index=False), flush=True)
    (OUT / "config.json").write_text(json.dumps(dict(
        deliverable="D4_epi_arc", backbone=BACKBONE, phase="rest_post",
        marker="seeded heat-kernel, leave-self-out, strength-residualised", n_scales=N_S,
        single_scale="s=10 (heat_t5, tau=10/lambda_max)", multiscale="mean residual marker over grid",
        null="strength-matched fake-SOZ (R=200, 4 quantile bins)", cohort=COHORT, bands=BANDS,
        beat_targets="audit_132 AUC delta.80/gl.74/beta.69; audit_117 AUC.81 prec@5 .60; SOTA .70-.86",
        precision_note="foreground precision@k/lift/PR-AUC, not just ROC-AUC",
    ), indent=2))
    print(f"[D4-epi] {len(df)} cells in {time.time()-t0:.0f}s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
