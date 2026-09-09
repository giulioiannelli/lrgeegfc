#!/usr/bin/env python3
"""FEASIBILITY PROBE (2026-09-09) -- NOT A TEST. No permutation null, no LOPO fusion, no CI.
Regenerates the SOZ-marker numbers quoted in
.agents/plans/active/2026-09-08_new-paths-to-objectives.md (objective 4).

What it asks: on the rest_pre / rest_post |ImCoh| graphs (canonical backbone), do heat-kernel node
scores separate clinically marked SOZ contacts from the rest better than (a) strength/degree/
centrality baselines and (b) an electrode-density spatial score (Conrad 2022 confound)?
Scores per node: strength (dense, backbone), degree, eigenvector centrality, PageRank,
heat-kernel K_ii and same-shaft-excluded off-diagonal row sum at s in S_LIST, each also
strength-residualised (degree-1 polyfit, audit_101 convention), plus the July seeded affinity
(mean K_ij over the SOZ seed set, self excluded) which USES the labels and is reported only as
a transductive reference, never as a detector.
AUC is P[score_soz > score_non] with 0.5 for ties, per patient; cohort median and count > 0.7;
pooled AUC on within-patient z-scores is labelled pooled. Patients with < 3 SOZ contacts or
all-SOZ are skipped and listed.
Run from the worktree: PYTHONPATH=src <lapbrain python> -u <this file>.
"""
import time, numpy as np, pandas as pd
from scipy.stats import spearmanr
from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig
from lrg_eegfc.utils.io.patient import build_epi_masks
from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.probe import build_probe_mask
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.phase_graphs import load_phase_fc
from lrg_eegfc.workflow.substrate import CANONICAL
ROOT = setup_script_env()
OUT = ROOT / "data" / "paper_final" / "feasibility"; OUT.mkdir(parents=True, exist_ok=True)
S_LIST = [0.5, 1.0, 2.0, 4.72, 10.0, 15.9, 50.0, 180.0]
MIN_EPI = 3

def auc(case, ctrl):
    if case.size == 0 or ctrl.size == 0: return np.nan
    gt = (case[:, None] > ctrl[None, :]).sum(); eq = (case[:, None] == ctrl[None, :]).sum()
    return (gt + 0.5 * eq) / (case.size * ctrl.size)

def resid(score, strength):
    b = np.polyfit(strength, score, 1); return score - np.polyval(b, strength)

def heat(ev, V, s):
    w = np.exp(-(s / ev[-1]) * ev); K = (V * w) @ V.T; return K / np.trace(K)

def pagerank(B, d=0.85, it=200):
    N = B.shape[0]; M = B / np.maximum(B.sum(0, keepdims=True), 1e-12); p = np.ones(N) / N
    for _ in range(it): p = (1 - d) / N + d * M @ p
    return p

def density(xyz, same):
    """Electrode-density scores from contact coordinates (mm). NaN if no coordinates."""
    if xyz is None: return {}
    D = np.sqrt(((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(-1)); np.fill_diagonal(D, np.inf)
    out = {"dens10": (D < 10).sum(1).astype(float), "dens20": (D < 20).sum(1).astype(float),
           "nn5_negdist": -np.sort(D, 1)[:, :5].mean(1)}
    Dx = D.copy(); Dx[same] = np.inf
    out["dens20_xshaft"] = (Dx < 20).sum(1).astype(float)
    return out

rows, skipped, t0 = [], [], time.time()
for ip, pat in enumerate(PATIENTS_4PHASE, 1):
    pm = build_epi_masks(pat); y = pm.epi_mask; N = y.size
    if y.sum() < MIN_EPI or y.all():
        skipped.append((pat, int(y.sum()), N)); continue
    same = build_probe_mask(list(pm.channels))
    xyz = None
    try:
        reg = load_channel_regions(pat)
        if len(reg) == N: xyz = reg[["x", "y", "z"]].to_numpy(float) / 1000.0
    except Exception as e:
        print(f"  {pat}: no coordinates ({e})", flush=True)
    dens = density(xyz, same)
    for band in BRAIN_BANDS_NAMES:
        for ph in ("rest_pre", "rest_post"):
            W = load_phase_fc(pat, ph, band)
            if W.shape[0] != N: print(f"skip {pat}/{band}/{ph}: N mismatch", flush=True); continue
            B = select_backbone(W, CANONICAL.backbone, frac=CANONICAL.frac); ev, V = laplacian_eig(B)
            st = W.sum(1); stb = B.sum(1)
            sc = {"strength_dense": st, "strength_bb": stb, "degree_bb": (B > 0).sum(1).astype(float),
                  "evc_bb": np.abs(V[:, -1]) if False else np.abs(np.linalg.eigh(B)[1][:, -1]),
                  "pagerank_bb": pagerank(B)}
            seed = np.where(y)[0]
            for s in S_LIST:
                K = heat(ev, V, s); Ko = K.copy(); np.fill_diagonal(Ko, 0.0); Kx = Ko.copy(); Kx[same] = 0.0
                sc[f"Kii_s{s:g}"] = np.diag(K); sc[f"Krow_xshaft_s{s:g}"] = Kx.sum(1)
                m = K[:, seed].sum(1) - np.where(y, np.diag(K), 0.0); sc[f"seeded_s{s:g}"] = m / (len(seed) - y)
            sc.update(dens)
            names = list(sc)
            for nm in names:
                v = sc[nm]
                variants = {"raw": v}
                if not nm.startswith(("strength", "dens", "nn5")): variants["stresid"] = resid(v, stb)
                for var, vv in variants.items():
                    rows.append(dict(patient=pat, band=band, phase=ph, score=nm, variant=var, N=N, n_soz=int(y.sum()),
                                     auc=auc(vv[y], vv[~y]), rho_strength=spearmanr(vv, stb)[0],
                                     rho_dens20=spearmanr(vv, dens["dens20"])[0] if dens else np.nan,
                                     z=(vv - vv.mean()) / (vv.std() + 1e-12), y=y))
    print(f"[{ip}/{len(PATIENTS_4PHASE)}] {pat} N={N} soz={int(y.sum())} coords={'yes' if xyz is not None else 'no'} {time.time()-t0:.0f}s", flush=True)

df = pd.DataFrame(rows)
per = df.drop(columns=["z", "y"]); per.to_csv(OUT / "soz_marker_probe.csv", index=False)
# cohort summary: median AUC, count > 0.7, pooled z-AUC (labelled pooled)
summ = []
for (band, ph, sc, var), g in df.groupby(["band", "phase", "score", "variant"]):
    a = g["auc"].to_numpy(float); zz = np.concatenate(g["z"].to_list()); yy = np.concatenate(g["y"].to_list())
    summ.append(dict(band=band, phase=ph, score=sc, variant=var, n_pat=len(g), median_auc=np.nanmedian(a),
                     n_gt07=int((a > 0.7).sum()), n_lt03=int((a < 0.3).sum()), pooled_z_auc=auc(zz[yy], zz[~yy]),
                     median_rho_strength=g["rho_strength"].median(), median_rho_dens20=g["rho_dens20"].median()))
S = pd.DataFrame(summ); S.to_csv(OUT / "soz_marker_probe_summary.csv", index=False)
print("skipped:", skipped)
pd.set_option("display.width", 200, "display.max_rows", 400)
top = S[S.phase == "rest_post"].sort_values("median_auc", ascending=False)
print("\n== rest_post: top 25 by cohort median AUC ==\n", top.head(25).to_string(index=False, float_format=lambda x: f"{x:.3f}"))
base = S[(S.score.str.startswith(("strength", "degree", "evc", "pagerank", "dens", "nn5"))) & (S.variant == "raw")]
print("\n== baselines (all bands, both phases; density is band-independent) ==\n",
      base.groupby(["score", "phase"])[["median_auc", "n_gt07", "pooled_z_auc"]].median().to_string(float_format=lambda x: f"{x:.3f}"))
best = S[(S.score.str.startswith("Krow_xshaft")) & (S.variant == "stresid")].pivot_table(index="score", columns=["phase", "band"], values="median_auc")
print("\n== Krow_xshaft strength-residualised: median AUC by s (rows) x band ==\n", best.to_string(float_format=lambda x: f"{x:.2f}"))
print(f"\nwall {time.time()-t0:.0f}s -> {OUT}")
