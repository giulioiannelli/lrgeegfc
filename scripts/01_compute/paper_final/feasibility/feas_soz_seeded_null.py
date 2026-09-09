#!/usr/bin/env python3
"""FEASIBILITY PROBE (2026-09-09) -- NOT A TEST. Companion of feas_soz_marker.py.
The seeded affinity (mean heat-kernel K_ij from node i to the SOZ seed set, self excluded) is the
only score that reaches AUC ~0.85, and it USES the labels. Two controls decide what that number means:
 (1) same-shaft exclusion: K_ij with same-shaft pairs zeroed before seeding (is it just shaft
     adjacency of SOZ contacts?);
 (2) random-seed null: draw seed sets of the same size (unmatched, and matched on the per-shaft
     count of SOZ contacts) and compute the AUC with which affinity-to-set recovers the set itself.
     A coherent random set also scores > 0.5 by construction; the SOZ claim is the excess over that.
Per patient x band, rest_post, canonical backbone, s in S_LIST. Reports observed AUC, null median,
null 95th percentile, per-patient z, and cohort counts. 200 draws per null.
"""
import time, numpy as np, pandas as pd
from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig
from lrg_eegfc.utils.io.patient import build_epi_masks
from lrg_eegfc.utils.probe import build_probe_mask
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.phase_graphs import load_phase_fc
from lrg_eegfc.workflow.substrate import CANONICAL
ROOT = setup_script_env(); OUT = ROOT / "data" / "paper_final" / "feasibility"
S_LIST = [1.0, 4.72, 10.0, 50.0]; NDRAW = 200; rng = np.random.default_rng(0)

def auc(case, ctrl):
    gt = (case[:, None] > ctrl[None, :]).sum(); eq = (case[:, None] == ctrl[None, :]).sum()
    return (gt + 0.5 * eq) / (case.size * ctrl.size)

def seeded_auc(K, y):
    seed = np.where(y)[0]; m = K[:, seed].sum(1) - np.where(y, np.diag(K), 0.0); m = m / (len(seed) - y)
    return auc(m[y], m[~y])

rows, t0 = [], time.time()
for ip, pat in enumerate(PATIENTS_4PHASE, 1):
    pm = build_epi_masks(pat); y = pm.epi_mask; N = y.size; probes = pm.probes
    same = build_probe_mask(list(pm.channels)); np.fill_diagonal(same, False)
    n_seed = int(y.sum()); shaft_ids, shaft_counts = np.unique(probes[y], return_counts=True)
    unmatched = [rng.choice(N, n_seed, replace=False) for _ in range(NDRAW)]
    matched = []
    for _ in range(NDRAW):
        pick = [rng.choice(np.where(probes == sh)[0], c, replace=False) for sh, c in zip(shaft_ids, shaft_counts)]
        matched.append(np.concatenate(pick))
    for band in BRAIN_BANDS_NAMES:
        W = load_phase_fc(pat, "rest_post", band)
        B = select_backbone(W, CANONICAL.backbone, frac=CANONICAL.frac); ev, V = laplacian_eig(B)
        for s in S_LIST:
            w = np.exp(-(s / ev[-1]) * ev); K = (V * w) @ V.T; K /= np.trace(K)
            for var, KK in (("all", K), ("xshaft", np.where(same, 0.0, K))):
                obs = seeded_auc(KK, y)
                for null_name, draws in (("unmatched", unmatched), ("shaftmatched", matched)):
                    nul = np.array([seeded_auc(KK, np.isin(np.arange(N), d)) for d in draws])
                    rows.append(dict(patient=pat, band=band, s=s, variant=var, null=null_name, n_soz=n_seed, N=N,
                                     obs=obs, null_med=np.median(nul), null_p95=np.percentile(nul, 95),
                                     z=(obs - nul.mean()) / (nul.std() + 1e-12), p=float(np.mean(nul >= obs))))
    print(f"[{ip}/10] {pat} {time.time()-t0:.0f}s", flush=True)

df = pd.DataFrame(rows); df.to_csv(OUT / "soz_seeded_null.csv", index=False)
g = df.groupby(["band", "s", "variant", "null"]).agg(obs_med=("obs", "median"), null_med=("null_med", "median"),
        null_p95=("null_p95", "median"), n_p05=("p", lambda p: int((p < 0.05).sum())), z_med=("z", "median")).reset_index()
pd.set_option("display.width", 200, "display.max_rows", 400)
print(g.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
print(f"\nwall {time.time()-t0:.0f}s")
