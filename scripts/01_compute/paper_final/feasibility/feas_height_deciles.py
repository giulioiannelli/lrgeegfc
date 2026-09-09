#!/usr/bin/env python3
"""FEASIBILITY PROBE (2026-09-09) -- NOT A TEST. No null.
Why: the k-level binning of lane_m1 is outlier-dominated at the top of the tree (k = 2..4 split off 1-7
peripheral contacts in 9/10 patients), so "coarse levels" mostly means "pairs touching peripheral contacts".
Here the scale axis is the cophenetic merge height itself (diffusion distance 1/rho at s = 1 on the rest_pre
backbone), binned into DECILES over the cross-shaft pairs of each patient (equal populations, no outlier
dominance). Profile: rho_sym per decile (heat-kernel s = 1 and raw readouts, canonical full phases), plus
split-half reliability per decile. Question: is the trace flat, monotone, or peaked along diffusion distance?
Output: data/paper_final/feasibility/height_deciles.csv. Run: PYTHONPATH=src <lapbrain python> -u <this file>
"""
import time, numpy as np, pandas as pd
from scipy.cluster.hierarchy import cophenet
from scipy.stats import spearmanr, wilcoxon
from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, linkage_at_scale, rho_sym
from lrg_eegfc.utils.io.patient import load_channel_labels
from lrg_eegfc.utils.probe import build_probe_mask
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.phase_graphs import load_phase_fc
from lrg_eegfc.workflow.substrate import CANONICAL
ROOT = setup_script_env(); OUT = ROOT / "data" / "paper_final" / "feasibility"
PH4 = ("A", "B", "task_test", "rest_post"); ND = 10; READ = ("rawW", "heat1")

def readout(W, iu):
    B = select_backbone(W, CANONICAL.backbone, frac=CANONICAL.frac); ev, V = laplacian_eig(B)
    w = np.exp(-(1.0 / ev[-1]) * ev); K = (V * w) @ V.T; return {"rawW": W[iu], "heat1": (K / np.trace(K))[iu]}

rows, t0 = [], time.time()
for ip, pat in enumerate(PATIENTS_4PHASE, 1):
    labels = load_channel_labels(pat)
    for band in BRAIN_BANDS_NAMES:
        Wpre = load_phase_fc(pat, "rest_pre", band); N = Wpre.shape[0]; iu = np.triu_indices(N, 1)
        keep = ~build_probe_mask(list(labels))[iu]
        B = select_backbone(Wpre, CANONICAL.backbone, frac=CANONICAL.frac); ev, V = laplacian_eig(B)
        h = cophenet(linkage_at_scale(ev, V, 1.0))
        q = np.quantile(h[keep], np.linspace(0, 1, ND + 1)); dec = np.clip(np.searchsorted(q, h, side="right") - 1, 0, ND - 1)
        R = {ph: readout(load_phase_fc(pat, ph, band), iu) for ph in PH4}
        for d in range(ND):
            m = keep & (dec == d)
            for r in READ:
                X = {ph: R[ph][r][m] for ph in PH4}
                rows.append(dict(patient=pat, band=band, readout=r, decile=d, n=int(m.sum()), h_med=float(np.median(h[m])),
                                 rho=rho_sym(X["A"], X["B"], X["task_test"], X["rest_post"])[0], rel=spearmanr(X["A"], X["B"])[0]))
    print(f"[{ip}/10] {pat} {time.time()-t0:.0f}s", flush=True)
df = pd.DataFrame(rows); df.to_csv(OUT / "height_deciles.csv", index=False)
pd.set_option("display.width", 220)
for r in READ:
    print(f"\n== {r}: median rho_sym by band x height decile (0 = closest in diffusion distance) ==")
    print(df[df.readout == r].pivot_table(index="band", columns="decile", values="rho", aggfunc="median").to_string(float_format=lambda x: f"{x:.2f}"))
print("\n== heat1: median split-half rel by band x decile ==")
print(df[df.readout == "heat1"].pivot_table(index="band", columns="decile", values="rel", aggfunc="median").to_string(float_format=lambda x: f"{x:.2f}"))
print("\n== heat1 beta: per-patient Spearman(decile, rho) and Wilcoxon across patients; and peak decile per patient ==")
for band in BRAIN_BANDS_NAMES:
    g = df[(df.readout == "heat1") & (df.band == band)]
    tr = g.groupby("patient").apply(lambda x: spearmanr(x.decile, x.rho)[0]).to_numpy(); pk = g.groupby("patient").apply(lambda x: int(x.decile.iloc[np.nanargmax(x.rho.to_numpy())])).to_numpy()
    print(f"  {band:11s} trend mean {tr.mean():+.2f} n_neg {(tr<0).sum()}/10 p {wilcoxon(tr).pvalue:.3f} | peak deciles {sorted(pk)}")
print(f"\nwall {time.time()-t0:.0f}s")
