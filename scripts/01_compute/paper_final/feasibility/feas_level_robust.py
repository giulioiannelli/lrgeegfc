#!/usr/bin/env python3
"""FEASIBILITY PROBE (2026-09-08) -- NOT A TEST. No null, no knob integration, no LOO.
Establishes nothing; it exists so the numbers quoted in
.agents/plans/active/2026-09-08_new-paths-to-objectives.md can be regenerated.
The methodology it probes is scoped (5-point preamble) in
.agents/guides/task-persistence-investigation/2026-09-08_hierarchy-level-decomposition.md.
Run from the worktree with PYTHONPATH=src and the lapbrain python. Output lands under
data/paper_final/feasibility/ (copy from the ROOT the script resolves).
"""
"""Robustness of the hierarchy-level gradient to the backbone fraction and the
reference-tree scale. Same construction as feas_level_decomp.py; only the
within-patient trend statistics are reported (beta, theta, alpha, low_gamma)."""
import time, sys, numpy as np, pandas as pd
from scipy.cluster.hierarchy import cophenet
from scipy.stats import spearmanr, wilcoxon
from lrg_eegfc.config.const import PATIENTS_4PHASE
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, linkage_at_scale, rho_sym
from lrg_eegfc.utils.io.patient import load_channel_labels
from lrg_eegfc.utils.probe import build_probe_mask
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.phase_graphs import load_phase_fc
ROOT = setup_script_env()
PH = ("A", "B", "task_learn", "task_test", "rest_post")
EDGES = [1, 2, 4, 8, 16, 32, 64, 10**6]
BANDS = ("alpha", "beta", "theta", "low_gamma")
FRACS = (0.07, 0.14, 0.20); SREF = (1.0, 4.72)

def heat(ev, V, s):
    w = np.exp(-(s / ev[-1]) * ev); r = (V * w) @ V.T; return r / np.trace(r)

rows = []; t0 = time.time()
for ip, pat in enumerate(PATIENTS_4PHASE, 1):
    labels = load_channel_labels(pat)
    for band in BANDS:
        dense = {ph: load_phase_fc(pat, ph, band) for ph in PH}
        Wpre = load_phase_fc(pat, "rest_pre", band)
        N = Wpre.shape[0]; iu = np.triu_indices(N, 1)
        same = build_probe_mask(labels)[iu] if (labels is not None and len(labels) == N) else np.zeros(iu[0].size, bool)
        for frac in FRACS:
            evp, Vp = laplacian_eig(select_backbone(Wpre, "mst", frac=frac))
            eig = {ph: laplacian_eig(select_backbone(dense[ph], "mst", frac=frac)) for ph in PH}
            R = {"rawW": {ph: dense[ph][iu] for ph in PH}, "heat1": {ph: heat(*eig[ph], 1.0)[iu] for ph in PH}}
            for sref in SREF:
                Z = linkage_at_scale(evp, Vp, sref); h = cophenet(Z)
                k = np.maximum(N - np.searchsorted(np.sort(Z[:, 2]), h, side="right"), 1)
                lev = np.digitize(k, EDGES[1:], right=True)
                for rname, X in R.items():
                    prof = []
                    for li in range(len(EDGES) - 1):
                        m = (~same) & (lev == li)
                        prof.append(rho_sym(X["A"][m], X["B"][m], X["task_test"][m], X["rest_post"][m])[0] if m.sum() >= 30 else np.nan)
                    prof = np.array(prof); ok = np.isfinite(prof)
                    slope = spearmanr(np.arange(len(prof))[ok], prof[ok])[0] if ok.sum() >= 5 else np.nan
                    cf = np.nanmean(prof[4:]) - np.nanmean(prof[:2])
                    rows.append(dict(patient=pat, band=band, frac=frac, sref=sref, readout=rname, slope=slope, fine_minus_coarse=cf))
    el = time.time() - t0
    print(f"[{ip}/10] {pat} {el:5.1f}s ETA {el/ip*(10-ip):5.1f}s", flush=True)

df = pd.DataFrame(rows); df.to_csv(ROOT / "data" / "paper_final" / "feasibility" / "level_robustness.csv", index=False)
pd.set_option("display.width", 220)
def summ(g):
    s = g.slope.dropna().to_numpy(); c = g.fine_minus_coarse.dropna().to_numpy()
    return pd.Series(dict(med_slope=np.median(s), n_pos=int((s > 0).sum()), p_slope=wilcoxon(s, alternative="greater").pvalue,
                          med_cf=np.median(c), p_cf=wilcoxon(c, alternative="greater").pvalue))
print("\n=== level-gradient robustness: rows = band x readout, cols = (frac, sref) ===")
out = df.groupby(["band", "readout", "frac", "sref"]).apply(summ).round(3)
for val in ("p_slope", "med_cf", "p_cf"):
    print(f"\n--- {val} ---"); print(out[val].unstack(["frac", "sref"]))
print(f"\nwall {time.time()-t0:.1f}s")
