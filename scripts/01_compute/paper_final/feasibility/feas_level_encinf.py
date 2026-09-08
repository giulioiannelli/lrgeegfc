#!/usr/bin/env python3
"""FEASIBILITY PROBE (2026-09-08) -- NOT A TEST. No null, no knob integration, no LOO.
Establishes nothing; it exists so the numbers quoted in
.agents/plans/active/2026-09-08_new-paths-to-objectives.md can be regenerated.
The methodology it probes is scoped (5-point preamble) in
.agents/guides/task-persistence-investigation/2026-09-08_hierarchy-level-decomposition.md.
Run from the worktree with PYTHONPATH=src and the lapbrain python. Output lands under
data/paper_final/feasibility/ (copy from the ROOT the script resolves).
"""
"""Feasibility: encoding vs inference BY HIERARCHICAL LEVEL (disjoint pair sets).
Per level: T_test = rho(g,p), T_learn = rho(e,p), C = T_learn - T_test, with
e = X_learn - X_A, g = X_test - X_A, p = X_post - X_B (Spearman within the level's pairs,
symmetrised over the A/B role swap like rho_sym). Descriptive, no null, order-confounded:
the question is only whether the level PROFILE of C is non-flat / sign-changing.
"""
import time, numpy as np, pandas as pd
from scipy.cluster.hierarchy import cophenet
from scipy.stats import spearmanr
from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, linkage_at_scale
from lrg_eegfc.utils.io.patient import load_channel_labels
from lrg_eegfc.utils.probe import build_probe_mask
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.phase_graphs import load_phase_fc
from lrg_eegfc.workflow.substrate import CANONICAL
ROOT = setup_script_env()
PH = ("A", "B", "task_learn", "task_test", "rest_post")
EDGES = [1, 2, 4, 8, 16, 32, 64, 10**6]
LEV = [f"k{EDGES[i]+1}-{EDGES[i+1]}" if EDGES[i+1] < 10**6 else f"k>{EDGES[i]}" for i in range(len(EDGES)-1)]

def heat(ev, V, s):
    w = np.exp(-(s / ev[-1]) * ev); r = (V * w) @ V.T; return r / np.trace(r)

def sym(fa, fb):                        # average over the A/B role swap
    return 0.5 * (fa + fb)

rows = []; t0 = time.time()
for ip, pat in enumerate(PATIENTS_4PHASE, 1):
    labels = load_channel_labels(pat)
    for band in BRAIN_BANDS_NAMES:
        dense = {ph: load_phase_fc(pat, ph, band) for ph in PH}
        Wpre = load_phase_fc(pat, "rest_pre", band)
        N = Wpre.shape[0]; iu = np.triu_indices(N, 1)
        same = build_probe_mask(labels)[iu] if (labels is not None and len(labels) == N) else np.zeros(iu[0].size, bool)
        evp, Vp = laplacian_eig(select_backbone(Wpre, CANONICAL.backbone, frac=CANONICAL.frac))
        Z = linkage_at_scale(evp, Vp, 1.0); h = cophenet(Z)
        k = np.maximum(N - np.searchsorted(np.sort(Z[:, 2]), h, side="right"), 1)
        lev = np.digitize(k, EDGES[1:], right=True)
        eig = {ph: laplacian_eig(select_backbone(dense[ph], CANONICAL.backbone, frac=CANONICAL.frac)) for ph in PH}
        R = {"rawW": {ph: dense[ph][iu] for ph in PH}, "heat1": {ph: heat(*eig[ph], 1.0)[iu] for ph in PH}}
        for rname, X in R.items():
            for li, lname in enumerate(LEV):
                m = (~same) & (lev == li)
                if m.sum() < 30: continue
                out = {}
                for (a, b) in (("A", "B"), ("B", "A")):
                    e = X["task_learn"][m] - X[a][m]; g = X["task_test"][m] - X[a][m]; p = X["rest_post"][m] - X[b][m]
                    out.setdefault("T_test", []).append(spearmanr(g, p)[0])
                    out.setdefault("T_learn", []).append(spearmanr(e, p)[0])
                Tt = float(np.mean(out["T_test"])); Tl = float(np.mean(out["T_learn"]))
                rows.append(dict(patient=pat, band=band, readout=rname, level=lname, li=li,
                                 T_test=Tt, T_learn=Tl, C=Tl - Tt))
    el = time.time() - t0
    print(f"[{ip}/{len(PATIENTS_4PHASE)}] {pat} {el:5.1f}s ETA {el/ip*(len(PATIENTS_4PHASE)-ip):5.1f}s", flush=True)

df = pd.DataFrame(rows)
out = ROOT / "data" / "paper_final" / "feasibility" / "level_encinf.csv"; out.parent.mkdir(parents=True, exist_ok=True); df.to_csv(out, index=False)
pd.set_option("display.width", 220)
for val in ("T_learn", "T_test", "C"):
    print(f"\n=== {val} by level, cohort MEDIAN, cross-shaft pairs ===")
    print(df.pivot_table(index=["band", "readout"], columns="level", values=val, aggfunc="median").reindex(columns=LEV).round(2))
print("\n=== C = T_learn - T_test: n patients with C>0 (of 10), readout=heat1 ===")
print(df[df.readout == "heat1"].groupby(["band", "level"])["C"].apply(lambda v: int((v > 0).sum())).unstack("level").reindex(columns=LEV))
print(f"\nwall {time.time()-t0:.1f}s -> {out}")
