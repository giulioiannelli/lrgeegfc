#!/usr/bin/env python3
"""FEASIBILITY PROBE (2026-09-08) -- NOT A TEST. No null, no knob integration, no LOO.
Establishes nothing; it exists so the numbers quoted in
.agents/plans/active/2026-09-08_new-paths-to-objectives.md can be regenerated.
The methodology it probes is scoped (5-point preamble) in
.agents/guides/task-persistence-investigation/2026-09-08_hierarchy-level-decomposition.md.
Run from the worktree with PYTHONPATH=src and the lapbrain python. Output lands under
data/paper_final/feasibility/ (copy from the ROOT the script resolves).
"""
"""Feasibility: decompose the cross-phase trace by the HIERARCHICAL LEVEL of the
rest_pre diffusion tree. Pairs are binned by the number of clusters k remaining when
they merge in the UPGMA tree of 1/rho(s=1) built on the full rest_pre backbone.
Bins use DISJOINT pair sets, so per-level statistics are independent by construction
(unlike whole-vector rho_sym at different s, which share all pairs).

Readouts (non-ultrametric, so within-level variance exists):
  rawW : dense |ImCoh| adjacency
  heat1/heat5/heat16 : rho(s) = e^{-tau L}/Tr on the backbone at s = 1, 4.72, 15.9
Statistic per level: rho_sym restricted to the level's pairs. Descriptive feasibility:
NO null here; what matters is whether the level-profile is non-flat and consistent
across patients (that is the scale-specific question), not the absolute level.
Also reported with same-shaft pairs removed (the finest levels are mostly same-shaft).
"""
import time, numpy as np, pandas as pd
from scipy.cluster.hierarchy import cophenet
from scipy.stats import spearmanr
from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, linkage_at_scale, rho_sym
from lrg_eegfc.utils.io.patient import load_channel_labels
from lrg_eegfc.utils.probe import build_probe_mask
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.phase_graphs import load_phase_fc
from lrg_eegfc.workflow.substrate import CANONICAL
ROOT = setup_script_env()
PH = ("A", "B", "task_learn", "task_test", "rest_post")
EDGES = [1, 2, 4, 8, 16, 32, 64, 10**6]     # k in (1,2],(2,4],(4,8],...
LEV = [f"k{EDGES[i]+1}-{EDGES[i+1]}" if EDGES[i+1] < 10**6 else f"k>{EDGES[i]}" for i in range(len(EDGES)-1)]
SCALES = {"heat1": 1.0, "heat5": 4.72, "heat16": 15.9}

def heat(ev, V, s):
    w = np.exp(-(s / ev[-1]) * ev); r = (V * w) @ V.T; return r / np.trace(r)

rows = []; t0 = time.time()
for ip, pat in enumerate(PATIENTS_4PHASE, 1):
    labels = load_channel_labels(pat)
    for band in BRAIN_BANDS_NAMES:
        try:
            dense = {ph: load_phase_fc(pat, ph, band) for ph in PH}
            Wpre = load_phase_fc(pat, "rest_pre", band)
        except Exception as e:
            print(f"skip {pat}/{band}: {e}", flush=True); continue
        N = Wpre.shape[0]; iu = np.triu_indices(N, 1)
        same = build_probe_mask(labels)[iu] if (labels is not None and len(labels) == N) else np.zeros(iu[0].size, bool)
        # reference tree: full rest_pre backbone at s=1
        Bpre = select_backbone(Wpre, CANONICAL.backbone, frac=CANONICAL.frac)
        evp, Vp = laplacian_eig(Bpre)
        Z = linkage_at_scale(evp, Vp, 1.0)
        h = cophenet(Z)                                   # condensed heights
        merges_le = np.searchsorted(np.sort(Z[:, 2]), h, side="right")   # merges with height <= h_ij
        k = N - merges_le                                 # clusters remaining after the pair merges
        k = np.maximum(k, 1)
        lev = np.digitize(k, EDGES[1:], right=True)       # 0..len(LEV)-1
        # readouts
        bb = {ph: select_backbone(dense[ph], CANONICAL.backbone, frac=CANONICAL.frac) for ph in PH}
        eig = {ph: laplacian_eig(bb[ph]) for ph in PH}
        R = {"rawW": {ph: dense[ph][iu] for ph in PH}}
        for name, s in SCALES.items():
            R[name] = {ph: heat(*eig[ph], s)[iu] for ph in PH}
        for rname, X in R.items():
            for mask_name, keep in (("all", np.ones(iu[0].size, bool)), ("xshaft", ~same)):
                for li, lname in enumerate(LEV):
                    m = keep & (lev == li)
                    n = int(m.sum())
                    if n < 30:
                        val = np.nan; rel = np.nan
                    else:
                        val = rho_sym(X["A"][m], X["B"][m], X["task_test"][m], X["rest_post"][m])[0]
                        rel = spearmanr(X["A"][m], X["B"][m])[0]
                    rows.append(dict(patient=pat, band=band, readout=rname, mask=mask_name,
                                     level=lname, li=li, n_pairs=n, frac_same=float(same[lev == li].mean()) if (lev == li).any() else np.nan,
                                     rho_sym=val, rel_AB=rel))
    el = time.time() - t0
    print(f"[{ip}/{len(PATIENTS_4PHASE)}] {pat} {el:5.1f}s ETA {el/ip*(len(PATIENTS_4PHASE)-ip):5.1f}s", flush=True)

df = pd.DataFrame(rows)
out = ROOT / "data" / "paper_final" / "feasibility" / "level_decomposition.csv"
out.parent.mkdir(parents=True, exist_ok=True); df.to_csv(out, index=False)
pd.set_option("display.width", 220); pd.set_option("display.max_columns", 30)
print("\n=== pairs per level (cohort median) and same-shaft fraction ===")
g0 = df[(df.readout == "rawW") & (df["mask"] == "all")].groupby("level")[["n_pairs", "frac_same"]].median().reindex(LEV).round(2)
print(g0.T)
for mask_name in ("xshaft",):
    sub = df[df["mask"] == mask_name]
    for val, title in (("rho_sym", "rho_sym"), ("rel_AB", "split-half reliability Spearman(A,B)")):
        print(f"\n=== {title} by level, cohort MEDIAN, mask={mask_name} ===")
        print(sub.pivot_table(index=["band", "readout"], columns="level", values=val, aggfunc="median").reindex(columns=LEV).round(2))
    sub = sub.assign(ratio=sub.rho_sym / sub.rel_AB)
    print(f"\n=== rho_sym / reliability by level, cohort MEDIAN, mask={mask_name} ===")
    print(sub.pivot_table(index=["band", "readout"], columns="level", values="ratio", aggfunc="median").reindex(columns=LEV).round(2))
print(f"\nwall {time.time()-t0:.1f}s -> {out}")
