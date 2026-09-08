#!/usr/bin/env python3
"""FEASIBILITY PROBE (2026-09-08) -- NOT A TEST. No null, no knob integration, no LOO.
Establishes nothing; it exists so the numbers quoted in
.agents/plans/active/2026-09-08_new-paths-to-objectives.md can be regenerated.
The methodology it probes is scoped (5-point preamble) in
.agents/guides/task-persistence-investigation/2026-09-08_hierarchy-level-decomposition.md.
Run from the worktree with PYTHONPATH=src and the lapbrain python. Output lands under
data/paper_final/feasibility/ (copy from the ROOT the script resolves).
"""
"""Feasibility: does a scale-LOCAL readout of the LRG heat kernel carry independent,
reliable content across s where the cophenetic vector does not?

Readouts per phase X, scale s (tau = s/lambda_max(X), same convention as cophenetic_at_scale):
  coph : cophenetic vector of the UPGMA tree of 1/rho_hat(tau)      [cumulative]
  heat : rho_hat(tau) = e^{-tau L}/Tr, upper triangle                 [cumulative low-pass]
  wave : -d rho_hat / d ln tau = tau * sum_k (lam_k - <lam>_tau) e^{-tau lam_k} phi_k phi_k^T / Z   [band-pass]
Quantities per patient x band x readout:
  n_eff_e : participation ratio of the cross-scale Spearman matrix of e_s = X_task - X_A
  n_eff_p : same for p_s = X_post - X_B
  rel(s)  : Spearman(X_A(s), X_B(s))  split-half reliability of the readout itself
  rho_sym(s) : the trace statistic (descriptive only here, no null)
"""
import time, numpy as np, pandas as pd
from scipy.stats import spearmanr
from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, cophenetic_at_scale, rho_sym
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.phase_graphs import load_phase_fc
from lrg_eegfc.workflow.substrate import CANONICAL
ROOT = setup_script_env()
S = np.logspace(np.log10(0.02), np.log10(180.0), 16)
PH = ("A", "B", "task_learn", "task_test", "rest_post")
iu = None

def heat_and_wave(ev, V, s):
    tau = s / ev[-1]
    w = np.exp(-tau * ev); Z = w.sum()
    rho = (V * w) @ V.T / Z
    lam_mean = (ev * w).sum() / Z
    wave = tau * ((V * ((ev - lam_mean) * w)) @ V.T) / Z
    return rho[iu], wave[iu]

def pr(C):
    lam = np.linalg.eigvalsh(C); lam = np.clip(lam, 0, None)
    return float(lam.sum()**2 / (lam**2).sum())

def xscale_neff(vecs):                     # vecs: list over s of 1-D arrays
    n = len(vecs); C = np.eye(n)
    for i in range(n):
        for j in range(i+1, n):
            C[i, j] = C[j, i] = spearmanr(vecs[i], vecs[j])[0]
    return pr(np.nan_to_num(C))

rows = []; t0 = time.time()
for ip, pat in enumerate(PATIENTS_4PHASE, 1):
    for band in BRAIN_BANDS_NAMES:
        try:
            Ws = {ph: select_backbone(load_phase_fc(pat, ph, band), CANONICAL.backbone, frac=CANONICAL.frac) for ph in PH}
        except Exception as e:
            print(f"skip {pat}/{band}: {e}", flush=True); continue
        N = Ws["A"].shape[0]; iu = np.triu_indices(N, 1)
        eig = {ph: laplacian_eig(W) for ph, W in Ws.items()}
        X = {r: {ph: [] for ph in PH} for r in ("coph", "heat", "wave")}
        for s in S:
            for ph in PH:
                ev, V = eig[ph]
                X["coph"][ph].append(cophenetic_at_scale(ev, V, s))
                h, w = heat_and_wave(ev, V, s)
                X["heat"][ph].append(h); X["wave"][ph].append(w)
        for r in X:
            e = [X[r]["task_test"][k] - X[r]["A"][k] for k in range(len(S))]
            p = [X[r]["rest_post"][k] - X[r]["B"][k] for k in range(len(S))]
            rel = [spearmanr(X[r]["A"][k], X[r]["B"][k])[0] for k in range(len(S))]
            rs = [rho_sym(X[r]["A"][k], X[r]["B"][k], X[r]["task_test"][k], X[r]["rest_post"][k])[0] for k in range(len(S))]
            rows.append(dict(patient=pat, band=band, readout=r, N=N,
                             n_eff_e=xscale_neff(e), n_eff_p=xscale_neff(p),
                             n_eff_rs=None,
                             rel_med=float(np.median(rel)), rel_min=float(np.min(rel)), rel_max=float(np.max(rel)),
                             **{f"rel_s{k}": rel[k] for k in range(len(S))},
                             **{f"rs_s{k}": rs[k] for k in range(len(S))}))
    el = time.time() - t0
    print(f"[{ip}/{len(PATIENTS_4PHASE)}] {pat} {el:5.1f}s ETA {el/ip*(len(PATIENTS_4PHASE)-ip):5.1f}s", flush=True)

df = pd.DataFrame(rows)
out = ROOT / "data" / "paper_final" / "feasibility_wavelet_scale.csv"
out.parent.mkdir(parents=True, exist_ok=True); df.to_csv(out, index=False)
pd.set_option("display.width", 200)
print("\n=== cohort MEDIAN over patients: effective number of independent scales (16-point grid, 4 decades) ===")
g = df.groupby(["band", "readout"])[["n_eff_e", "n_eff_p", "rel_med", "rel_min"]].median().round(2)
print(g.unstack("readout"))
print("\n=== reliability profile rel(s), cohort median, alpha & beta ===")
for band in ("alpha", "beta"):
    for r in ("coph", "heat", "wave"):
        v = df[(df.band == band) & (df.readout == r)][[f"rel_s{k}" for k in range(len(S))]].median().to_numpy()
        print(f"{band:5s} {r:5s} " + " ".join(f"{x:5.2f}" for x in v))
print("s grid: " + " ".join(f"{x:5.2f}" for x in S))
print("\n=== rho_sym(s) cohort median, alpha & beta (descriptive, no null) ===")
for band in ("alpha", "beta"):
    for r in ("coph", "heat", "wave"):
        v = df[(df.band == band) & (df.readout == r)][[f"rs_s{k}" for k in range(len(S))]].median().to_numpy()
        print(f"{band:5s} {r:5s} " + " ".join(f"{x:5.2f}" for x in v))
print(f"\nwall {time.time()-t0:.1f}s -> {out}")
