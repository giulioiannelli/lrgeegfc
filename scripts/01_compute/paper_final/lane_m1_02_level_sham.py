#!/usr/bin/env python3
"""Lane M1, sham rung — per-level ORDERED SHAM on the windowed |ImCoh| cache, plus the time-equated
learn-vs-test reinstatement profile (objective 3, windowed construction).
Scope: .agents/guides/task-persistence-investigation/2026-09-08_hierarchy-level-decomposition.md (caveat: reliability gradient).

Preamble.
(1) Claim: the fine-level (k >= 5) task_test -> rest_post rho_sym, and the fine-minus-coarse step, exceed what
    the same construction gives when the four "phases" are four consecutive chunks of ONE resting recording.
(2) Null (sham): pseudo-phases A, B, task, post := consecutive chunks Q1..Q4 of rest_pre (sham_pre) or of rest_post
    (sham_post), each nq windows of 30 s, the level assignment and the readouts computed exactly as for the data.
    Observed is built from the SAME window counts (A = pre[0:nq], B = pre[nq:2nq], task = task_test last nq or
    first nq windows, post = rest_post first nq windows), so duration/SNR is equated (no series truncation of the
    canonical full-phase result; this is a separate, windowed construction).
(3) Strongest alternative: fine-level pairs are simply more reproducible (within-module), so any two later chunks
    correlate more at fine levels than coarse ones -> a level step with no task content.
(4) The sham carries exactly that reproducibility gradient and the slow drift of a resting recording; it does not
    carry task-locked arousal or time-of-day effects that differ between rest_pre and rest_post (the matched-
    strength and tree-only nulls of lane_m1 do not either). It cannot reject "any state change, not this task".
(5) Falsifiers: sham profile shows the same step (paired Wilcoxon obs - sham on the contrast p > 0.05); fine-level
    excess not positive in >= 8/10 patients; nq too small for a stable reference (rel < 0.3 at fine levels).
Objective 3 add-on: rho_sym per level for task_learn (last nq) vs task_test (last nq, first nq) against rest_post
chunk j = 0..3, on the same level assignment: which block does rest_post reinstate, and does it decay with j.
Outputs: data/paper_final/lane_m1_level/sham_profiles.csv, sham_tests.csv, reinstatement_lag.csv.
Run: PYTHONPATH=src <lapbrain python> -u <this file>
"""
import time, numpy as np, pandas as pd
from scipy.cluster.hierarchy import cophenet
from scipy.stats import spearmanr, wilcoxon
from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, linkage_at_scale, rho_sym
from lrg_eegfc.utils.io.patient import load_channel_labels
from lrg_eegfc.utils.probe import build_probe_mask
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.phase_graphs import load_phase_fc
from lrg_eegfc.workflow.substrate import CANONICAL
ROOT = setup_script_env(); OUT = ROOT / "data" / "paper_final" / "lane_m1_level"; WC = CACHE_ROOT / "windowed_imcoh"
EDGES = [1, 2, 4, 8, 16, 32, 64, 10**6]; NLEV = 7; FINE0 = 2; LEV = ["k2", "k3-4", "k5-8", "k9-16", "k17-32", "k33-64", "k>64"]
READ = ("rawW", "heat1"); MIN_PAIRS = 30; S_READ = 1.0; WIN = 30

def levels_from_tree(Wpre):
    B = select_backbone(Wpre, CANONICAL.backbone, frac=CANONICAL.frac); ev, V = laplacian_eig(B); N = B.shape[0]
    Z = linkage_at_scale(ev, V, 1.0); h = cophenet(Z)
    return np.digitize(np.maximum(N - np.searchsorted(np.sort(Z[:, 2]), h, side="right"), 1), EDGES[1:], right=True)

def readout(W, iu):
    B = select_backbone(W, CANONICAL.backbone, frac=CANONICAL.frac); ev, V = laplacian_eig(B)
    w = np.exp(-(S_READ / ev[-1]) * ev); K = (V * w) @ V.T
    return {"rawW": W[iu], "heat1": (K / np.trace(K))[iu]}

def profile(R, lev, keep):
    """R: dict phase->dict readout->condensed. -> rho[readout](NLEV,), rel[readout](NLEV,)"""
    rho = {r: np.full(NLEV, np.nan) for r in READ}; rel = {r: np.full(NLEV, np.nan) for r in READ}
    for li in range(NLEV):
        m = keep & (lev == li)
        if m.sum() < MIN_PAIRS: continue
        for r in READ:
            rho[r][li] = rho_sym(R["A"][r][m], R["B"][r][m], R["task"][r][m], R["post"][r][m])[0]
            rel[r][li] = spearmanr(R["A"][r][m], R["B"][r][m])[0]
    return rho, rel

def stats(p):
    ok = ~np.isnan(p); fine = np.nanmean(p[FINE0:]) if ok[FINE0:].any() else np.nan
    coarse = np.nanmean(p[:FINE0]) if ok[:FINE0].any() else np.nan
    return fine, fine - coarse

rows, lag, t0 = [], [], time.time()
for ip, pat in enumerate(PATIENTS_4PHASE, 1):
    Wc = {ph: np.load(WC / f"{pat}_{ph}_w{WIN}s.npz")["W"] for ph in ("rest_pre", "task_learn", "task_test", "rest_post")}
    n = {ph: Wc[ph].shape[0] for ph in Wc}; nq = min(n["rest_pre"] // 4, n["rest_post"] // 4, n["task_learn"], n["task_test"])
    N = Wc["rest_pre"].shape[-1]; iu = np.triu_indices(N, 1); keep = ~build_probe_mask(list(load_channel_labels(pat)))[iu]
    chunk = lambda ph, bi, j: Wc[ph][j * nq:(j + 1) * nq, bi].mean(0).astype(float)
    for bi, band in enumerate(BRAIN_BANDS_NAMES):
        lev = levels_from_tree(load_phase_fc(pat, "rest_pre", band))          # reference tree = canonical full rest_pre
        A, B = readout(chunk("rest_pre", bi, 0), iu), readout(chunk("rest_pre", bi, 1), iu)
        postQ = [readout(chunk("rest_post", bi, j), iu) for j in range(n["rest_post"] // nq)]
        preQ = [readout(chunk("rest_pre", bi, j), iu) for j in range(4)]
        test_last = readout(Wc["task_test"][-nq:, bi].mean(0).astype(float), iu); test_first = readout(chunk("task_test", bi, 0), iu)
        learn_last = readout(Wc["task_learn"][-nq:, bi].mean(0).astype(float), iu)
        cons = {"obs_test_last": dict(A=A, B=B, task=test_last, post=postQ[0]), "obs_test_first": dict(A=A, B=B, task=test_first, post=postQ[0]),
                "obs_learn_last": dict(A=A, B=B, task=learn_last, post=postQ[0]),
                "sham_pre": dict(A=preQ[0], B=preQ[1], task=preQ[2], post=preQ[3]),
                "sham_post": dict(A=postQ[0], B=postQ[1], task=postQ[2], post=postQ[3])}
        for cname, R in cons.items():
            rho, rel = profile(R, lev, keep)
            for r in READ:
                f, c = stats(rho[r])
                rows.append(dict(patient=pat, band=band, construction=cname, readout=r, nq=nq, fine=f, contrast=c,
                                 **{f"rho_{LEV[l]}": rho[r][l] for l in range(NLEV)}, **{f"rel_{LEV[l]}": rel[r][l] for l in range(NLEV)}))
        for j, P in enumerate(postQ):                                             # reinstatement lag profile
            for tname, T in (("learn_last", learn_last), ("test_last", test_last), ("test_first", test_first)):
                rho, _ = profile(dict(A=A, B=B, task=T, post=P), lev, keep)
                for r in READ:
                    f, c = stats(rho[r]); lag.append(dict(patient=pat, band=band, block=tname, post_chunk=j, readout=r, nq=nq, fine=f, contrast=c, whole=np.nan))
    print(f"[{ip}/10] {pat} nq={nq} n_win={n} {time.time()-t0:.0f}s", flush=True)

df = pd.DataFrame(rows); df.to_csv(OUT / "sham_profiles.csv", index=False)
L = pd.DataFrame(lag); L.to_csv(OUT / "reinstatement_lag.csv", index=False)
tests = []
for (band, r), g in df.groupby(["band", "readout"]):
    piv = {st: g.pivot(index="patient", columns="construction", values=st) for st in ("fine", "contrast")}
    for obs in ("obs_test_last", "obs_test_first", "obs_learn_last"):
        for sham in ("sham_pre", "sham_post"):
            for st in ("fine", "contrast"):
                d = (piv[st][obs] - piv[st][sham]).dropna().to_numpy()
                p = wilcoxon(d).pvalue if d.size >= 5 else np.nan
                loo = max(wilcoxon(np.delete(d, j)).pvalue for j in range(d.size)) if d.size >= 6 else np.nan
                tests.append(dict(band=band, readout=r, obs=obs, sham=sham, stat=st, n=d.size, obs_mean=piv[st][obs].mean(), sham_mean=piv[st][sham].mean(),
                                  diff_mean=d.mean(), n_pos=int((d > 0).sum()), p_wilcoxon=p, p_loo_max=loo))
T = pd.DataFrame(tests); T.to_csv(OUT / "sham_tests.csv", index=False)
pd.set_option("display.width", 220, "display.max_rows", 500)
print("\n== obs_test_last vs shams (heat1, xshaft): fine-level rho and fine-coarse contrast ==")
print(T[(T.readout == "heat1") & (T.obs == "obs_test_last")].to_string(index=False, float_format=lambda x: f"{x:.3f}"))
print("\n== beta, all constructions: cohort median level profile ==")
print(df[df.band == "beta"].groupby(["readout", "construction"])[[f"rho_{l}" for l in LEV] + ["fine", "contrast"]].median().to_string(float_format=lambda x: f"{x:.2f}"))
print("\n== reinstatement lag (heat1): median fine-level rho by block x post chunk ==")
print(L[L.readout == "heat1"].pivot_table(index=["band", "block"], columns="post_chunk", values="fine", aggfunc="median").to_string(float_format=lambda x: f"{x:.2f}"))
print("\n== learn_last vs test_last, chunk 0 (heat1): paired Wilcoxon on fine-level rho ==")
for band in BRAIN_BANDS_NAMES:
    g = L[(L.readout == "heat1") & (L.post_chunk == 0) & (L.band == band)].pivot(index="patient", columns="block", values="fine")
    d = (g["test_last"] - g["learn_last"]).dropna().to_numpy()
    print(f"  {band:11s} test-learn mean {d.mean():+.3f} n_pos {(d>0).sum()}/{d.size} p {wilcoxon(d).pvalue:.3f}")
print(f"\nwall {time.time()-t0:.0f}s -> {OUT}")
