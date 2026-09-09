#!/usr/bin/env python3
"""Lane M1 / objective 3 — half-scale windowed construction (higher SNR than the quarter-chunk sham of
lane_m1_02, where 2-min chunks put every rho_sym at ~0.05 and the sham test was uninformative).
Not a null; a duration-equated learn-vs-test comparison on the level axis.
Windows (30 s cache): A = rest_pre first half, B = rest_pre second half; task blocks equated to m = min(n_learn,
n_test) windows: learn_last (last m of task_learn), test_last (last m of task_test), test_first (first m);
post = rest_post first half (H1) or second half (H2). Level assignment = canonical full rest_pre tree at s = 1,
cross-shaft pairs, heat-kernel (s = 1) and raw readouts. Reports per-level profiles, the fine-level (k >= 5)
rho_sym and fine-coarse contrast, paired Wilcoxon test_last - learn_last, and H1 vs H2 (decay vs plateau).
Recency caveat: test_last is temporally closest to rest_post; test_first and learn_last are both farther, so
test_first vs learn_last is the fairer block comparison at equal duration and unequal (but reversed) recency.
Outputs: data/paper_final/lane_m1_level/halves_profiles.csv. Run: PYTHONPATH=src <lapbrain python> -u <this file>
"""
import time, numpy as np, pandas as pd
from scipy.stats import wilcoxon
from lrg_eegfc.config.const import PATIENTS_4PHASE, BRAIN_BANDS_NAMES
from lrg_eegfc.config.paths import CACHE_ROOT
from lrg_eegfc.utils.io.patient import load_channel_labels
from lrg_eegfc.utils.probe import build_probe_mask
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.workflow.phase_graphs import load_phase_fc
import importlib.util, sys
ROOT = setup_script_env(); OUT = ROOT / "data" / "paper_final" / "lane_m1_level"; WC = CACHE_ROOT / "windowed_imcoh"
spec = importlib.util.spec_from_file_location("m102", ROOT / "scripts/01_compute/paper_final/lane_m1_02_level_sham.py")
# reuse helpers without executing the sham script body
src = open(spec.origin).read().split("rows, lag, t0 = [], [], time.time()")[0]; ns = {}; exec(compile(src, spec.origin, "exec"), ns)
levels_from_tree, readout, profile, stats, LEV, READ, WIN = (ns[k] for k in ("levels_from_tree", "readout", "profile", "stats", "LEV", "READ", "WIN"))

rows, t0 = [], time.time()
for ip, pat in enumerate(PATIENTS_4PHASE, 1):
    Wc = {ph: np.load(WC / f"{pat}_{ph}_w{WIN}s.npz")["W"] for ph in ("rest_pre", "task_learn", "task_test", "rest_post")}
    n = {ph: Wc[ph].shape[0] for ph in Wc}; m = min(n["task_learn"], n["task_test"]); hp, hq = n["rest_pre"] // 2, n["rest_post"] // 2
    N = Wc["rest_pre"].shape[-1]; iu = np.triu_indices(N, 1); keep = ~build_probe_mask(list(load_channel_labels(pat)))[iu]
    mean = lambda ph, bi, sl: Wc[ph][sl, bi].mean(0).astype(float)
    for bi, band in enumerate(BRAIN_BANDS_NAMES):
        lev = levels_from_tree(load_phase_fc(pat, "rest_pre", band))
        A, B = readout(mean("rest_pre", bi, slice(0, hp)), iu), readout(mean("rest_pre", bi, slice(hp, 2 * hp)), iu)
        post = {"H1": readout(mean("rest_post", bi, slice(0, hq)), iu), "H2": readout(mean("rest_post", bi, slice(hq, 2 * hq)), iu)}
        m2 = m // 2   # m = n_learn in every patient (learn is the shorter block), so first/last-m of learn coincide; m2 halves give a within-block position check
        blocks = {"learn_first2": readout(mean("task_learn", bi, slice(0, m2)), iu), "learn_last2": readout(mean("task_learn", bi, slice(n["task_learn"] - m2, n["task_learn"])), iu),
                  "test_first2": readout(mean("task_test", bi, slice(0, m2)), iu), "test_last2": readout(mean("task_test", bi, slice(n["task_test"] - m2, n["task_test"])), iu),
                  "learn_last": readout(mean("task_learn", bi, slice(n["task_learn"] - m, n["task_learn"])), iu),
                  "test_last": readout(mean("task_test", bi, slice(n["task_test"] - m, n["task_test"])), iu),
                  "test_first": readout(mean("task_test", bi, slice(0, m)), iu)}
        for bn, T in blocks.items():
            for hn, P in post.items():
                rho, rel = profile(dict(A=A, B=B, task=T, post=P), lev, keep)
                for r in READ:
                    f, c = stats(rho[r])
                    rows.append(dict(patient=pat, band=band, block=bn, post=hn, readout=r, m=m, fine=f, contrast=c,
                                     **{f"rho_{l}": rho[r][i] for i, l in enumerate(LEV)}, **{f"rel_{l}": rel[r][i] for i, l in enumerate(LEV)}))
    print(f"[{ip}/10] {pat} m={m} half_pre={hp} half_post={hq} {time.time()-t0:.0f}s", flush=True)
df = pd.DataFrame(rows); df.to_csv(OUT / "halves_profiles.csv", index=False)
pd.set_option("display.width", 220, "display.max_rows", 500)
h = df[df.readout == "heat1"]
print("\n== heat1: median fine-level rho_sym (k>=5) by band x block x post half ==")
print(h.pivot_table(index="band", columns=["block", "post"], values="fine", aggfunc="median").to_string(float_format=lambda x: f"{x:.2f}"))
print("\n== heat1: median fine-coarse contrast ==")
print(h.pivot_table(index="band", columns=["block", "post"], values="contrast", aggfunc="median").to_string(float_format=lambda x: f"{x:.2f}"))
print("\n== paired Wilcoxon (heat1, fine-level rho): test_first - learn_last | test_last - learn_last | H2 - H1 (test_last) ==")
for band in BRAIN_BANDS_NAMES:
    g = h[h.band == band]
    p1 = g[g.post == "H1"].pivot(index="patient", columns="block", values="fine"); ph = g[g.block == "test_last"].pivot(index="patient", columns="post", values="fine")
    def pw(a, b):
        d = (p1[a] - p1[b]).dropna().to_numpy(); return f"{d.mean():+.3f} {int((d>0).sum())}/10 p {wilcoxon(d).pvalue:.3f}"
    print(f"  {band:11s} m: testF-learn {pw('test_first','learn_last')} | testF-testL {pw('test_first','test_last')} || m/2: testF-learnF {pw('test_first2','learn_first2')} | testF-learnL {pw('test_first2','learn_last2')} | learnF-learnL {pw('learn_first2','learn_last2')} | testF-testL {pw('test_first2','test_last2')}")
print("\n== beta heat1 median profile by block (post H1) ==")
print(h[(h.band == "beta") & (h.post == "H1")].groupby("block")[[f"rho_{l}" for l in LEV] + [f"rel_{l}" for l in LEV]].median().to_string(float_format=lambda x: f"{x:.2f}"))
print(f"\nwall {time.time()-t0:.0f}s")
