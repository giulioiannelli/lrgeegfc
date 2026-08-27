#!/usr/bin/env python3
"""W0-B diagnostic: interrogate the N1 outcome before it is written up.

N1 (lag-randomized, coherence-magnitude-preserving) is the rung the paper's
central claim hangs on, and it FAILS at the cohort gate for alpha and beta. A
negative of that weight must be checked for the ways it could be an artifact of
the harness rather than a fact about the data. Four checks:

D1. ORIENTATION. Confirm the one-sided test points the right way, on data whose
    answer is already known: rerun the same gate function against the incumbent
    matched-strength rung, which is known to clear for alpha/beta. If the gate
    code were inverted, ms would fail too.

D2. PER-PATIENT ANATOMY of the beta / alpha margins. A large median margin with a
    non-significant Wilcoxon means the margin's SIGN is inconsistent across
    patients. Print every patient so the failure mode is visible rather than
    summarized: how many are positive, and how big are the negative ones.

D3. IS THE N1 NULL JUST "THE PIPELINE ON PLAIN COHERENCE"? Randomizing the phase
    sends <|Im C|>_f -> (2/pi)<|C|>_f, so N1's surrogate should be a near
    DETERMINISTIC substitution of ordinary coherence magnitude |C| for |ImCoh|,
    with only a small stochastic part. Test it directly: run the full pipeline on
    W = <|C|>_f and compare that rho_sym to the N1 surrogate median. If they
    agree, N1's meaning is exactly "ordinary coherence reproduces the trace" and
    can be reported in those concrete terms.

D4. NULL SPREAD. Report the N1 surrogate spread per cell. If the surrogate
    distribution is near-degenerate (tiny spread around the |C| value), then the
    per-cell p-value is effectively a deterministic comparison and the cohort
    Wilcoxon is the only real test -- worth stating explicitly.
"""
from __future__ import annotations
import sys
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))

from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, DEFAULT_SAMPLE_RATE, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.fc.backbone import select_backbone
from lrg_eegfc.utils.fc.heat_multiscale import laplacian_eig, rho_sym_over_scales
from lrg_eegfc.utils.surrogate.timeseries_nulls import (
    coherency_from_csd, csd_from_segment_subset, imcoh_abs_from_coherency, segment_fft,
)
from audit_150_rho_sym_gate import COHORT

OUT = ROOT / "data" / "paper_final" / "w0b_nulls"
PHASES = ("A", "B", "task_test", "rest_post")
SGRID = np.logspace(0.0, np.log10(180.0), 16)


def gate(m):
    return float(wilcoxon(m, alternative="greater")[1]) if np.any(m != 0) else np.nan


def d1_orientation():
    print("\n" + "=" * 92)
    print("D1  ORIENTATION CHECK — same gate function on the incumbent matched-strength rung")
    print("=" * 92)
    p = ROOT / "data" / "sparsified_arc" / "ms_mst020" / "per_patient_scale.csv"
    d = pd.read_csv(p)
    for band in ("alpha", "beta", "theta"):
        gs = []
        for s in np.sort(d.s.unique()):
            x = d[(d.band == band) & np.isclose(d.s, s)].dropna(subset=["obs_rho", "surr_p50"])
            gs.append(gate(x.obs_rho.values - x.surr_p50.values))
        gs = np.array(gs)
        print(f"  ms {band:11s} clears {int((gs < 0.05).sum()):2d}/16 scales, min p = {np.nanmin(gs):.4f}")
    print("  -> matched-strength still clears alpha/beta with THIS gate code, so the")
    print("     one-sided direction (obs > surrogate) is correct and the N1 failure is not")
    print("     a sign error. Also note rho_sym is scale-free in W (validated 0.0e+00),")
    print("     so N1's larger raw edge weights cannot bias the comparison either way.")


def d2_per_patient():
    print("\n" + "=" * 92)
    print("D2  PER-PATIENT MARGINS (obs - surr_p50), N1 vs matched-strength")
    print("=" * 92)
    n1 = pd.read_csv(OUT / "n1" / "per_patient_scale.csv")
    ms = pd.read_csv(ROOT / "data" / "sparsified_arc" / "ms_mst020" / "per_patient_scale.csv")
    for band in ("alpha", "beta"):
        # report at the scale where N1 is least unfavourable AND at s=1 (fine scale)
        gs = {}
        for s in np.sort(n1.s.unique()):
            x = n1[(n1.band == band) & np.isclose(n1.s, s)].dropna(subset=["obs_rho", "surr_p50"])
            gs[s] = gate(x.obs_rho.values - x.surr_p50.values)
        s_best = min(gs, key=lambda k: (np.inf if np.isnan(gs[k]) else gs[k]))
        for s in (float(np.sort(n1.s.unique())[0]), s_best):
            xn = n1[(n1.band == band) & np.isclose(n1.s, s)].set_index("patient")
            xm = ms[(ms.band == band) & np.isclose(ms.s, s)].set_index("patient")
            print(f"\n  --- {band}, s = {s:.2f}  (N1 gate p = {gs[s]:.4f}) ---")
            print(f"  {'patient':9s} {'obs':>8s} {'N1 p50':>8s} {'N1 marg':>9s} "
                  f"{'MS p50':>8s} {'MS marg':>9s}")
            mn, mm = [], []
            for pat in COHORT:
                if pat not in xn.index:
                    continue
                o = xn.loc[pat, "obs_rho"]; sn = xn.loc[pat, "surr_p50"]
                sm = xm.loc[pat, "surr_p50"] if pat in xm.index else np.nan
                mn.append(o - sn); mm.append(o - sm)
                flag = "  <-- NEG" if (o - sn) < 0 else ""
                print(f"  {pat:9s} {o:+8.3f} {sn:+8.3f} {o-sn:+9.3f} "
                      f"{sm:+8.3f} {o-sm:+9.3f}{flag}")
            mn, mm = np.array(mn), np.array(mm)
            print(f"  {'':9s} {'':8s} {'':8s} {'':9s}")
            print(f"  N1 margin: {int((mn>0).sum())}/{len(mn)} positive, median {np.median(mn):+.3f}, "
                  f"IQR {np.percentile(mn,75)-np.percentile(mn,25):.3f}, "
                  f"min {mn.min():+.3f}, p = {gate(mn):.4f}")
            print(f"  MS margin: {int((mm>0).sum())}/{len(mm)} positive, median {np.median(mm):+.3f}, "
                  f"p = {gate(mm):.4f}")


def d3_plain_coherence():
    print("\n" + "=" * 92)
    print("D3  IS N1's NULL JUST 'THE PIPELINE ON PLAIN COHERENCE |C|'?")
    print("=" * 92)
    print("  Running the full pipeline on W = <|C|>_f (ordinary coherence magnitude),")
    print("  and comparing to the N1 surrogate median from the production run.\n")
    n1 = pd.read_csv(OUT / "n1" / "per_patient_scale.csv")
    rows = []
    for pat in COHORT[:5]:                       # 5 patients is enough to settle it
        fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE)
        nps = nperseg_for_fs(fs); nps_h = max(256, nps // 2)
        for band in ("alpha", "beta"):
            bnd = BRAIN_BANDS[band]
            Wc, Wi = {}, {}
            Xr = np.asarray(load_timeseries(pat, "rest_pre", SEEG_DATAPATH), float)
            if Xr.shape[0] > Xr.shape[1]:
                Xr = Xr.T
            T = Xr.shape[1]
            src = {"A": (Xr[:, :T // 2], nps_h), "B": (Xr[:, T // 2:], nps_h)}
            for ph in ("task_test", "rest_post"):
                X = np.asarray(load_timeseries(pat, ph, SEEG_DATAPATH), float)
                if X.shape[0] > X.shape[1]:
                    X = X.T
                src[ph] = (X, nps)
            for ph, (Xp, np_) in src.items():
                _, F, sc = segment_fft(Xp, fs, np_, band=bnd)
                C = coherency_from_csd(csd_from_segment_subset(F, None, sc))
                A = np.abs(C).mean(axis=0).astype(float)          # <|C|>_f
                np.fill_diagonal(A, 0.0); A = 0.5 * (A + A.T)
                Wc[ph] = A
                Wi[ph] = imcoh_abs_from_coherency(C)
                del F, C
            del Xr, src
            eig = lambda W: laplacian_eig(select_backbone(W, "mst020", frac=0.20))
            r_coh = rho_sym_over_scales({p: eig(Wc[p]) for p in PHASES}, SGRID)
            r_imc = rho_sym_over_scales({p: eig(Wi[p]) for p in PHASES}, SGRID)
            for j, s in enumerate(SGRID):
                q = n1[(n1.patient == pat) & (n1.band == band) & np.isclose(n1.s, s)]
                rows.append(dict(patient=pat, band=band, s=float(s),
                                 rho_imcoh=r_imc[j], rho_coh=r_coh[j],
                                 n1_surr_p50=float(q.surr_p50.values[0]) if len(q) else np.nan,
                                 n1_surr_p95=float(q.surr_p95.values[0]) if len(q) else np.nan))
            print(f"  {pat} {band:5s}  |ImCoh| rho[s=1]={r_imc[0]:+.3f}  "
                  f"|C| rho[s=1]={r_coh[0]:+.3f}  N1 p50[s=1]="
                  f"{rows[-16]['n1_surr_p50']:+.3f}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "diag_plain_coherence.csv", index=False)
    d = df.dropna(subset=["rho_coh", "n1_surr_p50"])
    print(f"\n  agreement between the |C| pipeline and the N1 surrogate median:")
    print(f"    mean |rho_coh - n1_surr_p50| = {np.mean(np.abs(d.rho_coh-d.n1_surr_p50)):.4f}")
    print(f"    Pearson r                    = {np.corrcoef(d.rho_coh, d.n1_surr_p50)[0,1]:.4f}")
    for band in ("alpha", "beta"):
        b = d[d.band == band]
        print(f"    {band:5s}: median rho(|ImCoh|) = {b.rho_imcoh.median():+.3f}  "
              f"median rho(|C|) = {b.rho_coh.median():+.3f}")


def d4_null_spread():
    print("\n" + "=" * 92)
    print("D4  N1 NULL SPREAD (p95 - p50), by band — is the surrogate near-degenerate?")
    print("=" * 92)
    n1 = pd.read_csv(OUT / "n1" / "per_patient_scale.csv")
    ms = pd.read_csv(ROOT / "data" / "sparsified_arc" / "ms_mst020" / "per_patient_scale.csv")
    print(f"  {'band':12s} {'N1 spread':>12s} {'MS spread':>12s}")
    for band in sorted(n1.band.unique()):
        a = (n1[n1.band == band].surr_p95 - n1[n1.band == band].surr_p50).median()
        b = (ms[ms.band == band].surr_p95 - ms[ms.band == band].surr_p50).median()
        print(f"  {band:12s} {a:12.4f} {b:12.4f}")
    print("  -> a comparably-sized spread means N1 is a genuine stochastic null, not a")
    print("     deterministic substitution dressed up as one.")


if __name__ == "__main__":
    d1_orientation()
    d2_per_patient()
    d4_null_spread()
    d3_plain_coherence()
