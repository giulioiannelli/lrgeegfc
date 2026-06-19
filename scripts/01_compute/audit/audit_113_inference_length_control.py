#!/usr/bin/env python3
"""audit_113 — length-equalized (duration) control for the β inference-mark.

The inference-specific per-pair vector is ``f = D_taskTest − D_taskLearn``.
``task_test`` is 1.35–2.53× LONGER than ``task_learn`` in ALL 10 patients
(median 1.57×), so ``D_TT`` is estimated from ~1.5–2.5× more Welch segments than
``D_TL`` — a systematically SHARPER (less-noisy) cophenetic. This control asks
whether the cingulate inference hotspot (audit_110, R=1000 q=0.050) and the OFC
encoding/standard hotspots are an artifact of that data-amount asymmetry rather
than of task content.

5-point critical preamble
-------------------------
1. CLAIM. The β localization is content/inference-specific: encoding & standard
   over-accumulate in OFC; the inference-specific component (``f``, controlling
   encoding) over-accumulates in the CINGULATE.
2. NULL / CONFOUND. ``task_test`` has more data than ``task_learn`` → ``D_TT``
   less noisy than ``D_TL`` → ``f = D_TT − D_TL`` carries a duration-driven
   "sharpening" direction. If that direction rank-aligns with the persistent
   rest_post trace ``p`` on cingulate-incident edges, the cingulate inference
   hotspot is a duration artifact, not inference.
3. STRONGEST ALTERNATIVE THE CONTROL MUST ADDRESS. "The cingulate inference
   hotspot is driven by ``task_test`` simply having ~1.6× more samples than
   ``task_learn``."
4. DOES IT ADDRESS IT (by mechanism). YES for the data-amount channel: recompute
   ``D_TT`` from a ``task_test`` segment TRUNCATED to ``task_learn``'s sample
   count, so ``D_TT'`` and ``D_TL`` use the same number of Welch segments
   (matched estimation precision), then rerun the EXACT audit_110 localization
   with ONLY ``D_task_test`` swapped. Two built-in checks make it airtight:
   (a) the ENCODING target ``e = D_TL − D_preA`` has NO task_test dependence, so
   its localization MUST be bit-identical full-vs-matched (any change ⇒ bug);
   (b) a full-length recompute of ``D_TT`` must reproduce the cached ``D_TT``
   (pipeline == canonical cache). What it CANNOT do: regenerate the
   matched-strength surrogate on truncated data (this is an observed-statistic
   robustness stage, NOT a new null); exclude genuine learn/test CONTENT
   differences (those are the signal); fully rule out within-test
   non-stationarity (mitigated by 3 window placements head/center/tail).
5. FALSIFICATION / LIMITS. If the per-patient demeaned cingulate inference
   unit-mean flips sign or collapses toward 0 under length-matching across
   patients and placements, the localization is duration-confounded → withdraw /
   heavily caveat. DIRECTION-OF-BIAS NOTE: truncation DISCARDS data, so ``D_TT'``
   is NOISIER than the cached ``D_TT`` — this biases the concordance toward ZERO.
   Survival under matched (noisier) length is therefore CONSERVATIVE evidence for
   content-specificity.

Output: data/audit/inference_localization/length_control.csv
"""
from __future__ import annotations

import argparse
import gc
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr  # noqa: F401  (parity import)
from lrg_eegfc.utils.metrics.node_localization import epi_keep_mask

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_110_inference_mark_localization import (  # type: ignore
    COHORT, PHASES5, _targets_from_D,
)
from audit_63_split_baseline_surrogate import (  # type: ignore
    ensure_half_fcs, load_phase_fc,
)
from audit_83_localization_matched_strength import (  # type: ignore
    DROP, _canon_cophenet, unit_means_from_s,
)

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "hypothesis_tests"))
from _fc_split_half import _band_abs  # type: ignore  # same band-abs as the rsPre halves

OUT = ROOT / "data/audit/inference_localization"
BANDS = ["beta", "alpha", "low_gamma"]
TARGETS = ("standard", "encoding", "inference_pe")
PLACEMENTS = ("head", "center", "tail")
# the system each target localizes to (audit_110); we track that cell full-vs-matched
KEY = {"standard": "OFC", "encoding": "OFC", "inference_pe": "cingulate"}
# full-length recompute validation (proves pipeline == cache); a few spanning patients
VALIDATE_PATS = {"Pat_02", "Pat_06", "Pat_15"}


# --------------------------------------------------------------------------- #
# recompute task_test cophenetic from a length-matched window (all bands, 1 CSD)
# --------------------------------------------------------------------------- #
def _window(X: np.ndarray, n_keep: int | None, placement: str) -> np.ndarray:
    T = X.shape[1]
    if n_keep is None or n_keep >= T:
        return X
    if placement == "head":
        return X[:, :n_keep]
    if placement == "tail":
        return X[:, T - n_keep:]
    s = (T - n_keep) // 2          # center
    return X[:, s:s + n_keep]


def _tt_cophenets(X: np.ndarray, fs: float, n_keep: int | None,
                  placement: str) -> dict:
    """task_test canonical cophenetic per band, from ONE imcoh CSD of the window.

    Reproduces the canonical imcoh_abs → load_phase_fc cleanup → ``_canon_cophenet``
    path exactly (band-abs via ``_band_abs``, fill-diag 0, clip[0,1], symmetrize)."""
    Xw = np.ascontiguousarray(_window(X, n_keep, placement))
    freqs, Coh = compute_msc_welch(Xw, fs, nperseg=nperseg_for_fs(fs),
                                   metric="imcoh")
    out = {}
    for band in BANDS:
        flo, fhi = BRAIN_BANDS[band]
        W = np.asarray(_band_abs(Coh, freqs, flo, fhi), dtype=float)
        np.fill_diagonal(W, 0.0)
        W = np.clip(W, 0.0, 1.0)
        W = 0.5 * (W + W.T)
        out[band] = _canon_cophenet(W)
    del Coh
    return out


# --------------------------------------------------------------------------- #
# per-patient
# --------------------------------------------------------------------------- #
def run_patient(pat: str, epis, placements, verbose=True) -> list[dict]:
    """Recompute the expensive task_test CSDs ONCE, evaluate both epi modes off
    them (cophenetics are epi-independent; only the unit-mean keep-mask differs)."""
    rows: list[dict] = []
    # cached per-band cophenetics for the 5 phases (full)
    try:
        D_full = {b: {ph: _canon_cophenet(load_phase_fc(pat, ph, b))
                      for ph in PHASES5} for b in BANDS}
    except FileNotFoundError:
        if verbose:
            print(f"  [skip] {pat}: cached phase FC missing")
        return rows

    rdf = load_channel_regions(pat)
    uv = rdf["system"].to_numpy().astype(str)

    # geometry + full-length per-pair targets per band (epi-independent)
    geom, t_full_by_b = {}, {}
    for b in BANDS:
        N = int((1 + np.sqrt(1 + 8 * D_full[b]["task_test"].size)) / 2)
        iu = np.triu_indices(N, k=1)
        geom[b] = (iu[0].astype(int), iu[1].astype(int))
        t_full_by_b[b] = _targets_from_D(D_full[b])

    # timeseries + matched length
    X = np.asarray(load_timeseries(pat, "task_test", SEEG_DATAPATH), float)
    if X.shape[0] > X.shape[1]:
        X = X.T
    Xl = load_timeseries(pat, "task_learn", SEEG_DATAPATH)
    n_TL, n_TT = int(max(Xl.shape)), int(X.shape[1])
    del Xl
    fs = float(FS_OVERRIDES.get(pat, 2048.0))

    # validation: full-length recompute reproduces cached D_TT (pipeline == cache)
    val = {b: (np.nan, np.nan) for b in BANDS}
    if pat in VALIDATE_PATS:
        ttc_full = _tt_cophenets(X, fs, None, "head")
        for b in BANDS:
            if ttc_full[b].size == D_full[b]["task_test"].size:
                val[b] = (float(spearmanr(ttc_full[b],
                                          D_full[b]["task_test"]).correlation),
                          float(np.max(np.abs(ttc_full[b]
                                              - D_full[b]["task_test"]))))
        if verbose:
            vb = val["beta"]
            print(f"  [{pat}] full-recompute vs cache (β): rho={vb[0]:.6f} "
                  f"max|Δ|={vb[1]:.2e}  (ratio TT/TL={n_TT / n_TL:.2f})")

    # matched-length windows: recompute task_test cophenetics + targets ONCE
    t_m_by_place = {}
    for placement in placements:
        ttc = _tt_cophenets(X, fs, n_TL, placement)
        t_m = {}
        for b in BANDS:
            if ttc[b].size != D_full[b]["task_test"].size:
                if verbose:
                    print(f"  [warn] {pat}/{b}/{placement}: N mismatch")
                continue
            D_m = dict(D_full[b])
            D_m["task_test"] = ttc[b]
            t_m[b] = _targets_from_D(D_m)
        t_m_by_place[placement] = t_m
    del X
    gc.collect()

    # evaluate per epi mode (cheap unit-mean aggregation off shared targets)
    for epi_excl in epis:
        keep = epi_keep_mask(rdf, pat) if epi_excl else None
        etag = "exclude" if epi_excl else "include"
        full_um = {b: {t: unit_means_from_s(t_full_by_b[b][t], *geom[b], uv,
                                            keep, DROP["system"])
                       for t in TARGETS} for b in BANDS}
        for placement, t_m in t_m_by_place.items():
            for b in t_m:
                for t in TARGETS:
                    u = KEY[t]
                    mm = unit_means_from_s(t_m[b][t], *geom[b], uv, keep,
                                           DROP["system"])
                    rows.append({
                        "band": b, "epi": etag, "pat": pat,
                        "placement": placement, "target": t, "system": u,
                        "full": float(full_um[b][t].get(u, np.nan)),
                        "matched": float(mm.get(u, np.nan)),
                        "n_TL": n_TL, "n_TT": n_TT, "ratio": n_TT / n_TL,
                        "val_corr_beta": val["beta"][0],
                    })
    return rows


# --------------------------------------------------------------------------- #
# cohort summary
# --------------------------------------------------------------------------- #
def summarize(df: pd.DataFrame) -> None:
    for b in BANDS:
        print(f"\n==== length-control cohort summary : {b} "
              f"(matched = median over {len(PLACEMENTS)} placements) ====")
        print("  target/system    epi      med_full    med_match   ratio  "
              "pos(m)  sign-agree  W(m>0)  W(m≠f)  rho_fm")
        for t in TARGETS:
            for epi in sorted(df.epi.unique()):
                sub = df[(df.band == b) & (df.target == t) & (df.epi == epi)]
                if sub.empty:
                    continue
                g = (sub.groupby("pat")
                     .agg(full=("full", "first"), matched=("matched", "median"))
                     .dropna())
                if len(g) < 3:
                    continue
                f_, m_ = g.full.to_numpy(), g.matched.to_numpy()
                med_f, med_m = np.median(f_), np.median(m_)
                ratio = med_m / med_f if med_f != 0 else np.nan
                n_pos_m = int((m_ > 0).sum())
                n_sign = int((np.sign(f_) == np.sign(m_)).sum())
                try:
                    w_pos = wilcoxon(m_, alternative="greater").pvalue
                except ValueError:
                    w_pos = np.nan
                try:
                    w_diff = wilcoxon(m_, f_).pvalue
                except ValueError:
                    w_diff = np.nan
                rho = (float(spearmanr(f_, m_).correlation)
                       if len(f_) > 2 else np.nan)
                flag = ""
                if t == "encoding":  # MUST be identical (no task_test dependence)
                    flag = "  [enc⊥TT: full==matched expected]" \
                        if np.max(np.abs(f_ - m_)) < 1e-9 else "  [!! ENC CHANGED — BUG]"
                print(f"  {t:13s} {KEY[t]:9s} {epi:7s}: "
                      f"{med_f:+.3e} {med_m:+.3e}  {ratio:5.2f}  "
                      f"{n_pos_m:2d}/{len(m_):<2d}   {n_sign:2d}/{len(m_):<2d}     "
                      f"{w_pos:.3f}  {w_diff:.3f}  {rho:+.2f}{flag}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--epi-mode", choices=["include", "exclude"], default=None)
    ap.add_argument("--placements", nargs="+", default=list(PLACEMENTS))
    ap.add_argument("--patients", nargs="+", default=COHORT)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    print("[audit_113] pre-flight: ensure rsPre half FCs cached")
    for pat in args.patients:
        ensure_half_fcs(pat, BANDS)

    epis = ([args.epi_mode == "exclude"] if args.epi_mode else [False, True])
    print(f"\n#### length-equalized control : epi={[('exclude' if e else 'include') for e in epis]} "
          f"placements={args.placements} ####")
    allrows: list[dict] = []
    for pat in args.patients:
        allrows += run_patient(pat, epis, args.placements)

    df = pd.DataFrame(allrows)
    out = OUT / "length_control.csv"
    df.to_csv(out, index=False)
    summarize(df)
    print(f"\n  -> {out} ({len(df)} rows)")


if __name__ == "__main__":
    main()
