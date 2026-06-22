#!/usr/bin/env python3
"""audit_103b — STAGE-2-ARC: duration control on the WHOLE-BRAIN consolidation arc.

The arc band-dissociation (audit_103: alpha consolidates ENCODING only; beta
consolidates ENCODING + an INFERENCE-SPECIFIC reorganization) is built on the same
inference-specific functional ``f = D_taskTest - D_taskLearn`` whose LOCALIZATION
just failed the duration control (audit_113b: inference->cingulate downgraded once
`task_test` is length-matched to `task_learn`). Because `task_test` is 1.35-2.53x
LONGER than `task_learn` in every patient, `D_taskTest` is built from ~1.6x more
Welch segments than `D_taskLearn`, so `f` could inflate from data-amount alone. This
script closes that gap for the cohort-level arc claim: it re-runs the EXACT audit_103
matched-strength null with ONLY the `task_test` view (observed cophenetic AND its
matched-strength surrogate ensemble) swapped for the length-matched version, then asks
whether the beta inference-specific persistence ``T_infspec_pe = partial(f, p | e)``
still beats its strength-matched null after the duration advantage is removed.

5-point critical preamble
-------------------------
1. CLAIM. After length-matching `task_test` to `task_learn`, the beta rhythm STILL
   consolidates an inference-specific reorganization beyond encoding: cohort
   ``T_infspec_pe`` clears the matched-strength null in beta (one-sided paired
   Wilcoxon obs > surrogate-median p<0.05), and the alpha-vs-beta dissociation
   (alpha = encoding-only) survives.
2. NULL. The audit_103 matched-strength null (independent per-phase 4-cycle +/-delta
   rewiring, swap 20, canonical seed) with the `task_test` ensemble REPLACED by its
   length-matched ensemble (the very surrogates audit_113b generated on the truncated
   `task_test`); the four unchanged phases (rest_pre_A/B, task_learn, rest_post) keep
   their canonical full-length surrogates (cache hits).
3. STRONGEST ALTERNATIVE. The full-length beta ``T_infspec_pe`` (p=0.0068) was carried
   by `task_test`'s extra data sharpening `D_taskTest`; at matched length both the
   observed `f` and its surrogate widen, and the cohort no longer clears the null.
4. DOES IT ADDRESS IT (mechanism). YES: observed AND surrogate are both rebuilt on the
   truncated `task_test`, so the data-amount advantage is removed from BOTH sides of
   the inference functional. Built-in gates: (a) the truncated `task_test` surrogate
   path is the SAME cache audit_113b wrote (reuse, not regenerate); (b) the ENCODING
   functional ``T_learn = rho(e, p)`` has NO `task_test` dependence -> the matched-length
   observed must reproduce the canonical audit_103 ``T_learn`` BIT-FOR-BIT (swap-logic
   bug-check); (c) the standard trace ``T_test`` stays a live positive control. CANNOT
   do: rescue power lost to truncation (fewer Welch segments -> noisier on both sides);
   test window placements other than head (audit_113 showed placement-invariance).
5. FALSIFICATION / LIMITS. If beta ``T_infspec_pe`` fails the matched-length null, the
   arc's inference-specific headline is duration-assisted and must be downgraded to an
   encoding-consolidation result (exactly as the localization was). A modest p-rise from
   truncation alone is a power statement; a rise past 0.05 with the standard trace still
   clearing is a confound statement. Encoding (T_learn) is the duration-invariant anchor.

Output
------
data/audit/consolidation_arc/arc_lenmatched_null_R{R}_{band}.csv
+ printed cohort matched-strength verdict (canonical-vs-matched obs side by side).
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.surrogate.matched_strength import surrogate_cache_path

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
# Reuse the EXACT audit_103 arc machinery (no fork): loaders, ultrametric, the
# matched-strength stack builder, the functional definitions, the per-cell RNG,
# the cohort verdict, the cohort list and output dir.
from audit_103_cophenetic_consolidation_arc import (  # type: ignore
    ARC_PHASES, BANDS, COHORT, FUNCTIONALS, OUT,
    _arc_cell_rng, _obs_functionals, _surr_functionals,
    ensure_half_fcs, load_phase_fc, lrg_ultrametric_condensed, ws,
)
from audit_83_wm_stratified_cophenetic import _surr_stacks  # type: ignore
from audit_113_inference_length_control import _window  # type: ignore

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "hypothesis_tests"))
from _fc_split_half import _band_abs  # type: ignore

SEED_CANON = ws.config_seed("full")
ARC_CSV = OUT / "arc_per_patient.csv"   # canonical (full-length) observed, for the gate


def _tt_fc_lm(pat: str, band: str, placement: str) -> np.ndarray:
    """Truncated (length-matched to task_learn) task_test imcoh_abs FC at a given
    window placement, cleaned exactly as load_phase_fc would. Mirrors audit_113b's
    `_tt_fc_lenmatched` but with the placement parameterized."""
    X = np.asarray(load_timeseries(pat, "task_test", SEEG_DATAPATH), float)
    if X.shape[0] > X.shape[1]:
        X = X.T
    n_TL = int(max(load_timeseries(pat, "task_learn", SEEG_DATAPATH).shape))
    fs = float(FS_OVERRIDES.get(pat, 2048.0))
    Xw = np.ascontiguousarray(_window(X, n_TL, placement))
    freqs, Coh = compute_msc_welch(Xw, fs, nperseg=nperseg_for_fs(fs), metric="imcoh")
    flo, fhi = BRAIN_BANDS[band]
    W = np.asarray(_band_abs(Coh, freqs, flo, fhi), dtype=float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    W = 0.5 * (W + W.T)
    del Coh, X
    return W


def _surr_stack_for(pat: str, band: str, phase: str, W: np.ndarray, R: int,
                    placement: str) -> tuple[np.ndarray, bool]:
    """Cophenetic surrogate stack (R, n_pairs) for one phase. For `task_test` use the
    length-matched ensemble at the placement's non-canonical cache path (cache HIT for
    beta-head from audit_113b); every other phase uses the canonical full-length
    ensemble (cache HIT)."""
    sf = ws.SWAP_FACTOR
    if phase == "task_test":
        path = surrogate_cache_path(pat, band, f"task_test_lm_{placement}", R, sf,
                                    SEED_CANON, "imcoh_abs")
        cache_hit = path.exists()
        coph, _ = _surr_stacks(W, path, np.random.default_rng(SEED_CANON), R, sf, False)
        return coph, cache_hit
    coph, _ = _surr_stacks(W, ws.surr_eig_path("full", pat, band, phase, R, sf),
                           _arc_cell_rng(pat, band, phase), R, sf, False)
    return coph, True


def run_band(band: str, R: int, placement: str = "head",
             verbose: bool = True) -> pd.DataFrame:
    canon = None
    if ARC_CSV.exists():
        c = pd.read_csv(ARC_CSV)
        canon = c[c["band"] == band].set_index("patient")
    rows = []
    for pat in COHORT:
        ensure_half_fcs(pat, [band])
        try:
            # FCs: canonical for 4 phases, truncated for task_test.
            Ws = {ph: load_phase_fc(pat, ph, band) for ph in ARC_PHASES}
            W_ttm = _tt_fc_lm(pat, band, placement)
            if W_ttm.shape[0] != Ws["task_test"].shape[0]:
                print(f"  [skip] {pat}: matched-length N mismatch")
                continue
            Ws["task_test"] = W_ttm

            # Observed: 4 canonical cophenetics + length-matched task_test cophenetic.
            D = {ph: lrg_ultrametric_condensed(Ws[ph]) for ph in ARC_PHASES}
            obs = _obs_functionals(D)

            # Surrogate stacks: canonical 4 phases + length-matched task_test ensemble.
            sc, tt_hit = {}, True
            for ph in ARC_PHASES:
                sc[ph], hit = _surr_stack_for(pat, band, ph, Ws[ph], R, placement)
                if ph == "task_test":
                    tt_hit = hit
            surr = _surr_functionals(sc, R)

            row = {"patient": pat, "band": band, "tt_surr_cache_hit": tt_hit}
            for k in FUNCTIONALS:
                s = surr[k][np.isfinite(surr[k])]
                row[f"{k}_obs"] = obs[k]
                row[f"{k}_surr_p50"] = float(np.median(s)) if s.size else np.nan
                row[f"{k}_p"] = float(np.mean(s >= obs[k])) if s.size else np.nan
                row[f"{k}_nsurr"] = int(s.size)
            # Gate (b): encoding T_learn must be bit-identical to canonical (no TT dep).
            if canon is not None and pat in canon.index:
                row["T_learn_canon"] = float(canon.loc[pat, "T_learn"])
                row["T_infspec_pe_canon"] = float(canon.loc[pat, "T_infspec_pe"])
            rows.append(row)
            if verbose:
                dl = (row.get("T_learn_canon", np.nan) - obs["T_learn"])
                print(f"  [{pat}] tt_hit={tt_hit}  "
                      f"T_learn={obs['T_learn']:+.3f}(dCanon={dl:+.0e})  "
                      f"T_infspec_pe={obs['T_infspec_pe']:+.3f}"
                      f"(p={row['T_infspec_pe_p']:.3f})  "
                      f"T_test={obs['T_test']:+.3f}(p={row['T_test_p']:.3f})")
        except FileNotFoundError as exc:
            print(f"  [skip] {pat}: {exc}")
        except Exception as exc:  # noqa: BLE001
            print(f"  [audit_103b] FAIL {pat} {band}: {type(exc).__name__}: {exc}")
    return pd.DataFrame(rows)


def _verdict(df: pd.DataFrame, band: str) -> None:
    print(f"\n=== STAGE-2-ARC MATCHED-LENGTH VERDICT  band={band}  "
          f"(positive = trace; obs vs surr-median) ===")
    # Gate reports
    if "T_learn_canon" in df:
        dmax = float(np.nanmax(np.abs(df["T_learn_obs"] - df["T_learn_canon"])))
        gate = "PASS (encoding bit-identical, swap-logic correct)" if dmax < 1e-9 \
            else f"FAIL (max|d|={dmax:.2e} -- task_test leaked into encoding!)"
        print(f"[gate b] T_learn matched vs canonical: max|d|={dmax:.2e} -> {gate}")
    nhit = int(df["tt_surr_cache_hit"].sum()) if "tt_surr_cache_hit" in df else 0
    print(f"[gate a] truncated task_test surrogate cache hits (reuse audit_113b): "
          f"{nhit}/{len(df)}")

    print(f"\n{'functional':<16}{'obs_med':>9}{'surr_med':>9}{'n<.05':>7}"
          f"{'wilcox_p':>10}{'LOp15_p':>9}{'canon_obs':>10}")
    for k in FUNCTIONALS:
        d = df.dropna(subset=[f"{k}_obs", f"{k}_surr_p50"])
        if d.empty:
            continue
        n_above = int((d[f"{k}_p"] < 0.05).sum())
        v = ws.cohort_verdict(d[f"{k}_obs"].values, d[f"{k}_surr_p50"].values, n_above)
        d2 = d[d["patient"] != "Pat_15"]
        v2 = ws.cohort_verdict(d2[f"{k}_obs"].values, d2[f"{k}_surr_p50"].values,
                               int((d2[f"{k}_p"] < 0.05).sum()))
        canon_obs = (float(np.median(df[f"{k}_canon"])) if f"{k}_canon" in df
                     else np.nan)
        print(f"{k:<16}{v['med_obs']:>+9.3f}{v['med_surr']:>+9.3f}"
              f"{n_above:>4}/{v['n_defined']:<2}{v['wilcoxon_p']:>10.4f}"
              f"{v2['wilcoxon_p']:>9.4f}{canon_obs:>+10.3f}")
    print("\n[read] T_infspec_pe = inference-specific persistence controlling "
          "encoding -- the claim under duration test. T_learn = encoding anchor "
          "(must be duration-invariant). T_test = standard trace (live positive "
          "control). Compare wilcox_p to the full-length audit_103 null.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", default="beta")
    ap.add_argument("--R", type=int, default=200)
    ap.add_argument("--placement", default="head", choices=["head", "center", "tail"])
    args = ap.parse_args()
    print(f"#### STAGE-2-ARC duration control : band={args.band} R={args.R} "
          f"placement={args.placement} ####\n")
    df = run_band(args.band, args.R, args.placement)
    suffix = "" if args.placement == "head" else f"_{args.placement}"
    out = OUT / f"arc_lenmatched_null_R{args.R}_{args.band}{suffix}.csv"
    df.to_csv(out, index=False)
    print(f"\n[audit_103b] wrote {out} ({len(df)} rows)")
    _verdict(df, args.band)


if __name__ == "__main__":
    main()
