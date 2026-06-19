#!/usr/bin/env python3
"""audit_113b — STAGE 2 of the duration control: matched-strength null at matched length.

audit_113 (stage 1) showed the β cingulate inference hotspot does not MOVE when
`task_test` is truncated to `task_learn`'s length (8/8 sign-preserved, matched
indistinguishable from full). That is observed-statistic robustness. This stage
closes the question with the proper null: regenerate the matched-strength surrogate
ON THE TRUNCATED `task_test` and recompute the per-system matched-strength p + BH q,
so the verdict is "the hotspot still beats its strength-matched null after
length-matching", not just "the point estimate didn't move".

5-point critical preamble
-------------------------
1. CLAIM. After removing the `task_test`-vs-`task_learn` data-amount asymmetry, the
   β inference-specific component still over-accumulates in the CINGULATE above a
   strength-preserving null (BH q<0.05 across the 9 a-priori systems within β).
2. NULL. Matched-strength (4-cycle ±δ, swap 20, canonical seed) surrogates of the
   TRUNCATED `task_test` FC, combined with the cached full-length surrogates of the
   four unchanged phases (rest_pre_A/B, task_learn, rest_post) — i.e. the EXACT
   audit_110 null with ONLY the `task_test` ensemble swapped for its length-matched
   version. Independent per-phase rewiring (audit_63 convention).
3. STRONGEST ALTERNATIVE. The full-length cingulate q≈0.05 was carried by
   `task_test`'s extra data sharpening `D_TT`; at matched length the surrogate tail
   widens and the observed no longer clears BH.
4. DOES IT ADDRESS IT (mechanism). YES: observed AND surrogate are both built on the
   truncated `task_test`, so the data-amount advantage is removed from BOTH sides of
   the comparison — the only thing left is whether the truncated observed concordance
   still sits in the tail of its truncated-matched-strength null. Built-in checks:
   the observed matched-length cingulate value reproduces audit_113's head-placement
   value; the canonical full-length surrogate for the 4 unchanged phases is reused
   verbatim (cache hit). CANNOT do: rescue power lost to truncation (fewer Welch
   segments → noisier observed AND surrogate); test placements other than head
   (stage 1 already showed placement-invariance).
5. FALSIFICATION / LIMITS. If the cingulate matched-length p fails BH (q≥0.05), the
   honest verdict becomes "the inference hotspot is directionally preserved (stage 1,
   8/8 sign) but attenuated below the multiple-comparison threshold once the duration
   advantage is removed" — i.e. real but duration-assisted, NOT a pure artifact.
   Truncation lowers power on both sides, so a modest q-rise is expected and is a
   power statement, not necessarily a confound statement.

Output: data/audit/inference_localization/lenmatched_null_R{R}_{include,exclude}.csv
"""
from __future__ import annotations

import argparse
import gc
import sys

import numpy as np
import pandas as pd

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS, FS_OVERRIDES, nperseg_for_fs
from lrg_eegfc.config.paths import SEEG_DATAPATH
from lrg_eegfc.utils.fc.msc.msc import compute_msc_welch
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.io.regions import load_channel_regions
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.metrics.node_localization import epi_keep_mask
from lrg_eegfc.utils.surrogate import cophenetic_condensed_from_eigs
from lrg_eegfc.utils.surrogate.matched_strength import (
    load_or_compute_eigs_at_path, surrogate_cache_path,
)

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_110_inference_mark_localization import (  # type: ignore
    COHORT, PHASES5, SEED, SWAP, TARGETS, _targets_from_D,
    load_surrogate_cophenet5,
)
from audit_63_split_baseline_surrogate import (  # type: ignore
    ensure_half_fcs, load_phase_fc,
)
from audit_83_localization_matched_strength import (  # type: ignore
    DROP, _canon_cophenet, unit_means_from_s,
)
from audit_113_inference_length_control import _window  # type: ignore

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "hypothesis_tests"))
from _fc_split_half import _band_abs  # type: ignore

OUT = ROOT / "data/audit/inference_localization"
SYSTEMS_TARGETS = ("standard", "encoding", "inference_pe")
KEY = {"standard": "OFC", "encoding": "OFC", "inference_pe": "cingulate"}
PLACEMENT = "head"   # stage-1 showed placement-invariance; one window for the null


def _tt_fc_lenmatched(pat: str, band: str) -> np.ndarray:
    """Truncated (length-matched to task_learn) task_test imcoh_abs FC, cleaned
    exactly as load_phase_fc would (fill-diag 0, clip[0,1], symmetrize)."""
    X = np.asarray(load_timeseries(pat, "task_test", SEEG_DATAPATH), float)
    if X.shape[0] > X.shape[1]:
        X = X.T
    n_TL = int(max(load_timeseries(pat, "task_learn", SEEG_DATAPATH).shape))
    fs = float(FS_OVERRIDES.get(pat, 2048.0))
    Xw = np.ascontiguousarray(_window(X, n_TL, PLACEMENT))
    freqs, Coh = compute_msc_welch(Xw, fs, nperseg=nperseg_for_fs(fs),
                                   metric="imcoh")
    flo, fhi = BRAIN_BANDS[band]
    W = np.asarray(_band_abs(Coh, freqs, flo, fhi), dtype=float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    W = 0.5 * (W + W.T)
    del Coh, X
    return W


def run_band(band: str, epis, R: int, verbose=True) -> list[dict]:
    # accumulators: [epi][target][system][pat]
    obs_um = {e: {t: {} for t in SYSTEMS_TARGETS} for e in epis}
    surr_um = {e: {t: {} for t in SYSTEMS_TARGETS} for e in epis}

    for pat in COHORT:
        # full-length canonical surrogate cophenetics for all 5 phases
        surco, miss = load_surrogate_cophenet5(pat, band, R)
        if surco is None:
            if verbose:
                print(f"  [skip] {pat}: full-length surrogate missing ({miss}, R={R})")
            continue
        try:
            W_ttm = _tt_fc_lenmatched(pat, band)
        except FileNotFoundError:
            if verbose:
                print(f"  [skip] {pat}: timeseries missing")
            continue

        # observed cophenetics: cached 4 phases + length-matched task_test
        D_obs = {ph: _canon_cophenet(load_phase_fc(pat, ph, band))
                 for ph in PHASES5}
        if W_ttm.shape[0] != int((1 + np.sqrt(1 + 8 * D_obs["task_test"].size)) / 2):
            if verbose:
                print(f"  [skip] {pat}: matched-length N mismatch")
            continue
        D_obs["task_test"] = _canon_cophenet(W_ttm)
        t_obs = _targets_from_D(D_obs)
        N = W_ttm.shape[0]
        iu = np.triu_indices(N, k=1)
        iu_i, iu_j = iu[0].astype(int), iu[1].astype(int)

        # override task_test surrogate with the LENGTH-MATCHED matched-strength
        # ensemble (non-canonical cache path so the canonical cache is untouched)
        path = surrogate_cache_path(pat, band, f"task_test_lm_{PLACEMENT}", R,
                                    SWAP, SEED, "imcoh_abs")
        rng = np.random.default_rng(SEED)
        evals, evecs = load_or_compute_eigs_at_path(path, W_ttm, R, SWAP, rng,
                                                    verbose=verbose)
        surco["task_test"] = [
            None if not np.all(np.isfinite(evals[r]))
            else cophenetic_condensed_from_eigs(evals[r], evecs[r])
            for r in range(R)
        ]

        # per-realization surrogate targets (shared across epi modes)
        sr = []
        for r in range(R):
            Dr = {ph: surco[ph][r] for ph in PHASES5}
            sr.append(None if any(Dr[ph] is None for ph in PHASES5)
                      else _targets_from_D(Dr))

        rdf = load_channel_regions(pat)
        uv = rdf["system"].to_numpy().astype(str)

        for epi_excl in epis:
            keep = epi_keep_mask(rdf, pat) if epi_excl else None
            for t in SYSTEMS_TARGETS:
                om = unit_means_from_s(t_obs[t], iu_i, iu_j, uv, keep, DROP["system"])
                for u, m in om.items():
                    obs_um[epi_excl][t].setdefault(u, {})[pat] = m
                per = {u: [] for u in om}
                for r in range(R):
                    if sr[r] is None:
                        continue
                    sm = unit_means_from_s(sr[r][t], iu_i, iu_j, uv, keep,
                                           DROP["system"])
                    for u in per:
                        if u in sm:
                            per[u].append(sm[u])
                for u in om:
                    surr_um[epi_excl][t].setdefault(u, {})[pat] = np.array(per[u])
        del sr, evals, evecs, surco
        gc.collect()
        if verbose:
            print(f"  [{pat}] matched-length null done (N={N}, R={R})")

    # cohort verdict per epi × target × system
    rows = []
    for epi_excl in epis:
        etag = "exclude" if epi_excl else "include"
        for t in SYSTEMS_TARGETS:
            sys_rows = []
            for u, pm in obs_um[epi_excl][t].items():
                samplers = sorted(pm)
                K = len(samplers)
                if K == 0:
                    continue
                M_obs = float(np.median([pm[p] for p in samplers]))
                arrs = [surr_um[epi_excl][t][u].get(p) for p in samplers]
                arrs = [a for a in arrs if a is not None and len(a) > 0]
                if not arrs:
                    continue
                Lmin = min(len(a) for a in arrs)
                Rmat = np.vstack([a[:Lmin] for a in arrs])
                M_surr = np.median(Rmat, axis=0)
                nval = M_surr.size
                p_ms = (1 + int(np.sum(M_surr >= M_obs))) / (nval + 1)
                row = {"band": band, "epi": etag, "target": t, "unit": u,
                       "K_implanted": K, "M_obs": M_obs,
                       "surr_median": float(np.median(M_surr)),
                       "matched_strength_p": p_ms, "n_surr": nval,
                       "bh_q_system": np.nan}
                rows.append(row)
                sys_rows.append(row)
            if sys_rows:
                qs = bh_fdr([r["matched_strength_p"] for r in sys_rows])
                for r, q in zip(sys_rows, qs):
                    r["bh_q_system"] = q
            if verbose:
                print(f"\n  == {band} {t} epi-{etag} (matched length, R={R}) ==")
                for r in sorted(sys_rows, key=lambda r: r["matched_strength_p"]):
                    star = "*" if r["bh_q_system"] < 0.05 else " "
                    tag = " <<< KEY" if r["unit"] == KEY[t] else ""
                    print(f"    {r['unit']:16s} M_obs={r['M_obs']:+.3e} "
                          f"p={r['matched_strength_p']:.4f} q={r['bh_q_system']:.3f}"
                          f"{star} (K={r['K_implanted']}){tag}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--band", default="beta")
    ap.add_argument("--R", type=int, default=200)
    ap.add_argument("--epi-mode", choices=["include", "exclude"], default=None)
    args = ap.parse_args()
    epis = ([args.epi_mode == "exclude"] if args.epi_mode else [False, True])

    OUT.mkdir(parents=True, exist_ok=True)
    print(f"[audit_113b] pre-flight: ensure rsPre half FCs cached ({args.band})")
    for pat in COHORT:
        ensure_half_fcs(pat, [args.band])

    print(f"\n#### matched-strength null at MATCHED LENGTH : {args.band} "
          f"R={args.R} placement={PLACEMENT} ####")
    rows = run_band(args.band, epis, args.R)
    df = pd.DataFrame(rows)
    for epi_excl in epis:
        etag = "exclude" if epi_excl else "include"
        sub = df[df.epi == etag]
        out = OUT / f"lenmatched_null_R{args.R}_{etag}.csv"
        sub.to_csv(out, index=False)
        print(f"  -> {out} ({len(sub)} rows)")


if __name__ == "__main__":
    main()
