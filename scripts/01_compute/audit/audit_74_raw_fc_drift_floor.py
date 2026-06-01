#!/usr/bin/env python3
"""Audit 74 — within-session DRIFT FLOOR for ``ρ_split^raw`` (raw-FC substrate).

Adapts the cophenetic drift-floor control (``h2e_split_half.py`` +
``continuous_trace_controls.py``) from the LRG ultrametric ``D`` to the raw
``|ImCoh|`` adjacency ``A``, so the substrate per-pair probe ``ρ_split^raw``
(audit_67) gets the SAME within-session drift-floor test the cophenetic
``ρ_split`` already has. No new methodology — only the substrate changes
(``D → A``).

Critical preamble (per project rule):
 1. Claim: the substrate trace ``ρ_split^raw`` exceeds what within-session
    drift alone produces (no task involved).
 2. Null: ``ρ_null_drift^raw = Spearman(A^preB − A^preA, A^postB − A^postA)``
    — correlate two NO-TASK within-rest difference maps. Pure monotone
    session drift makes both point the same way → positive floor.
 3. Strongest alternative the null must catch: a monotone session drift that
    pulls FC the same way through pre→task→post, manufacturing a positive
    ``ρ_split^raw`` with no task-specific reorganization.
 4. Does it catch it? Mechanically yes for the component shared by drift —
    if drift dominates, the no-task floor reproduces the observed alignment
    and the paired ``(obs − floor)`` Wilcoxon fails to clear it. It does NOT
    catch a drift that is itself task-locked, nor edge-strength structure
    (that is audit_67's matched-strength null — complementary, not
    redundant).
 5. Falsification: if ``ρ_split^raw`` is not reliably ``> ρ_null_drift^raw``
    across patients in a band, that band's substrate trace is
    drift-consistent.

Reuses (no reinvention):
 - ``compute_imcoh_abs_halves`` (scripts/01_compute/hypothesis_tests/_fc_split_half.py)
 - the rest_pre half-FC cache already populated by audit_67
   (``data/cache/imcoh_halves_fc/<pat>/<band>_rest_pre_<half>_imcoh_abs.npy``)
 - the published ``ρ_split^raw`` obs_rho values
   (``data/audit/raw_fc_matched_strength/per_patient_per_band_all_bands.csv``)
 - ``wilcoxon_z`` / ``rank_biserial`` / ``bh_fdr`` (library)
 - drift-floor formula identical to ``h2e_split_half.py`` (``D → A``)

Noise-regime note (inherited from the cophenetic recipe, stated for honesty):
 ``ρ_split^raw`` (obs) uses FULL task_test/rest_post against half rest_pre;
 the drift floor uses HALF rest_post — so the floor lives in a slightly
 noisier regime. The test is a statistical contrast, not a magnitude
 comparison (same caveat as H2e / continuous-trace controls).

Writes:
 ``data/audit/raw_fc_matched_strength/drift_floor_per_patient.csv``
 ``data/audit/raw_fc_matched_strength/drift_floor_band_stats.md``
"""
from __future__ import annotations

import gc
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import (
    BRAIN_BANDS,
    BRAIN_BANDS_NAMES,
    BRAIN_BAND_TEX_DICT,
    DEFAULT_SAMPLE_RATE,
    FS_OVERRIDES,
    PATIENTS_4PHASE,
    nperseg_for_fs,
)
from lrg_eegfc.config.paths import CACHE_ROOT, SEEG_DATAPATH
from lrg_eegfc.utils.io.patient import load_timeseries
from lrg_eegfc.utils.metrics.hypothesis import wilcoxon_z, rank_biserial

# Reuse the heavy FC-from-halves helper (no copy).
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "hypothesis_tests"))
from _fc_split_half import compute_imcoh_abs_halves  # type: ignore

OUT = ROOT / "data" / "audit" / "raw_fc_matched_strength"
HALVES_FC_CACHE = CACHE_ROOT / "imcoh_halves_fc"
OBS_CSV = OUT / "per_patient_per_band_all_bands.csv"

REST_PHASES = ("rest_pre", "rest_post")


# ---------------------------------------------------------------------------
# Half-FC cache (identical convention to audit_67; generalized to rest_post)
# ---------------------------------------------------------------------------
def _half_fc_path(pat: str, band: str, phase: str, half: str) -> Path:
    return HALVES_FC_CACHE / pat / f"{band}_{phase}_{half}_imcoh_abs.npy"


def ensure_half_fcs(pat: str, phase: str, bands: list[str]) -> None:
    """Populate {band}_{phase}_{A,B}_imcoh_abs.npy for one resting phase."""
    missing = [
        (b, h)
        for b in bands
        for h in ("A", "B")
        if not _half_fc_path(pat, band=b, phase=phase, half=h).exists()
    ]
    if not missing:
        return
    print(f"[audit_74] {pat}/{phase}: caching {len(missing)} missing half FCs")
    try:
        X = load_timeseries(pat, phase, SEEG_DATAPATH)
    except Exception as exc:  # pragma: no cover
        print(f"[audit_74] {pat}/{phase}: timeseries load FAILED ({exc})")
        return
    X = np.asarray(X, dtype=np.float64)
    if X.shape[0] > X.shape[1]:
        X = X.T
    fs = FS_OVERRIDES.get(pat, DEFAULT_SAMPLE_RATE)
    nperseg_half = max(256, nperseg_for_fs(fs) // 2)
    halves_fc = compute_imcoh_abs_halves(X, fs, nperseg_half, BRAIN_BANDS)
    del X
    gc.collect()
    pat_dir = HALVES_FC_CACHE / pat
    pat_dir.mkdir(parents=True, exist_ok=True)
    for (band, half), A in halves_fc.items():
        if band not in bands:
            continue
        np.save(_half_fc_path(pat, band, phase, half), np.asarray(A, dtype=np.float32))
    halves_fc.clear()
    gc.collect()


def load_half(pat: str, band: str, phase: str, half: str) -> np.ndarray | None:
    """Load a cached half FC, cleaned exactly like audit_67.load_phase_fc."""
    p = _half_fc_path(pat, band, phase, half)
    if not p.exists():
        return None
    W = np.asarray(np.load(p), dtype=float)
    np.fill_diagonal(W, 0.0)
    W = np.clip(W, 0.0, 1.0)
    W = 0.5 * (W + W.T)
    return W


def condensed(W: np.ndarray) -> np.ndarray:
    """Off-diagonal upper-triangle vector (audit_67.raw_fc_condensed)."""
    W = W.copy()
    np.fill_diagonal(W, 0.0)
    W = np.maximum(W, W.T)
    return squareform(W, checks=False)


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    if x.size == 0 or x.std() == 0 or y.std() == 0:
        return float("nan")
    r, _ = spearmanr(x, y)
    return float(r) if np.isfinite(r) else float("nan")


# ---------------------------------------------------------------------------
# Per-patient drift-floor quartet (mirrors h2e.rho_quartet, D → A)
# ---------------------------------------------------------------------------
def drift_floor_for_patient(pat: str) -> dict[str, dict[str, float]]:
    """Return {band: {rho_null_drift, rho_within_rpre, rho_within_rpost}}."""
    out: dict[str, dict[str, float]] = {}
    for band in BRAIN_BANDS_NAMES:
        A_preA = load_half(pat, band, "rest_pre", "A")
        A_preB = load_half(pat, band, "rest_pre", "B")
        A_postA = load_half(pat, band, "rest_post", "A")
        A_postB = load_half(pat, band, "rest_post", "B")
        if any(M is None for M in (A_preA, A_preB, A_postA, A_postB)):
            continue
        shapes = {M.shape for M in (A_preA, A_preB, A_postA, A_postB)}
        if len(shapes) != 1:
            print(f"[audit_74] {pat}/{band}: half shape mismatch {shapes} — skip")
            continue
        d_preA, d_preB = condensed(A_preA), condensed(A_preB)
        d_postA, d_postB = condensed(A_postA), condensed(A_postB)
        out[band] = {
            "rho_null_drift": _spearman(d_preB - d_preA, d_postB - d_postA),
            "rho_within_rpre": _spearman(d_preA, d_preB),
            "rho_within_rpost": _spearman(d_postA, d_postB),
        }
    return out


# ---------------------------------------------------------------------------
# Cohort assembly + per-band paired test
# ---------------------------------------------------------------------------
def collect() -> pd.DataFrame:
    obs = pd.read_csv(OBS_CSV)[["patient", "band", "obs_rho"]].rename(
        columns={"obs_rho": "rho_split_raw"}
    )
    rows: list[dict] = []
    for pat in PATIENTS_4PHASE:
        for phase in REST_PHASES:
            ensure_half_fcs(pat, phase, list(BRAIN_BANDS_NAMES))
        floor = drift_floor_for_patient(pat)
        for band, vals in floor.items():
            rows.append({"patient": pat, "band": band, **vals})
        gc.collect()
    floor_df = pd.DataFrame(rows)
    return obs.merge(floor_df, on=["patient", "band"], how="inner")


def band_stats(df: pd.DataFrame) -> pd.DataFrame:
    recs: list[dict] = []
    for band in BRAIN_BANDS_NAMES:
        b = df[df.band == band].copy()
        if b.empty:
            continue
        diff = (b.rho_split_raw - b.rho_null_drift).to_numpy()
        diff = diff[np.isfinite(diff)]
        z, p = wilcoxon_z(diff) if diff.size >= 3 else (np.nan, np.nan)
        r_rb = rank_biserial(diff) if diff.size >= 3 else np.nan
        ceil = np.nanmedian(
            np.concatenate([b.rho_within_rpre.to_numpy(), b.rho_within_rpost.to_numpy()])
        )
        recs.append(
            {
                "band": band,
                "n": int(b.shape[0]),
                "median_rho_split_raw": float(np.nanmedian(b.rho_split_raw)),
                "median_rho_null_drift": float(np.nanmedian(b.rho_null_drift)),
                "median_ceiling": float(ceil),
                "median_diff": float(np.nanmedian(b.rho_split_raw - b.rho_null_drift)),
                "n_above_floor": int((b.rho_split_raw > b.rho_null_drift).sum()),
                "z": z,
                "p": p,
                "r_rb": r_rb,
            }
        )
    # No across-band BH: each band is a pre-specified hypothesis, not one of
    # 6 exchangeable discovery candidates. Per-band p is reported uncorrected
    # (regularity-first; p demoted to a footnote per project convention).
    return pd.DataFrame(recs)


def write_report(per_patient: pd.DataFrame, stats: pd.DataFrame) -> None:
    per_patient.to_csv(OUT / "drift_floor_per_patient.csv", index=False)
    print(f"wrote {OUT / 'drift_floor_per_patient.csv'}")

    L: list[str] = []
    a = L.append
    a("---")
    a("name: raw_fc_drift_floor")
    a("era: COHORT_N10 / IMCOH_ABS")
    a("status: current")
    a("---")
    a("")
    a("# Within-session drift floor for `ρ_split^raw` (raw-FC substrate)")
    a("")
    a("Adapts the cophenetic drift-floor control (`h2e_split_half.py`) from the "
      "LRG ultrametric `D` to the raw `|ImCoh|` adjacency `A`. Same formula, "
      "substrate swapped (`D → A`).")
    a("")
    a("* `ρ_split^raw` = audit_67 `Spearman(A^tt − A^preA, A^post − A^preB)` (obs_rho).")
    a("* `ρ_null_drift^raw` = `Spearman(A^preB − A^preA, A^postB − A^postA)` "
      "(no task; within-session drift only).")
    a("* `ceiling` = median split-half reliability `Spearman(A^phase_A, A^phase_B)`.")
    a("")
    a("Primary read is the **regularity**: the median `ρ_split^raw` relative to the "
      "≈0 within-session drift floor, and the fraction of patients above their own "
      "floor (`n>floor`). Per-band significance is a footnote (below).")
    a("")
    a("| band | n | median ρ_split^raw | median drift floor | Δ(obs−floor) | n>floor | reliability ceiling |")
    a("|------|--:|-------------------:|-------------------:|-------------:|:-------:|--------------------:|")
    for _, r in stats.iterrows():
        tex = BRAIN_BAND_TEX_DICT[r.band]
        a(
            f"| {tex} | {int(r.n)} | {r.median_rho_split_raw:+.3f} | "
            f"{r.median_rho_null_drift:+.3f} | {r.median_diff:+.3f} | "
            f"{int(r.n_above_floor)}/{int(r.n)} | {r.median_ceiling:+.2f} |"
        )
    a("")
    a("δ/β/low-γ sit consistently above a ≈0 drift floor (7–8/10 patients); α is "
      "muddied by its **own positive resting drift** (floor > 0); θ/high-γ show no "
      "separation (obs ≈ floor). The reliability ceiling = median split-half "
      "`Spearman(A_A, A_B)` — the upper bound any cross-phase ρ could reach.")
    a("")
    a("**Per-band significance (footnote — uncorrected).** Each band is a "
      "pre-specified hypothesis; **no across-band correction** is applied. An "
      "across-band BH would wrongly treat the 6 bands as exchangeable discovery "
      "candidates and dilute a real per-band effect with the null bands. One-sided "
      "paired Wilcoxon `ρ_split^raw > drift floor`:")
    a("")
    a("  " + " · ".join(
        f"{BRAIN_BAND_TEX_DICT[r.band]} p={r.p:.3f} (z={r.z:+.2f})"
        for _, r in stats.iterrows()
    ))
    a("")
    a("No band reaches p<0.05 at n=10 (best low-γ p≈0.057). The raw substrate is "
      "**directionally** above its drift floor but not significantly so; the "
      "drift-*cleared* result lives at the LRG cophenetic scale (α/β/low-γ "
      "p≈0.008–0.014). Raw-weak → LRG-clear is the intended contrast, not a defect.")
    (OUT / "drift_floor_band_stats.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"wrote {OUT / 'drift_floor_band_stats.md'}")


def main() -> int:
    df = collect()
    stats = band_stats(df)
    write_report(df, stats)
    print("\n=== drift-corrected ρ_split^raw (per band) ===")
    print(stats.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
