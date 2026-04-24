#!/usr/bin/env python3
"""H2d on |ImCoh|² — direct comparison with |ImCoh| at n=9.

The user's testable claim: "with |ImCoh|² band-selectivity was clearer
than with |ImCoh|." This script runs the same H2d pipeline on both FC
metrics and produces a side-by-side figure.

Steps:
  1. For each (patient, phase, band), derive |ImCoh|² FC from the
     freq-resolved ImCoh cache (via `load_fc_matrix`) and compute LRG,
     caching under `fc_method="imcoh_sq"`.
  2. Run the H2d per-k block-persistence computation on |ImCoh|² using
     the SAME `compute_persistence` helper as the canonical H2d.
  3. Load the existing |ImCoh| H2d raw CSV for comparison.
  4. Build the band-unanimity curves for both metrics, side-by-side.

Outputs:
    data/reports/imcoh_vi/h2d_imcoh_sq_raw.csv
    data/reports/imcoh_vi/figures/imcoh_sq_vs_abs.{pdf,png}

Memory: FC loads are one band at a time, LRG caches persist, peak RSS
typically < 2 GB. Runs ~30-60 min on first call; instant on cache hit.
"""
from __future__ import annotations

import argparse
import gc
import resource

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PATIENTS_LIST
from lrg_eegfc.config.paths import IMCOH_LRG_CACHE, REPORTS_ROOT
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import compute_lrg_analysis, load_lrg_result

from h2d_coactivation_persistence import compute_persistence, K_RANGE


REST_PHASES = ("rest_pre", "rest_post")
TASK_PHASES = ("task_learn", "task_test")
ALL_PHASES = REST_PHASES + TASK_PHASES


def _set_mem_cap(gb: float) -> None:
    try:
        resource.setrlimit(resource.RLIMIT_AS,
                            (int(gb * 1024 ** 3), int(gb * 1024 ** 3)))
    except (ValueError, OSError):
        pass


# ─────────────────────────── step 1: populate caches ───────────────────────────


def ensure_imcoh_sq_lrg(verbose: bool = False) -> None:
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            for phase in ALL_PHASES:
                try:
                    A = load_fc_matrix(pat, phase, band, fc_method="imcoh_sq")
                except (FileNotFoundError, OSError):
                    continue
                if A is None:
                    continue
                compute_lrg_analysis(
                    np.asarray(A), pat, phase, band,
                    fc_method="imcoh_sq",
                    cache_root=IMCOH_LRG_CACHE,
                    use_cache=True,
                )
                del A
            gc.collect()
        if verbose:
            print(f"  {pat}: imcoh_sq LRG caches populated")


# ─────────────────────────── step 2: H2d on imcoh_sq ───────────────────────────


def load_Z(pat: str, phase: str, band: str, fc_method: str) -> np.ndarray | None:
    try:
        r = load_lrg_result(pat, phase, band, fc_method=fc_method,
                            cache_root=IMCOH_LRG_CACHE)
    except Exception:
        return None
    return np.asarray(r.linkage_matrix) if r is not None else None


def collect_h2d(fc_method: str) -> pd.DataFrame:
    rows = []
    for pat in PATIENTS_LIST:
        for band in BRAIN_BANDS_NAMES:
            Z_rpre  = load_Z(pat, "rest_pre",  band, fc_method)
            Z_ttest = load_Z(pat, "task_test", band, fc_method)
            Z_rpost = load_Z(pat, "rest_post", band, fc_method)
            for k in K_RANGE:
                out = compute_persistence(Z_rpre, Z_ttest, Z_rpost, k)
                if out is None:
                    continue
                rows.append({"patient": pat, "band": band, "k": k,
                             "rho_task": out[0], "rho_inert": out[1],
                             "delta_rho": out[0] - out[1]})
    return pd.DataFrame(rows)


# ─────────────────────────── step 3: figure ───────────────────────────


BAND_COLORS = {
    "delta":      "#3b6e9c",
    "theta":      "#c44e4e",
    "alpha":      "#5ea85e",
    "beta":       "#b28ad1",
    "low_gamma":  "#e5a24b",
    "high_gamma": "#7d5a50",
}


def _frac_strong(rho_task: np.ndarray, rho_inert: np.ndarray,
                 thresh: float = 2.0) -> np.ndarray:
    finite = np.isfinite(rho_task) & np.isfinite(rho_inert)
    strong = finite & (rho_task > thresh * rho_inert)
    n = finite.sum(axis=0)
    s = strong.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, s / n, np.nan)


def _frac_positive(mat: np.ndarray) -> np.ndarray:
    finite = np.isfinite(mat)
    pos = finite & (mat > 0)
    n = finite.sum(axis=0)
    p = pos.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, p / n, np.nan)


def _build_mats(df: pd.DataFrame, patients: list[str], ks: list[int]):
    per_band: dict[str, dict[str, np.ndarray]] = {}
    for band in BRAIN_BANDS_NAMES:
        bdf = df[df["band"] == band]
        d  = np.full((len(patients), len(ks)), np.nan)
        rt = np.full_like(d, np.nan)
        ri = np.full_like(d, np.nan)
        for ip, pat in enumerate(patients):
            sub = bdf[bdf["patient"] == pat].set_index("k")
            for jk, k in enumerate(ks):
                if k in sub.index:
                    d[ip, jk]  = sub.at[k, "delta_rho"]
                    rt[ip, jk] = sub.at[k, "rho_task"]
                    ri[ip, jk] = sub.at[k, "rho_inert"]
        per_band[band] = {"delta": d, "rho_task": rt, "rho_inert": ri}
    return per_band


def draw_panel(ax, per_band, ks, title: str) -> dict[str, float]:
    strength: dict[str, float] = {}
    for band in BRAIN_BANDS_NAMES:
        curve = _frac_strong(per_band[band]["rho_task"],
                             per_band[band]["rho_inert"])
        ax.plot(ks, curve, color=BAND_COLORS[band], lw=2.2,
                label=BRAIN_BAND_TEX_DICT[band])
        strength[band] = float(np.nanmean(curve))
    ax.axhline(0.5, color="#888", lw=0.8, linestyle="--")
    ax.axhline(7/9, color="#000", lw=0.5, linestyle=":")
    ax.set_ylim(0.0, 1.05)
    ax.set_title(title, fontsize=11, loc="left")
    ax.set_ylabel(r"fraction with $\rho_{task} > 2\,\rho_{inert}$", fontsize=10)
    return strength


def figure_comparison(df_abs: pd.DataFrame, df_sq: pd.DataFrame) -> None:
    patients = sorted(set(df_abs["patient"].unique()) |
                      set(df_sq["patient"].unique()))
    ks = sorted(set(df_abs["k"].unique()) | set(df_sq["k"].unique()))

    per_abs = _build_mats(df_abs, patients, ks)
    per_sq  = _build_mats(df_sq,  patients, ks)

    fig, axes = plt.subplots(2, 1, figsize=(12.0, 8.5), dpi=160, sharex=True)
    s_abs = draw_panel(axes[0], per_abs, ks,
                        "|ImCoh| — band-selective H2d block-persistence curves")
    s_sq  = draw_panel(axes[1], per_sq, ks,
                        "|ImCoh|² — same analysis on the squared metric")
    axes[1].set_xlabel("k (dendrogram cut scale)", fontsize=11)
    axes[0].legend(ncol=6, loc="lower right", frameon=False, fontsize=10)

    # Textbox with band ranking comparison
    rows = ["band   |ImCoh|  |ImCoh|²"]
    for band in BRAIN_BANDS_NAMES:
        rows.append(f"{BRAIN_BAND_TEX_DICT[band]:>8}"
                    f"  {s_abs[band]:.3f}  {s_sq[band]:.3f}")
    textbox = "\n".join(rows)
    fig.text(0.82, 0.55, textbox, fontsize=9, family="monospace",
             bbox=dict(facecolor="white", edgecolor="#888", alpha=0.95))

    fig.suptitle(
        "H2d band-selectivity: |ImCoh| vs |ImCoh|² on n=9.\n"
        "Higher curves = stronger band-specific multiscale persistence; "
        "θ (red) is expected to be lowest if band-selectivity holds.",
        fontsize=12, y=1.00,
    )
    fig.tight_layout()
    out_dir = REPORTS_ROOT / "imcoh_vi" / "figures"
    for ext in ("pdf", "png"):
        out = out_dir / f"imcoh_sq_vs_abs.{ext}"
        fig.savefig(out, bbox_inches="tight")
        print(f"saved {out}")
    plt.close(fig)


# ─────────────────────────── main ───────────────────────────


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--mem-gb", type=float, default=10.0)
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("--skip-compute", action="store_true",
                   help="Assume imcoh_sq LRG caches already populated.")
    args = p.parse_args()

    _set_mem_cap(args.mem_gb)

    if not args.skip_compute:
        print("Populating imcoh_sq LRG caches...")
        ensure_imcoh_sq_lrg(verbose=args.verbose)

    print("Running H2d on imcoh_sq...")
    df_sq = collect_h2d("imcoh_sq")
    df_sq.to_csv(REPORTS_ROOT / "imcoh_vi" / "h2d_imcoh_sq_raw.csv", index=False)
    print(f"  {len(df_sq)} rows for imcoh_sq")

    print("Loading H2d on imcoh_abs from existing CSV...")
    df_abs = pd.read_csv(REPORTS_ROOT / "imcoh_vi" / "h2d_persistence_raw.csv")
    print(f"  {len(df_abs)} rows for imcoh_abs")

    figure_comparison(df_abs, df_sq)


if __name__ == "__main__":
    main()
