#!/usr/bin/env python3
"""Band × k unanimity landscape for H1/H2a/H2b/H2d/H3.

Produces:
  - data/outputs/figures/metric_exploration/band_k_landscape/{hyp}_imcoh_abs.pdf
    Panel (6 bands × 1 k-axis) showing at each k the fraction of contributing
    patients with positive contrast.  Contour lines at 50%, 71% (≥5/7),
    90% (≥9/10), and 100% (strict unanimity).  The dark ridges above 90%
    identify the scales where the hypothesis carries a measured cross-patient
    effect.
  - data/reports/imcoh_vi/band_k_landscape.md
    Per-(hypothesis, band) numerical summary: longest contiguous runs at
    three agreement levels (≥n-1, strict, all-but-two).

Reads:
  - ``data/reports/imcoh_vi/hypothesis_contrasts.csv`` (H1, H2a, H2b, H3 —
    produced by ``compute_imcoh_vi.py``).
  - ``data/reports/imcoh_vi/h2d_persistence_raw.csv`` (H2d — produced by
    ``h2d_coactivation_persistence.py``).

FC method is imcoh_abs by construction. Note that H2c (continuous drift
direction) has no per-k resolution and is therefore not in the landscape;
see ``data/reports/imcoh_vi/figures/h2c_unanimity.pdf`` for the H2c
patient × band heatmap.
"""
from __future__ import annotations

from itertools import groupby
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT
from lrg_eegfc.config.paths import FIGURES_ROOT, REPORTS_ROOT

HYPS = ["H1", "H2a", "H2b", "H2d", "H3"]
HYP_LABEL = {
    "H1":  "H1 — task stability   (VI(TL,TT) < others)",
    "H2a": "H2a — task trace     (VI(rest_pre,rest_post) > VI(task_test,rest_post))",
    "H2b": "H2b — approach/exit  (VI(rest_pre,task_test) > VI(task_test,rest_post))",
    "H2d": "H2d — block persistence   (Δρ = ρ_task − ρ_inert > 0)",
    "H3":  "H3 — within < cross",
}


def runs(ks: list[int]) -> list[tuple[int, int]]:
    ks = sorted(set(int(k) for k in ks))
    out = []
    for _, grp in groupby(enumerate(ks), key=lambda ix: ix[1] - ix[0]):
        vals = [v for _, v in grp]
        out.append((vals[0], vals[-1]))
    return out


def main() -> None:
    src = REPORTS_ROOT / "imcoh_vi" / "hypothesis_contrasts.csv"
    if not src.exists():
        raise FileNotFoundError(
            f"Missing {src}.  Run compute_imcoh_vi.py first."
        )
    df = pd.read_csv(src)

    # H2d lives in a separate raw CSV (k-resolved Δρ rather than VI).
    h2d_src = REPORTS_ROOT / "imcoh_vi" / "h2d_persistence_raw.csv"
    if h2d_src.exists():
        h2d_raw = pd.read_csv(h2d_src)
        # Convert to the same (hypothesis, band, k, patient, contrast, sign) schema
        h2d_df = pd.DataFrame({
            "hypothesis": "H2d",
            "band":       h2d_raw["band"],
            "k":          h2d_raw["k"].astype(int),
            "patient":    h2d_raw["patient"],
            "contrast":   h2d_raw["delta_rho"],
            "sign":       np.sign(h2d_raw["delta_rho"]).astype(int),
        })
        df = pd.concat([df, h2d_df], ignore_index=True)
    else:
        print(f"  (skipping H2d landscape — {h2d_src} not found)")

    # Per (hyp, band, k) fraction of contributing patients with positive sign.
    frac = (
        df.assign(pos=(df["sign"] > 0).astype(int))
          .groupby(["hypothesis", "band", "k"])
          .agg(n=("sign", "size"), pos=("pos", "sum"))
          .reset_index()
    )
    frac["frac_pos"] = frac["pos"] / frac["n"]

    fig_dir = FIGURES_ROOT / "metric_exploration" / "band_k_landscape"
    fig_dir.mkdir(parents=True, exist_ok=True)

    cmap = LinearSegmentedColormap.from_list(
        "ridge",
        [(1.0, 1.0, 1.0), (1.0, 0.95, 0.7),
         (1.0, 0.6, 0.2), (0.7, 0.05, 0.05)],
        N=256,
    )

    # ── Figure per hypothesis ─────────────────────────────────────────
    for hyp in HYPS:
        hdf = frac[frac["hypothesis"] == hyp]
        if hdf.empty:
            continue
        k_min, k_max = int(hdf["k"].min()), int(hdf["k"].max())
        k_grid = np.arange(k_min, k_max + 1)
        # Matrix (bands × k) of fraction-positive; NaN where no data.
        M = np.full((len(BRAIN_BANDS_NAMES), len(k_grid)), np.nan)
        N = np.full_like(M, np.nan)
        for ib, band in enumerate(BRAIN_BANDS_NAMES):
            bdf = hdf[hdf["band"] == band].set_index("k")
            for ik, k in enumerate(k_grid):
                if k in bdf.index:
                    M[ib, ik] = bdf.at[k, "frac_pos"]
                    N[ib, ik] = bdf.at[k, "n"]

        fig, ax = plt.subplots(figsize=(11.5, 3.6))
        im = ax.imshow(
            M, aspect="auto", cmap=cmap, vmin=0.0, vmax=1.0,
            extent=[k_grid[0] - 0.5, k_grid[-1] + 0.5,
                    len(BRAIN_BANDS_NAMES) - 0.5, -0.5],
            interpolation="nearest",
        )

        # Overlay black dots for strict unanimity (frac_pos == 1.0)
        unan_yx = np.argwhere(M >= 0.9999)
        if unan_yx.size:
            ys, xs = unan_yx[:, 0], unan_yx[:, 1]
            ax.scatter(k_grid[xs], ys, s=3, c="black", marker="s")

        ax.set_yticks(range(len(BRAIN_BANDS_NAMES)))
        ax.set_yticklabels([BRAIN_BAND_TEX_DICT[b] for b in BRAIN_BANDS_NAMES])
        ax.set_xlabel("k (number of communities)")
        ax.set_title(HYP_LABEL[hyp], loc="left", fontsize=11)
        cb = fig.colorbar(im, ax=ax, shrink=0.9, pad=0.01)
        cb.set_label("fraction of contributing patients with positive contrast")
        # Levels annotated on colorbar
        for lev, lab in [(0.5, "majority"), (0.71, "≥5/7"),
                         (0.90, "≥9/10"), (1.00, "unanimous")]:
            cb.ax.axhline(lev, color="black", lw=0.6)
            cb.ax.text(1.2, lev, lab, transform=cb.ax.transData,
                       fontsize=8, va="center")
        fig.tight_layout()

        out = fig_dir / f"{hyp}_imcoh_abs.pdf"
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
        print(f"  saved {out}")

    # ── Numerical summary markdown ────────────────────────────────────
    md_path = REPORTS_ROOT / "imcoh_vi" / "band_k_landscape.md"
    lines: list[str] = []
    ap = lines.append
    ap("# Band × k unanimity landscape on H1 / H2a / H2b / H2d / H3")
    ap("")
    ap("FC method: `imcoh_abs`. Contributing patients per hypothesis depend on "
       "which phases each patient has (Pat_10 skipped on rest↔task pairs, "
       "Pat_13 no rest_pre, Pat_14 no task_test).")
    ap("")
    ap("`unan+` = all contributing patients positive (strict).  "
       "`≥n-1` = at most one dissenter.  "
       "`≥n-2` = at most two dissenters.  "
       "Runs shown are the contiguous k ranges where the threshold holds "
       "(length L in k steps; total k span is 2..N-1 per patient).")
    ap("")
    for hyp in HYPS:
        hdf = frac[frac["hypothesis"] == hyp]
        if hdf.empty:
            continue
        ap(f"## {HYP_LABEL[hyp]}")
        ap("")
        ap("| band | k cells | n (med) | unan+ | ≥n-1 + | ≥n-2 + | "
           "longest unan+ | longest ≥n-1 + |  ≥n-1 contiguous runs (L≥5) |")
        ap("|------|--------:|--------:|------:|------:|------:|--------------:|---------------:|-----------------------------|")
        for band in BRAIN_BANDS_NAMES:
            bdf = hdf[hdf["band"] == band]
            if bdf.empty:
                continue
            unan_k, near1_k, near2_k = [], [], []
            n_list = []
            for _, row in bdf.iterrows():
                k, n, pos = int(row["k"]), int(row["n"]), int(row["pos"])
                n_list.append(n)
                if pos == n:       unan_k.append(k)
                if pos >= n - 1:   near1_k.append(k)
                if pos >= n - 2:   near2_k.append(k)
            n_med = int(np.median(n_list)) if n_list else 0
            r_unan  = runs(unan_k)
            r_near1 = runs(near1_k)
            long_unan  = max((b - a + 1 for a, b in r_unan),  default=0)
            long_near1 = max((b - a + 1 for a, b in r_near1), default=0)
            big = [f"{a}-{b}(L={b-a+1})" for a, b in r_near1 if (b - a + 1) >= 5]
            ap(f"| {BRAIN_BAND_TEX_DICT[band]} | {len(bdf):>8} | {n_med:>7} | "
               f"{len(unan_k):>5} | {len(near1_k):>5} | {len(near2_k):>5} | "
               f"{long_unan:>13} | {long_near1:>14} | "
               f"{', '.join(big) if big else '—'} |")
        ap("")

    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  saved {md_path}")


if __name__ == "__main__":
    main()
