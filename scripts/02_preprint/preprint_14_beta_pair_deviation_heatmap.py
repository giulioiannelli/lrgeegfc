#!/usr/bin/env python3
r"""β cohort joint-rank deviation heatmap (standalone single-panel figure).

Pools the within-patient rank pairs ``(u, v) = (rank Δ_task, rank Δ_rest)``
of cophenet-distance differences across all ten cohort patients and colours
each cell on a ``n_bins × n_bins`` grid of ``[0, 1]^2`` by the relative
deviation from the uniform null

.. math::

    z(u_{\rm bin}, v_{\rm bin}) =
        \frac{P_{\rm obs}(u_{\rm bin}, v_{\rm bin})}{P_{\rm null}} - 1,
    \qquad P_{\rm null} = 1\ \text{on}\ [0,1]^2.

Reframes the noisy 60-bin density (which dissolves the +0.22 cohort
ρ_split^coph into per-cell sampling fluctuations) into a coarse-bin
view where each cell pools enough pairs to surface the Gaussian-
copula-predicted diagonal enrichment / anti-diagonal depletion.

Usage
-----
    python preprint_14_beta_pair_deviation_heatmap.py
    python preprint_14_beta_pair_deviation_heatmap.py --band alpha --bins 30
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import rankdata, spearmanr

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()


COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]

PAIR_SPLIT_DIR = ROOT / "data" / "reports" / "imcoh_continuous_trace" / "per_pair_split"

# Palette anchored to the project's trace/anti convention.
CLR_TRACE = "#5fa844"   # grass green
CLR_ANTI = "#c0392b"    # brick red


def load_pair_data(pat: str, band: str) -> tuple[np.ndarray, np.ndarray]:
    npz_path = PAIR_SPLIT_DIR / f"{pat}_{band}.npz"
    d = np.load(npz_path)
    return (np.asarray(d["dD_task"], dtype=float),
            np.asarray(d["dD_rest"], dtype=float))


def pool_within_patient_ranks(band: str) -> tuple[np.ndarray, np.ndarray]:
    us, vs = [], []
    for pat in COHORT:
        dt, dr = load_pair_data(pat, band)
        n = len(dt)
        if n < 2:
            continue
        u = (rankdata(dt, method="ordinal") - 0.5) / n
        v = (rankdata(dr, method="ordinal") - 0.5) / n
        us.append(u); vs.append(v)
    return np.concatenate(us), np.concatenate(vs)


def render(band: str, n_bins: int, vmax: float) -> Path:
    band_tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")
    u_pool, v_pool = pool_within_patient_ranks(band)
    n_total = u_pool.size

    edges = np.linspace(0.0, 1.0, n_bins + 1)
    H, _, _ = np.histogram2d(u_pool, v_pool, bins=[edges, edges])
    expected = n_total / (n_bins * n_bins)
    z = H / expected - 1.0

    cmap_dev = LinearSegmentedColormap.from_list(
        "trace_dev",
        [CLR_ANTI, "#f7f3eb", CLR_TRACE],
        N=256,
    )

    fig = plt.figure(figsize=(7.6, 7.0))
    gs = fig.add_gridspec(
        1, 1, left=0.11, right=0.87, top=0.93, bottom=0.10,
    )
    ax = fig.add_subplot(gs[0, 0])

    im = ax.imshow(z.T, origin="lower", extent=[0.0, 1.0, 0.0, 1.0],
                   cmap=cmap_dev, vmin=-vmax, vmax=+vmax,
                   aspect="equal", interpolation="nearest", zorder=1)

    ax.plot([0, 1], [0, 1], color="0.18", lw=1.4, ls="--",
            alpha=0.55, zorder=3)
    ax.plot([0, 1], [1, 0], color="0.18", lw=0.9, ls=":",
            alpha=0.35, zorder=3)

    ax.set_xlim(0.0, 1.0); ax.set_ylim(0.0, 1.0)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_xticklabels(["0", r"$\frac{1}{2}$", "1"])
    ax.set_yticklabels(["0", r"$\frac{1}{2}$", "1"])
    ax.set_xlabel(r"within-patient rank of "
                  r"$\Delta_{\mathrm{task}}(i,j)$  ($u$)",
                  fontsize=10.5, labelpad=4)
    ax.set_ylabel(r"within-patient rank of "
                  r"$\Delta_{\mathrm{rest}}(i,j)$  ($v$)",
                  fontsize=10.5, labelpad=4)
    ax.tick_params(labelsize=9.5)

    cbar_ax = fig.add_axes([0.89, 0.10, 0.025, 0.83])
    cb = fig.colorbar(im, cax=cbar_ax, orientation="vertical")
    cb.set_label(r"$(P_{\rm obs} / P_{\rm null}) - 1$",
                 rotation=270, labelpad=18, fontsize=10)
    cb.ax.tick_params(labelsize=9)
    cb.set_ticks([-vmax, -vmax/2, 0.0, +vmax/2, +vmax])
    cb.set_ticklabels([rf"$-{vmax:.2f}$",
                       rf"$-{vmax/2:.2f}$", "$0$",
                       rf"$+{vmax/2:.2f}$",
                       rf"$+{vmax:.2f}$"])

    fig.text(0.04, 0.965,
             rf"{band_tex} cohort joint-rank density vs uniform null  "
             rf"($n_{{\rm bins}}={n_bins}$, $n_{{\rm pairs}}={n_total:,}$, "
             r"$n_{\rm pat}=10$)",
             ha="left", va="top", fontsize=11.0, fontweight="bold")

    pooled_rho, _ = spearmanr(u_pool, v_pool)
    print(f"  pooled ρ_split^coph = {pooled_rho:+.4f}")
    print(f"  n_pairs = {n_total:,}  (10 patients)")
    print(f"  diagonal-corner cell z = {float(z[-1, -1]):+.3f}")
    print(f"  anti-diagonal-corner cell z = {float(z[-1, 0]):+.3f}")

    out_dir = ROOT / "data" / "preprint" / "figures" / band / "per_pair_trace"
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"fig_{band}_deviation_heatmap.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


def main() -> Path:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--band", default="beta",
                        choices=["delta", "theta", "alpha", "beta",
                                 "low_gamma", "high_gamma"])
    parser.add_argument("--bins", type=int, default=20,
                        help="Bin count per side (default 20).")
    parser.add_argument("--vmax", type=float, default=0.60,
                        help="Symmetric cmap clip on relative deviation "
                             "(default 0.60).")
    args = parser.parse_args()
    return render(args.band, args.bins, args.vmax)


if __name__ == "__main__":
    main()
