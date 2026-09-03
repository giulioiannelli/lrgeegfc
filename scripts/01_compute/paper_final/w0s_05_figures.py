#!/usr/bin/env python3
"""Figures for lane W0-S -- what the scale axis looks like from the readout's side.

Three panels-of-record, in the order the argument is made.

fig_scale_axis_degeneracy
    Does the readout's INPUT change along the scale axis at all? Left: cross-scale
    Spearman of the baseline cophenetic vector. Middle: cross-scale overlap of the
    per-pair contributions -- the carriers of the trace itself -- with the
    matched-strength reference. Right: the locked scale units (N_eff, the
    communication neighbourhood, and the number of distinct cophenetic values the
    rank statistic actually sees). READING RULE: a similarity matrix that stays
    near +1 everywhere is a statistic scored on the same information at every
    scale, whatever the hierarchy is doing.

fig_scale_stratum_margin_surface
    The (scale x tree level) surface of the cohort margin, for the additive
    contribution decomposition and for the within-stratum local trace, with the
    stratum the diffusion is currently resolving overdrawn. READING RULE: colour
    in horizontal bands = no scale structure; colour tracking the overdrawn line =
    the trace is carried by whichever pairs the diffusion is resolving.

fig_scale_variability_vs_null
    Does the effect vary along the axis more than its own noise? Observed
    cross-patient shape agreement and effective number of independent scales,
    each against the held-out-surrogate null band for that same readout.
    READING RULE: a bar that does not stand clear of its null band is noise, not
    scale structure.
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals import imshow_colorbar_caxdivider
from lrg_eegfc.visuals.styles import band_color, use_lrg_style

ROOT = setup_script_env()

from w0s_00_artifacts import load_cells                       # noqa: E402

BASE = Path(os.environ.get(
    "W0S_OUT_ANALYSIS", ROOT / "data" / "paper_final" / "lane_s_scale"))
FIGS = BASE / "figures"
REF_BAND = "beta"
NQ = 5


def _tex(b: str) -> str:
    return BRAIN_BAND_TEX_DICT.get(b, b)


def _logticks(ax, s, axis="both"):
    idx = [0, len(s) // 3, 2 * len(s) // 3, len(s) - 1]
    lab = [f"{s[i]:.2g}" for i in idx]
    if axis in ("both", "x"):
        ax.set_xticks(idx)
        ax.set_xticklabels(lab)
    if axis in ("both", "y"):
        ax.set_yticks(idx)
        ax.set_yticklabels(lab)


def fig_degeneracy(cs) -> None:
    s = cs.s
    xs = np.nanmedian(cs.diag(REF_BAND, "xs_coph"), axis=0)
    car = BASE / "carrier" / f"{cs.patients[0]}__{REF_BAND}.npz"
    cos = None
    if (BASE / "carrier").exists():
        stack, nulls = [], []
        for p in cs.patients:
            f = BASE / "carrier" / f"{p}__{REF_BAND}.npz"
            if f.exists():
                z = np.load(f)
                stack.append(z["cos"])
                nulls.append(np.nanmedian(z["cos_null"], axis=0))
        if stack:
            cos = np.nanmedian(np.stack(stack), axis=0)
            cos_null = np.nanmedian(np.stack(nulls), axis=0)

    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.3), constrained_layout=True)
    im = axes[0].imshow(xs, origin="lower", vmin=0, vmax=1, cmap="viridis")
    axes[0].set_xlabel(r"scale $s$")
    axes[0].set_ylabel(r"scale $s$")
    _logticks(axes[0], s)
    imshow_colorbar_caxdivider(im, axes[0], size="4%", pad=0.08)
    axes[0].set_title(r"input: $\rho_S\,[\,C_A(s_i),\,C_A(s_j)\,]$", fontsize=9)

    if cos is not None:
        im = axes[1].imshow(cos, origin="lower", vmin=0, vmax=1, cmap="viridis")
        _logticks(axes[1], s)
        axes[1].set_xlabel(r"scale $s$")
        imshow_colorbar_caxdivider(im, axes[1], size="4%", pad=0.08)
        axes[1].set_title("carriers: overlap of per-pair contributions",
                          fontsize=9)
    else:
        axes[1].axis("off")

    ax = axes[2]
    for key, lab, ls in (("n_eff", r"$N_{\mathrm{eff}}(s)$", "-"),
                         ("n_distinct", "distinct cophenetic values", "--"),
                         ("m_comm", r"$m(s)$ neighbourhood", ":")):
        y = np.nanmedian(cs.diag(REF_BAND, key), axis=0)
        ax.plot(s, y, ls, color=band_color(REF_BAND), label=lab)
    if cos is not None:
        ax.plot(s, np.diagonal(np.ones_like(cos)) * np.nan, alpha=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"scale $s=\tau\lambda_{\max}$")
    ax.set_ylabel("count")
    ax.legend(frameon=False, fontsize=7, loc="lower left")
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "fig_scale_axis_degeneracy.pdf", transparent=True)
    plt.close(fig)


def fig_surface(cs) -> None:
    s = cs.s
    bands = cs.bands
    fig, axes = plt.subplots(2, len(bands), figsize=(2.05 * len(bands), 4.4),
                             constrained_layout=True, sharex=True, sharey=True)
    for row, (pre, lab) in enumerate((("Qcon_q", "contribution"),
                                      ("Qloc_q", "local trace"))):
        surf = []
        for b in bands:
            M = np.stack([cs.margins(b, f"{pre}{q+1}")[0] for q in range(NQ)])
            surf.append(np.nanmedian(M, axis=1))          # (NQ, nS)
        vmax = np.nanmax(np.abs(np.stack(surf)))
        vmax = float(vmax) if np.isfinite(vmax) and vmax > 0 else 1.0
        for k, b in enumerate(bands):
            ax = axes[row, k]
            im = ax.pcolormesh(np.arange(len(s) + 1) - 0.5,
                               np.arange(NQ + 2) - 0.5,
                               np.vstack([surf[k], np.full((1, len(s)), np.nan)]),
                               cmap="RdBu_r",
                               norm=TwoSlopeNorm(0.0, -vmax, vmax),
                               shading="flat")
            qs = np.nanmedian(cs.diag(b, "qstar"), axis=0) - 1
            ax.plot(np.arange(len(s)), qs, color="k", lw=1.0)
            ax.set_ylim(-0.5, NQ - 0.5)
            if row == 0:
                ax.set_title(_tex(b), fontsize=9)
            if row == 1:
                ax.set_xlabel(r"scale $s$")
                _logticks(ax, s, axis="x")
            if k == 0:
                ax.set_ylabel(f"{lab}\ntree level (quintile)")
                ax.set_yticks(range(NQ))
                ax.set_yticklabels([str(q + 1) for q in range(NQ)])
        imshow_colorbar_caxdivider(im, axes[row, -1], size="4%", pad=0.08)
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "fig_scale_stratum_margin_surface.pdf", transparent=True)
    plt.close(fig)


def fig_variability(ax_df) -> None:
    head = ["T", "Ccon_diag", "Tloc_diag", "Qcon_diag", "Qloc_diag", "Thei"]
    d = ax_df[ax_df.readout.isin(head)]
    bands = list(dict.fromkeys(d.band))
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.2), constrained_layout=True)
    w = 0.8 / len(head)
    for ai, (col, nullcol, lab) in enumerate((
            ("shape_agreement", "shape_agreement_null_med",
             "cross-patient agreement on profile shape"),
            ("n_eff_pr", "n_eff_null_med",
             r"$n_{\mathrm{eff}}$ (independent scales)"))):
        ax = axes[ai]
        for i, m in enumerate(head):
            sub = d[d.readout == m].set_index("band").reindex(bands)
            x = np.arange(len(bands)) + i * w - 0.4 + w / 2
            ax.bar(x, sub[col].to_numpy(), width=w * 0.9,
                   color=[band_color(b) for b in bands],
                   alpha=0.35 + 0.65 * (m == "T"),
                   edgecolor="k", linewidth=0.4, label=m if ai == 0 else None)
            ax.plot(x, sub[nullcol].to_numpy(), "_", color="k", ms=6, mew=1.2)
        ax.set_xticks(np.arange(len(bands)))
        ax.set_xticklabels([_tex(b) for b in bands])
        ax.set_ylabel(lab)
        ax.axhline(0.0, color="0.6", lw=0.6)
    axes[0].legend(frameon=False, fontsize=7, ncol=3)
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / "fig_scale_variability_vs_null.pdf", transparent=True)
    plt.close(fig)


def main() -> None:
    use_lrg_style()
    cs = load_cells()
    fig_degeneracy(cs)
    fig_surface(cs)
    p = BASE / "axis_and_neff.csv"
    if p.exists():
        fig_variability(pd.read_csv(p))
    print(f"[w0s-figures] -> {FIGS}", flush=True)


if __name__ == "__main__":
    main()
