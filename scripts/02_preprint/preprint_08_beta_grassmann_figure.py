#!/usr/bin/env python3
"""Global subspace persistence (Grassmann `d_G(k)`) preprint figure (per band).

Two-panel composite for any band's Grassmann paragraph of the preprint:

Panel (a)  cohort-median trace scalar T_G(k) over the full k range
           with the matched-strength surrogate cohort-median curve
           overlaid. Wilcoxon p<0.05 cells marked uniformly (no
           pre-framed window). Peak annotation at the cohort min-p k.

Panel (b)  same plot under epi-zone exclusion (audit_67).

Usage
-----
    python preprint_08_beta_grassmann_figure.py             # default: beta
    python preprint_08_beta_grassmann_figure.py --band alpha
    python preprint_08_beta_grassmann_figure.py --band low_gamma

Inputs
------
data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv
data/audit/grassmann_matched_strength_surrogate/per_patient_per_band_per_k.csv
data/audit/grassmann_epi_exclusion/cohort_summary.csv
data/audit/grassmann_epi_exclusion/per_patient_per_band_per_k.csv

Output (PDF only, no PNG sibling)
---------------------------------
data/preprint/figures/<band>/grassmann/fig_<band>_grassmann.pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()


GRASS_DIR = ROOT / "data" / "audit" / "grassmann_matched_strength_surrogate"
GRASS_EPIX_DIR = ROOT / "data" / "audit" / "grassmann_epi_exclusion"

CLR_OBS = "#1f3d6e"
CLR_OBS_EPIX = "#7f3b2c"
CLR_SURR = "#7d7d7d"
CLR_SURR_FILL = "#dcdcdc"
CLR_SIG = "#1f7a1f"
CLR_WIN = "#d9efd2"        # pale window fill
CLR_PEAK = "#2c2c2c"


def load_cohort(band: str, cohort_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(cohort_csv)
    df = df[df.band == band].sort_values("k").reset_index(drop=True)
    return df


def load_per_patient(band: str, per_pat_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(per_pat_csv)
    df = df[df.band == band].sort_values(["k", "patient"]).reset_index(drop=True)
    return df


def cohort_surrogate_envelope(per_pat: pd.DataFrame) -> pd.DataFrame:
    """Per-k cohort summary of patient surrogate distributions.

    Returns columns: k, surr_p50_cohort, surr_p25_cohort, surr_p75_cohort
    where the cohort statistics are taken across patients of each
    patient's own surrogate quantile.
    """
    g = per_pat.groupby("k")
    out = pd.DataFrame({
        "k": np.array(sorted(per_pat.k.unique())),
    })
    out["surr_p50_cohort"] = [
        np.nanmedian(g.get_group(k).surr_T_G_p50.values) for k in out.k
    ]
    out["surr_p25_cohort"] = [
        np.nanmedian(g.get_group(k).surr_T_G_p25.values) for k in out.k
    ]
    out["surr_p75_cohort"] = [
        np.nanmedian(g.get_group(k).surr_T_G_p75.values) for k in out.k
    ]
    return out


def longest_contiguous(mask: np.ndarray) -> tuple[int, int, int]:
    best = 0
    bs = -1
    be = -1
    cur = 0
    cs = -1
    for i, m in enumerate(mask):
        if m:
            if cur == 0:
                cs = i
            cur += 1
            if cur > best:
                best = cur
                bs = cs
                be = i
        else:
            cur = 0
    return best, bs, be


def draw_panel(ax: plt.Axes, cohort: pd.DataFrame, env: pd.DataFrame,
               label_obs: str, clr_obs: str, title: str,
               annotate_peak: bool, show_xticks: bool) -> None:
    ks = cohort.k.values.astype(int)
    obs_med = cohort.obs_median_T_G.values.astype(float)
    surr_med = cohort.surr_median_T_G_per_patient_median.values.astype(float)
    p_vals = cohort.paired_wilcoxon_p.values.astype(float)
    n_below = cohort.n_patients_below_own_surrogate.values.astype(int)

    env_k = env.k.values.astype(int)
    env_p25 = env.surr_p25_cohort.values.astype(float)
    env_p75 = env.surr_p75_cohort.values.astype(float)

    sig_mask = p_vals < 0.05
    n_sig = int(sig_mask.sum())
    run_len, run_s, run_e = longest_contiguous(sig_mask)
    if run_s >= 0:
        k_run_lo = int(ks[run_s])
        k_run_hi = int(ks[run_e])
    else:
        k_run_lo = k_run_hi = -1

    # Mark every Wilcoxon-significant cell uniformly (no pre-framed
    # contiguous window). The cohort gate is Wilcoxon p<0.05 only.
    for k_sig in ks[sig_mask]:
        ax.axvspan(k_sig - 0.5, k_sig + 0.5, facecolor=CLR_WIN,
                    alpha=0.15, zorder=0)

    # Surrogate IQR envelope (per-patient p25/p75 medians across cohort)
    ax.fill_between(env_k, env_p25, env_p75, color=CLR_SURR_FILL,
                    alpha=0.95, zorder=1,
                    label="matched-strength surrogate IQR")
    ax.plot(env_k, np.array([
        np.nanmedian([env_p25[i], env_p75[i]]) for i in range(env_k.size)
    ]), color=CLR_SURR, lw=0.0, zorder=1)
    ax.plot(ks, surr_med, color=CLR_SURR, lw=1.0, zorder=2,
            label="surrogate cohort median")

    # Observed cohort median
    ax.plot(ks, obs_med, color=clr_obs, lw=1.6, zorder=4,
            label=label_obs)

    # Zero line
    ax.axhline(0, color="0.4", lw=0.6, ls="--", zorder=1)

    # Wilcoxon p<0.05 tick markers along the top
    if sig_mask.any():
        ax.scatter(ks[sig_mask], np.full(sig_mask.sum(),
                                         ax.get_ylim()[1] if False else 0.18),
                   marker="|", color=CLR_SIG, s=18, zorder=5,
                   label=r"$p_{\mathrm{cohort}}<0.05$")

    # Peak annotation: read the actual cohort-paired-Wilcoxon minimum k
    # from the data rather than pre-pinning to a chosen "k=40". This
    # avoids privileging a single k in the figure.
    if annotate_peak and sig_mask.any():
        peak_k = int(ks[int(np.argmin(p_vals))])
        idx = int(np.where(ks == peak_k)[0][0])
        o = obs_med[idx]
        s = surr_med[idx]
        ratio = abs(o) / max(abs(s), 1e-12)
        ax.annotate(
            rf"$k = {peak_k}$  (cohort min-$p$)" "\n"
            rf"ratio $= {ratio:.1f}\times$" "\n"
            rf"Wilcoxon $p = {p_vals[idx]:.3f}$",
            xy=(peak_k, o), xytext=(peak_k - 14, o - 0.18),
            fontsize=8.5, color=CLR_PEAK,
            arrowprops=dict(arrowstyle="->", color=CLR_PEAK, lw=0.8),
            bbox=dict(facecolor="white", edgecolor="#bbb",
                      boxstyle="round,pad=0.30"),
        )

    ax.set_title(title, loc="left", fontsize=10.5, pad=4)
    # Sign convention locked 2026-05-26: T_G > 0 = TRACE.
    ax.set_ylabel(r"$T_G(k) = d_{\mathrm{chord}}^{\mathrm{rsPre, taskT}}"
                  r" - d_{\mathrm{chord}}^{\mathrm{taskT, rsPost}}$",
                  fontsize=10)
    if show_xticks:
        ax.set_xlabel(r"subspace dimension $k$", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(ks.min() - 0.5, ks.max() + 0.5)

    # y-limits: leave headroom for the top-strip tick marks
    ymin = float(np.nanmin([obs_med.min(), env_p25.min(), surr_med.min(),
                            -0.7]))
    ymax = float(np.nanmax([obs_med.max(), env_p75.max(), surr_med.max(),
                            0.15]))
    ax.set_ylim(ymin - 0.10, ymax + 0.20)

    # Per-panel summary text (lower right) — report the count of
    # significant cells, the range over which they spread, and the
    # longest contiguous run as one descriptor of distribution shape
    # rather than as the result itself.
    if n_sig > 0:
        sig_ks = ks[sig_mask]
        summary = (
            rf"sig $k$ cells (Wilcoxon $p<0.05$): {n_sig} / {len(ks)}"
            "\n"
            rf"spread: $k \in [{int(sig_ks.min())}, {int(sig_ks.max())}]$"
            "\n"
            rf"longest contiguous run: $\ell = {run_len}$ "
            rf"(at $k\in[{k_run_lo},{k_run_hi}]$)"
        )
    else:
        summary = rf"sig $k$ cells (Wilcoxon $p<0.05$): 0 / {len(ks)}"
    ax.text(0.99, 0.03, summary, transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8.5,
            bbox=dict(facecolor="white", edgecolor="#bbb",
                      boxstyle="round,pad=0.30"))


def main(band: str = "beta") -> Path:
    band_tex = BRAIN_BAND_TEX_DICT.get(band, rf"${band}$")
    out_dir = ROOT / "data" / "preprint" / "figures" / band / "grassmann"
    out_dir.mkdir(parents=True, exist_ok=True)

    cohort_full = load_cohort(band, GRASS_DIR / "cohort_summary.csv")
    per_pat_full = load_per_patient(band, GRASS_DIR / "per_patient_per_band_per_k.csv")
    env_full = cohort_surrogate_envelope(per_pat_full)

    cohort_epiX = load_cohort(band, GRASS_EPIX_DIR / "cohort_summary.csv")
    per_pat_epiX = load_per_patient(band, GRASS_EPIX_DIR / "per_patient_per_band_per_k.csv")
    env_epiX = cohort_surrogate_envelope(per_pat_epiX)

    fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(11.5, 7.4), sharex=False)

    draw_panel(
        ax_a, cohort_full, env_full,
        label_obs=rf"observed cohort median, {band_tex}",
        clr_obs=CLR_OBS,
        title=rf"(a) {band_tex} Grassmann $d_G(k)$ — full FC vs matched-strength surrogate",
        annotate_peak=True,
        show_xticks=False,
    )

    draw_panel(
        ax_b, cohort_epiX, env_epiX,
        label_obs=rf"observed cohort median, {band_tex} (epi-zone excluded)",
        clr_obs=CLR_OBS_EPIX,
        title=rf"(b) {band_tex} Grassmann $d_G(k)$ — non-epileptic cortex only",
        annotate_peak=False,
        show_xticks=True,
    )

    handles = [
        mlines.Line2D([], [], color=CLR_OBS, lw=1.8,
                      label=r"observed cohort median $T_G(k)$  (full FC)"),
        mlines.Line2D([], [], color=CLR_OBS_EPIX, lw=1.8,
                      label=r"observed cohort median $T_G(k)$  (epi-excluded)"),
        mlines.Line2D([], [], color=CLR_SURR, lw=1.0,
                      label="surrogate cohort median"),
        mpatches.Patch(facecolor=CLR_SURR_FILL,
                       label="surrogate IQR  (cohort-median of per-patient p25/p75)"),
        mpatches.Patch(facecolor=CLR_WIN, alpha=0.30,
                       label=r"per-$k$ Wilcoxon $p_{\mathrm{cohort}} < 0.05$  "
                              r"(every sig cell, no fixed window)"),
        mlines.Line2D([], [], color=CLR_SIG, marker="|", lw=0, markersize=10,
                      label=r"sig tick"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=8.5, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout(rect=(0, 0.06, 1, 1))

    out = out_dir / f"fig_{band}_grassmann.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--band", default="beta",
                        choices=["delta", "theta", "alpha", "beta",
                                 "low_gamma", "high_gamma"])
    args = parser.parse_args()
    main(args.band)
