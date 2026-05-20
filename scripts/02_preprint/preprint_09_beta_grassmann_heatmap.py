#!/usr/bin/env python3
"""Grassmann principal-angle heatmap (mode i × cutoff k) — per band.

Default: loops over all six bands and writes one PDF per band to
``data/preprint/figures/<band>/grassmann/fig_<band>_grassmann_heatmap.pdf``.
Run with one or more band names on the command line to restrict, e.g.

    python preprint_09_beta_grassmann_heatmap.py beta
    python preprint_09_beta_grassmann_heatmap.py alpha beta low_gamma

Per-cell principal-angle decomposition of the matched-strength-controlled
Grassmann probe. The cohort-median per-mode persistence

    Δθ_i(k) = θ_i^{rsPre, taskT} − θ_i^{taskT, rsPost}   (rad)

is rendered as a diverging RdBu_r heatmap on the triangular region i ≤ k
(positive Δθ ⇒ rsPost subspace closer to taskT than rsPre ⇒ trace ⇒ red).

Per-`k` cohort-paired matched-strength Wilcoxon `p<0.05` is drawn as a
two-row sig-bar above the heatmap (full FC + epi-X) and as Wilcoxon-tick
strips below. No 7/10 reference line, no contiguous-window pre-frame —
the cohort gate is per-`k` Wilcoxon p<0.05 only.

Right panel: per-mode profiles at six representative cutoffs k spanning
2..112, plasma-coloured small→large; each profile naturally truncates at
i = k.

Inputs
------
data/audit/grassmann_matched_strength_surrogate/cohort_summary.csv
data/audit/grassmann_epi_exclusion/cohort_summary.csv
data/cache/imcoh_lrg/Pat_NN/<band>_<phase>_lrg_imcoh-abs.npz
  (loaded via workflow.lrg.load_lrg_result)

Outputs (PDF only, no PNG sibling; one per band)
------------------------------------------------
data/preprint/figures/<band>/grassmann/fig_<band>_grassmann_heatmap.pdf
data/preprint/cache/grassmann_principal_angles_<band>_full_k.csv  (cache)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.lrg import load_lrg_result

ROOT = setup_script_env()
use_lrg_style()

# ---------------------------------------------------------------------------
# Centralised label sizing. The lrg_eegfc mplstyle baseline is tuned for
# compact 2x2 / 1x4 publication figures (font.size=9). This figure is a
# single full-page composite; bump every text element by LABEL_SCALE so
# labels remain legible at the intended print/screen size. Edit only this
# block to retune the whole figure — never hardcode `fontsize=N` below.
# ---------------------------------------------------------------------------
LABEL_SCALE = 2.0
plt.rcParams.update({
    "font.size":             9 * LABEL_SCALE,
    "axes.titlesize":        10 * LABEL_SCALE,
    "axes.labelsize":        9 * LABEL_SCALE,
    "xtick.labelsize":       7 * LABEL_SCALE,
    "ytick.labelsize":       7 * LABEL_SCALE,
    "legend.fontsize":       7 * LABEL_SCALE,
    "legend.title_fontsize": 8 * LABEL_SCALE,
    "figure.titlesize":      10 * LABEL_SCALE,
})


ALL_BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
COHORT = [
    "Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
    "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15",
]
PHASES = ("rest_pre", "task_test", "rest_post")

GRASS_DIR = ROOT / "data" / "audit" / "grassmann_matched_strength_surrogate"
GRASS_EPIX_DIR = ROOT / "data" / "audit" / "grassmann_epi_exclusion"
CLUSTER_DIR = ROOT / "data" / "audit" / "grassmann_cluster_extent"

GRASS_PER_PAT_CSV = "per_patient_per_band_per_k.csv"

CACHE_DIR = ROOT / "data" / "preprint" / "cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# Full k axis: matched-strength cohort summary covers k=2..112.
K_MIN, K_MAX = 2, 112
I_MAX = K_MAX

CLR_FULL = "#1d1d1d"
CLR_EPIX = "#7f3b2c"
CLR_SIG = "#1f7a1f"
CLR_SURR = "#7d7d7d"           # T_G surrogate cohort-median line
CLR_SURR_FILL = "#dcdcdc"      # T_G surrogate IQR envelope fill
CLR_TG_OBS = "#1f3d6e"         # T_G observed (full FC)
CLR_TG_EPIX = CLR_EPIX         # T_G observed (epi-X)

# Cluster-extent verdict colour scheme (locked 2026-05-19, VERDICT_LEDGER.md).
CLR_VERDICT_STRONG = "#0d5c1c"     # deep green: cluster_p < 0.01
CLR_VERDICT_WEAK = "#a48b22"        # ochre:     0.01 ≤ cluster_p < 0.05
CLR_VERDICT_NONE = "#9a9a9a"        # grey:      cluster_p ≥ 0.05

VERDICT_COLOR = {
    "strong": CLR_VERDICT_STRONG,
    "weak": CLR_VERDICT_WEAK,
    "no_trace": CLR_VERDICT_NONE,
}

K_PROFILES = [10, 30, 50, 70, 90, 110]


# ---------------------------------------------------------------------------
# Compute helpers
# ---------------------------------------------------------------------------
def _topk_basis(eigvecs: np.ndarray, k: int) -> np.ndarray:
    return np.ascontiguousarray(eigvecs[:, 1 : k + 1])


def _principal_angles(V_a: np.ndarray, V_b: np.ndarray) -> np.ndarray:
    sigma = np.linalg.svd(V_a.T @ V_b, compute_uv=False)
    sigma = np.clip(sigma, 0.0, 1.0)
    return np.arccos(sigma)


def compute_principal_angles_full_k(band: str) -> pd.DataFrame:
    """Cohort principal-angle Δθ_i(k) over k=2..112 for the given band.

    Convention: Δθ_i = θ_i^{rsPre,taskT} − θ_i^{taskT,rsPost}. Positive
    ⇒ rsPost subspace sits closer to taskT than rsPre does along mode i
    ⇒ trace direction.
    """
    cache = CACHE_DIR / f"grassmann_principal_angles_{band}_full_k.csv"
    if cache.exists():
        return pd.read_csv(cache)

    per_pat = []
    for pat in COHORT:
        EV = {}
        for ph in PHASES:
            res = load_lrg_result(pat, ph, band, fc_method="imcoh_abs")
            EV[ph] = res.eigenvectors
        N = EV[PHASES[0]].shape[0]
        k_hi = min(K_MAX, N - 1)
        for k in range(K_MIN, k_hi + 1):
            Vpre = _topk_basis(EV["rest_pre"], k)
            Vtt = _topk_basis(EV["task_test"], k)
            Vpost = _topk_basis(EV["rest_post"], k)
            th_pre_tt = _principal_angles(Vpre, Vtt)
            th_tt_post = _principal_angles(Vtt, Vpost)
            for i, (a, b) in enumerate(zip(th_pre_tt, th_tt_post), start=1):
                per_pat.append(dict(
                    patient=pat, k=k, mode_index=i,
                    delta_theta=float(a - b),
                ))
    pp = pd.DataFrame(per_pat)
    rows = []
    for (k, i), g in pp.groupby(["k", "mode_index"]):
        vals = g["delta_theta"].values
        rows.append(dict(
            k=int(k),
            mode_index=int(i),
            n_patients=int(g["patient"].nunique()),
            delta_theta_median_radians=float(np.median(vals)),
            delta_theta_iqr_radians=float(
                np.subtract(*np.percentile(vals, [75, 25]))),
        ))
    cohort = pd.DataFrame(rows).sort_values(["k", "mode_index"])
    cohort.to_csv(cache, index=False)
    print(f"  [{band}] cached principal-angle decomposition → {cache.name}")
    return cohort


def build_heatmap(pa: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    ks = list(range(K_MIN, K_MAX + 1))
    H = np.full((I_MAX, len(ks)), np.nan)
    for _, row in pa.iterrows():
        k = int(row.k)
        i = int(row.mode_index) - 1
        if k < K_MIN or k > K_MAX or i >= I_MAX or i < 0:
            continue
        col = ks.index(k)
        H[i, col] = float(row.delta_theta_median_radians)
    return H, np.array(ks)


def cohort_strip(ms: pd.DataFrame, band: str
                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sub = ms[(ms.band == band) & (ms.k >= K_MIN) & (ms.k <= K_MAX)]
    sub = sub.sort_values("k")
    return (sub.k.values.astype(int),
            sub.n_patients_below_own_surrogate.values.astype(float) / 10.0,
            sub.paired_wilcoxon_p.values.astype(float))


def longest_run_span(ks: np.ndarray, p_vals: np.ndarray
                     ) -> tuple[int, int, int]:
    """Return (length, k_start, k_end) of the longest contiguous p<0.05 run."""
    sig = (p_vals < 0.05).astype(int)
    best_len = 0
    best_s = -1
    best_e = -1
    cur_len = 0
    cur_s = -1
    for idx, s in enumerate(sig):
        if s:
            if cur_len == 0:
                cur_s = idx
            cur_len += 1
            if cur_len > best_len:
                best_len = cur_len
                best_s = cur_s
                best_e = idx
        else:
            cur_len = 0
    if best_len == 0:
        return 0, -1, -1
    return best_len, int(ks[best_s]), int(ks[best_e])


def cohort_surrogate_envelope(per_pat: pd.DataFrame, band: str
                               ) -> pd.DataFrame:
    """Per-k cohort-median of per-patient surrogate quantiles for T_G.

    Returns columns: k, surr_p25, surr_p50, surr_p75.
    """
    sub = per_pat[per_pat.band == band]
    rows = []
    for k, g in sub.groupby("k"):
        rows.append(dict(
            k=int(k),
            surr_p25=float(np.nanmedian(g.surr_T_G_p25.values)),
            surr_p50=float(np.nanmedian(g.surr_T_G_p50.values)),
            surr_p75=float(np.nanmedian(g.surr_T_G_p75.values)),
        ))
    return pd.DataFrame(rows).sort_values("k").reset_index(drop=True)


def cluster_verdict_row(band: str) -> dict:
    """Read the locked cluster-extent verdict for *band* from audit_70."""
    df = pd.read_csv(CLUSTER_DIR / "cohort_summary.csv")
    row = df[df.band == band]
    if row.empty:
        raise RuntimeError(f"no cluster-extent row for band={band}")
    r = row.iloc[0]
    return dict(
        verdict=str(r.verdict_cluster_extent),
        cluster_p_LR=float(r.cluster_p_longest_run),
        cluster_p_mass=float(r.cluster_p_cluster_mass),
        obs_LR=int(r.obs_longest_run),
        null_mean=float(r.null_mean_LR),
        null_p95=float(r.null_p95_LR),
        null_max=float(r.null_max_LR),
    )


# ---------------------------------------------------------------------------
# Plot driver
# ---------------------------------------------------------------------------
def build_figure(band: str, vmax_global: float | None = None,
                 visual_only: bool = False) -> Path:
    """Render the per-band Grassmann heatmap figure.

    visual_only=False (default): full annotated figure with titles,
    verdict box, legends, reading note. Saved as
    `fig_<band>_grassmann_heatmap.pdf`.

    visual_only=True: same plot elements (heatmap, sig-bar with
    verdict-coloured outline, bottom strips, right per-mode profiles,
    colorbar) but every title / label / legend / text-box stripped.
    Saved as `fig_<band>_grassmann_heatmap_visual.pdf`.
    """
    band_tex = BRAIN_BAND_TEX_DICT[band]
    out_dir = ROOT / "data" / "preprint" / "figures" / band / "grassmann"
    out_dir.mkdir(parents=True, exist_ok=True)

    pa = compute_principal_angles_full_k(band)
    ms_full = pd.read_csv(GRASS_DIR / "cohort_summary.csv")
    ms_epix = pd.read_csv(GRASS_EPIX_DIR / "cohort_summary.csv")

    H, ks = build_heatmap(pa)
    ks_full, n_full, p_full = cohort_strip(ms_full, band)
    ks_epix, n_epix, p_epix = cohort_strip(ms_epix, band)

    # Locked cluster-extent verdict for this band (VERDICT_LEDGER.md).
    verdict = cluster_verdict_row(band)
    run_len, run_k_lo, run_k_hi = longest_run_span(ks_full, p_full)
    verdict_clr = VERDICT_COLOR[verdict["verdict"]]

    finite = H[np.isfinite(H)]
    if vmax_global is not None:
        vmax = vmax_global
    else:
        vmax = float(np.percentile(np.abs(finite), 95))
        vmax = max(vmax, 0.05)

    fig = plt.figure(figsize=(15.0, 11.5))
    # 3-column × 3-row layout:
    #   row 0 (heatmap):    ax_h  | ax_p  | cax
    #   row 1 (T_G line):   ax_tg | ax_legend_tg     | (blank)
    #   row 2 (n_<surr):    ax_n  | ax_legend_nsurr  | (blank)
    # Colorbar lives in its own column so its height matches ax_h.
    # The right column below the colorbar is white space — used for
    # the T_G and n_<surr legends rather than figure-level legends.
    outer = GridSpec(
        3, 3, figure=fig,
        width_ratios=[3.4, 1.0, 0.05],
        height_ratios=[3.6, 1.0, 0.7],
        wspace=0.10, hspace=0.22,
        left=0.06, right=0.95, top=0.95, bottom=0.07,
    )
    ax_h = fig.add_subplot(outer[0, 0])
    ax_p = fig.add_subplot(outer[0, 1])
    cax = fig.add_subplot(outer[0, 2])
    ax_tg = fig.add_subplot(outer[1, 0], sharex=ax_h)
    ax_legend_tg = fig.add_subplot(outer[1, 1])
    ax_legend_tg.axis("off")
    ax_n = fig.add_subplot(outer[2, 0], sharex=ax_h)
    ax_legend_n = fig.add_subplot(outer[2, 1])
    ax_legend_n.axis("off")

    # Pre-compute per-k significance masks for the overlay (drawn AFTER
    # the heatmap so the stripes sit above the data).
    sig_full_ks = np.array([int(kk) for kk, pp in zip(ks_full, p_full)
                             if K_MIN <= kk <= K_MAX and pp < 0.05])
    sig_epix_ks = np.array([int(kk) for kk, pp in zip(ks_epix, p_epix)
                             if K_MIN <= kk <= K_MAX and pp < 0.05])
    if len(ks_epix):
        epix_k_max = int(ks_epix.max())
    else:
        epix_k_max = K_MIN - 1

    # ---- heatmap -----------------------------------------------------------
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#f2f2f2")
    H_masked = np.ma.masked_invalid(H)
    im = ax_h.imshow(
        H_masked, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax,
        origin="upper", interpolation="nearest",
        extent=(ks.min() - 0.5, ks.max() + 0.5, I_MAX + 0.5, 0.5),
    )
    ax_h.plot([ks.min() - 0.5, ks.max() + 0.5],
              [ks.min() - 0.5, ks.max() + 0.5],
              color="0.40", lw=0.7, ls="--", zorder=4)

    # ---- per-k Wilcoxon p<0.05 overlay above heatmap data -----------------
    # Extend the y-axis upward so the overlay band sits ABOVE the data
    # rather than on top of i=1..3 cells. Two lanes in data coords:
    # full FC (top) + epi-X (just below). Y-ticks restricted to i ≥ 1
    # so the overlay region is "out of frame" in tick numbering.
    OVERLAY_HEIGHT = 5.0
    FULL_Y0, FULL_Y1 = -4.7, -2.7
    EPIX_Y0, EPIX_Y1 = -2.5, -0.5

    from matplotlib.patches import Rectangle

    # Background lanes (light grey).
    ax_h.add_patch(Rectangle(
        (K_MIN - 0.5, FULL_Y0), K_MAX - K_MIN + 1, FULL_Y1 - FULL_Y0,
        facecolor="#ededed", edgecolor="none", zorder=8))
    if epix_k_max >= K_MIN:
        ax_h.add_patch(Rectangle(
            (K_MIN - 0.5, EPIX_Y0), epix_k_max - K_MIN + 1,
            EPIX_Y1 - EPIX_Y0,
            facecolor="#ededed", edgecolor="none", zorder=8))

    # Sig-cell stripes.
    for k_sig in sig_full_ks:
        ax_h.add_patch(Rectangle(
            (k_sig - 0.5, FULL_Y0), 1.0, FULL_Y1 - FULL_Y0,
            facecolor=CLR_SIG, edgecolor="none", zorder=9))
    for k_sig in sig_epix_ks:
        ax_h.add_patch(Rectangle(
            (k_sig - 0.5, EPIX_Y0), 1.0, EPIX_Y1 - EPIX_Y0,
            facecolor=CLR_SIG, edgecolor="none", zorder=9))

    # Longest-contiguous-run outline (verdict colour).
    if run_len > 0:
        ax_h.add_patch(Rectangle(
            (run_k_lo - 0.5, FULL_Y0),
            run_k_hi - run_k_lo + 1, FULL_Y1 - FULL_Y0,
            fill=False, edgecolor=verdict_clr, lw=2.0, zorder=10))

    # Row labels for the overlay (data coords, in the left margin).
    if not visual_only:
        ax_h.text(K_MIN - 1.8, (FULL_Y0 + FULL_Y1) / 2,
                   "full FC", ha="right", va="center", color="0.30")
        ax_h.text(K_MIN - 1.8, (EPIX_Y0 + EPIX_Y1) / 2,
                   "epi-X", ha="right", va="center", color="0.30")

    if not visual_only:
        ax_h.text(K_MAX * 0.97, K_MAX * 0.97 - 1.8, r"$i = k$",
                  color="0.30", ha="right", va="bottom",
                  rotation=-45)
        ax_h.text(0.92, 0.92,
                  "grey: $i > k$\nno principal angle defined",
                  transform=ax_h.transAxes, ha="right", va="top",
                  color="0.30",
                  bbox=dict(facecolor="white", edgecolor="#d8d8d8",
                            boxstyle="round,pad=0.30", alpha=0.92))

        # ---- locked cluster-extent verdict box (top-right of heatmap) -----
        verdict_label = {
            "strong": "STRONG TRACE",
            "weak": "WEAK TRACE",
            "no_trace": "NO TRACE",
        }[verdict["verdict"]]
        if run_len > 0:
            run_span_text = (rf"$\ell_{{\mathrm{{obs}}}} = {verdict['obs_LR']}$  "
                              rf"at $k \in [{run_k_lo}, {run_k_hi}]$")
        else:
            run_span_text = r"$\ell_{\mathrm{obs}} = 0$"
        verdict_text = (
            rf"$\bf{{Cluster\text{{-}}extent\ verdict}}$:  {verdict_label}"
            "\n"
            rf"{run_span_text}"
            "\n"
            rf"null $\bar{{\ell}} = {verdict['null_mean']:.2f}$,  "
            rf"$\ell^{{(95)}} = {verdict['null_p95']:.1f}$,  "
            rf"$\ell^{{\max}} = {verdict['null_max']:.0f}$"
            "\n"
            rf"cluster $p_{{\ell}} = {verdict['cluster_p_LR']:.4f}$"
        )
        ax_h.text(
            0.02, 0.98, verdict_text,
            transform=ax_h.transAxes, ha="left", va="top",
            color=verdict_clr,
            bbox=dict(facecolor="white", edgecolor=verdict_clr,
                      boxstyle="round,pad=0.45", linewidth=1.6, alpha=0.97),
            zorder=10,
        )

    # Extend y-axis upward to hold the overlay band above the data.
    ax_h.set_ylim(I_MAX + 0.5, -OVERLAY_HEIGHT)

    x_ticks = [2, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 112]
    y_ticks = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 112]
    ax_h.set_xticks(x_ticks)
    ax_h.set_xticklabels([str(t) for t in x_ticks])
    ax_h.set_yticks(y_ticks)
    ax_h.set_yticklabels([str(t) for t in y_ticks])
    if visual_only:
        ax_h.set_ylabel(r"$i$")
    else:
        ax_h.set_ylabel(r"mode index  $i$  (1 = slowest non-trivial mode)")
        ax_h.set_title(
            rf"(a) {band_tex} Grassmann — cohort-median per-mode persistence"
            rf"  $\Delta\theta_i(k) = \theta_i^{{\mathrm{{rsPre,taskT}}}}"
            rf" - \theta_i^{{\mathrm{{taskT,rsPost}}}}$ (rad).  "
            r"Red cell $\Leftrightarrow$ rsPost closer to taskT than rsPre is, along mode $i$.",
            loc="left", pad=4)
    plt.setp(ax_h.get_xticklabels(), visible=False)
    ax_h.spines[["top", "right"]].set_visible(False)

    # ---- shared colorbar (GridSpec cell — height matches ax_h) ------------
    cb = fig.colorbar(im, cax=cax)
    if visual_only:
        cb.set_label(r"$\Delta\theta_i(k)$  (rad)")
    else:
        cb.set_label(
            r"$\Delta\theta_i(k)$  (rad)" "\n"
            r"red $> 0$: rsPost closer to taskT than rsPre is" "\n"
            r"            $\Leftrightarrow$ trace direction along mode $i$" "\n"
            r"blue $< 0$: anti direction (rsPre closer to taskT)")

    # ---- T_G(k) panel: observed cohort median vs matched-strength ---------
    # Load per-patient T_G to build the surrogate IQR envelope (median
    # over patients of each patient's own p25/p50/p75 across R=200
    # surrogate replicates). Plot observed cohort-median T_G(k) for both
    # full FC and epi-X overlaid; shade the surrogate IQR (full FC) so
    # the cohort signal vs strength-controlled null is visible by eye.
    per_pat_full = pd.read_csv(GRASS_DIR / GRASS_PER_PAT_CSV)
    per_pat_epix = pd.read_csv(GRASS_EPIX_DIR / GRASS_PER_PAT_CSV)
    env_full = cohort_surrogate_envelope(per_pat_full, band)
    env_epix = cohort_surrogate_envelope(per_pat_epix, band)

    cohort_full_df = ms_full[ms_full.band == band].sort_values("k")
    cohort_epix_df = ms_epix[ms_epix.band == band].sort_values("k")

    tg_obs_full = cohort_full_df.obs_median_T_G.values.astype(float)
    tg_obs_epix = cohort_epix_df.obs_median_T_G.values.astype(float)
    tg_surr_full = cohort_full_df.surr_median_T_G_per_patient_median.values.astype(float)

    h_iqr_full = ax_tg.fill_between(
        env_full.k.values, env_full.surr_p25.values, env_full.surr_p75.values,
        color=CLR_SURR_FILL, alpha=0.95, zorder=1)
    h_iqr_epix = ax_tg.fill_between(
        env_epix.k.values, env_epix.surr_p25.values, env_epix.surr_p75.values,
        color=CLR_EPIX, alpha=0.12, zorder=1)
    h_surr_full, = ax_tg.plot(ks_full, tg_surr_full,
                               color=CLR_SURR, lw=1.0, zorder=2)
    h_obs_full, = ax_tg.plot(ks_full, tg_obs_full,
                              color=CLR_TG_OBS, lw=1.6, zorder=4)
    h_obs_epix, = ax_tg.plot(ks_epix, tg_obs_epix,
                              color=CLR_TG_EPIX, lw=1.6, zorder=4)
    ax_tg.axhline(0, color="0.4", lw=0.6, ls="--", zorder=1)
    ax_tg.set_xlim(ks.min() - 0.5, ks.max() + 0.5)
    ax_tg.set_xticks(x_ticks)
    ax_tg.set_xticklabels([str(t) for t in x_ticks])
    ax_tg.set_ylabel(r"$T_G(k)$")
    ax_tg.spines[["top", "right"]].set_visible(False)
    plt.setp(ax_tg.get_xticklabels(), visible=False)

    # Right-of-axes legend for the T_G panel.
    ax_legend_tg.legend(
        handles=[
            mlines.Line2D([], [], color=CLR_TG_OBS, lw=1.8,
                          label=r"observed cohort median (full FC)"),
            mlines.Line2D([], [], color=CLR_TG_EPIX, lw=1.8,
                          label=r"observed cohort median (epi-X)"),
            mlines.Line2D([], [], color=CLR_SURR, lw=1.0,
                          label="surrogate cohort median"),
            mpatches.Patch(facecolor=CLR_SURR_FILL,
                           label="surrogate IQR (full FC)"),
            mpatches.Patch(facecolor=CLR_EPIX, alpha=0.18,
                           label="surrogate IQR (epi-X)"),
        ],
        loc="center left", frameon=False, handlelength=2.0,
    )

    # ---- combined n_<surr panel: both probes on a single axes -------------
    h_n_full, = ax_n.plot(ks_full, n_full, color=CLR_FULL, lw=1.4, zorder=3)
    h_n_epix, = ax_n.plot(ks_epix, n_epix, color=CLR_EPIX, lw=1.4, zorder=3)
    # Wilcoxon-sig ticks at the top of the strip (full FC only — that's
    # the locked verdict layer).
    sig_full_mask = p_full < 0.05
    if sig_full_mask.any():
        ax_n.scatter(ks_full[sig_full_mask],
                      np.full(sig_full_mask.sum(), 1.06),
                      marker="|", color=CLR_SIG, s=22, zorder=4)
    ax_n.set_xlim(ks.min() - 0.5, ks.max() + 0.5)
    ax_n.set_ylim(-0.05, 1.18)
    ax_n.set_yticks([0.0, 0.5, 1.0])
    ax_n.set_yticklabels(["0", "5/10", "10/10"])
    ax_n.set_xticks(x_ticks)
    ax_n.set_xticklabels([str(t) for t in x_ticks])
    ax_n.set_ylabel(r"$n_{<\mathrm{surr}}$")
    ax_n.set_xlabel(r"$k$")
    ax_n.spines[["top", "right"]].set_visible(False)

    # Right-of-axes legend for the n_<surr panel.
    ax_legend_n.legend(
        handles=[
            mlines.Line2D([], [], color=CLR_FULL, lw=1.6, label="full FC"),
            mlines.Line2D([], [], color=CLR_EPIX, lw=1.6,
                          label="epi-zone excluded"),
            mlines.Line2D([], [], color=CLR_SIG, marker="|", lw=0,
                          markersize=14,
                          label=r"Wilcoxon $p_{\mathrm{cohort}}<0.05$"),
        ],
        loc="center left", frameon=False, handlelength=2.0,
    )

    # ---- right-side companion: per-mode profiles at representative k -----
    cmap_profiles = plt.get_cmap("plasma")
    norm_profiles = matplotlib.colors.Normalize(
        vmin=min(K_PROFILES), vmax=max(K_PROFILES))

    def _profile(ax, k, color):
        col = pa[pa.k == k].sort_values("mode_index")
        if col.empty:
            return
        y = col.mode_index.values.astype(int)
        dth = col.delta_theta_median_radians.values.astype(float)
        ax.plot(dth, y, color=color, lw=1.2, zorder=2,
                label=f"{k}")

    for k in K_PROFILES:
        clr = cmap_profiles(norm_profiles(k))
        _profile(ax_p, k, clr)

    ax_p.axvline(0, color="0.4", lw=0.6, ls="--", zorder=1)
    ax_p.invert_yaxis()
    ax_p.set_ylim(I_MAX + 0.5, 0.5)
    ax_p.set_yticks(y_ticks)
    ax_p.set_yticklabels([])
    ax_p.spines[["top", "right"]].set_visible(False)
    ax_p.set_xlabel(r"$\Delta\theta_i(k)$  (rad)")
    ax_p.legend(loc="lower right", frameon=False,
                 ncol=2, handlelength=1.4, title=r"$k$")
    if not visual_only:
        ax_p.set_title(r"(b) per-mode profiles  $\Delta\theta_i(k)$  at "
                        r"representative $k$",
                        loc="left", pad=4)
        ax_p.text(0.02, 0.02,
                  "each profile truncates at $i = k$\n"
                  "(no principal angle beyond)",
                  transform=ax_p.transAxes, ha="left", va="bottom",
                  color="0.35")

    # The per-panel legends (T_G, n_<surr) live in ax_legend_tg /
    # ax_legend_n at the right of the figure; no figure-level legend in
    # this layout. The verdict box on the heatmap (annotated mode) and
    # the verdict-coloured outline rectangle on the top-of-heatmap
    # overlay together carry the cluster-extent reading.

    suffix = "_visual" if visual_only else ""
    out = out_dir / f"fig_{band}_grassmann_heatmap{suffix}.pdf"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"  [{band}{' visual' if visual_only else ''}] saved → {out}")
    return out


def main() -> list[Path]:
    bands = sys.argv[1:] if len(sys.argv) > 1 else ALL_BANDS
    bad = [b for b in bands if b not in ALL_BANDS]
    if bad:
        raise SystemExit(
            f"unknown band(s): {bad}; choose from {ALL_BANDS}")

    # Precompute principal-angle caches first so we can use a common
    # v-max across bands if desired. For now keep per-band v-max so the
    # weakest bands aren't washed out by β's stronger signal.
    outs = []
    for b in bands:
        outs.append(build_figure(b, visual_only=False))
        outs.append(build_figure(b, visual_only=True))
    return outs


if __name__ == "__main__":
    main()
