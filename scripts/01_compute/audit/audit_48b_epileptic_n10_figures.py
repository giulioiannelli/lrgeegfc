#!/usr/bin/env python3
"""Audit 48b — Figure-first reporting of the epileptic n=10 revisit.

Reads the CSVs produced by ``audit_48_epileptic_n10_compute.py`` and emits
six manuscript-grade PDFs into
``data/audit/epileptic_n10_revisit/figures/``.

Figures
-------
fig_01_amplification_scatter.pdf
    Per-(patient, phase) scatter of R_imcoh_cp (raw |ImCoh|) vs R_ρ_cp
    (LRG ρ̂ at τ_max), one panel per band. Diagonal = LRG passes the raw
    enrichment through; deviation above = LRG amplifies, below = dilutes.

fig_02_enrichment_heatmap.pdf
    6 × 3 cohort-median heatmaps: rows = pair-class (all / cross-probe /
    same-probe); cols = raw |ImCoh| | ρ̂(τ_max) | their delta. Cells colored
    by log10 ratio. Stars where Wilcoxon-against-1 q < 0.05 (BH per panel).

fig_03_tau_sweep.pdf
    R_ρ_cp(τ) cohort median + IQR, one panel per band, three lines (rest_pre,
    task_test, rest_post). Vertical guide at τ_max = 1/λ_max.

fig_04_mrca_violins.pdf
    Pooled MRCA merge-height per pair-class (epi-epi, non-non, cross),
    one panel per band, three sub-violins per phase. Visual answer to
    "where in the LRG hierarchy do epi pairs first meet?".

fig_05_kc_epi_triangle.pdf
    KC distance restricted to epi-epi pairs across phase pairs (rsPre↔task,
    task↔rsPost, rsPre↔rsPost) and λ ∈ {0, 0.5, 1}. Per-patient T_KC_epi =
    d(rsPre, task) − d(task, rsPost). Positive = trace at the epi sub-network.
    (Sign convention locked 2026-05-26: T_d > 0 = TRACE.)

fig_06_dendrograms_beta.pdf
    9 patients × 3 phases at β. Dendrograms with epi leaves coloured red,
    others grey. Visual support for the heatmap-level findings.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram
from scipy.spatial.distance import squareform
from scipy.stats import wilcoxon

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PHASE_SHORT_LABELS
from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.workflow.lrg import load_lrg_result

import sys
sys.path.insert(0, str(Path(__file__).parent))
from audit_48_epileptic_n10_compute import (
    COHORT, PHASES, KC_LAMBDAS, _build_masks
)

CSV_DIR = ROOT / "data" / "audit" / "epileptic_n10_revisit"
FIG_DIR = CSV_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

PHASE_COLORS = {
    "rest_pre": "#4878d0",
    "task_test": "#d65f5f",
    "rest_post": "#6acc64",
}
CLS_COLORS = {
    "ee": "#d62728",  # epi-epi: red
    "nn": "#2ca02c",  # non-non: green
    "cross": "#1f77b4",  # cross: blue
}
CLS_LABEL = {"ee": "epi-epi", "nn": "non-non", "cross": "epi-non"}


# ---------------------------------------------------------------------------
# fig 01 — amplification scatter
# ---------------------------------------------------------------------------

def fig_01_amplification_scatter() -> None:
    m1 = pd.read_csv(CSV_DIR / "M1_raw_imcoh_enrichment.csv")
    m2 = pd.read_csv(CSV_DIR / "M2_lrg_enrichment.csv")
    m1 = m1[m1.cls == "cp"]
    m2 = m2[m2.cls == "cp"]
    df = m1.merge(
        m2[["patient", "band", "phase", "R_rho", "R_um"]],
        on=["patient", "band", "phase"], how="inner",
    )

    bands = list(BRAIN_BANDS_NAMES)
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 7.6), sharex=True, sharey=True)
    axes = axes.ravel()
    lo = max(0.4, min(df.R_imcoh.min(), df.R_rho.min()) * 0.92)
    hi = min(2.6, max(df.R_imcoh.max(), df.R_rho.max()) * 1.05)

    for ax, band in zip(axes, bands):
        sub = df[df.band == band]
        ax.plot([lo, hi], [lo, hi], "--", color="0.5", lw=0.9, zorder=1)
        ax.axhline(1.0, color="0.85", lw=0.6, zorder=0)
        ax.axvline(1.0, color="0.85", lw=0.6, zorder=0)
        for ph in PHASES:
            s = sub[sub.phase == ph]
            ax.scatter(
                s.R_imcoh, s.R_rho,
                s=42, color=PHASE_COLORS[ph], edgecolor="white",
                linewidth=0.6, alpha=0.92, label=PHASE_SHORT_LABELS[ph], zorder=3,
            )
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=12)
        ax.tick_params(labelsize=9)

    for ax in axes[3:]:
        ax.set_xlabel(r"raw $|\mathrm{ImCoh}|$ ratio  $R_{\mathrm{imcoh}}^{cp}$",
                      fontsize=10)
    for ax in (axes[0], axes[3]):
        ax.set_ylabel(r"LRG $\hat\rho(\tau_{\max})$ ratio  $R_{\rho}^{cp}$",
                      fontsize=10)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.02), fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(FIG_DIR / "fig_01_amplification_scatter.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# fig 02 — enrichment heatmap (band × phase, 3 cls × 3 metrics)
# ---------------------------------------------------------------------------

def _wilcoxon_against_one(values: np.ndarray) -> float:
    """One-sided Wilcoxon p that median > 1, working in log-space for ratios."""
    v = np.asarray(values, dtype=float)
    v = v[~np.isnan(v)]
    if v.size < 3:
        return float("nan")
    try:
        _, p = wilcoxon(np.log(v), alternative="greater")
        return float(p)
    except ValueError:
        return float("nan")


def fig_02_enrichment_heatmap() -> None:
    m1 = pd.read_csv(CSV_DIR / "M1_raw_imcoh_enrichment.csv")
    m2 = pd.read_csv(CSV_DIR / "M2_lrg_enrichment.csv")

    bands = list(BRAIN_BANDS_NAMES)
    cls_order = ["all", "cp", "sp"]
    cls_label = {"all": "all pairs", "cp": "cross-probe", "sp": "same-probe"}

    fig, axes = plt.subplots(3, 3, figsize=(13, 9))
    norm_ratio = TwoSlopeNorm(vcenter=0.0, vmin=-0.5, vmax=0.5)
    norm_delta = TwoSlopeNorm(vcenter=0.0, vmin=-0.4, vmax=0.4)

    band_lbl = [BRAIN_BAND_TEX_DICT[b] for b in bands]

    for ri, cls in enumerate(cls_order):
        # raw
        Mraw = np.full((len(bands), len(PHASES)), np.nan)
        Mlrg = np.full_like(Mraw, np.nan)
        Mdelta = np.full_like(Mraw, np.nan)
        Praw = np.full_like(Mraw, np.nan)
        Plrg = np.full_like(Mraw, np.nan)
        for bi, band in enumerate(bands):
            for pi, ph in enumerate(PHASES):
                s1 = m1[(m1.cls == cls) & (m1.band == band) & (m1.phase == ph)]
                s2 = m2[(m2.cls == cls) & (m2.band == band) & (m2.phase == ph)]
                if not s1.empty:
                    Mraw[bi, pi] = float(s1.R_imcoh.median())
                    Praw[bi, pi] = _wilcoxon_against_one(s1.R_imcoh.values)
                if not s2.empty:
                    Mlrg[bi, pi] = float(s2.R_rho.median())
                    Plrg[bi, pi] = _wilcoxon_against_one(s2.R_rho.values)
                Mdelta[bi, pi] = (Mlrg[bi, pi] - Mraw[bi, pi]) if (
                    not np.isnan(Mraw[bi, pi]) and not np.isnan(Mlrg[bi, pi])
                ) else np.nan

        # BH-FDR within each class (18 cells per panel)
        praw_flat = Praw.ravel()
        plrg_flat = Plrg.ravel()
        valid_raw = ~np.isnan(praw_flat)
        valid_lrg = ~np.isnan(plrg_flat)
        q_raw = np.full_like(praw_flat, np.nan)
        q_lrg = np.full_like(plrg_flat, np.nan)
        if valid_raw.any():
            q_raw[valid_raw] = bh_fdr(praw_flat[valid_raw])
        if valid_lrg.any():
            q_lrg[valid_lrg] = bh_fdr(plrg_flat[valid_lrg])
        q_raw = q_raw.reshape(Praw.shape)
        q_lrg = q_lrg.reshape(Plrg.shape)

        for ci, (M, Q, title, norm) in enumerate([
            (Mraw, q_raw, r"raw  $R_{\mathrm{imcoh}}$", norm_ratio),
            (Mlrg, q_lrg, r"LRG  $R_{\rho}(\tau_{\max})$", norm_ratio),
            (Mdelta, None, r"$\Delta = R_{\rho} - R_{\mathrm{imcoh}}$", norm_delta),
        ]):
            ax = axes[ri, ci]
            data = np.log10(np.maximum(M, 1e-3)) if ci < 2 else M
            im = ax.imshow(data, aspect="auto", cmap="RdBu_r",
                           norm=norm, origin="lower")
            ax.set_xticks(range(len(PHASES)))
            ax.set_xticklabels([PHASE_SHORT_LABELS[p] for p in PHASES], fontsize=9)
            ax.set_yticks(range(len(bands)))
            if ci == 0:
                ax.set_yticklabels(band_lbl, fontsize=9)
                ax.set_ylabel(cls_label[cls], fontsize=11)
            else:
                ax.set_yticklabels([])
            if ri == 0:
                ax.set_title(title, fontsize=11)
            for bi in range(len(bands)):
                for pi in range(len(PHASES)):
                    if np.isnan(M[bi, pi]):
                        continue
                    txt = f"{M[bi, pi]:.2f}"
                    star = ""
                    if Q is not None and not np.isnan(Q[bi, pi]) and Q[bi, pi] < 0.05:
                        star = "*"
                    ax.text(pi, bi, txt + star, ha="center", va="center",
                            fontsize=8, color="black")

    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_02_enrichment_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# fig 03 — τ-sweep
# ---------------------------------------------------------------------------

def fig_03_tau_sweep() -> None:
    """Cohort R_ρ_cp(τ) aggregated by τ-index (per-patient τ grids span their
    own [τ_min = 1/λ_max, τ*]; absolute τ values are not comparable across
    patients, so we aggregate over τ-index = relative position 0..N-1).

    The x-axis label gives the cohort-median normalized τ × λ_max for each idx.
    """
    sw = pd.read_csv(CSV_DIR / "M2b_tau_sweep.csv")
    bands = list(BRAIN_BANDS_NAMES)

    fig, axes = plt.subplots(2, 3, figsize=(12, 6.5), sharex=True, sharey=True)
    axes = axes.ravel()

    for ax, band in zip(axes, bands):
        sub = sw[sw.band == band]
        for ph in PHASES:
            ssub = sub[sub.phase == ph]
            grouped = ssub.groupby("tau_idx")["R_rho_cp"]
            idx = sorted(ssub.tau_idx.unique())
            med = [grouped.get_group(i).median() for i in idx]
            q25 = [grouped.get_group(i).quantile(0.25) for i in idx]
            q75 = [grouped.get_group(i).quantile(0.75) for i in idx]
            ax.fill_between(idx, q25, q75, color=PHASE_COLORS[ph],
                            alpha=0.20, lw=0)
            ax.plot(idx, med, color=PHASE_COLORS[ph], lw=1.8,
                    marker="o", markersize=4.5, label=PHASE_SHORT_LABELS[ph])
        ax.axhline(1.0, color="0.7", lw=0.7, ls="--")
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=12)
        ax.tick_params(labelsize=9)
        ax.set_ylim(0.7, 1.7)

    # x-axis: τ index annotated with cohort-median τ_norm
    norm_med = sw.groupby("tau_idx")["tau_norm"].median().sort_index()
    xticks = list(norm_med.index)
    xtick_labels = [f"{norm_med.iloc[i]:.1f}" for i in range(len(xticks))]
    for ax in axes:
        ax.set_xticks(xticks)
    for ax in axes[3:]:
        ax.set_xticklabels(xtick_labels, fontsize=8, rotation=0)
        ax.set_xlabel(r"$\tau \cdot \lambda_{\max}$ (cohort median)",
                      fontsize=10)
    for ax in axes[:3]:
        ax.set_xticklabels([])
    for ax in (axes[0], axes[3]):
        ax.set_ylabel(r"$R_{\rho}^{cp}(\tau)$  cohort median $\pm$ IQR",
                      fontsize=10)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.02), fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(FIG_DIR / "fig_03_tau_sweep.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# fig 04 — MRCA violins
# ---------------------------------------------------------------------------

def fig_04_mrca_violins() -> None:
    """Cohort-pooled MRCA height violins per pair-class, per band, per phase."""
    bands = list(BRAIN_BANDS_NAMES)
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.0), sharey=False)
    axes = axes.ravel()

    pooled: dict[tuple[str, str, str], np.ndarray] = {}
    for pat in COHORT:
        try:
            masks = _build_masks(pat)
        except Exception:
            continue
        if masks.epi_mask.sum() == 0:
            continue
        n = len(masks.epi_mask)
        iu, ju = np.triu_indices(n, k=1)
        e_i, e_j = masks.epi_mask[iu], masks.epi_mask[ju]
        cls_arr = np.where(
            e_i & e_j, "ee",
            np.where((~e_i) & (~e_j), "nn", "cross"),
        )
        for band in bands:
            for ph in PHASES:
                r = load_lrg_result(pat, ph, band, fc_method="imcoh_abs")
                if r is None:
                    continue
                from lrg_eegfc.utils.metrics.tree_distance import kc_vectors
                _, M = kc_vectors(np.asarray(r.linkage_matrix))
                for cls in ("ee", "nn", "cross"):
                    msk = cls_arr == cls
                    vals = M[msk]
                    if vals.size == 0:
                        continue
                    key = (band, ph, cls)
                    pooled.setdefault(key, []).append(vals)

    for ax, band in zip(axes, bands):
        positions = []
        data = []
        colors = []
        labels = []
        for pi, ph in enumerate(PHASES):
            for ci, cls in enumerate(("ee", "nn", "cross")):
                key = (band, ph, cls)
                vals_list = pooled.get(key, [])
                if not vals_list:
                    continue
                vals = np.concatenate(vals_list)
                positions.append(pi * 4 + ci)
                data.append(vals)
                colors.append(CLS_COLORS[cls])
                labels.append(f"{ph}-{cls}")

        if not data:
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center")
            continue
        parts = ax.violinplot(
            data, positions=positions, showmedians=True, widths=0.85,
            showextrema=False,
        )
        for body, c in zip(parts["bodies"], colors):
            body.set_facecolor(c)
            body.set_alpha(0.7)
            body.set_edgecolor("black")
            body.set_linewidth(0.4)
        if "cmedians" in parts:
            parts["cmedians"].set_color("black")
            parts["cmedians"].set_linewidth(1.0)

        ax.set_xticks([0.5 + 4 * i for i in range(len(PHASES))])
        ax.set_xticklabels([PHASE_SHORT_LABELS[p] for p in PHASES], fontsize=10)
        ax.set_title(BRAIN_BAND_TEX_DICT[band], fontsize=12)
        ax.tick_params(labelsize=9)

    for ax in (axes[0], axes[3]):
        ax.set_ylabel("MRCA merge height\n(LRG dendrogram)", fontsize=10)

    proxy = [
        plt.Rectangle((0, 0), 1, 1, fc=CLS_COLORS["ee"], ec="black",
                      label="epi-epi"),
        plt.Rectangle((0, 0), 1, 1, fc=CLS_COLORS["nn"], ec="black",
                      label="non-non"),
        plt.Rectangle((0, 0), 1, 1, fc=CLS_COLORS["cross"], ec="black",
                      label="epi-non"),
    ]
    fig.legend(handles=proxy, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.02), fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(FIG_DIR / "fig_04_mrca_violins.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# fig 05 — KC epi-subtree triangle
# ---------------------------------------------------------------------------

def fig_05_kc_epi_triangle() -> None:
    kc = pd.read_csv(CSV_DIR / "M3_kc_epi_subtree.csv")

    # Triangle scalar T_KC_epi = d(rsPre, task) - d(task, rsPost)
    # Project convention: T_KC_epi > 0 = trace at epi sub-network.
    rows = []
    for (pat, band, lam), g in kc.groupby(["patient", "band", "lam"]):
        d_pre_tt = g.loc[(g.phase_a == "rest_pre") & (g.phase_b == "task_test"),
                          "d_kc_epi"]
        d_tt_post = g.loc[(g.phase_a == "task_test") & (g.phase_b == "rest_post"),
                           "d_kc_epi"]
        if d_pre_tt.empty or d_tt_post.empty:
            continue
        rows.append(dict(
            patient=pat, band=band, lam=lam,
            d_pre_tt=float(d_pre_tt.iloc[0]),
            d_tt_post=float(d_tt_post.iloc[0]),
            T_KC_epi=float(d_pre_tt.iloc[0] - d_tt_post.iloc[0]),
        ))
    triangle = pd.DataFrame(rows)
    triangle.to_csv(CSV_DIR / "M3_kc_epi_subtree_triangle.csv", index=False)

    bands = list(BRAIN_BANDS_NAMES)
    fig, axes = plt.subplots(1, len(KC_LAMBDAS), figsize=(13, 4.0), sharey=True)

    for ax, lam in zip(axes, KC_LAMBDAS):
        sub = triangle[triangle.lam == lam]
        # boxplot per band
        data = [sub[sub.band == b]["T_KC_epi"].values for b in bands]
        bp = ax.boxplot(
            data, positions=range(len(bands)), widths=0.6,
            patch_artist=True, showmeans=False, showfliers=False,
            medianprops=dict(color="black", linewidth=1.2),
            boxprops=dict(facecolor="#e6e6e6", edgecolor="0.3", linewidth=0.6),
            whiskerprops=dict(color="0.3", linewidth=0.6),
            capprops=dict(color="0.3", linewidth=0.6),
        )
        # per-patient dots
        rng = np.random.default_rng(7)
        for bi, b in enumerate(bands):
            v = sub[sub.band == b]["T_KC_epi"].values
            jitter = (rng.random(len(v)) - 0.5) * 0.18
            ax.scatter(np.full(len(v), bi) + jitter, v, s=22,
                       color="#1f77b4", alpha=0.85, edgecolor="white",
                       linewidth=0.4, zorder=3)
            # per-band Wilcoxon > 0 (trace direction)
            if len(v) >= 3:
                try:
                    _, p = wilcoxon(v, alternative="greater")
                    if p < 0.05:
                        ax.text(bi, ax.get_ylim()[1] * 0.92 if ax.get_ylim()[1] > 0
                                else ax.get_ylim()[0] * 0.92,
                                "*", ha="center", fontsize=14, color="#d62728")
                except ValueError:
                    pass
        ax.axhline(0, color="0.4", lw=0.7, ls="--")
        ax.set_xticks(range(len(bands)))
        ax.set_xticklabels([BRAIN_BAND_TEX_DICT[b] for b in bands],
                           fontsize=9, rotation=20)
        kind = {0.0: "topology", 0.5: "balanced", 1.0: "heights"}[lam]
        ax.set_title(rf"$\lambda = {lam}$ ({kind})", fontsize=11)
        ax.tick_params(labelsize=9)
    axes[0].set_ylabel(r"$T_{KC}^{\mathrm{epi}}$ = $d(\mathrm{rsPre},\mathrm{task}) - "
                       r"d(\mathrm{task},\mathrm{rsPost})$"
                       "\n  (positive = trace at epi sub-network)", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig_05_kc_epi_triangle.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# fig 06 — dendrograms with epi leaves coloured (β only)
# ---------------------------------------------------------------------------

def fig_06_dendrograms_beta() -> None:
    band = "beta"
    fig, axes = plt.subplots(len(COHORT), len(PHASES),
                              figsize=(11.5, 1.7 * len(COHORT)),
                              sharex=False, sharey=False)
    if axes.ndim == 1:
        axes = axes[None, :]

    for ri, pat in enumerate(COHORT):
        try:
            masks = _build_masks(pat)
        except Exception:
            continue
        if masks.epi_mask.sum() == 0:
            for ax in axes[ri]:
                ax.set_visible(False)
            continue
        for ci, ph in enumerate(PHASES):
            ax = axes[ri, ci]
            r = load_lrg_result(pat, ph, band, fc_method="imcoh_abs")
            if r is None:
                ax.set_visible(False)
                continue
            Z = np.asarray(r.linkage_matrix)
            heights = np.sort(Z[:, 2])
            tmin = max(heights[0] * 0.8, 1e-6)
            tmax = heights[-1] * 1.05
            n_leaves = Z.shape[0] + 1

            # Build link colors based on whether the descendant subtree is
            # entirely epi-leaves (highlight) or not (grey).
            from scipy.cluster.hierarchy import to_tree

            tree = to_tree(Z)

            def gather_leaves(node):
                """Return frozenset of leaf-ids under `node`."""
                if node.is_leaf():
                    return frozenset([int(node.id)])
                left = gather_leaves(node.get_left())
                right = gather_leaves(node.get_right())
                cache[node.id] = left | right
                return cache[node.id]

            cache: dict[int, frozenset[int]] = {}
            gather_leaves(tree)

            link_colors: dict[int, str] = {}
            for k in range(Z.shape[0]):
                cid = n_leaves + k
                leafset = cache.get(cid, frozenset())
                if all(masks.epi_mask[i] for i in leafset):
                    link_colors[cid] = "#d62728"
                elif any(masks.epi_mask[i] for i in leafset):
                    link_colors[cid] = "#cc8a8a"
                else:
                    link_colors[cid] = "#bdbdbd"

            def link_color_func(c):
                return link_colors.get(c, "#bdbdbd")

            with np.errstate(divide="ignore"):
                dendrogram(
                    Z, ax=ax, no_labels=True,
                    color_threshold=None,
                    link_color_func=link_color_func,
                    above_threshold_color="#bdbdbd",
                )
            ax.set_yscale("log")
            ax.set_ylim(tmin, tmax)
            ax.set_xticks([])
            ax.tick_params(labelsize=7, axis="y")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            if ri == 0:
                ax.set_title(PHASE_SHORT_LABELS[ph], fontsize=10)
            if ci == 0:
                ax.set_ylabel(pat, fontsize=9, rotation=0, ha="right",
                              va="center", labelpad=18)

    fig.tight_layout(w_pad=0.4, h_pad=0.4)
    fig.savefig(FIG_DIR / "fig_06_dendrograms_beta.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("[audit_48b] fig_01 amplification scatter…")
    fig_01_amplification_scatter()
    print("[audit_48b] fig_02 enrichment heatmap…")
    fig_02_enrichment_heatmap()
    print("[audit_48b] fig_03 τ-sweep…")
    fig_03_tau_sweep()
    print("[audit_48b] fig_04 MRCA violins…")
    fig_04_mrca_violins()
    print("[audit_48b] fig_05 KC epi-subtree triangle…")
    fig_05_kc_epi_triangle()
    print("[audit_48b] fig_06 dendrograms (β)…")
    fig_06_dendrograms_beta()
    print(f"\n[audit_48b] all figures at {FIG_DIR}")


if __name__ == "__main__":
    main()
