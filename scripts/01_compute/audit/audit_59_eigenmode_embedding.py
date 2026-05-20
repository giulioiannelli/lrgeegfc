#!/usr/bin/env python3
"""Audit 59 — Eigenmode-embedding motion across phases.

Scope at .agents/guides/task-persistence-investigation/2026-05-08_eigenmode_embedding.md.

For each (patient, band) at fc_method='imcoh_abs', the leading k Laplacian
eigenvectors of the cached LRG result span an N x k embedding Y(k); each
contact i has a point Y_i(k) in R^k. Cross-phase motion is the norm of
the embedding shift after orthogonal Procrustes alignment of the eigenbases:

    Y_pre R_pre_to_tt  ~  Y_tt
    Y_pre R_pre_to_post ~  Y_post

Per-contact motion magnitudes:
    dist_rspre_taskt(i)  = || Y_tt - Y_pre R_pre_to_tt ||_2  on row i
    dist_taskt_rspost(i) = || Y_post - Y_tt R_tt_to_post ||_2 on row i
    dist_rspre_rspost(i) = || Y_post - Y_pre R_pre_to_post ||_2 on row i

Library reuse:
    load_lrg_result               cached eigvals/eigvecs
    load_channel_labels           per-patient labels (FC ordering)
    chordal_distance              sanity-only Grassmann reference
    orthogonal_procrustes         scipy.linalg

First-pass cell: Pat_02 beta at k=3. Outputs the per-contact CSV for the
full n=10 cohort at k=3 (cheap) for the universality teaser. Figure shows
2D scatters for Pat_02 beta + per-contact motion bar-plots + cohort teaser.

Outputs
-------
    data/audit/eigenmode_embedding/per_contact_motion.csv
    data/reports/notes_verification_2026-05-08/figures/eigenmode_motion_pat02_beta.pdf
"""
from __future__ import annotations

import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.linalg import orthogonal_procrustes

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.io import load_channel_labels
from lrg_eegfc.workflow.lrg import load_lrg_result


# n=10 cohort locked 2026-04-25.
COHORT = ("Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15")
PHASES = ("rest_pre", "task_test", "rest_post")
PHASE_TEX = {
    "rest_pre":  r"$\mathrm{R}_{\mathrm{pre}}$",
    "task_test": r"$\mathrm{T}_{\mathrm{T}}$",
    "rest_post": r"$\mathrm{R}_{\mathrm{post}}$",
}
PHASE_COLOR = {
    "rest_pre":  "#1f77b4",  # blue
    "task_test": "#d62728",  # red
    "rest_post": "#2ca02c",  # green
}
K_FIRST = 3
COHORT_BAND = "beta"

OUT_CSV = ROOT / "data" / "audit" / "eigenmode_embedding" / "per_contact_motion.csv"
OUT_PDF = (ROOT / "data" / "reports" / "notes_verification_2026-05-08"
           / "figures" / "eigenmode_motion_pat02_beta.pdf")


# ─────────────────────────── core math ───────────────────────────


def topk_embedding(eigenvectors: np.ndarray, k: int) -> np.ndarray:
    """Return the N x k embedding from columns 1..k of an ascending-eigenvalue
    eigenvector matrix. Column 0 is the trivial mode (lambda_1 = 0).
    """
    if k + 1 > eigenvectors.shape[1]:
        raise ValueError(f"k={k} too large for N={eigenvectors.shape[1]}")
    return np.ascontiguousarray(eigenvectors[:, 1: k + 1])


def aligned_embedding(Y_ref: np.ndarray, Y_other: np.ndarray) -> np.ndarray:
    """Orthogonal-Procrustes-align Y_ref to Y_other: returns Y_ref @ R minimising
    ||Y_ref @ R - Y_other||_F. Uses scipy.linalg.orthogonal_procrustes which
    works on tall N x k matrices (returns R of shape k x k).
    """
    R, _ = orthogonal_procrustes(Y_ref, Y_other)
    return Y_ref @ R


def per_contact_distance(Y_a: np.ndarray, Y_b: np.ndarray) -> np.ndarray:
    """Row-wise Euclidean distance ||Y_a[i] - Y_b[i]||_2."""
    return np.linalg.norm(Y_a - Y_b, axis=1)


def gap_min(eigvals: np.ndarray, k_max: int) -> float:
    """Smallest non-zero gap between the leading k_max+1 eigenvalues. Used
    as a degeneracy diagnostic — gap_min < 1e-6 flags non-unique eigenvectors
    in the Procrustes-aligned subspace.
    """
    sub = eigvals[: k_max + 1]
    diffs = np.diff(sub)
    return float(np.min(diffs[diffs > 0])) if np.any(diffs > 0) else 0.0


# ─────────────────────────── per-patient computation ───────────────────────────


def _load_eigvecs_three_phases(pat: str, band: str
                                ) -> tuple[dict[str, np.ndarray] | None,
                                           dict[str, np.ndarray] | None]:
    eigvecs: dict[str, np.ndarray] = {}
    eigvals: dict[str, np.ndarray] = {}
    for phase in PHASES:
        r = load_lrg_result(pat, phase, band, "imcoh_abs")
        if r is None or r.eigenvectors is None or r.eigenvalues is None:
            return None, None
        eigvecs[phase] = np.asarray(r.eigenvectors)
        eigvals[phase] = np.asarray(r.eigenvalues)
    Ns = {V.shape[0] for V in eigvecs.values()}
    if len(Ns) != 1:
        return None, None
    return eigvecs, eigvals


def compute_motions(pat: str, band: str, k: int) -> list[dict] | None:
    """Returns one row per contact with per-pair motion magnitudes."""
    eigvecs, eigvals = _load_eigvecs_three_phases(pat, band)
    if eigvecs is None:
        return None
    try:
        labels = load_channel_labels(pat)
    except Exception:
        labels = [f"i{n}" for n in range(eigvecs["rest_pre"].shape[0])]
    if len(labels) != eigvecs["rest_pre"].shape[0]:
        # Defensive: fall back to integer labels, log row's count.
        labels = [f"i{n}" for n in range(eigvecs["rest_pre"].shape[0])]

    Y_pre = topk_embedding(eigvecs["rest_pre"], k)
    Y_tt = topk_embedding(eigvecs["task_test"], k)
    Y_post = topk_embedding(eigvecs["rest_post"], k)

    # Anchor on rest_pre. Align the other two phases to rest_pre's embedding.
    Y_tt_a = aligned_embedding(Y_tt, Y_pre)
    Y_post_a = aligned_embedding(Y_post, Y_pre)

    # The third pair (tt -> post) is computed via the chained alignment: we
    # measure motion in the rest_pre reference frame, so consistency demands
    # |Y_post_a - Y_tt_a| (both in the same frame), NOT a fresh tt -> post
    # Procrustes (which would rotate them away from the shared frame).
    d_pre_tt = per_contact_distance(Y_pre, Y_tt_a)
    d_tt_post = per_contact_distance(Y_tt_a, Y_post_a)
    d_pre_post = per_contact_distance(Y_pre, Y_post_a)

    gap_pre = gap_min(eigvals["rest_pre"], k)
    gap_tt = gap_min(eigvals["task_test"], k)
    gap_post = gap_min(eigvals["rest_post"], k)

    rows: list[dict] = []
    for i, lab in enumerate(labels):
        rows.append(dict(
            patient=pat, band=band, k=k, contact=i, channel_label=lab,
            dist_rspre_taskt=float(d_pre_tt[i]),
            dist_taskt_rspost=float(d_tt_post[i]),
            dist_rspre_rspost=float(d_pre_post[i]),
            gap_min_pre=gap_pre,
            gap_min_tt=gap_tt,
            gap_min_post=gap_post,
        ))
    return rows


# ─────────────────────────── figure ───────────────────────────


def _scatter_pair(ax, Y_pre, Y_tt, Y_post, dim_a: int, dim_b: int,
                   d_pre_tt: np.ndarray, d_tt_post: np.ndarray):
    """Single eigenmode-pair scatter with 3 phases overlaid."""
    # Marker size proportional to (max(d_pre_tt, d_tt_post)) per contact;
    # keeps the eye on contacts with large motion.
    motion = np.maximum(d_pre_tt, d_tt_post)
    s = 8 + 60 * (motion / max(motion.max(), 1e-9))
    ax.scatter(Y_pre[:, dim_a], Y_pre[:, dim_b],
               c=PHASE_COLOR["rest_pre"], s=s, alpha=0.6,
               edgecolor="none", label=PHASE_TEX["rest_pre"])
    ax.scatter(Y_tt[:, dim_a], Y_tt[:, dim_b],
               c=PHASE_COLOR["task_test"], s=s, alpha=0.6,
               edgecolor="none", label=PHASE_TEX["task_test"])
    ax.scatter(Y_post[:, dim_a], Y_post[:, dim_b],
               c=PHASE_COLOR["rest_post"], s=s, alpha=0.6,
               edgecolor="none", label=PHASE_TEX["rest_post"])
    ax.set_xlabel(f"$v_{{{dim_a + 2}}}$", fontsize=9)
    ax.set_ylabel(f"$v_{{{dim_b + 2}}}$", fontsize=9)
    ax.tick_params(labelsize=7)
    ax.grid(True, alpha=0.25, linewidth=0.4)


def _bar_top_movers(ax, labels: list[str], d_pre_tt: np.ndarray,
                     d_tt_post: np.ndarray, n_top: int = 15):
    """Bar plot of top-n_top contacts by max(motion). Two bars per contact."""
    motion = np.maximum(d_pre_tt, d_tt_post)
    order = np.argsort(motion)[::-1][:n_top]
    pos = np.arange(n_top)
    width = 0.4
    ax.bar(pos - width / 2, d_pre_tt[order], width=width,
           color=PHASE_COLOR["task_test"], alpha=0.85,
           edgecolor="none",
           label=r"$\Delta_i^{R_{\mathrm{pre}}\to T_T}$")
    ax.bar(pos + width / 2, d_tt_post[order], width=width,
           color=PHASE_COLOR["rest_post"], alpha=0.85,
           edgecolor="none",
           label=r"$\Delta_i^{T_T\to R_{\mathrm{post}}}$")
    ax.set_xticks(pos)
    ax.set_xticklabels([labels[i] for i in order],
                       rotation=70, fontsize=6.5, ha="right")
    ax.set_ylabel("embedding shift", fontsize=9)
    ax.tick_params(axis="y", labelsize=7)
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.4)
    ax.legend(loc="upper right", fontsize=7, frameon=False)
    ax.set_xlim(-0.7, n_top - 0.3)


def _cohort_teaser(ax, df_all: pd.DataFrame, band: str):
    """One bar per patient: cohort-mean motion for pre->tt and tt->post."""
    sub = df_all[df_all.band == band]
    pats = sorted(sub.patient.unique())
    means_a = []
    means_b = []
    for p in pats:
        psub = sub[sub.patient == p]
        means_a.append(psub.dist_rspre_taskt.mean())
        means_b.append(psub.dist_taskt_rspost.mean())
    pos = np.arange(len(pats))
    width = 0.4
    ax.bar(pos - width / 2, means_a, width=width,
           color=PHASE_COLOR["task_test"], alpha=0.85,
           edgecolor="none",
           label=r"$\overline{\Delta}^{R_{\mathrm{pre}}\to T_T}$")
    ax.bar(pos + width / 2, means_b, width=width,
           color=PHASE_COLOR["rest_post"], alpha=0.85,
           edgecolor="none",
           label=r"$\overline{\Delta}^{T_T\to R_{\mathrm{post}}}$")
    ax.set_xticks(pos)
    short = [p.replace("Pat_", "P") for p in pats]
    ax.set_xticklabels(short, fontsize=6.5)
    ax.set_ylabel("contact-mean shift", fontsize=9)
    ax.tick_params(axis="y", labelsize=7)
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.4)
    ax.legend(loc="upper right", fontsize=7, frameon=False)


def make_figure(df_all: pd.DataFrame, k: int) -> None:
    """Render the Pat_02 beta first-pass figure plus cohort teaser."""
    pat = "Pat_02"
    band = COHORT_BAND
    eigvecs, _ = _load_eigvecs_three_phases(pat, band)
    if eigvecs is None:
        raise RuntimeError(f"{pat} {band}: no cached eigenvectors")
    Y_pre = topk_embedding(eigvecs["rest_pre"], k)
    Y_tt = aligned_embedding(topk_embedding(eigvecs["task_test"], k), Y_pre)
    Y_post = aligned_embedding(topk_embedding(eigvecs["rest_post"], k), Y_pre)
    d_pre_tt = per_contact_distance(Y_pre, Y_tt)
    d_tt_post = per_contact_distance(Y_tt, Y_post)

    sub = df_all[(df_all.patient == pat) & (df_all.band == band)]
    labels = sub.sort_values("contact")["channel_label"].tolist()
    if len(labels) != Y_pre.shape[0]:
        labels = [f"i{n}" for n in range(Y_pre.shape[0])]

    fig = plt.figure(figsize=(11.0, 7.4))
    gs = fig.add_gridspec(
        2, 3,
        height_ratios=[1.0, 1.0],
        hspace=0.45, wspace=0.30,
        left=0.07, right=0.985, top=0.92, bottom=0.10,
    )
    ax_12 = fig.add_subplot(gs[0, 0])
    ax_13 = fig.add_subplot(gs[0, 1])
    ax_23 = fig.add_subplot(gs[0, 2])
    ax_bar = fig.add_subplot(gs[1, :2])
    ax_cohort = fig.add_subplot(gs[1, 2])

    _scatter_pair(ax_12, Y_pre, Y_tt, Y_post, 0, 1, d_pre_tt, d_tt_post)
    _scatter_pair(ax_13, Y_pre, Y_tt, Y_post, 0, 2, d_pre_tt, d_tt_post)
    _scatter_pair(ax_23, Y_pre, Y_tt, Y_post, 1, 2, d_pre_tt, d_tt_post)

    # Single shared phase legend on the figure (top centre); per-rule, no
    # axis-level legend on the scatter panels.
    handles, lbls = ax_12.get_legend_handles_labels()
    fig.legend(handles, lbls, loc="upper center",
               bbox_to_anchor=(0.5, 0.985),
               ncol=len(handles), frameon=False, fontsize=9)

    _bar_top_movers(ax_bar, labels, d_pre_tt, d_tt_post, n_top=15)
    _cohort_teaser(ax_cohort, df_all, band)

    # Header text via fig.text (no fig.suptitle per project rule).
    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    fig.text(
        0.07, 0.955,
        rf"Pat_02 {band_tex} eigenmode-embedding motion (k=3)",
        fontsize=10, fontweight="bold",
    )
    fig.text(
        0.985, 0.955,
        f"cohort teaser (n={len(df_all.patient.unique())}, {band_tex})",
        fontsize=9, ha="right", color="0.35",
    )

    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PDF, format="pdf")
    plt.close(fig)


# ─────────────────────────── orchestrator ───────────────────────────


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("--k", type=int, default=K_FIRST,
                   help=f"Embedding order (default: {K_FIRST})")
    p.add_argument("--band", default=COHORT_BAND,
                   help=f"Band for cohort sweep (default: {COHORT_BAND})")
    args = p.parse_args(argv)

    rows: list[dict] = []
    for pat in COHORT:
        out = compute_motions(pat, args.band, args.k)
        if out is None:
            if args.verbose:
                print(f"[audit_59] skip {pat} {args.band}: missing eigenvectors")
            continue
        rows.extend(out)
        if args.verbose:
            n = len(out)
            print(f"[audit_59] {pat} {args.band}: k={args.k}, {n} contacts")

    if not rows:
        print("[audit_59] no rows computed — eigenvector cache missing?",
              file=sys.stderr)
        return 1

    df = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"[audit_59] CSV: {OUT_CSV}  ({len(df)} rows, "
          f"{df.patient.nunique()} patients)")

    make_figure(df, args.k)
    print(f"[audit_59] PDF: {OUT_PDF}")

    # Cohort summary line.
    sub = df[df.band == args.band].groupby("patient").agg(
        mean_pre_tt=("dist_rspre_taskt", "mean"),
        mean_tt_post=("dist_taskt_rspost", "mean"),
        max_pre_tt=("dist_rspre_taskt", "max"),
    ).round(4)
    print("\n[audit_59] cohort summary (k=3, band="
          f"{args.band}):\n{sub.to_string()}")

    # Pat_02 beta: top-3 movers by max(motion).
    p02 = df[(df.patient == "Pat_02") & (df.band == args.band)].copy()
    p02["motion_max"] = p02[["dist_rspre_taskt",
                             "dist_taskt_rspost"]].max(axis=1)
    top = p02.nlargest(5, "motion_max")[
        ["channel_label", "dist_rspre_taskt",
         "dist_taskt_rspost", "dist_rspre_rspost"]]
    print(f"\n[audit_59] Pat_02 {args.band} top-5 movers (k={args.k}):")
    print(top.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
