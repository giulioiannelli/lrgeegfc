#!/usr/bin/env python3
r"""Talk slide — the "network-analysis pipeline" as FIVE separate clean panels.

One script, five standalone high-res transparent PNGs (Canva import; + QA PNGs). Left-to-right
in Canva and draws the arrows between them, so each panel is a self-contained icon with no
cross-panel decoration. Every panel uses ONE (patient, band) — default Pat_05 / beta — so the
sequence tells a single coherent story: raw signal -> connectivity -> diffusion -> hierarchy ->
cross-phase similarity.

    1  pipeline_seq_1_timeseries.pdf        a handful of raw sEEG contacts (rest_pre window)
    2  pipeline_seq_2_imcoh_matrix.pdf      the |ImCoh| functional-connectivity matrix A
    3  pipeline_seq_3_distance.pdf          the communication distance D(tau)=(1-d_ij)/K, K=exp(-tau L)
    4  pipeline_seq_4_dendrogram.pdf        the cophenetic (LRG) dendrogram
    5  pipeline_seq_5_rhocoph_tanglegram.pdf   rho^coph explainer: rest_pre vs rest_post trees

Panel 5 is an EXPLAINER (no cohort statistic): rho^coph = Spearman of the two phases'
cophenetic-distance vectors. To stay legible we crop to one mid-size rest_post clade (~18
contacts), untangle both induced sub-trees by topology-preserving branch rotation, and
HIGHLIGHT a sub-clade that keeps its grouping across the two phases (parallel amber ribbons =
retained; crossing grey ribbons = reorganised). The annotated rho is the LOCAL subset value
(labelled "illustrative subset"); the honest full-matrix rho^coph is printed underneath.

Nothing is recomputed that a cache already holds; if a required cache is missing the script
reports which one and stops (no fabricated data). House rules: use_lrg_style(), PDF-only +
transparent, no suptitle, no rasterisation, plt.close per panel.

Usage:
    /home/giulio/Documents/miniconda3/envs/lapbrain/bin/python \
        scripts/07_figures/gen_pipeline_sequence.py --patient Pat_05 --band beta
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.cluster.hierarchy import dendrogram, optimal_leaf_ordering
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.config.paths import SEEG_DATAPATH, IMCOH_LRG_CACHE, FIGURES_ROOT
from lrg_eegfc.config.const import FS_OVERRIDES, DEFAULT_SAMPLE_RATE
from lrg_eegfc.utils.io import load_timeseries
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.visuals.fc_templates import plot_fc_adjacency

# ---- tanglegram machinery: reuse the cophenetic-tanglegram helpers (DRY) ----
# These modules live under scripts/ and self-configure sys.path at import; we
# only need their directories importable first.
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "figures_embedded"))
from fig_cophenetic_tanglegram import (            # type: ignore  # noqa: E402
    induced_linkage, untangle, ribbon, subset_rho, crossings, clades, draw_subtree,
)
from fig_branch_origin_tree import coph_square      # type: ignore  # noqa: E402
from audit_63_split_baseline_surrogate import load_phase_fc   # type: ignore  # noqa: E402

OUTDIR = FIGURES_ROOT / "talk"
QA = Path(
    "/tmp/claude-1000/-home-giulio-Documents-research-neural-networks-lrgeegfc/"
    "20c35621-890b-48ca-8413-660e8b128ed1/scratchpad"
)

# palette
DARK = "#1f2a37"          # raw sEEG traces / dendrogram ink
BRACKET_INK = "#2b2f36"   # tree brackets + phase headers
HL = "#e08a1e"            # retained-clade highlight (amber)
GREY_RIB = "#b9bdc4"      # reorganised (faint) connectors
BEAD_GREY = "#cfd3d9"     # non-highlighted leaf beads


def _save(fig, name: str) -> Path:
    """Write a panel to a high-res TRANSPARENT PNG for Canva import (user opt-in 2026-07-11:
    Canva imports PNG far more reliably than PDF) + a white-matte QA PNG, then close it."""
    OUTDIR.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    png = OUTDIR / f"{name}.png"
    fig.savefig(png, bbox_inches="tight", transparent=True, pad_inches=0.03, dpi=300)
    fig.savefig(QA / f"{name}.png", bbox_inches="tight", facecolor="white", dpi=150)
    plt.close(fig)
    print(f"      wrote {png}", flush=True)
    return png


# ============================ Panel 1 — raw sEEG ============================
def panel1_timeseries(patient: str, *, n_contacts: int, win_s: float, start_s: float) -> Path:
    """A few stacked raw sEEG contacts over a short rest_pre window (common gain, seconds axis)."""
    fs = FS_OVERRIDES.get(patient, DEFAULT_SAMPLE_RATE)
    ts = np.asarray(load_timeseries(patient, "rest_pre", SEEG_DATAPATH), dtype=float)
    if ts.shape[0] < ts.shape[1]:                    # -> (n_samples, n_channels)
        ts = ts.T
    n_samp, n_chan = ts.shape
    n = int(round(win_s * fs))
    s0 = int(round(min(start_s, max(0.0, n_samp / fs - win_s)) * fs))
    seg = ts[s0:s0 + n]
    seg = seg - seg.mean(0, keepdims=True)
    sel = np.linspace(0, n_chan - 1, n_contacts).round().astype(int)
    # common gain (median across-channel std * 6) keeps relative amplitudes honest
    spacing = 6.0 * float(np.median(seg[:, sel].std(0))) or 1.0
    t = np.arange(n) / fs

    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    for k, ci in enumerate(sel):
        y0 = (n_contacts - 1 - k) * spacing
        ax.plot(t, seg[:, ci] + y0, lw=0.6, color=DARK, solid_capstyle="round")
    ax.set_yticks([])
    ax.set_xlim(float(t[0]), float(t[-1]))
    ax.set_xlabel("time (s)")
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    print(f"      window {s0 / fs:.1f}-{s0 / fs + win_s:.1f}s @ {int(fs)} Hz  "
          f"contacts={list(map(int, sel))}", flush=True)
    return _save(fig, "pipeline_seq_1_timeseries")


# ======================= Panel 2 — |ImCoh| FC matrix =======================
def _load_fc_or_die(patient: str, phase: str, band: str) -> np.ndarray:
    try:
        A = load_fc_matrix(patient, phase, band, "imcoh_abs")
    except FileNotFoundError as e:
        sys.exit(f"!! missing FC cache for {patient}/{phase}/{band} (imcoh_abs): {e}")
    if A is None:
        sys.exit(f"!! missing FC cache for {patient}/{phase}/{band} (imcoh_abs)")
    return np.asarray(A, dtype=float)


def panel2_fc(patient: str, band: str, A: np.ndarray) -> Path:
    """The |ImCoh| functional-connectivity matrix, log colour, generic FC_{ij} colorbar."""
    fig, ax = plt.subplots(figsize=(3.4, 3.0))
    plot_fc_adjacency(
        A, ax=ax, fc_method="imcoh_abs", band=None,
        colorbar_label=r"$\mathrm{FC}_{ij}$", tick_labels="generic", log_scale=True,
    )
    return _save(fig, "pipeline_seq_2_imcoh_matrix")


# ===================== Panel 3 — communication distance D=1/K ======================
def panel3_distance(patient: str, band: str, A: np.ndarray) -> tuple[Path, float]:
    """LRG communication distance D_ij(tau) = (1-delta_ij)/K_ij, K=exp(-tau L), tau=1/lambda_max.

    This is the matrix that actually FEEDS the hierarchical clustering (UPGMA runs on
    D, not on K). It is the photographic negative of the FC matrix — strong coupling ->
    SMALL distance (dark), weak coupling -> LARGE distance (bright) — so it reads as a
    distinct object from panel 2 rather than a near-copy of it. Off-diagonal, log colour,
    robust vmax (98th pct) so a few huge weak-pair distances don't wash out the structure."""
    L = np.diag(A.sum(1)) - A
    lam = np.linalg.eigvalsh(L)
    tau = 1.0 / float(lam.max())
    K = expm(-tau * L)
    if not np.allclose(K, K.T, atol=1e-9):
        sys.exit("!! propagator K is not symmetric — refusing to plot")
    with np.errstate(divide="ignore"):
        Dist = 1.0 / K                                   # communication distance
    np.fill_diagonal(Dist, 0.0)                          # D_ii = 0 (Villegas 1-delta_ij)
    Dist = np.maximum(Dist, Dist.T)
    off = Dist[np.triu_indices_from(Dist, 1)]
    vmin = float(off[off > 0].min())
    vmax = float(np.percentile(off, 98))                 # robust: clip the weak-pair tail
    fig, ax = plt.subplots(figsize=(3.4, 3.0))
    plot_fc_adjacency(
        Dist, ax=ax, fc_method="imcoh_abs", band=None,
        colorbar_label=r"$D_{ij}(\tau)$", tick_labels="generic", log_scale=True,
        vmin=vmin, vmax=vmax,
    )
    return _save(fig, "pipeline_seq_3_distance"), tau


# ===================== Panel 4 — cophenetic dendrogram =====================
def panel4_dendrogram(patient: str, band: str) -> Path:
    """The LRG cophenetic dendrogram: clean monochrome tree, log merge-height axis, no leaves."""
    r = load_lrg_result(patient, "rest_pre", band, "imcoh_abs", IMCOH_LRG_CACHE)
    if r is None:
        sys.exit(f"!! missing LRG cache for {patient}/rest_pre/{band} (imcoh_abs)")
    Z = r.linkage_matrix
    try:
        Z = optimal_leaf_ordering(Z, r.ultrametric_matrix)
    except Exception:                                # noqa: BLE001
        pass
    h = r.linkage_matrix[:, 2]

    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    dendrogram(Z, ax=ax, no_labels=True, color_threshold=0.0,
               above_threshold_color=BRACKET_INK)
    ax.set_yscale("log")
    ax.set_ylim(float(h.min()) * 0.8, float(h.max()) * 1.05)   # house convention
    ax.set_ylabel("ultrametric distance")
    ax.set_xticks([])
    for sp in ("top", "right", "bottom"):
        ax.spines[sp].set_visible(False)
    return _save(fig, "pipeline_seq_4_dendrogram")


# ==================== Panel 5 — rho^coph tanglegram ========================
def _subclades_local(Z: np.ndarray, m: int, lo: int, hi: int) -> list[list[int]]:
    """Local-index clades (leaf lists) of an induced tree of ``m`` leaves, size in [lo, hi]."""
    members = {i: [i] for i in range(m)}
    out: list[list[int]] = []
    for k in range(m - 1):
        mm = members[int(Z[k, 0])] + members[int(Z[k, 1])]
        members[m + k] = mm
        if lo <= len(mm) <= hi:
            out.append(sorted(mm))
    return out


def _layout_S(Cpre: np.ndarray, Cpost: np.ndarray, S: list[int]):
    """Untangled two-tree layout for a subset S: fix rest_post order, branch-rotate rest_pre to it."""
    m = len(S)
    Zp = induced_linkage(Cpre, S)
    Zr = induced_linkage(Cpost, S)
    Cs = Cpost[np.ix_(S, S)].astype(float)
    Zr = optimal_leaf_ordering(
        Zr, squareform(0.5 * (Cs + Cs.T) - np.diag(np.diag(Cs)), checks=False)
    )
    pos_post = {l: i for i, l in enumerate(dendrogram(Zr, no_plot=True)["leaves"])}
    pos_pre = {l: i for i, l in enumerate(untangle(Zp, m, pos_post))}
    return m, Zp, Zr, pos_pre, pos_post, subset_rho(Cpre, Cpost, S)


def _select_retained(Cpre, Cpost, Zpost, N, *, size_lo=16, size_hi=22):
    """Pick (subset S, retained clade G) for the explainer.

    A retained clade = a rest_post sub-clade that is ALSO a contiguous block in the untangled
    rest_pre order (so its grouping survived the phase change). Composite score rewards a
    crossing-free clade (gc=0), a larger clade, and a subset rho near 0.6 so the surrounding
    ribbons still show some reorganisation for contrast. Deterministic single pass.
    """
    best = None
    for S in clades(Zpost, N, size_lo, size_hi):
        m, Zp, Zr, pos_pre, pos_post, rhoS = _layout_S(Cpre, Cpost, S)
        for G in _subclades_local(Zr, m, 4, 8):
            PL = sorted(pos_pre[l] for l in G)
            spread = PL[-1] - PL[0]
            blockL = (len(G) - 1) / spread if spread > 0 else 1.0
            if blockL < 0.999:                       # must be contiguous in rest_pre too
                continue
            gc = crossings([pos_pre[l] for l in G], [pos_post[l] for l in G], len(G))
            in_band = 0.40 <= rhoS <= 0.85
            score = len(G) - 3.0 * gc - 2.0 * abs(rhoS - 0.6) - (0.0 if in_band else 5.0)
            key = (score, len(G), -gc)
            if best is None or key > best[0]:
                best = (key, dict(S=S, G=set(G), m=m, Zp=Zp, Zr=Zr,
                                  pos_pre=pos_pre, pos_post=pos_post,
                                  rhoS=rhoS, blockL=blockL, gc=gc))
    return None if best is None else best[1]


def panel5_tanglegram(patient: str, band: str) -> tuple[Path, dict]:
    """rho^coph explainer tanglegram: rest_pre (left) vs rest_post (right), one retained clade lit."""
    Wpre = load_phase_fc(patient, "rest_pre", band)
    Wpost = load_phase_fc(patient, "rest_post", band)
    if Wpre.shape != Wpost.shape:
        sys.exit(f"!! phase N mismatch {Wpre.shape} vs {Wpost.shape} — cannot align pairs")
    Cpre, _ = coph_square(Wpre)
    Cpost, Zpost = coph_square(Wpost)
    N = Cpre.shape[0]
    iu = np.triu_indices(N, 1)
    rho_full = float(spearmanr(Cpre[iu], Cpost[iu]).statistic)

    pick = _select_retained(Cpre, Cpost, Zpost, N)
    if pick is None:
        sys.exit("!! no clean retained clade found for the tanglegram (relax size window)")
    m, Zp, Zr = pick["m"], pick["Zp"], pick["Zr"]
    pos_pre, pos_post, G, rhoS = pick["pos_pre"], pick["pos_post"], pick["G"], pick["rhoS"]

    RAIL_L, RAIL_R, EXT = 0.0, 1.0, 0.34
    allh = np.concatenate([Zp[:, 2], Zr[:, 2]])
    pos = allh[allh > 0]
    hlo, hhi = float(np.log(pos.min())), float(np.log(pos.max()))

    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    # left tree opens left (leaves face centre), right tree opens right
    draw_subtree(ax, Zp, m, pos_pre, RAIL_L, -1, EXT, hlo, hhi)
    draw_subtree(ax, Zr, m, pos_post, RAIL_R, +1, EXT, hlo, hhi)
    # reorganised connectors first (faint), retained on top (bold amber)
    for l in range(m):
        if l in G:
            continue
        ribbon(ax, RAIL_L, pos_pre[l], RAIL_R, pos_post[l], GREY_RIB, alpha=0.55, lw=1.3, z=2)
    for l in G:
        ribbon(ax, RAIL_L, pos_pre[l], RAIL_R, pos_post[l], HL, alpha=0.95, lw=2.4, z=4)
    for l in range(m):
        in_g = l in G
        ax.scatter([RAIL_L, RAIL_R], [pos_pre[l], pos_post[l]],
                   s=46 if in_g else 30, color=HL if in_g else BEAD_GREY,
                   ec="#20242a" if in_g else "#9aa0a8", lw=0.8, zorder=6 if in_g else 5)

    gy = float(np.mean([pos_post[l] for l in G]))
    ax.text(RAIL_R + EXT + 0.06, gy, "retained\nclade", ha="left", va="center",
            fontsize=9, color=HL, fontweight="bold")
    ymax = m - 1
    ax.text(RAIL_L, ymax + 1.1, "rest-pre", ha="center", va="bottom",
            fontsize=10.5, fontweight="bold", color=BRACKET_INK)
    ax.text(RAIL_R, ymax + 1.1, "rest-post", ha="center", va="bottom",
            fontsize=10.5, fontweight="bold", color=BRACKET_INK)
    ax.text(0.5, ymax + 1.9,
            r"$\rho^{\mathrm{coph}}\ \mathrm{(illustrative\ subset)}=%.2f$" % rhoS,
            ha="center", va="bottom", fontsize=9.5, color=BRACKET_INK)
    ax.text(0.5, -1.5, r"full-matrix $\rho^{\mathrm{coph}}=%.2f$" % rho_full,
            ha="center", va="top", fontsize=8, color="#6a7078")
    ax.set_xlim(RAIL_L - EXT - 0.1, RAIL_R + EXT + 0.5)
    ax.set_ylim(-2.0, ymax + 3.0)
    ax.axis("off")

    info = dict(rho_full=rho_full, rho_subset=rhoS, subset_size=m,
                clade_size=len(G), subset=sorted(pick["S"]), clade_local=sorted(G))
    return _save(fig, "pipeline_seq_5_rhocoph_tanglegram"), info


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Five-panel network-analysis pipeline sequence for a talk slide."
    )
    ap.add_argument("--patient", default="Pat_05")
    ap.add_argument("--band", default="beta")
    ap.add_argument("--n-contacts", type=int, default=6, help="panel 1: stacked sEEG contacts")
    ap.add_argument("--win-s", type=float, default=2.0, help="panel 1: window length (s)")
    ap.add_argument("--start-s", type=float, default=60.0, help="panel 1: window start (s)")
    a = ap.parse_args()

    t0 = time.perf_counter()
    use_lrg_style()
    print(f"[pipeline-seq] patient={a.patient} band={a.band}  ->  {OUTDIR}", flush=True)

    print("[1/5] raw sEEG timeseries", flush=True)
    panel1_timeseries(a.patient, n_contacts=a.n_contacts, win_s=a.win_s, start_s=a.start_s)

    print("[2/5] |ImCoh| FC matrix", flush=True)
    A = _load_fc_or_die(a.patient, "rest_pre", a.band)
    panel2_fc(a.patient, a.band, A)

    print("[3/5] communication distance D(tau)=1/K", flush=True)
    _, tau = panel3_distance(a.patient, a.band, A)
    print(f"      tau = 1/lambda_max = {tau:.4g}", flush=True)

    print("[4/5] cophenetic dendrogram", flush=True)
    panel4_dendrogram(a.patient, a.band)

    print("[5/5] rho^coph tanglegram (explainer)", flush=True)
    _, info = panel5_tanglegram(a.patient, a.band)
    print(f"      full rho^coph={info['rho_full']:.4f}   subset rho^coph={info['rho_subset']:.4f}"
          f"   |subset|={info['subset_size']}   retained clade |G|={info['clade_size']}", flush=True)

    print(f"[done] {time.perf_counter() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
