#!/usr/bin/env python3
r"""fig:arc_c — the encoding component anchors in orbitofrontal cortex (Results §2, R2.3).

CORE MESSAGE: what the network retains from being SHOWN the premises (the encoding
component e = D_task_learn - D_rest_pre, before any inference) concentrates, in beta, in
ORBITOFRONTAL cortex and nowhere else, while prefrontal cortex is significantly depleted.
This is the same orbitofrontal hotspot the test-phase trace occupies (R1.2), so OFC is
anchored by BOTH task phases: learning already writes the trace there.

Per-system enrichment of the encoding cophenetic trace, beta, against the R=1000
strength-matched surrogate, BH-corrected across the 9 a-priori systems within each
direction (enrichment tail q_up, depletion tail q_lo):
  OFC  enriched  q_up = 0.040 (include) / 0.010 (exclude)  -- the only enriched carrier
  PFC  depleted  q_lo = 0.010 (both modes)
OFC sits in nodes of BELOW-average coupling (strength_dev < 0) -- a concentration, not a
hub artifact. OFC is sampled in 5/10 patients (coverage limit, not a null).

Reads : data/audit/inference_localization_rhosym/encoding_localization_rhosym_{include,exclude}.csv
Writes: data/reports/results_section2/fig_arc_c_encoding_ofc_anchor.pdf
"""
from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from lrg_eegfc.utils.metrics.hypothesis import bh_fdr
from lrg_eegfc.utils.scripting import setup_script_env
from lrg_eegfc.visuals.styles import use_lrg_style

ROOT = setup_script_env()
use_lrg_style()

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "figures_embedded"))
try:
    from fig_trace_b_ofc_localization import pooled_coords_systems  # type: ignore
    from nilearn.plotting import plot_glass_brain
    _BRAIN_OK = True
except Exception as exc:  # pragma: no cover
    print(f"[warn] glass-brain machinery unavailable ({exc}); shipping panel a only")
    _BRAIN_OK = False

BASE = ROOT / "data/audit/inference_localization_rhosym"
OUT = ROOT / "data/reports/results_section2/fig_arc_c_encoding_ofc_anchor.pdf"

DROP = {"other", "non_anatomical"}
SYS_ORDER = ["OFC", "cingulate", "MTL", "insula", "lateral_temporal",
             "PFC", "parietal", "sensorimotor", "occipital"]
WINNER = "OFC"
C_ENR, C_DEP, C_NS = "#2a9d5c", "#d1791f", "#9a9a9a"
CAP = 3.0
Q_SIG = 0.05


def load_sys(inc: str) -> dict:
    """system -> dict(enr, p_dir, q_dir, signed) for beta encoding, direction-aware
    one-sided BH (enrichment tail vs depletion tail), matching the R2.3 prose."""
    df = pd.read_csv(BASE / f"encoding_localization_rhosym_{inc}.csv")
    d = df[(df.band == "beta") & (df.granularity == "system") &
           (~df.unit.isin(DROP))].copy()
    r = int(d.n_surr.iloc[0])
    enr = (d.M_obs > d.surr_median).to_numpy()
    up = np.clip(d.matched_strength_p.astype(float), 1.0 / (r + 1), 1.0)
    lo = np.clip(d.matched_strength_p_lower.astype(float), 1.0 / (r + 1), 1.0)
    q_up = np.asarray(bh_fdr(up.tolist()))
    q_lo = np.asarray(bh_fdr(lo.tolist()))
    p_dir = np.where(enr, up, lo)
    q_dir = np.where(enr, q_up, q_lo)
    signed = np.where(enr, 1.0, -1.0) * np.minimum(-np.log10(p_dir), CAP)
    out = {}
    for u, e, pv, qq, sv, sd in zip(d.unit, enr, p_dir, q_dir, signed,
                                    d.strength_dev_median):
        out[str(u)] = dict(enr=bool(e), p=float(pv), q=float(qq),
                           signed=float(sv), strdev=float(sd))
    return out


def draw_lollipop(ax, inc, exc):
    y0 = np.arange(len(SYS_ORDER))[::-1]
    ypos = {s: y for s, y in zip(SYS_ORDER, y0)}
    ax.axvline(0, color="0.55", lw=0.9, zorder=1)
    for s in SYS_ORDER:
        rec = inc.get(s)
        if rec is None:
            continue
        y, x = ypos[s], rec["signed"]
        col = C_ENR if rec["enr"] else C_DEP
        sig = rec["q"] < Q_SIG
        ax.plot([0, x], [y, y], color=col, lw=1.9, alpha=0.9, zorder=2)
        # epi-excluded replicate: pale shadow marker (uniform, all systems)
        xe = exc.get(s, {}).get("signed")
        if xe is not None:
            ax.scatter([xe], [y + 0.16], s=16, marker="o", facecolors="none",
                       edgecolors=col, linewidths=0.8, alpha=0.55, zorder=3)
        if s == WINNER and sig:
            ax.scatter([x], [y], s=180, marker="o", facecolors=col,
                       edgecolors="black", linewidths=1.6, zorder=5)
        elif sig:
            ax.scatter([x], [y], s=105, marker="o", facecolors=col,
                       edgecolors="black", linewidths=1.1, zorder=4)
        else:
            ax.scatter([x], [y], s=70, marker="o", facecolors="white",
                       edgecolors=col, linewidths=1.4, zorder=4)
        if sig:
            ax.annotate(r"$\ast$", (x, y), textcoords="offset points",
                        xytext=(8 if x >= 0 else -8, 3), ha="center",
                        fontsize=13, color=col, zorder=6)
    # headline q's on the two carrier systems (enrichment / depletion)
    for s, side, lbl in ((WINNER, +1, "enriched"), ("PFC", -1, "depleted")):
        rec = inc.get(s)
        if rec and rec["q"] < Q_SIG:
            qi, qe = rec["q"], exc.get(s, {}).get("q", float("nan"))
            ax.annotate(f"$q={qi:.3f}$ / ${qe:.3f}$", (rec["signed"], ypos[s]),
                        textcoords="offset points",
                        xytext=(14 * side, 15), ha="center", fontsize=9.5,
                        color=(C_ENR if side > 0 else C_DEP))

    ax.set_yticks(y0)
    ax.set_yticklabels(SYS_ORDER)
    for lab in ax.get_yticklabels():
        if lab.get_text() in (WINNER, "PFC"):
            lab.set_fontweight("bold")
    ax.set_ylim(y0.min() - 0.7, y0.max() + 0.7)
    ax.set_xlim(-CAP - 0.6, CAP + 0.6)
    ax.set_xlabel(r"$\leftarrow$ depleted      signed enrichment  "
                  r"$\mathrm{sign}(M_{\mathrm{obs}}-\tilde{M}_0)\cdot"
                  r"-\log_{10}p_{\mathrm{MS}}$      enriched $\rightarrow$",
                  fontsize=11.5)
    ax.spines[["top", "right"]].set_visible(False)


def draw_brain(fig, rect, coords, systems, sysd):
    disp = plot_glass_brain(None, display_mode="lzr", axes=rect, figure=fig)

    def mask(names):
        return np.isin(systems, list(names))

    enr_sig = {s for s, r in sysd.items() if r["enr"] and r["q"] < Q_SIG}
    dep_sig = {s for s, r in sysd.items() if (not r["enr"]) and r["q"] < Q_SIG}
    ns = set(SYS_ORDER) - enr_sig - dep_sig
    m_ns = mask(ns | (set(np.unique(systems)) - set(SYS_ORDER)))
    if m_ns.any():
        disp.add_markers(coords[m_ns], marker_color=C_NS, marker_size=6, alpha=0.20)
    if mask(dep_sig).any():
        disp.add_markers(coords[mask(dep_sig)], marker_color=C_DEP, marker_size=26,
                         alpha=0.9)
    if mask(enr_sig - {WINNER}).any():
        disp.add_markers(coords[mask(enr_sig - {WINNER})], marker_color=C_ENR,
                         marker_size=26, alpha=0.9)
    if WINNER in enr_sig and mask({WINNER}).any():
        disp.add_markers(coords[mask({WINNER})], marker_color=C_ENR, marker_size=72,
                         alpha=0.98, edgecolors="black", linewidths=1.4)
    return disp


def main():
    inc, exc = load_sys("include"), load_sys("exclude")

    fig = plt.figure(figsize=(11.6, 5.7))
    axa = fig.add_axes([0.085, 0.17, 0.44, 0.76])
    draw_lollipop(axa, inc, exc)
    fig.text(0.02, 0.96, r"$\mathbf{a}$", fontsize=17, va="top", fontweight="bold")

    if _BRAIN_OK:
        try:
            coords, systems = pooled_coords_systems()
            draw_brain(fig, (0.56, 0.10, 0.42, 0.82), coords, systems, inc)
            fig.text(0.565, 0.96, r"$\mathbf{b}$", fontsize=17, va="top",
                     fontweight="bold")
            fig.text(0.775, 0.14, r"encoding $\rightarrow$ OFC", ha="center",
                     fontsize=12, color=C_ENR, fontweight="bold")
            print(f"  brain: pooled {len(coords)} contacts")
        except Exception as exc_b:  # pragma: no cover
            print(f"[warn] brain panel failed ({exc_b}); shipping panel a only")

    handles = [
        Line2D([], [], color=C_ENR, marker="o", ls="", mfc=C_ENR, mec="black",
               mew=1.4, ms=11, label=r"enriched, BH $q<0.05$ ($\ast$)"),
        Line2D([], [], color=C_DEP, marker="o", ls="", mfc=C_DEP, mec="black",
               mew=1.1, ms=10, label=r"depleted, BH $q<0.05$ ($\ast$)"),
        Line2D([], [], color="0.4", marker="o", ls="", mfc="white", mec="0.4",
               mew=1.4, ms=9, label=r"n.s. ($q\geq0.05$)"),
        Line2D([], [], color=C_ENR, marker="o", ls="", mfc="none", mec=C_ENR,
               mew=0.9, ms=6, alpha=0.6, label="epilepsy-excluded replicate"),
    ]
    if _BRAIN_OK:
        handles.append(Patch(facecolor=C_ENR, edgecolor="black",
                             label="OFC (encoding carrier)"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.03),
               ncol=5, frameon=False, fontsize=9.2, handletextpad=0.5,
               columnspacing=1.3)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight")
    plt.close(fig)

    print("fig:arc_c — encoding -> OFC anchor (beta)")
    for s in SYS_ORDER:
        r = inc.get(s, {})
        print(f"  {s:16s} {'ENR' if r.get('enr') else 'dep'} "
              f"signed={r.get('signed', float('nan')):+.2f} q={r.get('q', float('nan')):.3f} "
              f"strdev={r.get('strdev', float('nan')):+.2f}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
