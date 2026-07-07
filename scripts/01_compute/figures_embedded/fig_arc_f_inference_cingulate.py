#!/usr/bin/env python3
r"""fig:arc_f — the inference component leans to the cingulate, away from OFC (§2, R2.6).

CORE MESSAGE: encoding anchors in OFC (R2.3); the inference-specific component
(f = D_test - D_learn, encoding partialled out) concentrates elsewhere. In beta its one
cortical concentration at FULL LENGTH is the cingulate -- the only system clearing
correction (BH q = 0.030 include / 0.040 exclude), in nodes of near-average coupling --
while OFC, the encoding anchor, carries none of it. But this is held to a DIRECTION, not a
fixed location: f carries the test/learn length asymmetry, and the cingulate concentration
does NOT survive matching the two recordings' lengths (length-matched q = 0.15 include /
0.30 exclude). So the cingulate is the DIRECTIONAL LEAD for the inferred relations
(consistent sign, away from the orbitofrontal encoding anchor), length-assisted, not an
established location. The provisional cell is drawn with a distinct hatched glyph and its
own legend entry -- never at equal status to the encoding->OFC result.

Reads : data/audit/inference_localization_rhosym/inference_pe_localization_rhosym_{include,exclude}.csv
        data/audit/inference_localization_rhosym/lenmatched_null_rhosym_R200_{include,exclude}.csv
Writes: data/reports/results_section2/fig_arc_f_inference_cingulate.pdf
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
OUT = ROOT / "data/reports/results_section2/fig_arc_f_inference_cingulate.pdf"

DROP = {"other", "non_anatomical"}
SYS_ORDER = ["OFC", "cingulate", "MTL", "insula", "lateral_temporal",
             "PFC", "parietal", "sensorimotor", "occipital"]
WINNER, ANCHOR = "cingulate", "OFC"
C_INF, C_DEP, C_NS = "#d1491c", "#d1791f", "#9a9a9a"
C_INF_PALE = "#f0a688"                    # length-assisted cingulate (provisional)
CAP = 3.0
Q_SIG = 0.05


def load_sys(inc: str) -> dict:
    df = pd.read_csv(BASE / f"inference_pe_localization_rhosym_{inc}.csv")
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


def lenmatched_q(inc: str) -> float:
    df = pd.read_csv(BASE / f"lenmatched_null_rhosym_R200_{inc}.csv")
    row = df[df.unit == WINNER]
    return float(row.bh_q_system.iloc[0]) if not row.empty else float("nan")


def draw_lollipop(ax, inc, exc, lm):
    y0 = np.arange(len(SYS_ORDER))[::-1]
    ypos = {s: y for s, y in zip(SYS_ORDER, y0)}
    ax.axvline(0, color="0.55", lw=0.9, zorder=1)
    for s in SYS_ORDER:
        rec = inc.get(s)
        if rec is None:
            continue
        y, x = ypos[s], rec["signed"]
        prov = (s == WINNER)                          # the length-assisted lead
        col = C_INF if rec["enr"] else C_DEP
        sig = rec["q"] < Q_SIG
        ax.plot([0, x], [y, y], color=col, lw=1.9, alpha=0.9, zorder=2)
        xe = exc.get(s, {}).get("signed")
        if xe is not None:
            ax.scatter([xe], [y + 0.16], s=16, marker="o", facecolors="none",
                       edgecolors=col, linewidths=0.8, alpha=0.55, zorder=3)
        if prov and sig:
            # distinct provisional glyph: hatched diamond + dashed-look outer ring
            ax.scatter([x], [y], s=150, marker="D", facecolors=C_INF_PALE,
                       edgecolors=C_INF, linewidths=1.7, hatch="////", zorder=5)
            ax.scatter([x], [y], s=290, marker="D", facecolors="none",
                       edgecolors=C_INF, linewidths=1.1, alpha=0.55, zorder=5)
        elif sig:
            ax.scatter([x], [y], s=105, marker="o", facecolors=col,
                       edgecolors="black", linewidths=1.1, zorder=4)
            ax.annotate(r"$\ast$", (x, y), textcoords="offset points",
                        xytext=(8 if x >= 0 else -8, 3), ha="center",
                        fontsize=13, color=col, zorder=6)
        else:
            ax.scatter([x], [y], s=70, marker="o", facecolors="white",
                       edgecolors=col, linewidths=1.4, zorder=4)
    # cingulate: full-length vs length-matched q, spelling out the demotion
    wr = inc[WINNER]
    ax.annotate(rf"full-length $q={wr['q']:.2f}$/{exc[WINNER]['q']:.2f}"
                "\n"
                rf"length-matched $q={lm['include']:.2f}$/{lm['exclude']:.2f} (fails)",
                (wr["signed"], ypos[WINNER]), textcoords="offset points",
                xytext=(4, 20), ha="center", fontsize=9.0, color=C_INF)

    ax.set_yticks(y0)
    ax.set_yticklabels(SYS_ORDER)
    for lab in ax.get_yticklabels():
        if lab.get_text() in (WINNER, ANCHOR):
            lab.set_fontweight("bold")
    ax.set_ylim(y0.min() - 0.7, y0.max() + 0.85)
    ax.set_xlim(-CAP - 0.6, CAP + 0.6)
    ax.set_xlabel(r"$\leftarrow$ depleted      signed enrichment  "
                  r"$\mathrm{sign}(M_{\mathrm{obs}}-\tilde{M}_0)\cdot"
                  r"-\log_{10}p_{\mathrm{MS}}$      enriched $\rightarrow$",
                  fontsize=11.5)
    ax.spines[["top", "right"]].set_visible(False)


def draw_brain(fig, rect, coords, systems, sysd):
    # Isolate the cingulate lead: enrichment and depletion share the orange family
    # here (unlike encoding's green/orange), so painting depletions would blur the
    # single message. Cingulate prominent (provisional glyph); everything else faint.
    disp = plot_glass_brain(None, display_mode="lzr", axes=rect, figure=fig)
    m_win = np.isin(systems, [WINNER])
    if (~m_win).any():
        disp.add_markers(coords[~m_win], marker_color=C_NS, marker_size=6, alpha=0.16)
    if WINNER in sysd and sysd[WINNER]["enr"] and sysd[WINNER]["q"] < Q_SIG and m_win.any():
        # provisional: pale face + heavy ring (never the saturated winner treatment)
        disp.add_markers(coords[m_win], marker_color=C_INF_PALE, marker_size=64,
                         alpha=0.98, edgecolors=C_INF, linewidths=1.9)
    return disp


def main():
    inc, exc = load_sys("include"), load_sys("exclude")
    lm = {"include": lenmatched_q("include"), "exclude": lenmatched_q("exclude")}

    fig = plt.figure(figsize=(11.8, 5.8))
    axa = fig.add_axes([0.085, 0.17, 0.44, 0.72])
    draw_lollipop(axa, inc, exc, lm)
    fig.text(0.02, 0.96, r"$\mathbf{a}$", fontsize=17, va="top", fontweight="bold")

    if _BRAIN_OK:
        try:
            coords, systems = pooled_coords_systems()
            draw_brain(fig, (0.56, 0.10, 0.42, 0.82), coords, systems, inc)
            fig.text(0.565, 0.96, r"$\mathbf{b}$", fontsize=17, va="top",
                     fontweight="bold")
            fig.text(0.775, 0.13, r"inference $\rightarrow$ cingulate (directional lead)",
                     ha="center", fontsize=11.5, color=C_INF, fontweight="bold")
        except Exception as exc_b:  # pragma: no cover
            print(f"[warn] brain panel failed ({exc_b}); shipping panel a only")

    handles = [
        Line2D([], [], color=C_INF, marker="D", ls="", mfc=C_INF_PALE, mec=C_INF,
               mew=1.7, ms=12,
               label="inference lead — length-assisted (fails length-matched null)"),
        Line2D([], [], color=C_DEP, marker="o", ls="", mfc=C_DEP, mec="black",
               mew=1.1, ms=10, label=r"depleted, BH $q<0.05$ ($\ast$)"),
        Line2D([], [], color="0.4", marker="o", ls="", mfc="white", mec="0.4",
               mew=1.4, ms=9, label=r"n.s. ($q\geq0.05$) — incl. OFC (encoding anchor)"),
        Line2D([], [], color=C_INF, marker="o", ls="", mfc="none", mec=C_INF,
               mew=0.9, ms=6, alpha=0.6, label="epilepsy-excluded replicate"),
    ]
    if _BRAIN_OK:
        handles.append(Patch(facecolor=C_INF_PALE, edgecolor=C_INF,
                             label="cingulate (provisional lead)"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.04),
               ncol=3, frameon=False, fontsize=9.0, handletextpad=0.5,
               columnspacing=1.2)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight")
    plt.close(fig)

    print("fig:arc_f — inference -> cingulate directional lead (beta)")
    for s in SYS_ORDER:
        r = inc.get(s, {})
        print(f"  {s:16s} {'ENR' if r.get('enr') else 'dep'} "
              f"signed={r.get('signed', float('nan')):+.2f} q={r.get('q', float('nan')):.3f} "
              f"strdev={r.get('strdev', float('nan')):+.2f}")
    print(f"  cingulate length-matched q: include={lm['include']:.3f} "
          f"exclude={lm['exclude']:.3f} (fails -> directional lead)")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
