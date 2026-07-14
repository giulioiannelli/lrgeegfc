#!/usr/bin/env python3
r"""fig:cophenetic-tanglegram — what ρ^coph actually *sees*: the split hierarchy of one small
clade, drawn as a three-panel tanglegram across rest_pre → task → rest_post, for an ANCHOR,
a TRACE and a RESET exemplar in ONE β patient (§5.3 explainer / methods-intuition figure).

CORE MESSAGE. ρ^coph = Spearman of the cophenetic-distance vectors between two phases, i.e. it
scores whether the RANK ORDER OF SPLITS (the height at which each contact-pair first joins) is
preserved. This figure makes that literal: crop the huge full-patient dendrogram down to one
~8-contact clade and draw its induced sub-tree in each phase, linking the same contact across
adjacent panels. Lines that stay parallel = the grouping was kept; lines that cross = it
reorganised. Contacts are coloured by their REST_POST sub-clade, so the three fates read as a
colour-block pattern before you count a single line:

    ANCHOR : blocks contiguous in all three rails         → nothing moved
    TRACE  : blocks scrambled in rest_pre, contiguous in task AND rest_post
                                                          → written at task, held offline
    RESET  : blocks contiguous in rest_pre AND rest_post, scrambled in task
                                                          → moved at task, reverted

────────────────────────────────────────────────────────────────────────────────────────────
CRITICAL PREAMBLE (5 points, per house rule — this is an EXPLAINER, read the honesty).
  1. Claim (illustrative). On a chosen clade, ρ^coph's three regimes (high/high, low/high,
     "outer-high/inner-low") correspond to anchor / trace / reset; the tanglegram visualises
     the split-rank (dis)agreement that ρ^coph summarises.
  2. Null. NONE — no gate here. The cohort claim stays the matched-strength gates elsewhere
     (β cophenetic trace 23.7× MS, q .009–.013 → OFC; inference gate T_infspec_pe β p .0098,
     6/10). This figure explains the measure; it does not test it.
  3. Strongest way it could mislead. Dendrogram leaf order is arbitrary up to branch flips, so
     naive panels could FABRICATE or HIDE crossings. Mitigation: every tree is UNTANGLED —
     greedy branch-rotation (topology-preserving) to a rest_pre→task→rest_post reference chain
     — so residual crossings are genuine topological disagreement, not a drawing artefact.
  4. What it cannot show. Subset ρ^coph (annotated per gap) is LOCAL and ≠ the whole-tree
     ρ^coph that the paper reports; it is a single patient; clades are picked by a DISCLOSED
     objective score on that subset ρ. It is an exemplar, not cohort evidence.
  5. Falsify / limits. A "trace" clade whose task→rest_post gap untangles to ~0 crossings AND
     high ρ is consistent; if it did not, the label would be wrong. Crossing counts and ρ are
     logged together so the two must agree. Reset "reverts" is shown by ρ(pre,post) high with
     both inner gaps crossing — annotated, not asserted.
────────────────────────────────────────────────────────────────────────────────────────────

Clade selection (disclosed). Enumerate every clade of the REST_POST tree with 6–10 contacts;
with a=ρ(pre,task) b=ρ(task,post) c=ρ(pre,post) on that subset, score
  anchor = min(a,b,c)        trace = b − max(a,c)        reset = c − max(a,b).
Best clade per fate; the patient is the one maximising min(anchor,trace,reset) so all three
exemplars are decent (override with --patient).

Reuses coph_square / load_phase_fc / lrg_linkage (DRY with fig_branch_origin_tree).
Writes: data/preprint/figures/_drafts/fig_cophenetic_tanglegram_{pat}_{band}_DRAFT.pdf
"""
from __future__ import annotations

import argparse
import sys
import time

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Patch
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgb
from scipy.cluster.hierarchy import linkage, dendrogram, optimal_leaf_ordering
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.utils.probe import contact_labels, split_label
from lrg_eegfc.visuals.network_templates import load_probe_labels

sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "figures_embedded"))
from audit_63_split_baseline_surrogate import load_phase_fc              # type: ignore
from audit_65_kc_matched_strength_surrogate import lrg_linkage           # type: ignore  (noqa)
from fig_branch_origin_tree import coph_square, BAND_SYM                 # type: ignore

OUTDIR = ROOT / "data/preprint/figures/_drafts"
QA = ("/tmp/claude-1000/-home-giulio-Documents-research-neural-networks-lrgeegfc/"
      "c2197e94-c496-4e4e-97cc-89e43564e7da/scratchpad")
COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]

RAIL = (0.0, 1.0, 2.0)                 # x of the three phase leaf-rails
EXT = 0.34                             # max horizontal extent of a dendrogram (outer phases)
EXT_MID = 0.55 * EXT                   # shorter for the middle tree (shares gap 1)
GROUP_COLORS = ["#3b6db3", "#e0842b", "#4d9221", "#8073ac",   # rest_post sub-clade tints
                "#c44e52", "#4bb2b6"]
BRACKET = "#454a52"


def n_groups(m):
    """Number of colour groups: 3 for small clades, up to 5 as the clade grows (~1 per 5 leaves)."""
    return int(np.clip(m // 5, 3, 5))
FATES = ("anchor", "trace", "reset", "reorganized")
ROW_META = {                           # fate -> (title, one-line descriptor)
    "anchor":      ("ANCHOR",      "grouping unchanged — held across every phase"),
    "trace":       ("TRACE",       "reorganised at task, then retained after"),
    "reset":       ("RESET",       "reorganised at task, then reverted to baseline"),
    "reorganized": ("REORGANISED", "reorganised at every transition — never settles"),
}

# labeled-bead layout: a fixed physical bead + a fixed inches-per-y-unit so beads never
# overlap regardless of clade size N (figure height scales with N via GridSpec ratios).
UNIT_IN = 0.30                         # inches per dendrogram y-unit
NODE_S = 178.0                         # bead marker area (pt²) ≈ 15 pt diameter
NODE_LW = 1.0                          # bead border width
LBL_FS = 5.4                           # in-bead contact-label font size
FIG_W = 9.6
SPAN_TOP, SPAN_MID, BOTTOM_IN = 3.6, 2.4, 1.35   # extra y-units per row / bottom margin
_SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")


def bead_label(name):
    """Contact name as a compact in-bead string: 'A10' -> 'A¹⁰' (shaft + superscript index)."""
    shaft, idx = split_label(name)
    return f"{shaft}{idx.translate(_SUP)}" if idx else shaft


def text_ink(color):
    """Black or white ink for a bead fill, by perceived luminance."""
    r, g, b = to_rgb(color)
    return "#0f1216" if 0.299 * r + 0.587 * g + 0.114 * b > 0.62 else "#ffffff"


def contact_names(pat, n):
    """FC-aligned contact names for a patient, or numeric fallback on any load failure."""
    try:
        names = contact_labels(load_probe_labels(pat), n)
        if len(names) >= n:
            return names
    except Exception:                                            # noqa: BLE001
        pass
    return [str(i) for i in range(n)]


# ───────────────────────── cophenetic / clade machinery ─────────────────────────
def load_cophs(pat, band, task_phase):
    """Per-phase N×N cophenetic matrices (aligned leaves) + rest_post linkage. None on mismatch."""
    phases = ["rest_pre", task_phase, "rest_post"]
    C, Zpost, ns = {}, None, []
    for ph in phases:
        W = load_phase_fc(pat, ph, band)
        Cq, Z = coph_square(W)
        C[ph] = Cq
        ns.append(W.shape[0])
        if ph == "rest_post":
            Zpost = Z
    if len(set(ns)) != 1:
        return None
    return C, Zpost


def clades(Z, N, lo, hi):
    """Leaf-index lists of every rest_post clade with size in [lo, hi]."""
    members = {i: [i] for i in range(N)}
    out = []
    for k in range(N - 1):
        m = members[int(Z[k, 0])] + members[int(Z[k, 1])]
        members[N + k] = m
        if lo <= len(m) <= hi:
            out.append(sorted(m))
    return out


def subset_rho(Ca, Cb, S):
    """Spearman of the upper-triangular induced cophenetic distances on clade S."""
    iu = np.triu_indices(len(S), 1)
    a = Ca[np.ix_(S, S)][iu]
    b = Cb[np.ix_(S, S)][iu]
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return np.nan
    return float(spearmanr(a, b).statistic)


def abc(C, S, task_phase):
    a = subset_rho(C["rest_pre"], C[task_phase], S)
    b = subset_rho(C[task_phase], C["rest_post"], S)
    c = subset_rho(C["rest_pre"], C["rest_post"], S)
    return a, b, c


def scores(a, b, c):
    """a=ρ(pre,task) b=ρ(task,post) c=ρ(pre,post). Higher score = cleaner exemplar of the fate.
    anchor=all high; trace=held by task not baseline; reset=reverted to baseline; reorganized=all
    low (changes at every transition)."""
    return dict(anchor=min(a, b, c), trace=b - max(a, c),
                reset=c - max(a, b), reorganized=-max(a, b, c))


# ───────────────────────── induced sub-tree + untangling ─────────────────────────
def induced_linkage(C, S):
    """Average-linkage on the induced cophenetic sub-matrix == the full tree's induced sub-tree
    (a sub-matrix of an ultrametric is ultrametric, so UPGMA recovers it exactly)."""
    Cs = C[np.ix_(S, S)].astype(float)
    Cs = 0.5 * (Cs + Cs.T)
    np.fill_diagonal(Cs, 0.0)
    return linkage(squareform(np.clip(Cs, 0, None), checks=False), method="average")


def untangle(Z, m, refpos):
    """Rotation-only leaf order of Z closest to reference positions refpos[leaf] (greedy
    branch flips: at each merge the child with the smaller mean reference position goes first)."""
    def rec(node):
        if node < m:
            return [node]
        k = node - m
        L = rec(int(Z[k, 0])); R = rec(int(Z[k, 1]))
        return L + R if np.mean([refpos[i] for i in L]) <= np.mean([refpos[i] for i in R]) else R + L
    return rec(2 * m - 2)


def subtree_leaves(Z, m, node):
    if node < m:
        return [node]
    k = node - m
    return subtree_leaves(Z, m, int(Z[k, 0])) + subtree_leaves(Z, m, int(Z[k, 1]))


def top_k_groups(Z, m, K):
    """Colour groups = the K highest sub-clades of the rest_post induced tree (successive split
    of the current tallest cluster). Cosmetic grouping for the tanglegram, NOT a metric."""
    cl = [2 * m - 2]
    while len(cl) < K:
        internal = [c for c in cl if c >= m]
        if not internal:
            break
        c = max(internal, key=lambda nd: Z[nd - m, 2])
        cl.remove(c)
        k = c - m
        cl += [int(Z[k, 0]), int(Z[k, 1])]
    grp = {}
    for gi, c in enumerate(cl):
        for l in subtree_leaves(Z, m, c):
            grp[l] = gi
    return grp


def crossings(pa, pb, m):
    return sum(1 for i in range(m) for j in range(i + 1, m)
               if (pa[i] - pa[j]) * (pb[i] - pb[j]) < 0)


# ───────────────────────── drawing ─────────────────────────
def ribbon(ax, x0, y0, x1, y1, color, alpha=0.72, lw=1.7, z=2):
    dx = x1 - x0
    verts = [(x0, y0), (x0 + 0.42 * dx, y0), (x1 - 0.42 * dx, y1), (x1, y1)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
    ax.add_patch(PathPatch(Path(verts, codes), fc="none", ec=color, lw=lw,
                           alpha=alpha, zorder=z, capstyle="round"))


def draw_subtree(ax, Z, m, ypos, x_rail, direction, ext, hlo, hhi):
    def ext_of(h):
        if h <= 0:
            return 0.0
        t = (np.log(h) - hlo) / (hhi - hlo) if hhi > hlo else 1.0
        return direction * ext * float(np.clip(t, 0.06, 1.0))

    def rec(node):
        if node < m:
            return ypos[node], 0.0
        k = node - m
        y1, h1 = rec(int(Z[k, 0]))
        y2, h2 = rec(int(Z[k, 1]))
        xp = x_rail + ext_of(Z[k, 2])
        x1 = x_rail + ext_of(h1)
        x2 = x_rail + ext_of(h2)
        for xa, xb, ya, yb in ((xp, xp, y1, y2), (xp, x1, y1, y1), (xp, x2, y2, y2)):
            ax.plot([xa, xb], [ya, yb], color=BRACKET, lw=1.2, solid_capstyle="round", zorder=3)
        return (y1 + y2) / 2.0, Z[k, 2]
    rec(2 * m - 2)


_PHASE_LABEL = {"rest_pre": "rest-pre", "rest_post": "rest-post"}
_PHASE_SHORT = {"rest_pre": "pre", "rest_post": "post"}


def _ph_label(ph, task_phase):
    return f"task ({task_phase.split('_')[-1]})" if ph == task_phase else _PHASE_LABEL.get(ph, ph)


def _ph_short(ph, task_phase):
    return "task" if ph == task_phase else _PHASE_SHORT.get(ph, ph)


def build_row(ax, C, S, task_phase, is_top, labels=None, phase_order=None,
              color_phase="rest_post", show_outer_rho=True,
              rail=None, ext=None, ext_mid=None):
    """Three-panel tanglegram of clade S. `phase_order` = (left, mid, right) phase keys;
    default is temporal (rest_pre, task, rest_post). `color_phase`'s sub-clades set the bead
    colours — its blocks stay contiguous under any branch rotation, so THAT column reads as
    clean colour blocks; put the reference phase there. For a RESET we pass
    (rest_pre, rest_post, task): the reverted pair sits ADJACENT so gap 1 = pre↔post is the
    similar/parallel one and the task is the lone scrambled departure on the right.

    `rail` (3 x-positions of the leaf-rails), `ext` (outer-tree horizontal extent) and `ext_mid`
    (middle-tree extent) default to the module constants; pass a tighter `rail` + larger `ext` to
    give the dendrograms more width and the cross-phase ribbons less (e.g. the side-by-side quad)."""
    m = len(S)
    rail = RAIL if rail is None else rail
    ext = EXT if ext is None else ext
    ext_mid = EXT_MID if ext_mid is None else ext_mid
    if phase_order is None:
        phase_order = ("rest_pre", task_phase, "rest_post")
    p0, p1, p2 = phase_order
    Z = {ph: induced_linkage(C[ph], S) for ph in phase_order}

    # reference chain: left panel keeps its own optimal order, comb rightward
    Cs0 = C[p0][np.ix_(S, S)].astype(float)
    Z[p0] = optimal_leaf_ordering(Z[p0], squareform(0.5 * (Cs0 + Cs0.T)
                                                    - np.diag(np.diag(Cs0)), checks=False))
    pos0 = {l: i for i, l in enumerate(dendrogram(Z[p0], no_plot=True)["leaves"])}
    pos1 = {l: i for i, l in enumerate(untangle(Z[p1], m, pos0))}
    pos2 = {l: i for i, l in enumerate(untangle(Z[p2], m, pos1))}
    groups = top_k_groups(Z[color_phase], m, n_groups(m))

    allh = np.concatenate([Z[p][:, 2] for p in phase_order])
    posh = allh[allh > 0]
    hlo, hhi = float(np.log(posh.min())), float(np.log(posh.max()))

    draw_subtree(ax, Z[p0], m, pos0, rail[0], -1, ext, hlo, hhi)
    draw_subtree(ax, Z[p1], m, pos1, rail[1], -1, ext_mid, hlo, hhi)
    draw_subtree(ax, Z[p2], m, pos2, rail[2], +1, ext, hlo, hhi)

    for l in range(m):
        col = GROUP_COLORS[groups[l] % len(GROUP_COLORS)]
        ribbon(ax, rail[0], pos0[l], rail[1], pos1[l], col)
        ribbon(ax, rail[1], pos1[l], rail[2], pos2[l], col)
    for l in range(m):
        col = GROUP_COLORS[groups[l] % len(GROUP_COLORS)]
        ys = [pos0[l], pos1[l], pos2[l]]
        ax.scatter(list(rail), ys, s=NODE_S, color=col, ec="#20242a", lw=NODE_LW, zorder=5)
        if labels is not None:
            txt, ink = bead_label(labels[l]), text_ink(col)
            for x, y in zip(rail, ys):
                ax.text(x, y, txt, ha="center", va="center", fontsize=LBL_FS,
                        color=ink, fontweight="bold", zorder=6)

    a, b, c = abc(C, S, task_phase)                    # canonical (pre,task / task,post / pre,post)
    rho01 = subset_rho(C[p0], C[p1], S)                # displayed gaps follow the column order
    rho12 = subset_rho(C[p1], C[p2], S)
    rho02 = subset_rho(C[p0], C[p2], S)
    x1c = crossings([pos0[l] for l in range(m)], [pos1[l] for l in range(m)], m)
    x2c = crossings([pos1[l] for l in range(m)], [pos2[l] for l in range(m)], m)
    for xm, val, cross in ((0.5 * (rail[0] + rail[1]), rho01, x1c),
                           (0.5 * (rail[1] + rail[2]), rho12, x2c)):
        held = val >= 0.6
        ax.text(xm, m + 0.42, rf"$\rho={val:.2f}$", ha="center", va="bottom", fontsize=9.8,
                color=("#2f6b34" if held else "#9a3b3b"), fontweight="bold")
        ax.text(xm, m + 0.30, f"{cross} crossing{'s' if cross != 1 else ''}", ha="center",
                va="top", fontsize=7.0, color=("#2f6b34" if held else "#9a3b3b"))
    if show_outer_rho:
        ax.text(rail[2] + ext + 0.12, -0.82,
                rf"$\rho_{{\rm {_ph_short(p0, task_phase)},{_ph_short(p2, task_phase)}}}={rho02:.2f}$",
                ha="right", va="top", fontsize=8, color="#555b63")

    if is_top:
        for x, ph in zip(rail, phase_order):
            ax.text(x, m + 1.35, _ph_label(ph, task_phase), ha="center", va="bottom",
                    fontsize=11, color="#2b2f36", fontweight="bold")

    ax.set_xlim(rail[0] - ext - 0.18, rail[2] + ext + 0.55)
    ax.set_ylim(-1.15, m + (2.05 if is_top else 1.05))
    ax.axis("off")
    return a, b, c, x1c, x2c, m


# ───────────────────────── patient / clade selection ─────────────────────────
def best_clades(C, Zpost, task_phase, lo, hi):
    N = Zpost.shape[0] + 1
    best = {k: (-np.inf, None, None) for k in FATES}
    for S in clades(Zpost, N, lo, hi):
        a, b, c = abc(C, S, task_phase)
        if not np.all(np.isfinite([a, b, c])):
            continue
        sc = scores(a, b, c)
        for fate, v in sc.items():
            if v > best[fate][0]:
                best[fate] = (v, S, (a, b, c))
    return best


def pick_patient(patients, band, task_phase, lo, hi):
    print(f"[scan] ranking {len(patients)} patients for a clean anchor+trace+reset triad "
          f"({band}, task={task_phase}, clade {lo}-{hi})", flush=True)
    rows = []
    for i, pat in enumerate(patients, 1):
        try:
            got = load_cophs(pat, band, task_phase)
        except Exception as e:                                   # noqa: BLE001
            print(f"  [{i}/{len(patients)}] {pat}: skip ({e})", flush=True)
            continue
        if got is None:
            print(f"  [{i}/{len(patients)}] {pat}: skip (phase N mismatch)", flush=True)
            continue
        C, Zpost = got
        best = best_clades(C, Zpost, task_phase, lo, hi)
        bal = min(best[k][0] for k in best)
        rows.append((bal, pat, best))
        print(f"  [{i}/{len(patients)}] {pat}: anchor={best['anchor'][0]:.2f} "
              f"trace={best['trace'][0]:.2f} reset={best['reset'][0]:.2f} "
              f"reorg={best['reorganized'][0]:.2f}  min={bal:.2f}", flush=True)
    rows.sort(reverse=True)
    return rows


def legend_footnotes(fig, b_frac, n=3):
    """Shared bottom-band decoration: sub-clade colour key (n patches) + the two honesty footnotes."""
    handles = [Patch(fc=GROUP_COLORS[i], label=f"rest_post sub-clade {i + 1}") for i in range(n)]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, b_frac * 0.62), ncol=n,
               frameon=False, fontsize=8.6, handlelength=1.4, columnspacing=1.7)
    fig.text(0.5, b_frac * 0.42, "colour = rest_post sub-clade (grouping key)   ·   lines link "
             "the same contact across phases:  parallel = retained,  crossing = reorganised",
             ha="center", va="top", fontsize=7.6, color="#6a7078")
    fig.text(0.5, b_frac * 0.18, r"$\rho$ = local (subset) cophenetic Spearman — illustrative; "
             r"trees untangled by branch rotation.   Cohort claim: matched-strength gate.",
             ha="center", va="top", fontsize=7.4, color="#9aa0a8")


def draw(pat, band, task_phase, best):
    C = load_cophs.cache[pat]
    names = contact_names(pat, C["rest_post"].shape[0])
    ms = [len(best[f][1]) for f in FATES]
    spans = [m + (SPAN_TOP if i == 0 else SPAN_MID) for i, m in enumerate(ms)]
    fig_h = sum(spans) * UNIT_IN + BOTTOM_IN

    use_lrg_style()
    fig = plt.figure(figsize=(FIG_W, fig_h))
    b_frac = BOTTOM_IN / fig_h                          # reserved bottom band for legend/footnotes
    gs = GridSpec(len(FATES), 1, height_ratios=spans, hspace=0.16,
                  top=0.995, bottom=b_frac, figure=fig)
    for r, fate in enumerate(FATES):
        ax = fig.add_subplot(gs[r])
        _, S, _ = best[fate]
        lab = [names[l] for l in S]
        a, b, c, x1, x2, m = build_row(ax, C, S, task_phase, is_top=(r == 0), labels=lab)
        title, desc = ROW_META[fate]
        ax.text(-0.012, 0.5, title, transform=ax.transAxes, rotation=90, ha="right", va="center",
                fontsize=12, fontweight="bold", color="#2b2f36")
        ax.text(RAIL[0], -0.82, desc, ha="left", va="top", fontsize=8.3,
                color="#6a7078", style="italic")
        print(f"  {fate:11s} S={S}  a={a:.2f} b={b:.2f} c={c:.2f}  "
              f"crossings pre→task={x1} task→post={x2}", flush=True)

    legend_footnotes(fig, b_frac)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / f"fig_cophenetic_tanglegram_{pat}_{band}_DRAFT.pdf"
    fig.savefig(out, bbox_inches="tight", transparent=True, pad_inches=0.04)
    fig.savefig(f"{QA}/tanglegram_{pat}_{band}.png", dpi=155, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print(f"wrote {out.name}", flush=True)


def main():
    ap = argparse.ArgumentParser(description="cophenetic tanglegram taxonomy explainer")
    ap.add_argument("--patient", default=None, help="force a patient; default = scan cohort")
    ap.add_argument("--band", default="beta")
    ap.add_argument("--task-phase", default="task_test")
    ap.add_argument("--size-min", type=int, default=6)
    ap.add_argument("--size-max", type=int, default=10)
    a = ap.parse_args()

    t0 = time.perf_counter()
    load_cophs.cache = {}                                        # memoise per-patient cophenetics

    def cached(pat):
        if pat not in load_cophs.cache:
            load_cophs.cache[pat] = load_cophs(pat, a.band, a.task_phase)
        return load_cophs.cache[pat]

    if a.patient:
        got = cached(a.patient)
        if got is None:
            print(f"!! {a.patient}: phase N mismatch — cannot draw")
            return
        best = best_clades(got[0], got[1], a.task_phase, a.size_min, a.size_max)
        load_cophs.cache[a.patient] = got[0]                     # draw() wants C only
        draw(a.patient, a.band, a.task_phase, best)
    else:
        rows = pick_patient(COHORT, a.band, a.task_phase, a.size_min, a.size_max)
        if not rows:
            print("!! no usable patient")
            return
        bal, pat, best = rows[0]
        print(f"[pick] {pat} (balanced min score {bal:.2f})", flush=True)
        load_cophs.cache = {pat: cached(pat)[0]}
        draw(pat, a.band, a.task_phase, best)
    print(f"[done] {time.perf_counter() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
