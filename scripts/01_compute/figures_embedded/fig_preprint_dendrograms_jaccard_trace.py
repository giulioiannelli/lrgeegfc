"""Per-subtree Jaccard trace coloring on the 4-phase dendrogram.

For each internal node S of any phase's dendrogram, score the subtree by
how well it's preserved as a clade across phases — using a phase-specific
formula that asks the trace question directly:

    J_X(S) = max over internal-node subtrees S' of tree X of  |S∩S'|/|S∪S'|

Per-phase score:
    tTest, tLearn :  s(S) =  J_post(S)  -  J_pre(S)              ∈ [-1, +1]
    rPost         :  s(S) =  J_test(S)  -  J_pre(S)              ∈ [-1, +1]
    rPre          :  s(S) = (J_test(S) + J_post(S))/2  -  1      ∈ [-1,  0]
                       0  = anchor (subtree fully preserved into both tasks)
                      -1  = reset (subtree absent in both task phases)

Coloured RdBu_r centred at 0, symmetric vlim across all four phases. Red
on tTest *and* on rPost at the *same* leaf set ⇒ trace clade preserved
between task and post-rest. Blue on rPre at the same leaves ⇒ those
leaves belonged to a different clade in baseline = "freshly assembled by
the task and held".

Single patient × single band (Pat_06 × β) for first iteration.
Output:
    data/reports/preprint/figure1_beta_trace/dendrograms_kc_trace/
        Pat_06_beta_jaccard_trace.pdf
"""
import argparse
from pathlib import Path

import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrgsglib.plotlib import imshow_colorbar_caxdivider

use_lrg_style()


PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
PHASE_TITLES = {
    "rest_pre": r"$\mathrm{rPre}$",
    "task_learn": r"$\mathrm{tLearn}$",
    "task_test": r"$\mathrm{tTest}$",
    "rest_post": r"$\mathrm{rPost}$",
}
FC = "imcoh_abs"

ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "data/reports/preprint/figure1_beta_trace/dendrograms_kc_trace"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def leaves_under_each_node(Z: np.ndarray) -> dict[int, frozenset[int]]:
    n = Z.shape[0] + 1
    out: dict[int, frozenset[int]] = {i: frozenset({i}) for i in range(n)}
    for i in range(n - 1):
        a = int(Z[i, 0])
        b = int(Z[i, 1])
        out[n + i] = out[a] | out[b]
    return out


def best_jaccard(S: frozenset[int], leaves_other: dict[int, frozenset[int]]) -> float:
    """Max Jaccard of S vs any non-singleton subtree in `leaves_other`."""
    best = 0.0
    for T in leaves_other.values():
        if len(T) < 2:
            continue
        inter = len(S & T)
        if inter == 0:
            continue
        j = inter / len(S | T)
        if j > best:
            best = j
    return best


def per_subtree_jaccard(
    leaves_X: dict[int, frozenset[int]],
    leaves_Y: dict[int, frozenset[int]],
) -> dict[int, float]:
    return {
        nid: best_jaccard(S, leaves_Y)
        for nid, S in leaves_X.items()
        if len(S) >= 2
    }


def score_for_phase(
    phase: str,
    leaves_X: dict[int, frozenset[int]],
    leaves_pre: dict[int, frozenset[int]],
    leaves_test: dict[int, frozenset[int]],
    leaves_post: dict[int, frozenset[int]],
) -> dict[int, float]:
    if phase == "rest_pre":
        Jt = per_subtree_jaccard(leaves_X, leaves_test)
        Jp = per_subtree_jaccard(leaves_X, leaves_post)
        return {nid: 0.5 * (Jt[nid] + Jp[nid]) - 1.0 for nid in Jt}
    elif phase == "rest_post":
        Jt = per_subtree_jaccard(leaves_X, leaves_test)
        Je = per_subtree_jaccard(leaves_X, leaves_pre)
        return {nid: Jt[nid] - Je[nid] for nid in Jt}
    else:  # task_test, task_learn — both compared as: post-side − pre
        Jp = per_subtree_jaccard(leaves_X, leaves_post)
        Je = per_subtree_jaccard(leaves_X, leaves_pre)
        return {nid: Jp[nid] - Je[nid] for nid in Jp}


def main(patient: str, band: str) -> None:
    out_pdf = OUT_DIR / f"{patient}_{band}_jaccard_trace.pdf"
    Zs = {p: load_lrg_result(patient, p, band, FC).linkage_matrix for p in PHASES}
    n = Zs[PHASES[0]].shape[0] + 1
    leaves = {p: leaves_under_each_node(Z) for p, Z in Zs.items()}

    scores = {
        p: score_for_phase(
            p,
            leaves[p],
            leaves["rest_pre"],
            leaves["task_test"],
            leaves["rest_post"],
        )
        for p in PHASES
    }

    # Sanity prints
    for p in PHASES:
        v = np.fromiter(scores[p].values(), dtype=float)
        print(f"{p:12s}  n_subtrees={v.size:3d}  score min={v.min():+.3f}  "
              f"max={v.max():+.3f}  median={np.median(v):+.3f}")

    # Symmetric vlim across all phases
    all_v = np.concatenate([np.fromiter(s.values(), dtype=float) for s in scores.values()])
    vmax = max(float(np.max(np.abs(all_v))), 0.05)
    vmin = -vmax
    print(f"\nsymmetric vlim = ±{vmax:.3f}")

    cmap = plt.get_cmap("RdBu_r")
    norm = Normalize(vmin=vmin, vmax=vmax)

    def color_for_phase(phase: str):
        sc = scores[phase]

        def cf(node_id: int) -> str:
            v = sc.get(node_id)
            if v is None or not np.isfinite(v):
                return "#bbbbbb"
            return mcolors.to_hex(cmap(norm(v)))

        return cf

    tmin = min(Z[0, 2] for Z in Zs.values()) * 0.8
    tmax = max(Z[-1, 2] for Z in Zs.values()) * 1.05

    fig, axes = plt.subplots(
        1, 4, figsize=(13.5, 3.8),
        sharey=True,
        gridspec_kw={"width_ratios": [1, 1, 1, 1.10]},
    )
    fig.subplots_adjust(left=0.06, right=0.94, top=0.84, bottom=0.06, wspace=0.10)

    for col, (ax, phase) in enumerate(zip(axes, PHASES)):
        ax.set_title(PHASE_TITLES[phase])
        Z = Zs[phase]
        dendrogram(
            Z, ax=ax,
            no_labels=True,
            link_color_func=color_for_phase(phase),
            above_threshold_color="#bbbbbb",
        )
        ax.set_ylim(tmin, tmax)
        ax.set_xticks([])
        ax.set_xlabel("")
        if col == 0:
            ax.set_ylabel(r"merge height  $\hat{D}(\tau)$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    _, _, cb = imshow_colorbar_caxdivider(sm, axes[-1], size="4%", pad=0.10)
    cb.set_label(
        r"per-subtree Jaccard trace score"
        "\n"
        r"$+1$: perfectly preserved in partner; $-1$: absent in partner",
        fontsize=7,
    )

    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    fig.text(
        0.06, 0.94,
        rf"{patient}  |  {band_tex}  |  Jaccard trace coloring  |  $N = {n}$",
        fontweight="bold",
    )
    fig.text(
        0.06, 0.91,
        "red on tTest and rPost on the same leaves = trace clade  |  "
        "red on rPre = anchor  |  blue on rPre = reset",
        fontsize=7, color="#444",
    )

    fig.savefig(out_pdf)
    plt.close(fig)
    print(f"\nDONE: {out_pdf}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--patient", default="Pat_06")
    parser.add_argument("--band", default="beta")
    args = parser.parse_args()
    main(args.patient, args.band)
