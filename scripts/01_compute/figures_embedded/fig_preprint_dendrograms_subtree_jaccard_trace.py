"""Subtree clade-preservation trace measure — direct answer to:

    "Which subtrees of each phase's dendrogram are more similar to
     tTest/rPost than to rPre, as clades?"

For every internal node S of each phase X's tree (with leaf set L(S)):

    J_{X→Y}(S) = max over internal subtrees T of phase Y of
                 |L(S) ∩ L(T)| / |L(S) ∪ L(T)|

The trace question is asymmetric in the PHASE the subtree comes from.
Each phase's formula excludes the trivial self-match:

    tTest  →  trace = J_{test→post}(S)  −  J_{test→pre}(S)   ∈ [-1, 1]
    rPost  →  trace = J_{post→test}(S)  −  J_{post→pre}(S)   ∈ [-1, 1]
    tLearn →  trace = ½(J_{learn→test} + J_{learn→post})
                     −  J_{learn→pre}                         ∈ [-1, 1]
    rPre   →  not applicable — rPre IS the baseline → drawn gray.

No self-match (so no `[-1, 0]` floor on rPre or `[-½, 1]` ceiling on
tTest/rPost). All three colored phases share the same `[-1, +1]` range.

Reading:
- **red U-shape:** S's leaf-set forms a clade in the partner phase but
  not in rPre  →  trace clade.
- **blue U-shape:** S's leaf-set forms a clade in rPre but not in the
  partner phase.
- **near-white:** clade preserved in both pre and partner (anchor), or
  in neither.
- **gray rPre panel:** baseline reference; trace measure does not apply
  to subtrees of the baseline by construction.

Difference from `kc_lambda0_trace`:

- KC trace measures **pair-wise tree-depth** similarity averaged over
  the subtree (per-pair s on the m-vector). A subtree can score red
  even if its exact clade isn't preserved, as long as its leaves sit
  at similar MRCA depths in test/post.
- This script measures **clade preservation itself**: the subtree's
  leaf-set must (approximately) form a clade in the partner tree to
  score red. Stricter, more structural.

Per-leaf identity stubs: same probe-shaft coloring as in
`fig_preprint_dendrograms_kc_trace.py` — leaves on the same sEEG shaft
share a stub color across all four panels.

Topology-only by construction (Jaccard is set-overlap). A height-aware
variant for KC λ=0.5 analog is not produced by this script.

One PDF per band, one page per patient × 4 phases.

Output:
    data/reports/preprint/figure1_beta_trace/dendrograms_subtree_jaccard_trace/
        subtree_jaccard_trace_{band}_all_patients.pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.config.const import BRAIN_BANDS, BRAIN_BAND_TEX_DICT
from lrg_eegfc.utils.io.patient import load_channel_labels
from lrg_eegfc.utils.probe import extract_probe_labels
from lrg_eegfc.visuals.styles import use_lrg_style
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrgsglib.plotlib import imshow_colorbar_caxdivider

use_lrg_style()


COHORT = ["Pat_02", "Pat_03", "Pat_05", "Pat_06", "Pat_07",
          "Pat_08", "Pat_10", "Pat_13", "Pat_14", "Pat_15"]
PHASES = ["rest_pre", "task_learn", "task_test", "rest_post"]
PHASE_TITLES = {
    "rest_pre": r"$\mathrm{rPre}$",
    "task_learn": r"$\mathrm{tLearn}$",
    "task_test": r"$\mathrm{tTest}$",
    "rest_post": r"$\mathrm{rPost}$",
}
FC = "imcoh_abs"

ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "data/reports/preprint/figure1_beta_trace/dendrograms_subtree_jaccard_trace"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def leaves_under_each_node(Z: np.ndarray) -> dict[int, frozenset[int]]:
    n = Z.shape[0] + 1
    out: dict[int, frozenset[int]] = {i: frozenset({i}) for i in range(n)}
    for i in range(n - 1):
        a = int(Z[i, 0])
        b = int(Z[i, 1])
        out[n + i] = out[a] | out[b]
    return out


def best_jaccard(S: frozenset[int],
                 leaves_other: dict[int, frozenset[int]]) -> float:
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


def render_patient_page(patient: str, band: str, pdf: PdfPages) -> None:
    try:
        Zs = {p: load_lrg_result(patient, p, band, FC).linkage_matrix for p in PHASES}
    except Exception as err:
        print(f"  [{patient} {band}] skip ({err})")
        return
    n = Zs[PHASES[0]].shape[0] + 1
    leaves_per_phase = {p: leaves_under_each_node(Z) for p, Z in Zs.items()}

    # J[px][py] = per-subtree Jaccard of phase px subtrees against
    # phase py. Self-matches (px == py) are NOT computed: every formula
    # below uses only cross-phase comparisons, so the trivial self-match
    # never enters the trace measure.
    J: dict[str, dict[str, dict[int, float]]] = {}
    for px in PHASES:
        J[px] = {}
        for py in PHASES:
            if px == py:
                continue
            J[px][py] = per_subtree_jaccard(
                leaves_per_phase[px], leaves_per_phase[py]
            )

    # Per-phase trace formula (cross-phase only — no self-match bias):
    #   tTest  : J_test→post  − J_test→pre
    #   rPost  : J_post→test  − J_post→pre
    #   tLearn : ½(J_learn→test + J_learn→post) − J_learn→pre
    #   rPre   : baseline reference, no trace value.
    trace_per_phase: dict[str, dict[int, float]] = {}
    valid_nids_test = J["task_test"]["rest_post"].keys()
    trace_per_phase["task_test"] = {
        nid: J["task_test"]["rest_post"][nid] - J["task_test"]["rest_pre"][nid]
        for nid in valid_nids_test
    }
    valid_nids_post = J["rest_post"]["task_test"].keys()
    trace_per_phase["rest_post"] = {
        nid: J["rest_post"]["task_test"][nid] - J["rest_post"]["rest_pre"][nid]
        for nid in valid_nids_post
    }
    valid_nids_learn = J["task_learn"]["task_test"].keys()
    trace_per_phase["task_learn"] = {
        nid: 0.5 * (J["task_learn"]["task_test"][nid]
                    + J["task_learn"]["rest_post"][nid])
        - J["task_learn"]["rest_pre"][nid]
        for nid in valid_nids_learn
    }
    # rPre intentionally absent — colored gray below.

    all_traces = np.fromiter(
        (v for d in trace_per_phase.values() for v in d.values()),
        dtype=float,
    )
    vmax = max(float(np.max(np.abs(all_traces))), 0.1)
    cmap = plt.get_cmap("RdBu_r")
    norm = Normalize(vmin=-vmax, vmax=vmax)

    print(f"  [{patient} {band}] trace ranges (cross-phase only): "
          + ", ".join(
              f"{p}=[{min(trace_per_phase[p].values()):+.2f},{max(trace_per_phase[p].values()):+.2f}]"
              for p in trace_per_phase
          )
          + f"  sym vlim ±{vmax:.2f}  (rPre = baseline, gray)")

    def color_for_phase(phase: str):
        if phase == "rest_pre":
            def cf_gray(_node_id: int) -> str:
                return "#bbbbbb"
            return cf_gray

        scores = trace_per_phase[phase]

        def cf(node_id: int) -> str:
            sc = scores.get(node_id)
            if sc is None or not np.isfinite(sc):
                return "#888888"
            return mcolors.to_hex(cmap(norm(sc)))

        return cf

    # Probe-shaft leaf identity stubs (consistent across panels).
    channel_labels = load_channel_labels(patient)
    probe_per_leaf = extract_probe_labels(channel_labels)
    unique_probes = sorted(set(probe_per_leaf))
    palette = plt.get_cmap("tab20")
    probe_to_color = {
        pr: mcolors.to_hex(palette(i % palette.N))
        for i, pr in enumerate(unique_probes)
    }
    leaf_anchor_color = [probe_to_color[pr] for pr in probe_per_leaf]

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
        d = dendrogram(
            Z, ax=ax,
            no_labels=True,
            link_color_func=color_for_phase(phase),
            above_threshold_color="#888888",
        )
        leaf_order = d["leaves"]
        leaf_xs = [5 + 10 * i for i in range(len(leaf_order))]
        stub_top = tmin + 0.025 * (tmax - tmin)
        for x, L in zip(leaf_xs, leaf_order):
            ax.plot(
                [x, x], [tmin, stub_top],
                color=leaf_anchor_color[L], linewidth=2.0,
                solid_capstyle="butt", zorder=20,
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
        "per-subtree clade-preservation trace\n"
        r"tTest: $J_{\to\mathrm{post}}-J_{\to\mathrm{pre}}$;  "
        r"rPost: $J_{\to\mathrm{test}}-J_{\to\mathrm{pre}}$;  "
        r"tLearn: $\bar{J}_{\to\{\mathrm{test,post}\}}-J_{\to\mathrm{pre}}$",
        fontsize=6,
    )

    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    fig.text(
        0.06, 0.94,
        rf"{patient}  |  {band_tex}  |  subtree-Jaccard trace coloring  |  "
        rf"$\tau = 1/\lambda_{{\max}}$  |  $N = {n}$  |  sym vlim $\pm{vmax:.2f}$",
        fontweight="bold",
    )
    fig.text(
        0.06, 0.91,
        r"red U-shape = clade preserved in the partner phase but not in "
        r"rPre (no self-match in any formula).  rPre = baseline (gray).  "
        "Bottom-of-leg stubs: probe (shaft) identity.",
        fontsize=7, color="#444",
    )

    pdf.savefig(fig)
    plt.close(fig)


def build_band_pdf(band: str) -> Path:
    out_pdf = OUT_DIR / f"subtree_jaccard_trace_{band}_all_patients.pdf"
    print(f"\n=== {band} → {out_pdf.name} ===")
    with PdfPages(out_pdf) as pdf:
        for pat in COHORT:
            render_patient_page(pat, band, pdf)
    print(f"DONE: {out_pdf}")
    return out_pdf


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--band", default=None,
                        help="single band (default: all 6 bands)")
    args = parser.parse_args()
    bands = [args.band] if args.band else list(BRAIN_BANDS.keys())
    for band in bands:
        build_band_pdf(band)


if __name__ == "__main__":
    main()
