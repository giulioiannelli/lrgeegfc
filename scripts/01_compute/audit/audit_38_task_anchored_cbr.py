#!/usr/bin/env python3
"""Audit 38 — Task-anchored Cohesion CBR.

Per scope: ``.agents/guides/task-persistence-investigation/2026-04-25_task-anchored-cbr.md``
Per Section-5 plan §3.1: localization probe complementary to the existing
rest_post-anchored ``audit_12_cohesion_cbr.py``.

Enumerates internal subtrees S of the task_test dendrogram (T^TT). For each
S, computes:

    COH_RPre  = robust_cohesion(S, T^RPre)
    COH_RPost = robust_cohesion(S, T^RPost)

Soft 4-corner affinity in the (J_RPre, J_RPost) unit-square:

    TRACE corner    = (0, 1)  -- task-induced (low RPre) AND persists into RPost
    PERSIST corner  = (1, 1)  -- subtree always existed (anchor)
    RESET corner    = (1, 0)  -- existed in RPre, in TT, but gone in RPost
    REARRANGE corner= (0, 0)  -- subtree only exists in TT, gone everywhere else

A subtree is classified TRACE iff `argmax_p a_p = trace`. The dominant
classification + confidence flow into the per-leaf projection following the
audit_12 smallest-containing-anchor rule.

Outputs
-------
``data/audit/per_patient_hierarchy_task_anchored/{Pat_XX}/{band}_hierarchy-crossphase.pdf``
``data/audit/per_patient_hierarchy_task_anchored/leaf_assignment.csv``
``data/audit/per_patient_hierarchy_task_anchored/anchor_classification.csv``
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

from lrg_eegfc.utils.scripting import setup_script_env

ROOT = setup_script_env()

# Reuse audit_12 helpers + audit_08 layout primitives
sys.path.insert(0, str(ROOT / "scripts" / "01_compute" / "audit"))
from audit_08_per_patient_hierarchy import (  # noqa: E402
    BAND_ORDER, PHASE_ORDER, PHASE_TEX, BRAIN_BAND_TEX_DICT,
    PATIENTS_4PHASE, PATTERN_COLOR, BACKGROUND_GRAY,
    SIZE_MIN, SIZE_MAX,
    nodes_with_leaves, _load_Z,
)
from audit_12_cohesion_cbr import (  # noqa: E402
    robust_cohesion, best_jaccard_local, corner_affinities, project_to_leaves,
    link_colors_for_phase,
    PATTERNS, ALPHA_CLASSIFIED, ALPHA_UNCLASSIFIED, NEUTRAL_GRAY,
)

OUT_DIR = ROOT / "data" / "audit" / "per_patient_hierarchy_task_anchored"


def per_anchor_task(Zs: dict[str, np.ndarray]) -> tuple[list[dict], dict[str, list[dict]]]:
    """Enumerate T^TT subtrees; classify in (COH_rpre, COH_rpost) plane."""
    nodes = {ph: nodes_with_leaves(Z) if Z is not None else None for ph, Z in Zs.items()}
    if any(nodes[ph] is None for ph in ("task_test", "rest_pre", "rest_post")):
        return [], nodes

    Z_tt = Zs["task_test"]
    dmax_tt = float(Z_tt[-1, 2])
    rows = []
    for η in nodes["task_test"]:
        if η["size"] < SIZE_MIN or η["size"] > SIZE_MAX:
            continue
        S = η["leaves"]
        COH_rpre = robust_cohesion(S, nodes["rest_pre"])
        COH_rpost = robust_cohesion(S, nodes["rest_post"])
        # TRACE corner = (low rpre, high rpost) on the (COH_rpre, COH_rpost) plane.
        aff = corner_affinities(COH_rpre, COH_rpost)
        # Jaccard sidecar
        J_rpre = best_jaccard_local(S, nodes["rest_pre"])
        J_rpost = best_jaccard_local(S, nodes["rest_post"])
        rows.append({
            "anchor_row": int(η["row"]),
            "size": int(η["size"]),
            "h_rel": float(η["h"]) / dmax_tt if dmax_tt > 0 else 0.0,
            "COH_rpre": COH_rpre, "COH_rpost": COH_rpost,
            "J_rpre": J_rpre, "J_rpost": J_rpost,
            "a_trace": aff["trace"], "a_persist": aff["persist"],
            "a_reset": aff["reset"], "a_rearrange": aff["rearrange"],
            "leaves": S,
        })
    return rows, nodes


def plot_panel(ax, ax_strip, Z, dominant, alpha, title=""):
    if Z is None:
        ax.text(0.5, 0.5, "missing", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        ax_strip.set_axis_off()
        return
    n_leaves = Z.shape[0] + 1
    link_colors = link_colors_for_phase(Z, dominant, alpha)
    dendrogram(
        Z, ax=ax, no_labels=True, color_threshold=0.0,
        link_color_func=lambda nid: link_colors.get(nid, BACKGROUND_GRAY),
        above_threshold_color=BACKGROUND_GRAY,
    )
    ax.set_title(title)
    leaf_order = [int(t.get_text()) for t in ax.get_xticklabels() if t.get_text()]
    if not leaf_order:
        leaf_order = list(range(n_leaves))
    leaf_order = leaf_order[:n_leaves]
    # Colour codebar
    bar = np.zeros((n_leaves, 4))
    for i, ℓ in enumerate(leaf_order):
        bar[i] = _leaf_rgba(dominant[ℓ], alpha[ℓ])
    ax_strip.imshow(bar.reshape(1, n_leaves, 4), aspect="auto", interpolation="nearest")
    ax_strip.set_yticks([])
    ax_strip.set_xticks([])


def _leaf_rgba(label, alpha):
    if label == "unclassified":
        rgba = list(to_rgba(NEUTRAL_GRAY))
    else:
        rgba = list(to_rgba(PATTERN_COLOR[label]))
    rgba[3] = alpha
    return tuple(rgba)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    leaf_rows, anchor_rows_all = [], []
    for pat in PATIENTS_4PHASE:
        for band in BAND_ORDER:
            Zs = {ph: _load_Z(pat, ph, band) for ph in PHASE_ORDER}
            if Zs["task_test"] is None:
                continue
            rows, nodes = per_anchor_task(Zs)
            n_leaves = Zs["task_test"].shape[0] + 1
            dom, alp, A_tab, classified = project_to_leaves(n_leaves, rows)
            for r in rows:
                anchor_rows_all.append({
                    k: v if k != "leaves" else "|".join(str(x) for x in sorted(v))
                    for k, v in {**r, "patient": pat, "band": band}.items()
                })
            for ℓ in range(n_leaves):
                leaf_rows.append({
                    "patient": pat, "band": band, "leaf_id": ℓ,
                    "dominant": dom[ℓ], "confidence": float(alp[ℓ]),
                    "A_trace": float(A_tab[ℓ, 0]), "A_persist": float(A_tab[ℓ, 1]),
                    "A_reset": float(A_tab[ℓ, 2]), "A_rearrange": float(A_tab[ℓ, 3]),
                    "classified": bool(classified[ℓ]),
                })
            # Per-patient-band plot
            (OUT_DIR / pat).mkdir(parents=True, exist_ok=True)
            fig, axes = plt.subplots(2, 4, figsize=(16, 5),
                                     gridspec_kw={"height_ratios": [4, 0.4], "hspace": 0.05})
            for j, ph in enumerate(PHASE_ORDER):
                plot_panel(axes[0, j], axes[1, j], Zs[ph], dom, alp,
                           title=PHASE_TEX.get(ph, ph))
            fig.tight_layout()
            fig.savefig(OUT_DIR / pat / f"{band}_hierarchy-crossphase.pdf")
            plt.close(fig)

    pd.DataFrame(leaf_rows).to_csv(OUT_DIR / "leaf_assignment.csv", index=False)
    pd.DataFrame(anchor_rows_all).to_csv(OUT_DIR / "anchor_classification.csv", index=False)

    df_leaf = pd.DataFrame(leaf_rows)
    summary = (
        df_leaf.assign(is_trace_leaf=(df_leaf["dominant"] == "trace").astype(int))
        .groupby(["patient", "band"])["is_trace_leaf"]
        .agg(["count", "sum"])
        .reset_index()
        .rename(columns={"count": "n_leaves", "sum": "n_trace_leaves"})
    )
    summary["frac_trace_leaves"] = summary["n_trace_leaves"] / summary["n_leaves"]
    summary.to_csv(OUT_DIR / "per_patient_summary.csv", index=False)

    cohort = (
        summary.groupby("band")["frac_trace_leaves"]
        .agg(["mean", "median", "min", "max", "count"])
        .reset_index()
    )
    cohort["n_pats_ge_10pct"] = summary.groupby("band").apply(
        lambda g: int((g["frac_trace_leaves"] >= 0.10).sum())
    ).values
    cohort["n_pats_ge_20pct"] = summary.groupby("band").apply(
        lambda g: int((g["frac_trace_leaves"] >= 0.20).sum())
    ).values
    cohort.to_csv(OUT_DIR / "cohort_summary.csv", index=False)
    print(cohort.to_string(index=False))
    print(f"[audit_38] outputs at {OUT_DIR}")


if __name__ == "__main__":
    main()
