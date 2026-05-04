"""Build the d_P vs d_S example-pairs figure.

Two example (patient, band, phase_A, phase_B) cells are shown:

1. **High d_S / low d_P**: Pat_14, low_gamma, task_learn ↔ rest_post.
   Rank ordering of edges shuffled (d_S = 0.67) but the linear fit on
   raw values is still tight (d_P = 0.08). Topological reorganization
   without amplitude redistribution.

2. **Low d_S / high d_P**: Pat_02, high_gamma, task_test ↔ rest_post.
   Linear fit poor (d_P = 0.26) but ranks largely preserved
   (d_S = 0.13). Volume / amplitude redistribution without topological
   reorganization.

For each example, four panels: adjacency heatmap A, adjacency heatmap
B, scatter of `triu(B)` vs `triu(A)` (with Pearson + Spearman lines),
network drawing of (A, B) overlaid on a circular layout with edges
above the 90th percentile.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import networkx as nx
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import pearsonr, spearmanr

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

from lrg_eegfc.config.paths import DATA_ROOT
from lrg_eegfc.workflow.fc import load_fc_matrix
from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT

OUT_BASE = DATA_ROOT / "audit" / "raw_fc_phase_distance"

# Two examples, picked from the d_S vs d_P z-score disagreement table
EXAMPLES = [
    {
        "label": "high $d_S$, low $d_P$  →  topological reorganization",
        "patient": "Pat_14",
        "band": "low_gamma",
        "phase_A": "task_learn",
        "phase_B": "rest_post",
    },
    {
        "label": "low $d_S$, high $d_P$  →  volume / amplitude redistribution",
        "patient": "Pat_02",
        "band": "high_gamma",
        "phase_A": "task_test",
        "phase_B": "rest_post",
    },
]


def _common_giant_indices(adj_phases):
    sets = []
    for adj in adj_phases.values():
        G = nx.from_numpy_array(np.abs(adj))
        comps = list(nx.connected_components(G))
        sets.append(set(max(comps, key=len)) if comps else set())
    common = set.intersection(*sets) if sets else set()
    return np.array(sorted(common), dtype=int)


def _triu(A):
    return A[np.triu_indices_from(A, k=1)]


def _network_panel(ax, A, title, q=0.90):
    """Draw an undirected graph with edges above quantile q on a circular layout."""
    N = A.shape[0]
    ut = _triu(A)
    thr = float(np.quantile(ut, q))
    G = nx.Graph()
    G.add_nodes_from(range(N))
    edges = []
    weights = []
    for i in range(N):
        for j in range(i + 1, N):
            if A[i, j] >= thr:
                G.add_edge(i, j, weight=float(A[i, j]))
                edges.append((i, j))
                weights.append(float(A[i, j]))
    pos = nx.circular_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=12, node_color="#222",
                           edgecolors="white", linewidths=0.3, ax=ax)
    if edges:
        weights = np.asarray(weights)
        # Map weight to alpha and width
        wnorm = (weights - weights.min()) / (weights.max() - weights.min() + 1e-9)
        for (u, v), w in zip(edges, wnorm):
            ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                    color="tab:red", lw=0.3 + 1.4 * w, alpha=0.2 + 0.6 * w)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=8)


pdf_path = OUT_BASE / "dP_dS_example_pairs.pdf"
with PdfPages(pdf_path) as pdf:
    fig, axes = plt.subplots(2, 5, figsize=(15, 6.4),
                              gridspec_kw={"width_ratios": [1, 1, 1.05, 1, 1]})
    for r, ex in enumerate(EXAMPLES):
        # Load both phase matrices, restrict to V*
        adj_phases = {}
        for ph in (ex["phase_A"], ex["phase_B"]):
            A = load_fc_matrix(ex["patient"], ph, ex["band"], "imcoh_abs")
            adj_phases[ph] = np.asarray(A, dtype=np.float64)
        common_idx = _common_giant_indices(adj_phases)
        A_a = adj_phases[ex["phase_A"]][np.ix_(common_idx, common_idx)]
        A_b = adj_phases[ex["phase_B"]][np.ix_(common_idx, common_idx)]
        ut_a = _triu(A_a); ut_b = _triu(A_b)

        d_P = float(1.0 - pearsonr(ut_a, ut_b).statistic)
        d_S = float(1.0 - spearmanr(ut_a, ut_b).statistic)
        d_F = float(np.linalg.norm(A_a - A_b, ord="fro") /
                    np.sqrt(np.linalg.norm(A_a, ord="fro")
                            * np.linalg.norm(A_b, ord="fro")))

        vmax = max(A_a.max(), A_b.max())
        # Panel 1: heatmap A
        ax = axes[r, 0]
        ax.imshow(A_a, cmap="magma", vmin=0, vmax=vmax, aspect="equal")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"$A_{{{ex['phase_A']}}}$", fontsize=9)

        # Panel 2: heatmap B
        ax = axes[r, 1]
        ax.imshow(A_b, cmap="magma", vmin=0, vmax=vmax, aspect="equal")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"$A_{{{ex['phase_B']}}}$", fontsize=9)

        # Panel 3: scatter triu(B) vs triu(A)
        ax = axes[r, 2]
        ax.scatter(ut_a, ut_b, s=2, alpha=0.35, color="#3a3a3a")
        # Pearson best-fit line (linear regression)
        slope, intercept = np.polyfit(ut_a, ut_b, 1)
        xline = np.linspace(ut_a.min(), ut_a.max(), 100)
        ax.plot(xline, intercept + slope * xline, color="tab:blue", lw=1.4,
                label="Pearson fit")
        ax.plot([ut_a.min(), ut_a.max()], [ut_a.min(), ut_a.max()],
                color="k", lw=0.6, ls="--", alpha=0.6, label="identity")
        ax.set_xlabel(f"triu $A_{{{ex['phase_A'][:4]}}}$", fontsize=8)
        ax.set_ylabel(f"triu $A_{{{ex['phase_B'][:4]}}}$", fontsize=8)
        ax.tick_params(labelsize=6)
        ax.legend(fontsize=6, loc="upper left")
        ax.text(0.05, 0.65,
                f"$d_P = {d_P:.3f}$\n$d_S = {d_S:.3f}$\n$d_F = {d_F:.3f}$",
                transform=ax.transAxes, fontsize=8, va="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="gray", alpha=0.85))

        # Panel 4: network A
        _network_panel(axes[r, 3], A_a,
                       f"network $A_{{{ex['phase_A']}}}$ (top 10%)")
        # Panel 5: network B
        _network_panel(axes[r, 4], A_b,
                       f"network $A_{{{ex['phase_B']}}}$ (top 10%)")

        # Row title — left of all panels
        axes[r, 0].annotate(
            f"{ex['label']}\n{ex['patient']} · {BRAIN_BAND_TEX_DICT[ex['band']]} · "
            f"{ex['phase_A']} ↔ {ex['phase_B']}",
            xy=(-0.35, 0.5), xycoords="axes fraction",
            fontsize=9, ha="right", va="center", rotation=90,
            wrap=True,
        )

    fig.text(0.99, 0.01, "dP_dS_example_pairs — d_P (volume) vs d_S (topology)",
             ha="right", fontsize=6, color="gray")
    fig.tight_layout()
    pdf.savefig(fig, dpi=200); plt.close(fig)
print(f"wrote {pdf_path}")
