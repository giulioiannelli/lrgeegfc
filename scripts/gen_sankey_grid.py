#!/usr/bin/env python3
"""Generate Sankey community evolution diagrams for all phase/band combos.

Produces interactive Plotly HTML + static PDF for each (phase, band) pair.
Uses fixed (τ, n) pairs for the alluvial visualization (communities must
decrease monotonically for a clean Sankey layout).

Run: python scripts/gen_sankey_grid.py
"""
import sys
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import networkx as nx
import plotly.graph_objects as go
from plotly.express import colors as px_colors

from lrg_eegfc.config.const import BRAIN_BANDS_NAMES, BRAIN_BAND_TEX_DICT, PHASE_LABELS
from lrg_eegfc.workflow.msc import load_msc_matrix
from lrg_eegfc.visuals.lrg import _load_channel_labels
from lrg_eegfc.visuals.metastable import compute_clustering_across_tau

# ── Config ────────────────────────────────────────────────────────────
PATIENT = "Pat_02"
MSC_CACHE = ROOT / "data" / "msc_cache"
DATASET_ROOT = ROOT / "data" / "stereoeeg_patients"
OUTPUT_DIR = ROOT / "data" / "figures" / "report_mslcd_section" / PATIENT / "sankey"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

COMBOS = [(phase, band) for phase in PHASE_LABELS for band in BRAIN_BANDS_NAMES]

TAU_NCLUST_PAIRS = [
    (0.1, 15), (0.3, 12), (0.5, 8), (0.8, 6),
    (1.0, 5), (2.0, 3), (5.0, 2),
]
ALLUVIAL_TAUS = [t for t, _ in TAU_NCLUST_PAIRS]
ALLUVIAL_NCLUST = [n for _, n in TAU_NCLUST_PAIRS]

COLOR_PALETTE = (px_colors.qualitative.Set3 + px_colors.qualitative.Dark2
                 + px_colors.qualitative.Pastel)

# ── Load channel labels (shared) ─────────────────────────────────────
channel_labels = _load_channel_labels(PATIENT, DATASET_ROOT)
N = len(channel_labels)
print(f"{PATIENT}: {N} channels, {len(COMBOS)} combos to generate\n")


def make_sankey(phase, band, partitions):
    """Build and save Plotly Sankey for one (phase, band) combo."""
    band_tex = BRAIN_BAND_TEX_DICT[band]

    sankey_labels = []
    sankey_hovers = []
    sankey_colors = []
    tau_to_offset = {}
    node_count = 0

    for i, tau in enumerate(ALLUVIAL_TAUS):
        tau_to_offset[tau] = node_count
        clusters_arr = partitions[tau]
        unique_c = np.unique(clusters_arr)
        for c in unique_c:
            members = sorted([channel_labels[j] for j in range(N)
                              if clusters_arr[j] == c])
            n_mem = len(members)
            sankey_labels.append(f"C{c} ({n_mem})")
            hover_pins = "<br>".join(members) if n_mem <= 20 else (
                "<br>".join(members[:20]) + f"<br>... +{n_mem - 20} more")
            sankey_hovers.append(
                f"τ={tau:.2f} C{c}: {n_mem} pins<br>{hover_pins}")
            sankey_colors.append(
                COLOR_PALETTE[(c - 1) % len(COLOR_PALETTE)])
        node_count += len(unique_c)

    sources, targets, values, link_colors = [], [], [], []
    for i in range(len(ALLUVIAL_TAUS) - 1):
        tau1, tau2 = ALLUVIAL_TAUS[i], ALLUVIAL_TAUS[i + 1]
        c1_arr, c2_arr = partitions[tau1], partitions[tau2]
        trans = Counter(zip(c1_arr, c2_arr))
        uc1 = list(np.unique(c1_arr))
        uc2 = list(np.unique(c2_arr))
        for (cl, cr), cnt in trans.items():
            sources.append(tau_to_offset[tau1] + uc1.index(cl))
            targets.append(tau_to_offset[tau2] + uc2.index(cr))
            values.append(cnt)
            tc = sankey_colors[tau_to_offset[tau2] + uc2.index(cr)]
            if tc.startswith("rgb"):
                link_colors.append(
                    tc.replace("rgb", "rgba").replace(")", ",0.4)"))
            elif tc.startswith("#"):
                hx = tc.lstrip("#")
                r, g, b = (int(hx[j:j+2], 16) for j in (0, 2, 4))
                link_colors.append(f"rgba({r},{g},{b},0.4)")
            else:
                link_colors.append("rgba(180,180,180,0.4)")

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=20, thickness=30,
            line=dict(color="black", width=0.5),
            label=sankey_labels, color=sankey_colors,
            hovertemplate="%{customdata}<extra></extra>",
            customdata=sankey_hovers,
        ),
        link=dict(source=sources, target=targets,
                  value=values, color=link_colors),
    )])
    fig.update_layout(
        title_text=f"{PATIENT} {phase} {band_tex} — Cluster Evolution",
        font_size=10, width=1400, height=700,
    )

    out_html = OUTPUT_DIR / f"fig_sankey_{phase}_{band}.html"
    fig.write_html(str(out_html))

    out_pdf = OUTPUT_DIR / f"fig_sankey_{phase}_{band}.pdf"
    try:
        fig.write_image(str(out_pdf), format="pdf",
                        width=1400, height=700, scale=2)
    except Exception:
        try:
            out_png = out_pdf.with_suffix(".png")
            fig.write_image(str(out_png), format="png",
                            width=1400, height=700, scale=2)
        except Exception:
            pass

    return out_html.name


# ── Generate ──────────────────────────────────────────────────────────
for phase, band in COMBOS:
    print(f"  {phase}/{band}...", end=" ")

    A = load_msc_matrix(PATIENT, phase, band,
                        cache_root=MSC_CACHE, sparsify="none",
                        n_surrogates=0, nperseg=4096)
    if A is None:
        print("SKIP (no MSC cache)")
        continue
    np.fill_diagonal(A, 0)
    assert A.shape[0] == N

    G = nx.from_numpy_array(A)
    Gcc_nodes = max(nx.connected_components(G), key=len)
    Gcc = G.subgraph(Gcc_nodes).copy()

    partitions, _ = compute_clustering_across_tau(
        None, np.array(ALLUVIAL_TAUS), Gcc,
        n_clusters_list=ALLUVIAL_NCLUST,
    )

    fname = make_sankey(phase, band, partitions)
    print(f"saved: {fname}")

print(f"\nDone! Figures in {OUTPUT_DIR}")
