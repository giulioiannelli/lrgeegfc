#!/usr/bin/env python3
"""Generate 3D brain connectome HTML plots at multiple community scales.

For each patient, phase, and band, generates interactive Plotly 3D brain
figures with nodes colored by community partition at scales:
  n = 3, 5, 8, 10, 15, 20, 30, and n* (optimal from LRG).

Output:
  data/figures/brain_connectome_multiscale/{patient}/
      fig_brain_connectome_n{N}_{phase}_{band}.html

Run:
  python scripts/gen_brain_connectome_multiscale.py
"""
import re
import warnings
from pathlib import Path

from lrg_eegfc.utils.scripting import setup_script_env
ROOT = setup_script_env()

import matplotlib.colors
import networkx as nx
import numpy as np
import plotly.graph_objects as go
from scipy.cluster.hierarchy import fcluster
from scipy.spatial import ConvexHull

from lrg_eegfc.config.const import BRAIN_BAND_TEX_DICT, nperseg_for_fs
from lrg_eegfc.visuals.lrg import _load_channel_labels
from lrg_eegfc.visuals.spatial import load_spatial_metadata, prepare_spatial_coordinates
from lrg_eegfc.workflow.lrg import load_lrg_result
from lrg_eegfc.workflow.msc import load_msc_matrix

# ── Configuration ────────────────────────────────────────────────────
PATIENTS = ["Pat_02", "Pat_03", "Pat_05", "Pat_07", "Pat_08"]
PHASES = ["rest_pre", "rest_post"]
BANDS = ["delta", "theta", "alpha", "beta"]
FC_METHOD = "imcoh_abs"
FIXED_N_VALUES = [3, 5, 8, 10, 15, 20, 30]

from lrg_eegfc.config.paths import IMCOH_CACHE, IMCOH_LRG_CACHE, SEEG_DATAPATH, FIGURES_ROOT

MSC_CACHE = IMCOH_CACHE
LRG_CACHE = IMCOH_LRG_CACHE
DATASET_ROOT = SEEG_DATAPATH
OUTPUT_ROOT = FIGURES_ROOT / "brain_connectome_multiscale_imcoh"

# Patient-specific nperseg values (Pat_03 is 1024 Hz)
PATIENT_NPERSEG = {
    "Pat_02": 4096,
    "Pat_03": 2048,
    "Pat_05": 4096,
    "Pat_07": 4096,
    "Pat_08": 4096,
}

# Edge display parameters
EDGE_PERCENTILE = 85   # show top 15%
ARC_NPTS = 14          # points per arc
ARC_BULGE = 6.0        # mm offset at midpoint
EDGE_POWER = 1.8
WIDTH_MIN, WIDTH_MAX = 0.5, 16.0
N_BUCKETS = 15
NODE_RADIUS = 0.9      # mm

# ── Curated distinct colors ─────────────────────────────────────────
_DISTINCT_COLORS = [
    "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
    "#d63b00", "#f032e6", "#469990", "#9A6324", "#800000",
    "#808000", "#000075", "#e03080", "#1a7820", "#7030c0",
    "#b05010", "#206060", "#c04080", "#305090", "#704020",
]


def get_distinct_colors(n):
    """Generate n maximally perceptually distinct colors."""
    if n <= len(_DISTINCT_COLORS):
        return [matplotlib.colors.to_rgba(c) for c in _DISTINCT_COLORS[:n]]
    colors = []
    for i in range(n):
        hue = (i * 0.618033988749895) % 1.0
        colors.append(matplotlib.colors.hsv_to_rgb([hue, 0.75, 0.85]))
    return colors


# ── Brain surface (fsaverage pial) — loaded once ────────────────────
print("Loading fsaverage brain meshes...")
from nilearn.datasets import fetch_surf_fsaverage
from nilearn.surface import load_surf_mesh

fsaverage = fetch_surf_fsaverage()
brain_meshes = []
for hemi_name, mesh_key in [("L", "pial_left"), ("R", "pial_right")]:
    mesh = load_surf_mesh(fsaverage[mesh_key])
    brain_meshes.append((hemi_name, mesh.coordinates, mesh.faces))
print("  Done.")


# ── Unit sphere mesh (reused per node) ──────────────────────────────
def _unit_sphere(n_lat=8, n_lon=12):
    """Create unit sphere vertices and triangle faces."""
    verts = [[0, 0, 1]]  # top pole
    for i in range(1, n_lat):
        lat = np.pi * i / n_lat
        for j in range(n_lon):
            lon = 2 * np.pi * j / n_lon
            verts.append([
                np.sin(lat) * np.cos(lon),
                np.sin(lat) * np.sin(lon),
                np.cos(lat),
            ])
    verts.append([0, 0, -1])  # bottom pole
    verts = np.array(verts)

    faces = []
    # Top cap
    for j in range(n_lon):
        faces.append([0, 1 + j, 1 + (j + 1) % n_lon])
    # Middle strips
    for i in range(n_lat - 2):
        ring = 1 + i * n_lon
        nxt = ring + n_lon
        for j in range(n_lon):
            j1 = (j + 1) % n_lon
            faces.append([ring + j, nxt + j, ring + j1])
            faces.append([ring + j1, nxt + j, nxt + j1])
    # Bottom cap
    bottom = len(verts) - 1
    last_ring = 1 + (n_lat - 2) * n_lon
    for j in range(n_lon):
        faces.append([bottom, last_ring + (j + 1) % n_lon, last_ring + j])
    return verts, np.array(faces)


_sph_v, _sph_f = _unit_sphere()
_n_sph_v = len(_sph_v)


# ── Sphere offsets for convex hull inflation ─────────────────────────
def _sphere_offsets(radius, n_pts=26):
    pts = []
    golden = (1 + np.sqrt(5)) / 2
    for i in range(n_pts):
        theta = np.arccos(1 - 2 * (i + 0.5) / n_pts)
        phi = 2 * np.pi * i / golden
        pts.append([
            radius * np.sin(theta) * np.cos(phi),
            radius * np.sin(theta) * np.sin(phi),
            radius * np.cos(theta),
        ])
    return np.array(pts)


_offsets_tube = _sphere_offsets(2.0, 26)     # per-electrode tubes
_offsets_bridge = _sphere_offsets(3.0, 26)   # cross-electrode bridges


# ── Arc interpolation ───────────────────────────────────────────────
def _make_arc(pi, pj, centroid):
    """Create a curved arc from pi to pj, bulging outward from centroid."""
    mid = (pi + pj) / 2
    outward = mid - centroid
    norm = np.linalg.norm(outward)
    if norm < 1e-6:
        outward = np.array([0, 0, 1.0])
    else:
        outward = outward / norm
    edge_dir = pj - pi
    edge_len = np.linalg.norm(edge_dir)
    if edge_len > 1e-6:
        perp = np.cross(edge_dir / edge_len, outward)
        perp_norm = np.linalg.norm(perp)
        if perp_norm > 1e-6:
            perp = perp / perp_norm
            outward = 0.7 * outward + 0.3 * perp
            outward = outward / np.linalg.norm(outward)
    control = mid + ARC_BULGE * outward
    t = np.linspace(0, 1, ARC_NPTS)
    pts = (
        np.outer((1 - t) ** 2, pi)
        + np.outer(2 * t * (1 - t), control)
        + np.outer(t ** 2, pj)
    )
    return pts


def _edge_rgba(w_n):
    gray = int(80 * (1 - w_n))
    alpha = 0.15 + 0.75 * np.power(w_n, 1.0)
    return f"rgba({gray},{gray},{gray},{alpha:.3f})"


# ── Build common traces for a given patient/phase/band ──────────────
def build_shared_traces(coords_mni, channel_labels, A, N):
    """Build electrode shaft traces and edge traces (shared across n* values)."""
    # Electrode groups
    electrode_groups = {}
    for i, lbl in enumerate(channel_labels):
        m = re.match(r"([A-Za-z']+)", lbl)
        if m:
            electrode_groups.setdefault(m.group(1), []).append(i)

    # Shaft traces
    shaft_traces = []
    for prefix, indices in electrode_groups.items():
        if len(indices) < 2:
            continue
        indices_sorted = sorted(
            indices,
            key=lambda idx: int(re.search(r"(\d+)", channel_labels[idx]).group(1))
            if re.search(r"(\d+)", channel_labels[idx])
            else 0,
        )
        x = [coords_mni[i, 0] for i in indices_sorted]
        y = [coords_mni[i, 1] for i in indices_sorted]
        z = [coords_mni[i, 2] for i in indices_sorted]
        shaft_traces.append(go.Scatter3d(
            x=x, y=y, z=z, mode="lines",
            line=dict(width=3.5, color="rgba(60,60,60,0.35)"),
            hoverinfo="none", showlegend=False,
        ))

    # Edge selection and bucketing
    triu_idx = np.triu_indices(N, k=1)
    edge_w_all = A[triu_idx]
    edge_threshold = np.percentile(edge_w_all, EDGE_PERCENTILE)
    strong_mask = edge_w_all >= edge_threshold
    edge_i = triu_idx[0][strong_mask]
    edge_j = triu_idx[1][strong_mask]
    edge_w = edge_w_all[strong_mask]

    if len(edge_w) == 0:
        return shaft_traces, [], edge_w_all.min(), edge_w_all.max()

    w_min_e, w_max_e = edge_w.min(), edge_w.max()
    w_norm = (edge_w - w_min_e) / (w_max_e - w_min_e + 1e-10)

    centroid = coords_mni.mean(axis=0)
    bucket_bounds = np.linspace(0, 1.001, N_BUCKETS + 1)

    edge_traces = []
    for b in range(N_BUCKETS):
        bmask = (w_norm >= bucket_bounds[b]) & (w_norm < bucket_bounds[b + 1])
        if not np.any(bmask):
            continue
        bw = w_norm[bmask].mean()
        bwidth = WIDTH_MIN + np.power(bw, EDGE_POWER) * (WIDTH_MAX - WIDTH_MIN)
        bcolor = _edge_rgba(bw)
        x, y, z = [], [], []
        for ei, ej in zip(edge_i[bmask], edge_j[bmask]):
            arc = _make_arc(coords_mni[ei], coords_mni[ej], centroid)
            x.extend(arc[:, 0].tolist() + [None])
            y.extend(arc[:, 1].tolist() + [None])
            z.extend(arc[:, 2].tolist() + [None])
        edge_traces.append(go.Scatter3d(
            x=x, y=y, z=z, mode="lines",
            line=dict(width=bwidth, color=bcolor),
            hoverinfo="none", showlegend=False,
        ))

    return shaft_traces, edge_traces, w_min_e, w_max_e


def build_figure_for_n_precomputed(
    n_cut, cl, coords_mni, channel_labels, N,
    shaft_traces, edge_traces, w_min_e, w_max_e,
    patient, phase, band,
):
    """Build a complete Plotly figure with pre-computed cluster labels."""
    palette = get_distinct_colors(n_cut)
    unique_cl = sorted(np.unique(cl))
    cl_cmap = {c: palette[i % len(palette)] for i, c in enumerate(unique_cl)}
    node_hex = [matplotlib.colors.to_hex(cl_cmap[c]) for c in cl]

    # Brain surface traces
    btrs = []
    for hemi_name, vtx, tri in brain_meshes:
        btrs.append(go.Mesh3d(
            x=vtx[:, 0], y=vtx[:, 1], z=vtx[:, 2],
            i=tri[:, 0], j=tri[:, 1], k=tri[:, 2],
            color="lightgray", opacity=0.18,
            hoverinfo="none", showlegend=False,
            lighting=dict(ambient=0.9, diffuse=0.1),
        ))

    # Map each node to electrode prefix
    node_electrode = {}
    for i, lbl in enumerate(channel_labels):
        m = re.match(r"([A-Za-z']+)", lbl)
        node_electrode[i] = m.group(1) if m else "?"

    # Cluster hull traces
    hull_trs = []
    for c in unique_cl:
        c_indices = np.where(cl == c)[0]
        if len(c_indices) < 2:
            continue
        elec_groups = {}
        for idx in c_indices:
            elec_groups.setdefault(node_electrode[idx], []).append(idx)

        hex_c = matplotlib.colors.to_hex(cl_cmap[c])

        # Layer 1: per-electrode tubes
        for elec, indices in elec_groups.items():
            pts = coords_mni[indices]
            if pts.shape[0] < 2:
                continue
            inflated = [pts]
            for off in _offsets_tube:
                inflated.append(pts + off[None, :])
            inflated = np.vstack(inflated)
            try:
                hull = ConvexHull(inflated)
            except Exception:
                continue
            hull_trs.append(go.Mesh3d(
                x=inflated[hull.vertices, 0],
                y=inflated[hull.vertices, 1],
                z=inflated[hull.vertices, 2],
                alphahull=0, color=hex_c, opacity=0.30,
                hoverinfo="none", showlegend=False,
                lighting=dict(ambient=0.7, diffuse=0.3, specular=0.15, fresnel=0.1),
                flatshading=False,
            ))

        # Layer 2: cross-electrode bridge
        if len(elec_groups) >= 2:
            all_pts = coords_mni[c_indices]
            inflated = [all_pts]
            for off in _offsets_bridge:
                inflated.append(all_pts + off[None, :])
            inflated = np.vstack(inflated)
            try:
                hull = ConvexHull(inflated)
            except Exception:
                pass
            else:
                hull_trs.append(go.Mesh3d(
                    x=inflated[hull.vertices, 0],
                    y=inflated[hull.vertices, 1],
                    z=inflated[hull.vertices, 2],
                    alphahull=0, color=hex_c, opacity=0.14,
                    hoverinfo="none", showlegend=False,
                    lighting=dict(ambient=0.8, diffuse=0.2, specular=0.05, fresnel=0.1),
                    flatshading=False,
                ))

    # Node spheres (Mesh3d)
    node_sphere_trs = []
    for c in unique_cl:
        c_idxs = np.where(cl == c)[0]
        all_v, all_fi, all_fj, all_fk = [], [], [], []
        for k_node, idx in enumerate(c_idxs):
            sv = _sph_v * NODE_RADIUS + coords_mni[idx]
            offset = k_node * _n_sph_v
            all_v.append(sv)
            all_fi.append(_sph_f[:, 0] + offset)
            all_fj.append(_sph_f[:, 1] + offset)
            all_fk.append(_sph_f[:, 2] + offset)
        all_v = np.vstack(all_v)
        all_fi = np.concatenate(all_fi)
        all_fj = np.concatenate(all_fj)
        all_fk = np.concatenate(all_fk)
        hex_c = matplotlib.colors.to_hex(cl_cmap[c])
        node_sphere_trs.append(go.Mesh3d(
            x=all_v[:, 0], y=all_v[:, 1], z=all_v[:, 2],
            i=all_fi, j=all_fj, k=all_fk,
            color=hex_c, opacity=1.0,
            hoverinfo="none", showlegend=False,
            lighting=dict(ambient=0.6, diffuse=0.4, specular=0.3, fresnel=0.2),
        ))

    # Node labels + hover
    LABEL_OFFSET_Z = NODE_RADIUS + 1.5
    label_coords = coords_mni.copy()
    label_coords[:, 2] += LABEL_OFFSET_Z
    node_label_tr = go.Scatter3d(
        x=label_coords[:, 0], y=label_coords[:, 1], z=label_coords[:, 2],
        mode="text",
        text=[channel_labels[i] for i in range(N)],
        textposition="middle center",
        textfont=dict(size=11, color="rgba(30,30,30,0.85)"),
        hoverinfo="text",
        hovertext=[f"{channel_labels[i]}<br>C{cl[i]}" for i in range(N)],
        showlegend=False,
    )

    # Colorbar
    cb_tr = go.Scatter3d(
        x=[None], y=[None], z=[None], mode="markers",
        marker=dict(
            size=0.1,
            colorscale=[[0, "rgba(100,100,100,0.08)"], [1, "rgba(0,0,0,0.85)"]],
            cmin=w_min_e, cmax=w_max_e,
            colorbar=dict(
                title="MSC", thickness=20, len=0.6, x=0.98,
                titlefont=dict(size=16), tickfont=dict(size=13),
            ),
            color=[w_min_e],
        ),
        hoverinfo="none", showlegend=False,
    )

    # Assemble figure
    band_tex = BRAIN_BAND_TEX_DICT.get(band, band)
    title_text = f"{patient} | {phase} | {band_tex} | n={n_cut} communities"

    fig = go.Figure(
        data=(
            btrs + hull_trs + edge_traces + shaft_traces
            + [cb_tr] + node_sphere_trs + [node_label_tr]
        )
    )
    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            aspectmode="data",
            bgcolor="white",
            camera=dict(
                eye=dict(x=-1.8, y=0.0, z=0.3),
                up=dict(x=0, y=0, z=1),
            ),
        ),
        showlegend=False,
        width=1400, height=900,
        paper_bgcolor="white",
        margin=dict(l=0, r=0, t=40, b=0),
        title=dict(text=title_text, x=0.5, font=dict(size=18)),
    )

    return fig


# ══════════════════════════════════════════════════════════════════════
# Main generation loop
# ══════════════════════════════════════════════════════════════════════
total_generated = 0
total_skipped = 0
total_errors = 0

for patient in PATIENTS:
    nperseg = PATIENT_NPERSEG[patient]
    print(f"\n{'='*60}")
    print(f"Patient: {patient} (nperseg={nperseg})")
    print(f"{'='*60}")

    # Load channel labels once per patient
    try:
        channel_labels = _load_channel_labels(patient, DATASET_ROOT)
    except Exception as e:
        print(f"  ERROR loading channel labels: {e}")
        total_errors += 1
        continue

    # Load spatial metadata once per patient
    try:
        metadata = load_spatial_metadata(patient, DATASET_ROOT)
        # Check for NaN coordinates and track which channels are valid
        valid_coord_mask = metadata[["x", "y", "z"]].notna().all(axis=1).values
        n_missing = (~valid_coord_mask).sum()
        if n_missing > 0:
            print(f"  WARNING: {n_missing} channels have missing coordinates, "
                  f"will be excluded from visualization")
        coords_mni = prepare_spatial_coordinates(
            metadata[valid_coord_mask].reset_index(drop=True),
            scale="mm", to_mni=True,
        )
    except Exception as e:
        print(f"  ERROR loading spatial coordinates: {e}")
        total_errors += 1
        continue

    for phase in PHASES:
        for band in BANDS:
            print(f"\n  {patient} / {phase} / {band}")

            # Load LRG result
            lrg = load_lrg_result(patient, phase, band, FC_METHOD, LRG_CACHE)
            if lrg is None:
                print(f"    SKIP: no LRG result cached")
                total_skipped += 1
                continue

            # Load FC matrix directly (ImCoh)
            fc_path = MSC_CACHE / patient / f"{band}_{phase}_imcoh_sparsify-none_nperseg-{nperseg}.npy"
            if not fc_path.exists():
                print(f"    SKIP: no FC matrix at {fc_path}")
                total_skipped += 1
                continue
            A = np.load(fc_path)
            if A is None:
                total_skipped += 1
                continue

            np.fill_diagonal(A, 0)
            N_full = A.shape[0]

            # Verify base dimension matches labels
            if N_full != len(channel_labels):
                print(f"    SKIP: dimension mismatch N={N_full} vs labels={len(channel_labels)}")
                total_skipped += 1
                continue

            # If some channels lack coordinates, filter the matrix and labels
            if n_missing > 0:
                A_viz = A[np.ix_(valid_coord_mask, valid_coord_mask)]
                channel_labels_eff = [channel_labels[i] for i in range(N_full) if valid_coord_mask[i]]
            else:
                A_viz = A
                channel_labels_eff = channel_labels

            N = A_viz.shape[0]

            if N != coords_mni.shape[0]:
                print(f"    SKIP: dimension mismatch N={N} vs coords={coords_mni.shape[0]}")
                total_skipped += 1
                continue

            # The LRG linkage matrix was computed on the giant component of the
            # FULL matrix, so fcluster returns labels for lrg.n_nodes nodes.
            # We need to map those back to the full index, then filter to valid coords.
            # For most patients giant_component == N_full, but handle the general case.
            lrg_n = lrg.n_nodes

            # Compute n* (optimal) from LRG
            cl_optimal = fcluster(
                lrg.linkage_matrix,
                t=lrg.optimal_threshold,
                criterion="distance",
            )
            n_star = len(np.unique(cl_optimal))

            # Build giant component node mapping (same order as LRG computation)
            G_full = nx.from_numpy_array(A)
            gcc_nodes = sorted(max(nx.connected_components(G_full), key=len))
            gcc_set = set(gcc_nodes)

            # Verify LRG node count matches giant component
            if lrg_n != len(gcc_nodes):
                print(f"    SKIP: LRG n_nodes={lrg_n} != giant component={len(gcc_nodes)}")
                total_skipped += 1
                continue

            # Build complete n-value list (fixed + n*)
            n_values = sorted(set(FIXED_N_VALUES + [n_star]))
            # Filter out n values larger than lrg_n-1 (can't have more clusters than nodes)
            n_values = [n for n in n_values if 2 <= n <= lrg_n - 1]

            print(f"    N={N} (full={N_full}, giant={lrg_n}), n*={n_star}, "
                  f"generating n={n_values}")

            # Build shared traces (edges + shafts) once per band/phase/patient
            try:
                shaft_traces, edge_traces, w_min_e, w_max_e = build_shared_traces(
                    coords_mni, channel_labels_eff, A_viz, N,
                )
            except Exception as e:
                print(f"    ERROR building shared traces: {e}")
                total_errors += 1
                continue

            # Generate figure for each n value
            output_dir = OUTPUT_ROOT / patient
            output_dir.mkdir(parents=True, exist_ok=True)

            for n_cut in n_values:
                out_path = output_dir / f"fig_brain_connectome_n{n_cut}_{phase}_{band}.html"
                if out_path.exists():
                    print(f"    n={n_cut}: already exists, skipping")
                    total_skipped += 1
                    continue

                try:
                    # fcluster on the LRG linkage gives labels for giant component nodes
                    cl_full = fcluster(lrg.linkage_matrix, n_cut, criterion="maxclust")
                    # Map back to full matrix indices, then filter to valid coords
                    cl_all = np.zeros(N_full, dtype=int)
                    for gi, ni in enumerate(gcc_nodes):
                        cl_all[ni] = cl_full[gi]
                    # Assign non-giant-component nodes to cluster 0
                    # (they won't appear if we filter by valid_coord_mask)
                    if n_missing > 0:
                        cl_viz = cl_all[valid_coord_mask]
                    else:
                        cl_viz = cl_all

                    fig = build_figure_for_n_precomputed(
                        n_cut, cl_viz, coords_mni, channel_labels_eff, N,
                        shaft_traces, edge_traces, w_min_e, w_max_e,
                        patient, phase, band,
                    )
                    fig.write_html(str(out_path))
                    total_generated += 1
                    suffix = " [n*]" if n_cut == n_star else ""
                    print(f"    n={n_cut}: saved{suffix}")
                except Exception as e:
                    print(f"    n={n_cut}: ERROR - {e}")
                    import traceback
                    traceback.print_exc()
                    total_errors += 1

print(f"\n{'='*60}")
print(f"SUMMARY")
print(f"  Generated: {total_generated}")
print(f"  Skipped:   {total_skipped}")
print(f"  Errors:    {total_errors}")
print(f"  Output:    {OUTPUT_ROOT}")
print(f"{'='*60}")
