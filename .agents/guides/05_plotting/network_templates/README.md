---
name: 05_plotting / network_templates / README
type: figure_template_index
era: IMCOH_ABS × COHORT_N10
status: round_1_iterating
created: 2026-05-05
updated: 2026-05-05
pointers:
  - ../README.md
  - ../fc_templates/README.md
  - ../colormaps-and-styles.md
  - src/lrg_eegfc/visuals/correlation.py
  - src/lrg_eegfc/visuals/msc.py
  - src/lrg_eegfc/visuals/spatial.py
  - src/lrg_eegfc/visuals/lrg.py
  - src/lrg_eegfc/visuals/plotting.py
  - scripts/archive/2026-02_imcoh-dev-notes/fig_helper_layouts.py
  - scripts/archive/2026-02_imcoh-dev-notes/fig_helper_lrg_layout.py
  - scripts/archive/2026-02_imcoh-dev-notes/fig_helper_lrg_sfdp.py
  - scripts/archive/2026-02_imcoh-dev-notes/fig_helper_rank_transform.py
  - scripts/archive/2026-02_imcoh-dev-notes/fig_helper_weight_transform.py
  - scripts/archive/2026-02_imcoh-dev-notes/fig_helper_edge_style.py
  - scripts/archive/2026-02_imcoh-dev-notes/fig_helper_backbone_layout.py
  - scripts/03_analysis/compare_network_variants.py
---

# Network-drawing templates

Library + per-template style sheets + thin demo scripts for FC network
drawings.  The network analogue of `fc_templates/`.  Library lives at
`src/lrg_eegfc/visuals/network_templates.py`; templates dispatched
here.

## Templates in this folder

| Template | Layout | Status | Script |
|---|---|---|---|
| `single_network` | one patient × one band × one phase, one panel | round_1_iterating | [`single_network.py`](single_network.py) |
| `row_per_phase` | one patient × one band × 4 phases | round_1_iterating | [`row_per_phase.py`](row_per_phase.py) |
| `row_per_band` | one patient × one phase × 6 bands | round_1_iterating | [`row_per_band.py`](row_per_band.py) |
| `mosaic_band_phase` | 4 phases × 6 bands (one patient) | round_1_iterating | [`mosaic_band_phase.py`](mosaic_band_phase.py) |
| `mosaic_layout_compare` | every layout side-by-side, one case | round_1_iterating | [`mosaic_layout_compare.py`](mosaic_layout_compare.py) |
| `lrg_seeded` | LRG-driven layout (gt SFDP / nx KK) | round_1_iterating | [`lrg_seeded.py`](lrg_seeded.py) |
| `circular_dendrogram_network` | gt curvy chord + full circular-dendrogram overlay | round_1_iterating | [`circular_dendrogram_network.py`](circular_dendrogram_network.py) |
| `chord_highlighted_edges` | gt curvy chord with highlighted edge subset on a faint background | round_1_iterating | [`chord_highlighted_edges.py`](chord_highlighted_edges.py) |
| `matrix_plus_network` | FC matrix imshow + network with cmap-linked edges; phases as rows | round_1_iterating | [`matrix_plus_network.py`](matrix_plus_network.py) |
| `metastable_sankey` | Plotly interactive Sankey of LRG cluster evolution across τ | round_1_iterating | [`metastable_sankey.py`](metastable_sankey.py) |
| `metastable_sankey_mpl` | matplotlib static-PDF alluvial twin of `metastable_sankey` | round_1_iterating | [`metastable_sankey_mpl.py`](metastable_sankey_mpl.py) |

Generated PDFs go to `data/outputs/figures/network_templates/<template_name>/`.

## Hierarchy-bundled chord family (extracted 2026-05-26 from preprint_07)

The last two templates share a single render pipeline that crystallises
the two visual primitives in panel (c) of
``preprint_07_beta_rho_split_figure_test2.py``:

  * `circular_dendrogram_network` — the "backbone" pass: every edge
    drawn under the canonical same-probe / cross-probe rule, with the
    full LRG circular dendrogram (every merge in
    ``lrg.linkage_matrix``) overlaid as radial + arc polylines.  Use
    when the LRG hierarchy itself is the figure's subject.
  * `chord_highlighted_edges` — the "highlight" pass: a faint grey
    "null backbone" of every pair, with user-named highlight classes
    (colour + per-pair magnitude + width/alpha range) layered on top.
    Use when a small subset of edges is the figure's subject.

Both call into the same library helpers
(`build_chord_depth2_layout`, `render_hierarchy_chord`,
`draw_circular_dendrogram_overlay`, `draw_radial_leaf_labels`) at
`src/lrg_eegfc/visuals/network_templates.py`.  Picking the right
template is a choice of subject: hierarchy structure (Template 1) or
selected edges (Template 2).  Set ``show_dendrogram=True`` on
Template 2 for a Template-1+2 hybrid when you genuinely need both.

## Matrix + network family (extracted 2026-05-26 from MSC-era presentation_figures)

  * `matrix_plus_network` — matrix imshow on the left + network on the
    right with edges coloured by the same imshow `cmap + norm`.  Dark
    cells map to dark edges, bright cells to bright edges — the reader
    reads both panels in the same colour language.  Phases as rows
    (``--phases rest_pre rest_post`` reproduces the reference figure
    layout).  Modular: layout, edge γ-power recipe, coloring mode
    (`cmap` / `probe` / `shaft` / `signed` / `weight`), node colouring
    (`shaft` / `community` / `uniform`), backend (NetworkX / graph-tool
    via the layout name) are independent knobs.  Spring auto-tunes
    `k = k_base/√N` by default; pass `--spring-k-scale fixed` or
    `--k <val>` for diagnostic sweeps.

The cmap-linked edge mode is a new `coloring="cmap"` option on
`draw_gamma_edges` (alongside the existing `probe / shaft / signed /
weight` modes); pass `cmap=` and `norm=` to lock the network's edges
to any external imshow scale.

## Sankey / alluvial family (extracted 2026-05-26 from MSC-era `05_metastable_sankey.ipynb`)

  * `metastable_sankey` — interactive Plotly **HTML** Sankey of LRG
    cluster evolution across a τ-schedule.  Columns = τ values, stacked
    bars = clusters within a column, ribbons = ``c_t → c_{t+1}``
    transition counts.  Cluster counts must decrease monotonically for
    the alluvial to read as hierarchical merging.  Hovering a bar
    reveals the contact-label list and cluster size.
  * `metastable_sankey_mpl` — pure-matplotlib **PDF** alluvial twin
    using cubic-bezier ribbon polygons.  No Plotly, no ``kaleido``, no
    rasterisation; fully vector at ~20 KB for a 5-column figure.  Use
    this for paper figures; use the Plotly variant when interactive
    hover is needed.

Both backends share ``compute_sankey_flows(partitions, tau_values,
node_labels)`` which returns a ``SankeyFlowData`` dataclass (cluster
sizes per τ, sources / targets / values, per-bar colours, hover
strings).  Plotly + matplotlib renderers are thin wrappers around this
dataclass.  Library lives at
``src/lrg_eegfc/visuals/metastable.py``.  ``create_sankey_diagram`` is
preserved as a back-compat shim around the new two-step path.

## What follows is the original audit + design (kept for reference)

Skip ahead to "§13 Where this lives" if you just want the build
layout; everything above describes the lessons-learned that went into
the templates.

## TL;DR

A network drawing is **three independent knobs**: *positions* (layout),
*edge style* (width + alpha + color), *node decoration* (color + size).
Decoupling them is what makes a template reusable. Our archive already
contains the core recipe — a γ-power edge-renderer with same-probe red
highlighting, eight layout variants (NetworkX + graph-tool), and an
LRG-seeded module-aware layout — but every helper still lives as
`_private` code in `scripts/archive/2026-02_imcoh-dev-notes/` with
three rule violations baked in (rasterized edges, `suptitle`, PNG
output). Promoting them to `src/lrg_eegfc/visuals/network_templates.py`
+ thin demo scripts is the work.

## Read me first

1. The repo already has a working network-drawing recipe (γ=6 power
   law on edge widths, same-probe red, cross-probe gray, sort-by-weight
   draw order). It just lives in archived dev scripts, not the library.
2. Eight layout strategies have been tried; for our dense weighted
   ImCoh case **only deterministic ones survive** (KK, spectral,
   Laplacian-PCA, LRG-seeded KK, graph-tool SFDP). Plain spring is a
   hairball without rank- or weight-preprocessing.
3. The "LRG metric" layouts already exist in two flavours — KK with
   inflated inter-community distances (NetworkX), and SFDP with
   `groups=` set from the LRG dendrogram cut (graph-tool). Both work;
   SFDP scales better.
4. Three rules from the FC sweep apply here too: PDF-only, no
   `suptitle`, no `set_rasterized(True)`. The archive code violates
   all three; porting must strip them.
5. Probe-shaft coloring is non-negotiable — it is how we visually
   verify that an apparent "module" isn't just a single sEEG shaft
   (same-probe MSC bias). Same-probe edges go red, cross-probe gray.
6. The proposed template family is **`single_network`,
   `row_per_phase`, `mosaic_band_phase`, `mosaic_layout_compare`,
   `lrg_seeded`, `glass_brain_3d`** — six templates mirroring the FC
   structure plus two network-specific (layout-compare, LRG-seeded).
7. Open questions block the build: graph-tool dependency, default
   layout choice, MSC vs ImCoh edge-color rule, 3D scope.

---

## §1 What's in the repo today

### Active library helpers (`src/lrg_eegfc/visuals/`)

| File | Function | Layout | Edge style | Notes |
|---|---|---|---|---|
| `correlation.py:370` | `plot_correlation_and_network` | spring, `seed=42` | `width = w` (raw) | side-by-side heatmap + net; PNG output |
| `msc.py:103` | `plot_msc_and_network` | spring, `k=0.1, iter=50` | `width = w` (raw) | accepts `sparsify={"none","soft"}`; PNG |
| `msc.py:370` | `plot_msc_summary` | spring, `k=0.5, iter=50` | `width = w` (raw) | + percolation curves panel |
| `lrg.py:1004` | `plot_lrg_full_panel` (panel d) | `nx.spectral_layout` | power law `widths = 4·w²` | embedded in 5-panel LRG figure |
| `spatial.py:545` | `plot_spatial_network_3d` (Plotly) | electrode MNI coords | weight-modulated alpha; threshold + max_edges | interactive HTML |
| `spatial.py:714` | `plot_spatial_network_3d_mpl` | electrode MNI coords | per-edge `ax.plot3D`, alpha = `0.2+0.3·w` | static 3D, PNG |
| `spatial.py:1050` | `view_brain_connectome` | nilearn glass-brain | colormap by weight | interactive HTML |
| `spatial.py:1181` | `plot_brain_connectome` | nilearn glass-brain | percentile threshold | static glass-brain |
| `spatial.py:1320` | `plot_brain_connectome_at_n` | nilearn glass-brain | top-`edge_top_pct` | drawn into existing axes |
| `plotting.py:131` | `plot_graph` | caller-provided layout | `width = w` (raw) | low-level CLI helper |

**Verdict:** every active helper is a one-shot wrapper around
`nx.draw(...)`-style primitives. None decompose layout / edge /
decoration as separate concerns. None handle the dense-uniform ImCoh
case (no preprocessing, no γ-power, no same-probe highlighting).

### Active scripts that draw networks

- `scripts/03_analysis/compare_network_variants.py` — corr vs MSC
  dense vs MSC validated. Uses `nx.spring_layout(k=0.2, iter=50)` and
  raw `nx.draw`. PNG. **Rule violations: `suptitle`, PNG, `fig.colorbar(ax=...)`,
  uniform lightblue nodes (no probe coloring).**

That's it for active scripts. Everything else is in `wp1_spatial.py`
(glass-brain orchestrator) or in archive.

### Archive dev scripts (`scripts/archive/2026-02_imcoh-dev-notes/`) — the goldmine

| Script | What it sweeps | Why it matters |
|---|---|---|
| `fig_helper_layouts.py` | 8 layouts × 6 cases | Canonical gallery — the "one of these layouts will work" reference |
| `fig_helper_edge_style.py` | width_range × γ on locked layout | Where γ=6, width=(0.15, 4.0), alpha=(0.03, 0.9) defaults come from |
| `fig_helper_lrg_layout.py` | LRG cut × inter-community separation factor | NetworkX KK driven by LRG dendrogram |
| `fig_helper_lrg_sfdp.py` | LRG cut × SFDP γ | graph-tool SFDP with LRG groups (scales better than NetworkX KK) |
| `fig_helper_rank_transform.py` | (rank/E)^α × {KK, spring k} | Uniform→sharpened weight distribution for layout only |
| `fig_helper_weight_transform.py` | spring k × iterations on w^α | Direct power preprocessing alternative |
| `fig_helper_backbone_layout.py` | top-k% / disparity backbone × {FR, ARF, SFDP} | Layout on sparse backbone, draw full graph |
| `fig_helper_three_layouts.py` | KK / spring / spectral side-by-side | Quick comparison harness |

All eight share the same `_draw_gamma_edges` recipe and `_shaft_colors`
node-decoration. That recipe is what we promote.

---

## §2 The three-layer model (the design contract)

A network template parameterizes three orthogonal concerns:

1. **Layout** — `(N,) → (N, 2)` or `(N, 3)`. Pure positional. Reads
   `A` and (optionally) prior labels (LRG community, electrode
   coords). Returns coordinates only.
2. **Edge representation** — `(pos, A, probes) → drawn LineCollection`.
   Reads positions + the *original* weights. Decides width, alpha,
   color per edge. **Same-probe red, cross-probe gray** is the project
   default. γ-power on `w/w_max` controls dynamic range.
3. **Node decoration** — `(pos, labels) → ax.scatter(...)`. Color by
   electrode shaft (default), LRG community (cluster figures), or
   anatomical region. Size constant 14–20pt unless explicitly
   degree-modulated.

The template `.py` is just a thin wrapper:

```python
A, probes = load_inputs(patient, band, phase, fc_method)
pos = LAYOUT_FN(A, ...)        # layer 1
draw_edges(ax, pos, A, probes)  # layer 2
ax.scatter(pos[:, 0], pos[:, 1], c=node_colors(...))  # layer 3
```

Anything fancier (preprocessing, multiple panels, mosaic) composes
these three primitives. Layout is the only knob that varies between
templates.

---

## §3 Layout strategies — catalog

| Layout | Determinism | Edge weights | Use for | Source |
|---|---|---|---|---|
| Spring (FR) `k=0.01–0.1` | stochastic | yes | only after rank/weight preprocessing | nx, all archive scripts |
| Spring `k=k_base/√N` | stochastic | yes | dense-aware variant; still hairball without preproc | `fig_helper_layouts.py:101` |
| Kamada–Kawai (KK) | deterministic | distance=1/w attribute | **default candidate**; needs distance edge attr | `fig_helper_layouts.py:109` |
| Spectral (λ₂, λ₃) | deterministic | yes (Laplacian) | reveals bipartition; flat for our nearly-fully-connected case | `fig_helper_layouts.py:124` |
| Laplacian-PCA (5→2) | deterministic | yes | richer than spectral when there's no clear bipartition | `fig_helper_layouts.py:193` |
| Community-grouped ring | stochastic (sub-spring) | yes | separates Louvain modules visually but sub-spring inside is unstable | `fig_helper_layouts.py:153` |
| Backbone-guided (disparity α=0.05) | stochastic | layout on sparse, draw full | exposes statistically-significant skeleton | `fig_helper_layouts.py:212` |
| Circular by shaft | deterministic | none | sanity-check baseline; use to confirm modules aren't shaft-trivial | `fig_helper_layouts.py:139` |
| **LRG-seeded KK** | deterministic given dendrogram | distance×separation on inter-comm | explicit LRG-driven module separation (NetworkX) | `fig_helper_lrg_layout.py:79` |
| **LRG-seeded SFDP** | deterministic | yes + groups | graph-tool, scales to dense; **best candidate for the cluster figure** | `fig_helper_lrg_sfdp.py:95` |
| ARF (Geipel 2007) | deterministic post-init | yes | graph-tool alternative to SFDP for very dense graphs | `fig_helper_backbone_layout.py:126` |
| Electrode MNI coords (3D) | fixed | n/a | anatomical reference; for spatial network panels | `spatial.py:427` |

### Lessons applied to our dense-uniform ImCoh case

- **Plain `nx.spring_layout` produces a hairball** because ImCoh
  weights are within ~1 order of magnitude of each other. All edges
  pull with similar force; nothing separates.
- **Three fixes that work, in order of preference:**
  1. **LRG-seeded SFDP** with `groups=` from the LRG cut — modules
     are forced apart by the algorithm, not by a layout trick.
  2. **Rank-transform preprocessing** — replace weights with
     `(rank/E)^α` (α=4–8 works) before the layout. Drawing uses the
     original weights so the visual stays comparable.
  3. **Disparity-filter backbone-guided** — layout on the
     null-model-significant subset, draw the full graph on top.

The recipe is "preprocess for layout, draw with original weights" —
this is documented in `fig_helper_rank_transform.py` and is what we
codify in the template helpers.

---

## §4 Edge-representation lessons (the canonical recipe)

The recipe is reproduced almost verbatim in five archive scripts:

```python
# Inputs: pos (N, 2), A (N, N) original weights, probes list[str]
r, c = np.triu_indices(N, k=1)
w = A[r, c]
active = np.where(w > 0)[0]
w_act = w[active]

t = (w_act / w_act.max()) ** gamma          # γ = 6 default
widths = width_range[0] + t * (width_range[1] - width_range[0])  # (0.15, 4.0)
alphas = alpha_range[0] + t * (alpha_range[1] - alpha_range[0])  # (0.03, 0.9)

order = np.argsort(w_act)                    # ascending → strongest drawn last
segments, colors, linew = [], [], []
for oi in order:
    idx = active[oi]
    i, j = r[idx], c[idx]
    segments.append([pos[i], pos[j]])
    linew.append(widths[oi])
    same_probe = probes[i] == probes[j]
    rgb = (0.8, 0.2, 0.2) if same_probe else (0.35, 0.35, 0.35)
    colors.append(rgb + (alphas[oi],))
ax.add_collection(LineCollection(segments, colors=colors, linewidths=linew, zorder=1))
```

### What it gets right

- **γ=6 power law** compresses the long tail of weak edges so they
  don't fill the plot with gray noise; the strong edges still pop.
- **Sort-ascending draw order** puts strong edges on top of weak ones
  visually, no matter the layout.
- **Same-probe red** is the cohort-specific check: any apparent
  module that's also red-dominated is a same-probe artefact (the MSC
  bias documented in `.agents/guides/02_methods/probe-bias-guide.md`).
- **`LineCollection`** instead of looped `nx.draw_networkx_edges` is
  ~50× faster for ~6.5k edges.

### What needs porting fixes

- The archive code passes `rasterized=True` to `LineCollection`. The
  current rule (`feedback_no_rasterization.md` + the never/always
  list) is **never rasterize** — vector PDFs handle 6.5k thin lines
  fine. **Strip `rasterized=True` on port.**
- Hard-coded `(0.8, 0.2, 0.2)` and `(0.35, 0.35, 0.35)` should
  become module-level constants `SAME_PROBE_RGB` / `CROSS_PROBE_RGB`
  exported from the helpers module so `/figure` can re-skin if needed.
- The recipe assumes magnitude weights `≥ 0`. Signed FC (raw
  correlation, `imcoh`) needs a sign-aware variant — diverging
  colormap on weight, no probe red/gray. **Open question for design.**

### Default parameters I recommend locking

| param | default | rationale |
|---|---|---|
| `gamma` | `6.0` | empirically best across `imcoh_abs` from `fig_helper_edge_style.py` |
| `width_range` | `(0.15, 4.0)` | thinnest still visible; thickest still proportional |
| `alpha_range` | `(0.03, 0.9)` | weak edges nearly invisible; strong edges nearly opaque |
| `same_probe_rgb` | `(0.8, 0.2, 0.2)` | red, calls attention to bias |
| `cross_probe_rgb` | `(0.35, 0.35, 0.35)` | mid-gray; doesn't compete with shaft node colors |

---

## §5 Node-decoration patterns

| Mode | Function | Used for | Source |
|---|---|---|---|
| **Electrode shaft (default)** | `_shaft_colors(probes)` — `tab20` over unique shafts | sanity check; reveals whether modules are probe-trivial | all archive scripts |
| **LRG community** | `_comm_colors(labels)` — `tab20` over cluster ids | the cluster figure; pairs with LRG-seeded layout | `fig_helper_lrg_layout.py:45`, `fig_helper_lrg_sfdp.py:49` |
| **Uniform `lightblue`** | hard-coded | nothing — current `correlation.py` / `msc.py` default; **drop in templates** | active library only |
| **Size by degree (strength)** | `s = 60 + 200·(strength - min)/(max-min)` | glass-brain emphasis on hubs | `spatial.py:1432` |
| **Anatomical region** | (not yet) — would join with `implant_pat_NN.csv` Desikan-Killiany | future cohort-anatomy figure | — |
| **Epileptic-zone highlight** | (not yet) — would mark seizure-onset contacts | future clinical figure | — |

**Default for the templates: shaft coloring.** The cluster-network
template is the only one where LRG-community coloring is the default,
because that figure's whole point is to show LRG modules.

Constant decoration params:
- `node_size`: `15–20` (smaller for n×m mosaics, larger for single
  panels)
- `edgecolors="white"`, `linewidths=0.3`
- `zorder=5` — above edges (zorder=1)

---

## §6 Weight / rank preprocessing — when to apply it

**Rule of thumb:** preprocess only the matrix that goes into the
layout. Always draw with the original weights so panels stay
comparable.

| Transform | Formula | Use when |
|---|---|---|
| **Identity** | `A` | MSC validated (heavy-tailed already) |
| **Power** | `(A / A.max())^α`, α=2 | mild sharpening; doesn't change ordering |
| **Rank** | `(rank(w) / E)^α`, α=4–8 | ImCoh dense (near-uniform weights) |
| **Backbone (disparity)** | `disparity_filter(A, α=0.05)` | when modules are statistically significant but visually faint |
| **Top-k%** | `A * (A ≥ percentile(A, 100−k))` | quick-and-dirty backbone; less principled than disparity |

`fig_helper_rank_transform.py` already implements rank transform.
`fig_helper_weight_transform.py` implements power transform with
spring k × iterations sweep. Both keep `_draw_edges(ax, pos, A, ...)`
on the original `A`.

---

## §7 LRG-driven layouts — the answer to "metric given by the LRG"

Two implementations exist, with different tradeoffs.

### LRG-seeded Kamada–Kawai (NetworkX, `fig_helper_lrg_layout.py:79`)

```python
def lrg_seeded_layout(A, community_labels, separation=5.0):
    D = 1.0 / (A + 1e-6); np.fill_diagonal(D, 0.0)
    inter = community_labels[:, None] != community_labels[None, :]
    D[inter] *= separation
    G = build_graph_with_distance_attr(A, D)
    return np.array([nx.kamada_kawai_layout(G, weight="distance")[i] for i in range(N)])
```

- **Knobs:** `n_communities` (LRG dendrogram cut level), `separation`
  (multiplicative factor on inter-community distances, 1.5 → 10).
- **Pro:** runs anywhere NetworkX runs; deterministic given dendrogram.
- **Con:** O(N²) memory for distance matrix; KK is slow at N=115.

### LRG-seeded SFDP (graph-tool, `fig_helper_lrg_sfdp.py:95`)

```python
def sfdp_lrg(A, groups, gamma=0.1):
    g, ew = _gt_graph(A)
    grp = g.new_vertex_property("int"); grp.a = np.asarray(groups, dtype=int)
    pos = gt.sfdp_layout(g, eweight=ew, groups=grp, gamma=gamma, ..., max_iter=0)
    return np.array([list(pos[v]) for v in range(A.shape[0])], dtype=float)
```

- **Knobs:** `n_communities`, `gamma` (SFDP module-separation; 0.03 →
  0.6 spans "weak" → "very strong" separation).
- **Pro:** SFDP is multilevel, scales to thousands of nodes; `groups`
  argument is the native API for module-aware layout.
- **Con:** graph-tool is a heavy dep (C++), not in `pip install
  lapbrain`'s base. **Open question: do we make graph-tool a hard
  dependency for the network templates, or keep SFDP as the optional
  premium path with KK as the default?**

### Recommendation

- **Default LRG-aware layout: NetworkX KK with separation factor.**
  It's already wired up, deterministic, and avoids the graph-tool
  dependency.
- **Keep an SFDP variant as `lrg_seeded(..., backend="gt")` for users
  who have graph-tool installed.** Same template, different backend.

---

## §8 graph-tool vs NetworkX

graph-tool exists in the codebase only in the three archive scripts
listed in §1 (`fig_helper_gt_layout.py`, `fig_helper_lrg_sfdp.py`,
`fig_helper_backbone_layout.py`). No active library or script imports
it.

| | NetworkX | graph-tool |
|---|---|---|
| Install | pure Python, in lapbrain env | C++, not currently in lapbrain env |
| Speed (N=115, dense) | KK ~5 s | SFDP ~0.3 s |
| Edge-weighted layouts | spring, KK | SFDP, ARF, FR (proper) |
| Module-aware layouts | not native (we hack via distance matrix) | native (`groups=`) |
| SBM / hierarchical SBM | none | first-class |

**Verdict for templates:** ship NetworkX-only by default. Add
graph-tool as `optional_backend="gt"` on layouts that benefit
(SFDP, ARF). If user later wants SBM-driven layouts (the user
explicitly mentioned "SBM-inspired"), we promote graph-tool to a
core dep at that point.

---

## §9 3D variants — already mature

Three independent 3D paths exist and they don't need a template
overhaul, only documentation:

- `plot_spatial_network_3d` (Plotly, interactive HTML) — for
  notebooks and presentations.
- `plot_spatial_network_3d_mpl` (matplotlib 3D, static PNG/PDF) —
  for paper figures, when 3D angle is informative.
- `view_brain_connectome` / `plot_brain_connectome` (nilearn glass
  brain) — for anatomical context. Always zoom to electrode bbox
  per `feedback_brain_connectome_zoom.md`.

**Template scope decision:** keep 3D out of the v0 template family.
Add a `glass_brain_3d` template later only if we generate >5 such
figures for the paper. The 2D dense-network problem is the
load-bearing one.

---

## §10 Anti-patterns observed (to forbid in templates)

| Anti-pattern | Where seen | Why it's wrong |
|---|---|---|
| `nx.draw(G, ...)` directly | `correlation.py`, `msc.py`, `compare_network_variants.py` | no width/alpha modulation, no zorder, no LineCollection speedup, no probe coloring |
| `node_color="lightblue"` | `correlation.py:359`, `msc.py:222` | erases probe / community structure |
| `set_rasterized(True)` on edge LineCollection | every archive `_draw_gamma_edges` | violates `feedback_no_rasterization.md` |
| `fig.suptitle(...)` | every archive `run_*` function | violates `feedback_no_suptitles.md` |
| `fig.colorbar(im, ax=ax, fraction=...)` | `compare_network_variants.py:68` | use `imshow_colorbar_caxdivider` |
| PNG output | `correlation.py`, `msc.py`, `compare_network_variants.py`, `fig_helper_*.py` (via `save_helper_fig`) | violates PDF-only rule |
| Plain spring without preprocessing | active `correlation.py` / `msc.py` | hairball on dense ImCoh |
| Threshold-by-quantile for the *layout* matrix | `compare_network_variants.py:32` | drops the layout's structural information; threshold the *backbone* (disparity) instead |
| Mixing layout and drawing in one `nx.draw` call | active library | breaks the three-layer contract |

---

## §11 Proposed template family (mirror of `fc_templates/`)

| Template | Layout | When to use | New library function |
|---|---|---|---|
| `single_network` | configurable (`spring`, `kk`, `spectral`, `lrg_kk`) | one patient × one band × one phase | `plot_fc_network` |
| `row_per_phase` | fixed across panels | one patient × one band × {RPre, TL, TT, RPost} | `plot_fc_network_row` |
| `row_per_band` | fixed across panels | one patient × one phase × 6 bands | reuse `plot_fc_network_row` |
| `mosaic_band_phase` | fixed across panels | one patient × 4 phases × 6 bands | `plot_fc_network_grid` |
| `mosaic_layout_compare` | one row per layout | layout sanity check (KK / spectral / LRG-KK / SFDP) | `plot_layout_gallery` |
| `lrg_seeded` | LRG KK or SFDP | the cluster-figure (shows LRG modules in space) | `plot_fc_network_lrg` |

Each template:
- ships a `<name>.md` style sheet (inputs, defaults, anti-patterns)
- ships a `<name>.py` thin demo (calls the library function with
  argparse over patient/band/phase/`fc_method`)
- writes a vector PDF to
  `data/outputs/figures/network_templates/<name>/`

### Library port plan (concrete)

Create `src/lrg_eegfc/visuals/network_templates.py` exporting:

```python
# layer 1: layouts (return (N, 2) np.ndarray)
def layout_spring(A, k_base=5.0, seed=42)
def layout_kk(A, weight="distance")
def layout_spectral(A)
def layout_laplacian_pca(A, n_components=2)
def layout_circular_by_shaft(probes)
def layout_community_grouped(A, seed=42)
def layout_backbone_guided(A, alpha=0.05)
def layout_lrg_kk(A, community_labels, separation=5.0)
def layout_lrg_sfdp(A, community_labels, gamma=0.1)  # optional, graph-tool

# preprocessing (return modified A; never used by drawing)
def rank_transform(A, alpha=8.0)
def power_transform(A, alpha=2.0)
def disparity_backbone(A, alpha=0.05)

# layer 2: edges
SAME_PROBE_RGB = (0.8, 0.2, 0.2)
CROSS_PROBE_RGB = (0.35, 0.35, 0.35)
def draw_gamma_edges(ax, pos, A, probes, *, gamma=6.0, width_range=(0.15, 4.0), alpha_range=(0.03, 0.9))

# layer 3: nodes
def shaft_colors(probes)
def community_colors(labels)
def draw_nodes(ax, pos, *, c, s=18, edgecolors="white", linewidths=0.3, zorder=5)

# top-level templates
def plot_fc_network(A, probes, *, ax, layout="kk", node_color="shaft", ...) -> ax
def plot_fc_network_row(patient, band, fc_method, *, phases=PHASE_LABELS, layout="kk", ...) -> fig
def plot_fc_network_grid(patient, *, fc_method, layout="kk", ...) -> fig
def plot_fc_network_lrg(patient, band, phase, fc_method, *, n_communities=10, separation=5.0, backend="nx") -> fig
def plot_layout_gallery(A, probes, *, layouts=DEFAULT_GALLERY) -> fig
```

`__init__.py` re-exports the user-facing entry points.

**Estimated diff:** ~600 LOC library (porting + cleanup of archive
helpers), ~80 LOC per demo script × 6 templates, 6 `.md` style
sheets. Two days of focused work.

---

## §12 Design decisions (resolved 2026-05-05)

1. **graph-tool: REQUIRED.** Verified `graph_tool 2.63` in `lapbrain`.
   Library ships a `nx_to_gt(G)` and `matrix_to_gt(A)` adaptor so any
   layout can be called via either backend without caller code change.
2. **Default layout: `spring`.** All layouts must be testable via
   `mosaic_layout_compare` — we don't pick a winner up front. Spring
   stays the default for `single_network`.
3. **Signed-FC: SUPPORTED in v0.** Edge renderer takes a
   `coloring` kwarg:
   - `coloring="probe"` — same-probe red / cross-probe gray (default
     for magnitude FC: `imcoh_abs`, `msc`, `|corr|`)
   - `coloring="signed"` — diverging `RdBu_r` colormap on signed
     weight, `TwoSlopeNorm` centered at 0 (for raw `imcoh`, raw
     `corr`). Width/alpha still drive off `|w|`.
4. **No sparsification by default.** Full dense graph drawn; the
   γ-power on width/alpha is what carries the visual structure.
   `disparity_backbone` and `top_k_backbone` remain available as
   opt-in preprocessing for users who want them.
5. **Keep legacy `plot_correlation_and_network` /
   `plot_msc_and_network`** as the CLI back-ends. New templates do
   not replace them; they live alongside in
   `lrg_eegfc.visuals.network_templates`. The CLI may later be wired
   to the new helpers, but that's a separate change.
6. **Same-probe rule for non-sEEG cohorts.** Open. Pat_06's scalp
   EEG `f3`, `c3` will trivially form their own "probe family". In
   v0 we accept this; flag added to docstring; revisit if it becomes
   visually noisy.

---

## §13 Where this lives once approved

```
.agents/guides/05_plotting/network_templates/
├── README.md              # this file (review + design)
├── single_network.md      # style sheet
├── single_network.py      # thin demo
├── row_per_phase.md
├── row_per_phase.py
├── mosaic_band_phase.md
├── mosaic_band_phase.py
├── mosaic_layout_compare.md
├── mosaic_layout_compare.py
└── lrg_seeded.md
    lrg_seeded.py

src/lrg_eegfc/visuals/network_templates.py   # library

data/outputs/figures/network_templates/
├── single_network/
├── row_per_phase/
├── mosaic_band_phase/
├── mosaic_layout_compare/
└── lrg_seeded/
```

After approval, the work order is: (a) port helpers to
`network_templates.py` with rule fixes, (b) build `single_network`
end-to-end as the proof-of-template, (c) extend to `row_per_phase`
+ `mosaic_band_phase`, (d) wire `lrg_seeded` once we settle the
graph-tool question, (e) update `fc_templates/README.md`'s routing
table to point here for "network", "graph", "spring layout",
"connectome" requests.
