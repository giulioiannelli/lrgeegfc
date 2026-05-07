# Archive — 2026-02 imcoh-dev-notes

**Development-era figure drafts and layout experiments written while
the ImCoh framework was still being prototyped (Feb–Mar 2026). Preserved
because they contain useful visualization helpers and layout recipes
that informed the canonical Section 2 figures.**

## Contents

- `fig_AB/C/C2/D/E/F/G/G2_*.py` — exploratory figure drafts
  (adjacency, spectral distribution, distributions, enrichment, network,
  metrics-on-graph, eigenvalue, susceptibility).
- `fig_helper_*.py` — layout helpers (backbone, edge style, GT layout,
  layouts, sweep, LRG SFDP, rank/weight transform, etc.).
- `fig_01_probe_bias.py` — early probe-bias visualization.
- `_shared.py` — old copy of the stats helpers (now in library).
- `section3/gen_fig_H/I/J_*.py` — Section 3 figure drafts
  (dendrograms, brain communities, concordance comparison).
- `data/` — local figure data (gitignored, kept on disk for reference).

## Status

Don't re-run as a batch. Canonical Section 2 visualization recipe lives
in `data-layout.md` + `imcoh_network_drawing_recipe` memory. Individual
helpers may be worth porting to `src/lrg_eegfc/visuals/` if reused
elsewhere — see `.agents/guides/04_rules/coding-rules.md` (≥2 callers
→ library).
