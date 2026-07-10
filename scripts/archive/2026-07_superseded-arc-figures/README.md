# Superseded R2 encoding/inference figures (2026-07-09)

Three rejected attempts at the R2 "disentangling / dynamics" figure, kept for history.
All SUPERSEDED by `scripts/01_compute/figures_embedded/fig_arc_attractor_3d.py`
(live windowed state-space trajectory; see the scope report
`.agents/guides/task-persistence-investigation/2026-07-08_windowed-consolidation-flow.md`).

- `fig_arc_h1_gated_loop.py` / `fig_arc_h2_windowed_flow.py` — flat 2-D axis panels
  (user: "horrible ... flat 2d panels ... no cool interpretation").
- `fig_arc_ribbons_3d.py` — six phase-scale ribbons in 3-D; effectively 1-D curves,
  and the huge task_test inference peak buried the small rest_post residue
  (user: "really 1d curves, no real need for 3d").

Fix that stuck: a dynamical-systems state-space portrait — amplitude of the task
excursion is a transient; the result is where rest_post SETTLES (limit set).
