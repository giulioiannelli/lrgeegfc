---
name: writing_directive_2026-07-13_sparsified-recovery-OVERVIEW
kind: writing-directive
era: IMCOH_ABS_COHORT_N10
status: current
created: 2026-07-13
scope: shared context for the three parallel Results-rewrite agents — the mst@0.20 sparsified-recovery story, the new τ-swept pipeline, the master number tables, the honesty rules, and how to weave the role of τ into every section. READ THIS FIRST, then open your section's own diff file.
pointers:
  - .agents/reports/2026-07-12_mst020-recovery-arc.md            # full source of truth
  - .agents/preprint/overleaf/results_sec_1.tex                  # §Trace (agent 1)
  - .agents/preprint/overleaf/results_sec_2.tex                  # §Abstraction / enc–inf (agent 2)
  - .agents/preprint/overleaf/results_sec_3.tex                  # §Epileptogenic marker (agent 3)
  - .agents/preprint/directives/writing_directive_2026-07-13_sec1-trace-conceptual-diff.md
  - .agents/preprint/directives/writing_directive_2026-07-13_sec2-encinf-conceptual-diff.md
  - .agents/preprint/directives/writing_directive_2026-07-13_sec3-epi-conceptual-diff.md
---

# Sparsified-recovery Results rewrite — shared overview

## Head

The Results were written for a **dense, fully-connected** diffusion propagator
read at a **single** time \(\tau = 1/\lambda_{\max}\). We have since shown that
on the fully-connected graph the propagator is **degenerate** — at \(\tau_{\min}\)
the heat-kernel geometry is indistinguishable from the raw connectivity, so
"diffusion" was buying nothing (audit_174). The fix is a **minimal
sparsification**: keep each network's **maximum spanning tree plus its strongest
20 % of edges** (`mst@0.20`), which makes the propagator non-degenerate and turns
\(\tau\) into a genuine **multiscale scanner**. On this backbone the old trace
**reproduces exactly at \(\tau_{\min}\)** and then **sharpens and resolves by band
as \(\tau\) is swept up**. Three things change in the paper: (1) the pipeline gains
a sparsification step and \(\tau\) becomes a swept axis, not a fixed operating
point; (2) a **controls ladder** now shows the band-selective trace is unique to
the cophenetic (multiscale) read-out — simple pairwise, network, geodesic, and
even a **spectral** method all fail; (3) the inter-patient spread is **explained**
(electrode laterality), not just controlled. One result **weakens honestly**
(inference is no longer "\(\beta\) alone"); the epilepsy marker is **unchanged**
(robust to the scheme). This file is the shared context; your section's diff file
has the specifics.

---

## 1. What changed in the method (one paragraph you can lift)

We read functional connectivity through a **diffusion propagator** — the heat
kernel \(e^{-\tau L}\) on the Laplacian of the connectivity graph — and cluster the
resulting communication geometry into a hierarchy whose **cophenetic** distances
we compare across task phases (\(\rho_{\mathrm{sym}}\), the split-half-symmetrized
cross-phase correlation; Methods). Two things are new relative to the earlier
draft. First, the graph is **sparsified to a backbone** — the maximum spanning
tree together with the strongest 20 % of remaining edges — because on the fully
connected graph the propagator is degenerate at the fine scale and adds nothing to
the raw weights. Second, the diffusion time \(\tau\) is **swept** across scales
(dimensionless \(s \equiv \tau\lambda_{\max}\), from \(s=1\) at the fine operating
point up to \(s\approx 180\)) rather than fixed, so each band's trace is read as a
**function of scale** rather than at a single point. The only null is the
**strength-matched surrogate**; the drift/timeshift nulls are dropped because the
task is directional and a genuine trace is monotonic by construction, which makes
those nulls degenerate with the alternative they were meant to exclude.

**Why `mst@0.20` and not another backbone (the robustness point).** The \(\beta\)
trace is **sparsification-invariant** — it survives on the dense graph, on every
`mst@f`, on the percolation backbone, and on TMFG (best-scale \(\rho\) ≈
0.20–0.30 throughout). The \(\alpha\) trace is **not** scheme-free: it needs a
backbone that keeps the **strongest** edges. `mst@0.20` retains \(\alpha\)
(best-scale \(\rho\) ≈ 0.20); the parameter-free **percolation** backbone at the
*same density* **destroys** it (\(\rho\) ≈ 0.08), because percolation keeps a
spanning skeleton rather than the heaviest links. So the reported backbone is the
minimal one that preserves **both** carrier bands. Present this as *robustness of
the finding*, not as a tuned choice — the headline is scheme-independent for
\(\beta\) and requires only "keep the strong edges" for \(\alpha\).

---

## 2. The role of \(\tau\) — the through-line to weave into every section

\(\tau\) is the **scale knob** of the analysis, and it is the single thread that
ties the three sections together. Wherever a section currently reads the trace "at
\(\tau=1/\lambda_{\max}\)", it should now read it "across scales, with \(\tau\) the
diffusion depth." The concrete facts to carry:

- **At the fine scale (\(s=1\), \(\tau_{\min}\)) the earlier numbers are recovered
  exactly** — \(\beta\) \(p=0.014\), \(\alpha\) \(p=0.024\) (the \(\alpha\) value is
  the *same* as the published full-graph gate). So nothing already claimed is lost;
  \(\tau\)-sweeping is an *addition*.
- **Sweeping \(\tau\) up separates the bands by their scale signature.** \(\beta\)
  is **scale-invariant** (clears at all 16 scales, deepening to \(p=0.001\) at the
  mesoscale). \(\alpha\) is **mesoscale-tuned** (clears 12/16, peaks near \(s\approx
  5\), \(p=0.007\)). Non-tracing bands stay silent at *every* scale (\(\theta\),
  \(\mathrm{low}\text{-}\gamma\): 0/16). This band-differentiated \(\tau\)-signature
  **is** the multiscale content of the trace.
- **\(\tau\) is what a single-scale or non-diffusion method cannot supply.** The
  band-selectivity lives in *how the trace depends on scale*; read at one scale, or
  through raw edges, or through a spectral distance, the bands are no longer
  separable (§controls ladder). This is the mechanistic reason the hierarchy
  succeeds where simpler measures fail — put \(\tau\) at the centre of that
  explanation.
- **Read \(\tau\) per-scale, never at its best value.** Reporting the single most
  significant scale is a forking path that falsely lights up \(\delta\) and
  high-\(\gamma\) at coarse \(\tau\) where the surrogate median collapses. Those are
  **artifacts, not traces**, and must not be claimed. The honest statement is
  "clears at N of 16 scales including the fine and mesoscale operating points."

---

## 3. Master numbers (fresh, `mst@0.20`, matched-strength, per-scale)

**Trace gate** — `13_matched_strength_mst020` / `05` T\_test (two runs agree):

| band | \(p\) at \(s=1\) (\(\tau_{\min}\)) | \(p\) at mesoscale (\(s\approx5.6\)) | scales cleared | verdict |
|---|---|---|---|---|
| \(\beta\) | 0.014 | **0.001** | **16/16** | scale-invariant trace (flagship) |
| \(\alpha\) | **0.024** | 0.007 | 12/16 | mesoscale trace (peak \(s\approx5\)) |
| \(\theta\) | ns | ns | 0/16 | clean null |
| low-\(\gamma\) | ns | ns | 0/16 | clean null |
| \(\delta\) | fail | fail | 8/16 (non-contiguous) | **scale-artifact — not claimed** |
| high-\(\gamma\) | fail | ns | 4/16 (coarse only) | **scale-artifact — not claimed** |

**Controls ladder** — gate \(p\) for T\_test (the trace), `14_controls_ladder_mst020`.
Bold = clears (\(p<0.05\)):

| descriptor | \(\delta\) | \(\theta\) | \(\alpha\) | \(\beta\) | low-\(\gamma\) | high-\(\gamma\) |
|---|---|---|---|---|---|---|
| raw FC (edges) | **.024** | .053 | **.019** | **.024** | **.032** | .065 |
| node strength | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| clustering | .14 | .54 | **.019** | .22 | .14 | .25 |
| geodesic | **.042** | .053 | .065 | .065 | **.014** | .25 |
| **resistance (spectral)** | .19 | .50 | .54 | .12 | .50 | .58 |
| cophenetic (\(\tau_{\min}\)) | .19 | .65 | **.019** | **.014** | .39 | .16 |
| **cophenetic (mesoscale)** | .12 | .50 | **.007** | **.001** | .22 | .16 |

Reading: **only the cophenetic read-out is band-selective (\(\alpha\)/\(\beta\),
nothing else).** Raw edges fire four bands (non-selective). Strength is degenerate
(the surrogate preserves it exactly → \(p=1\)). Clustering catches \(\alpha\) only.
Geodesic fires the **wrong** bands (\(\delta\)/low-\(\gamma\), misses \(\alpha\)/\(\beta\)).
The **spectral resistance distance is dead everywhere** — a multiscale-blind
Laplacian read-out that *should* be the strongest foil and finds nothing.

**Encoding / inference** — `05_enc_inf_arc` (mst020, T\_learn bug fixed), per-scale:

| component | band | \(p\) at mesoscale | scales | note |
|---|---|---|---|---|
| encoding \(T_{\mathrm{learn}}\) | \(\beta\) | 0.003 | 16/16 | scale-invariant; **cophenetic-UNIQUE** (raw/clustering/geodesic/spectral all fail) |
| encoding \(T_{\mathrm{learn}}\) | \(\alpha\) | 0.032 | 3/16 | marginal (0.053 in ladder run); **readable from raw FC + clustering = basic metric** |
| inference \(T_{\mathrm{inf}\cdot e}\) | \(\beta\) | 0.005 | 7/16 | mesoscale; **also readable from raw FC (0.010)** → value-add = localization |
| inference \(T_{\mathrm{inf}\cdot e}\) | \(\alpha\) | 0.024 | 7/16 | mesoscale |
| inference \(T_{\mathrm{inf}\cdot e}\) | \(\delta\) | 0.019 | 10/16 | mesoscale |

> **Encoding is NOT washed out** (an earlier bug faked that). The band-split is the point:
> **β encoding is the multiscale-exclusive result** (only the hierarchy sees it), while
> **α encoding lives in a basic metric** (raw FC/clustering read it directly). For β the
> "encoding=basic / inference=multiscale" story is in fact **inverted** — β *encoding* is
> hierarchy-only; β *inference* is raw-visible (its coph value-add is localization, not
> detection). Do not carry the stale "encoding strength-explained / null in multiscale".

**Epilepsy marker** — `06_epi_arc` (mst020), strength-residual seed affinity AUC:

| band | AUC single-\(\tau\) | AUC multi-\(\tau\) | beats null | multi>single | prec@5 (single band) |
|---|---|---|---|---|---|
| \(\delta\) | 0.81 | 0.83 | 7/10 | 7/10 | 0.40 |
| low-\(\gamma\) | 0.76 | 0.82 | 8/10 | 7/10 | 0.40 |
| \(\beta\) | 0.72 | 0.745 | 8/10 | 8/10 | 0.30 |
| \(\alpha\) | 0.61 | 0.61 | 5/10 | 5/10 | 0.20 |

**Inter-patient variation (pillar)** — `15_patient_variation`: the \(\beta\) trace
is **left-lateralized** — its per-patient strength tracks the left-hemisphere
contact fraction (Spearman \(\rho = +0.685\), \(p = 0.029\)); left-only implants
trace, the right-only implant (Pat_15) runs null. \(\alpha\) is **not** lateralized
(\(\rho = 0.16\), ns) — a trace with no anatomical address, consistent with its
whole-brain-but-placeless localization.

---

## 4. Honesty rules (non-negotiable — this is a brutal-honesty project)

1. **Per-scale, never best-scale.** Any claim reads at the fine and/or mesoscale
   operating points and reports "N/16 scales." \(\delta\) and high-\(\gamma\) clear
   only at coarse \(\tau\) where the surrogate collapses → artifacts, **not** traces.
2. **Matched-strength is the only null.** Do not reintroduce drift / timeshift /
   within-baseline framing as verification.
3. **State what weakened up front.** The "inference in \(\beta\) alone" headline does
   **not** survive (§sec 2). Say so plainly; do not launder it.
4. **The multiscale ladder is a fraction, not a count.** On the dense graph 100 % of
   patient×band cells have a single specific-heat peak (degenerate); at `mst@0.20`
   ~38 % develop a second peak (percolation ~45 %), and ~0 % above density 0.5. Do
   **not** write "the C(τ) curve has several peaks" as if universal — the defensible
   multiscale evidence is the **band-differentiated \(\tau\)-signature** (β broad, α
   mesoscale) plus the controls ladder, not a multi-peak thermodynamic curve.
5. **Don't invent numbers.** Everything citable is in §3 or your section diff; if a
   value isn't there, mark it a placeholder and flag it, don't guess.

---

## 5. What is NOT touched by this rewrite

- **The localization results** (β→OFC, encoding→OFC, inference→cingulate,
  low-γ→cingulate memory, SOZ-divergence, tissue-class carrier) are downstream of
  the same cophenetic trace and **stand as written** — the anatomy did not change,
  only the pipeline that produces the per-pair cophenetic distances feeding it.
  Localization is, in fact, now the *headline value-add* for inference (§sec 2).
- **The Grassmann table** (sec 1) is a complementary spectral read on the same
  operator; it was matched-strength-verified on the dense/whole-task scheme
  (audit_66) and is **not re-derived** on the sparsified backbone. Keep it, but do
  not re-attribute it to `mst@0.20`; if in doubt, treat it as a methods-level
  companion, not a sparsified-pipeline result.
- **Reinstatement / held-not-replayed** (sec 1) inherits the static trace's control
  and is unchanged.
