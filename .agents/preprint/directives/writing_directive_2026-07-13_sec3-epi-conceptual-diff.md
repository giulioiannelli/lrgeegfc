---
name: writing_directive_2026-07-13_sec3-epi-conceptual-diff
kind: writing-directive
era: IMCOH_ABS_COHORT_N10
status: current
created: 2026-07-13
scope: conceptual diff for the agent rewriting results_sec_3.tex (§"Propagator-inspired markers of epileptogenic tissue"). The marker is ROBUST to the sparsification — smallest diff of the three. Add the honest τ-contribution (ranking, not precision) and the robustness note. Read the OVERVIEW first.
pointers:
  - .agents/preprint/directives/writing_directive_2026-07-13_sparsified-recovery-OVERVIEW.md
  - .agents/preprint/overleaf/results_sec_3.tex
  - .agents/reports/2026-07-12_mst020-recovery-arc.md
---

# §3 (epilepsy marker) — conceptual diff

## Head

This is the **smallest** change of the three, and the change is good news: the
epileptogenic-tissue marker **reproduces on the sparsified `mst@0.20` backbone**, so
it is **robust to the pipeline choice** rather than an artifact of the dense graph.
Every headline stands — \(\delta\) is the strongest band, the band ordering is
\(\delta \gtrsim\) low-\(\gamma \gtrsim \beta\), the cohort splits into a co-diffusing
"community" majority and a right-hemisphere "hub" minority, the fused leave-one-out
detector is calibrated for triage, and the top-five precision is ~60 %. Your only real
edits: (1) add one honest sentence on **what the multiscale \(\tau\)-scan contributes**
— it lifts the ranking/AUC but **not** the top-of-list precision (that 60 % is a
band-**fusion** effect, orthogonal to \(\tau\)); and (2) note the marker is
**scheme-robust** (same numbers on dense and sparsified). Keep everything else.

## The framing (why this section barely moves)

Sections 1–2 needed sparsification because the *cross-phase cophenetic trace* was
degenerate on the dense graph. This section reads the propagator **differently** — as
a **within-recording seed affinity** (heat placed on known SOZ contacts, where does it
flow), not a cross-phase hierarchy comparison — and that read-out was never degenerate.
So re-deriving it on `mst@0.20` mostly **confirms** the published numbers. Present that
as a **feature**: the same operator that carries the cognitive trace localizes
epileptogenic tissue, and it does so **regardless of how the graph is sparsified** —
one more piece of scheme-independence.

## What STAYS (essentially the whole section)

- **Relational, not contact-by-contact**: SOZ contacts form a strength-independent
  co-diffusing group; strength-residual seed affinity separates SOZ from healthy above
  node strength. Unchanged.
- **Band ordering and AUCs**: \(\delta\) strongest, then low-\(\gamma\), then \(\beta\);
  \(\alpha\) marginal. On `mst@0.20`: AUC \(\delta\) 0.83, low-\(\gamma\) 0.82, \(\beta\)
  0.745, \(\alpha\) 0.61 — **matching or slightly above** the published
  0.80/0.74/0.69/0.60. Keep the published values or update to these; either is honest,
  but say they reproduce.
- **Distant-seed discovery** (seed one electrode's SOZ, rank another's; \(\delta\) AUC
  0.72, 8/10). Unchanged.
- **Calibrated LOO detector**: fused across bands, median AUC 0.87, prob 0.38 vs 0.08,
  9/10 above chance, label-shuffle null 0.48; interpretable logistic over nonlinear.
  Unchanged.
- **Two populations**: 8 community (AUC up to 0.99, 5/8 > 0.90) + 2 right-hemisphere hub
  (fail, AUC 0.49/0.57; strength rescues them); per-patient reporting recovers 9/10.
  Unchanged — and note this dovetails with §1's laterality result (the right-hemisphere
  implants are the same ones that fail the β trace).
- **Triage framing**: prec@5 = 60 % (~7× base rate), spans gray/white (41 % of SOZ in
  white matter), occult candidates for prospective testing. Unchanged.

## What CHANGES (two sentences, really)

**(A) The honest \(\tau\) contribution — ranking, not precision.** Sweeping the
diffusion time \(\tau\) (using the multiscale seed affinity rather than a single scale)
**improves the AUC** in the majority of patients — multiscale beats single-scale in
7/10 (\(\delta\)), 8/10 (\(\beta\)), 7/10 (low-\(\gamma\)) — so \(\tau\) sharpens the
**overall ranking**. But it does **not** improve the **top-of-list precision**: the
single-band prec@5 stays ~0.40, and the reported **60 %** comes from **fusing the six
bands**, which is independent of the \(\tau\)-scan. State this cleanly: *the multiscale
scan helps discrimination (AUC) but the precision headline is a band-fusion effect, not
a \(\tau\) effect.* This keeps the section a **marker** (ranking) rather than
overselling \(\tau\) as a precision booster.

**(B) Scheme-robustness note.** One clause: the marker's numbers reproduce on the
sparsified backbone, so the localization of epileptogenic tissue is not an artifact of
leaving the graph dense — it is a property of the diffusion operator itself. This is
the same "robust to the sparsification" thread that runs through §1's β result.

## τ — how to weave it

Physical and short: heat placed on the seizure contacts spreads through the network;
**letting it diffuse deeper (larger \(\tau\)) improves how well the seizure zone
separates as a co-diffusing community** (multiscale > single-scale AUC in the majority),
consistent with \(\delta\) — the slowest, most spatially extended rhythm — carrying the
marker. But the *shortlist* precision is set by fusing bands, not by depth. That's the
whole \(\tau\) story here; don't inflate it.

## Numbers you may cite (fresh, mst@0.20)

- Seed-affinity AUC (single→multi-\(\tau\)): δ 0.81→0.83, low-γ 0.76→0.82, β 0.72→0.745,
  α 0.61→0.61. Beats matched-strength null: δ 7/10, low-γ 8/10, β 8/10, α 5/10.
  Multi>single: δ 7/10, β 8/10, low-γ 7/10.
- prec@5 (single band): δ 0.40, low-γ 0.40, β 0.30, high-γ 0.40 — none reach 0.60;
  the 0.60 is the **fused** LOO detector (unchanged from the published pipeline).
- Everything else (distant-seed, calibrated LOO, two-populations, tissue spread, occult
  candidates): as published; reproduced, not revised.

## Do-not

- Do not claim \(\tau\) / multiscale raises the precision — it raises AUC/ranking only.
- Do not drop the two-populations honesty or the "marker not detector" ceiling.
- Do not best-scale the AUC; report single-vs-multi and beats-null counts.
