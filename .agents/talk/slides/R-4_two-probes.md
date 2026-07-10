---
name: talk-slide-R4-two-probes
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-07-10
updated: 2026-07-10
slide: R-4
part: III — Results
duration: ~90 s
pointers:
  - .agents/reports/2026-07-07_talk-structure-20min.md
  - scripts/01_compute/audit/audit_169_grassmann_drift_null.py
---

# R-4 · It's real — and it survives the strongest null

## Slide placeholder (copy into Canva)

**Main point.** Three ways of reading the same connectivity — **raw edges**, the
**cophenetic tree** (fine/multiscale), the **Grassmann subspace** (coarse/global) —
against **two nulls at once**: matched-strength *and* drift. Under matched-strength the
tree and the subspace **both** clear β (corroboration). Under the **stronger drift
null**, **raw drowns**, **Grassmann is drift-vulnerable**, and **only the cophenetic β
survives both** — the fine multiscale tree is the one drift-robust probe.

**Concepts to land.**
- **A band must overcome BOTH nulls to count** — matched-strength (is it more than a
  strength change?) *and* drift (is it more than slow session non-stationarity?).
- **Matched-strength — the corroboration:** cophenetic clears **β + α**, Grassmann
  clears **β + γ_low + δ**; the only band **both** independent strength-immune probes
  agree on is **β** (a clean double dissociation — α is tree-only, γ_low/δ subspace-only).
- **Drift — the discriminator:** **raw clears *nothing*** (its big per-patient values are
  matched by equally big drift; paired p ≥ 0.19 every band → it **drowns**). **Cophenetic
  clears β** (real 0.198 vs drift 0.051, **p = 0.042**) — drift *removes α*. **Grassmann
  clears nothing** — β only **trends** (real 0.425 vs drift 0.235, **p = 0.116**).
- **So cophenetic β is the sole survivor of the dual gate.** The **fine multiscale tree
  is drift-robust; the coarse global subspace and the raw edge are drift-vulnerable** —
  the same message as "the trace is fine-grained, not coarse-spectral."

**Figures / visuals.**
- **Dual-null bloom triptych** — one radial plot per measure (raw · cophenetic ·
  Grassmann). Per-patient **triangles**; **grey annulus = matched-strength null**;
  **red dashed waterline = drift null**; a **wedge glows only if the band clears BOTH**.
  Reads at a glance: raw = *no glow* (petals impaled on the red drift line), cophenetic =
  *β glows alone*, Grassmann = *no glow* (β trends, doesn't clear). Paths:
  `data/preprint/figures/_drafts/fig_dual_bloom_raw_DRAFT.pdf` ·
  `…/fig_dual_bloom_coph_DRAFT.pdf` · `…/fig_dual_bloom_grassmann_DRAFT.pdf`
  (builder `scripts/01_compute/figures_embedded/fig_reasoning_bloom.py --dual`).
- *(MS-only alternative, if you want the plain double dissociation:*
  `data/outputs/figures/talk/fig_coph_grassmann_confirmation_matrix.pdf` — β the only
  band both probes clear under matched-strength.*)*

**References.** None external — the Grassmann subspace distance and the drift null are
our own controls (audit_66 / audit_169). The LRG operator is Villegas 2023/2025 (M-3).

---

## Keep honest (content constraints, not styling)

- **Grassmann β is drift-UNVERIFIED, not drift-clean and not "confounded."** It trends
  positive (real 0.425 > drift 0.235) but does **not** clear the windowed drift
  (p = 0.116). Never say "Grassmann survives drift" (false) — say "Grassmann corroborates
  β under matched-strength but is not drift-robust / drift-untested."
- **Which drift matters — be precise.** The dual gate here uses the **windowed** drift
  (audit_167/169), which is the **conservative** construction: it also fails cophenetic
  **α** (p = 0.188), whereas the fair full-duration **C2** rescues α (drift-clean 0.0068).
  There is **no fair full-duration Grassmann drift** (would need a C2-style build). So
  "β only" for cophenetic is under the *windowed* drift; the flagship β itself is
  drift-clean under **both** windowed (0.042) and fair C2 (0.0137).
- **Matched-strength "two probes agree on β" still holds** — that is the corroboration,
  and it is a real double dissociation (α tree-only, γ_low/δ subspace-only). Don't let the
  drift result erase it; state both layers.
- **The drift gate is PAIRED** (each patient's real vs its own drift). The red arc is the
  *cohort* drift level (illustrative); the **glow uses the correct paired Wilcoxon** — so
  read the verdict from the glow, not from petals-vs-arc.
- **Never call any of this "nonlinear."** |ImCoh| is second-order; tree/subspace are
  higher-order **graph** features emergent from pairwise connectivity.
- New this session: `audit_169_grassmann_drift_null.py` (whole-task Grassmann windowed
  drift; valid because whole-task T_G is a simple, non-conditional contrast).
