---
name: talk-slide-14-detect-vs-discriminate
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
slide: 14
status: draft
updated: 2026-07-14
canva: REWRITTEN 2026-07-14 around the ACTUAL message — the trace is a genuine MULTISCALE
  reorganization. (1) the whole standard network toolkit (pairwise edges, strength, clustering,
  geodesic, resistance) cannot resolve it; (2) the FIELD-STANDARD spectral-clustering measure
  (Grassmann subspace distance) is MISLED — it flags γ_low, a global/whole-graph spectral shift
  that is NOT a multiscale reorganization (tell: the only other method flagging γ_low is the graph
  geodesic, the other global metric), and it MISSES α; (3) only the multiscale hierarchy isolates
  the true carriers α+β and is not fooled by the global γ_low. Figures: fig_grassmann_contrast (the
  full read-out ladder, centerpiece), fig_pedestal_vs_peak (magnitude intuition), and optionally
  fig_descriptor_bloom_contrast (per-patient spread). Sole null = matched-strength.
---

# Slide 14 — Detect ≠ discriminate: the trace is a multiscale reorganization

1. TITLE
Detect ≠ discriminate — the trace is a *multiscale* reorganization the standard toolkit can't
resolve and spectral clustering mis-reads

2. MAIN CONCEPT
- **HEAD (the juice):** the trace lives in **α and β** — and it is a genuinely *multiscale*
  reorganization, not something any network measure would catch. On the *same* data and the *same*
  matched-strength null: the entire **standard network toolkit is blind or non-selective**, the
  **field-standard spectral-clustering method is actively misled** (it flags a *global* γ_low
  structure that is **not** a multiscale reorganization, and it **misses α**), and **only the
  multiscale hierarchy names the true carriers α+β** without being fooled.
- **We tested the standard toolkit (Fig B, top five rows).** Every everyday network-science
  descriptor on the same mst@0.20 backbone against the same null, representativeness-gated
  (cohort p<0.05 AND leave-one-out AND ≥3/10 patients above own null):
  - **raw FC edges** (pairwise) → β only (α, γ_high are cohort-significant but LOO-fragile). Raw
    detects β because β is a strong pairwise convergence — but it cannot name α.
  - **node strength / clustering coef. / effective resistance** → **blind** (nothing representative).
  - **graph geodesic** (shortest paths — a *global* metric) → **γ_low only**, a non-carrier band.
  - *(Read as white = no trace.)* Most of this block is white: the standard toolkit does not
    resolve the reorganization.
- **The field standard is MISLED (Fig B, Grassmann row).** The go-to method for comparing two
  graphs is **spectral clustering** — the **Grassmann chordal distance** between the k-dimensional
  Laplacian eigen-subspaces of the two phases. Same combinatorial Laplacian `L = D − W` as our
  hierarchy, but it **freezes one subspace dimension k** (a single, essentially *global* spectral
  scale) instead of diffusing across **all** scales `e^{−τL}`. Its verdict: **β** (agrees) **+
  γ_low**, and it **misses α**.
- **γ_low is the tell (Fig B, the γ_low column).** γ_low is flagged by **exactly the two global
  methods** — the graph geodesic and Grassmann — and by **nothing else**, and **never by the
  multiscale read**. So γ_low is a **whole-graph / global spectral shift, not a multiscale
  reorganization**: single-scale spectral clustering conflates a global change with structure.
  That is an argument *for* multiscale, not a neutral "different answer."
- **Only multiscale resolves the carriers (Fig B, bottom row).** The cophenetic hierarchy lands on
  **α and β and nothing else**: it catches **α**, which lives at a scale a fixed subspace cannot
  isolate (α is multiscale-*exclusive*), and it is **not fooled by the global γ_low**. β is the
  convergence carrier all real methods agree on.
- **The through-line:** the multiscale read is the tool *matched to the phenomenon*. Raw FC stays a
  robust, complementary detector of β; the standard toolkit is blind to the reorganization; the
  field-standard spectral method mistakes a global shift for structure and misses α. Multiscale is
  what turns detection into *resolution*.

3. ON-SLIDE TEXT

same data · one null (matched-strength) — is the α/β trace really MULTISCALE?

we tested the standard network toolkit — and the field standard:
  raw edges → β (pairwise) · strength / clustering / resistance → blind
  graph geodesic → γ_low (a global path metric) · **Grassmann (spectral clustering) → β, γ_low — misses α**
  **cophenetic (multiscale · ours) → α, β — and nothing else**

γ_low is flagged ONLY by the two GLOBAL methods (geodesic, Grassmann) — never by the multiscale read:
  ⇒ γ_low = a whole-graph / global spectral shift, NOT a multiscale reorganization
  ⇒ α is multiscale-exclusive · β is the convergence carrier

the field standard — Grassmann chordal distance between Laplacian eigen-subspaces (spectral clustering):

  L = D - W,\qquad U_X=[\,u_1,\dots,u_k\,]\ \ (k\ \text{slowest modes of phase }X)

  d_G(X,Y)=\sqrt{\,k-\textstyle\sum_i \sigma_i^{2}\,}=\sqrt{\textstyle\sum_i \sin^{2}\theta_i},\qquad
  \sigma_i=\cos\theta_i=\mathrm{svd}\!\left(U_X^{\top}U_Y\right)

  one FIXED subspace k (a single, global scale) — not the full diffusion e^{-τL} across all scales

---

4. SPEECH  (~95 s)

We found a reorganization in alpha and beta. The real question isn't whether it's there — it's
whether it's genuinely multiscale, or just something any network measure would pick up. So we ran
the whole standard toolkit on the same data, against the same matched-strength null. [Fig B]

Look at the top five rows — the everyday network-science descriptors. Pairwise edges, node strength,
clustering, shortest paths, effective resistance. Look how much of it is white: most are blind. Raw
connectivity sees beta — beta is a strong pairwise convergence, so of course it does — but it can't
name alpha.

Now the sophisticated one, the measure the field actually reaches for when it compares two graphs:
spectral clustering, formalized as the Grassmann distance between the Laplacian eigen-subspaces. Same
Laplacian we build our hierarchy on — the only difference is it freezes one subspace dimension, a
single global scale, instead of diffusing across all scales. It agrees with us on beta. But watch the
low-gamma column. Spectral clustering flags low-gamma — and the only other method that flags it is
the graph geodesic, the other global metric. Two global measures, one band. Low-gamma is a
whole-graph spectral shift; it is not a multiscale reorganization, and single-scale spectral
clustering can't tell the difference — it reads a global change as structure. And it misses alpha
entirely, because alpha lives at a scale a fixed subspace can't isolate.

Only the multiscale hierarchy, the bottom row, lands on both true carriers — alpha and beta — and
nothing else. It isn't fooled by the global low-gamma, and it reaches the scale where alpha lives.
That's the whole argument for multiscale: not that we alone see a trace — the edges see beta, spectral
clustering sees beta — but that we resolve the real reorganization while the standard tools are either
blind or misled.

**Careful:**
- Matched-strength is the SOLE null (drift RETIRED). Figs B/C use the strict representativeness gate
  (cohort p<0.05 + LOO + ≥3/10); Fig A uses the naive cohort gate (magnitude view). State which gate
  when quoting a count.
- The γ_low reading is the POINT: γ_low is flagged only by the two GLOBAL methods (geodesic +
  Grassmann) and never by the multiscale read ⇒ it is a global/whole-graph shift, NOT a multiscale
  reorganization. Do NOT soften this to "Grassmann sees an extra band" — the interpretation is that
  spectral clustering mistakes a global change for structure.
- Do NOT say "only multiscale sees a trace" / "invisible to simple methods" — raw DETECTS β and so
  does spectral clustering. raw = NON-selective but COMPLEMENTARY, never "blind/fragile."
- **Grassmann is WHOLE-GRAPH (dense imcoh_abs FC), NOT the mst@0.20 backbone**; its per-patient T_G
  aggregates over the subspace-dimension grid k (no single "s"). Honest claim: "same Laplacian
  operator, different construction."
- **Grassmann does NOT confirm the α+β carrier set.** Its matched-strength cluster-extent verdict
  (+ LOO): β strong & LOO-stable, γ_low strong, δ significant-but-LOO-fragile, θ/α/γ_high none. It
  AGREES on β, MISSES α, FLAGS γ_low. Never "Grassmann confirms our carriers."
- α is cophenetic-EXCLUSIVE, β a CONVERGENCE band. Never "coph → β only."
- Marker grammar: band colour = clears (bold/fill = representative, ring = sig-but-LOO-fragile),
  **white = no trace**. Only genuine null bands are white.
- Do NOT name β's anatomy here (β is DELOCALIZED; OFC is the *encoding* anchor only). → later slide.
- Higher-order β (coph-beyond-raw) is RETIRED (BH-marginal, q≈.064; clean cell was δ, not β).

---

5. FIGURES  (talk PNGs, transparent white-ink, dark-slide ready, in data/outputs/figures/talk/)

- **Fig B — the full read-out ladder (CENTERPIECE).** `fig_grassmann_contrast.png`
  (gen `scripts/07_figures/talk_fig_grassmann_contrast.py`; standard + cophenetic verdicts from
  `controls_ladder_apples/per_cell.csv` via representativeness gate, Grassmann from
  `data/audit/grassmann_cluster_extent/cohort_summary.csv`). Seven rows in three tiers (dashed
  separators): **standard toolkit** (raw FC · node strength · clustering · graph geodesic ·
  effective resistance) → **spectral clustering** (Grassmann chordal · whole graph) → **multiscale**
  (cophenetic · ours). × 6 bands. FILL = robust · band-colour RING = sig-but-LOO-fragile · white =
  no trace. Two highlighted columns: **β** (convergence carrier) and **γ_low** ("global — not
  multiscale"; lit only at geodesic + Grassmann). α callout "multiscale-only" (only cophenetic).
  Grassmann FORMULA in §3 on-slide text. This one figure carries the whole slide.
- **Fig A — pedestal vs peak (magnitude intuition, optional).** `fig_pedestal_vs_peak.png`
  (gen `talk_fig_pedestal_vs_peak.py`; `controls_ladder/cohort_gate.csv` dense raw + backbone
  control). raw = broad positive pedestal (all 6 positive, clears 4/6: δ .024, α .019, β .024,
  γ_low .032; naive gate), cophenetic = sharp α/β peak (α .007, β .001) with θ negative. θ boxed as
  internal control. Good as a warm-up "detect ≠ discriminate" glance BEFORE the ladder.
- **Fig C — per-patient blooms (spread, optional).** `fig_descriptor_bloom_contrast.png`
  (gen `talk_fig_descriptor_bloom.py`). Three radial blooms — raw · Grassmann · cophenetic — ten
  patient petals per band, wedge glowing where the band clears its cohort gate. Grassmann panel on
  its OWN radial scale (whole-graph T_G, mean over k). Shows the spread behind the ladder's dots.
- **RETIRED:** `fig_coph_beyond_raw.png` (higher-order β — BH-marginal q≈.064). Do NOT resurrect.

⚠ COMPOSITION NOTE: Fig B is now the standalone centerpiece — it carries "toolkit blind + spectral
clustering misled + multiscale resolves." Fig A (magnitude) and Fig C (per-patient spread) are
optional supports; the ladder subsumes their "raw non-selective / we tested many metrics" point. If
the slide is crowded, run **Fig B alone**, or **Fig B + one** of A/C. Decide on deck.

6. REFERENCES
- Grassmann chordal distance / principal angles between subspaces — standard subspace geometry
  (spectral clustering / spectral embedding). Graph geodesic (shortest-path) and effective
  resistance (Laplacian pseudo-inverse; Klein & Randić 1993) are the global path/spectral foils.
  Controls-ladder + matched-strength null (`controls_ladder_apples`) and the Grassmann
  matched-strength surrogate + cluster-extent gate (audit_66 / `grassmann_cluster_extent`) are ours.

7. CANVA STATUS
REWRITTEN 2026-07-14 to the correct message. The empirical "the trace is MULTISCALE" slide: the
standard network toolkit can't resolve it; the field-standard spectral-clustering measure is misled
(flags γ_low = a global shift, not a multiscale reorganization — tell: only geodesic, the other
global metric, agrees; and it misses α); only the multiscale hierarchy names α+β. Centerpiece =
`fig_grassmann_contrast` (full 7-row ladder). Grassmann formula in on-slide text. Canva job: (1)
place Fig B (± Fig A/C, see COMPOSITION NOTE); (2) paste §3 on-slide text incl. the Grassmann LaTeX;
(3) speech → presenter notes. No β→OFC; matched-strength only. What s means + localization → later.
