---
name: talk-slide-15-detect-vs-discriminate
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20 recovery)
slide: 15
status: draft
updated: 2026-07-16
canva: REWRITTEN 2026-07-14 around the ACTUAL message — the trace is a genuine MULTISCALE reorganization. (1) the whole standard network toolkit (pairwise edges, clustering, geodesic, resistance; node strength excluded — it is the matched-strength control) cannot resolve it; (2) the FIELD-STANDARD spectral-clustering measure (Grassmann subspace distance) is MISLED — it flags γ_low, a global/whole-graph spectral shift that is NOT a multiscale reorganization (tell: the only other method flagging γ_low is the graph geodesic, the other global metric), and it MISSES α; (3) only the multiscale hierarchy isolates the true carriers α+β and is not fooled by the global γ_low. Figures: fig_grassmann_contrast (the full read-out ladder, centerpiece), fig_pedestal_vs_peak (magnitude intuition), and optionally fig_descriptor_bloom_contrast (per-patient spread). Sole null = matched-strength.
---

# Slide 15 — Detection and discrimination

1. TITLE
Detection and discrimination

2. MAIN CONCEPT
- **HEAD (the juice):** the trace lives in **α and β** — and **α in particular is a genuinely *multiscale* reorganization no single-scale measure resolves** (β is a robust low-order convergence every real detector catches). On the *same* data and the *same* matched-strength null: the entire **standard network toolkit is blind or non-selective**, the **field-standard spectral-clustering method is actively misled** (it flags a *global* γ_low structure that is **not** a multiscale reorganization, and it **misses α**), and **only the multiscale hierarchy names the true carriers α+β** without being fooled.
- **We tested the standard toolkit (Fig B, top four rows).** Every everyday network-science descriptor on the same mst@0.20 backbone against the same null, representativeness-gated (cohort p<0.05 AND leave-one-out AND ≥3/10 patients above own null):
  - **raw FC edges** (pairwise) → β only (α, γ_high are cohort-significant but LOO-fragile). Raw detects β because β is a strong pairwise convergence — but it cannot name α.
  - **clustering coef. / effective resistance** → **blind** (nothing representative). (Node strength is omitted — the matched-strength null fixes it by construction, so it is the control, not a probe.)
  - **graph geodesic** (shortest paths — a *global* metric) → **γ_low only**, a non-carrier band.
  - *(Read as white = no trace.)* Most of this block is white: the standard toolkit does not resolve the reorganization.
- **The field standard is MISLED (Fig B, Grassmann row).** The go-to method for comparing two graphs is **spectral clustering** — the **Grassmann chordal distance** between the k-dimensional Laplacian eigen-subspaces of the two phases. Same combinatorial Laplacian `L = D − W` as our hierarchy, but it **freezes one subspace dimension k** (a single, essentially *global* spectral scale) instead of diffusing across **all** scales `e^{−τL}`. Its verdict: **β** (agrees) **+ γ_low**, and it **misses α**.
- **γ_low is the tell (Fig B, the γ_low column).** γ_low is flagged by **exactly the two global methods** — the graph geodesic and Grassmann — and by **nothing else**, and **never by the multiscale read**. So γ_low is a **whole-graph / global spectral shift, not a multiscale reorganization**: single-scale spectral clustering conflates a global change with structure. That is an argument *for* multiscale, not a neutral "different answer."
- **Only multiscale resolves the carriers (Fig B, bottom row).** The cophenetic hierarchy lands on **α and β and nothing else**: it catches **α**, which lives at a scale a fixed subspace cannot isolate (α is multiscale-*exclusive*). **Why it does *not* also flag the global γ_low — the key logical point:** the tree is *built* from diffusion across all scales, but two phases are *compared* by **ρ_sym, the rank (Spearman) correlation of cophenetic distances**, which is **invariant to an across-the-board global drift by construction** (shift/scale every pair equally → ranks unchanged). So the hierarchy is **not a superset** of the global metrics — it is a **relational lens**: it reads *which nodes merge before which*, resolving the relational α while a roughly-uniform global shift (γ_low) cancels. "Doesn't see γ_low" is by design, not a coverage gap. β is the convergence carrier all real methods agree on.
- **The through-line:** the multiscale read is the tool *matched to the phenomenon*. Raw FC stays a robust, complementary detector of β; the standard toolkit is blind to the reorganization; the field-standard spectral method mistakes a global shift for structure and misses α. Multiscale is what turns detection into *resolution*.

3. ON-SLIDE TEXT

different bands reorganize at different SCALES — no single metric catches them all

  • raw fires broad (every band positive, even θ) — coph suppresses θ:  detect ≠ discriminate
  • β  → pairwise · spectral · hierarchy agree     a strong low-order convergence (robust)
  • α  → the multiscale hierarchy ONLY             a mesoscale change no single-scale read reaches
  • γl → geodesic + Grassmann only                 global geometry, not a relational reorganization
  • a relational lens — tree RANKS, not distances
  • global drift cancels → resolves α, ignores γ_low
  • only the multiscale read catches α, and the global γ_low doesn't fool it

Grassmann distance:

  L = D - W,\qquad U_X=[\,u_1,\dots,u_k\,]\ \ (k\ \text{slowest modes of phase }X)

  d_G(X,Y)=\sqrt{\,k-\textstyle\sum_i \sigma_i^{2}\,}=\sqrt{\textstyle\sum_i \sin^{2}\theta_i},\qquad
  \sigma_i=\cos\theta_i=\mathrm{svd}\!\left(U_X^{\top}U_Y\right)

  one FIXED scale k — not all scales (e^{-τL})

---

4. SPEECH  (~110 s)

So we have a lasting trace, in alpha and beta — but everything we've said about it, that it's band-selective and multiscale, we read through one lens: our hierarchy. The fair question is whether that lens is doing any real work. Is this genuinely a *multiscale* reorganization — or would any ordinary network measure have found the same thing? Detection is easy; the real question is discrimination. So on the same data, against the same matched-strength null, we put three lenses side by side [top radials] — raw pairwise edges, spectral clustering (the Grassmann distance between the two graphs' eigen-subspaces), and our multiscale hierarchy — with the whole standard toolkit behind them [the matrix].

Beta first — it's the easy one, and it makes the opposite point too. Raw pairwise connectivity, the left radial, catches beta — but look how broadly it fires. The stem plot below says it in magnitude: raw is a broad pedestal, positive on every band, even theta; our cophenetic read is a sharp alpha–beta peak, and theta actually drops *negative*. That's detection versus discrimination in a single band — raw's pedestal leans positive on theta, the hierarchy suppresses it. Raw's breadth isn't sensitivity, it's non-selectivity: beta is a genuine low-order convergence everyone agrees on, but raw lights up bands the hierarchy correctly zeroes. Raw detects; it doesn't discriminate.

Alpha is the real test, and the spine of the slide. In the matrix, read the alpha column. Clustering, shortest paths, effective resistance — all blind. Raw does flag alpha, but only fragilely: it's cohort-significant, then it collapses the moment you leave one patient out — it hinges on a single subject. Only the bottom row, the multiscale hierarchy, resolves alpha *robustly*, surviving leave-one-out. Alpha is a mesoscale reorganization no single-scale read holds onto.

The field's sophisticated tool misses it too — that spectral-clustering measure. It runs on the same Laplacian we build our tree on, but it freezes one subspace dimension, a single global scale, instead of diffusing across all of them. It fires instead on low-gamma — lit by exactly the two global methods, Grassmann and the graph geodesic, and nothing else. Low-gamma is a whole-graph shift, not a multiscale reorganization; single-scale spectral clustering mistakes that global change for structure.

So why doesn't our all-scales tree flag that global low-gamma too? Because we compare phases by the tree's *ranks* — who merges before whom — not by distances, and ranks ignore any across-the-board drift. The hierarchy isn't a superset that sees everything plus more; it's a relational *lens* — it resolves the relational alpha and isn't fooled by the global low-gamma. Not "only we see a trace" — the edges see beta, spectral clustering sees beta — but we *resolve* the reorganization the standard tools miss and the field standard mis-reads.

**Careful:**
- **★ SETTLEMENT — what does raw FC fire? (recomputed fresh 2026-07-16 from `controls_ladder{,_apples}/per_cell.csv`, T_test).** The "raw → β only" (Fig B/C) vs "raw → 4/6" (Fig A) clash is NOT a contradiction — it is two graphs × two gates. **Canonical (same mst@0.20 backbone as the cophenetic read, matched-strength cohort-p + LOO):** raw = **β robust (FILL, p=.024, LOO .049)**, **α + γ_high fragile RINGS (p=.042, LOO .082 — cohort-significant, fail leave-one-out)**, δ/θ/γ_low null. So "raw → β only" means β is raw's only *LOO-robust* carrier; α/γ_high are fragile rings, NOT white. **Fig A's "4/6 (δ,α,β,γ_low)"** is a DIFFERENT object: DENSE raw + cohort-p only (no LOO) — the non-selectivity *pedestal*, and note dense raw fires γ_low+δ (null on the backbone) and NOT γ_high. NEVER put the two counts side by side as if comparable. For the α spine: raw is NOT "blind to α" — raw flags α but fragilely (single-patient); only the hierarchy resolves α *robustly*. Say "robustly only in the hierarchy", never "raw misses α".
- Matched-strength is the SOLE null (drift RETIRED). Fig B/C gate = cohort-p + **LOO** (the ≥3/10 count does not decide any raw cell here — LOO does; keep the honest gate = p+LOO, house rule counts-never-the-gate). Fig A = naive cohort-p, dense, magnitude view. State which construction+gate when quoting any count.
- The γ_low reading is the POINT: γ_low is flagged only by the two GLOBAL methods (geodesic + Grassmann) and never by the multiscale read ⇒ it is a global/whole-graph shift, NOT a multiscale reorganization. Do NOT soften this to "Grassmann sees an extra band" — the interpretation is that spectral clustering mistakes a global change for structure.
- Do NOT say "only multiscale sees a trace" / "invisible to simple methods" — raw DETECTS β and so does spectral clustering. raw = NON-selective but COMPLEMENTARY, never "blind/fragile."
- **θ is the detect≠discriminate emblem (internal control).** Frame it as raw's NON-SELECTIVITY: raw's pedestal is positive on *every* band incl. θ, while the cophenetic read pushes θ to zero / negative (Fig A boxes θ; coph θ is the negative bar). This is a statement about the MEASURES' selectivity — raw fires broad, coph suppresses θ — NOT an ontological "θ has no trace" claim (banned: no-band-is-trace-free). Say "coph suppresses θ" / "raw leans positive even on θ", never "θ is trace-free". Accuracy: raw θ does NOT clear the strict gate (radial p≈.05, misses Fig A's naive gate) — the load-bearing fact is the *sign contrast* (raw θ positive pedestal vs coph θ negative), so lead with the pedestal-vs-peak sign, not a "raw significantly fires on θ" overclaim.
- **Grassmann is WHOLE-GRAPH (dense imcoh_abs FC), NOT the mst@0.20 backbone**; its per-patient T_G aggregates over the subspace-dimension grid k (no single "s"). Honest claim: "same Laplacian operator, different construction." ⇒ **the γ_low "misled" contrast is partly CONFOUNDED by graph construction** (dense Grassmann vs backbone hierarchy), so it is NOT airtight. The unconfounded, same-backbone result is the spine: **α is resolved ONLY by the hierarchy** while raw / clustering / geodesic / resistance (all on the same mst@0.20 backbone) miss it. Lead with α; keep γ_low secondary.
- **Why the hierarchy doesn't flag the global γ_low is NOT "it misses that scale."** The cross-phase trace is ρ_sym = rank (Spearman) correlation of cophenetic distances → invariant to across-the-board global drift by construction. Frame as a **relational lens**, never as "multiscale spans all scales so it sees everything the global metrics see" (that is false — it would then flag γ_low).
- **Grassmann does NOT confirm the α+β carrier set.** Its matched-strength cluster-extent verdict (+ LOO): β strong & LOO-stable, γ_low strong, δ significant-but-LOO-fragile, θ/α/γ_high none. It AGREES on β, MISSES α, FLAGS γ_low. Never "Grassmann confirms our carriers."
- α is cophenetic-EXCLUSIVE, β a CONVERGENCE band. Never "coph → β only."
- Marker grammar: band colour = clears (bold/fill = representative, ring = sig-but-LOO-fragile), **white = no trace**. Only genuine null bands are white.
- Do NOT name β's anatomy here (β is DELOCALIZED; OFC is the *encoding* anchor only). → later slide.
- Higher-order β (coph-beyond-raw) is RETIRED (BH-marginal, q≈.064; clean cell was δ, not β).

---

5. FIGURES  (talk PNGs, transparent white-ink, dark-slide ready, in data/outputs/figures/talk/)

- **Fig B — the full read-out ladder (CENTERPIECE).** `fig_grassmann_contrast.png` (gen `scripts/07_figures/talk_fig_grassmann_contrast.py`; standard + cophenetic verdicts from `controls_ladder_apples/per_cell.csv` via representativeness gate, Grassmann from `data/audit/grassmann_cluster_extent/cohort_summary.csv`). Six rows in three tiers (dashed separators): **standard toolkit** (raw FC · clustering · graph geodesic · effective resistance) → **spectral clustering** (Grassmann chordal · whole graph) → **multiscale** (cophenetic · ours). × 6 bands. FILL = robust · band-colour RING = sig-but-LOO-fragile · white = no trace. Two highlighted columns: **β** (convergence carrier) and **γ_low** ("global — not multiscale"; lit only at geodesic + Grassmann). α callout "multiscale-only" (only cophenetic). Grassmann FORMULA in §3 on-slide text. This one figure carries the whole slide.
- **Fig A — pedestal vs peak (magnitude intuition, optional).** `fig_pedestal_vs_peak.png` (gen `talk_fig_pedestal_vs_peak.py`; `controls_ladder/cohort_gate.csv` — ⚠ **DENSE** raw, cohort-p ONLY, no LOO — a DIFFERENT graph and a WEAKER gate than the backbone matrix Fig B). raw = broad positive pedestal (all 6 positive, clears 4/6: δ .024, α .019, β .024, γ_low .032; naive gate), cophenetic = sharp α/β peak (α .007, β .001) with θ negative. θ boxed as internal control. ⚠ **DO NOT quote the "4/6" as a carrier count that competes with Fig B's "raw → β robust."** They are not the same measurement (dense/no-LOO here vs backbone/+LOO there), and the dense raw even fires **γ_low + δ** (null on the backbone). Fig A is ONLY the *non-selectivity pedestal* (raw fires broad → magnitude), NOT a gated carrier verdict. Best used as a warm-up "detect ≠ discriminate" magnitude glance BEFORE the ladder — OR dropped, since Fig B alone carries the gated story and mixing the two is exactly what makes the raw count look self-contradictory. See the SETTLEMENT bullet in §4.
- **Fig C — per-patient blooms (spread, optional).** `fig_descriptor_bloom_contrast.png` (gen `talk_fig_descriptor_bloom.py`). Three radial blooms — raw · Grassmann · cophenetic — ten patient petals per band, wedge glowing where the band clears its cohort gate. Grassmann panel on its OWN radial scale (whole-graph T_G, mean over k). Shows the spread behind the ladder's dots.
- **RETIRED:** `fig_coph_beyond_raw.png` (higher-order β — BH-marginal q≈.064). Do NOT resurrect.

⚠ COMPOSITION NOTE: Fig B is now the standalone centerpiece — it carries "toolkit blind + spectral clustering misled + multiscale resolves." Fig A (magnitude) and Fig C (per-patient spread) are optional supports; the ladder subsumes their "raw non-selective / we tested many metrics" point. If the slide is crowded, run **Fig B alone**, or **Fig B + one** of A/C. Decide on deck.

6. REFERENCES
- Grassmann chordal distance / principal angles between subspaces — standard subspace geometry (spectral clustering / spectral embedding). Graph geodesic (shortest-path) and effective resistance (Laplacian pseudo-inverse; Klein & Randić 1993) are the global path/spectral foils. Controls-ladder + matched-strength null (`controls_ladder_apples`) and the Grassmann matched-strength surrogate + cluster-extent gate (audit_66 / `grassmann_cluster_extent`) are ours.

7. CANVA STATUS
REWRITTEN 2026-07-14 to the correct message. The empirical "the trace is MULTISCALE" slide: the standard network toolkit can't resolve it; the field-standard spectral-clustering measure is misled (flags γ_low = a global shift, not a multiscale reorganization — tell: only geodesic, the other global metric, agrees; and it misses α); only the multiscale hierarchy names α+β. Centerpiece = `fig_grassmann_contrast` (full 7-row ladder). Grassmann formula in on-slide text. Canva job: (1) place Fig B (± Fig A/C, see COMPOSITION NOTE); (2) paste §3 on-slide text incl. the Grassmann LaTeX; (3) speech → presenter notes. No β→OFC; matched-strength only. What s means + localization → later.
