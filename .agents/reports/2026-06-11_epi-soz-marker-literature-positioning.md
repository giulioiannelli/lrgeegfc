---
name: epi-soz-marker-literature-positioning
era: IMCOH_ABS_COHORT_N10
date: 2026-06-11
status: current
kind: reference
scope: Literature positioning for the LRG-propagator SOZ marker — where a heat-kernel diffusion-community marker on |ImCoh| sits in the interictal-iEEG SOZ-localization state of the art. Deep-research-backed, cited, for the manuscript related-work / discussion.
companion_of: .agents/reports/2026-06-08_epi-propagator-soz-marker.md
provenance: deep-research workflow 2026-06-11 (18 primary sources fetched, 87 claims → 25 adversarially verified, 24 confirmed / 1 killed)
---

# Where our SOZ marker sits in the literature

**Head.** Our LRG heat-kernel propagator marker on |ImCoh| connectivity sits in a
**crowded, well-validated** field of interictal-iEEG SOZ node-markers, but occupies a
**genuinely distinct niche**: the multiscale diffusion-community read-off of the
propagator ρ(τ)=e^{−τL̂}/Z has, on the evidence gathered, **not previously been
applied to SOZ localization**. Our LOPO AUC **0.62–0.73** (n=10) is **squarely
within the field-standard range (≈0.70–0.86)** for interictal node-markers — and we
reach it with a **strength-orthogonal** feature where node strength is near chance.
**Honest framing note:** this leave-one-patient-out holds out the classifier
*coefficients* only — the test patient's features still use that patient's own SOZ
labels, so it shows the community structure is **cross-patient consistent**, not that we
can predict SOZ in a *label-free* new patient (that fails, audit_94). The genuine
predictive evidence is the within-patient held-out masking recovery + few-seed curve.
The honest gap vs the benchmarks: **they validate against resection +
Engel/ILAE outcome; we do not have resection ground truth**, so we position as a
**triage/shortlist** tool, not an outcome-validated predictor.

---

## TOC (plain-English)
1. The field has two marker families — directed "inward/sink" and undirected "hyperconnectivity/centrality".
2. The hyperconnected-vs-isolated debate is real, unresolved, and **metric-driven** — our strength-orthogonal undirected community sits between the poles.
3. |ImCoh| is already field-standard for SOZ mapping — our FC substrate is not exotic.
4. The LRG heat-kernel propagator is exactly e^{−τL̂}/Z — and is **new to SOZ**.
5. Our AUC is in-range; our **honest deficit is ground truth**, our honest edge is the strength null.
6. What to cite where, and the four open questions deep research flagged.

---

## 1. Two validated marker families (the benchmark)

**(a) Directed "inward/sink" markers — the dominant interictal characterization.**
The epileptogenic zone at rest behaves as an inward-influenced *sink* (high inward,
low outward connectivity), **not** an outward-driving hyperconnected hub:
- **Interictal Suppression Hypothesis** — Doss et al., *Brain* 2024;147(9):3009
  (n=81 resting SEEG, partial directed coherence + PCA). The high-inward/low-outward
  motif was **PC1 in 64/81 (79%) of patients** and the strongest SOZ differentiator
  (ANOVA SOZ/PZ/NIZ P≈2.6e-9). No common evidence for the competing driver-hub pattern.
- **Source-Sink Index** — Gunnarsdottir et al., *Brain* 2022;145(11):3901 (n=65).
  EZ nodes are "sinks strongly inhibited by neighbouring sources". Localizes EZ at
  **73%** accuracy; predicts surgical outcome at **AUC 0.86±0.07** (vs HFO 0.72±0.07);
  Ps-values track Engel class; failed surgeries had untreated high-metric regions.
- **Directed ImCoh/PDC** — Narasimhan/Park et al., *Epilepsia* 2020 (PMC7899016, n=25,
  alpha): inward PDC **AUC 0.84**; combined directed+nondirected **AUC 0.88** (0.93 in
  favorable-outcome subgroup). Outward "source" reversal appears only at *ictal* onset.

**(b) Undirected hyperconnectivity / centrality / abnormality markers.**
- **Within-resection-zone hyperconnectivity** — Shah et al., *NeuroImage:Clinical*
  2019;23:101908 (n=27, β 15–25 Hz). RZ–RZ connections stronger than RZ–OUT/OUT–OUT,
  and stronger in good- vs poor-outcome patients — **and this survived a
  spatially-constrained null model**. This is the **closest precedent to our result**:
  internally-hyperconnected EZ that survives a *space* control.
- **Eigenvector-centrality biomarker** — Roy et al., *Front. Netw. Physiol.* 2025
  (PMC12129916, n=65, 6 centers). EZ localization **AUC 0.70±0.01 under 26-fold
  leave-one-patient-out**; outcome prediction AUC 0.75. Comparable to source-sink.
- **Normative band-power abnormality** — Taylor/Wang et al., *Brain* 2022;145(3):939
  (atlas from 234 participants). Distinguishability of resected-vs-spared (DRS)
  **AUC 0.75, P=0.0003** (62 patients). *Cite for the validation paradigm
  (DRS / resected-vs-spared / normative atlas), not as a connectivity precedent — it
  is a band-power marker.*

**Takeaway for us:** the canonical AUC band is **≈0.70–0.86**, and the canonical
validation is **LOPO + resected-vs-spared + outcome-stratified**. Our LOPO AUC
0.62–0.73 (n=10) lands inside that band — its lower edge.

## 2. The hyperconnected-vs-isolated debate is metric-driven

The debate our own framing flagged is genuinely unresolved, and the disagreement
**traces to methodology**: directed studies (ISH, source-sink) find *inward*
hyperconnectivity/sink; undirected studies (Shah 2019; Goodale/Narasimhan 2020,
PMC7225010, alpha ImCoh — internal **and** external hyperconnectivity, P<.01 corrected)
find hyperconnectivity. Doss explicitly attributes the discrepancy with high-outward
studies to **measuring strength vs counting significant connections**. So "hyperconnected"
and "inward-sink" are partly reconcilable (inward hyperconnectivity is still
hyperconnectivity, just asymmetric).

**Our position:** a **strength-orthogonal, undirected, path-summed diffusion community**
sits *between* the poles. It aligns most with Shah's internally-hyperconnected-and-
space-robust view, but is built on diffusion *affinity* (all-paths heat-kernel), not raw
edge strength — which is exactly why it can dissociate from the strength confound that
drives the debate. **Frame as sitting between the poles, not as adjudicating them.**

## 3. |ImCoh| is field-standard for SOZ

- **Nolte et al. 2004**, *Clin Neurophysiol* 115:2292 (>3000 cites): the coherency of
  volume-conducted/non-interacting sources is necessarily real, so the **imaginary part
  is immune to zero-lag volume conduction**.
- **Park et al. 2022**, *Neurology* (DOI 10.1212/WNL.0000000000200386, n=32 resting
  SEEG): between-region **ImCoh higher in EZ** than non-EZ every day (p=0.007–0.023).
- **Goodale/Narasimhan et al. 2020**, *Neurosurgery* (PMC7225010, n=15): alpha-band
  imaginary coherence chosen to "ignore zero-time-lag signals and minimize volume
  conduction"; FC model predicts epileptogenicity at **AUC 0.78**.

**Caveat to honor:** ImCoh is immune only to *instantaneous/zero-lag* mixing — it does
not remove all spurious connectivity and it discards genuine zero-lag physiological
coupling. The "100% of patients" MEG figure (Elisevich/Xiang 2014) is within-patient
*sensitivity* only, no controls. State this; do not overclaim ImCoh as a cure-all.

## 4. The LRG propagator is e^{−τL̂}/Z — and is new to SOZ

Villegas et al., *Nature Physics* 2023;19:445 (arXiv:2203.07230); open-access review
arXiv:2406.02337. The network propagator **K̂ = e^{−τL̂}** is "the discrete counterpart
of the path-integral formulation of diffusion"; the density operator **ρ̂(τ) = K̂/Tr(K̂)
= e^{−τL̂}/Z**, L̂ = D̂ − Â. Each element **K̂_ij sums diffusion trajectories over all
paths** between i and j at time τ — the formal justification for reading node-pair
affinity and diffusion communities off the propagator (consistent with spectral-graph
theory: exp(−τL) admits the Poisson random-walk-mixture expansion Σ_k(τ^k e^{−τ}/k!)P^k).

The nearest prior graph-spectral node-marker in spirit is **communicability** exp(A)
(a walk-counting affinity). But **the heat-kernel / LRG multiscale-diffusion formulation
applied to SOZ localization is, per this corpus, unprecedented.** Our novelty is the
**multiscale band-specific diffusion community** (segregated block surviving strength +
space nulls) with a **strength-orthogonal LOPO marker**.

⚠ **Novelty is an absence-of-evidence claim** (caveat, deep research): the corpus has
the LRG definition and many SOZ markers but no LRG-for-SOZ. Before asserting outright
priority in the manuscript, run a **dedicated negative search** for
"communicability / heat-kernel / random-walk diffusion + SOZ/epilepsy localization",
since communicability node-markers and diffusion parcellation exist in adjacent
brain-network literature.

## 5. Our edge and our deficit, stated honestly

- **Edge:** the field has repeatedly had to *control for* degree/strength (Doss'
  strength-vs-count discrepancy; Shah's spatial null). Our marker is **strength-
  orthogonal by construction** (matched-strength-z features; node strength at chance
  cross-patient) **and** space-irreducible (C5 spatial null). That is the cleanest
  possible answer to the confound the field keeps fighting.
- **Deficit:** every benchmark validates against **resection margin + Engel/ILAE
  outcome**. We have **clinically-labelled SOZ, not resection/outcome ground truth**.
  So: (i) we cannot claim outcome-meaningfulness; (ii) we position as a **seed-based
  triage/shortlist** (precision@k, enrichment/lift) — which the field has **no canonical
  semi-supervised benchmark to compare against** (open question 2), so our few-seed
  curve is reported on its own honest terms (δ/low-γ/β 3–5× enrichment, ~30% precision).

## 6. Citation map + open questions

**Cite where:**
- *Related work — node markers / benchmark AUC band:* Gunnarsdottir 2022 (source-sink,
  AUC 0.86), Roy 2025 (EVC, LOPO AUC 0.70), Taylor 2022 (DRS paradigm, AUC 0.75).
- *Characterization debate:* Doss 2024 (ISH), Shah 2019 (internal hyperconnectivity +
  spatial null — our closest precedent), Goodale 2020 (ImCoh hyperconnectivity).
- *FC substrate:* Nolte 2004 (ImCoh foundational), Park 2022 + Goodale 2020 (ImCoh→SOZ).
- *Method:* Villegas 2023 + review arXiv:2406.02337 (LRG propagator).

**Open questions deep research flagged (verbatim priorities):**
1. Dedicated negative search for diffusion/communicability+SOZ before claiming priority.
2. No canonical semi-supervised "mark-few-discover-many" SOZ benchmark surfaced — our
   triage framing may lack a direct precedent to cite; confirm or fill the gap.
3. Cross-walk our undirected diffusion community to the directed inward/sink (ISH,
   source-sink) family on the same cohort would strengthen positioning.
4. Open benchmark datasets (Taylor normative cohort; Roy/Gunnarsdottir multi-center)
   could permit external LOPO validation against resection + outcome — the path from
   triage-shortlist to outcome-validated.

**Verification provenance.** Deep-research workflow 2026-06-11: 5 search angles, 18
primary sources fetched, 87 falsifiable claims extracted, 25 adversarially verified
(3-vote, need 2/3 to kill), **24 confirmed / 1 killed**. The killed claim was an
inverse-framing of a real Taylor 2022 result (spared-more-abnormal direction), correctly
dropped. Two findings rest on 2-1 votes (Roy EVC "across folds" vs held-out; Taylor
band-power-not-connectivity) — reproduce those nuances precisely when citing. AUC
landscape current as of 2024–2025; re-check before submission.
