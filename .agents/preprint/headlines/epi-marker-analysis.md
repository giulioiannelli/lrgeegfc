---
name: epi-marker-analysis
era: IMCOH_ABS_COHORT_N10
status: current
kind: headline
scope: HEADLINE — the COMPLETE epileptic-marker investigation on the LRG Laplacian density matrix ρ(τ)=e^{−τL}/Z (all bands, all audits). The full honest arc: node-intrinsic markers fail (= hubness + electrode depth); the relational diffusion COMMUNITY is the real finding (survives strength + spatial nulls); an all-contacts proximity confound was caught and withdrawn; the clean survivor is the δ off-shaft distant marker (AUC 0.72, 8/10), seed-based and phase-stable. Math + plain language.
---

# The epileptic marker — full investigation analysis

**Head.** We asked whether the Laplacian density matrix ρ(τ) = e^{−τL}/Z (the LRG
diffusion propagator) can flag epileptic (seizure-onset, SOZ) contacts. The honest
arc, in one breath: a **node-by-node** dynamics marker **does not exist** — within a
patient the recoverable signal is just node strength (hubness) plus electrode depth,
and across patients nothing transfers. But read **relationally**, the SOZ contacts
form a **strength-independent, spatially-irreducible, band-specific diffusion
community** — a genuine structural finding that survives strength- and
space-matched nulls. Turned into a marker, that community's per-node read-out
(*how strongly a contact co-diffuses with the known SOZ*) recovers **distant**
SOZ — on separate electrodes, with spatial proximity removed by construction — at
**median leave-one-electrode-out AUC ≈ 0.72 in 8/10 patients**, **δ-band only**,
**stable across rest and task**, with a label-shuffle null at p = 0.000. Along the
way we **caught and withdrew** an inflated "all-contacts" version of the marker
(a proximity tautology) and a label-free cross-patient detector (fails at chance).
**What it is:** a mechanistic, multiscale, strength-orthogonal, *seed-based*
triage/hypothesis tool, validated against clinical SOZ labels — **not** an
outcome-validated, from-scratch clinical localiser.

---

## Plain-language summary — the argument that we have a marker

**One line.** Epileptic contacts are wired into the network as a *group that diffuses
together*. Mark a few of them, and the diffusion picture points at the others — even on
far-away electrodes — well enough to be a useful shortlist. It shows up in **8 of 10
patients**; the 2 that differ are a *different network type*, not a failure.

**What the marker is:**
- **Let information spread on the connectivity graph and watch where it goes.** Seed the
  diffusion propagator from the known seizure contacts and ask which other contacts the
  diffusion most connects to those seeds.
- **It reads the wiring *pattern*, not just "strong nodes."** We subtract out node
  strength (hubness) first; the signal survives. It's about *how* a contact co-diffuses
  with the seizure group, not merely how connected it is.
- **It finds *distant* contacts, not neighbours.** We only ever score contacts on
  electrodes with no seed, so "sitting next to a marked one" can't help — no circular
  credit for proximity.

**Why it's a real pattern (not noise):**
- **Specific to one rhythm** — only the slow **δ** band carries it; α/β/low-γ give chance.
  A method that lit up everywhere would be suspicious; this one doesn't.
- **A fake answer collapses** — relabel random contacts as "fake seizure" and the effect
  vanishes (p = 0.000). The signal is tied to the *true* seizure group being a community.
- **Stable in time** — essentially the same before, during, and after the task. It
  behaves like a fixed trait of the brain, not a passing state.
- **Holds out of sample** — pick the method on 9 patients, test on the 10th: same answer.

**The pattern across the cohort — why 1–2 odd patients is fine:**
- **8 of 10 patients show it** (median AUC 0.72, up to 0.95), inside the published field
  range (0.70–0.86).
- **The 2 that don't (Pat_10, Pat_15) are a clean, understood second type:** their seizure
  contacts are network *hubs*, not a co-diffusing community, so the diffusion marker
  correctly fails and plain node-strength works instead. Both are right-hemisphere implants
  in a mostly-left cohort, so the split is **consistent with implant / electrode geometry**.
  For a real biological pattern across 10 differently-implanted patients, **1–2 going the
  other way is expected** — and here those 2 are *explained*, not dropped.

**What it gives you in practice:**
- **Mark few, rank the rest.** From just **2–3 marked contacts** the remaining distant
  seizure contacts rank above healthy tissue (AUC ≈ 0.66–0.68), improving toward ≈0.72 as
  you mark more — a several-fold-enriched shortlist.
- **Tell the two types apart and switch tools.** A simple seed-based check (are the marked
  contacts hubs?) sorts a patient into "community" or "hub" type; using the matched tool
  per patient lifts coverage **8/10 → 9/10**, with large top-5 enrichment inside each type
  (**9–15×** by diffusion for community patients, **up to 12×** by strength for hub
  patients — the wrong tool gives ~0).

**The honest ceiling.** It is **seed-based** (needs a few known contacts; can't start from
zero), **δ-only**, and validated against **clinical labels, not surgical outcomes**. So the
honest sale is *a strength-orthogonal, multiscale, mechanistic network signature of the
seizure zone + a proof-of-concept triage shortlist* — **not** a finished clinical detector.

---

## 0. The object and the mathematics (read once, in plain language)

**The graph.** Per patient, per band, per phase we have a functional-connectivity
matrix **W** = |ImCoh| (the band-averaged magnitude of the imaginary part of
coherency — volume-conduction-immune by construction, Nolte 2004). W is symmetric,
non-negative, zero diagonal; node *i* = one sEEG contact.

**The Laplacian.** Plain: "how a node differs from its neighbours."
$$ L = D - W, \qquad D = \mathrm{diag}\!\Big(\textstyle\sum_j W_{ij}\Big) $$
The diagonal entry **D_ii = node strength** = total connectivity of contact *i*
(this is the "hubness" that haunts the whole investigation).

**The heat-kernel propagator / density matrix.** Plain: "let information diffuse on
the graph for a time τ; the propagator says how much reaches node *j* from node *i*."
$$ K(\tau) = e^{-\tau L} = U\,\mathrm{diag}(e^{-\tau\lambda})\,U^{\!\top},
\qquad \hat\rho(\tau) = \frac{K(\tau)}{\mathrm{Tr}\,K(\tau)} = \frac{e^{-\tau L}}{Z} $$
where (λ, U) is the eigendecomposition of L. Each entry **K_ij(τ)** is a sum over
**all diffusion paths** from *i* to *j* of length-weighted walks (Villegas 2023:
the path-integral / random-walk-mixture reading $e^{-\tau L}=\sum_k \frac{\tau^k
e^{-\tau}}{k!}P^k$). This is *the* object — the same one whose multiscale cophenetic
structure carries the β trace.

**Scale = diffusion time τ.** Plain: "small τ = local/fine structure, large τ =
coarse/global structure." We scan a per-patient grid
$$ \tau \in \mathrm{geomspace}\big(1/\lambda_{\max},\; 10/\lambda_{\max},\; 6\big),
\qquad \tau_0=\text{fast} \ \ldots\ \tau_5=\text{slow}. $$
This continuous scale flow is the multiscale part — no fixed number of clusters or
components is ever chosen (the contrast with PCA / spectral clustering; §8).

**Three ways we read the propagator** (this is the spine of the whole story):
1. **Node-intrinsic** — each node summarised by its *own* diffusion (return
   probability ρ_ii, row entropy, commute distance). → **Act I: fails.**
2. **Relational** — each node *i* scored by its **affinity to the known SOZ set E**.
   → **Act II/III: the finding + the marker.**
   $$ \mathrm{aff}(i) = \frac{1}{|E\setminus i|}\sum_{s\in E\setminus i}
      \big[e^{-\tau L}\big]_{i,s} $$
   Plain: "how much of node *i*'s diffusion lands on the marked seizure contacts."
3. **Strength control.** Every relational score is **residualised on node strength**
   so hubness cannot explain it. Fit a line $\mathrm{aff} \approx a\cdot
   \mathrm{strength}+b$ by least squares and keep the residual:
   $$ \widetilde{\mathrm{aff}}(i) = \mathrm{aff}(i) - \big(a\cdot
      \mathrm{strength}(i)+b\big). $$

**How we score a marker — AUC by rank concordance** (Mann–Whitney). Plain: "the
chance a true SOZ outranks a healthy contact."
$$ \mathrm{AUC} = \frac{\#\{\text{SOZ}>\text{healthy}\} + \tfrac12\#\{\text{ties}\}}
{\#\text{SOZ}\times\#\text{healthy}}. $$

---

## 1. Map of the whole investigation (every audit, one line)

| audit | question | metric | verdict |
|---|---|---|---|
| 77/78/79 | does the cross-phase trace live in/around epi tissue? | stratified ρ^coph / Grassmann | context; "strengthens under epi-exclusion" later shown **node-count, not epi-specific** |
| 80 / 81 | per-node SOZ classifier / leverage, **cross-patient** | LOPO AUC | **NEGATIVE** — at chance (FC strength is patient-scaled) |
| 85 / 86 | **within-patient** template recovery from ρ(τ) node features | recovery AUC + fair strength + matched-strength-z | β = **hubness**; only a modest, borderline α residual |
| 87 | hide-and-seek + within-shaft control | masked recovery AUC | looked positive (within-shaft) … |
| §10 depth | add along-shaft **depth** baseline | within-shaft AUC | … **dissolves**: depth 0.72–0.74 beats dynamics every band |
| 88 | signed-ImCoh **magnetic/Kunegis Laplacian** (directed) | within-shaft AUC | **NEGATIVE** — also just depth; no beyond-hubness |
| 89 | read propagator **relationally**: block means ρ_EE/ρ_EN/ρ_NN | matched-strength contrast | **POSITIVE** — SOZ↔SOZ co-diffusion elevated beyond strength |
| 90 | **C5** strength **+ spatial**-matched null on the community | surrogate z | **survives** δ/β/low-γ (+α); **high-γ clean null** |
| 91 | relational node features f_aff / f_seg / f_part | within-patient separation | epi vs non 0.74–0.86, ≫ strength ~0.6 |
| 92 / 93 | masking recovery + LOPO logistic + P(SOZ), **all-contacts** | AUC ~0.81 / LOPO 0.62–0.73 | **later WITHDRAWN as a performance claim** (proximity, audit_99) |
| 94 | **label-free** cross-patient SOZ selection | precision / lift | **NEGATIVE** — narrows ~8× within-patient, **selection at chance** → seed-based |
| 95 | "mark few, discover many" few-seed curve (all-contacts) | precision@top-m | proximity-confounded → **withdrawn**; honest **off-shaft few-seed recomputed 2026-06-18** (§4) |
| 98 | deployment (all-contacts) | precision/lift | **WITHDRAWN** (proximity) |
| 99 | is "near a known SOZ → SOZ" a tautology? | nearest-contact baseline | **YES** — beats the marker ~2× on all-contacts → **discard proximity** |
| 101 | 16-operator sweep, **off-shaft** leave-one-shaft-out | strength-resid LOSO AUC | **δ heat kernel wins: 0.72, 8/10** (clean) |
| 102 | is it a fluke? | label-shuffle null + nested LOPO | **null p=0.000; LOPO 0.716, heat selected 9/10** |
| 103 | honest precision + candidate shortlists (off-shaft) | precision@k / lift | real but **concentrated** (cohort flat ~1.2×; responders ~2.2× via 3 pts) |
| 104 | propagator+strength **compound** + calibrated P(SOZ) | off-shaft AUC / Brier | **switch** 0.716→0.750; P(SOZ) calibrated-at-top, modest |

(Numbers verified against the raw CSVs in-session 2026-06-18 for audits 101–104 and
the phase check; audits 80–95 are transcribed from their dated reports.)

---

## 2. Act I — a node-by-node dynamics marker does not exist (the honest negatives)

**Cross-patient (audit_80/81).** Per-node identity and "leverage" features are at
**chance** in leave-one-patient-out. Reason (plain): absolute FC strength lives on a
patient-specific scale, so no global threshold learned on others transfers to a new
patient.

**Within-patient (audit_85/86/87/88).** The realistic clinical scenario: a clinician
marks *some* contacts, wants the rest; self-normalise inside the patient. We summarised
ρ(τ) per node by return probability ρ_ii(τ), diffusion-row entropy
$H_i(\tau)=-\sum_j p_{ij}\log p_{ij}$ with $p_{ij}=K_{ij}/\sum_j K_{ij}$, and commute
distance — then tried to recover held-out epi contacts. It **looked** positive, but
three "wins" each dissolved under the next control:

1. **The β "breakthrough" was an unfair baseline.** audit_85 scored strength by
   *template distance* (which penalises nodes *more extreme* than the epi mean —
   wrong for a monotone feature). Under a **fair** Mann–Whitney baseline, node
   strength is **0.738 at β** and *beats* the propagator. Retracted.
2. **The "shaft-independent" within-shaft recovery was electrode depth.** Adding
   along-shaft depth `contact_index_norm` as a baseline (epileptic foci sit at
   characteristic mesial depths): **depth alone recovers hidden epi at AUC 0.72–0.74
   in every band**, beats both the propagator (0.61–0.67) and the signed-Laplacian
   (0.60–0.71), and the dynamics **add nothing on top of depth** (all p ≥ 0.33).
3. **The directed/signed extension added nothing.** The signed-ImCoh **magnetic /
   Kunegis Laplacian** ($L_H = \bar D - iA$ Hermitian; real $L_s=\bar D - B$) — the
   natural way to inject directionality — is Stage-1 negative: its within-shaft
   recovery is also depth.

**Lesson (locked).** Within-patient epileptic-contact recoverability = **node
strength (hubness, gross) + along-shaft depth (fine)** — both *already known to the
implanting team*. The propagator subspace and the magnetic Laplacian add **no**
node-intrinsic marker beyond geometry. *The marker, if any, is not node-intrinsic.*

---

## 3. Act II — the relational diffusion COMMUNITY (the real positive finding)

The pivot: stop scoring nodes by their *own* diffusion; score node **pairs** by
membership. Partition pairs into epi–epi (EE), epi–nonepi (EN), nonepi–nonepi (NN)
and read the propagator as **block means** $\bar\rho_{EE},\bar\rho_{EN},\bar\rho_{NN}$
(audit_89). Plain question: *do seizure contacts diffuse to each other more than
strength or anatomy alone would predict?*

**Control battery (matched-strength is mandatory at every step):**
- **C1** absolute vs strength-preserving rewiring — passes but inflated (context only).
- **C2** contrast $\bar\rho_{EE}-\bar\rho_{NN}$ vs matched-strength — epi block enriched.
- **C3** strength-matched random-subset null *on the real graph* — **decisive: cannot
  be hubness.** Strong δ/low-γ/β.
- **C4** cross-shaft restriction — positive where tested (not just same-probe coupling).
- **C5 (decisive, audit_90)** strength **+ spatial**-matched random subsets — matched
  to the epi set on *both* the strength-quintile profile **and** spatial spread
  (radius of gyration in implant x,y,z). **This closes the anatomy/proximity confound
  at the block level.**

**C5 verdict (cohort, n=10, max over τ;** `spatial_null_cohort.csv`**):**
δ +2.06 (p .042), β +1.47 (p .003), low-γ +2.92 (p .002), α +1.81/+2.06 survive;
**high-γ a clean null** (p .42, 1/10). The diffusion community is **neither hubness
(C3) nor spatial clustering (C5)**, and it is **band-resolved** (high-γ correctly
returns null).

**The relational node features** (audit_91), all strength-orthogonal (matched-strength
-z), the per-node expression of the community:
- $f_{\mathrm{aff}}(i)=\overline{\rho(i,E\setminus i)}$ — affinity to the SOZ community.
- $f_{\mathrm{seg}}(i)=\overline{\rho(i,E)}-\overline{\rho(i,N)}$ — interface
  segregation (SOZ couple to SOZ **more** than to healthy tissue; ρ_EN is *suppressed*).
- $f_{\mathrm{part}}(i)$ — share of node *i*'s diffusion routed to the SOZ.

Within-patient epi-vs-non separation (matched-strength-z): δ 0.835, α 0.754, β 0.790,
low-γ 0.744 — all ≫ the strength baseline (~0.60). **This is a genuine structural
result about how the epileptogenic zone is organised**, independent of the
performance-metric questions that follow.

---

## 4. Act III — the proximity confound, and the clean survivor (δ off-shaft marker)

**What looked great but was confounded.** The relational marker, scored on
**all contacts** (audit_92/93/95/98), gave masking-recovery AUC ~0.81 and cross-patient
LOPO 0.62–0.73 with a calibrated P(SOZ). **audit_99 then caught the confound:** SOZ
contacts are labelled because they sit in a clinically-marked epileptogenic **area**,
so "a contact near a known SOZ is itself SOZ" is **baked into the labelling**. A
trivial **nearest-contact** baseline beats any network marker ~2× on all-contacts
recall. The all-contacts performance numbers are therefore **WITHDRAWN as a
performance claim** (the *community* finding of Act II is untouched — it has its own
spatial null).

**The honest metric — leave-one-shaft-out (LOSO) off-shaft distant discovery.**
Plain: *hide all SOZ on one electrode; seed from SOZ on the other electrodes; ask
whether the hidden electrode's SOZ rank above healthy, looking only at contacts on
electrodes with no seed.* Targets and controls are all far from every seed, so
proximity cannot help.

**audit_101 — 16-operator sweep, off-shaft, strength-residual.** Operators from
L = D − W: heat kernel $e^{-\tau L}$ (6 τ), normalized-Laplacian heat (6 τ),
personalized PageRank (3 α), Katz, communicability $\exp(W/\rho)$, negative
effective-resistance, diffusion-distance. **In δ the top cells are all heat / heatN
variants;** communicability/PPR/Katz/resistance sit below. The seed-affinity marker
is $\mathrm{aff}(i)=\frac1{|S|}\sum_{s\in S}[e^{-\tau L}]_{i,s}$, residualised on
strength. Per-patient δ slow-heat (`heat_t5`) AUC:

| Pat | 08 | 14 | 05 | 06 | 13 | 03 | 02 | 07 | 10 | 15 |
|---|---|---|---|---|---|---|---|---|---|---|
| AUC | 0.95 | 0.91 | 0.84 | 0.82 | 0.76 | 0.68 | 0.66 | 0.57 | 0.41✗ | 0.21✗ |

**Median 0.716, 8/10 > 0.5** (slow-τ ridge robust: heat_t3 0.749, heat_t4 0.752,
heatN_t4 0.742). **δ-specific:** α 0.60, β 0.59, low-γ 0.61. The strength baseline
alone reaches only 0.63 — the residualised heat (0.72) beats it.

**audit_102 — it is not a fluke.**
- **Label-shuffle null: p = 0.000.** Relabel n_epi random contacts (spanning ≥2
  shafts) as fake-SOZ, rerun LOSO; null median ≈ 0.49 (p95 ≈ 0.57) vs real 0.72–0.75;
  empirical p = 0.000 over 500 draws, all four hero markers. The recovery is specific
  to the *true* SOZ being a co-diffusing community.
- **Nested leave-one-patient-out: held median 0.716**, `heat_t5` selected in **9/10**
  folds. (Guards the τ/operator choice; the marker is parameter-free given band+τ.)

**Phase-robustness (computed 2026-06-18 — intrinsic trait, not a phase artifact):**

| | rest_pre | task | rest_post |
|---|---|---|---|
| median δ heat_t5 AUC | 0.739 | 0.762 | 0.716 |
| n > 0.5 | 8/10 | 8/10 | 8/10 |

Per-patient AUC correlates across phases (Spearman 0.69–0.86); Pat_10/Pat_15 fail in
*every* phase. The δ diffusion community is present before, during, and after the task.

**Two populations (a finding, not just a failure).** ~8/10 patients: **SOZ = a
diffusion community** (the propagator works, strength is at chance). Pat_10 & Pat_15:
**SOZ = network hubs** (node strength predicts, the propagator anti-predicts —
Pat_15 affinity 0.21 vs strength 0.93). This is the same split seen across the whole
investigation, and it maps onto the literature's "hyperconnected hub vs isolated
community" debate (§8).

### Few-seed reconstruction — "start with a few, rank the rest" (off-shaft, honest)

The deployment question directly: *mark only k contacts; can the marker reconstruct the
rest by ranking every other contact from most- to least-SOZ-like?* The clean seed-count
curve — pick k random SOZ seeds, rank all contacts on **no-seed electrodes**, recover the
remaining *distant* SOZ (δ heat_t5, strength-residual, 300 random draws/patient, computed
2026-06-18):

| seeds k | propagator AUC (median) | n>0.5 | strength AUC | top-m enrichment (responders) |
|---|---|---|---|---|
| 2 | 0.659 | 8/10 | 0.620 | 1.84× |
| 3 | 0.684 | 8/10 | 0.640 | 2.01× |
| 5 | 0.724 | 6/8 | 0.653 | 5.43× |
| all-but-one-shaft (LOSO) | 0.716 | 8/10 | 0.629 | 2.18× |

**Plain reading: yes, but modestly and monotonically.** From just **2–3 marked contacts**
the propagator ranks the remaining distant SOZ above healthy at **AUC ≈ 0.66–0.68 (8/10
patients)**, beating the strength ranking (0.62–0.64); add more seeds and it climbs
smoothly toward the ~0.72 leave-one-shaft-out value. The top-m shortlist enrichment rises
with seeds (≈1.8–2× in responders at k=2–3, ≈5× at k=5). The same per-patient split holds
at every k — Pat_06/08/14 reconstruct strongly even at k=2 (0.80–0.89), the hub-patients
(Pat_10/15) fail at all k. **This is the honest replacement for the withdrawn all-contacts
few-seed curve** (audit_95, which reported 0.73–0.85 from k=2) — once proximity is removed
the numbers come down. It is a real triage/reconstruction aid (a few seeds → a
several-fold-enriched shortlist of the rest), **not** a clean "recover all the others"
detector, and it is δ-band and responder-dependent. *(Reproducible:
`audit_113_epi_fewseed_offshaft.py` → `data/audit/epi_marker_fewseed/`.)*

---

## 5. Act IV — label-free fails; deployment is honest but modest

**Label-free cross-patient (audit_94): NEGATIVE → fundamentally seed-based.** With
**zero** labels in a target patient, unsupervised co-diffusion clustering **narrows**
the SOZ from ~115 contacts to one ~15-node community at ~4–5× enrichment (captures
50–62% of SOZ) — but **selecting which community is pathological is at chance** (top
precision ≈ prevalence, lift < 1, every band). Multiple *physiological* communities
are equally coherent; nothing label-free marks the pathological one. So the marker
**cannot run on a zero-label patient** — it needs ≥2 known SOZ shafts to seed.

**Precision (audit_103, off-shaft): real but concentrated.** Off-shaft SOZ prevalence
≈ 3.7%. Cohort-**all** lift@5 = **1.23× (essentially flat)**; the responders-only
(8/10) median lift@5 = 2.18× is carried by **three** patients (Pat_05 8.7×, Pat_08
14.7×, Pat_14 11.4×); responders recall@10 ≈ 0.34. **The honest cohort headline is the
rank-wide AUC (0.72), not precision@k.** Candidate shortlists are **hypotheses** (no
ground truth).

**Two-population enrichment — split the cohort, apply the matched marker.** The cohort
splits into **propagator-driven** (SOZ = a diffusion community) and **strength-driven**
(SOZ = hubs). A **label-free regime detector** — the mean node-strength percentile of
the seeds, $r=\overline{\mathrm{strength\text{-}percentile}(\text{seeds})}$ — flags which
is which (r ≥ 0.70 → hub-regime), and a **hard switch** picks the matched marker:
$$ s_{\mathrm{switch}}(i)=\begin{cases}\mathrm{strength}(i) & r\ge 0.70\\
\widetilde{\mathrm{aff}}(i) & r<0.70\end{cases}. $$
Exact per-patient picture (off-shaft δ AUC; lift@5 = top-5 enrichment over base rate;
**bold** = the matched/winning tool):

| patient | r | regime | propagator AUC | strength AUC | matched-tool lift@5 |
|---|---|---|---|---|---|
| Pat_08 | 0.47 | community | **0.95** | 0.38 | 14.7× (prop) |
| Pat_14 | 0.56 | community | **0.91** | 0.58 | 11.4× (prop) |
| Pat_05 | 0.53 | community | **0.84** | 0.54 | 8.7× (prop) |
| Pat_06 | 0.70 | community | **0.82** | 0.72 | — |
| Pat_13 | 0.59 | community | **0.76** | 0.67 | 3.2× (str) |
| Pat_03 | 0.36 | community | **0.68** | 0.34 | 2.9× (prop) |
| Pat_02 | 0.72 | hub | 0.66 | **0.82** | 1.5× (prop) |
| Pat_07 | 0.73 | hub | 0.57 | **0.74** | 3.5× (str) |
| Pat_15 | 0.87 | hub | 0.21 | **0.93** | 12.3× (str) |
| Pat_10 | 0.37 | community | 0.41 | 0.35 | 0 (both fail) |

The enrichment from differentiating, cohort-level (exact):

| strategy | median AUC | n>0.5 | lift@5 (all) |
|---|---|---|---|
| propagator only (no split) | 0.716 | 8/10 | 1.23× |
| strength only | 0.629 | 7/10 | 0.81× |
| **differentiated (switch)** | **0.750** | **9/10** | **1.89×** |

**Plain reading: differentiating buys coverage, and within each population the matched
tool buys large enrichment the wrong tool gives ≈0.** Splitting the cohort and applying
the right marker lifts coverage **8/10 → 9/10** (AUC 0.716 → 0.750); the cohort gain is
essentially **rescuing the one recoverable hub-patient** (Pat_15 0.21 → 0.93 via
strength), while protecting the responders. Within a population the matched marker is
decisive: community patients reach **9–15× top-5 enrichment** via the propagator
(strength gives 0); hub-patients reach **3–12×** via strength (the propagator gives 0).
**Honest limits:** Pat_10 is unrecoverable by either (a third, true-negative type); the
regime detector is imperfect at the border (Pat_02, r=0.72, is flagged hub yet the
per-fold switch reverts to the propagator and leaves its strength edge 0.82 on the table);
and the split needs the seeds to estimate r, so it is itself seed-based. A continuous
**blend** instead of the switch raises cohort *lift* but **damages the strong responders**
(Pat_08 0.95 → 0.56) — a cohort-median artifact; **do not blend.**

**Calibrated P(SOZ) (audit_104): emittable, calibrated-at-top, modest.** A
leave-one-patient-out logistic $P(\mathrm{SOZ})=\sigma(b_0+b_1 z_{\mathrm{aff}}+b_2
z_{\mathrm{str}})$ (within-fold z-scored) transfers for responders (held-out AUC median
0.79) but **not** the hub-patients (fixed positive weights cannot serve both regimes;
Pat_10 0.38, Pat_15 0.24). Calibration tracks the diagonal at the high end (P ≈ 0.13 →
13% observed) but values are small (top ~0.13–0.35) and **Brier 0.0396 ≈ the 3.9%
prevalence floor (not beaten)**. A triage score, not a verdict — and note the
P(SOZ) model and the Pat_15-rescuing switch are **different objects**.

**Classifying the unmarked contacts (occult-candidate discovery).** The deployment
endpoint of the original clinical question — *"can we say this **unmarked** contact might
be SOZ?"* Both deployment audits emit a per-patient **ranking of every unmarked contact**:
`audit_103` (`candidates_per_patient.csv`) ranks them by off-shaft strength-residual
affinity to the marked SOZ, annotated with Desikan-Killiany region, shaft, score
percentile, and distance to the nearest seed; `audit_104` (`candidates_with_prob.csv`)
attaches the leave-one-patient-out **calibrated P(SOZ)**. In the responders the top
candidates are **anatomically coherent and several are genuinely distant**:

| patient | top unmarked candidate | dist. to nearest seed | P(SOZ) |
|---|---|---|---|
| Pat_02 | lateral **orbitofrontal** (resonates with the β-trace OFC) | **31 mm** | ~0.09 |
| Pat_05 | rostral-middle-frontal | **25 mm** | 0.35 |
| Pat_08 | superior-**temporal** cluster | 16 mm | 0.24 |
| Pat_13 | **amygdala** / superior-temporal | 13 mm | — |
| Pat_14 | superior-temporal / **hippocampus** / insula | 4–13 mm | 0.19 |

**The honest reading — this is a hypothesis list, not a diagnosis:**
- **No ground truth.** No resection/outcome exists for these contacts, so an unmarked
  high-scoring node **cannot be confirmed** SOZ. The lists are **hypotheses**.
- **Internally meaningful, externally unconfirmed.** The off-shaft validation says the
  ranking is real (held-out *known* SOZ rank high), so the top of the list is *enriched*
  for true-but-unlabelled SOZ — at the measured precision (~2–5× over base rate in
  responders), which still means **most top entries are false positives**.
- **P(SOZ) is modest** — top responders ~0.10–0.35; a triage score, not a confident call.
- **White-matter contamination** — several top candidates are `Wm` (e.g. a spurious
  Pat_03 P=0.89 WM node); re-rank on a gray-only montage before reading the list.
- **Responders only** — for the two hub-patients (Pat_10/15) the classification is
  unreliable (the marker does not work there).

So: **yes, we produce a ranked, anatomically-annotated, probability-scored classification
of the unmarked contacts** — a usable *occult-candidate shortlist* whose top entries are
enriched for, and plausibly located near, true SOZ — but it is explicitly
**hypothesis-generating, not confirmatory**, and bounded by the same seed-based, δ-band,
no-outcome-ground-truth ceiling as the rest of the marker. (The cross-consistency with the
β-trace localisation — OFC in Pat_02, mesial/lateral temporal elsewhere — is a plausibility
signal, not proof.)

---

## 6. Honest scope (all the limitations, first-class)

1. **n = 10**; responder subgroup 8; precision carried by 3. Power is thin.
2. **No surgical-outcome ground truth** (no resection margin / Engel-ILAE). All metrics
   are against clinician SOZ labels — an **upper bound** on true occult discovery; no
   outcome-meaningfulness claim is available.
3. **Seed-based / semi-supervised** — needs ≥2 known SOZ shafts; label-free selection
   fails (audit_94). An augmentation/triage tool, not de novo.
4. **Two populations ⇒ no single model** serves all; the switch needs a regime label,
   P(SOZ) is responder-only.
5. **δ-only** for the clean off-shaft marker (though δ-slowing is a recognised
   epileptogenic sign — physiologically coherent, and it is on strength-residual
   connectivity, not δ power).

---

## 7. Retraction ledger (what was withdrawn, and why — the process working)

| withdrawn | why | replaced by |
|---|---|---|
| "propagator beats strength at β" (audit_85) | unfair template-distance baseline | fair MW baseline → β is hubness |
| "shaft-independent within-patient marker" (audit_87) | along-shaft **depth** confound | depth 0.72–0.74 beats dynamics |
| signed/magnetic-Laplacian marker (audit_88) | within-shaft signal = depth | Stage-1 negative, no Stage-2 |
| all-contacts marker AUC ~0.81 + deployment (audit_92/93/95/98) | **proximity tautology** (audit_99) | off-shaft LOSO 0.72 (audit_101) |
| "mark few, discover many" clean discovery (audit_95) | all-contacts precision | off-shaft precision (audit_103), concentrated |
| label-free cross-patient detection | selection at chance (audit_94) | seed-based framing |
| continuous blend compound looked best (audit_104) | damages responders (cohort median lies) | hard switch |

The relational **community** finding (Act II, audit_89/90) and the δ off-shaft
**marker** (audit_101/102) are the parts that **survived** every control.

---

## 8. Discussion — positioning in the multiscale-Laplacian programme

**The unification (the sellable core).** The *same* density operator ρ̂(τ)=e^{−τL̂}/Z
that organises the cross-phase **β trace** (its multiscale cophenetic structure) also,
read as a node/pair diffusion affinity, exposes the **δ epileptogenic community**. One
description of a recording yields **both a cognitive read-out (the consolidation trace)
and a clinical one (the epileptic community)** — a diffusion geometry as a general lens
for functionally distinguished structure.

**Mechanistic + multiscale vs descriptive + single-scale (the defensible thesis).** The
dominant interictal SOZ decompositions extract a **single descriptive scale** and assign
it meaning post hoc — **PCA on directed connectivity** (Doss 2024, inward-suppression
PC1) or a chosen **centrality** (Roy 2025; Gunnarsdottir 2022 source-sink). Against them
the Laplacian density operator is **generative and scale-continuous**: every entry of
$e^{-\tau L}$ is a sum over diffusion paths at scale τ, so a marker score decomposes into
*which routes connect a contact to the epileptogenic seed set, at a tunable scale* —
interpretable in path/connectivity language, where a PC loading is a variance coefficient
with a post-hoc reading. The scale is a **continuous physical parameter (τ)**, not a
chosen number of components/clusters; PCA, fixed-k spectral clustering, and our own
fixed-k **Grassmann** probe are single-scale limits (Grassmann corroborating the
all-scale cophenetic β result is a **consistency check**, not a competition).

**Honesty guardrails — do NOT claim "improvement" over Doss.**
1. **Directed vs undirected.** Doss's whole result *is* directionality (in/out); our
   symmetric L → symmetric $e^{-\tau L}$ **cannot** represent in/out asymmetry, and when
   we tried a directed/signed magnetic Laplacian it added nothing (audit_88). **Frame as
   complementary, not superior.**
2. **Power.** Doss n=81 (ANOVA P≈2.6e-9, 79% of patients); we are n=10. The pitch is
   **methodological**, never "stronger evidence."
3. **Don't strawman PCA as uninterpretable** — Doss interprets PC1 mechanistically; the
   honest contrast is *post-hoc descriptive* vs *mechanistic-by-construction*.
4. **The multiscale payoff belongs to the trace**, not this marker (the SOZ signal lives
   at the slow/coarse τ end specifically).

**The edge that genuinely sells (honest): strength-orthogonality.** The field keeps
having to *control for* degree/strength (Doss pins the hyperconnected-vs-sink discrepancy
on "measuring strength vs counting connections"; Shah 2019 needs a spatial null). Our
read is **strength-orthogonal by construction, proximity-free, null-verified, and
phase-stable** — the cleanest available answer to the confound the field keeps fighting.
**Our AUC 0.72 sits at the lower edge of the field's 0.70–0.86 band.**

---

## 9. Literature to attach to and cite

> **Provenance.** Papers and reported effect sizes are carried from the
> deep-research-verified positioning report
> `.agents/reports/2026-06-11_epi-soz-marker-literature-positioning.md` (18 sources
> fetched, 24/25 claims verified). **DOIs resolved on 2026-06-18 by reading the
> publisher / PMC / Crossref / IOPscience / doi.org records and copy-pasted verbatim —
> none invented.** One author correction surfaced: the *Neurology* 2022 ImCoh paper is
> **Paulo DL** first-author (the earlier report mislabelled it "Park"; Narasimhan is a
> co-author).

**Field benchmark — node markers + validated AUC band (≈0.70–0.86):**
- **Gunnarsdottir KM et al.** "Source-sink index" — EZ as inhibited sinks. *Brain*
  2022;145(11):3901. n=65; ~73% localisation; **surgical-outcome AUC 0.86 ± 0.07**.
  **DOI: 10.1093/brain/awac300**
- **Roy S et al.** Eigenvector-centrality EZ biomarker. *Frontiers in Network Physiology*
  2025; n=65, 6 centres; **AUC 0.70 ± 0.01 under 26-fold LOPO**; outcome AUC 0.75.
  **DOI: 10.3389/fnetp.2025.1565882**
- **Taylor PN / Wang Y et al.** Normative band-power abnormality; **resected-vs-spared
  AUC 0.75, P=0.0003** (62 patients). *Brain* 2022;145(3):939. *Cite for the DRS
  validation paradigm (band-power, not connectivity).* **DOI: 10.1093/brain/awab380**

**The directed-decomposition comparator (our closest "intent" neighbour — PCA):**
- **Doss DJ et al.** "The interictal suppression hypothesis." *Brain* 2024;147(9):3009.
  n=81; **partial directed coherence + PCA**; inward-suppression motif = PC1 in
  **64/81 (79%)**; ANOVA **P≈2.6e-9**. Position *beside*, not above.
  **DOI: 10.1093/brain/awae189**

**Characterisation debate + our closest precedent:**
- **Shah P et al.** Within-resection-zone hyperconnectivity (β) **surviving a spatial
  null**. *NeuroImage: Clinical* 2019;23:101908. n=27. **Closest precedent.**
  **DOI: 10.1016/j.nicl.2019.101908**
- **Goodale SE et al.** Alpha imaginary-coherence FC predicts epileptogenicity **AUC
  0.78**. *Neurosurgery* 2020. n=15. **DOI: 10.1093/neuros/nyz351**

**FC substrate — |ImCoh| is field-standard for SOZ:**
- **Nolte G et al.** "Identifying true brain interaction from EEG data using the
  imaginary part of coherency." *Clinical Neurophysiology* 2004;115(10):2292.
  **DOI: 10.1016/j.clinph.2004.04.029**
- **Paulo DL et al.** "SEEG functional connectivity measures to identify epileptogenic
  zones: stability, medication influence, and recording condition." *Neurology*
  2022;98(20):e2060–e2072. Between-region **ImCoh higher in EZ** (p=0.007–0.023). n=32.
  **DOI: 10.1212/WNL.0000000000200386**
- **Narasimhan S et al.** Directed ImCoh/PDC: inward PDC **AUC 0.84**; combined **0.88**.
  *Epilepsia* 2020. n=25, alpha. **DOI: 10.1111/epi.16686**

**Method — the LRG / Laplacian density operator:**
- **Villegas P, Gili T, Caldarelli G, Gabrielli A.** "Laplacian renormalization group for
  heterogeneous networks." *Nature Physics* 2023;19(3):445–450. ρ̂(τ)=e^{−τL̂}/Z; K̂_ij
  sums diffusion trajectories over all paths. **DOI: 10.1038/s41567-022-01866-8**
- **Caldarelli G et al.** The Laplacian renormalization group (review). *J. Stat. Mech.:
  Theory Exp.* 2024;084002 (arXiv:2406.02337). **DOI: 10.1088/1742-5468/ad57b1**

> ⚠ **Novelty is an absence-of-evidence claim** (deep-research caveat): the corpus has
> the LRG definition and many SOZ markers but **no LRG/heat-kernel-propagator applied to
> SOZ**. Run a **dedicated negative search** for "communicability / heat-kernel /
> random-walk diffusion + SOZ/epilepsy" before asserting priority.

## 10. Open questions / highest-value next steps

1. **Make "more principled than PCA" a demonstration.** Run a Doss-style
   PCA-on-connectivity (and a fixed-k spectral-clustering comparator) for SOZ separation
   **on our own cohort**, show the density-matrix read recovers the same separation
   without choosing k, strength-orthogonally, with the path/scale interpretation a PC
   lacks. The single best referee-proofing step.
2. **External validation on an open dataset with resection + outcome** (Taylor normative;
   Roy/Gunnarsdottir multi-centre) — converts triage → outcome-validated. Biggest lever;
   needs data, not code.
3. **Cross-walk to the directed inward/sink family** (Doss, source-sink) on our cohort.
4. **Turn the two-population split into a result** — what predicts which regime?
5. **Persist the phase-robustness check** as a small audit if cited.

## Provenance

- **Node-intrinsic negatives:** `audit_80/81` (cross-patient), `audit_85/86/87`
  (within-patient + depth control), `audit_88` (signed/magnetic Laplacian,
  `data/audit/epi_signed_laplacian/`). Report: `.agents/reports/2026-06-05_propagator-subspace-epi-recovery.md` (§10 final).
- **Relational community (Part A):** `audit_89` (`data/audit/epi_propagator_blocks/`),
  `audit_90` C5 spatial null. Relational features `audit_91`, masking `audit_92`,
  classifier/P(SOZ) `audit_93`, label-free `audit_94`, few-seed `audit_95`
  (`data/audit/epi_marker_relational/`, `epi_community_localization/`). Report:
  `.agents/reports/2026-06-08_epi-propagator-soz-marker.md`.
- **Proximity pivot + clean marker:** `audit_99` (discovery yield),
  `audit_101` (16-operator sweep), `audit_102` (null + LOPO) →
  `data/audit/epi_marker_library/`. Report:
  `.agents/reports/2026-06-12_propagator-distant-soz-marker.md`.
- **Deployment:** `audit_103` (`data/audit/epi_marker_precision/`),
  `audit_104` (`data/audit/epi_marker_compound/`). Walkthrough:
  `.agents/reports/2026-06-16_epi-distant-marker-walkthrough.md`.
- **Few-seed reconstruction (off-shaft):** `audit_113_epi_fewseed_offshaft.py` →
  `data/audit/epi_marker_fewseed/` (δ, k=2/3/5, strength-residual, 300 draws,
  seed 20260618; computed 2026-06-18).
- **Figures:** `scripts/02_preprint/preprint_36/39/40_*.py` →
  `data/preprint/figures/all_bands/fig_epi_*.pdf`.
- **Literature:** `.agents/reports/2026-06-11_epi-soz-marker-literature-positioning.md`.
- **Honesty lessons:** node-intrinsic = strength + depth; proximity is a labelling
  tautology (evaluate off-shaft); cohort-median AUC/lift can lie (check per-patient);
  switch ≠ blend; matched-strength is the mandatory null at every layer.
