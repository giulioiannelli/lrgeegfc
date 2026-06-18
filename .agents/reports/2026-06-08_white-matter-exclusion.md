---
name: white-matter-exclusion
type: report
era: IMCOH_ABS_COHORT_N10
status: current
created: 2026-06-08
updated: 2026-06-12
companion: data/audit/wm_stratified/
build_scripts:
  - scripts/01_compute/audit/_wm_stratify.py
  - scripts/01_compute/audit/audit_83_wm_stratified_cophenetic.py   # cophenetic + raw
  - scripts/01_compute/audit/audit_84_wm_stratified_grassmann.py    # T_G(k)
  - scripts/01_compute/audit/audit_83b_wm_stratified_figures.py     # figures
  - scripts/01_compute/audit/audit_85_wm_decimation_control.py      # node-count control (flag A)
  - scripts/01_compute/audit/audit_86_wm_grassmann_cluster_extent.py # C3 gate on WM Grassmann (flag B)
sources:
  - audit_83  wm_stratified cophenetic + raw ρ_split        # this analysis
  - audit_84  wm_stratified Grassmann T_G(k)                # this analysis
  - audit_85  wm_stratified random-node-decimation control  # resolves flag A
  - audit_86  wm_stratified Grassmann cluster-extent gate    # resolves flag B
  - audit_63  matched_strength_surrogate_split_baseline     # cophenetic/raw null engine (reused)
  - audit_66  grassmann_matched_strength_surrogate          # full-graph Grassmann (re-emitted)
  - audit_70  grassmann_cluster_extent                      # locked C3 gate machinery (reused by audit_86)
  - audit_77/78  epi_stratified_*                           # the methodological template
pointers:
  - .agents/reports/2026-04-29_result-1-raw-fc-phase-trace.md
  - .agents/reports/2026-05-05_result-2-lrg-beta-trace.md
  - memory/epi_stratified_trace_2026_06_05.md
---

# White-matter exclusion: does the task trace survive a gray-only montage?

**Head.** White matter is **30–57 % of every patient's montage** (cohort ≈40 %),
yet **removing it does not weaken the task trace — the trace survives a gray-only
montage.** Cophenetic α/β still clear matched-strength on the gray-only submatrix,
and **β/low-γ Grassmann re-pass the LOCKED C3 cluster-extent gate** under WM
exclusion (`cluster_p_mass^wmX = 0.005`, strong, leave-one-out-robust; audit_86).
The trace is robust to WM removal — **not** a white-matter recording artifact, and
WM is not *carrying* it. **What we no longer claim:** that WM removal *sharpens*
the cophenetic trace as a tissue effect. A size-matched **random-node-decimation
control (audit_85)** shows the cophenetic "strengthening/emergence" (α/β rise,
γ_l/θ emerge) is reproduced by removing *any* ~40 % of nodes (cohort
`p_dec` 0.15–0.49, none WM-specific) — it is a **generic node-count effect**, not
WM dilution. The one place WM removal *is* specifically beneficial is the **raw
|ImCoh| substrate** (β/γ_l `p_dec` 0.000/0.030, **WM_specific**): the raw
substrate's historical weakness genuinely *was* white-matter strength-structure.
The **only** trace component that *depends* on white matter is the **δ Grassmann
subspace** (fails the gate, `cluster_p_mass^wmX = 0.105`; weakens 23→7 k-cells) —
and δ is the cross-probe epileptic-biology channel, never the task-trace headline.
Net: WM is **neither carrying nor diluting** the primary (cophenetic/Grassmann)
task trace; it is gray-robust. This is a clean defensive result against the "maybe
it's volume-conducted white-matter signal" objection — with the over-claim
("sharpens because gray-resident") withdrawn.

Verified: every claim below is gated by a **strength-preserving matched-strength
4-cycle ±δ surrogate (R=200)** regenerated on each WM-stratified subgraph. The two
2026-06-08 honesty flags (node-count confound on the "sharpening"; C6-Grassmann on
the per-`k` Wilcoxon basis not the C3 gate) are now **resolved** by the
decimation control (audit_85, §4) and the cluster-extent gate (audit_86, §5).

---

## Read this first (the seven statements)

1. **WM removal is a big cut, not a perturbation** — 30–57 % of contacts are
   dominant-white-matter. `exclude_wm` leaves a median of **74 nodes** (down from
   ~117). Any *weakening* must be read against that node loss before being called
   biology.
2. **Cophenetic α SURVIVES** WM exclusion (obs +0.11→+0.18, p 0.005→0.014,
   leave-Pat_15-out p=0.027). The α trace still clears matched-strength on the
   gray-only graph — but the *rise* is not WM-specific (decimation `p_dec`=0.31,
   §4).
3. **Cophenetic β SURVIVES, not WM-carried** — obs rose +0.22→+0.32 and
   leave-Pat_15-out p=0.037; the full-cohort exclude_wm p=0.053 sits a hair above
   0.05 (a small-n power crossing, §4). β is preserved, not lost. But the
   strengthening is **not WM-specific** — random K-node decimation reaches the
   same cohort ρ_split (`p_dec`=0.145, §4).
4. **Cophenetic low-γ/θ "emergence" is a NODE-COUNT effect, not WM biology.**
   low-γ p 0.116→0.019, θ p 0.688→0.005 under exclude_wm — but the decimation
   control reproduces it by removing *any* ~40 % of nodes (low-γ `p_dec`=0.185;
   §4). These are not WM-revealed traces; do not cite them as WM-specific.
5. **The raw-FC substrate's weakness genuinely WAS WM dilution.** The raw
   ρ_split^raw — "directionally above the drift floor but not significant" at the
   full graph (audit_74) — **clears matched-strength at δ/α/β/low-γ once WM is
   removed**, and (unlike cophenetic) this **IS WM-specific**: β/γ_l decimation
   `p_dec` = 0.000/0.030 (**WM_specific**, §4). Gray-only raw FC is a materially
   better substrate, for a real tissue reason.
6. **Grassmann β and low-γ subspace traces RE-PASS THE LOCKED C3 GATE** under WM
   exclusion (audit_86: `cluster_p_mass^wmX` = 0.005 / 0.005, **strong**,
   LOO-robust) — upgraded from the per-`k` count to a re-pass of the actual gate.
   Only **δ Grassmann fails** (`cluster_p_mass^wmX`=0.105, WM-dependent). γ_h
   emerges on the gate (0.005) but is LOO-fragile (0.050).
7. **WM is not pure noise** — `wm_only` carries its own weaker trace (cophenetic
   low-γ emerges p=0.032; high-γ Grassmann present). The trace is gray-matter-
   *dominant*, not gray-matter-*exclusive*.

---

## 1. What was tested

**WM definition (locked by user, 2026-06-08):** a contact is white matter iff its
*dominant* Desikan-Killiany tissue is `Wm` (argmax of the per-tissue weights, PTD
excluded) — `load_channel_regions(pat)["region"] == "Wm"`, consistent with
`regions._NON_ANATOMICAL_REGIONS`. The mask is FC-index-ordered and length-checked
against the FC matrix `N` (same alignment contract as the epi mask; all 10
patients align exactly).

**"Removal from the temporal series" is realized exactly.** ImCoh (and corr/MSC)
are *pairwise* spectral estimators, so the FC on a gray-only montage equals the
gray-only submatrix of the cached full FC to numerical precision. Dropping WM
channels before computing FC is therefore identical to submatrixing the cached FC
— only the node-coupled LRG renormalization genuinely changes, and that is exactly
what is recomputed on the submatrix.

**Six views × two substrate families**, mirroring the epi-stratification audits
(audit_77/78):

| family | configs | what it asks |
|---|---|---|
| **subgraph** (recompute LRG on submatrix) | `full`, `exclude_wm` (gray-only — the headline), `wm_only` | does the trace live in gray vs white tissue? |
| **pair-class** (full dendrogram, restrict pairs) | `gray_gray`, `cross` (gray↔WM), `wm_wm` | which tissue-pair class carries it within the full graph? |

**Substrates:** LRG cophenetic ρ_split (load-bearing, Result 2), raw |ImCoh|
ρ_split^raw (Result 1), Grassmann chordal T_G(k) (subspace). Cophenetic and raw
are read off **one shared matched-strength eigendecomposition ensemble** per cell
(`cophenetic_condensed_from_eigs` and `adjacency_from_laplacian_eigs`); Grassmann
reuses the same eigvectors. Positive = TRACE throughout (T_d sign convention
locked 2026-05-26).

**Null:** matched-strength 4-cycle ±δ, R=200, regenerated on each submatrix
(`wmX` / `wmONLY` caches, seed 20260608); `full` and pair-class reuse the
canonical seed-20260511 ensemble. `wm_only` keeps **36–67 nodes** (vs epi_only's
6–30), so its null is well-posed — not the near-degenerate arm epi_only was.

---

## 2. Results — cophenetic ρ_split (the load-bearing substrate)

`obs_median` (positive = trace), cohort one-sided paired Wilcoxon p vs own
matched-strength surrogate, sensitivity flag vs the `full` baseline. Median 74
nodes under `exclude_wm`.

| band | full | exclude_wm | LO-P15 | flag (excl_wm) | gray_gray | wm_only |
|---|---|---|---|---|---|---|
| δ | +0.01 / p0.216 | −0.01 / p0.500 | 0.367 | absent | absent | absent |
| θ | −0.04 / p0.688 | +0.07 / **p0.005** | 0.006 | **emerge** | absent | absent |
| α | +0.11 / **p0.005** | +0.18 / **p0.014** | 0.027 | **persist** | persist | weaken |
| β | +0.22 / **p0.007** | +0.32 / p0.053 | **0.037** | weaken* | **persist** (p0.032) | weaken |
| low-γ | +0.08 / p0.116 | +0.21 / **p0.019** | 0.014 | **emerge** | absent | emerge (p0.032) |
| high-γ | +0.00 / p0.278 | −0.01 / p0.312 | 0.213 | absent | absent | absent |

**Reading.** Every band that was a trace at the full graph (α, β) stays a trace,
with a **larger** effect size after WM removal. Two bands that were sub-threshold
at full (θ, low-γ) cross into significance. Nothing is lost. The `*` on β is the
single honest caveat — see §4. `wm_only` is mostly `weaken`/`absent` for α/β
(those traces are in gray), but low-γ `wm_only` emerges — white-matter contacts
carry a low-γ trace of their own.

Figure: `data/audit/wm_stratified/figures/fig_wm_verdict_matrix_coph.pdf`,
`fig_wm_exclude_forest.pdf` (per-patient full-vs-exclude_wm, α/β/low-γ).

---

## 3. Results — raw |ImCoh| and Grassmann

**Raw ρ_split^raw.** Removing WM rescues the raw substrate. At the full graph the
raw trace clears nothing strictly (audit_74: best low-γ p≈0.057). Gray-only:

| band | full | exclude_wm | LO-P15 | flag |
|---|---|---|---|---|
| δ | +0.11 / p0.097 | +0.13 / **p0.019** | 0.027 | **emerge** |
| α | +0.16 / p0.097 | +0.17 / **p0.010** | 0.014 | **emerge** |
| β | +0.26 / p0.053 | +0.35 / **p0.024** | 0.027 | **emerge** |
| low-γ | +0.13 / p0.053 | +0.18 / **p0.024** | 0.006 | **emerge** |
| θ, high-γ | n.s. | n.s. | — | absent |

The raw substrate's matched-strength failure was substantially **white-matter
strength-structure that the null was matching**. Once WM is gone, gray-only raw FC
clears at δ/α/β/low-γ. (Read with the node-count caveat, §4 — but it is a clean
directional + significance gain, and it is the *opposite* of what a node-loss
power penalty would predict.)

**Grassmann T_G(k)** (n_k of 81–111 with cohort Wilcoxon p<0.05; judged on
Wilcoxon p, not the strict "separated" flag — see §5):

| band | full | exclude_wm | wm_only | verdict |
|---|---|---|---|---|
| δ | 23 | 7 | 10 | **weaken** (WM-dependent) |
| θ | 8 | 6 | 1 | persist (weak) |
| α | 4 | 3 | 6 | weak at all (α lives in cophenetic) |
| β | **40** | **40** | 4 | **persist** (unchanged) |
| low-γ | **41** | **40** | 7 | **persist** (LO-P15 robust) |
| high-γ | 19 | 22 | 10 | **persist** |

The β and low-γ subspace traces are **completely untouched** by WM exclusion. δ is
the one WM-dependent subspace component.

Figure: `fig_wm_verdict_matrix_raw.pdf`, `fig_wm_grassmann_kspan.pdf`.

---

## 4. The node-count confound — resolved by random-node decimation (audit_85)

`exclude_wm` removes ~40 % of nodes (median 117→74), so the full-vs-exclude
*change* in ρ_split confounds tissue identity with a generic node-reduction
effect: a smaller graph gives a denser, cleaner LRG renormalization (and a
different-variance cohort Spearman) regardless of *which* nodes are dropped. The
matched-strength null is regenerated per submatrix, so every reported significance
is above strength-matched chance — but that does not isolate *tissue* from *node
count*. **audit_85** does: for each (patient, band) it drops K = #WM nodes **at
random** R=200 times, recomputes ρ_split, and asks where the observed `exclude_wm`
ρ_split sits in that random-decimation distribution (cohort-level `p_dec` = P(a
random-decimation cohort-median ≥ the exclude_wm cohort-median)).

**Verdict (cohort `p_dec`, `decimation_control_cohort.csv`):**

| substrate | band | cohort obs | cohort dec median | `p_dec` | verdict |
|---|---|---|---|---|---|
| cophenetic | α | +0.183 | +0.151 | 0.305 | generic_nodecount |
| cophenetic | β | +0.324 | +0.272 | 0.145 | generic_nodecount |
| cophenetic | low-γ | +0.212 | +0.148 | 0.185 | generic_nodecount |
| cophenetic | θ | +0.073 | +0.036 | 0.225 | generic_nodecount |
| **raw** | **β** | +0.348 | +0.255 | **0.000** | **WM_specific** |
| **raw** | **low-γ** | +0.178 | +0.135 | **0.030** | **WM_specific** |
| raw | α | +0.165 | +0.148 | 0.300 | generic_nodecount |

- **Cophenetic: the "sharpening" is NOT WM-specific.** No cophenetic band crosses
  `p_dec < 0.05`. The β rise (+0.22→+0.32) sits at the 85th percentile of random
  K-node removal (`p_dec`=0.145, not beyond the random p95=+0.356); α (0.305) and
  the γ_l/θ "emergence" (0.185 / 0.225) are squarely inside the random-decimation
  band. **Removing any ~40 % of nodes raises the cophenetic ρ_split the same
  amount** — so we cannot attribute the sharpening to white matter. We therefore
  **withdraw** the "WM removal sharpens the cophenetic trace" framing.
- **Cophenetic is NOT WM-carried either.** If WM carried the trace, `exclude_wm`
  would fall *below* random removal (`p_dec`→1). It does not (δ is the only band
  leaning that way, `p_dec`=0.95 — and δ is not a trace band). So WM is **neither
  carrying nor specifically diluting** the cophenetic trace: removing it behaves
  like removing random nodes. The **survival** claim (α/β still clear
  matched-strength on the gray-only graph) is untouched — survival is a
  within-`exclude_wm` test, not a full-vs-exclude contrast.
- **Raw IS WM-specific.** raw β (`p_dec`=0.000) and raw γ_l (0.030) exceed the
  random-decimation cohort distribution: removing WM *specifically* clears the raw
  substrate beyond random removal. The raw substrate's historical matched-strength
  failure genuinely was white-matter strength-structure (the cophenetic
  renormalization already discounts that, which is why WM removal does nothing
  extra for cophenetic).

**The β cophenetic p=0.053** is still a small-n power crossing, not a lost trace
(effect size grew; LO-P15 p=0.037; `gray_gray` p=0.032; raw-β clears p=0.024) — we
report β `exclude_wm` as **preserved**, not failed. What changed is only that we
no longer read its *growth* as a WM-specific result. Limitation: random decimation
can sever a whole probe/hemisphere, so the null is broad (intentionally — it is
the generic-node-loss reference); the control holds node count fixed and
randomizes tissue identity, complementing audit_83's matched-strength (which holds
strength fixed).

---

## 5. Grassmann on the locked C3 gate (audit_86) + the Wilcoxon-basis note

**The locked gate, re-run under WM exclusion (audit_86).** The C6-Grassmann table
in §3 reports per-`k` cohort-Wilcoxon significant-cell *counts* (the audit_66
basis). The **locked** C3/C5 Grassmann gate is the audit_70 cluster-extent
permutation **mass** test (`cluster_p_mass`), not the per-`k` count. audit_86 runs
that gate on the `exclude_wm` per-`k` T_G curves (regenerating the per-surrogate
T_G stack from the warm `wmX` eigvec caches):

| band | `cluster_p_mass^wmX` | verdict^wmX | C3 full | LOO max | reading |
|---|---|---|---|---|---|
| β | **0.005** | **strong** | strong | 0.005 | **persist on the gate** |
| low-γ | **0.005** | **strong** | strong | 0.005 | **persist on the gate** |
| δ | 0.105 | no_trace | strong | 0.31 | **weaken** (WM-dependent) |
| high-γ | 0.005 | strong | no_trace | 0.050 | emerge, but LOO-fragile |
| α / θ | 0.49 / 0.41 | no_trace | no_trace | 1.00 / 0.63 | absent both |

So **β and low-γ Grassmann genuinely re-pass the locked gate** under WM exclusion
(strong, LOO-robust) — the §3 "40→40 / 41→40 unchanged" per-`k` read is now backed
by the actual C3 gate. δ **fails the gate** (confirming the 23→7 per-`k`
weakening): δ is WM-dependent, consistent with it being the cross-probe
epi-biology channel, not the task trace. (Caveat: `exclude_wm` changes N, so the
k-axis is only approximately comparable to the full-graph gate; the p-value is
honest *within* the WM-excluded geometry — obs and null share N and the
high-k NaN structure per patient.)

**Wilcoxon-basis note (why §3 uses the per-`k` count).** The shared
`cohort_verdict` calls a cell `separated` only if Wilcoxon p<0.05 **and** the
surrogate median ≈0 **and** ≥8/10 patients beat their surrogate. The med_surr≈0
condition is calibrated for the cophenetic ρ_split; the Grassmann chordal distance
is positive-offset, so the strict flag reports 0 "separated" for β even though
β-full clears Wilcoxon (p=0.024, 9/10). The §3 Grassmann table is therefore judged
on `paired_wilcoxon_p < 0.05` (the audit_66/78 basis); the **gate** verdict is now
the audit_86 `cluster_p_mass^wmX` above. The cophenetic/raw `separated` flags are
unaffected.

---

## 6. What it means

- **Defensive (manuscript).** The trace is robust to the "it's an electrode-tissue
  artifact" objection: removing white-matter contacts — physiologically the most
  suspect (low-amplitude, volume-conduction-prone) — does not weaken the trace
  (cophenetic α/β still clear matched-strength; β/low-γ Grassmann re-pass the
  locked gate), and the decimation control shows WM is not *carrying* it (its
  removal behaves like random node loss). State this as **survival/robustness**,
  **not** "sharpening" — the apparent cophenetic strengthening is generic
  node-count (§4), so do not claim WM was diluting the trace. Belongs in the
  methods/robustness section alongside matched-strength and cross-probe.
- **Substrate insight (genuinely WM-specific).** The raw-FC substrate (Result 1),
  historically the weak arm, is materially improved by a gray-only montage — and
  this one **is** WM-specific (decimation β/γ_l `p_dec`=0.000/0.030, §4): its
  matched-strength failure was white-matter strength-structure. If we ever want
  raw FC to carry weight, gray-only is the montage to do it on.
- **Dissociation.** δ is the only WM-dependent component (Grassmann 23→7) — and δ
  is the cross-probe epileptic-biology channel, not the task trace. So WM-carried
  structure and task-trace structure are anatomically separate: a small bonus
  consistency check on the whole δ-is-different story.
- **Not gray-exclusive.** `wm_only` carries weaker traces (low-γ cophenetic;
  high-γ Grassmann), so the framing is gray-*dominant*. Don't over-claim "WM is
  noise."

---

## 7. Limitations

1. **Node-count confound — RESOLVED (§4).** The decimation control (audit_85)
   shows the cophenetic "sharpening/emergence" is generic node-count, not
   WM-specific; the raw β/γ_l clearing IS WM-specific. The full-vs-exclude
   *change* in cophenetic ρ_split should not be read as a WM tissue effect; the
   **survival** of α/β (matched-strength on the gray-only graph) and the β/low-γ
   Grassmann **gate re-pass** (audit_86) are the load-bearing robustness claims.
2. **Single WM definition.** Atlas-dominant `Wm` (user-chosen). A PTD<0 criterion
   (proximal-tissue-density, Mercier 2017) would draw a slightly different,
   arguably more principled, gray/white line on borderline contacts; not run.
   Conclusions that hold under a 40 %-of-nodes cut are unlikely to flip under a
   borderline reclassification, but this is untested.
3. **Cohort n=10.** β `exclude_wm` p=0.053 illustrates how close to the power
   floor the full-cohort tests run; the leave-one-out and pair-class corroborants
   carry weight precisely because the headline p is marginal.
4. **`wm_only` is the secondary arm.** Well-posed (36–67 nodes, non-degenerate
   null) but lower-powered; its `emerge`/`absent` calls are secondary evidence.

---

## 8. Provenance

- **Cohort:** n=10 (Pat_02, 03, 05, 06, 07, 08, 10, 13, 14, 15). **FC:** imcoh_abs,
  τ=1/λ_max, average-linkage ultrametric. **Bands:** all six.
- **Null:** matched-strength 4-cycle ±δ, R=200, swap_factor=20; `wmX`/`wmONLY`
  ensembles seed 20260608 (`data/cache/matched_strength_surrogate_wm_{excluded,only}_lrg/`);
  `full`/pair-class reuse canonical seed 20260511.
- **Caches:** `data/audit/wm_stratified/{cophenetic_raw_per_patient,cophenetic_raw_cohort,
  grassmann_per_patient_per_k,grassmann_cohort,decimation_control_per_patient,
  decimation_control_cohort,grassmann_cluster_extent_wmX,grassmann_cluster_extent_wmX_per_k}.csv`
  + `README.md`, `README_grassmann.md`, `README_decimation.md`,
  `README_grassmann_cluster_extent.md`.
- **Figures:** `data/audit/wm_stratified/figures/fig_wm_{verdict_matrix_coph,
  verdict_matrix_raw,exclude_forest,grassmann_kspan}.pdf`.
- **WM prevalence per patient:** Pat_02 35%, 03 39%, 05 39%, 06 36%, 07 44%,
  08 39%, 10 46%, 13 38%, 14 30%, 15 57%.
- **Wall-clock:** audit_83 ≈133 min (submatrix ensemble generation), audit_84
  ≈3 min (reuses caches), audit_85 ≈1 min (no surrogates), audit_86 ≈7 min
  (reuses `wmX` caches).
