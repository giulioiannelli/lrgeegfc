---
name: writing-directive-methods-anatomy
era: IMCOH_ABS_COHORT_N10
status: directive
kind: writing-agent-directive
date: 2026-05-19
target: methods §sssec:methods_anatomy — complete the placeholder subsubsection
priority: high (final missing methods subsubsection blocking anatomy paragraph in §β results and downstream bands)
source_of_truth:
  - .agents/preprint/locked/ANATOMY_CONTROLS.md  # locked A1-A4 battery
  - .agents/preprint/locked/ANATOMY_LEDGER.md    # locked per-(band, probe) verdicts
  - .agents/preprint/locked/CONTROLS.md          # parent trace battery (C3 matched-strength = same surrogate as A3)
---

# Writing directive — complete `sssec:methods_anatomy` (2026-05-19; rev 2026-05-19 pm)

**Head.** The current LaTeX placeholder for §sssec:methods_anatomy mentions
two gates (hypergeometric + matched-strength surrogate) but does not formally
define the per-probe "trace-flagged unit", the surrogate construction, or the
joint-gate language. This directive supplies the exact prose + equations
needed to finalize the subsubsection. The locked battery in the **manuscript
is A1 (hypergeometric) + A2 (matched-strength surrogate)**. **A2 + A4 (the
deferred sampling-corrected bootstrap + implant-geometry regression) are NOT
disclosed in the manuscript methods** — they belong in the sensitivity
supplement with their own methodology description.

**Manuscript ↔ lab label mapping (locked 2026-05-19 pm).** In the manuscript:
**A1 = hypergeometric** (lab A1, unchanged) and **A2 = matched-strength
surrogate** (lab A3). In the lab artifacts (`ANATOMY_CONTROLS.md`,
`ANATOMY_LEDGER.md`, audit CSV column names, and the prose of this directive
outside the drop-in LaTeX block) the labels **A1 / A3 stay**, with **lab A2
(sampling-corrected bootstrap)** and **lab A4 (implant-geometry regression)**
reserved for the sensitivity supplement. The drop-in LaTeX block in this
directive uses **manuscript labels A1 / A2**; surrounding methodology prose
uses **lab labels A1 / A3** to keep traceability with the locked artifacts.

---

## What the placeholder currently has (correct fragments to preserve)

The existing paragraph correctly identifies:
- Two-gate structure (hypergeometric + matched-strength surrogate).
- BH-FDR correction across DK regions at `q_BH < 0.05`.
- Matched-strength empirical-null gate at `p_emp < 0.05`.
- Joint-gate semantics (a region passes only when it clears both gates).

These do not need to be rewritten. What is missing is the *what* each gate
operates on (the trace-flagged unit per probe), the *how* of the matched-
strength surrogate construction, the *separate* gate at Grassmann (where A1
is structurally sparse and A3 alone is the gate), and the explicit deferral
of A2 + A4.

---

## What the subsubsection must define (in order)

### 1. Trace-flagged unit per probe

The hypergeometric and matched-strength tests need a "trace-flagged unit" to
enrich over. The unit differs by probe:

- **Cophenet probe (`\rhosplit` on `\Dcoph`).** The trace-flagged unit is the
  set of pairs `(i, j)` whose per-pair contribution to the cohort `\rhosplit`
  signal sits in the **top decile** of `\lvert \Delta \rhosplit^{(i,j)} \rvert`,
  pooled across the cohort. Concretely, for each pair `(i, j)` compute
  `\Delta \rhosplit^{(i,j)} = \rhosplit^{coph}(i,j; \text{rspre_A}, \text{rspost})
  - \rhosplit^{coph}(i,j; \text{rspre_A}, \text{rspre_B})` (the per-pair contribution
  to the split-baseline minus split-drift difference), pool the absolute values
  across the `n = 10` patients, and select the top 10\% as the trace-flagged
  set.

- **Grassmann probe (`d_G(k)` on `U_k`).** The trace-flagged unit is the set of
  nodes whose **per-node participation** averaged across the significance-
  thresholded set `S(b) = {k : p_k(b) < α_k}` sits in the top decile, pooled
  across the cohort. For node `i` and cutoff `k`, participation is the squared
  loading of node `i` in the subspace `U_k = span{φ_2, …, φ_{k+1}}` of the
  combinatorial Laplacian (definition in §sssec:methods_lrg). `S(b)` is the
  **support of the cluster-mass statistic `T_G^*`** (the same set of k-cells
  whose per-k cohort-paired Wilcoxon clears `α_k`); there is no
  contiguous-significant window — under the locked all-clusters
  `T_G^*` paradigm, the retired `K*(b)` ("longest contiguous-significant
  window") is replaced by `S(b)` everywhere. The per-node statistic is the
  **unweighted average** across `S(b)`:
  ```
  part(i; b) = (1 / |S(b)|) · Σ_{k ∈ S(b)} ||U_k^T e_i||²_2
  ```
  Significance-weighted aggregation (`weight k by −log10 p_k(b)`) is **not**
  used — `T_G^*` already weights by `−log10 p_k`, so re-weighting here would
  double-count the same evidence. Unweighted aggregation reads the cleaner
  question "across significance-thresholded subspaces, where do the heavy
  contacts sit", independent of how strongly each k cell contributes to the
  cohort signal.

### 2. Hypergeometric test (A1)

For each DK region `r` with cohort-wide coverage, let `k_r` be the count of
trace-flagged units intersecting `r`, `K = \sum_r k_r` the total number of
trace-flagged units, `n_r` the number of pair- or node-positions in `r`, and
`N` the total number of pair- or node-positions in the cohort. Under the
one-sided Fisher's exact test, the per-region p-value is

```
p_{hyper}(r) = \mathbb{P}\bigl[\, K_r \ge k_r \mid \text{Hypergeometric}(N, n_r, K) \,\bigr].
```

Multiple-testing correction uses BH-FDR across the
`R = N_{DK regions with coverage}` regions per (band, probe), gating regions
at `q_{BH}(r) < 0.05`. The enrichment ratio cited in the ledger is
`k_r / \mathbb{E}[K_r] = k_r N / (n_r K)`.

### 3. Matched-strength surrogate test (A3, mandatory)

For each (patient, band, phase) the same matched-strength surrogate ensemble
used for control C3 (4-cycle ±`\delta` rewiring, `R = 200`, seed 20260511) is
rerun, and the entire A1 enrichment pipeline (recompute the per-pair `\Delta
\rhosplit^{(i,j)}` distribution or per-node participation across the cohort,
re-extract the top decile, and recompute `k_r` per DK region) is applied to
each surrogate. The matched-strength surrogate preserves per-node strength
while randomizing topology, so the A3 null isolates **topological**
enrichment from strength-driven enrichment.

For each DK region `r`, let `\mathrm{enr}_r^{obs}` denote the observed
enrichment ratio and `\{ \mathrm{enr}_r^{(s)} \}_{s=1}^{R}` the surrogate
ensemble. The matched-strength statistic is

```
z_r^{A3} = \bigl(\mathrm{enr}_r^{obs} - \overline{\mathrm{enr}_r^{(s)}}\bigr)
                 \big/ \mathrm{sd}\bigl(\mathrm{enr}_r^{(s)}\bigr),
```

with empirical p-value (`R + 1` plug-in)

```
p_{emp}(r) = \bigl(1 + \#\{s : \mathrm{enr}_r^{(s)} \ge \mathrm{enr}_r^{obs}\}\bigr) / (R + 1).
```

A region passes A3 iff **both** `p_{emp}(r) < 0.05` and `z_r^{A3} > 2.0`.

### 4. Joint-gate semantics (cophenet) and A3-only fallback (Grassmann)

- **Cophenet probe**: a region passes the **anatomical enrichment claim** iff
  it clears **both** A1 (`q_{BH}(r) < 0.05`) and A3 (`p_{emp}(r) < 0.05`,
  `z_r^{A3} > 2.0`).

- **Grassmann probe**: the per-node participation statistic does not partition
  the cohort into independent units the way pair counts do, and A1 enrichment
  ratios for individual regions are structurally compressed (often `<\!\!1\times`
  even for A3-significant regions). The **locked gate at the Grassmann probe
  is A3 alone** — `p_{emp}(r) < 0.05` AND `z_r^{A3} > 2.0`, with A1 reported
  descriptively. This A3-only convention is documented in `ANATOMY_LEDGER.md`
  and applies to β Grassmann, γ_l Grassmann, δ Grassmann (full), and δ
  Grassmann (epi-X).

### 5. Scope limit — A2 + A4 deferred to sensitivity supplement

Two further controls are part of the locked anatomy battery but **not gates**
in the present analysis:

- **A2 (sampling-corrected bootstrap)** — `N = 1000` patient-resampling
  bootstrap conditioned on each region's cohort-wide contact distribution;
  gates on 95\% CI lower bound `>\!\!1.0` for regions passing A1. Purpose: reject
  enrichment claims that arise solely from a region being sampled in many
  patients (e.g., hippocampus, heavily sampled in clinical sEEG practice).
- **A4 (implant-geometry cohort regression)** — per-contact regression of the
  trace-flag against DK label, distance to nearest epileptogenic-zone contact,
  centroid distance from the implant center of mass, and hemisphere. Flags
  geometric covariates that independently predict trace-flag at `q_{BH}<0.05`.

Both controls are deferred to the sensitivity supplement; their absence is
disclosed here so that the locked claim is **strong localized under A1 + A3
(or A3 alone at Grassmann), with A2 + A4 pending**. The methods subsubsection
must state this explicitly — concealing the deferral would be a verification
failure under `feedback_brutal_honesty_no_sycophancy.md`.

---

## Suggested LaTeX replacement (drop-in for the placeholder paragraph)

The writing agent may copy and edit; the structure below is the locked
content, in the prose register of the surrounding Methods.

```latex
\subsubsection{Anatomical enrichment of the per-pair and subspace probes}
    \label{sssec:methods_anatomy}
%
Where a trace exists, its anatomical localization is read by per-region
enrichment of a probe-specific \emph{trace-flagged} unit against a parcellation
of contacts into Desikan--Killiany (DK) cortical and subcortical regions. For
the per-pair probe (\(\rhosplit\) on \(\Dcoph\)), the trace-flagged unit is the
set of contact pairs \((i, j)\) whose per-pair contribution
\(\lvert \rhosplit^{\,(i,j)}(\text{rspre\_A}, \text{rspost})
  - \rhosplit^{\,(i,j)}(\text{rspre\_A}, \text{rspre\_B}) \rvert\)
sits in the cohort-pooled top decile. For the subspace probe (\(d_G(k)\) on
\(U_k\)), the trace-flagged unit is the set of contacts whose per-node
participation \(\norm{U_k^{\transp} \vec{e}_i}_2^2\) averaged across the
significance-thresholded set
\(\mathcal{S}(b) = \{\,k : p_k(b) < \alpha_k\,\}\) --- the support of the
cluster-mass statistic \(T_G^*\)
(Section~\ref{sssec:methods_compare_grassmann}) --- sits in the cohort-pooled
top decile. The aggregation is unweighted: every \(k \in \mathcal{S}(b)\)
contributes equally to the per-node statistic. The retired notion of a
``contiguous-significant \(k\)-window'' is not used under the locked
all-clusters \(T_G^*\) paradigm.
%
Two complementary tests gate the per-region enrichment claim. The
\textbf{hypergeometric test (A1)} counts trace-flagged units intersecting each
DK region \(r\) and tests the count against the marginal cohort coverage via
one-sided Fisher's exact, multiple-testing-corrected by Benjamini--Hochberg
across the \(R = N_{\mathrm{DK,covered}}\) regions per (band, probe), with
gate \(q_{\mathrm{BH}}(r) < 0.05\). The \textbf{matched-strength surrogate
test (A2)} replays the entire A1 pipeline (top-decile selection plus per-region
count) on \(R = 200\) strength-preserving (4-cycle \(\pm \delta\)) surrogate
adjacency ensembles --- the same surrogates used for control C3
(Section~\ref{sssec:methods_compare_ctm}) --- and gates the observed enrichment
against the surrogate distribution at
\(p_{\mathrm{emp}}(r) = (1 + \#\{s : \mathrm{enr}_r^{(s)} \ge \mathrm{enr}_r^{obs}\}) / (R + 1) < 0.05\)
jointly with \(z_r^{A2} > 2.0\). Because matched-strength rewiring preserves
per-node strength while randomizing topology, A2 isolates topological
enrichment from strength-driven enrichment.

For the per-pair probe, a region is reported as anatomically enriched only
when it clears both A1 and A2 (\textbf{joint gate}). For the subspace probe,
the per-node participation statistic compresses the A1 enrichment ratios
toward unity; the locked gate is A2 alone, with A1 reported descriptively. A
band passes the anatomical-localization claim at a probe when at least one
non-trivial (i.e., not white-matter, not unknown) DK region clears the
applicable gate, and is tagged \emph{strong localized} in the per-band brief.
```

The above LaTeX assumes the following macros are already defined elsewhere in
the manuscript: `\rhosplit`, `\Dcoph`, `\gls{bh}`, `\gls{lrg}`. If `\gls{bh}`
is not defined the text reads "Benjamini--Hochberg" verbatim (used above).

---

## What NOT to do

- Do **not** cite KC-era anatomy CSVs (`data/audit/lrg_localization_anatomy/`)
  or the retired "Hippocampus + left fusiform" claim from
  `memory/result_2_lrg_beta_trace.md`. These were retired by the 2026-05-18
  trace-side lockdown and further refined by the 2026-05-19 pm cluster-extent
  rerun; see `ANATOMY_LEDGER.md` revision-history note and the
  cross-band overlap table for the locked replacement (parahippocampal at β
  cophenet; Hippocampus at β Grassmann; **left fusiform appears nowhere
  under the locked cluster-extent paradigm `S(b)` — neither at β, γ_l, nor
  δ**).
- Do **not** introduce a per-band anatomy result in the Methods
  subsubsection. The methodology is band-agnostic; per-band verdicts go in
  the results sections, anchored to `ANATOMY_LEDGER.md` rows.
- Do **not** describe A2 or A4 as "future work" or "next steps" — they are
  **part of the locked battery** but deferred to the sensitivity
  supplement. The wording matters: "deferred to sensitivity supplement"
  preserves the load-bearing status of A1 + A3 without overclaiming.
- Do **not** silently equate A3 to C3. They use the **same surrogate ensemble**
  (`R = 200`, 4-cycle ±δ, seed 20260511), but A3 is a separate test (enrichment
  permutation per region) and the per-region empirical p-value is computed
  independently of C3's cohort-paired Wilcoxon. The Methods should say "the
  same matched-strength surrogate ensemble as for control C3" without
  conflating the tests themselves.

---

## Sanity checks for the writing agent before committing

1. The placeholder's two-gate structure is preserved.
2. The trace-flagged unit is defined per probe (top-decile pairs for cophenet;
   top-decile per-node participation for Grassmann).
3. The A3 surrogate construction is referenced as identical to C3's matched-
   strength ensemble (same `R`, same algorithm, same seed).
4. The A3-only fallback at Grassmann is stated explicitly.
5. A2 + A4 are named and deferred (not silently omitted).
6. No per-band verdict appears in the Methods subsubsection.
7. The four control names A1, A2, A3, A4 are stable across the manuscript and
   the locked `ANATOMY_CONTROLS.md` — do not rename mid-text.

---

## Cross-references downstream

After this subsubsection lands, the Results paragraphs that cite anatomical
localization at β (and α, γ_l, δ once those bands are drafted) inherit the
following hooks:

- "the joint-gate anatomical enrichment of Section~\ref{sssec:methods_anatomy}"
  (β cophenet, α cophenet)
- "the matched-strength-only anatomical enrichment of
  Section~\ref{sssec:methods_anatomy}" (β Grassmann, γ_l Grassmann, δ Grassmann)
- "the deferred A2 + A4 sensitivity controls
  (Section~\ref{sssec:methods_anatomy})" — for any anatomy paragraph that
  flags the deferral

All region citations in the Results sections must point to the corresponding
`data/audit/anatomy_<band>_<probe>[_epiX]/cohort_summary.csv` row, per the
`ANATOMY_LEDGER.md` per-band verdict-detail subsections. Verification
follows the `EVALUATION_PROTOCOL.md` checklist item 1 (Numbers — CSV-row
trace).

---

## Revision history

- **2026-05-19** — Initial directive. Supplies prose, equations, joint-gate
  semantics, A3-only Grassmann fallback, and A2/A4 deferral language for
  §sssec:methods_anatomy.
