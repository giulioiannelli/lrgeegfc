---
name: writing-directive-grassmann-subsection-results
era: IMCOH_ABS_COHORT_N10
status: current
kind: directive
scope: Results encapsulation for the Grassmann subsection (ssec title "Grassmann-analysis probe on communication subspaces: the global multiscale counterpart for per-pair measures"). What the whole-network subspace probe SHOWS at cohort level + how it clears/does-not-clear controls + per-band reading of the all-bands compound figure + how the cross-probe pattern illuminates the ρ^coph per-pair result. Measure definition (principal angles, U_k, T_G) belongs to Methods, NOT here. Changes no locked verdict.
created: 2026-06-04
companion: directives/writing_directive_2026-06-04_coph-subsection-results.md (the per-pair twin this subsection complements); directives/writing_directive_2026-06-03_raw-copula-figure-swap.md (raw). Figure: data/preprint/figures/all_bands/fig_bands_grassmann_heatmap.pdf (first figure, all-bands compound).
verified_against: .agents/preprint/locked/VERDICT_LEDGER.md (Grassmann verdict tiers, Decision 12 2026-05-28); data/audit/grassmann_cluster_extent/cohort_summary.csv (cluster_p_mass, longest run, T_G*, LOO max); data/audit/matched_strength_surrogate_split_baseline/cohort_summary.csv (ρ^coph gate, for the cross-probe table).
---

# Grassmann subsection — results encapsulation

**Head.** The Grassmann probe asks the per-pair question one level up: does the
**whole leading-mode communication subspace** rotate back toward the task
configuration in post-task rest, rather than do the individual pair distances.
Read against the per-pair `ρ^coph` result, it **disentangles the trace by
structural scale**. Only **β lights both probes** — its reorganization is
coherent at the per-pair *and* the low-rank-subspace level, which is what makes
β the **primary, multi-probe-convergent** trace band. **γ_l (and weakly δ)
light the subspace probe but are silent per-pair**; **α lights the per-pair
probe but is silent in the subspace**; **θ and γ_h are silent on both**. The
two probes are **complementary, not redundant**: `ρ^coph` alone would have
missed γ_l, the subspace probe alone would have missed α. All subspace
verdicts clear (or fail) the **matched-strength cluster-mass gate** — the same
mandatory control as the per-pair probe.

> **Do not describe the measure here** (principal angles, the U_k subspace, the
> per-mode Δθ_i(k), T_G, cluster-extent permutation) — that is Methods. This
> subsection reports *what the subspace probe shows*, *how it behaves under the
> controls*, and *what the cross-probe pattern means for the `ρ^coph` reading*.

---

## 1. The result to encapsulate (one short paragraph)

Frame it as **disentanglement by structural scale**, the global-multiscale
complement of the per-pair subsection. The subspace probe re-asks the trace
question for the *collective* leading-mode geometry instead of for individual
pairs, and the bands sort differently than they did per-pair: the subspace
trace is **strong in β and γ_l, weak in δ, absent in α, θ and γ_h**. The
single fact that organizes the subsection is the **overlap with the per-pair
result is exactly one band — β** — and the **two mismatches (α per-pair-only,
γ_l subspace-only) are informative**, not noise: they localize each band's
reorganization to a different level of the hierarchy (§4). This is the
"global multiscale counterpart" the subsection title promises.

## 2. First figure — `fig_bands_grassmann_heatmap.pdf` (per-band comments)

All-bands 2×3 compound (top row δ θ α; bottom row β γ_l γ_h), the global
counterpart of the per-pair joint-density compound. Each panel is the
matched-strength residual on the mode-index `i` × subspace-cutoff `k` triangle,
weighted by the per-k cohort evidence and a band-level verdict scale; **red
below the `i = k` diagonal = trace** (post-task subspace closer to task than
pre-task rest, in excess of the matched-strength expectation), **blue =
anti-trace**, **white/cream = no excess**. Short reading, band by band:

- **β** — strong, extended red trace band running below the diagonal across the
  intermediate cutoffs (k ≈ 20–60), with the complementary anti-trace blue
  beneath it at higher mode index; the textbook subspace trace.
- **γ_l** — the most saturated panel: a broad strong-red region across k ≈ 12–80.
  The subspace trace is at its largest here — but **only here**, not per-pair.
- **δ** — a moderate, less coherent red residual at low-to-mid k, rendered at
  **half intensity** to mark its **weak** tier (see §3, §5 — Pat_08-leveraged at
  full data).
- **α** — faint red wisps only, no coherent extended trace region: the subspace
  probe is **silent** at α. Contrast directly with α's clean green per-pair
  diagonal in the `ρ^coph` figure — same band, opposite probe verdict.
- **θ, γ_h** — washed out (θ a near-empty panel; γ_h a faint near-diagonal wisp
  at high k that sits at the null). No subspace trace.

The one-sentence contrast to make explicit: where the per-pair figure lit
**α and β**, the subspace figure lights **β, γ_l and (weakly) δ** — and **β is
the only band lit in both**.

> **Figure-vs-ledger note for the writing agent (do not paraphrase into the
> text — it is a provenance caveat).** The δ panel is at half intensity on
> purpose. The per-band `_visual_bootstrap` PDFs scale each band by the raw
> `cluster_p_mass` column of `grassmann_cluster_extent/cohort_summary.csv`,
> which still tags δ "strong" (pre-Decision-12). The locked ledger demoted δ to
> **weak** (Decision 12, LOO binds at Pat_08). This compound encodes the
> **locked tier** (β/γ_l strong → full, δ weak → half, θ/α/γ_h → muted), so it
> is the ledger-faithful figure. Cite δ as **weak**, never strong.

## 3. Controls — what clears, what does not (verified, n=10)

The subspace probe runs against the **same mandatory matched-strength gate** as
the per-pair probe (`feedback_matched_strength_mandatory`); for the subspace it
takes the form of a **cluster-mass permutation** over the matched-strength
surrogate ensemble (the gate quantity is `cluster_p_cluster_mass`). The
strong-tier verdict additionally requires the trace to be robust to dropping any
single patient (full-data leave-one-out max p < 0.05) — this is what separates β
and γ_l from δ.

| band | gate `cluster_p_mass` | LOO max p (full) | C5 epi-zone exclusion (secondary) | **subspace verdict** |
|---|---|---|---|---|
| **β** | **0.005** | **0.005** (Pat_02) — robust | strengthens (mass 70 → 89, p=0.005) | **strong** |
| **γ_l** | **0.005** | **0.040** (Pat_05) — robust | contracts (mass 66 → 33, p=0.030, LOO fragile) | **strong** |
| δ | 0.005 | **0.055** (Pat_08) — fails | resolves (mass 38 → 44, p=0.005, LOO→0.005) | **weak** (LOO binds) |
| γ_h | 0.060 | — | — | no trace |
| θ | 0.144 | — | — | no trace |
| α | 0.099 | — | — | no trace |

**Verdict reading.** β and γ_l clear the gate **and** are LOO-robust → strong.
δ clears the *cohort* gate at the empirical floor (`p_mass = 0.005`) but a single
patient (Pat_08) leverages it over the line, so it is **weak**, not strong; the
epi-zone-exclusion analysis is a **secondary mechanistic observation** showing
that leverage is epi-zone-driven (it never promotes the tier). γ_h sits just
outside the gate (0.060), θ and α well outside.

## 4. How the Grassmann result sheds light on the `ρ^coph` reading (the key payload)

This is the reason the subsection exists; give it the most space. The two probes
read **different geometric facets of the same LRG hierarchy** — `ρ^coph` the
displacement-and-recovery of *individual pair* communication distances, the
Grassmann probe the rotation of the *collective leading-mode subspace* — so
agreement and disagreement between them are both interpretable:

| band | `ρ^coph` per-pair (gate p) | Grassmann subspace (gate p_mass) | joint reading |
|---|---|---|---|
| **β** | **trace** (0.005) | **trace** (0.005) | **both probes — primary, convergent** |
| α | **trace** (0.002) | no trace (0.099) | per-pair only |
| γ_l | no trace (0.116) | **trace** (0.005) | subspace only — the inverse of β |
| δ | no trace (0.278) | weak (0.005 / LOO 0.055) | subspace only, weak |
| θ | **anti** (0.722, ρ<0) | no trace (0.144) | reversal, no trace |
| γ_h | no trace (0.246) | no trace (0.060) | silent on both |

Three points to make, in order:

1. **β is the convergent primary.** It is the *only* band whose reorganization is
   coherent at both levels — the individual pair distances are displaced and
   recovered, *and* the leading-mode subspace rotates back. Two probes sensitive
   to different geometry agreeing on one band is what licenses calling β the
   primary trace; this is stronger evidence than either probe alone, and it is
   why β carries the manuscript.

2. **α and γ_l dissociate in opposite directions — and that is the disentanglement.**
   α has a real per-pair trace but no subspace trace: its reorganization is
   **distributed across many pairs without coalescing into a coherent rotation of
   the dominant subspace** — a fine-grained, effectively high-rank change the
   low-rank subspace probe averages away. γ_l is the mirror image: a **coherent
   collective rotation of a band of ~20 slow modes** that does **not** survive as
   consistent per-pair cophenetic displacement — the cophenetic aggregation
   disperses the γ_l per-pair structure that exists at the raw substrate (γ_l is
   "the inverse of β": the cophenet step *demotes* its per-pair signal where it
   *amplifies* β's). δ behaves like a weak γ_l.

3. **Therefore `ρ^coph` silence is not "no reorganization."** The pair-resolution
   probe being quiet in a band (γ_l, δ) means *no reorganization at the per-pair
   level* — the band can still reorganize coherently in its leading subspace, as
   γ_l does. Symmetrically, the subspace probe being quiet where `ρ^coph` speaks
   (α) means α's reorganization is real but not low-rank-coherent. The two
   probes are **complementary windows**: running both is what disentangles, per
   band, **whether** the hierarchy reorganizes and **at what structural scale**.
   This is the concrete content of the "structural enrichment / disentanglement"
   framing for the LRG step — not a single sharper number, but a *resolution of
   the trace into its per-pair and collective components*.

## 5. Guards

- **β is the headline, primary band.** Do not let γ_l's larger panel intensity
  read as primacy: γ_l is **subspace-only**; β's primacy comes from **clearing
  both probes**. State this explicitly if the figure's saturation invites the
  wrong inference.
- **δ is `weak`, never strong** (Decision 12; the CSV/per-band figures are
  stale — see §2 note). Report it as a weak, Pat_08-leveraged subspace trace
  resolved under epi-zone exclusion.
- **No anatomy, no per-patient localization** — retracted and shown delocalized
  this era; keep the subsection cohort-level and structural, never spatial.
- **Do not pool the two probes into a consensus scalar** (forbidden). The point
  is the *pattern of agreement/disagreement across probes*, reported as such.
- **Sign convention**: positive residual / red = trace (post-task subspace
  closer to task than pre-task rest). Never "negative = trace".
- Measure mechanics (principal angles, U_k = span{φ_2..φ_{k+1}}, T_G,
  cluster-extent permutation) go to **Methods**; one back-reference suffices.
- The gate is the **matched-strength cluster-mass** permutation — name it as the
  same mandatory control family as the per-pair probe, so the reader sees one
  battery applied at two structural levels.

## 6. Verified numbers (sources)

Grassmann subspace — `grassmann_cluster_extent/cohort_summary.csv` (audit_70)
+ verdict tiers from `locked/VERDICT_LEDGER.md` (Decision 12, 2026-05-28):

| band | longest run (cells) | window k | `cluster_p_mass` | `cluster_p_LR` | T_G* (norm) | LOO max p | tier |
|------|--------------------|----------|------------------|----------------|-------------|-----------|------|
| δ    | 7  | low-k     | 0.005 | 0.025 | 38.07 (0.149) | 0.055 (Pat_08) | **weak** |
| θ    | 5  | —         | 0.144 | 0.099 | —             | —              | no trace |
| α    | 4  | k=11..14  | 0.099 | 0.159 | —             | —              | no trace |
| β    | 29 | k=27..55  | 0.005 | 0.005 | 69.76 (0.273) | 0.005 (Pat_02) | **strong** |
| γ_l  | 13 | k=12..23  | 0.005 | 0.015 | 66.14 (0.259) | 0.040 (Pat_05) | **strong** |
| γ_h  | 9  | —         | 0.060 | 0.055 | —             | —              | no trace |

Per-pair `ρ^coph` (for the cross-probe table) —
`matched_strength_surrogate_split_baseline/cohort_summary.csv`:
δ +0.008 (p=0.278) · θ −0.040 (p=0.722, anti) · α +0.105 (**p=0.002**) ·
β +0.221 (**p=0.005**) · γ_l +0.083 (p=0.116) · γ_h +0.000 (p=0.246).

C5 epi-zone exclusion (secondary, audit_72) — `grassmann_epi_exclusion/
c5_wilcoxon_cohort.csv`: β mass 69.76 → 89.04 (p=0.005, LOO 0.005); γ_l
66.14 → 32.75 (p=0.030, LOO 0.159 fragile); δ 38.07 → 43.99 (p=0.005, LOO
0.005). Report C5 as a mechanistic observation only — **never a tier-promoter**.

> Note on probe-number provenance: cite the **subspace** verdicts and stats from
> the cluster-extent file + ledger; cite the **per-pair** ρ and gate-p from the
> matched-strength split-baseline file (as the per-pair subsection does). The
> qualitative cross-probe pattern (β both, α per-pair-only, γ_l/δ subspace-only,
> θ/γ_h none) is identical across every source.
