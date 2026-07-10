---
name: methods-directive-TG-normalization
filename_note: "Filename uses `methods_directive_` prefix rather than the WRITING_GUIDE.md `writing_directive_` convention because this is a methods-side directive (changes a Methods-section equation), not a writing-agent task. Kept verbatim under `status: applied` 2026-05-28 — not renamed to avoid breaking inbound references."
era: IMCOH_ABS_COHORT_N10
status: applied
applied_date: 2026-05-28
kind: methods-agent-directive
date: 2026-05-19
target: methods §sssec:methods_compare_grassmann — normalize the cluster mass T_G^*
priority: high (changes a methods-equation definition and a primary scalar)
source_of_truth:
  - .agents/preprint/directives/methods_grassmann_cluster_extent.md
  - .agents/preprint/directives/writing_directive_2026-05-19_beta_post_methods_revision.md
  - data/audit/grassmann_cluster_extent/cohort_summary.csv
applied_in:
  - .agents/preprint/METHODS_AUDIT_ISSUES.md (C1 lock with normalized values table)
  - .agents/preprint/directives/methods_grassmann_cluster_extent.md (Head + §5b normalized formula + §6 + §9 checklist updated 2026-05-28)
  - .agents/preprint/bands/00_cohort.md §3 (normalized T_G^* column added 2026-05-28)
  - .agents/preprint/bands/03_gammalow.md (γ_l raw 66.14 / normalized 0.259 cited 2026-05-28)
  - .agents/preprint/bands/06_delta.md (δ raw 38.07 / normalized 0.149 cited 2026-05-28)
  - .agents/preprint/locked/VERDICT_LEDGER.md (Decision 12 cascade 2026-05-28; δ verdict flipped strong → weak)
notes_on_application: |
  The "Open item" (line 64 of body) — discrepancy between cohort_summary CSV
  (longest contiguous cluster only) and methods companion equation (sums over
  all sig cells) — was RESOLVED 2026-05-19 pm: the CSV is now all-clusters
  formula. Post-fix raw values: β 69.76 / γ_l 66.14 / δ 38.07 / γ_h 31.67 /
  θ 12.97 / α 7.98. Normalized values (denominator 255.65 for full-data
  cohort with n_k^cohort = 111 and R = 200): β 0.273 / γ_l 0.259 / δ 0.149 /
  γ_h 0.124 / θ 0.0507 / α 0.0312.

  The body text below references pre-fix values (β = 52.97 / normalized
  0.2073) and the un-resolved "Open item" — kept verbatim as historical
  record. **Live citations should use the post-fix values listed under
  `notes_on_application` above, or pull directly from
  `data/audit/grassmann_cluster_extent/cohort_summary.csv`.**
---

# Methods directive — normalize the Grassmann cluster mass T_G^* (2026-05-19)

**Head.** The cluster mass `T_G^*` defined in Methods Eq.~\eqref{eq:methods_TGstar} currently produces unbounded raw values that are hard to read across patients and bands. **User directive 2026-05-19**: redefine `T_G^*` as the **normalized form**, divided by the maximum possible mass given the empirical-null floor, so it lives in `[0, 1]` and is naturally reported as a percentage. **The normalized form is the only form**: no unnormalized version is reported anywhere in Methods, Results, captions, or supplementary material. The symbol `T_G^*` (no tilde, no superscript-tilde, no special decoration) refers to the normalized cluster mass throughout. The cohort gate (`p_mass < 0.05` on the empirical-null permutation) is **unchanged** — only the definition of the scalar changes.

The same normalization applies at both levels:
- **Band-level (cohort)** scalar `T_G^*(b) ∈ [0, 1]` (Eq.~\eqref{eq:methods_TGstar}).
- **Per-patient** scalar `T_G^{*, p}(b) ∈ [0, 1]` (new equation, see §3 below).

Both share the same denominator `K · log10(R+1)`, so they are directly comparable.

## What changes — one equation, one paragraph

### Eq.~\eqref{eq:methods_TGstar} → redefine as the normalized cluster mass (band level)

**Replacement (this IS the new definition of `T_G^*`; no separate "tilde" form):**
```
T_G^*(b) = ( Σ_{k : p_k(b) < α_k} (−log10 p_k(b)) ) / [ K · log10(R+1) ]
         ∈ [0, 1]
```

with `K = |k-grid| = 111` (`k ∈ {2, …, 112}`) and `R = 200` the matched-strength surrogate ensemble size. `−log10(1/(R+1)) = log10(R+1) ≈ 2.303` is the maximum contribution any single `k`-cell can make under the `1/(R+1)` empirical-null floor; multiplying by `K` gives the maximum possible (unnormalized) cluster mass, `K · log10(R+1) ≈ 255.65`. Dividing by this denominator forces `T_G^*(b) ∈ [0, 1]`. `T_G^*(b) = 1` ⇔ every `k`-cell in the grid hits the empirical-null floor (maximum trace evidence); `T_G^*(b) = 0` ⇔ no `k`-cell clears `α_k`.

**Notation lock (user directive 2026-05-19).** Use the symbol `T_G^*` (no tilde, no `\widetilde`, no `\check`) for the normalized form. The unnormalized log-sum is **not** a quantity the manuscript carries — neither as an intermediate symbol nor as a value in any table or figure. Eq.~\eqref{eq:methods_TGstar} replaces the prior unnormalized equation directly; any prior reference to `T_G^*` in the manuscript is automatically the normalized form after this update.

### Phantom-surrogate null — same normalization on the null distribution

The empirical null `mass^null(r)` in §5c (lines 246 of the methods companion) must be normalized **by the same denominator** for the permutation p-value to remain calibrated:

```
mass_normalized^null(r) = mass^null(r) / [ K · log10(R+1) ]
```

The empirical p-value `cluster_p_mass(b)` is **unchanged** under simultaneous normalization of observation and null — both sit in `[0, 1]` and the ranking is preserved. The methods text should state this explicitly:

> "Normalizing both `T_G^*` and `mass_normalized^{null}_r` by the same denominator `K · log10(R+1)` leaves the empirical p-value `p_mass(b) = (1 + #{r : mass_normalized^{null}_r ≥ T_G^*(b)}) / (R + 1)` unchanged."

### Verdict gate — unchanged

`p_mass(b) < 0.05` remains the gate. The β verdict `strong trace` at `p_mass = 0.005` survives unchanged because the normalization is monotone.

## Numerical anchors (post-normalization)

For β at the manuscript-window definition of cluster mass (sum over the longest contiguous-significant cluster `k ∈ [27, 55]`, 29 cells) the current `obs_cluster_mass_neglog10p = 52.97` stored in `data/audit/grassmann_cluster_extent/cohort_summary.csv`:

- Raw: `T_G^*(β) = 52.97`
- Normalized: `T_G^*(β) = 52.97 / 255.65 = 0.2073` ≈ **20.7%** of maximum possible cluster mass

If the methods equation is read literally (sum over **all** significant `k`-cells, not just the longest contiguous cluster), the raw cluster mass for β is 69.76 over 40 cells, giving `T_G^*(β) = 69.76 / 255.65 = 0.273` ≈ **27.3%**.

**Open item** — see `writing_directive_2026-05-19_beta_post_methods_revision.md` for the standing discrepancy between the cohort_summary CSV (sums over the longest contiguous cluster only) and the methods companion equation (sums over all sig cells). Methods agent should pick **one** and lock it; both numbers above are computed for cross-reference.

## Per-band cluster-mass values to update in any table the methods section carries

If the methods section cites raw cluster masses anywhere (e.g. comparison across bands), update to the normalized form. Using the cohort_summary.csv stored values (manuscript-window convention):

| Band | Raw `T_G^*` | Normalized `T_G^*` | Percentage | `p_mass` |
|---|---|---|---|---|
| **β** | **52.97** | **0.207** | **20.7%** | **0.005** |
| γ_l | 19.17 | 0.075 | 7.5% | 0.035 |
| δ | 12.78 | 0.050 | 5.0% | 0.025 |
| α | 7.98 | 0.031 | 3.1% | 0.099 |
| γ_h | 15.75 | 0.062 | 6.2% | 0.065 |
| θ | 7.11 | 0.028 | 2.8% | 0.144 |

The relative ordering of bands is preserved (β >> γ_l > γ_h > δ > α > θ). The normalized scalar now reads as "fraction of maximum possible trace evidence in the `k`-spectrum" rather than as an unbounded log-sum.

## What the methods agent updates

1. **Eq.~\eqref{eq:methods_TGstar}** — replace with the normalized form (or attach the normalization as a separate equation immediately after).
2. **Paragraph around Eq.~\eqref{eq:methods_TGstar}** — add one sentence defining the denominator `K · log10(R+1)` and noting that the empirical-null permutation remains calibrated under simultaneous normalization.
3. **Add the per-patient cluster mass `T_G^{*, p}` to the methods section.** Currently the per-patient cluster mass is defined only in the β table caption (`overleaf/tables/beta_per_patient.tex` + `.md`). Move that definition **out of the table caption and into the methods section**, immediately after the band-level `T_G^*` definition. Suggested placement: a new short subsubsection or a paragraph titled "Per-patient cluster mass" inside §sssec:methods_compare_grassmann (just after the cluster-extent / `T_G^*` paragraph). Equation to add:
   ```
   T_G^{*, p}(b) = ( Σ_{k : p_k^{p}(b) < α_k} (−log10 max(p_k^{p}(b), 1/(R+1))) )
                          / [ K · log10(R+1) ]
   ```
   with `p_k^{p}(b)` the per-patient one-sided empirical p-value at cutoff `k` against that patient's `R = 200` matched-strength surrogate distribution (definition in §sssec:methods_compare_stats), regularized at the same `1/(R+1)` floor so single-cell contributions are bounded by `\log_{10}(R+1)`. Sign convention: `T_G^{*, p}(b) ∈ [0, 1]`, monotone in trace strength.
   
   The per-patient cluster mass is **not a separate test** — it is the per-patient analogue of the band-level cluster mass, reported as descriptive per-patient evidence in the supplementary per-patient table (Table~\ref{tab:beta_per_patient}). The verdict gate remains the band-level `p_mass < 0.05` on `T_G^*`. The methods section should state this relationship explicitly so a reviewer does not read the per-patient values as a separate Wilcoxon family.
4. **`5d Mass-only verdict gate`** in the methods companion — state explicitly that the gate is on `T_G^*` (or `T_G^*` if the symbol is retained for the normalized form), unchanged in mechanics.
5. **Any table in the methods section that quotes `T_G^*` values** — replace with normalized values per the per-band table above. (β-only, this is just one value: `T_G^*(β) = 0.207` / 20.7%, `p_mass = 0.005`.)
6. **`L_obs` is unchanged** — the longest contiguous-significant run statistic is methods-side descriptive only (results paragraph does not cite it, per the 2026-05-19 writing directive); no normalization needed for `L_obs`.

### Definition that the methods agent inherits from the writing-side (so the β table caption can be slimmed)

The current β table caption (`overleaf/tables/beta_per_patient.tex`) carries a full definition of the per-patient cluster mass — this content **should move to Methods**. The table caption then references the methods definition only:

> "Cols.\ 6--7: per-patient normalized cluster mass \(T_G^{\ast, p}\) (see Methods Eq.~\eqref{eq:methods_TGstar_perpatient}) and the count \(n_{\rm sig, k}\) of \(k\)-cells passing the per-cell threshold \(\alpha_k = 0.05\) at the patient level."

This slims the caption to two lines and gives the reader a single canonical location (Methods) for the definition. The writing-side directive will be updated in parallel to slim the table caption once the methods definition lands.

## What does **not** change

- The cluster-extent permutation algorithm (phantom-surrogate construction, R = 200).
- The per-`k` cohort-paired Wilcoxon test (`p_k(b)`).
- The verdict gate `p_mass < 0.05`.
- The verdict for any band (β strong, γ_l + δ weak, α + γ_h + θ none).
- The methods companion's broader content (subspace definition, chordal distance, per-`k` statistic).

## Where to file the methods edit

`.agents/preprint/directives/methods_grassmann_cluster_extent.md` (the locked methods companion). The user's actual LaTeX methods section (`ssec:methods_compare`) currently quoted in the manuscript should be edited in step with the companion.

## Revision history

- **2026-05-19** — Initial directive. Normalize Grassmann cluster mass `T_G^*` to `[0, 1]` by dividing by `K · log10(R+1) ≈ 255.65`; preserve verdict gate and band ordering; β normalized value 0.207 (20.7% of max), `p_mass = 0.005` unchanged.
