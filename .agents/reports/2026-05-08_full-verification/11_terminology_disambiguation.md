---
type: report
status: current
date: 2026-05-08
era: IMCOH_ABS / COHORT_N10
scope: CRITICAL — the word "trace" carries three distinct senses in §5; rewrite each first occurrence
---

# Terminology disambiguation — "trace" overload

**Head.** The word **trace** is used in §5 (and back-references from §1 and §6) in three structurally distinct senses: (i) a directional cohort-aggregate verdict at substrate or per-pair level (`T_d < 0`, `ρ_split > 0`); (ii) a clade in the §5.6 four-class taxonomy passing the (rspre fragmented, taskt coherent, taskt-rspost coherent) gates; (iii) a contact whose per-leaf score exceeds the within-patient `p95` null in §5.5 anatomy, also called a "trace-positive contact" in some prose passages. The senses are not interchangeable. The proposed vocabulary fix is **trace direction** (the verdict), **trace module** (a §5.6 clade), **trace leaf** (a contact in the §5.5 sense). This is the **CRITICAL** rewrite of the manuscript: any reader who imports the §4 reading of "trace" into §5.5 or §5.6 will misread the per-pair / per-leaf / per-module distinction.

## The three senses, exactly

### Sense 1 — trace direction (a directional verdict)

**What it is.** A cohort-aggregate, scalar-statistic claim about whether the post-task resting state sits closer to the task state than the pre-task resting state does. Lives at the *cohort* level. Defined identically at substrate (§4) and at LRG per-pair correlation level (§5.3).

**Where it appears.**
- §4.2: `T_d := d(taskTest, rsPost) − d(rsPre, taskTest)`. `T_d < 0` is the trace direction.
- §4.4: "The β band reaches seven patients on `d_S` with cohort-median `T_d^{(d_S)} = −0.038`, and is the convergence band on which all three distances agree".
- §4.5 / §4.6: drift controls applied at this level.
- §5.3: `ρ_split(p, b) := ρ_S(Δ_task, Δ_rest)`. `ρ_split > 0` is the trace direction at the LRG per-pair layer.
- §5.3: "Three trace bands α, β, low-γ at within-probe BH-q = 0.027".
- §6.1 (band-resolved synthesis) consistently in this sense.

**Recommended disambiguating phrase.** "Trace direction" or "directional trace pattern" or just "the trace at the cohort level". Avoid the bare word "trace" in this sense unless the surrounding sentence makes the cohort-aggregate, statistical-claim reading unambiguous (e.g., "the trace at α survives the within-baseline null").

### Sense 2 — trace module (a clade in the §5.6 taxonomy)

**What it is.** A subtree of the dendrogram of `D(τ)` for one (patient, band) cell that satisfies the four-class taxonomy match rules:
- `match(rspre, taskt) = fragmented` (Jaccard < 0.5)
- `match(taskt, rspost) = coherent` (Jaccard ≥ 0.6)

Lives at the *per-patient, per-band* level. Counted in cells, not patients.

**Where it appears.**
- §5.6 introductory paragraph: "the four-class taxonomy of cross-phase modules — trace, reset, rearrange, anchor".
- §5.6 Tab 4 row "trace": match status (fragmented, coherent, —).
- §5.6 prose: "rearrange is the most abundant of the four module-level classes by a factor of three to five (560 rearrange clades against 144 anchor, 107 trace, 83 reset across the sixty (patient, band) cells)".
- Fig. 29 (trace class shown at Pat_07 β): "six clades of multiscale sizes (4, 4, 5, 5, 6, 24 leaves; 32 trace leaves of 116 total)".

**Recommended disambiguating phrase.** "Trace module" (or "trace clade", or the §5.6 explicit form "module-level trace class"). The first time §5.6 uses the bare word "trace" it should say "trace module" to anchor the reader.

### Sense 3 — trace leaf (a contact passing the §5.5 per-leaf score gate)

**What it is.** A single contact (electrode channel, leaf of the dendrogram) whose per-leaf trace score `f_i = Σ_j (matching predicate)/|partners|` exceeds the within-patient `p95` of the split-half null. Lives at the *per-contact* level. Counted in leaves (typically 30–75 per band across the cohort).

**Where it appears.**
- §5.5 introductory paragraph: "we resolve this question on the same diffusion substrate as §5.3, at the same `τ' = 1/λ_max`, through a per-leaf calibrated catalog that flags individual contacts by their participation in `ρ_split`".
- §5.5: "We flag a contact as a trace-leaf at band b when its score exceeds the within-patient `p95` of the null distribution".
- §5.5: "the cohort cortical contact pool is `N_total = 866`, with cortical baseline trace rates of 3.46% at α (K = 30 trace-leaves), 5.31% at β (K = 46), and 8.66% at γ_l (K = 75)".
- Also called "trace-positive contact" in the per-leaf decomposition of §5.2 ("A leaf-pair (i, j) is trace-positive when ...").

**Subtle gotcha.** §5.2 defines a per-leaf-**pair** trace-positive predicate `T^+(i, j)` for the KC λ=0 decomposition; this is per-pair-of-leaves, not per-leaf. The "trace-positive" terminology in §5.2 leaks back into §5.5 prose without the leaf-pair vs single-leaf distinction being explicit. §5.5 then aggregates per-pair contributions into a per-leaf score `f_i`, but the language treats "trace-positive contact" and "trace-leaf" as if they were the same object.

**Recommended disambiguating phrase.** "Trace leaf" (or "trace-positive contact", but only after Sense 3 has been pinned to "single contact passing the per-leaf p95 gate"). For the §5.2 per-pair predicate use the explicit phrase "trace-positive leaf-pair" or "trace-positive (i, j) pair", not bare "trace-positive".

## Where the senses get conflated

The conflation risk is highest in three places. These are the rewrite targets.

### Conflation 1 — §5 head paragraph

The §5 head ("Multiscale Laplacian phase-distance investigation") says:

> "the headline result is band-resolved. The per-pair imprint of the cognitive episode survives every control at α, β, and low-γ (q = 0.027). The β imprint registers at every geometric level the analysis tests: per-pair correlation, tree topology, tree heights, partition cut, and spectral subspace. The low-γ imprint registers at the per-pair level and the heights of the dendrogram, with anatomical localization to left fusiform cortex (`p_hyper = 5.2 × 10^−6`). The α imprint lives almost exclusively in the per-pair geometry: it does not restructure the dendrogram tree distance but it shifts the per-pair communication distances coherently, controlled at p = 0.007."

The word "imprint" here is a softer alias for Sense 1 (the directional cohort verdict). The §5 head does **not** use bare "trace" for these. Good. **Action**: keep "imprint" or "per-pair imprint" in the §5 head. Don't substitute "trace" without qualification.

### Conflation 2 — §5.5 anatomy prose

The §5.5 prose says:

> "**The β anatomy is uncorrected-only at every cell.** At full cohort, the three uncorrected candidates are Hippocampus (5/27, 3.49×, p_hyper = 0.011), left fusiform (5/38, 2.48×, p_hyper = 0.045), and left superior temporal (6/53, 2.13×, p_hyper = 0.055), none of which reaches the Bonferroni threshold."

These "trace-leaf" counts (5/27, 5/38, 6/53) are Sense 3. The same paragraph then talks about "the β anatomical headline" and later about "the β imprint" — here "imprint" is Sense 1 (cohort directional verdict from §5.3). The reader has to track that the *β imprint* (Sense 1) does **not** localize anatomically (no Bonferroni-surviving region), but β *trace-leaves* (Sense 3) exist and are spread across multiple regions. **Action**: rewrite the first occurrence of each in §5.5 to say "trace-leaf" and "β imprint at the per-pair correlation level (§5.3)" explicitly.

### Conflation 3 — §5.6 first paragraph

§5.6 first paragraph says:

> "the four-class taxonomy of cross-phase modules — trace, reset, rearrange, anchor".

A reader coming from §5.3 will read "trace" as Sense 1 (the directional cohort verdict). It is actually Sense 2 (a per-cell module-level class). **Action**: rewrite as:

> "the four-class **module-level** taxonomy of cross-phase modules — **trace module**, reset module, rearrange module, anchor module".

Then the rest of §5.6 can drop the "module" suffix once the disambiguation is established.

## Suggested rewrites — per first occurrence

### §1 closing paragraph (preview)

Currently: "controlled cohort claim at α, β, and low-γ via the per-pair correlation on `D(τ)`".

Keep "trace" out of this sentence; the sentence is already explicit ("controlled cohort claim", "per-pair correlation"). No rewrite needed at §1.

### §4.2 first occurrence

Currently:
> "A negative `T_d` indicates that the post-task resting state sits closer to the task state than the pre-task resting state does, which is the directional signature of a *trace*: the cognitive episode left a structural mark that the network carries into the post-task rest."

This is fine as-is — the italics and the "directional signature of a *trace*" frame it as Sense 1 explicitly. No rewrite needed.

### §5.2 KC decomposition

Currently:
> "A leaf-pair (i, j) is trace-positive when ..."

**Rewrite first occurrence**: "A leaf-pair (i, j) is **trace-positive at the leaf-pair level** when ... (the trace-positive leaf-pair predicate is per-pair, distinct from the trace-leaf predicate of §5.5)".

After this anchor, "trace-positive" can stay bare in §5.2.

### §5.3 first occurrence after the head

Currently:
> "Two further controls operate on the same per-(patient, band) cell: the drift-floor null `ρ_drift(p, b)` ... captures the within-session drift component of the per-pair correlation".

The §5.3 head already established "ρ_split > 0 is the trace direction" — keep "trace direction" here, not bare "trace".

### §5.5 first occurrence

Currently:
> "We flag a contact as a trace-leaf at band b when its score exceeds the within-patient p95 of the null distribution".

This is the canonical Sense 3 anchor. Keep as-is — "trace-leaf" with the hyphen is the right disambiguator.

But the *next* sentence:
> "After dropping white-matter (Wm) and parcellation-failure (Unk) labels, the cohort cortical contact pool is N_total = 866, with cortical baseline trace rates of 3.46% at α (K = 30 trace-leaves)..."

reads "trace rates" — keep "trace-leaf rates" (rate of trace-leaves in the cortical pool) not bare "trace rates". One word substitution.

### §5.5 anatomy headline prose

Currently:
> "The β anatomical reading is therefore not a single-patient artifact in disguise: Pat_02 is both the cohort's anatomically densest fusiform implant and its highest-rate fusiform tracer".

"Highest-rate fusiform tracer" — Sense 3 (Pat_02 has the highest per-patient rate of trace-leaves in fusiform). **Rewrite**: "highest per-patient rate of trace-leaves in fusiform". Avoid the gerund "tracer" which is undefined.

### §5.6 first paragraph

See **Conflation 3** above. Rewrite "the four-class taxonomy of cross-phase modules — trace, reset, rearrange, anchor" to make the module-level reading explicit on first occurrence.

### §5.6 trace-class prose

Currently (§5.6 paragraph after Tab 4):

> "Fig. 29 shows trace at Pat_07 β, where six clades of multiscale sizes (4, 4, 5, 5, 6, 24 leaves; 32 trace leaves of 116 total) fragment in rspre, crystallize coherently in taskt, and persist into rspost".

Here "trace" appears twice in different senses: "trace at Pat_07 β" (Sense 2 — module class) and "32 trace leaves of 116 total" (Sense 3 — leaf count, but here "leaves" of a clade not "trace-leaves" of the cortical pool — yet another local meaning).

**Rewrite**: "Fig. 29 shows the **trace module** at Pat_07 β, where six clades of multiscale sizes (4, 4, 5, 5, 6, 24 leaves; **32 leaves of 116 belong to a trace module** of total) fragment in rspre, crystallize coherently in taskt, and persist into rspost". The leaves-of-a-clade reading is now distinct from the "trace-leaf" of §5.5.

Note: Pat_07 β taxonomy "trace leaves" (32 leaves in trace modules) is **NOT** the same as Pat_07 β §5.5 "trace-leaves" (cortical contacts with `f_i > p95`). At Pat_07, focal leaf `i*=62` and `|P|=37` trace-positive partners (§5.2) — yet again a different per-pair count. These three numbers (32 / 37 / and any §5.5-style p95 trace-leaf count for Pat_07 β cortical contacts) coexist in the same paragraph block and must not be conflated.

### §6.1 synthesis

Currently:
> "Anatomically, no single region carries it under multi-region correction: the strongest full-cohort candidate (Hippocampus, 3.49×) does not survive Pat_03 dropout, and the only Pat_03 -robust candidate (left superior temporal, 2.13×) remains uncorrected."

"Carries it" — refers back to "the β imprint" (Sense 1). Acceptable as-is because §6.1 is the synthesis and the antecedent is unambiguous. But the phrase "the data test single-region concentration rather than distributional spread" two sentences later is hand-wavy; Pat_02 fusiform double-density disclosure should be repeated as a one-liner here for self-containment.

## Action — writing agent

1. Adopt the three-term vocabulary: **trace direction**, **trace module**, **trace leaf**.
2. On the **first** occurrence of each in **each** of §5.2, §5.3, §5.5, §5.6, qualify with the exact phrase ("trace direction", "trace module", "trace leaf"). Subsequent occurrences in the same subsection can drop the qualifier.
3. The bare word "trace" remains acceptable in the §1 / §4 prose because §4 always means Sense 1 and §1 only previews. But §5 onwards should never use bare "trace" without one of the three qualifiers.
4. Update the abbreviation glossary if there is one (current `acronym` list does not contain "trace" — add it as a three-row entry: trace direction / trace module / trace leaf).
5. Tab 4 (the §5.6 four-class table): the "trace" row should read "trace **module**" (or "module-level trace") for self-containment.
6. Figure captions of Figs. 29 / 32 / 34 (the §5.6 class composites): reads "trace" / "anchor" labels in the panel headers — these are class names, leave as-is; but the caption prose should explicitly say "trace module" on first reference.

## Connection — `.agents/guides/01_project/terminology.md`

The terminology guide already establishes the **trace / anchor / reset / emergent** taxonomy at the agent level (and the "anchor" / "reset" / "rearrange" terms used in §5.6 share the spirit of that taxonomy). The manuscript should be consistent with the guide. The §5.6 four classes are: trace, reset, rearrange, anchor. The agent terminology guide uses: trace, anchor, reset, emergent. Note that **emergent** (a clade present in rspost that did not exist in rspre) is **not** the same as **rearrange** (post-task coherence with no antecedent at any earlier phase). They are similar concepts; the manuscript uses "rearrange" (to keep "emergent" available for the LRG-only case where a community membership was created by the diffusion process). This distinction is acceptable but should be flagged in §5.6 with a one-liner: "we use 'rearrange' for the §5.6 module-level class to distinguish it from the LRG-emergent partition-membership concept of `terminology.md` §3."
