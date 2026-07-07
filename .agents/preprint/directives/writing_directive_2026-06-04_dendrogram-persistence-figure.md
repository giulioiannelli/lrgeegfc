---
name: writing-directive-dendrogram-persistence-figure
era: IMCOH_ABS_COHORT_N10
status: current
kind: directive
scope: Fill the single [PENDING] illustrative figure + one-paragraph slot in the cophenetic subsection — the circular-dendrogram "persistence" triptych. Tells the writing agent which file to use, how the colour works, how to read it, and the hard guardrails. ILLUSTRATIVE single-patient/single-band visual only; changes no locked verdict and makes no new claim.
created: 2026-06-04
companion: directives/writing_directive_2026-06-04_coph-subsection-results.md (the cohort result this figure illustrates — the real claim lives there, not here).
figure_source: data/preprint/figures/beta/dendrogram_persistence/fig_beta_dendrogram_persistence_Pat_05_taskref_cont_sp.pdf (generator: scripts/02_preprint/preprint_15_dendrogram_persistence_triptych.py --patient Pat_05 --band beta --reference task --gate-mode continuous)
---

# Dendrogram-persistence triptych — illustrative figure + paragraph

**Head.** Fill the `[PENDING]` figure/paragraph slot with the circular-dendrogram
triptych for \patient{05}, \(\beta\). It is a **visual intuition** for what the
cophenetic persistence (ρ^coph trace) *looks like* — one patient, one band — not
a per-patient result. The reading is a **panel-level contrast**: the \acrshort{rspost}
tree is broadly more task-coloured than the \acrshort{rspre} tree.

> **Two corrections to the LaTeX comment — apply both:**
> 1. **File changed.** Use
>    `figures/fig_beta_dendrogram_persistence_Pat_05_taskref_cont_sp.pdf`,
>    **not** `..._taskref_excl.pdf`.
> 2. **Framing changed.** The comment says "shared baseline removed; persistence =
>    branches retained in BOTH rest phases." That is the *old* (excl) plan and is
>    **wrong** for this figure. The chosen render does **no** baseline removal and
>    is **not** an intersection. Use the continuous-shade reading below.

---

## 1. What the figure is (keep it simple)

Three circular dendrograms side by side — \acrshort{rspre}, \acrshort{taskt}
(the **reference**, marked with a star ★), \acrshort{rspost} — for one patient,
one band.

- The **task** panel sets a continuous rainbow over its leaf order; both rest
  panels are re-laid-out to match that order, so angular position is comparable
  across panels.
- In each **rest** panel, every leaf is coloured by its **per-leaf cophenetic
  similarity to task**, ρ^coph_c = the rank (Spearman) agreement between that
  contact's row of communication-distances-to-all-others in task vs in that rest
  phase. High ρ^coph_c → that contact sits in the **same overall place** in the
  hierarchy as during task; low/negative → a different place.
- **Colour = saturation.** Most task-similar leaf → full hue; least similar →
  faded to light grey. The scale is **shared** across the two rest panels, so a
  broadly-less-similar phase simply reads greyer. No threshold, no cut.
- Branch colour is the size-weighted mean of its leaves' colours, so a coherent
  task-like cluster shows as a coloured wedge and a scrambled one greys out.

## 2. How to read it (the one sentence the paragraph must land)

**\acrshort{rspost} is broadly more saturated than \acrshort{rspre}**: the
\(\beta\) communication hierarchy present during \acrshort{taskt} is broadly
**retained into \acrshort{rspost}** and was fainter / less task-like in
\acrshort{rspre}. That brightness contrast *is* the visual sense of the
cophenetic persistence. (The task panel is fully coloured by construction — it
is the reference.)

Illustrative numbers you may quote for **this patient only**:
ρ^coph(\acrshort{taskt}, \acrshort{rspre}) ≈ **+0.30**,
ρ^coph(\acrshort{taskt}, \acrshort{rspost}) ≈ **+0.48** — the \acrshort{rspost}
tree is the more cophenetically task-like of the two.

## 3. Hard guardrails (do not violate)

- **Illustrative only.** Single patient (\patient{05}), single band (\(\beta\)).
  Frame as "what the persistence looks like." The actual claim is the **cohort**
  ρ_split^coph result in the companion subsection — point there for the result;
  this figure carries no claim of its own. Subsection stays cohort-level.
- **ρ^coph_c is an *overall*-position measure, not "this leaf stayed put."** It
  is dominated by the bulk of far-apart pairs, so a coloured leaf can still sit
  in a different local cluster. Therefore:
  - **No node-by-node narration**, no naming individual contacts as "persistent."
  - **No anatomy / localization** read off this figure (the multi-region DK lists are
    retracted; β's one audited localization — the OFC *system*, `ANATOMY_LEDGER.md` 2026-06-10 —
    is a separate system-scale analysis, not read off this dendrogram).
  - **Do not** describe a *sparse set* of persistent nodes. The colour is broad
    and graded; that breadth is exactly the delocalized nature of the trace.
- **Do not redefine the measure here** — cophenetic distance / ρ^coph belong to
  Methods. One back-reference suffices.

## 4. Gaps to fill

- **`\caption{[PENDING]}`** — plain-language: name the three panels and the ★
  reference; state colour = per-leaf cophenetic similarity to task, faded to
  light grey on a shared scale; give the reading (rspost more task-coloured than
  rspre). 2–3 sentences, no in-figure-claim language.
- **One short paragraph** before/after the float — the §2 reading, explicitly
  tagged illustrative, with the back-reference to the cohort result. Keep
  `\label{fig:beta_dendrogram_persistence}`.
