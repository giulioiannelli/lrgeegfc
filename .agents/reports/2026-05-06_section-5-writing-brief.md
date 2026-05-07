---
name: 2026-05-06_section-5-writing-brief
type: report
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-05-06
updated: 2026-05-06
pointers:
  - .agents/reports/2026-05-05_result-2-lrg-beta-trace.md
  - .agents/reports/2026-05-05_section-5-manuscript-draft.md
  - .agents/reports/2026-04-29_result-1-raw-fc-phase-trace.md
  - data/reports/section_5_lrg_trace/README.md
  - data/reports/section_5_lrg_trace/headline/README.md
---

# Section 5 — writing-agent brief

**Purpose.** This is what to give a writing agent so they can produce the
manuscript Section 5 on the LRG-layer β-band trace finding without
re-deriving the analysis or re-discovering the framing. Self-contained
read-order + final framing decisions + do/don't claims + provenance.

---

## 0. Read order (do not skip)

1. **Plain-language synthesis (start here):**
   `.agents/reports/2026-05-05_result-2-lrg-beta-trace.md` —
   "In plain words" section first, then the technical evidence chain.
2. **Drop-in manuscript prose (already written, edit-don't-redraft):**
   `.agents/reports/2026-05-05_section-5-manuscript-draft.md` — six
   paragraphs (head, methods, β results, low_γ + cross-validation,
   anatomy, limitations) plus four figure caption stubs.
3. **Section 4 reference (so Section 5 builds, not repeats):**
   `.agents/reports/2026-04-29_result-1-raw-fc-phase-trace.md`.
4. **Per-measure index (only when checking a specific claim):**
   `data/reports/section_5_lrg_trace/README.md`. Each row links to a
   self-contained measure report.
5. **Headline figures + figure-purpose guide:**
   `data/reports/section_5_lrg_trace/headline/README.md` describes
   what each of the four figures proves and which reviewer question
   it answers.

If the writing agent reads only one file, it's #1. If two, #1 + #2.

---

## 1. Final framing decisions (locked)

### 1a. Lead with the dendrogram-amplification finding

The most novel methodological contribution of Section 5 is **the LRG
dendrogram amplifies the same underlying signal ~67× over matrix
distance at β** (effect 5.3 vs 0.08 for the same Wilcoxon p=0.001). This
is what the volcano figure shows, and it justifies why the LRG layer is
needed. Without LRG, β would look like a small effect (0.08) easy to
dismiss as drift; with LRG, the effect is unmistakable (≈5).

→ The volcano figure (`headline/figures/volcano.pdf`) is Figure 5A and
should be cited in the methods + results.

### 1b. β is "load-bearing", low_γ is "secondary via different route"

Avoid "β is the only band". The honest framing:
- **β** has 4 cells surviving joint Bonferroni m=48 (all KC variants +
  d_F) — the *partition-structure* route.
- **low_γ** has 2 cells surviving (Grassmann k=13 global + per-leaf
  controlled localization) — the *eigenstructure* + *localization*
  route.
- **α / δ / θ / high_γ** show no global-probe trace under controls.

The two-band picture is more defensible than a one-band claim and
matches the data.

### 1c. Anatomy = Hippocampus + left fusiform (sampling-corrected)

After hypergeometric correction for cohort-wide contact density:
- β: **Hippocampus** (3.5×, p=0.011, 3/9 patients) + **left fusiform
  cortex** (2.5×, p=0.045, 2/9). Left superior temporal borderline.
- low_γ: **left fusiform** (4.0×, p<0.001) + **left inferior temporal**
  (2.3×, p=0.026).
- **Left middle temporal at β is NOT a finding** — at baseline rate
  (1.1×, p=0.48). The original raw-count interpretation conflated
  signal with implant-sampling density. This was caught and corrected
  on 2026-05-06; the writing agent must NOT revert to "left temporal
  lobe" generically.

→ The anatomy figure (`headline/figures/anatomy_2d.pdf`) is Figure 5D.

### 1d. Result 2 = structural enrichment of Result 1, not independent confirmation

Same data, different geometric layer. Section 5 should explicitly say
"this is not independent replication; it is a methodologically distinct
view of the same imaginary-coherency matrices." The reviewer-defense
move is in the *limitations paragraph* of the manuscript draft.

The novel content of Section 5 over Section 4: partition-structure
shifts (LRG-layer) vs edge-weight shifts (raw-FC). Different question,
related answer.

### 1e. Voice and tone

- Third person plural, past tense. NeuroImage / J Neurosci tier.
- Renormalization-style: lead each paragraph with a 1-sentence head, fill
  in technical detail below. Per the project never/always rule.
- Use the **trace** terminology, not "persistence" (per
  `feedback_trace_terminology.md` memory + `terminology.md` guide).
- Variable names: `T_d`, `n_trace`, never `n_persist`.

---

## 2. Defensible claims (cite-with-confidence list)

| claim | numerical anchor | source CSV |
|:--|:--|:--|
| β all 4 facets pass joint Bonferroni m=48 | p=0.000977, 10/10 patients each | `data/audit/lrg_global_probe_controls/cohort_controls_summary.csv` |
| β KC λ=0 Pat_03-dropout robust | p=0.019 → 0.037 (still significant) | `data/audit/lrg_global_probe_controls/pat03_dropout.csv` |
| Dendrogram amplifies effect 67× over D-rank at β | KC β λ=0 effect 5.29 / d_F effect 0.08 | (computed in volcano figure code) |
| 5/6 bands have NO uncorrected global-probe trace | only β + low_γ have any cell with p<0.05 | `cohort_controls_summary.csv` (filter by p<0.05) |
| low_γ has 1 joint-Bonferroni-surviving global cell | Grassmann k=13 p=0.000977 | same CSV |
| 9/10 patients have ≥3 calibrated trace-leaves at β | per-measure Bonferroni m=6 surviving (p=0.002) | `data/audit/per_leaf_rho_null/cohort_band_summary.csv` |
| Bidirectional cross-validation at β | enrichment p=0.004, concentration p=0.004 | `data/audit/lrg_cross_validation/cohort_band_summary.csv` |
| β trace enriched at Hippocampus 3.5× | hypergeometric p=0.011, 3/9 patients | `data/audit/lrg_localization_anatomy/region_enrichment.csv` |
| β trace enriched at left fusiform 2.5× | hypergeometric p=0.045, 2/9 patients | same |
| low_γ trace enriched at left fusiform 4.0× | hypergeometric p<0.001, 3/10 patients | same |

Every number above is in a CSV on disk. The writing agent should verify
each one before submission.

---

## 3. Do NOT write (claims that look right but aren't)

| claim | why it's wrong |
|:--|:--|
| "β is the only band with a trace" | low_γ has joint-MTC-surviving evidence too. |
| "Left temporal lobe + Hippocampus dominate β" | After sampling correction, left middletemporal is at baseline (p=0.48). Be specific: Hippocampus + left fusiform. |
| "Section 5 is independent replication of Section 4" | Same data, different geometric layer. Explicitly NOT independent. |
| "The β trace lasts hours / persists indefinitely" | We measured a fixed post-task rest window. Duration is unknown. |
| "Memory localizes to β-band" | We show network signatures, not behavior. No behavioral data is in this paper's analysis. |
| "Cohort is N=10 healthy adults" | Cohort is N=10 sEEG epilepsy patients. Clinical context matters. |
| "All probes agree" | Grassmann k=13 at β is uncorrected only (p=0.024 BH-FDR, not Bonferroni). β trace is partition-structure, NOT eigenstructure. State this as a methodological distinction, not a failure of agreement. |
| "Trace lateralizes to dominant hemisphere" | Implant cohort is 64.5 % LH by clinical sampling design. The LH skew of trace-leaves is consistent with that, not evidence of lateralization. Would need handedness/Wada-test data. |
| "Effect size 0.08 (D-rank)" or "Effect size 5 (KC)" without context | These two numbers describe the SAME underlying effect viewed through different probes. The 67× amplification is the headline; reporting either alone is misleading. |
| "Pat_03 drives the β finding" | KC λ=0 survives Pat_03 dropout (p=0.037). Pat_03 contributes the largest |T_d| (-19.7) but the cohort signal does not depend on it. |

---

## 4. Joint MTC numbers (cite carefully)

| scope | m | Bonferroni threshold | surviving cells |
|:--|:-:|:--:|:--|
| controls (measure 11) only | 42 | 0.00119 | 5 (β: KC λ=0/0.5/1 + d_F; low_γ: Grassmann k=13) |
| controls + per-leaf (measures 11+12) | **48** | **0.00104** | **6** (above + low_γ per-leaf) |
| controls + per-leaf + cross-validation (11+12+09) | 60 | 0.000833 | 5 (β cells + low_γ Grassmann); per-leaf low_γ misses by 0.00014 |

The manuscript headline number is **m=48, 6 cells surviving**. Cross-
validation is reported as methodological consistency, not statistical
power.

---

## 5. The four headline figures and their roles

| figure | reviewer question it answers | path |
|:--|:--|:--|
| 5A volcano (3-panel) | "Is β really unique vs other bands?" | `headline/figures/volcano.pdf` |
| 5B real-vs-null scatter | "Is the result driven by 2-3 patients?" | `headline/figures/per_patient_slopes.pdf` |
| 5C facet-agreement matrix | "Is it specific to one geometric measure?" | `headline/figures/cohort_fingerprint.pdf` |
| 5D anatomy enrichment | "Where in the brain?" | `headline/figures/anatomy_2d.pdf` |

Caption stubs already drafted in the manuscript draft (Section 5
manuscript draft, "Figure caption stubs" subsection).

---

## 6. Result 1 → Result 2 narrative arc

The two-paragraph version that should anchor the manuscript flow:

> Section 4 (raw-FC) showed: in β-band imaginary coherency, 7/10
> patients have d_S^β triangle T_d < 0 (controls pass) — task-shaped
> reorganization in resting-state edges, persisting into rest_post.
>
> Section 5 (LRG) shows: at the LRG dendrogram layer, β reorganization
> is much stronger and shows up across multiple geometric facets (KC
> tree distance λ=0/0.5/1, Frobenius matrix distance d_F), with effect
> sizes 67× larger than at the matrix-distance layer. The β trace is a
> *partition-structure* shift, not just an edge-weight shift. After
> sampling correction, the trace is anatomically enriched at the
> Hippocampus and left fusiform cortex — known memory-system regions
> consistent with the patients having performed a memory-encoding task.

Section 5 builds on Section 4 by saying "the same data also carries a
multiscale-partition-level signal that the edge-level test alone could
not have revealed."

---

## 7. Numbers to triple-check before submission

The writing agent should run this exact snippet at the bottom of the
manuscript draft and verify the numbers it prints against the cited
prose:

```python
import pandas as pd
ROOT = "data"

ctl = pd.read_csv(f"{ROOT}/audit/lrg_global_probe_controls/cohort_controls_summary.csv")
beta = ctl[(ctl.band == "beta") & (ctl.wilcoxon_p_real_lt_null < 0.05/len(ctl))]
print("β Bonferroni m=42 surviving:", beta[["probe", "variant", "n_real_below_null",
                                            "wilcoxon_p_real_lt_null"]].to_string(index=False))

p3 = pd.read_csv(f"{ROOT}/audit/lrg_global_probe_controls/pat03_dropout.csv")
print("\nPat_03 dropout:", p3[["probe", "variant", "wilcoxon_p_full",
                                "wilcoxon_p_drop"]].to_string(index=False))

anat = pd.read_csv(f"{ROOT}/audit/lrg_localization_anatomy/region_enrichment.csv")
beta_anat = anat[(anat.band == "beta") & (anat.enrichment >= 2.0) &
                 (anat.p_hyper < 0.05) & (anat.n_pat_trace >= 2)]
print("\nβ enriched regions:", beta_anat[["region", "enrichment", "n_pat_trace",
                                          "p_hyper"]].to_string(index=False))
```

Expected output: 4 β cells at m=42, KC β λ=0 dropout p_full=0.019 →
p_drop=0.037, enriched β regions Hippocampus + ctx-lh-fusiform.

---

## 8. What the writing agent should NOT do

- Do NOT recompute the analysis. Numbers are CSV-on-disk.
- Do NOT re-derive the framing decisions in §1. Those are locked.
- Do NOT add new figures. The 4 headline figures are the manuscript set.
- Do NOT change the joint MTC m=48 number. m=42 is controls only;
  m=60 includes cross-validation but loses the joint claim.
- Do NOT re-introduce "left temporal lobe" as the anatomical headline.
  Specific regions: Hippocampus + left fusiform, sampling-corrected.
- Do NOT cite measures 09 / 13 as "load-bearing" — they are
  consistency / descriptive layers. The load-bearing pair is
  measure 11 (β: 4 facets) + measure 12 (low_γ per-leaf).

---

## 9. Files the writing agent will produce

The expected deliverables:
1. Section 5 LaTeX/Markdown source — six paragraphs (head, methods, β
   results, low_γ + cross-validation, anatomy, limitations) ≈ 1500 words.
2. Figure caption text for 5A-5D using the stubs in the manuscript draft.
3. Supplementary table list pointing at the source CSVs.
4. A reproducibility appendix snippet (use §7 above).

Estimated effort: half a day of writing + revision, given that all
analysis, figures, and prose drafts already exist.
