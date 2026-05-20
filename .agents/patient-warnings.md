---
name: patient-warnings
type: warning_index
era: IMCOH_ABS × COHORT_N10
status: current
created: 2026-05-11
updated: 2026-05-11
pointers:
  - .agents/guides/02_methods/probe-bias-guide.md
  - .agents/guides/03_implementation/data-layout.md
---

# Patient-specific warnings and known biases

Per-patient pitfalls discovered during analysis that can mislead
cohort-level inference if not flagged. Append new findings here as
they surface; cite from any report/figure that touches the affected
patient × band × scale.

---

## Pat_05  ·  probe-shaft bias in exact-match LRG clades

**Surface:** β-band LRG dendrograms under `imcoh_abs`. Exact-equality
clade matching between tTest and rPost (script
`fig_preprint_dendrograms_matched_trace.py`, threshold = exact set
equality, MIN_CLADE_SIZE = 3, MAX_CLADE_SIZE = n // 2) returns several
large "preserved clades" that turn out to be **pure shaft bundles**:

- **18-leaf clade** = all leaves of probes `X'` (8) + `Y'` (10), no
  other probes touched.
- **10-leaf clade** = probes `G` (3) + `X` (3) + `Y` (4), three
  adjacent shafts only.

These bundles appear identically in tTest and rPost (and approximately
in rPre + tLearn) because adjacent contacts on the same sEEG shaft
retain residual coherence even under `imcoh_abs`. They are
anatomy-driven, not task-state-driven. CLAUDE.md flags this as the
"same-probe MSC bias" — `imcoh_abs` mitigates the zero-phase-lag
component but does not erase shaft-neighborhood clustering.

**Symptom:** When evaluating "preserved test↔post clade" measures
(matched_trace, KC, subtree Jaccard), Pat_05 produces big colored
blocks that are stable across all four phases and dominate the
visual signal.

**Reading rule:** Any clade with `len(set(probes_in_clade)) < 3`
should be flagged as probe-driven, not functional. The
matched_trace figure does *not* currently filter on probe diversity
— readers should manually check probe composition before claiming
trace.

**Inspection one-liner** (probe composition for a clade with leaves L):
```python
from lrg_eegfc.utils.io.patient import load_channel_labels
from lrg_eegfc.utils.probe import extract_probe_labels
labels = load_channel_labels("Pat_05")
probes = extract_probe_labels(labels)
print({p: sum(1 for L in clade_leaves if probes[L] == p)
       for p in sorted({probes[L] for L in clade_leaves})})
```

**Mitigation options** (not yet applied to the figure):
1. Filter — drop preserved clades with `n_unique_probes < 3`.
2. Annotate — show probe-count badge per preserved clade so reader
   can disregard probe-narrow matches.
3. Both — filter + annotate.

**Cross-band reach:** Probe-shaft bias is anatomy-driven, so it
generalises across bands. The X'+Y' bundle is likely preserved in
α / θ / γ as well; needs a quick check before generalising any
matched-trace finding from Pat_05 to the cohort.

**Source observation:** chat 2026-05-11, looking at
`data/reports/preprint/figure1_beta_trace/dendrograms_matched_trace/matched_trace_beta_all_patients.pdf`.
