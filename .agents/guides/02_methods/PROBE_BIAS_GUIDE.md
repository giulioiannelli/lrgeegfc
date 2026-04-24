---
name: probe-bias-guide
type: guide
era: IMCOH_ABS
status: current
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

# Same-Probe MSC Bias — Critical Guide

## The Problem

All sEEG data in this project uses **common-reference montage** (referenced
to contact G2). Contacts on the same electrode probe (e.g., A1-A12) share:

1. **Volume conduction** — physically adjacent contacts (~3-5mm spacing)
   record overlapping neural populations
2. **Common reference noise** — V_i - V_G2 injects G2 activity into every
   channel
3. **Shared far-field signals** — distant sources appear identically at
   nearby contacts

This produces **artificially high MSC between same-probe contacts** (2-8×
higher than cross-probe MSC), which dominates the network structure.

## Impact on LRG Analysis

The LRG hierarchy is built from the MSC-derived Laplacian. The inflated
same-probe MSC causes:

- **Community structure at coarse scales follows probe geometry.** At n=10
  communities, ~40% of same-community pairs are same-probe contacts (vs ~10%
  expected by chance). At n=30 this rises to 40-80%.
- **The dendrogram merges same-probe contacts first** because they have the
  shortest ultrametric distance.
- **Any metric derived from community assignment** (participation coefficient,
  community cosine, metastability, allegiance) is confounded at scales n<~15-20.
- **Intra-group coherence** for any spatially contiguous group (e.g., epileptic
  nodes on the same probe) is trivially inflated.

### Quantified same-probe enrichment in LRG communities

```
n_communities | Pat_02 | Pat_03 | Pat_05 | Pat_07 | Pat_08
-----------+--------+--------+--------+--------+--------
     3     |  0.9×  |  1.0×  |  1.0×  |  1.1×  |  1.0×
     5     |  1.6×  |  1.0×  |  1.1×  |  1.8×  |  1.3×
    10     |  3.9×  |  1.4×  |  1.4×  |  2.0×  |  1.5×
    15     |  5.6×  |  1.7×  |  3.7×  |  2.7×  |  2.1×
    20     |  6.2×  |  4.0×  |  3.7×  |  2.8×  |  2.4×
    30     |  8.2×  |  5.8×  |  6.5×  |  4.1×  |  4.0×
```

(Enrichment = fraction of same-community pairs that are same-probe,
divided by the expected fraction under random assignment.)

## How to Detect the Bias

**DO NOT use Variation of Information (VI)** to compare community labels
with probe labels. VI can appear "high" (suggesting communities ≠ probes)
even when communities clearly follow probe spatial geometry. VI misses
the many-to-one mapping where multiple probes map to one community.

**Instead use:**

```python
# Same-community/same-probe enrichment
same_comm = labels[i] == labels[j]  # for all i,j
same_probe = probe[i] == probe[j]
enrichment = (same_comm & same_probe).sum() / same_comm.sum()
expected = same_probe.sum() / (N * (N-1))
ratio = enrichment / expected  # >1 means probe-driven
```

**Or visually inspect** the 3D brain connectome HTML plots at
`data/figures/brain_connectome_multiscale/` — probe-driven communities
are immediately visible as coloured clusters along electrode trajectories.

## How to Debias

### THE CORRECT FIX: Use ImCoh instead of MSC

The proper solution is to use a functional connectivity metric that is
**immune to volume conduction by physics**, rather than trying to patch
a contaminated metric after the fact. The imaginary part of coherency
(ImCoh, Nolte 2004) achieves this: volume conduction is instantaneous,
so it produces zero phase lag, which means a purely real cross-spectral
density. The imaginary part of CSD is therefore exactly zero for any
volume-conducted signal, and ImCoh inherently ignores it.

**Results:** With ImCoh, the same-probe mean/cross-probe mean ratio drops
from 5.5x (MSC, alpha, Pat_02) to 1.35x. The bias is eliminated at
the source rather than patched post-hoc.

**Full same-probe ratio comparison:**

| Metric | Range across patients/bands | Typical |
|--------|---------------------------|---------|
| MSC | 1.6-6.0x | 3-5x |
| ImCoh | 0.6-2.0x | ~1x |

The ImCoh ratio of 0.6-2.0x is dramatically reduced from MSC's 1.6-6.0x.
The remaining ratio (especially ~2x in beta for Pat_07) reflects genuine
lagged neural coupling between nearby contacts, not artifact. This is
expected: physically adjacent sEEG contacts DO share genuine local circuit
interactions, and ImCoh correctly preserves these while removing the
volume-conduction component.

**LRG confirms ImCoh is not a rescaling:** The ImCoh-derived LRG
dendrograms are fundamentally different from MSC-derived ones (120/120
cases have different merge pairs and heights; community overlap at n=5 is
4-15%). This proves ImCoh captures a genuinely different network structure,
not just a rescaled version of MSC.

See `.agents/guides/IMCOH_GUIDE.md` for full details on the ImCoh
implementation, LRG results, and the caching bug.

---

### DEPRECATED APPROACHES (and why they fail)

The following approaches were investigated and found to be problematic.
They are documented here for completeness but **should not be used**.

#### Approach 1: Zero Same-Probe Edges (ad hoc, loses information)

```python
def debias_same_probe(A, channel_labels):
    """Zero out edges between contacts on the same probe."""
    import re
    probe_labels = [re.match(r"([A-Za-z]+'?)", l).group(1) for l in channel_labels]
    A_out = A.copy()
    for i in range(len(A)):
        for j in range(len(A)):
            if probe_labels[i] == probe_labels[j]:
                A_out[i, j] = 0.0
    return A_out
```

- Removes ~10% of edges (same-probe pairs)
- Graph remains connected (90% of edges are cross-probe)
- **Problem: Destroys real within-probe information** (local circuit
  dynamics between nearby contacts are genuine and scientifically
  interesting)
- **Problem: Ad hoc** -- there is no principled reason to zero these
  edges entirely rather than correct them
- **Problem: Would not survive peer review** -- reviewers would
  (correctly) ask why real connectivity was discarded

#### Approach 2: Percentile Rescaling (ad hoc, non-physical)

Rescale same-probe edges to match the distribution of cross-probe edges
by percentile mapping.

- **Problem: Entirely ad hoc** -- no physical or statistical justification
- **Problem: Non-physical** -- the rescaled values do not correspond to
  any meaningful coherence quantity
- **Problem: Would not survive peer review** -- arbitrary data
  manipulation with no theoretical backing

#### Approach 3: Bipolar Re-Referencing (MAKES THE BIAS WORSE)

Subtract consecutive contacts on the same probe BEFORE computing MSC:
`signal_bipolar_k = V_{i+1} - V_i`. This cancels the common reference
and shared far-field at the signal level.

- Standard in clinical sEEG
- Reduces channel count by ~10%
- Addresses the common-reference component
- **CRITICAL PROBLEM: Makes the same-probe bias WORSE, not better.**
  The same-probe ratio goes from ~5x (common reference) to **8-17x**
  (bipolar). This happens because adjacent bipolar channels share a
  physical contact: bipolar channel k = V_{k+1} - V_k and channel
  k+1 = V_{k+2} - V_{k+1} both contain V_{k+1} with opposite sign.
  This shared contact creates a new, stronger source of spurious
  coherence that was not present in the common-reference montage.
- **Do not use bipolar re-referencing to address same-probe bias.**

## Rules for Agents

1. **NEVER report community-level results without checking probe enrichment**
   at the relevant scale. If enrichment > 2×, the result is suspect.

2. **NEVER claim "epileptic nodes form a community"** or "cluster together
   in the dendrogram" without verifying on the debiased matrix. Same-probe
   contacts always cluster together — this is trivial.

3. **NEVER use intra-group MSC** (for any spatially contiguous group)
   without the probe-matched null: compare against random same-size
   groups drawn from the same probes with the same consecutive structure.

4. **When comparing epileptic vs non-epileptic contacts:**
   - Use within-probe comparisons (contacts on the same probe)
   - AND verify on the debiased matrix (same-probe edges zeroed)
   - Both must agree for the result to be trustworthy

5. **For LRG community analysis at any scale:**
   - Compute and report the same-community/same-probe enrichment
   - If enrichment > 2×, discount the finding or use the debiased matrix
   - The debiased LRG hierarchy (from same-probe-zeroed MSC) is the
     honest version for community-level claims

6. **The metastability index** from the LRG coarsening is ALSO affected.
   Metastable nodes may be "metastable" because they sit at probe
   boundaries (where the trivially-cohesive probe group meets another).
   Re-evaluate metastability on the debiased matrix before interpreting.

7. **Prefer ImCoh over MSC** for any analysis where same-probe bias
   matters. Use `fc_method="imcoh_abs"` in `workflow.fc.load_fc_matrix()`
   (the `|ImCoh|` convention — Ewald 2012 / Bastos & Schoffelen 2016).
   See `.agents/guides/02_methods/IMCOH_GUIDE.md` for the full three-method
   taxonomy (`imcoh` signed, `imcoh_abs`, `imcoh_sq`) introduced in the
   2026-04-15 reset.

## Probe Structure Reference

Each patient has 11-14 probes with 4-18 contacts per probe.
Same-probe pairs are ~9-10% of all pairs.

```
Pat_02: 117 channels, 11 probes (A:9, B:10, F:11, G:12, H:15, L:12, M:8, N:8, P:14, Q:14, T:4)
Pat_03: 122 channels, 12 probes (A:14, B:10, C:10, D:11, E:12, F:13, H:5, I:11, L:12, O:8, P:10, Y:7)
Pat_05: 118 channels, 14 probes (F:8, G:13, G':10, H:7, H':7, L:8, L':8, O:8, P':9, X:7, X':11, Y:10, Y':12)
Pat_07: 116 channels, 14 probes (A:10, F:15, I:8, J:10, K:10, L:8, M:10, P:18, U:8, X:8, Y:10)
Pat_08: 120 channels, 12 probes (A:11, B:11, C:12, D:12, F:12, P:15, R:10, S:9, T:8, U:10, W:11)
```

## Related Files

- Memory: `~/.claude/projects/.../memory/probe_bias_critical.md`
- Brain connectome plots: `data/figures/brain_connectome_multiscale/`
- Plot generator: `scripts/gen_brain_connectome_multiscale.py`
- Epileptic node loader: `src/lrg_eegfc/utils/io/patient.py:load_epileptic_nodes()`
- Epileptic investigation report: `data/figures/epileptic_analysis_v3/REPORT_FINAL.md`
- ImCoh guide: `.agents/guides/02_methods/IMCOH_GUIDE.md`
- ImCoh implementation: `src/lrg_eegfc/utils/fc/coherence/imcoh.py` (signed Nolte 2004)
- Loader + magnitude transforms: `src/lrg_eegfc/workflow/fc.py` (`fc_method` ∈ `{"imcoh", "imcoh_abs", "imcoh_sq"}`)
