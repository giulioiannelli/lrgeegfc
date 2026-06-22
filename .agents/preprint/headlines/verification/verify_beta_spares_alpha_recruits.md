---
name: verify-beta-spares-alpha-recruits
era: IMCOH_ABS_COHORT_N10
status: ready_to_run
kind: verification-brief
headline: N1.6 (bridge, .agents/preprint/headlines/01_trace.md)
owner_agent: white-matter / localization + null-model
updated: 2026-06-22
---

# Verify — β spares, α recruits the epileptic core (N1.6 bridge)

**Head.** Confirm the cross-band tissue dissociation that bridges the cognitive
trace (N1/N2) to the epileptogenic read-out (N3): the **β** cophenetic trace is
carried by healthy gray↔gray cortical pairs and the seizure core is **spared**,
while the **α** trace **recruits** the epileptic core (epi↔epi). The key thing to
verify is that this is a **genuine cross-band difference**, not a node-count
artifact or "α is just stronger everywhere."

## 5-point critical preamble

1. **Claim.** β epi↔epi pairs do **not** carry the trace (spared); α epi↔epi pairs
   **do** (recruited, matched-strength-significant). The **β-vs-α difference** on
   epi↔epi tissue is real.
2. **Null.** Per band, epi↔epi pair-class trace ≤ (a) matched-strength surrogate
   AND (b) a **random same-size pair subset** (pair-count null). For the
   interaction: β-vs-α epi↔epi difference ≤ a label/band-permutation null.
3. **Strongest alternative.** (A) **Node-count / subgraph confound** — rebuilding
   LRG on any node subset inflates ρ^coph. (B) **α-is-stronger** — α recruits the
   core only because α's trace is larger everywhere, not specifically in epi
   tissue. (C) **Epi-labeling** — SOZ labels are spatially clustered.
4. **Does the null control it.** vs (A): use **same-graph pair classes** (no
   rebuild) + the `audit_85 --mode pairclass` pair-count null (is epi↔epi more
   trace-bearing than a random equal-size pair set?). vs (B): test the
   **band × tissue interaction** (epi↔epi vs nonepi↔nonepi *difference*, β vs α) —
   "α stronger everywhere" predicts no interaction. vs (C): pair-class is
   tissue-membership, not proximity; |ImCoh| substrate removes zero-lag anyway.
5. **Falsification + limits.** Falsified if α epi↔epi fails its pair-count null, or
   the β-vs-α interaction is n.s., or β epi↔epi is actually carried. Limits:
   `epi_only` is underpowered (few patients with large SOZ); the dissociation is a
   *characterization* layer — neither band's verdict tag changes.

## Steps

1. Load the n=10 pair-class cophenetic table (epi↔epi / cross / nonepi↔nonepi per
   band) and the pair-count decimation table.
2. Confirm, per band: β epi↔epi matched-strength **n.s./depleted** (spared);
   α epi↔epi **concentrated + significant** (recruited).
3. Run the **pair-count null** (`audit_85 --mode pairclass`, R=1000): epi↔epi vs
   random equal-size pair subset, per band.
4. Test the **β × tissue vs α × tissue interaction** (paired across patients).
5. LOO-max per claim (no single-patient-p-driven).

## Data / scripts (no hardcoded numbers — read the CSVs)

- `data/audit/epi_stratified/cophenetic_cohort.csv` (n=10) ·
  `audit_77_epi_stratified_cophenetic.py` · 2026-06-16.
- `data/audit/epi_stratified/pairclass_decimation_cohort.csv` ·
  `audit_85_wm_decimation_control.py --mode pairclass` · 2026-06-12.
- Locked context: `../../locked/VERDICT_LEDGER.md` (α C5 epi↔epi), `../../bands/01_beta.md` §3.2.5.

## Pass / fail

- **Pass:** α epi↔epi beats its pair-count null AND the β-vs-α epi↔epi interaction
  is significant (LOO-robust). → N1.6 stands as a bridge.
- **Fail/soften:** interaction n.s. → report as a single-band α observation, drop
  the "dissociation" framing.

## Caveats

- `epi_only` underpowered; lead with the same-graph pair-class (clean) result.
- State both bands on the cophenetic object (raw |ImCoh| = baseline only).

## What to return

A one-paragraph verdict (pass/soften/fail), the CSV rows used, the interaction
test + LOO-max, and a figure (X-epi pair-class contrast; no C4). Fill below.

> **Result (agent fills):** …
