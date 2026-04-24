---
name: imcoh-recovery-plan
type: plan
era: IMCOH_ABS
status: current
created: 2026-04-24
updated: 2026-04-24
pointers: []
---

# ImCoh Recovery Plan — draft for user review

**Status:** draft — do not execute until user approves. Author: Claude Code session 2026-04-15.

---

## 1. What went wrong (plain)

- The project's computation code (`src/lrg_eegfc/utils/fc/msc/msc.py` `metric="imcoh"`) computes `Im(S_ij)² / (S_ii · S_jj)` — the **squared** imaginary coherence, range `[0, 1]`, non-negative.
- Across this session I (Claude) repeatedly referred to that quantity as `|ImCoh|`. That is wrong. `|ImCoh|` = `|Im(S)|/√(S_ii·S_jj)`, and `|ImCoh|² = Im(S)²/(S_ii·S_jj)` = what is cached. They differ by `sqrt`.
- Every user-facing number I quoted that said `|ImCoh| ≈ X` should have read `|ImCoh|² ≈ X`.
- Qualitative conclusions (rankings, community structure, MSC-vs-ImCoh trends, unanimity counts, enrichment ratios) are preserved under `sqrt` (monotonic transform on non-negatives). Quantitative absolute values differ.
- The canonical literature quantity (Nolte et al. 2004) is the **signed** `ImCoh = Im(S_ij)/√(S_ii·S_jj)`, range `[-1, 1]`. We have never computed that.

## 2. What is actually salvageable from 3 days of work

### Safe — stays correct under the sqrt rewrite
- **All LRG dendrograms** (linkage trees, ultrametric matrices, community assignments at every scale k). Identical up to monotonic transform — LRG uses ranks internally via the Laplacian.
- **Probe-bias enrichment ratios** (same-community / same-probe enrichment). Identical.
- **Unanimity counts** (H1, H2a, H2b, H3) per (band, k). Identical — driven by VI distances between partitions, which only depend on community structure.
- **Variation of Information profiles**. Identical.
- **Figure `fig_D1_enrichment_heatmap_redesigned`, `fig_D2`, `fig_C3`, `fig_E2` (networks from LRG ultrametric)**: visually identical.
- **Figure `fig_G1` (Laplacian eigenvalue spectrum)**: identical (Laplacian of squared vs magnitude differ only by sqrt of same vector — rank-preserving).

### Quantitatively changes — must be re-quoted
- **Absolute mean / median / percentiles of the FC matrices** (numbers like "same-probe ratio ≈ 5×").
- **Same-probe / cross-probe mean ratios**. Under `|·|`, ratios will be larger than under `|·|²` (because sqrt(a)/sqrt(b) > a/b when b < a < 1).
- **Spectral distribution plots (`fig_C2_spectral_distribution_*`)**. Numbers on the y-axis change.
- **Adjacency heatmap absolute color scales (`fig_A`, `fig_B`)**. Visual shape is the same, but vmax values change.
- **Every scalar "X = 0.XX" number I wrote into descriptions, handoff docs, or memory files**.

### Forever lost (only if we delete now without snapshotting)
- Nothing actually lost — signed `imcoh` can always be computed from scratch, and `|imcoh|` / `|imcoh|²` are trivial transforms of signed.

## 3. Proposed literature-backed taxonomy

| Name | Formula | Range | Source | When to use |
|---|---|---|---|---|
| `imcoh` | `Im(S_ij)/√(S_ii·S_jj)` | [-1, 1] | Nolte et al. 2004 | **canonical**. Any future directional / lag-direction analysis. |
| `absimcoh` | `|Im(S_ij)|/√(S_ii·S_jj)` | [0, 1] | Standard connectivity-strength convention (Ewald 2012; Bastos & Schoffelen 2016) | Undirected community analysis, LRG input, adjacency heatmaps, same-probe analysis. |
| ~~`imcoh_sq`~~ | `Im(S)²/(S_ii·S_jj)` | [0, 1] | Not a standard literature quantity | **DROP.** No scientific justification to cache or report this. |

**Caching policy proposal:** Cache only signed `imcoh` on disk. Compute `absimcoh = np.abs(signed)` on the fly in the loader. Saves disk, guarantees consistency, enables signed analysis with no extra cost. (Alternative: cache both explicitly. Costs 2× disk but avoids load-time `abs`.)

## 4. Open decisions for you

1. **Primary science quantity.** Use `absimcoh` for all current Section 2 analysis (matches literature "ImCoh" convention for connectivity strength)?  Yes / No.
2. **Signed `imcoh` computation priority.** Compute signed now (so both `imcoh` and `absimcoh` become available) or defer signed to a separate future work item?
3. **Disk caching policy.** Store only signed; derive abs on the fly. Yes / No. (If no: cache both.)
4. **Snapshot before wipe.** Before deleting `data/cache/imcoh*` (or renamed `absimcoh_sq*`), create a git tag `pre_imcoh_recovery` and / or copy to `data/cache/_archive_imcoh_sq/` so nothing is irreversible. Yes / No.
5. **Verification strategy.** After recompute, how do you want to verify correctness? Options: (a) compare one patient-band against an independent reference (e.g., MNE-Python's imaginary coherence), (b) analytic self-test (symmetry, diagonal zero, range), (c) visual sanity check only.
6. **Section 2 handoff state.** For `for_writing_agent/`: (a) leave as-is until recompute verified, then overwrite; (b) pull the figures now and replace with a note "ImCoh figures pending recomputation"; (c) rebuild now even if numbers shift a bit.
7. **Memory / docs / reports cleanup.** Every mention of numerical ImCoh values in `.agents/reports/` and memory files is potentially wrong (absolute values) or correct (rankings). Options: (a) mass-flag everything with a "QUANTITATIVE_INCORRECT" banner until updated, (b) update as we go, (c) leave and fix per-section.
8. **Failure safety.** Do I get to recompute without asking you again per step, or do you want to review each step (cache wipe → recompute signed → recompute LRG → regenerate figures)?

## 5. Proposed execution order (if/after plan is approved)

Step A — **Snapshot & freeze** (10 min)
- `git tag pre_imcoh_recovery HEAD`
- Create `data/cache/_archive_imcoh_sq/` if user chooses (4.4 above).
- Write a banner README in `for_writing_agent/` saying "figures reflect pre-recovery squared quantity; do not quote numbers pending recompute".

Step B — **Code layer fix** (30 min, low risk)
- `src/lrg_eegfc/utils/fc/msc/msc.py`: implement `metric="imcoh"` = signed (Nolte 2004); `metric="absimcoh"` = magnitude. Drop the squared branch.
- Update config, CLI, routing per the partially-completed rename agent.
- Unit test: small synthetic CSD, verify `imcoh ∈ [-1, 1]`, `absimcoh = |imcoh|`, `imcoh_symmetric_for_Im(S)_antisymmetric`.

Step C — **Compute signed `imcoh` for the cohort** (~70 min)
- 5 patients × 4 phases × 6 bands. Store signed only.
- Sanity: diagonal = 0, antisymmetric for ImCoh (Im part is antisymmetric on complex CSD).

Step D — **Recompute LRG from `absimcoh`** (~30 min)
- Wipe `lrg` caches, recompute from `|imcoh|`.
- Cross-check: community structure at k=5,10,20 should match the pre-recovery versions (we predict it will, since sqrt is monotonic).

Step E — **Regenerate figures & verify** (~30 min)
- Re-run every Section 2 figure script. Compare community figures to pre-recovery visually (expect identical structure). Compare absolute-number figures (adjacency, weight dist) and note any that look qualitatively different (these indicate a real discovery we missed in sq-version).
- Fix adjacency heatmap normalization (`PowerNorm(gamma=0.5)` or log-norm) to handle `|imcoh|`'s heavy tail.

Step F — **Update docs, memory, handoff** (~30 min)
- Mass-edit `.agents/` reports, memory files, CLAUDE.md to reflect the `imcoh` (signed) vs `absimcoh` (mag) split.
- Sync fresh figures into `for_writing_agent/`.

Total clock time: ~3.5 h, of which ~100 min is compute that runs unattended.

## 6. What I explicitly will NOT do without your go-ahead

- Touch `data/cache/`.
- Run any recompute.
- Delete any figure.
- Edit memory / docs / reports further.
- Spawn any more subagents that do destructive work.

## 7. Rollback strategy

If at any point the recomputed figures disagree qualitatively with the pre-recovery ones (indicating a bug in the new computation, not a discovery), roll back:
- `git checkout pre_imcoh_recovery -- src/lrg_eegfc/`
- Restore `data/cache/_archive_imcoh_sq/` → `data/cache/absimcoh_sq/`.
- Return to the (wrong-labeled but functional) prior state while we debug.

---

**Waiting for your call on items (1)–(8) in section 4 before any code moves.**
