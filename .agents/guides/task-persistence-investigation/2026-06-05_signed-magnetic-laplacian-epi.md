---
name: signed-magnetic-laplacian-epi
type: scope
era: COHORT_N10
status: stage1-negative-depth-confound
created: 2026-06-05
updated: 2026-06-05
pointers:
  - .agents/guides/task-persistence-investigation/2026-06-05_epi-node-trace-marker.md
  - scripts/01_compute/audit/audit_88_signed_magnetic_laplacian.py
  - .agents/guides/02_methods/imcoh-guide.md
---

# Signed/directional ImCoh + magnetic-Laplacian eigenmodes as an epi-node marker

**Head.** Everything so far used `imcoh_abs = ⟨|ImCoh|⟩` — it throws away the
**sign** of the imaginary coherence, which encodes the **lead/lag direction** of
coupling. The signed ImCoh matrix `A` is **antisymmetric** (`A = −Aᵀ`), so it is
a *directional/flow* object, not a Kunegis symmetric-signed graph. The principled
operator is therefore the **Hermitian (magnetic) Laplacian** `L_H = D̄ − iA`
with `D̄_ii = Σ_j |A_ij|` — Hermitian (real eigenvalues, complex eigenvectors),
positive semi-definite, whose eigenvector **phases** encode each node's position
in the lead/lag flow. The question: do epileptic contacts (putative drivers)
occupy a distinct region of this directional structure — by net flow, by flow
coherence, or by localising in specific magnetic modes — beyond what node
strength already says?

> **Note on the user's `L_s = |D| − A`.** For a *symmetric* signed graph that is
> the Kunegis signed Laplacian. ImCoh is *antisymmetric*, so `D̄ − A` is
> non-Hermitian (complex spectrum); `eigvalsh` would silently symmetrise it and
> give a meaningless answer. `L_H = D̄ − iA` is the correct Hermitian realisation
> of the same idea for a directional graph (`iA` is Hermitian when `A` is real
> antisymmetric). We implement `L_H`.

## 5-point critical preamble

1. **Claim.** Epileptic contacts have a distinctive *directional* signature in
   the signed-ImCoh magnetic Laplacian — distinctive beyond node strength
   (|ImCoh| degree) and beyond what `imcoh_abs` already captured.
2. **Null (staged).** Stage 1 (this pass): within-patient epi-vs-non-epi AUC vs a
   node-strength baseline + Spearman(feature, strength) confound. Stage 2 (gated
   on signal): a **sign/direction-randomising matched-strength surrogate**
   (preserve `|A|` row sums, randomise edge signs / antisymmetric orientation) —
   the mandatory cohort control; the existing non-negative matched-strength
   surrogate does NOT apply to a directional graph, so this is a NEW null to
   build only if Stage 1 shows signal.
3. **Strongest alternative.** It is (a) just `|ImCoh|` strength again, or (b)
   sEEG-shaft geometry (nearby contacts share lead/lag), or (c) volume-conduction
   residue — though ImCoh is built to suppress zero-lag (Nolte 2004), the *sign*
   is exactly the lag structure VC does not produce, which is the upside.
4. **Mechanical reach.** Stage-1 strength-correlation + the within-shaft control
   (reuse audit_87's design) address (a)/(b); the sign-randomising surrogate
   (Stage 2) is the decisive directional null.
5. **Falsification.** If no directional feature separates epi within-patient
   beyond strength and the within-shaft control, the sign carries no epi marker
   and `imcoh_abs` lost nothing — report that.

## Definitions (per node `i`, phase rest_post, band `b`)

`A` = signed band-averaged ImCoh (antisymmetric, `load_fc_matrix(...,"imcoh")`).

| name | formula | reads |
|---|---|---|
| `net_flow` | `Σ_j A_ij` | signed net lead−lag; a consistent driver/sink has large `|net_flow|` |
| `flow_imbalance` | `|Σ_j A_ij| / Σ_j |A_ij|` ∈ [0,1] | how *directional* the node is (1 = pure source/sink, 0 = balanced relay) |
| `mag_lowmode_loc` | `Σ_{k<K} |U_k(i)|²`, `U_k` = lowest-`λ` magnetic eigenvectors | localisation in the low (slow-flow) magnetic modes — flow bottleneck / frustration |
| `mag_phase_disp` | circular spread of `arg U_k(i)` over `k<K` | phase (in)consistency of the node across slow flow modes |
| (ref) `strength_abs` | `Σ_j |A_ij|` = `D̄_ii` | the `imcoh_abs` degree — the hubness baseline |

No new statistic is pooled into a scalar; the features stay a vector.

## Implementation

- `audit_88_signed_magnetic_laplacian.py` (Stage 1): load signed `A` (rest_post,
  6 bands, n=10), `L_H = D̄ − iA`, `eigh` (Hermitian), per-node features above,
  within-patient epi-vs-non-epi Mann-Whitney AUC + cohort Wilcoxon + strength
  Spearman + the audit_87 within-shaft recovery (does the directional signature
  find hidden epi within shaft). Output `data/audit/epi_signed_laplacian/`.
- Library: if a magnetic-Laplacian helper gains a 2nd caller, promote to
  `src/lrg_eegfc/utils/lrg/` with a general name (`magnetic_laplacian`,
  `hermitian_laplacian_eigs`) — never an epi/manuscript token.
- Stage 2 (deferred, gated): sign-randomising matched-strength surrogate +
  per-feature surrogate-z; only if Stage 1 separates.

## Open questions
- `K` for the low-mode features (start 5; sweep if signal).
- rest_post vs rest_pre (intrinsic) vs task — start rest_post (trace destination).
- Is `net_flow` sign meaningful cohort-wide, or only `|net_flow|`/imbalance
  (gauge/orientation of ImCoh sign is a fixed convention, so signed net_flow is
  comparable across nodes within a patient, but cross-patient sign needs care).

---

## RESULT (2026-06-05 — `audit_88` Stage 1; NEGATIVE, depth confound)

The Hermitian magnetic Laplacian `L_H = D̄ − iA` was built and is correct
(Hermitian, PSD). Directional features (`net_flow`, `flow_imbalance`,
`mag_lowmode_loc`) recover hidden epi within-shaft above chance in all bands
(p≤0.006; δ 0.717 striking) — BUT this is the **along-shaft DEPTH confound**, not
the ImCoh sign. Plain contact depth (`contact_index_norm`) recovers hidden epi
within-shaft at 0.72–0.74 every band, **beats** the directional subspace, and the
directional features **add nothing to depth** (paired Wilcoxon p ≥ 0.33). So the
ImCoh sign / magnetic-Laplacian carries **no epi marker beyond electrode depth**.
**Stage 2 (sign-randomising surrogate) NOT warranted.** Full audit trail +
depth-baseline table in `.agents/reports/2026-06-05_propagator-subspace-epi-recovery.md`
§10. The magnetic-Laplacian construction itself is reusable (correct operator for
directional ImCoh) and is kept for any future directional-FC question.
