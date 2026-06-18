---
name: epi-distant-marker-walkthrough
era: IMCOH_ABS_COHORT_N10
date: 2026-06-16
status: current
kind: report
scope: Step-by-step technical companion to fig_epi_distant_soz_marker.pdf (preprint_39) — every formula, score, and test behind the verified δ heat-kernel distant-SOZ marker, in plain language with math + code.
related: 2026-06-12_propagator-distant-soz-marker.md
---

# How the distant-SOZ marker works — formulas, scores, tests, step by step

**Head.** This walks through, one step at a time, exactly how we turn the connectivity
matrix into a number that points at distant epileptic contacts, and how we proved that
number is real. It pairs with the figure `fig_epi_distant_soz_marker.pdf`: each section
ends with **→ panel (x)** telling you which part of the picture it explains. Plain
language first in every section; the formula and the actual code sit just below it.

A quick glossary so nothing is opaque:
- **Contact** = one recording point on a depth electrode; **shaft** = one electrode (a line
  of contacts pushed into the brain).
- **SOZ** = *seizure-onset zone* = the contacts a clinician marked as where seizures start.
- **AUC** = a ranking score from 0 to 1. **0.5 = a coin flip, 1.0 = perfect.** It answers:
  "pick one true SOZ contact and one healthy contact at random — does my measure score the
  SOZ one higher?" 0.72 means "yes, 72% of the time."
- **τ (tau)** = diffusion time: how long we let signal spread through the network. Small τ =
  a quick spread that stays local; large τ = a slow spread that reaches far.

---

## Step 0 — The input: one connectivity matrix per patient, band, phase

**Plain.** For each patient we build a network of their contacts. The wire between contact
$i$ and contact $j$ is how strongly they are connected, measured by **imaginary coherence**
(a connectivity measure that is immune to volume conduction — it ignores the trivial
"two contacts pick up the same source" artefact). We use the **resting-after-task** recording
(`rest_post`) in a given frequency band.

This gives a symmetric matrix $W$ (size $N\times N$, $N$ = number of contacts), with
$W_{ij}\ge 0$ and zero on the diagonal.

```python
W = load_phase_fc(pat, "rest_post", band)   # |ImCoh| adjacency, N x N, symmetric, >= 0
```

The clinician's labels give us, for each contact, two things we use: is it SOZ (yes/no),
and which shaft it sits on.

```python
pm = build_epi_masks(pat)          # canonical loader
epi  = np.asarray(pm.epi_mask, bool)     # True = SOZ contact
probes = np.asarray(pm.probes, object)   # which shaft each contact is on
```

---

## Step 1 — The Laplacian and the heat-kernel propagator (the LRG object)

**Plain.** From the network we form the **graph Laplacian** $L$. The Laplacian is the
operator that governs diffusion on a graph: if you drop a unit of "heat" on one contact and
let it spread along the wires, $L$ tells you how. The amount of heat that has flowed from
contact $j$ to contact $i$ after time $\tau$ is one entry of the **heat-kernel propagator**
$e^{-\tau L}$. This is precisely the LRG object $\rho(\tau)=e^{-\tau L}/Z$ — the global
constant $Z$ (a normalisation) cancels out the moment we *rank* contacts, so we use the bare
$e^{-\tau L}$.

$$
L = D - W,\qquad D=\mathrm{diag}\!\Big(\textstyle\sum_j W_{ij}\Big),\qquad
\rho(\tau)=e^{-\tau L}=U\,\mathrm{diag}\!\big(e^{-\tau\lambda_k}\big)\,U^{\top}
$$

where $L = U\,\mathrm{diag}(\lambda_k)\,U^{\top}$ is the eigendecomposition ($\lambda_k\ge 0$).

```python
d = W.sum(1)
L = np.diag(d) - W
lam, U = np.linalg.eigh(L)                      # eigenvalues lam, eigenvectors U
heat_tau = U @ (np.exp(-tau * lam)[:, None] * U.T)   # = e^{-tau L}
```

**The diffusion time τ.** We don't pick one τ; we scan six, log-spaced from fast to slow,
set by the network's own fastest timescale $1/\lambda_{\max}$:

$$
\tau \in \Big[\tfrac{1}{\lambda_{\max}},\ \tfrac{10}{\lambda_{\max}}\Big],\quad
\text{6 values, geometric spacing}
$$

```python
taus = np.geomspace(1/lam[-1], 10/lam[-1], 6)   # t0=fast ... t5=slow
```

`heat_t0` = fastest (heat barely leaves home), `heat_t5` = slowest (heat reaches across the
brain). **The headline marker is `heat_t5`, the slowest** — because the task is to reach
*distant* contacts.

---

## Step 2 — The marker: a contact's diffusion affinity to the known SOZ

**Plain.** Given a set $S$ of *known* SOZ contacts (the "seeds"), the score of any other
contact $i$ is simply **how much heat reaches $i$ from the seeds**: average the propagator
over the seed columns.

$$
\mathrm{marker}_S(i)\;=\;\frac{1}{|S|}\sum_{s\in S}\big[e^{-\tau L}\big]_{i,s}
$$

```python
def marker(S):                 # S = indices of seed (known) SOZ
    return heat_tau[:, S].mean(axis=1)
```

A high score means "this contact is strongly tied, through slow diffusion, to the known
seizure zone." The whole claim is that true SOZ contacts get higher scores than healthy
ones — *even when they sit far from the seeds*.

---

## Step 3 — Why we throw spatial proximity away (the one honest choice)

**Plain.** Here is the trap that sinks naïve markers. A SOZ contact is labelled SOZ because
it sits inside a clinically-marked epileptogenic *region*. So "a contact near a known SOZ is
also SOZ" is **built into the labels** — it is circular, not a discovery. A dumb rule that
just flags the nearest contacts scores a great AUC and discovers nothing.

We confirmed this directly: on the easy "rank all contacts" task, plain distance-to-seed
**beats** the propagator about 2× (recall@10 ≈ 0.46–0.68 for distance vs 0.27–0.35 for the
marker — `audit_99`). That gap is the tautology talking, not skill. **So we discard proximity
entirely** and only ever ask the question proximity *cannot* answer: distant discovery.

→ This is the reason for the whole leave-one-shaft-out design in Step 4. (It is documented in
the figure's companion text rather than a panel, because it is a "why", not a result.)

---

## Step 4 — The honest score: leave-one-shaft-out distant discovery

**Plain.** Hide every SOZ contact on **one whole electrode shaft**. Use the SOZ on the
*other* shafts as the seeds. Now score contacts — but only contacts that live on shafts
containing **no seed**. Among those far-away contacts, do the hidden shaft's SOZ outrank the
healthy ones? Repeat, hiding each SOZ shaft in turn. Proximity is useless here by
construction: everything being scored is far from every seed.

```text
for each shaft p that carries SOZ:
    seeds   = SOZ contacts NOT on p
    targets = SOZ contacts ON p            (hidden — the "can we find them?" set)
    off     = contacts whose shaft has no seed   (the only contacts we score)
    cases   = targets within off           controls = healthy contacts within off
    accumulate concordance( marker[cases] > marker[controls] )
```

**The score itself is an AUC computed as rank concordance** (the Mann–Whitney form): over
every (case, control) pair, count how often the case outranks the control.

$$
\mathrm{AUC}=\frac{\#\{\,\mathrm{marker}(c)>\mathrm{marker}(h)\,\}+\tfrac12\#\{\text{ties}\}}
{\#\text{cases}\times\#\text{controls}},\qquad c\in\text{cases},\ h\in\text{controls}
$$

```python
def concordance(case, ctrl):
    c = sum((case[:,None] > ctrl[None,:]).sum() + 0.5*(case[:,None]==ctrl[None,:]).sum())
    return c / (case.size * ctrl.size)
```

(Full loop: `loso_auc` in `audit_101`.)

---

## Step 5 — Removing the second cheat: node strength (hubness)

**Plain.** The other trivial explanation is "SOZ are just the most-connected hubs." We remove
it. Before scoring, we **linearly regress the marker on each contact's total connection
strength and keep only the residual** — the part of the score that strength cannot explain.

$$
\text{strength}_i=\sum_j W_{ij},\qquad
\widetilde{\mathrm{marker}}_i=\mathrm{marker}_i-\big(a\cdot\text{strength}_i+b\big)
$$

where $a,b$ are the least-squares line fit of marker vs strength.

```python
a, b = np.polyfit(strength, marker, 1)        # best-fit line
resid = marker - (a*strength + b)             # AUC is computed on this
```

So the headline number — **"strength-residual leave-one-shaft-out AUC"** — has *neither*
proximity *nor* hubness propping it up. As a sanity check we also run plain node strength as
its own marker; on the residual scale it sits exactly at chance (0.5), as it must.

→ **panel (c)**: the red `strength (baseline)` bar sits at 0.5.

---

## Step 6 — The result, patient by patient

**Plain.** With the slow heat kernel (`heat_t5`) in the **δ band**, the strength-residual
distant-discovery AUC is:

| patient | AUC | reading |
|---|---|---|
| Pat_08 | 0.95 | recovered |
| Pat_14 | 0.91 | recovered |
| Pat_05 | 0.84 | recovered |
| Pat_06 | 0.82 | recovered |
| Pat_13 | 0.76 | recovered |
| Pat_03 | 0.68 | recovered |
| Pat_02 | 0.66 | recovered |
| Pat_07 | 0.57 | recovered (marginal) |
| Pat_10 | 0.41 | **fails — hub-patient** |
| Pat_15 | 0.21 | **fails — hub-patient** |

**Median 0.72, works in 8 of 10.** The two failures are the patients whose SOZ *are* the
network hubs (a different regime — there, plain strength predicts them and diffusion does
not). This is the same two-population split we see throughout the epilepsy investigation.

→ **panel (a)**: green = the 8 above the null, red = the 2 hub-patients; the green line is the
median 0.72.

---

## Step 7 — Why δ and why *slow*: the mechanism

**Plain.** Two pre-specified predictions, both confirmed in the data:

1. **Slow diffusion is required.** Fast τ keeps heat local, so it cannot reach a hidden
   shaft; only slow τ does. The cohort-median AUC climbs from ≈0.45 at the fastest τ to
   ≈0.75 at slow τ. This is not a knob we tuned — "slow reaches far" is exactly what you would
   write down *before* looking.
2. **It is δ-specific.** The slow δ rhythm carries the distant signal; α, β, low-γ do not
   (their AUCs hover at the null's upper edge, ≈0.59–0.61, not significant).

→ **panel (b)** shows the rise with τ; **panel (d)** shows δ separated from the other bands.

---

## Step 8 — The propagator wins a fair contest (16 operators)

**Plain.** We did not assume the heat kernel is special — we made it earn it. We built **16
different network operators** from the same Laplacian and scored them all on the identical
leave-one-shaft-out test. The heat kernel (our propagator) and its normalized-Laplacian
cousin came out on top; everything else trailed.

| operator | what it is | δ AUC (best) |
|---|---|---|
| **heat kernel** | $e^{-\tau L}$ — the LRG propagator | **0.75** |
| norm.-Lap. heat | $e^{-\tau L_{\mathrm{norm}}}$, $L_{\mathrm{norm}}=I-D^{-1/2}WD^{-1/2}$ | 0.74 |
| PageRank | $(1-\alpha)(I-\alpha D^{-1}W)^{-1}$ | 0.62 |
| communicability | $\exp(W/\rho(W))$ | 0.61 |
| diffusion-distance | distance built from heat columns | 0.59 |
| neg. resistance | $-R$, $R_{ij}=L^{+}_{ii}+L^{+}_{jj}-2L^{+}_{ij}$ | 0.56 |
| strength baseline | total connection strength | 0.50 |
| Katz | $(I-\beta W)^{-1}$ | 0.48 |

→ **panel (c)**.

---

## Step 9 — Proof it is real: two independent tests

**Plain.** A good-looking AUC could still be (i) a fluke that any contact-set would produce,
or (ii) the luckiest of 16 operators. We close both doors.

**(A) Fake-label (label-shuffle) null — is it the *true* SOZ that matters?**
Relabel a random, count-matched set of contacts (spread over ≥2 shafts) as *fake* SOZ and
rerun the whole pipeline. Repeat 500 times. If the real number is special, the fakes should
collapse to chance.

$$
p_{\text{empirical}}=\frac{\#\{\text{shuffles with median AUC}\ \ge\ \text{real median AUC}\}}{500}
$$

Result: fakes sit at **≈0.49** (95th percentile ≈0.57) versus the real **0.72–0.75**, so
**$p=0.000$** — not one of 500 random sets reached the real number. The recovery is specific
to the *true* SOZ being a tied-together diffusion community.

```python
# per shuffle: pick n_epi random contacts spanning >=2 shafts as fake SOZ, rerun LOSO
p_emp = np.mean([shuffle_median_auc >= real_median_auc for _ in range(500)])   # -> 0.0
```

**(B) Nested leave-one-patient-out — did we just pick the lucky operator?**
For each held-out patient, choose the best marker using the *other 9 patients only*, then
score it on the held-out 10th. This removes any "we chose the winner after seeing everyone".
The slow heat kernel (`heat_t5`) is selected in **9 of 10** folds, and the held-out median is
**0.716** — the same number, out of sample.

→ Both tests are baked into the figure as the grey **label-shuffle null band** (mean → 95th
percentile ≈ 0.49–0.57) behind panels (a)–(d): every green point clears it.

---

## Step 10 — Honest scope (the limits, up front)

- **δ band only** — the other rhythms carry no distant signal.
- **Seed-based** — it needs some *known* SOZ (on ≥2 shafts) to find more; it cannot start
  from zero labels or run a blank patient.
- **8 of 10** — it fails on the two hub-patients (Pat_10, Pat_15), where SOZ are hubs.
- **Clinical labels, not surgical outcome** — validated against what the clinician marked,
  not against which tissue, once removed, stopped the seizures.

Because of these, the honest framing is a **proof of principle**: *the LRG propagator carries
seizure-zone information beyond proximity and beyond hubness* — a true mechanistic claim — not
a finished diagnostic tool.

---

## Provenance — code and data

| Step | Script (`scripts/01_compute/audit/`) | Output |
|---|---|---|
| operators + LOSO score (Steps 1–8) | `audit_101_epi_marker_library.py` | `data/audit/epi_marker_library/marker_library_{per_patient,cohort}.csv` |
| null + nested-LOPO (Step 9) | `audit_102_epi_marker_verify.py` | `…/verify_{nulls,lopo}.csv` |
| why proximity is discarded (Step 3) | `audit_99_epi_discovery_yield.py` | `data/audit/epi_discovery_yield/` |
| figure | `scripts/02_preprint/preprint_39_epi_distant_marker.py` | `data/preprint/figures/all_bands/fig_epi_distant_soz_marker.pdf` |

Full verdict + supersession of the old proximity-confounded numbers:
`.agents/reports/2026-06-12_propagator-distant-soz-marker.md`.
