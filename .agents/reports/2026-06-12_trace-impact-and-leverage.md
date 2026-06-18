---
name: trace-impact-and-leverage
type: report
era: IMCOH_ABS / COHORT_N10
status: current
created: 2026-06-12
updated: 2026-06-12
pointers:
  - .agents/guides/task-persistence-investigation/README.md
  - .agents/reports/2026-06-11_epi-soz-marker-literature-positioning.md
  - memory/task_paradigm_transitive_inference.md
---

# What the cophenetic trace is *for*: impact, positioning, and leverage

**Head (the juice).** The β-band LRG cophenetic *trace* is best presented as a
**methodologically novel descriptor on a well-precedented substrate**: nobody
has applied a Laplacian-Renormalization-Group connectivity *hierarchy* to
intracranial EEG, and nobody tracks a hierarchy change that *persists into
post-task rest* — that persistence axis is a documented gap in every dominant
comparator (gradients, manifolds, flexibility). The finding sits at the
intersection of three literatures and the honest verdict is the same in all
three: **present it as a novel multiscale-hierarchy measure with a defensible
cognitive-map / consolidation MOTIVATION and a clinical-biomarker HYPOTHESIS —
never as a validated biomarker, never as a demonstrated consolidation
mechanism, and lead every claim with its missing control.** The single
highest-value, lowest-cost next experiment is to correlate per-patient trace
magnitude with transitive-inference behavior (the Bassett-2011
"reconfiguration-predicts-learning" template).

This report is the durable capture of a 2026-06-12 deep-research pass
(workflow `wf_474423d2-722`, 102 agents, 22 verified claims / 3 killed, all
primary sources) plus two targeted clinical top-ups. The task identity
(**transitive inference**) is recorded in
[[task_paradigm_transitive_inference]].

---

## 1. Positioning verdict

| Angle (user weighting) | Where the trace stands | The honest ceiling |
|---|---|---|
| **Clinical biomarker (#1)** | A *hypothesis-generating* candidate. Best clinical claim is **specificity** (trace concentrates in a cognitive hub and is **not** the epileptogenic network), framed as an **extension of the nociferous / remote-network hypothesis**. | No surgical outcome → cannot make any "predicts decline / localizes SOZ" claim. SOZ-*prediction* by the per-node cophenetic trace is already **negative** in-repo (audit_80/81). |
| **Methods positioning (#2)** | The **persistence-of-a-hierarchy** axis is genuinely novel; LRG-on-iEEG is unworked; cophenetic comparison of brains is precedented as a principled object. | Cophenetic/dendrogram comparison is **not** demonstrably *superior* to graph metrics (that claim was refuted) — claim "principled + multiscale", not "better". |
| **Consolidation (#3)** | The strongest *story* (Tambini–Davachi persistence; Liu–Behrens replay; OFC = cognitive map), tightly matched to the transitive-inference design. | It is a **motivation, not a result**: the data show a persistent FC-hierarchy change, not replay and not a causal memory link. The replay authors themselves decline the consolidation claim. |

---

## 2. The three angles, with citations

### 2.1 Methods positioning (#2)

- **LRG is prior art, but never on brains.** Villegas et al. 2023 (*Nature
  Physics* 19:445) defined ρ(τ)=e^{−τL}/Z on metabolic / interactome /
  Internet graphs — not connectomes. The conceptual ancestor *for brain FC* is
  the De Domenico spectral-entropy lineage (Nicolini 2020 *NeuroImage*;
  Benigni 2021 *Network Neuroscience*). Applying LRG to intracranial-EEG FC
  with a cross-phase comparison is novel.
  ⚠️ Cite Villegas **only** for the ρ(τ)/ultrametric/dendrogram primitive. The
  "Kadanoff supernodes / mesoscale communities at τ\*" reading was **refuted**
  in verification and is degenerate for our fully-connected weighted case (the
  project already retired the entropy-curve rung).
- **Ultrametric/cophenetic comparison of brains is precedented.** Chung & Lee
  et al. 2012 (*IEEE TMI*): single-linkage dendrograms + Gromov–Hausdorff =
  "first unified framework for measuring brain-network differences"
  (single-linkage minimax = cophenetic distance). ⚠️ Applied to *clinical
  groups*, not within-subject task-vs-rest; the "beats graph metrics" claim
  was **refuted** — claim *principled & precedented*, not *superior*.
- **The persistence axis is the gap.** The dominant comparators measure
  reorganization *within* task/learning epochs and report **no** offline/rest
  persistence: cortical gradients (Margulies 2016 *PNAS*), manifold
  eccentricity (Areshenkoff 2024 *eLife* — resting scan was *pre*-task; no
  post-learning rest analyzed), module-allegiance flexibility (Bassett 2011
  *PNAS*). An ultrametric tree is geometrically distinct from a 1-D continuous
  gradient. **What the hierarchy view adds: which node-sets block together at
  which scales, and whether that blocking persists.**
- **The leverage template.** Bassett 2011: network flexibility in one session
  predicts the *relative amount of learning in a future session* (task
  activation did not) — the precedent that a single per-subject reconfiguration
  scalar can be behaviorally predictive. Reliability caveat: Yang 2020
  (*NeuroImage*) flags under-characterized test-retest reliability of
  multilayer flexibility — any single trace scalar must show stability.

### 2.2 Clinical biomarker (#1)

- **The benchmark you'd be measured against** (interictal connectivity
  SOZ-markers, validated on resection + Engel/ILAE outcome, AUC ≈ 0.70–0.86):
  Gunnarsdottir 2022 (*Brain*, Source-Sink Index, n=65, outcome AUC 0.86);
  Doss 2024 (*Brain*, Interictal Suppression, n=81); Narasimhan/Park 2020
  (*Epilepsia*, directed+undirected AUC 0.88); Roy 2025 (eigenvector
  centrality, n=65, LOPO AUC 0.70); Taylor/Wang 2022 (*Brain*, normative
  band-power AUC 0.75); Jiang 2022 (*Advanced Science*, interictal **resting**
  sEEG, SOZ acc 88% / AUC 0.94, seizure-outcome AUC 0.93). **All
  resting-state/ictal — none task-persistence-based.**
- **The bar to be heard at all:** beat a **spatial-sampling null** (Conrad
  2022 *J Neural Eng*, n=110: electrode location/density alone classifies SOZ
  at AUC 0.70 — "a feature must outperform a spatial null model"); match
  **static-FC properties** (Laumann/Gratton 2024 *Imaging Neuroscience* — the
  timeseries-domain sibling of our graph-domain matched-strength null); then
  cross-patient label-free generalization and replication at n≫10.
- **Why SOZ-prediction is the wrong play for this measure:** the per-node
  cophenetic-trace → SOZ classifier was built in-repo and is **negative
  cross-patient** (audit_80/81, LOPO at chance). The positive SOZ marker in the
  project is a *different* object — the propagator diffusion-community
  (audit_89–98, [[epi_propagator_diffusion_community_2026_06_08]]). **Do not
  conflate them.**
- **The winning clinical claim is specificity.** The β-trace OFC concentration
  **survives epileptic-node exclusion** (audit_83 `--epi-mode exclude`),
  licensing: *"this task-induced reorganization concentrates in orbitofrontal
  cortex and is demonstrably not an epileptogenic-network artifact."*
- **Band confound to pre-empt:** Shah 2019 (*NeuroImage:Clinical*) found
  within-resection-zone hyperconnectivity in **15–25 Hz beta** — the same band
  as the trace. A referee will ask whether the β-trace is partly β
  epileptogenic hyperconnectivity. The epi-exclusion survival is the answer,
  but it must be run as a foregrounded specificity test (→ TC2).

### 2.3 Consolidation / offline reprocessing (#3)

- **The one true mechanistic anchor.** Tambini, Ketz & Davachi 2013 (*PNAS*):
  task-evoked hippocampal patterns that **persist into post-encoding rest**
  predict subsequent memory — i.e. *a task pattern that persists offline is
  behaviorally meaningful*. Healthy subjects, within-subject, single-session —
  but it is the precedent that legitimizes the whole "trace" concept's
  behavioral relevance.
- **The relational/abstract resonance.** Liu, Dolan, Kurth-Nelson & Behrens
  2019 (*Cell*): awake post-learning rest replay follows the **learned abstract
  order**, not the literal visual experience, and stitches together transitions
  never directly seen — exactly transitive inference. ⚠️ The authors
  **explicitly decline** the consolidation/causal claim; a 2025–26 reanalysis
  questions TDLM replay-detection reliability. Motivation, not result.
- **The OFC substrate.** Schuck, Cai, Wilson & Niv 2016 (*Neuron*): OFC encodes
  a cognitive map of **unobservable** task states; decoding accuracy tracks
  performance. Origin: Wilson/Takahashi/Schoenbaum/Niv 2014 (*Neuron*). A 2025
  *NeuroImage* social-transitive-inference study (n=25) finds HPC+OFC greater
  activation **and** stronger FC during inference — a learn-a-hierarchy-then-
  infer design parallel to ours. ⚠️ fMRI + social/spatial TI, within-task —
  corroborates substrate/framing, not the persistence axis or the iEEG context.

---

## 3. The crystallized clinical framing

The cleanest, honest clinical story is **not "biomarker"** — it is a three-layer
scaffold that the trace *extends*:

1. **Mechanistic precedent (the anchor):** Tambini & Davachi 2013 — offline
   persistence of task patterns is behaviorally meaningful.
2. **Disease scaffold (well-precedented):** the **nociferous-cortex /
   remote-network hypothesis** — epileptic tissue degrades cognition in
   *distant* regions more than in itself. Penfield & Jasper (1950s, concept);
   **Ung 2017 (*Brain*): interictal spikes *outside* the SOZ impair memory
   *more* than spikes inside**; Dahal 2019 (*Brain*, n=10 — our exact n):
   long-range (~6 cm) pathological IED–spindle coupling. Cognitive comorbidity
   is the rule (60–80% of chronic epilepsy; TLE→memory, FLE→executive — Novak
   2022; Stretton & Thompson 2012).
3. **Task/population precedent:** **Reber 2016 (*Hippocampus*)** — intracranial
   recordings in presurgical epilepsy patients performing *relational/transitive
   inference*. Validates our exact task in our exact population — but via
   **hippocampal ERPs**, not OFC, not a connectivity hierarchy: it backs the
   task and the MTL hint, not the OFC-β headline.

**The novelty (double-edged).** There is **no primary precedent** for OFC-
specific involvement in epilepsy cognition (that literature is *dorsolateral*
PFC / executive; OFC-in-epilepsy ties to OCD/affect), and **none** for a
task-induced offline-persisting FC-hierarchy used as a clinical marker. So the
OFC-β-hierarchy trace is a **novel extension to a non-canonical cognitive hub**
— do **not** write "consistent with known OFC pathology." Novel = no precedent
to lean on = **higher** evidentiary burden.

**The field you cite but do not compete with** (presurgical memory-network
prediction of post-resection decline; mature, AUC 0.7–0.9, n=50–101, all with
the surgical outcome we lack): Sidhu 2015 (*Neurology*, memory-fMRI
lateralization, n=50, R²=0.43); Balachandra/Bonilha 2020 (*Neurology*,
structural connectome + hippocampal volume, n=81 UCSD→UCSF held-out, AUC 0.90);
Busby/Gleichgerrcht/Bonilha 2023 (*Epilepsia*, white-matter MTL asymmetry,
n=101, 25–33% variance); **Audrain/McAndrews 2023 (*Epilepsia*, n=72:
preoperative resting-state integration of the *to-be-resected* region with the
memory network predicts verbal-memory decline, 44% variance — the closest
resting-state analog to a trace/persistence framing)**; Doucet/Sperling 2015
(*Epilepsia*, resting graph theory → neurocognitive outcome). Cite these as the
state of the art you are **not** matching; lead any biomarker sentence with
"no surgical outcome was available; conceptual analogy, not a validated marker."

---

## 4. Test-case menu (leverage)

Ranked by value × feasibility; each names the exact in-repo artifact it consumes.

### ⭐ TC1 — Trace magnitude vs transitive-inference behavior  *(highest value; needs behavioral CSVs)*
- **Data:** per-patient trace magnitude `data/audit/wm_stratified/cophenetic_raw_per_patient.csv`
  (`obs_stat` / `obs_z`, per patient × band) × per-patient TI performance
  (accuracy on novel non-adjacent pairs / RT / recall — obtainable, not yet in repo).
- **Design:** Spearman across n=10, per band, a-priori prediction **β carries it**.
- **Template:** Bassett 2011 (reconfiguration predicts learning); mechanism anchor Tambini & Davachi 2013.
- **Referee demands:** partial out node-strength + within-baseline drift; show test-retest stability of the scalar; pre-register the β-direction. **n=10 is the hard limit** — report effect size + CI, treat p as secondary, frame confirmatory not discovery.
- **Plug-in point:** add `load_task_performance(patient)` to `src/lrg_eegfc/utils/io/patient.py`; natural file `data/raw/stereoeeg_patients/{PAT}/task_performance.csv`.

### TC2 — Trace vs SOZ spatial dissociation (specificity)  *(runs today on repo data; recommended first build)*
- **Data:** build the merged per-contact table that does **not** yet exist —
  join `data/audit/lrg_localization_anatomy/per_trace_leaf.csv` (`rho_demeaned`,
  Desikan-Killiany `region`, MNI x/y/z) with SOZ flags (`load_epileptic_nodes`)
  and `load_channel_regions`.
- **Design:** does the OFC β-trace concentration sit **inside / adjacent to /
  away from** each patient's SOZ; does per-contact trace correlate with
  distance-to-SOZ?
- **Why:** operationalizes the §3 nociferous-extension story (trace in a
  cognitive hub, *spares* the seizure network) and kills the Shah-2019 β-confound.
- **Referee demands:** benchmark against the **Conrad spatial-sampling null
  (AUC 0.70 floor)** before attributing to physiology.

### TC3 — Does the hierarchy view see persistence the standard tools miss?  *(runs on repo data; defends methods novelty)*
- **Design:** compute cortical-gradient eccentricity, module-allegiance
  flexibility, and edgewise FC change **cross-phase**, and show they do **not**
  carry the rest_post-persisting signal that cophenetic distance does.
- **Why:** earns the methods-novelty claim empirically instead of asserting it.

### TC4 — Band-specificity + minimal external validity  *(controls referees force)*
- Why β; is α / low-γ independent or a β-spillover; and name the replication
  path (a second TI/relational-memory iEEG cohort — the only credible external
  validity absent Engel/ILAE here).

---

## 5. Brutal-honesty ledger (lead with these)

1. **No surgical outcome, no demographics** — biggest gap vs every clinical
   benchmark; trace = triage / hypothesis only.
2. **OFC localization flipped 3× in audit history.** Surviving claim: β
   cophenetic-trace hotspot in OFC (system scale, R=1000, BH q≈0.009–0.013,
   shaft- and LOO-robust). **MTL/hippocampus = sub-threshold, LOO-fragile,
   Pat_02-driven hint — never a co-headline** ([[localization_audit_plan_2026_05_29]]).
3. **n=10 + sEEG spatial sparsity** — location alone is partly predictive
   (spatial null AUC 0.70).
4. **Exact TI stimuli/timing undocumented** — weakens fine-grained mechanistic
   claims until the protocol is obtained.
5. **Consolidation is a motivation, not a result** — persistence of an
   FC-hierarchy change, not replay, not a causal memory link.
6. **Supporting cognitive literature is domain-mismatched** (fMRI, social/
   spatial TI; OFC-cognitive-map work is healthy subjects) — licenses framing
   and substrate, not direct prior art.
7. **Task-persistence as a clinical marker is unprecedented** — no precedent to
   lean on; higher burden, not lower.

---

## 6. Citation appendix (verified primary sources)

**Methods / hierarchy**
- Villegas et al. 2023, *Nature Physics* 19:445 — LRG primitive only. https://www.nature.com/articles/s41567-022-01866-8
- Nicolini et al. 2020, *NeuroImage* — scale-resolved spectral entropy. https://www.sciencedirect.com/science/article/pii/S1053811920300902
- Benigni et al. 2021, *Network Neuroscience* — multiscale spectral entropy. https://direct.mit.edu/netn/article/5/3/831/102505/
- Chung & Lee et al. 2012, *IEEE TMI* — single-linkage + Gromov–Hausdorff = cophenetic comparison. https://pages.stat.wisc.edu/~mchung/papers/tr_228.pdf
- Margulies et al. 2016, *PNAS* — principal cortical gradient. https://www.pnas.org/doi/10.1073/pnas.1608282113
- Areshenkoff et al. 2024, *eLife* — manifold eccentricity; no post-learning rest. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11198988/
- Bassett et al. 2011, *PNAS* — flexibility predicts future learning. https://arxiv.org/pdf/1010.3775

**Clinical — SOZ connectivity markers + controls**
- Gunnarsdottir et al. 2022, *Brain* — Source-Sink Index, n=65, outcome AUC 0.86.
- Doss et al. 2024, *Brain* — Interictal Suppression Hypothesis, n=81.
- Narasimhan/Park et al. 2020, *Epilepsia* — directed+undirected AUC 0.88.
- Shah et al. 2019, *NeuroImage:Clinical* — n=27, β 15–25 Hz within-zone hyperconnectivity surviving a spatial null.
- Roy et al. 2025, *Front Netw Physiol* — eigenvector centrality, n=65, LOPO AUC 0.70.
- Taylor/Wang et al. 2022, *Brain* — normative band-power, AUC 0.75.
- Jiang et al. 2022, *Advanced Science* — interictal resting sEEG, SOZ AUC 0.94, outcome AUC 0.93. https://pmc.ncbi.nlm.nih.gov/articles/PMC9218648/
- Conrad et al. 2022, *J Neural Eng*, n=110 — spatial-sampling null, AUC 0.70 floor. https://pmc.ncbi.nlm.nih.gov/articles/PMC9590099/
- Laumann/Snyder/Gratton 2024, *Imaging Neuroscience* — static-FC-matched null. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12315734/

**Epilepsy cognition / nociferous / remote-network**
- Ung et al. 2017, *Brain* — extra-SOZ spikes impair memory more than intra-SOZ. (via https://pmc.ncbi.nlm.nih.gov/articles/PMC8746065/)
- Dahal et al. 2019, *Brain* 142(11):3502 — n=10, long-range IED–spindle coupling. https://academic.oup.com/brain/article/142/11/3502/5566384
- Novak et al. 2022, *J Clin Med* — 60–80% cognitive impairment in chronic epilepsy. https://pmc.ncbi.nlm.nih.gov/articles/PMC8746065/
- Stretton & Thompson 2012, *Epilepsy Behav* — TLE executive dysfunction. https://pmc.ncbi.nlm.nih.gov/articles/PMC4251008/
- Stoub et al. 2019, *Hippocampus* — hippocampal disconnection predicts memory, n=50. https://pubmed.ncbi.nlm.nih.gov/28888031/
- Reber et al. 2016, *Hippocampus* 26(1):54 — iEEG epilepsy, relational/transitive inference (hippocampal ERPs). https://pubmed.ncbi.nlm.nih.gov/26136107/

**Consolidation / replay / OFC cognitive map**
- Tambini, Ketz & Davachi 2013, *PNAS* — task patterns persist into post-encoding rest, predict memory. https://www.pnas.org/content/110/48/19591
- Liu, Dolan, Kurth-Nelson & Behrens 2019, *Cell* — awake replay reorganizes toward learned abstract order. https://www.cell.com/cell/pdf/S0092-8674(19)30640-3.pdf
- Schuck, Cai, Wilson & Niv 2016, *Neuron* — OFC cognitive map of (unobservable) task space. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5044873/
- Wilson, Takahashi, Schoenbaum & Niv 2014, *Neuron* — OFC-as-cognitive-map origin. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4001869
- *NeuroImage* 2025 — social transitive inference, HPC+OFC activation + FC. https://www.sciencedirect.com/science/article/pii/S1053811925006834

**Presurgical memory-network prediction (engage, don't compete)**
- Sidhu et al. 2015, *Neurology* — memory-fMRI lateralization, n=50, R²=0.43. https://pmc.ncbi.nlm.nih.gov/articles/PMC4408284/
- Balachandra/Bonilha et al. 2020, *Neurology* — SC + hippocampal volume, n=81, AUC 0.90. https://pmc.ncbi.nlm.nih.gov/articles/PMC7455364
- Busby/Gleichgerrcht/Bonilha 2023, *Epilepsia* — WM MTL asymmetry, n=101, 25–33% variance. https://pmc.ncbi.nlm.nih.gov/articles/PMC10640663/
- Audrain/McAndrews 2023, *Epilepsia* — to-be-resected-region integration, n=72, 44% variance. https://pubmed.ncbi.nlm.nih.gov/37643922/
- Doucet/Sperling 2015, *Epilepsia* — resting graph theory → neurocognitive outcome. https://onlinelibrary.wiley.com/doi/10.1111/epi.12936

**Refuted in verification — DO NOT cite**
- LRG Kadanoff-supernode coarse-graining authorizing the dendrogram (0-3).
- GH-dendrogram comparison *beats* graph metrics (0-3).
- "FC adds minimal value beyond spatial geometry" (1-2, contested — cite neither direction).
