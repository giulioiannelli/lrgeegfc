---
name: talk-slide-03-seeg-dataset
type: report
era: IMCOH_ABS × COHORT_N10 (ρ_sym)
slide: 3
status: draft
updated: 2026-07-11
canva: page 3 · ~35%
---

# Slide 3 — The sEEG dataset

1. TITLE
The sEEG dataset

2. MAIN CONCEPT
- Intracranial sEEG from our Sapienza collaborators (Ferraina group) — epilepsy patients implanted for clinical monitoring.
- Crucially, this is primarily an epilepsy dataset — the electrodes are placed to localize seizures, not to answer our cognitive question; the task rides on top of clinical monitoring.
- They performed a transitive-inference task (many more pairs than the toy example).
- Protocol: rest_pre → task_learn → task_test → rest_post (the four phases we exploit).
- Cognitive protocol: in task_learn they study adjacent premise pairs; in task_test they judge the non-adjacent inferences — the disc-golf exercise, scaled to many items.
- Spatial-scale bridge (why this data suits a multiscale reading): each contact sits near the neuronal scale, recording local population activity; a group of neighbouring contacts stands for an area. So the contact network already spans spatial scales — from local activity up to area–area interaction. This is the spatial axis of the multiscale brain (Betzel), embodied in the electrode layout; it pays off in slides 5–6. NB: this is a property of the *data's spatial layout* — distinct from the *topological* coarse-graining our diffusion lens does later; don't conflate the two.
- Pros: dual resolution (2 kHz in time, mm in space), a direct high-SNR intracranial signal, a within-subject rest+task design, and — above all — a built-in reach across spatial scales. Cons: heterogeneous per-patient implants (no common montage), epilepsy-driven coverage (electrodes follow the seizures, not us), sparse sampling of only clinical sites, and a small cohort (n = 10) — together these make cross-patient pooling hard and force patient-level metrics + cohort nulls rather than a group-average map.
- Precedent: the same Sapienza group showed in monkeys that transitive inference recruits specific prefrontal rhythms (Di Bello, Ferraina 2024, Commun Biol). We extend to the multiscale network in humans.
- We approach it through network science: cognition is a property of collective organization.
- A double probe: the same recordings serve a cognitive question (the trace) and a clinical one (epileptogenic tissue) — one dataset, two payoffs (the coda closes this loop).

3. ON-SLIDE TEXT
Two-block PROS / CONS layout (keep the pros–cons structure). Context strip on top; protocol carried by the timeline figure. Trim each side to the strongest 2–3 on the Canva slide if it gets crowded.

Context strip: stereo-EEG (sEEG) · epilepsy patients · n = 10 · seizure-onset zones (SOZ) marked

PROS
- spans spatial scales — local → area
- resolution: 2 kHz + mm
- direct, clean intracranial signal
- within-subject rest + task

CONS
- heterogeneous implants
- coverage follows epilepsy, not us
- sparse — clinical sites only
- hard to pool across patients

(protocol shown by the timeline figure: rest_pre → task_learn → task_test → rest_post; cognitive task = learn adjacent premises → infer the non-adjacent pairs)

4. SPEECH
To answer that, we have a rare dataset: intracranial recordings from epilepsy patients doing a transitive-inference task — the same disc-golf exercise, scaled to many items: in the learning phase they study the adjacent premise pairs, then they're tested on the non-adjacent inferences. Four phases — rest, learn, test, rest again — recorded at two kilohertz with millimetre precision. The spatial side is what makes it special: each contact sits close to the neuronal scale, reading local population activity, and once you build a network between contacts, groups of neighbouring ones stand for whole areas — so the network already reaches across spatial scales, from local activity up to areas interacting. That built-in reach is exactly what we exploit. The price is a small, heterogeneous cohort of ten, a different implant in every patient, and — crucially — coverage driven by epilepsy surgery, not by our question: the electrodes follow the seizures, which makes pooling across patients genuinely hard. This isn't the first look at inference in the brain — the same Sapienza group has shown, in monkeys, that transitive inference recruits specific rhythms in prefrontal cortex. What's new here is reading it through a network-science lens — cognition as collective organization — and across scales.

5. FIGURES
- An sEEG image showing the shafts (external / clinical).
- Cohort implants overlaid on a brain — 3D pooled pial shell, contacts as mm-space spheres, one colour per patient, three angular views (left · superior · right); 1169 contacts, n=10. → data/outputs/figures/talk/cohort_implants_overlay_3d.png (static) + …_3d.html (interactive). Gen scripts/07_figures/gen_cohort_implants_3d.py. [replaces the old flat glass-brain overlay]
- The four-phase protocol timeline — rest_pre → task_learn → task_test → rest_post. BUILD (Lane F): no paper describes our protocol, so it must be a figure. → data/outputs/figures/talk/ti_protocol_timeline.pdf.
- Raw sEEG traces (wide horizontal banner) — what a contact records: 5 channels (3 non-SOZ dark + 2 SOZ red) over an 8 s window at 2048 Hz, common gain, a continuation arrow off the right end of each trace ("the recording goes on"), time + amplitude scale bars. → data/outputs/figures/talk/seeg_traces_stacked.pdf. Gen scripts/07_figures/gen_seeg_traces.py (Pat_05 rest_pre; --dur window length, --n-ctx/--n-soz/--seed random channel draw; amplitude in a.u. — pass --uv if the export is confirmed µV).

6. REFERENCES
- Dataset (self): Iannelli, G. et al. (2026), ongoing.
- Acquisition protocol: Ricci et al. (2023) — Bambino Gesù group, collaborators on data acquisition. [Lane R: pin exact citation.]
- Same-group TI precedent (macaque): Di Bello, F., Mione, V., Pani, P., Brunamonti, E., Ferraina, S. (2024), "Prefrontal cortex contribution in transitive inference task through the interplay of beta and gamma oscillations", Communications Biology 7(1), 1715. DOI 10.1038/s42003-024-07418-5 — TI recruits β/γ in PFC in monkeys (pre-echoes our β flagship; use as motivation, not identical evidence).
- sEEG methods / good practices (better placed in the methods section or refs appendix, not the intro): Mercier, M. et al. (2022), "Advances in human intracranial EEG research, guidelines and good practices", NeuroImage 260, 119438; Wang et al. (2023). [Lane R: pin Wang citation.]
- Acquisition context (dataset facts, for methods — not on-slide): 2048 Hz (Pat_03 1024 Hz); referenced to a common white-matter channel (G2 = 2nd row of channel_labels); notch + band-pass filtered only, no further preprocessing. → also propagate to the data-layout guide.

7. CANVA STATUS
Canva page TBD — deck shifted (~+5) after the slide-2 split; the sEEG slide is no longer on pg 3 (pg 3 is now a slide-2 disc-golf beat). Content: dataset text final; references clarified (Ricci acquisition + Di Bello/Ferraina monkey precedent + Reber; Mercier/Wang → methods). Figures: implants overlay ready (Lane F), sEEG photo (Lane R), protocol timeline (Canva). Re-sync the exact page in the one-pass Canva re-map.
