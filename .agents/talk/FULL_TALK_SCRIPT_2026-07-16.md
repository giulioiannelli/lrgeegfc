---
name: talk-full-script-2026-07-16
type: report
era: IMCOH_ABS × COHORT_N10 (mst@0.20)
status: active
updated: 2026-07-16
---

# Full talk script — 18 slides · 4774 spoken words · ~37–32 min (130–150 wpm)



## Slide 1 — Title  —  ~64s (139 w)

Good morning everyone, I'm Giulio Iannelli PostDoc at CREF and today I present to you some fresh results comingfrom last year work in the context of an European project about rvealing higher order features of Human Brain with Network neurocience. In paritcular today I am going to show you how cognition — and disease — reshape these neurological networks leaving traces that can only be detected with appropriate multiscale tools.
I beg pardon to the neuroscientific audience for my brutality in exposing the neurological side of this work, also because this is ongoing work with the Sapienza neuroscience group and i lead the methodology leveraging the expertise gathered on the network side of the field at my research center; it is not finished, so I'm genuinely after your questions and any solid neurophysiological interpretation you can suggest we might have overlooked.


## Slide 2 — Transitive inference: a puzzle of relational reasoning  —  ~93s (202 w)

Six paragraphs, one per sub-slide (copy one to each):

beat 1 — Let's do a quick transitive-inference exercise together. Take disc golf.

beat 2 — Here are our four players, ranked by ability. Federico — a world-champion ("iridato") Ultimate player, raised on pasta and good Italian principles. Andrea — the master of master equations. Diego — the pioppino-mushroom eater. And Sandro — our Italian Ulysses, home from Spain. In ability it goes: Federico beats Andrea, Andrea beats Diego, Diego beats Sandro. Memorise that.

beat 3 — Now a knockout tournament. Notice the first-round matchups pair players who never actually met — Battiston against Garlaschelli, Gabrielli against Meloni. Who wins?

beat 4 — You called it instantly — Battiston and Gabrielli — even though I never showed you those matches. Those first-round results are pure inference.

beat 5 — And the final, Battiston versus Gabrielli, you were effectively shown — so: Battiston, champion. Here's the point: without being shown every match, you inferred them all. You didn't answer from a lookup table of stored pairs — you answered from a structured cognitive map that organizes what you learned so you can resolve comparisons you were never explicitly trained on. That map is a kind of geometry of abstraction: the shape relational reasoning takes in the brain.

beat 6 — So where, and how, does the brain build that cognitive map? That's the question the rest of this talk chases.


## Slide 3 — The sEEG dataset  —  ~91s (197 w)

To answer that, we have a rare dataset: intracranial recordings from epilepsy patients doing a transitive-inference task — the same disc-golf exercise, scaled to many items: in the learning phase they study the adjacent premise pairs, then they're tested on the non-adjacent inferences. Four phases — rest, learn, test, rest again — recorded at two kilohertz with millimetre precision. The spatial side is what makes it special: each contact sits close to the neuronal scale, reading local population activity, and once you build a network between contacts, groups of neighbouring ones stand for whole areas — so the network already reaches across spatial scales, from local activity up to areas interacting. That built-in reach is exactly what we exploit. The price is a small, heterogeneous cohort of ten, a different implant in every patient, and — crucially — coverage driven by epilepsy surgery, not by our question: the electrodes follow the seizures, which makes pooling across patients genuinely hard. This isn't the first look at inference in the brain — the same Sapienza group has shown, in monkeys, that transitive inference recruits specific rhythms in prefrontal cortex. What's new here is reading it through a network-science lens — cognition as collective organization — and across scales.


## Slide 4 — Functional connectivity  —  ~57s (123 w)

How does a network scientist get at relations between brain areas? Functional connectivity. For every pair of contacts we take their time series and compute a coupling estimate — how much they fluctuate together — and that pairwise coupling is a weighted graph. Here it is for four of our patients, across the four phases of the experiment — the same functional object every time. It is genuinely functional: it does not just retrace the anatomy — strong coupling shows up even between regions with no direct structural link. On that graph, decades of network neuroscience have found real organization — small-world, a rich club of hubs, communities that trade off integration against segregation. A mature, productive toolkit. Our contribution is what we do next to that graph.


## Slide 5 — Multiscale organization: a higher-order feature  —  ~70s (151 w)

The brain is multiscale by nature, organised along three axes at once. [right figure] Time, through its rhythms; space, from single neurons up to whole lobes; and topology, from pairwise links up to communities and whole-network structure. Read it at a single scale and you throw away whatever lives across scales — and cognition is assumed to be exactly that kind of cross-scale emergence.

Now, the frontier past simple pairwise connectivity is higher-order. One road is explicit — hypergraphs, simplicial complexes — you add many-body structure by hand, and it already captures real phenomenology: decoding tasks, fingerprinting individuals, tracking behaviour. But there's a second, quieter higher-order that needs no extra modelling at all: multiscale organisation itself. Coarse-grain the graph you already have and collective, beyond-pairwise structure appears — groupings no single edge carries. That's an implicit higher-order feature of the topology. What's been missing is a simple, controllable way to read it — and that's what we bring.


## Slide 6 — Why it matters: four literature anchors, each → a facet of our work (text-forward)  —  ~150s (326 w)

Let me tell you why this matters. Humans are extraordinarily good at something machine learning, for all its data, still isn't: rich inference from very little — exactly the tournament you solved a moment ago. That few-shot relational reasoning is a hallmark of human intelligence, and it's worth knowing how the brain pulls it off.

Three ideas from the literature [gesture the header column, right] tell us where to look. First, cognitive maps. The brain doesn't store a list of facts — Behrens, in "What is a cognitive map?", calls it a structured representation of how things relate, organised across scales. And it's general: Son and FeldmanHall show the same map organises social worlds too [point: the social-network figure, bottom-left] — inferring ties you never saw from the few you did. So the object we're after is a multi-scale relational structure.

Second, offline consolidation. Ellenbogen and Walker, 2007 — the title says it, "human relational memory requires time and sleep": the inferred order isn't there right after learning, it's assembled offline, at rest. That forces the distinction this whole talk turns on — what you were shown, encoding, versus what you worked out, inference.

Third, replay. Liu and Behrens, 2019 — human replay spontaneously reorganises experience, learned structure reactivated and reshuffled at rest [point: the replay schematic, bottom-left]. So there's a resting trace to go and look for.

But all of that was shown through behaviour and theory. Nobody has read that offline, multi-scale reorganisation directly in the brain's network. That's our step — from detecting that something changed to characterising what kind of change [point: detect-vs-characterize, bottom-right]: not one did-it-change number, but a profile across scales — flat, peaked at a characteristic scale, or emerging only when you zoom out.

And it isn't only cognition. Kramer and Cash frame epilepsy as a disorder of cortical network organisation that spans scales — so the same lens reads the epileptogenic network, in the very same recording. One catch I'll repeat all talk: whatever we find counts only if it beats the boring explanation — that the change simply sat where the connections are already strongest.


## Slide 7 — Form mirrors form  —  ~102s (221 w)

Everything on the last slide — consolidation, replay, the cognitive map — is the brain doing its work offline, at rest. Here's the link that turns it into a bet you can test. When the cortex sorts out what it learned, it doesn't just tuck facts away — it lays them down as a *hierarchy*. That's McNaughton's point, in the line above [gesture the McNaughton subtitle]: knowledge is extracted from memory offline and laid down hierarchically.

Now look at our task — it *is* a relational hierarchy, a chain of ordered relations [gesture the small task tree, top-right]. So the two forms should meet: if consolidation lays knowledge down as a hierarchy, and the task already is one, the resting network afterwards should carry that same form. That's the wager — form mirrors form. We're not guessing it's hierarchical; the consolidation literature predicts it.

What we *add* is a single move: read that hierarchy as a *multiscale* object — the whole nested shape across scales, not a single edge and not one fixed partition. Picture the functional network on the brain giving rise to this tree of relations [gesture brain → dendrogram]. If we're right, we don't just learn that something moved at rest — we'll say which frequency bands carry the trace, and at what scale it lives, which a single number never could.

First, though, a methods problem: how do we even measure connectivity in these signals without fooling ourselves?


## Slide 8 — Imaginary Coherence for FC estimation  —  ~78s (170 w)

How do we measure connectivity between two contacts? We start from coherence — at a given frequency, do two signals hold a consistent phase relationship? We build it from the cross-spectrum, the Welch method, and normalize.

The trap, especially in intracranial recordings: volume conduction. One source bleeds into many contacts at the same instant — zero lag — and zero-lag mixing is purely *real*. Ordinary coherence counts that leakage as a connection that isn't there.

Nolte's fix, in 2004: keep only the *imaginary* part. Zero-phase leakage has no imaginary component, so it cancels — what's left is genuine, time-lagged coupling. And Bastos and Schoffelen confirm it directly [gesture]: ordinary coherence lights up with field-spread artefacts, imaginary coherence stays clean — robust even for contacts on the same shaft, where leakage is worst.

So our connectivity, in one line: per band, take the magnitude of imaginary coherence and average across the band. That non-negative number, W-i-j, is the edge weight — it keeps the strength, drops the direction, and being non-negative is exactly what the diffusion step needs.


## Slide 9 — The temporal scale: frequency bands  —  ~54s (117 w)

The temporal axis of the multiscale brain is its rhythms, so we split the signal into six frequency bands — from delta below four hertz, up through beta, and gamma all the way to three hundred hertz. That top end is unusual: scalp EEG can't see high-frequency activity that far, but our clean intracranial contacts can, so we keep it. Splitting into bands isn't bookkeeping — it's a control. A genuine oscillatory effect should be selective for a band; a broadband change — overall arousal, a drift in signal quality — hits every band the same. So band-specificity, when we see it later, is itself a piece of evidence that the effect is real. We build one functional network per band.


## Slide 10 — Revealing multiscale structure with diffusion  —  ~179s (388 w)

Our lens is diffusion — heat spreading along the network's connections, governed by the graph Laplacian. [gesture: equations, left] The propagator, the heat kernel, tells you how much reaches one contact from another after a diffusion time tau.

Let it flow [gesture: FC matrix / play], and how long you let it run sets how coarse the communities you read. On one implant [gesture: micro-meso-macro brain panel]: briefly, dozens of little single-population clusters; longer, they merge into a handful of mid-sized groups, then two big blocks — micro to meso to macro, which is what "multiscale" actually means (illustrative on the dense graph, not a data claim). And tau is the diffusion time, a scale parameter, not the neural time axis — that's the frequency bands, from a moment ago.

And the higher-order structure comes for free. [gesture: propagator series] The propagator sums walks of every length, so multi-step paths and communities emerge from the same pairwise graph, by linear diffusion — a higher-order graph feature, never nonlinear, no simplices built by hand. That's the "no extra model" road from earlier, cashed out.

From the propagator we read a communication distance between contacts, and linkage turns it into a hierarchy. [gesture: dendrogram] The tree is the coarse-graining written out — cut it low for many small communities, high for a few large ones.

Now, does that scale mean anything physical? [gesture: right-hand curve] Here's the honest answer, and it's read straight off community size: as the scale grows, the co-diffusing communities grow, and we can measure how much of the brain each one physically spans — that's this curve, the average extent of the community a contact belongs to, against the diffusion scale. To the left of this line you're still reading raw, single edges; cross it, where real diffusion communication begins, and the groups are already spread over the implant — a few centimetres across for most bands — not compact spots. Push it further and they merge out to the whole implant, about nine centimetres. So the scale does earn a physical footprint — local groups out to global — but through community size, as a communication ruler, never a millimetre resolution zoom. These are functional groups, not proximity clusters — precise from how we just measured connectivity, a couple of slides back.

So diffusion resolves the topology axis — from pairwise edges up to nested communities. [gesture: brain+tree 3D scene] With the bands giving time and the contacts giving space, our lens now covers all three axes of the multiscale brain.


## Slide 11 — The network-analysis pipeline  —  ~73s (158 w)

Here's the whole machine on one slide, on one example patient. From the raw sEEG we build the per-band imaginary-coherence network — and now the step you haven't seen yet. That network is fully connected: every contact couples to every other, and a fully-connected graph has essentially one scale — the diffusion has nothing multiscale to resolve. So we sparsify: keep the maximum spanning tree plus the strongest fifth of the edges — a connected, cycle-rich backbone. On THAT backbone the diffusion propagator, e to the minus tau L, genuinely telescopes: as we sweep the diffusion time the specific heat C of tau grows a mesoscale ladder — several peaks, several scales — instead of one collapse; the inset shows the propagator spreading as tau increases. The last panel reads that diffusion geometry into a cophenetic tree — the hierarchy of the network at scale. How we turn two of those trees, one per phase, into a single trace number — ρ — is the next slide.


## Slide 12 — The Measure  —  ~142s (308 w)

The pipeline turned each phase into a hierarchy; the measure asks one thing — how similar are two of those hierarchies, at a given scale. Three steps, all on the left. First, read a distance off the tree: for two contacts i and j, follow their branches up to where they first join — their lowest common ancestor — and take the height of that merge. That's the cophenetic distance, D-coph — one number per pair, and scale-aware, since the tree changes with the diffusion scale τ. Second, a displacement — how far a phase pushed that fingerprint from a reference X. For the trace X is rest-pre, so Δ-task is task minus rest-pre and Δ-rest is rest-post minus rest-pre. Third, the measure: rank-correlate those two displacements — Spearman, so it's the ordering of pairs that must agree — giving rho-coph, from minus one to one. If the task's reshaping predicts the rest-after reshaping, the trace survived. We take the baseline from two halves of rest-pre, a different half on each side so they share no noise, then swap and average — the symmetrized rho-coph. The colour on the right is that number made visible: each contact shaded by how well it kept its cophenetic place — green kept, grey lost — and the bar under each tree is the single number, read at every scale, against a matched-strength null. One thing about that bar: its top is not one. Split rest-pre into two halves and correlate the hierarchy with itself and you still get only about a half — each tree is estimated from finite, noisy data, so that self-reproducibility, not one, is the most any phase can score, and a perfect trace lands there. So a rho-coph near a half is read against that reproducibility and the null, never against one. Next, we watch it move: one patient, how each population's place holds or shifts across the phases.


## Slide 13 — Task-induced taxonomy of neuronal populations (trace / anchor / reset / reorganized)  —  ~67s (146 w)

Each of these contacts picks up a millimetre-scale population of neurons, and the hierarchy I've been showing you is a hierarchy of those populations. So the question is what happens to a population's place in that hierarchy. Before any statistics, here's what we're actually looking at. Take one patient, one band, and lay the hierarchy from the task next to the hierarchy from the rest afterwards; link the same contacts between them. Lines that stay parallel mean the order was kept; lines that cross mean it reorganized. Four things can happen. It can be reorganized by the task and hold into rest — that's a trace, what we're after. It can reorganize and then revert — a reset. It can stay put throughout — an anchor. Or it can keep changing in every phase, never settling — reorganized. These four words are the vocabulary for the rest of the talk.


## Slide 14 — A lasting, multiscale trace  —  ~234s (508 w)

Two slides ago we built one number — rho — how similar two phases' trees are at a given scale, read against a matched-strength null. Now we spend it, on the question this whole talk is about: after the task, does its reorganization survive in the resting tree? The short answer is yes — read the resting hierarchy afterwards and it doesn't snap back; its tree still looks like the task's. Here's that made quantitative, three ways: how strong the trace is, where it's shared across patients, and how it varies from patient to patient.

[Fig 3 — the curves: HOW STRONG.] This is rho itself — for every band, how much the rest tree still resembles the task tree, swept from fine to coarse scale. Read two things. The height: beta rides highest, its rest tree is the most task-like at every scale; alpha next; the others lower, some sitting right on the grey null floor where there is nothing above chance. And the dots: filled when the whole cohort clears the matched-strength null at that scale, open when it doesn't. Beta is filled all the way across; alpha across most of the range.

[Fig 1 — the map: WHERE IT'S SHARED.] Compress those dots into a map — band down the side, scale along the bottom, colour is how consistently the cohort agrees, grey is below the null — and the shape of the trace jumps out. Beta is a solid block across the entire width: the cohort agrees at every single scale, so the trace is scale-invariant — it doesn't live at one resolution, it's architectural. Alpha is broad too, strongest through the mesoscale. Delta and high-gamma clear only in patches at the extremes — the fine, near-raw end and the coarse-collapse tail, the least trustworthy ends of the axis. Theta and low-gamma stay grey. So the consistent, cohort-wide trace is selective — and, crucially, grey does not mean nothing happened.

[Fig 2 — the violins: HOW IT VARIES.] Here's what the grey hides. For each band, all ten patients, plotted as distance above their own null. Look at the widths — they're huge. Even in the grey bands individual patients sit far above their null; low-gamma reaches point-six, a very strong trace. The traces are there — they just don't point the same way from patient to patient, so the cohort median doesn't hold. That's the open ring: the patients disagree, never "nobody traced". Beta is the exception — filled even at its weakest scale. And the ring reads both ways: a median a hair above the line is still open, because being above your own null is not the same as the whole cohort agreeing.

So — the result, and what it means. The task leaves a lasting, offline trace in the tree: the reorganization outlives the task. And it is band-selective. Beta carries the shared, cohort-wide trace at every scale; alpha across most of them; the other bands reorganize too, sometimes strongly, but each patient in their own way, so they don't cohere into a cohort trace. The trace we can stand behind is consistency — the same reshaping across patients — not mere presence. What drives that per-patient heterogeneity, and where in the brain the trace sits, come next.


## Slide 15 — Detection and discrimination  —  ~206s (446 w)

So we have a lasting trace, in alpha and beta — but everything we've said about it, that it's band-selective and multiscale, we read through one lens: our hierarchy. The fair question is whether that lens is doing any real work. Is this genuinely a *multiscale* reorganization — or would any ordinary network measure have found the same thing? Detection is easy; the real question is discrimination. So on the same data, against the same matched-strength null, we put three lenses side by side [top radials] — raw pairwise edges, spectral clustering (the Grassmann distance between the two graphs' eigen-subspaces), and our multiscale hierarchy — with the whole standard toolkit behind them [the matrix].

Beta first — it's the easy one, and it makes the opposite point too. Raw pairwise connectivity, the left radial, catches beta — but look how broadly it fires. The stem plot below says it in magnitude: raw is a broad pedestal, positive on every band, even theta; our cophenetic read is a sharp alpha–beta peak, and theta actually drops *negative*. That's detection versus discrimination in a single band — raw's pedestal leans positive on theta, the hierarchy suppresses it. Raw's breadth isn't sensitivity, it's non-selectivity: beta is a genuine low-order convergence everyone agrees on, but raw lights up bands the hierarchy correctly zeroes. Raw detects; it doesn't discriminate.

Alpha is the real test, and the spine of the slide. In the matrix, read the alpha column. Clustering, shortest paths, effective resistance — all blind. Raw does flag alpha, but only fragilely: it's cohort-significant, then it collapses the moment you leave one patient out — it hinges on a single subject. Only the bottom row, the multiscale hierarchy, resolves alpha *robustly*, surviving leave-one-out. Alpha is a mesoscale reorganization no single-scale read holds onto.

The field's sophisticated tool misses it too — that spectral-clustering measure. It runs on the same Laplacian we build our tree on, but it freezes one subspace dimension, a single global scale, instead of diffusing across all of them. It fires instead on low-gamma — lit by exactly the two global methods, Grassmann and the graph geodesic, and nothing else. Low-gamma is a whole-graph shift, not a multiscale reorganization; single-scale spectral clustering mistakes that global change for structure.

So why doesn't our all-scales tree flag that global low-gamma too? Because we compare phases by the tree's *ranks* — who merges before whom — not by distances, and ranks ignore any across-the-board drift. The hierarchy isn't a superset that sees everything plus more; it's a relational *lens* — it resolves the relational alpha and isn't fooled by the global low-gamma. Not "only we see a trace" — the edges see beta, spectral clustering sees beta — but we *resolve* the reorganization the standard tools miss and the field standard mis-reads.


## Slide 16 — Encoding vs inference: the order it never saw, held offline  —  ~264s (572 w)

So the trace is real, it's multiscale, and it's carried by alpha and beta. But step back to what the task actually was. [callback: the disc-golf bracket, slide 2] You were shown a handful of head-to-head results — adjacent pairs — and from those you worked out an order you were never shown: the matchups between players who never met. Two different things happened in your head — the pairs you were *shown*, and the order you *inferred* — and right at the start we said the interesting one, the inference, is the part the brain is thought to build offline, at rest. So the sharp question is: this trace that rest holds onto — is it just the pairs the brain saw, or does it also keep the order it never saw but worked out? Encoding, or inference?

The four phases let us separate them, and here is the space they define. [gesture: the 3-D portraits] One axis is how much the hierarchy moved from rest to learning, while the pairs were on the screen — that is encoding. The other is the *further* move from learning to test, the structure the pairs never showed — that is inference. Each patient's network traces a path through that space, and we ask, by rank correlation across every pair of contacts: did the task's reshaping predict what rest kept? For inference we partial encoding out, so what's left is only what the inferred order added, beyond the pairs it saw.

And one honest detail that makes it trustworthy: we never measure a task-shift and a rest-shift against the *same* baseline — we split the pre-task rest into two independent halves, cross-pair them, and average. That symmetrisation is what stops a shared baseline from manufacturing a correlation; drop it and the number more than doubles.

And rest keeps both — you can read it straight off the held tree. [gesture: the tinted hierarchy] Some of its branches are encoding-origin, in amber; others are inference-origin, in teal — both kinds survive into rest. Encoding isn't a nuisance we scrub away: the pairs the brain saw genuinely stick — in beta at every one of the sixteen scales, and in alpha too, where even ordinary tools already see it.

Now the inferred order — the thing we actually care about — and this is the real result, read against the matched-strength null across diffusion scales. [gesture: the scale curves] The amber curve, encoding, rides high across the whole scale axis in beta. The teal curve, inference, does something different and telling: at the finest scale it sits right on the null — it isn't there — and it lifts clear of the null only as the network coarse-grains, in *both* bands, beta and alpha, seven of the sixteen scales each. That's the correction to the older version of this talk — it was never beta alone. And a signal that is absent at the fine scale and appears only once you zoom out is the signature of an abstraction — the very geometry of the cognitive map we opened with — not a replayed detail.

One honesty point: raw connectivity already detects the beta inference signal, because the hierarchy is built from the same edges — so detection is tied by construction. What the multiscale read adds is the *form* — that this persistence lives as a coarse re-grouping. And it's held, not replayed: a sustained state across the ten minutes of rest, spread across the network with no single hotspot, and it doesn't grow with how long the phases ran. Rest keeps not just what the brain saw, but the order it built.


## Slide 17 — From Cognition to Pathology  —  ~100s (216 w)

And the same tool pays a second dividend — a clinical one. Same propagator, same backbone that built the cognitive hierarchy — only now, instead of asking how the hierarchy persists across phases, we drop heat on a few known seizure-onset contacts in one recording and watch where it spreads. Score each contact by the seed-heat it receives, node strength removed so it's relational, read at a coarse scale so the heat reaches the onset contacts on *other* electrodes, not just neighbours — the distant discovery you see on the brain. The seizure zone lights up as its own co-diffusing community, above a null that scrambles the zone but keeps each node's connectivity — strongest in delta and low-gamma, AUC around point-eight.

Deploy delta alone, or fuse all six bands into a leave-one-patient-out detector, and it lands at median AUC point-eight-five, nine of ten patients above chance — a triage marker, not a diagnostic, its top-five shortlist about seven times base rate. The one honest exception, a single right-hemisphere implant at chance, is the *same* patient that's our cognitive outlier — a hub, not a community. And the two reads sit on different tissue: beta spares the seizure zone, delta marks it and carries no cognitive trace at all. One operator, one backbone — a cognitive lens that doubles as a clinical one.


## Slide 18 — Take-homes · ongoing · outlook  (the compound closer)  —  ~178s (386 w)

Let me pull it together. One diffusion operator, one null, reads the brain across its three axes at once — time in the frequency bands, space in the electrodes, topology in the hierarchy. And remember how we opened: cognitive maps, offline consolidation, replay, learning from little data — all of it established by behaviour and by theory, none of it read directly in the network. That's the move. We turned detecting a change into characterizing it: which band carries it — alpha and beta; at which scale — beta at every scale, alpha at the mesoscale; and a higher-order, coarse structure the edges simply miss. The trace is held — a sustained resting state, an offline consolidator, not a fleeting burst. And what rest keeps isn't just the pairs the brain saw, it's the order it inferred — a re-grouping that only appears when you zoom out, the shape of an abstraction. Then the same operator pays a clinical dividend: read inside one recording it localizes the seizure zone, triage-grade, and the bands split cleanly — cognition in alpha and beta, epilepsy in delta, low-gamma and beta, on different tissue. And none of it counts unless it beats the boring explanation — connection strength.

What are we working on now? Per-patient traces in the other bands, which are strong but don't yet line up across people; higher-order markers beyond pairwise edges; whether the low-gamma seizure signal is a trace or just stable tissue; exactly what the simpler metrics miss; and resolving whether beta has any spatial signature at all.

And where it goes — this is the part I'm most excited about. We want behaviour, to tie the characterization to how well someone actually reasons. But here's the twist: even without behaviour, a patient who sits far from the cohort — especially where dropping them flips a band's signature — might be telling us about their own pathology. The fluctuation could be the diagnostic. Then richer, nonlinear decompositions; bigger and more varied data; a version of the trace that doesn't depend on this substrate. And finally, full circle: understanding how the human brain builds so much structure from so little data is exactly the kind of mechanism that could inspire more data-efficient machines. Thank you — to our Sapienza collaborators for the data and the task, and to the LRG authors whose framework made all of this possible.
