Slides

## Slide 1: Title Slide

## Slide 2: Transitive Inference: a puzzle of relational reasoning
Here we need to introduce here the fundamental structure and question
- The starting point is a simple exercise that each of us is able to resolve: a transitive inference cognitive task
- Here we start with the example of what the task consist of (example of the best Ultimate (frisbee) player), we have four players Federico battiston > Andrea Gabrielli > Diego Garlaschelli > Sandro Meloni. If I show you this hierarchy and then ask you who would win the following tournament (FB vs DG and AG vs SM) you already know each the winner of each direct match: your brain has encoded the relation and then has performed an inference, implicitly building a cognitive map able to infer relations previously unseen
- Understanding how brain does this practically is an impotant question because ...
- Close with the final question where and how brain stores the relational map?
Figures: a picture of ultimate frisbee, the relations between the player, the tournament bracket

## Slide 3: The sEEG dataset
To answer the previous question the starting point is the dataset
- Our colleagues at sapienza had intracranial seeg data
- transitive inference task (way more of pairs)
- epilepsy patients
- The protocol: rspre → taskLearn → taskTest → rsPost
- strenght of dataset (pros): ultra accurate recording (2 kHz, add spatial res)
- challenges of the datasets(cons)  (heterogeneity, n=10)
- We approach it though network science, well known cognition is a propery of collective organization
Figures: a seeg picture showign the shafts, the actual implants, the protocol structure 

## Slide 4: Functional Connectivy
The network scientist path toward undertanding relational features of neuronal activity is functional connectivity
- In network neuroscience one distinguish structure and function
- Functional connectivity can be studied through netowkr properties and metrics
- No unit can be seen as disconnected from the network it belong 
- Allows for organization of complicated timeseries into matricial objects
- Maybe somethign i am missing? 
Figures: structure/function, oen functional connectivity from our patients general network analyssi pipeline, 3d brain with network

## Slide 5: Higher-order features and the forgotten Multiscale
While network neuroscience is amature iedl and basic graph metrics produced a wide variety of results, there is alor unexplored and the field is hunger for more results
- Part of the approach to detect higher order features beyon the simple network metric has moved to hypergraph and intrinsic higher order network analysis
- But one fo the most interesting higher order feature is the so called multiscale properties which requires no further modelling offer beyond pairwise
- Multiscale is known to be a fundamental concept hosting and explaining many phenomena occuring in human brain (add examples) 
- To be multiscale means that there is a fundamental resolution problem linked to the fact that any single scale analysis of the system woudl lose pehnomena intrinsically driven by this and congition is one of those phenomenon that is naturally assumed to be an emergence of intertwined scales in the system 
- Structure itself live on multiple scales if oen thinks about it from neuronal population to lobes and full brain behavior, 
- This is better represented form the three axis plot ad from the fundational paper of the field of multiscale neuroscience which is the thrre axis along which the brain is intrinsically a multiscale system: time (with multiple rythms regulating different processes), spatial, from the single neuron firing to the lobes communicating and topological, where relations betwteen parts of the brain range from the local pairwise interactions to the presence of symplexes whose emergence is a collective behavior of multiple entities
- because in network science mutliscale  emthods are often difficult to control this path has led not to being able to unify under a simple interpreting multiscale lens the extraction o multiscale features from the topological axis
- we will see that there is a particularly simple framework where multiscale features can be extracted and how this apply to the detection of patterns implied by cognition and in particular transitive inference tasks
Figures: betzel multiscale 3 axis, ona small hypergraph, a figure for some results from multiscale brain analysis (suggest the most important)

## Slide 6: Revealing multiscale structures with diffusion
The tool that we will use for revealing this multiscale structure in network science is as fundamental as one  of the most basic process in nature: the diffusion process.
- Recently developed by Pablo and andrea, which are surely better statistical phsyicsit than frisbee players
- The approach is as simple imagining that on the network there is a quantity that start spreading as heat does in real space
- In discrete geometry this has a pretty simple form because this is the process described by the dyanmical equation dxdt = -L x where L is nothing bt hte network laplacian, that is D-A
- Interestingly if you let the diffusion go it turns our that this reveal the structure linking the temporal scale to the spatial one, also allowing for the definition of intrinsic scales in the system
- indeed, if one consider the so called propagator of the process, this can in turn be used to establish a relational diastance between the nodes D_ij(tau) ~1/K_ij(tau)
- the resolution scale allows for a zoomin out and allows for resolving the multiscale structure of the network, by revealing nested organization 
- The sense of the higher order here can be intended by looking to the fact that the propagator encodes higher order steps (explain this better) terms of the e^(-tau L) = sum_n(-tau L)^n/n!
Figures: video

## Slide 7: Form Mirrors Form
the intrincic ability of the brain to perform inference is wriitten in the multiscale structure of functional relatiosn between brain areas
- What if the built relational map live in the higher order structure? 
- it is not a completely brand new idea that cognition resides in multiscale feature of brain (add examples and quote )
- inspired by the relational map beign a mltirelationsal inference we bet that the trace of such process could live in the multiscale structure of the relational network
- What we are hunting here is the presence of a trace in the structure that is invisible to regular network measure and spectral clusterign mehtods on the functional connectivity, and rely on multiscale feature of the topoloy of the network
- So the practical first step for doing this is define the protol of anlysis that we intend to use especially to extract fc from this kind of activity time series data
Figures: ... (what?)

## Slide 8: Imaginary Coherence for FC estimation
- Eliminates volume conduction
- No lag components
- Accurate connectivity assessment
- Enhanced signal interpretation
- What to say else here
- Which figure?

## Slide 9: The temporal scale: Frequency Bands
Along the temporal axis of the multiscale brain we split into frequency bands cause its known that dirrefenr rytim correspond to different ...
- what to say else here? 
- which figure?

## Slide 10: The Network analysis Pipeline 
here show the pipeline as a panel where we start from timeseries of brain and we arrive at the dendrogram, then we measure the spearman correlators for see the rank relation

## Slide 11: A lasting trace
- Show here that we see trace, 
- In raw fc no heterogeneity across bands, all seems to move together
- in grassmann (spectral analaysis) is different per band
- in cohp also is different per band, but disagree with sgrassmann
- mmmmm interesting why this heterog? What is going on
- here pur results without statistical control tests
  
## Slide 12: Statistical tests: the failing of standard pairwise and spectral clustering methods → it must be multiscale (find better title)
⚠ SUPERSEDED 2026-07-13 — this raw outline predates the current per-slide files. See slides 14 (the null) + 15 (controls ladder). DRIFT RETIRED: the "drift washes out grassmann/raw fc" bullet below is void — a directional task makes a drift null degenerate with the trace, so matched-strength is the SOLE null. The current "it must be multiscale" argument is the CONTROLS LADDER (raw / strength / clustering / geodesic / spectral-resistance all fail; only cophenetic carves the bands).
Put here the controls we performed, 
- describe controls and rationale
- matched strengt partly kills raw fc
- ~~drift completely washes out grassmann and raw fc~~ ← RETIRED (drift null discarded 2026-07-12)
- controls ladder: simple pairwise, network, geodesic, AND spectral (effective resistance) all fail; cophenetic (multiscale, τ-swept) is the only band-selective read
figures of trace pre and post tests 

## Slide 13: Localizaiton of beta to ofc

## Slide 14: encoding vs inference

## Slide 15: encoding vs inference

## Slide 16: Some more about variability and heterogeneity across bands
alpha, gamma_low have large other resutls, state them

## Slide 17: localization of encoding vs inference

## Slide 18: Propagator allows for an epilepsy marker
yet not completely detection tool

## Slide 19: Outlook and thanks
- bigger datasets
- different tasks
- indepdendent from substrate
- 