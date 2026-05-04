---
name: terminology
type: guide
era: CROSS_ERA
status: current
created: 2026-04-30
pointers:
  - CLAUDE.md
  - .agents/guides/04_rules/never-always-list.md
---

# Terminology — multiscale module behaviour across phases

**Single source of truth for the words we use to describe how a
structural feature behaves across the four phases `{rest_pre,
task_learn, task_test, rest_post}`. Read before writing any report,
figure caption, or chat reply that discusses cross-phase behaviour.**

The previous loose use of "persistence" for both "task-induced
fingerprint that lingered" and "module unchanged across phases"
muddied every result. Those are nearly opposite phenomena. The
canonical taxonomy below pins down four behaviours; pick the right
word every time.

---

## The four canonical behaviours

### 1. Trace (or "task-trace persistence")

**Definition.** Reorganized in `task_test` AND stays reorganized into
`rest_post`. The classical "task left a fingerprint" pattern.

| RPre | TT | RPost |
|---|---|---|
| baseline X | reorganized to Y | still close to Y (not X) |

**This is the headline phenomenon we are hunting.** Our raw-FC
triangle test (`T_d = d(TT, RPost) − d(RPre, TT) < 0`) is a *trace*
test. The previous loose "persistence in 4 of 6 bands" framing was a
misnomer.

**Synonyms (use these too).** "task-trace", "trace effect",
"reorganized persistence due to task", "task-induced persistence".
**Avoid:** the bare word "persistence" (ambiguous with anchor below).

### 2. Anchor (or "structural persistence")

**Definition.** Unchanged across all four phases. A module intrinsic
to the patient and indifferent to the cognitive state.

| RPre | TT | RPost |
|---|---|---|
| X | X | X |

**Don't conflate with trace.** When the literature says "persistent
module" it usually means *anchor*. Use **anchor** when you mean "never
changed".

**Synonyms.** "anchor module", "structural anchor", "phase-invariant
module".

### 3. Reset

**Definition.** Reorganized in `task_test`, then reverted to the
`rest_pre` baseline in `rest_post`. The brain "resets" after the task.

| RPre | TT | RPost |
|---|---|---|
| X | reorganized to Y | back to X |

The triangle test gives `T_d > 0` here (RPost is closer to RPre than
to TT). This is the **opposite** of trace, not a weaker version of it.

**Synonyms.** "transient reorganization", "reverted module".

### 4. Emergent (or "fully reorganized")

**Definition.** A module that did not exist in `rest_pre`, came
together in `task_test` (and possibly `rest_post`). Tracks the
*emergence* of new community structure rather than the modification of
an existing one.

| RPre | TT | RPost |
|---|---|---|
| absent / fragmented | new coherent module | new module persists or dissolves |

This is a *membership* statement (which nodes form a module), not a
*value* statement (how strong are the edges of an existing module).
Edge-rank distances cannot test for emergent modules — that is an
LRG-dendrogram (community membership) question.

**Synonyms.** "emergent module", "novel community", "fully
reorganized".

---

## Where the raw-FC verdict sits

The 2026-04-29 raw-FC verdict (`final_verdict_table.csv`) is a
**trace** result on edge-rank ordering at the cohort level:

- `T_d^(d_S) < 0` ⟺ rank-ordering of edges in `RPost` is closer to
  the rank-ordering in `TT` than to the rank-ordering in `RPre` —
  task reshuffled the rank ordering AND the reshuffle persisted into
  `RPost`.

It is **not** an anchor result, **not** a reset result, **not** an
emergent-module result. The "n_persist" column in the cohort CSVs is
literally an `n_trace` count and should be read that way.

---

## Naming rules (apply everywhere)

In any new report, figure caption, script docstring, or memory entry:

- **`d(TT, RPost) < d(RPre, TT)` cohort-wide** ⇒ call it a **trace**
  or "task-trace persistence" — never the bare "persistence".
- **A module unchanged across all phases** ⇒ call it an **anchor** —
  never call it "persistent" without a qualifier.
- **`T_d > 0` cohort-wide** ⇒ call it a **reset** pattern.
- **A module that exists only after task onset** ⇒ **emergent** /
  "fully reorganized".

When the qualifier matters across bands (e.g. anchor in θ, trace in
β), prepend the band: "anchor in θ, trace in β".

### Variable naming in code

| use | not |
|---|---|
| `n_trace` | `n_persist` |
| `n_anchor` | "persistence count" |
| `n_reset` | (no prior name) |
| `n_emergent` | (LRG / community-membership only) |

Existing CSVs (`final_verdict_table.csv` and friends) still use
`n_persist` for backwards compatibility but the writing bundle and
all new code use `n_trace`. CSV columns will be renamed at the next
artefact rebuild.

### Figure / caption phrases

| use | not |
|---|---|
| "trace scatter on `d_S`" | "persistence scatter on `d_S`" |
| "trace zone (below identity)" | "persistence zone" |
| "task-trace persistence" | bare "persistence" |
| "the trace pulls `RPost` toward `TT`" | "`RPost` persists toward `TT`" |

---

## Why this matters

A reviewer or co-author reading "persistence in 4 of 6 bands" can
reasonably interpret it as **anchor**, when we mean **trace**. The
two are nearly opposite (one says "nothing changed", the other says
"task changed it AND the change stuck"). Loose terminology costs us
the reading every time.

A second cost: the LRG step will eventually classify modules with all
four labels (some are anchors that never change, some are emergent in
task and dissolve in `rest_post`, some are traces). Without the
taxonomy we cannot describe the LRG result without reinventing the
vocabulary.

---

## How to apply (checklist)

1. When proposing a new task-trace measure, name it explicitly with
   the taxonomy term: "trace measure on …" / "anchor detector for …".
2. When writing figure captions, replace "persistence scatter" with
   "trace scatter"; `n_persist` annotation with `n_trace`.
3. When in doubt, call out the ambiguity inline: "we use 'persistence'
   in the trace sense — see `terminology.md`".
4. When a result mixes behaviours across bands, list them with the
   right label per band: "anchor in θ, trace in β, reset in high_γ".
