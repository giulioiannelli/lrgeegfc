---
name: overleaf-mirror-readme
era: IMCOH_ABS × COHORT_N10
status: current
kind: protocol
scope: the single copy-paste mirror of the manuscript sections as they live in the external Overleaf project
---

# Overleaf mirror — the copy-paste source of truth

**Head.** This folder holds the manuscript **sections** exactly as they live in
the external Overleaf project. These are the `.tex` files we edit together and
that you copy-paste to/from Overleaf. Nothing here is auto-synced (no Overleaf
API); the sync is the manual copy-paste protocol below. Everything that is *not*
a manuscript section — directives, draft paragraphs, notes — lives outside this
folder so the mirror stays a clean 1:1 image of the paper.

## Files (each mirrors one Overleaf file, same basename)

| File | Overleaf section |
|------|------------------|
| `methods.tex` | Materials and Methods |
| `results_sec_1.tex` | Results (section 1; figures embedded inline) |
| `results_sec_2.tex` | Results (section 2) |

Add `intro.tex`, `discussion.tex`, … here as they materialize. Figures are
embedded directly in each section `.tex` (`\includegraphics`), not `\input` floats.

## Copy-paste protocol

1. **We edit here → you push to Overleaf.** Copy the file's contents over the
   matching Overleaf file, compile, done.
2. **You edit on Overleaf → you pull back here.** Paste the changed section back
   into the matching file here, so agents always see the *current* manuscript
   state (never let the mirror drift).
3. **Preamble stays in Overleaf.** These files assume the macros / acronyms are
   defined in the Overleaf preamble (`math_commands.tex`, `acronyms.tex`,
   `\giampi`) — e.g. `\lapl \propag \Dcoph \rhocoph \txtacr \gls \num \qty`,
   `\acrshort{…}`. A file here is a *fragment* (`\section{…}` onward), not a
   standalone document.
4. **`\input` paths inside these files mirror Overleaf, not this repo** — do not
   rewrite them to repo paths.

## What is deliberately NOT here

- Find/replace + writing **directives** → `../directives/`
- **Per-result draft paragraphs** (`R1.1 … R3.4.tex`) that compose into
  `results_sec_1.tex` → `../directives/results_paragraphs/`
- **Methods companion notes** (Grassmann methodology, neurophysiological
  interpretation, displacement taxonomy) → `../directives/`
- **Superseded methods reviews** → `../directives/archive/2026-05/`
- Band briefs, locked ledgers, responses, table drafts + notes → their own folders.
