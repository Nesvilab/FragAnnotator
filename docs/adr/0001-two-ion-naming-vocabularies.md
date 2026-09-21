# 1. Two ion naming vocabularies, chosen per PSM

Date: 2026-09-21

## Status

Accepted

## Context

MSFragger declares labile fragmentation through `msfragger.mass_offsets_detailed`: per mass
offset, the diagnostic ions (`_d=`), peptide remainders (`_p=`) and fragment remainders
(`_f=`) it produces. Glycan fragmentation is declared through that *same* parameter — a
glyco search's offsets are its glycan compositions, its `_p=` list is that composition's Y
ladder, its `_d=` list is the oxonium markers, and `_f=203.07937` is the core-GlcNAc remnant.

So one generic annotator can cover both. The parameter names none of its ions, though, so a
glycan Y ion generated this way can only be labelled by its mass — `pep+1216` where the
previous glyco-specific code wrote `Y_N2H3`. Glycomics reads compositions, not masses, and
the mass label discards the very thing the analyst is looking for.

FragViz solved the overlap the other way: `workflow::is_glyco` detects a glyco search from
three separate parameters and excludes it from the generic path entirely, leaving glycans to
its own renderer. That heuristic exists because FragPipe's O-glyco workflows run in plain
`labile` mode — only N-glyco sets `labile_search_mode=nglycan` — so no single parameter
answers the question, and the three it reads have changed names across FragPipe versions.

We need both: generic coverage for ADPr and any future labile modification, and composition
labels for glycopeptides.

## Decision

Both vocabularies coexist, and **the PSM decides which one it gets**, not the search.

A PSM whose `Total Glycan Composition` parses is annotated with composition labels: peptide
remainders enumerated from its own glycan (`Y_N2H3`), diagnostic ions from the union of the
declared `_d=` list and the oxonium table. Every other PSM gets mass labels: `pep+203` from
the declared `_p=` list, `136.06` from the declared `_d=` list.

There is no search-level glyco flag.

## Consequences

The heuristic disappears. Nothing has to track how FragPipe spells its glyco switches,
because the question "is this a glycopeptide?" is answered by the PSM's own glycan
composition — written by the tool that actually assigned it.

The fallback stops being a special case. A glyco search whose PTM-Shepherd step did not run,
a result missing the glycan database, a PSM with no glycan assigned, and a file from a
mixed search all take the same path: no composition, so mass labels. We warn once per file
when the column is missing while the search declared `nglycan`, because that is the one case
where the analyst probably expected composition labels.

The cost is that a reader of one `psm.tsv` can meet two vocabularies in one column, and a
reader of the code will not find the `is_glyco` logic FragViz has. Hence this record.

The deeper cost is that composition labels depend on a *downstream* tool's output rather
than on the search parameters. If PTM-Shepherd's column is ever renamed, glyco results
silently degrade to mass labels rather than failing — annotated correctly, named less
usefully. The per-file warning is what makes that visible.
