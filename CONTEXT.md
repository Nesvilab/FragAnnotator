# Context

The domain language of FragAnnotator: the ion vocabulary it annotates spectra with, and
the MSFragger search settings it reads that vocabulary out of.

## Ion categories

Every annotated ion falls into exactly one of four categories. Each category is one group
of output columns in `psm.tsv`.

### Backbone ion

A peptide fragment produced by cleaving the backbone once, carrying **all modifications
intact**. The standard `a/b/c/x/y/z` series plus any custom series the search declared.
Always annotated, for every search.

Named by series letter, cleavage index and charge: `b5`, `y9++`. A custom series uses the
name the search gave it: `zOne5`.

### Fragment remainder ion

A backbone ion whose labile modification has been **partially or completely lost** during
fragmentation. Positional, like a backbone ion, but only defined for fragments that
**contain a labile modification site** — a fragment with no such site is just a backbone
ion, and must not be annotated twice.

Named by appending the signed remainder mass, rounded to the nearest integer: `b5+80`,
`y4+203`, `b4-18`. A remainder of `0` means the modification was lost entirely, leaving
the bare residue mass.

### Peptide remainder ion

The **intact peptide** carrying a partial modification loss. Not positional — there is one
per remainder mass, not one per cleavage site. Named `pep+203`, or by glycan composition
on a glyco search (see *Glycan Y ion*).

### Diagnostic ion

A low-m/z marker ion **characteristic of a modification**, detached from the peptide.
Not positional. Declared as *m/z*, not as a neutral mass.

## Modifications

### Labile modification

A modification that fragments in MS2 alongside the peptide backbone, and therefore
produces fragment remainder, peptide remainder and/or diagnostic ions. Declared by the
search as a **mass offset**, and only annotated as labile when the search ran in a labile
mode.

### Mass offset

One MSFragger `mass_offsets_detailed` entry: a mass, the **sites** it may occupy, and the
three ion lists it declares (`_d=` diagnostic, `_p=` peptide remainder, `_f=` fragment
remainder). Two offsets may share a mass and differ in sites and ion lists, so an offset
is identified by **mass and site together**, never mass alone.

### Site

The residues a mass offset may occupy. Either a plain residue list (`SKTYHDE`) or a
**sequon** — a bracketed motif (`{N[^P][ST]}`) describing the sequence context the offset
was searched in.

Only the **modified residue** of a sequon is a site. It is the residue named in
parentheses when the sequon has them, and the first residue otherwise. The rest of the
motif was the search's own constraint and is already satisfied by the search having
assigned the modification, so re-testing it can only reject what the search accepted.

## Searches

### Labile mode

Whether the search generated labile fragmentation at all. Decides whether the fragment
remainder, peptide remainder and diagnostic columns appear. A search that declares mass
offsets but ran with labile mode off is **not** a labile search — those offsets were never
searched.

### Glycopeptide PSM

A PSM whose labile modification is a glycan of known composition. Not a search-level
property and not a separate ion vocabulary: it annotates the same four categories, but
**names** its peptide remainder and diagnostic ions by glycan composition rather than by
mass. A PSM is a glycopeptide PSM exactly when its own glycan composition parses — the
search's parameters do not decide it.

### Glycan Y ion

A peptide remainder ion on a glycopeptide PSM, named by the composition it carries
(`Y_N2H3`) rather than by mass. Enumerated from that PSM's own glycan composition, so only
sub-compositions the identified glycan can actually produce are annotated.

## Remainder mass

The mass a labile modification leaves behind on an ion. `0` means the modification was lost
entirely; the full modification mass means it survived intact (which is a backbone ion, not
a remainder ion).

Where several labile modifications contribute to one ion, their remainders **sum**: the
ion's identity is its mass, so two combinations summing alike are one ion with one label.

## Localized and unlocalized

A modification the search placed on a specific residue is **localized**, and appears in the
PSM's assigned modifications. One it could not place is **unlocalized**, and survives only
as a delta mass.

Only a localized modification can produce fragment remainder ions, which need a position to
say which fragments contain it. Peptide remainder and diagnostic ions need no position and
are produced by both.

### Custom ion series

A user-defined backbone series (`ion_series_definitions`), declared as a name, a terminus
and a mass offset. Walks the whole backbone exactly like a standard series. Its offset is
stated against the **ordinary ion of its terminus** — `b` for N-terminal, `y` for
C-terminal.
