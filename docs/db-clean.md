# Database cleanup

This document describes how the metadata inside an input csv file is used to
remove unwanted sequences from the database. The objective of this step is to
improve the taxonomic resolution of the database and to increase taxonomic
recall of query sequences.

## Environmental sequences

Enironmental sequences were obtained by metabarcoding using amplicon sequencing.
The taxonomic origin of these sequence was usually determined by BLASTing the
respecitve sequence against the NCBI nucleotide database, which often results in
inorrect taxonomic assignment.

Environmental sequences are recognized by searching for the keywords
`uncultured` and `environmental` in the `organism_name` column and by searching
for `unverified` in the `title` column.

## Imprecise species annotation and invalid species

Some sequences are only annotated as `[Genus] sp.`, thus missing the information
to which species they belong. I think many of these sequences were also uploaded
as a result of metabarcoding studies and are thus unreliable for taxonomic
classification. Such sequences are removed by searching for `sp.` in the column `organism_name`.

Additionally, sequences with invalid species names are removed by searching the
`title` and `organism_name` columns for `nom. inval.`. See
[this Wikipedia page](https://en.wikipedia.org/wiki/Validly_published_name) for
more information.

## Missing COX2 gene

Some of the sequences don't contain the COX2 gene. Usually, these sequences are
in the database because the GenBank record contains the name COX2 somewhere.

Sequences not containing the COX2 gene are recognized because the values in the
columns `cox2_start` and/or `cox2_end` are missing.

## Too short sequences

Sequences that are considerably shorter than the targeted amplicon are removed
because they are missing information that allows for the distingtion between
closely related species. Keeping those sequences would result in a lower
resolution of the database.

## Prefer RefSeq sequences

RefSeq sequences are [curated by NCBI](https://www.ncbi.nlm.nih.gov/refseq/about/)
and therefore have a higher quality than other sequences. RefSeq sequences can
be distinguished by an `_` in their accession.
Therefore, if a species has at least one RefSeq sequence, all non-RefSeq
sequences are discarded.
This step has to be done before any clustering to avoid accidently removing
RefSeq sequences.

This step leads to a substantial higher resolution for the _Phythophthora_ genus.

## Fix organism name and taxid

Sometimes, the organism name is cluttered with additional information like
subspecies names, strain identifiers, or 'filler' terms like `var.` or `aff.`.
This results in errors during taxonomic classification. To fix this issue, only
the first two words of the species name, excluding any 'filler' terms are retained.
The fixed species names and tax_ids are saved to the columns
`organism_name_sanitized` and `tax_id_sanitized`.

## Trim sequences to COX2 region

Many sequences (mitochondiral genomes in particular) contain additional
sequences outside of the COX2 region. Therefore, any additional parts outside of
COX2 need to be removed before clustering/ dereplication. This is done using the
`cox2_start`/ `cox2_end` columns that specify the start/ end of the COX2 gene
relative to the entire sequence as annotated in the NCBI database.

> [!IMPORTANT]
> We found that the annotation in the database is not always correct:
> For some sequences, COX2 was annotated to start at the first 30-50 bp,
although in many cases the sequence of the forward primer was found before the
annotated COX2 start.

## Reverse complement

The `cox2_strand` column indicates weither COX2 is on the positive strand (1) or
on the reverse-complement strand (-1). Sequences on the reverse complement strand
are put in their correct orientation by reverse complementing them. This is
neccessarry to avoid issues with sequences on different strands for downstream applications.
