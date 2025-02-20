# Parsing Genbank records

File: [src/genbank_parser.py]()

The BioPython module `SeqIO` is used to parse genbank records. The following information is extracted directly from the sequence record:

- accession
- creation date
- authors
- sequence
- sequence length
- title of the record
- name of the organism

The taxid is a number that uniquely identifies a species in the NCBI Taxonomy database. It is obtained by searching for a `source` feature that has a `db_xref` qualifier. The taxid is rextracted from a `db_xref` qualified feature that starts with `taxon:`.

The position of the Cox2 gene within the sequence is extracted by searching for a feature of type `gene` where the first `gene`-qualified feature matches any of the Cox2 aliases. If no such feature exists, we search for a feature of type `CDS` where the `product`-qualifier matches any known Cox2 alias.

A record is discarded if the Cox2 position could not be found and the record title indicates that the record doesn't contain Cox2, the record is discarded.
