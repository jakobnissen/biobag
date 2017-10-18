# biobag

A collections of miscellaneous bioinformatics-related scripts

### Contents

__annotator.py__

A simple script to annotate metagenomic bins by Diamondblasting to a database. Uses a crude approach.

__deduplicatesam.py__

Parses a sorted SAM file and removes lines that correspond to the first pair of a read pair if both members map to the same gene.
This is useful when creating mapping statistics on shorter DNA pieces.

__misctools.py__

A collection of useful Python functions like a function for buffered printing to a file, iterating over a FASTA file,
creating directories which fail in a threadsafe way, and calculating statistics on assemblies.

__parsebins.py__

Parses a directory of FASTA files representing bins.
Can translate to amino acids, create a nonredundant FASTA file of all protein sequences in any bin,
and/or create an easily parse-able table of which contigs/genes belong in which bins.

__parsecanopy.py__

Parses Canopy clustering output, creating bin files and/or nonredundant FASTA file of all protein seqs in any bin.

__unifetch.py__

Script for searching a local or remote database for uniprot proteins given the primary citable accession number.
Also contains a function for creating a local database easily.
