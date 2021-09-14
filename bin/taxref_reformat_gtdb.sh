#!/bin/sh

# Reads the ar122 and bac120 SSU fasta files from GTDB (after first untarring)
# and outputs two new fasta files, one suitable for DADA2's assignTaxonomy()
# and addSpecies() functions.

# Untar any tar file in the working directory
for f in *.tar.gz; do
    tar xzf $f
done

# Write the assignTaxonomy() fasta file: assignTaxonomy.fna
cat ar122*.fna bac120*.fna | sed '/^>/s/>\([^ ]\+\) \([^[]\+\) \[.*/>\2(\1\)/' | sed '/^>/s/;s__.*//' | sed 's/[a-z]__//g' | sed 's/ /_/g' | sed '/^>/s/\(Archaea\)\|\(Bacteria\)/&;&/' > assignTaxonomy.fna

# Write the addSpecies() fasta file: addSpecies.fna
cat ar122*.fna bac120*.fna | sed '/^>/s/>\([^ ]\+\) .*;s__\([^[]\+\) \[.*/>\1 \2/' > addSpecies.fna
