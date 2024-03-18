#!/bin/sh

# Reads the ar* and bac* SSU fasta files from GTDB (after first untarring)
# and outputs two new fasta files, one suitable for DADA2's assignTaxonomy()
# and addSpecies() functions.

# Untar any tar file in the working directory
for f in *.tar.gz; do
    tar xzf $f
done

# Write the assignTaxonomy() fasta file: assignTaxonomy.fna
cat ar*.fna bac*.fna | sed '/^>/s/>[^ ]\+ \([^[]\+\) \[.*/>\1/' | sed '/^>/s/ \[.*//' | sed 's/[a-z]__//g' > assignTaxonomy.fna

# Write the addSpecies() fasta file: addSpecies.fna
cat ar*.fna bac*.fna | sed '/^>/s/>\([^ ]\+\) .*;s__\([^[]\+\) \[.*/>\1 \2/' > addSpecies.fna
