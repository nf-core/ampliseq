#!/bin/sh

# Reads the ar* and bac* SSU fasta files from GTDB (after first untarring/unzipping)
# and outputs two new fasta files, one suitable for DADA2's assignTaxonomy()
# and addSpecies() functions.

# Unzip any .fna.gz file in the working directory - versions 220 and newer
for f in *.fna.gz; do
    gunzip -c $f > $(basename "$f" .gz)
done

# Untar any tar file in the working directory - versions 214.1 and older
for f in *.tar.gz; do
    tar xzf $f
done

# Write the assignTaxonomy() fasta file: assignTaxonomy.fna
cat ar*.fna bac*.fna | sed '/^>/s/>[^ ]\+ \([^[]\+\) \[.*/>\1/' | sed '/^>/s/ \[.*//' | sed 's/[a-z]__//g' > assignTaxonomy.fna

# Write the addSpecies() fasta file: addSpecies.fna
cat ar*.fna bac*.fna | sed '/^>/s/>\([^ ]\+\) .*;s__\([^[]\+\) \[.*/>\1 \2/' > addSpecies.fna
