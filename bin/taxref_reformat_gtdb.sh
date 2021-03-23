#!/bin/sh

# Reads the ar122 and bac120 SSU fasta files from GTDB (after first untarring)
# and outputs two new fasta files, one suitable for DADA2's assignTaxonomy()
# and addSpecies() functions.

# Untar any tar file in the working directory
for f in *.tar.gz; do
  tar xzf $f
done

# Write the assignTaxonomy() fasta file: gtdb_assignTaxonomy.fna
sed '/^>/s/>\([^ ]\+\) \([^[]\+\) \[.*/>\2(\1\)/' ar122*.fna bac120*.fna | sed 's/[a-z]__//g' | sed 's/ /_/g' > assignTaxonomy.fna

# Write the assignTaxonomy() fasta file: addSpecies.fna
sed '/^>/s/>\([^ ]\+\) .*;s__\([^[]\+\) \[.*/>\1 \2/' ar122*.fna bac120*.fna > addSpecies.fna
