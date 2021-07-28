#!/bin/sh

# Untars a file ending with .gz downloaded from Unite, reformats a bit to
# assignTaxonomy.fna and reformats a copy to addSpecies.fna.

# Untar the Unite file
tar xzf *gz

# Remove leading "k__" and the like, remove ranks classified as "unknown",
# and replace space with underscore to create assignTaxonomy.fna
cat */*.fasta | sed '/^>/s/;[ks]__.*//' | sed '/^>/s/[a-z]__unidentified//g' | sed '/^>/s/[a-z]__//g' | sed '/^>/s/ /_/g' | sed 's/>.*|/&Eukaryota;/' > assignTaxonomy.fna

# Reformat to addSpecies format
sed 's/>\([^|]\+\)|\([^|]\+|[^|]\+\)|.*/>\2 \1/' assignTaxonomy.fna | sed '/^>/s/_/ /g' > addSpecies.fna
