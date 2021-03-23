#!/bin/sh

# Untars a file ending with .gz downloaded from Unite, reformats a bit to 
# assignTaxonomy.fna and reformats a copy to addSpecies.fna.

# Untar the Unite file
tar xzf *.gz

# Remove leading "k__" and the like and replace space with underscore to create assignTaxonomy.fna
sed 's/[a-z]__//g' */*.fasta | sed 's/ /_/g' > assignTaxonomy.fna

# Create a copy 
sed 's/>\([^|]\+\)|\([^|]\+|[^|]\+\)|.*/>\2 \1/' assignTaxonomy.fna | sed 's/_/ /g' > addSpecies.fna
