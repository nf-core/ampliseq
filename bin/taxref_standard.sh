#!/bin/sh

# Uses preformatted databses from DADA2 (https://benjjneb.github.io/dada2/training.html)
# The file for taxonomy assignment, identified by containing "train" in the name,
# and the file for add species, identified by containing "species" in the name, are renamed 

mv *train*gz assignTaxonomy.fna.gz
mv *species*gz addSpecies.fna.gz
