#!/bin/sh

# Uses preformatted databases from DADA2 (https://benjjneb.github.io/dada2/training.html)
# The file for taxonomy assignment, identified by containing "train" in the name,
# gets the first field duplicated:
gunzip -c *train*gz | sed 's/>\([^;]*\)/>\1;\1/' > assignTaxonomy.fna

# and the file for add species, identified by containing "species" in the name, is renamed
mv *species*gz addSpecies.fna.gz
