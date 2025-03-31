#!/bin/sh

# Handles the Greengenes2 DADA2 database.

# Write the assignTaxonomy() fasta file: assignTaxonomy.fna
zcat *.gz | sed '/^>/s/>[^ ]\+ \([^[]\+\) \[.*/>\1/' | sed '/^>/s/ \[.*//' | sed 's/[a-z]__//g' | gzip -c > assignTaxonomy.fna.gz

# Write the addSpecies() fasta file: addSpecies.fna
# remove last ";", remove last ";" to preserve genus, take last field separated by ";", add line number after '>' as fake ID
zcat assignTaxonomy.fna.gz | sed 's/;*$//g' | sed 's/\(.*\);/\1 /' | sed 's/.*;/> /' | awk '$1 == ">" {print $1 FNR, $2, $3; next} {print}' | gzip -c > addSpecies.fna.gz
