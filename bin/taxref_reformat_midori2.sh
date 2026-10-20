#!/bin/bash

# Handles the MIDORI database.

set -euo pipefail

# The Zenodo zip holds the files for every gene, take the preformatted DADA2 file for CO1
unzip -p *.zip 'DADA2_sp/uniq/*_CO1_DADA2.fasta.gz' > midori2_co1.fasta.gz
test -s midori2_co1.fasta.gz

# There are preformatted DADA2 files for assignTaxonomy() and addSpecies()
gunzip -c midori2_co1.fasta.gz | sed 's/\(phylum_\)\|\(class_\)\|\(order_\)\|\(family_\)\|\(genus_\)//g' | sed 's/_[0-9]\+;/;/g' > assignTaxonomy.fna

#unzip, remove last ";", take last field separated by ";" and add fake SeqID
gunzip -c midori2_co1.fasta.gz | sed 's/;*$//g' | sed 's/.*;/> /' | sed 's/> \(.*\)_\([0-9]\+$\)/>\2 \1/' > addSpecies.fna

rm midori2_co1.fasta.gz
