#!/bin/bash

# Handles the MIDORI database.

# There are preformatted DADA2 files for assignTaxonomy() and addSpecies()
gunzip -c *.fasta.gz | sed 's/\(phylum_\)\|\(class_\)\|\(order_\)\|\(family_\)\|\(genus_\)//g' | sed 's/_[0-9]\+;/;/g' > assignTaxonomy.fna

#unzip, remove last ";", take last field separated by ";" and add fake SeqID
gunzip -c *.fasta.gz | sed 's/;*$//g' | sed 's/.*;/> /' | sed 's/> \(.*\)_\([0-9]\+$\)/>\2 \1/' > addSpecies.fna

