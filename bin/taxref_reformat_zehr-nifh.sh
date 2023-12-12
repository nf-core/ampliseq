#!/bin/sh

# Write the assignTaxonomy() fasta file: assignTaxonomy.fna
cp *.fasta assignTaxonomy.fna

# Write the addSpecies() fasta file: addSpecies.fna
cut -d, -f 2,6,7 *.csv  | grep -v '^sequence,' | sed 's/\(.*\),[0-9]* \(.*\),\(.*\)/>\3 \2\n\1/' > addSpecies.fna
