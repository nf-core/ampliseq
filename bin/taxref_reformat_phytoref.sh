#!/bin/sh

# Write the assignTaxonomy() fasta file: assignTaxonomy.fna
cat PhytoRef_with_taxonomy.fasta | sed '/>/s/>[^|]*|/>/' | sed '/>/s/|/;/g' > assignTaxonomy.fna

# Write the addSpecies() fasta file: addSpecies.fna
cat PhytoRef_with_taxonomy.fasta | sed '/^>/s/>\([^|]\+\)|.*|\([^|]\+\)/>\1 \2/' > addSpecies.fna
