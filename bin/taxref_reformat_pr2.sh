#!/bin/sh

# Handles the PR2 database.

# There's a preformatted DADA2 file for assignTaxonomy() -- this is just ungzipped
gunzip -c *dada2.fasta.gz > assignTaxonomy.fna

# For addSpecies(), the UTAX file is downloaded and reformated to only contain the id and species.
# The second two sed calls are to replace "_" with space only in the species name and not the last part of the id (overdoing it a bit, as I don't the id actually matters as long as it's unique).
gunzip -c *UTAX.fasta.gz | sed '/^>/s/>\([^;]*\);.*,s:\(.*\)/>\1 \2/' | sed 's/_/ /g' | sed 's/ \([A-Z]\) /_\1 /' > addSpecies.fna
