#!/bin/bash

# Handles the MIDORI database.

# There are preformatted DADA2 files for assignTaxonomy() and addSpecies() -- this is just ungzipped
for f in MIDORI2_UNIQ_NUC_*_DADA2.fasta.gz; do
    if [[ $f == *"MIDORI2_UNIQ_NUC_SP_"* ]]; then
        gunzip -c $f > addSpecies.fna
    elif [[ $f == *"MIDORI2_UNIQ_NUC_"* ]]; then
        # gets the first field duplicated:
        gunzip -c $f | sed 's/>\([^;]*\)/>\1;\1/' > assignTaxonomy.fna
    fi
done
