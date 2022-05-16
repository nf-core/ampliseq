#!/bin/bash

# Handles the MIDORI database.

# There are preformatted DADA2 files for assignTaxonomy() and addSpecies() -- this is just ungzipped
for f in MIDORI_UNIQ_NUC_*_DADA2.fasta.gz; do
    if [[ $f == *"MIDORI_UNIQ_NUC_SP_"* ]]; then
        gunzip -c $f > addSpecies.fna
    elif [[ $f == *"MIDORI_UNIQ_NUC_"* ]]; then
        gunzip -c $f > assignTaxonomy.fna
    fi
done
