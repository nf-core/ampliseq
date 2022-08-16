#!/bin/bash

# Handles the MIDORI database.

# There are preformatted DADA2 files for assignTaxonomy() and addSpecies() -- this is just ungzipped
for f in MIDORI2_UNIQ_NUC_*_DADA2.fasta.gz; do
    if [[ $f == *"MIDORI2_UNIQ_NUC_SP_"* ]]; then
        #unzip, remove last ";", take last field separated by ";" and add fake SeqID
        gunzip -c $f | sed 's/;*$//g' | sed 's/.*;/>SeqID /' > addSpecies.fna
    elif [[ $f == *"MIDORI2_UNIQ_NUC_"* ]]; then
        gunzip -c $f > assignTaxonomy.fna
    fi
done