#!/bin/sh

# Unzip the qza files
unzip *seqs.qza
unzip *tax.qza

# Select and rename dynamic files
cat */data/taxonomy.tsv > silva.tax
cat */data/*sequences.fasta > silva.fna
