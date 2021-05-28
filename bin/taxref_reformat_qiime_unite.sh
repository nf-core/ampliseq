#!/bin/sh

# Untar the Unite file
tar xzf *.gz

# Select and rename dynamic files
cat */*_dynamic_*.fasta > unite.fna
cat */*_dynamic_*.txt > unite.tax
