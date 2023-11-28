#!/bin/sh

# Decompress files.
gzip -c -d 2022.10.seqs.fna.gz > 2022.10.seqs.fna
gzip -c -d 2022.10.taxonomy.md5.tsv.gz > 2022.10.taxonomy.md5.tsv

# Select and rename files
mv *.fna greengenes2022.fna
mv *.tsv greengenes2022.tax
