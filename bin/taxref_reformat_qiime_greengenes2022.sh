#!/bin/sh

# Decompress files.
gzip -c -d *.seqs.fna.gz > greengenes2.fna
gzip -c -d *.taxonomy.md5.tsv.gz > greengenes2.tax
