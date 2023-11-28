#!/bin/sh

# Decompress files.
gzip -d *.gz

# Select and rename files
mv *.fna greengenes2022.fna
mv *.tsv greengenes2022.tax
