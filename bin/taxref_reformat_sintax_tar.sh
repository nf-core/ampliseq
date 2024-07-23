#!/bin/sh

# Handles preformatted database files suitable for sintax

# Just rename the preformatted file
# Assumes only one (gzipped) file

# Extract the fasta file without _dev in its name
f=$(tar tfz *.tgz | grep fasta | grep -v '_dev')
tar xzf *.tgz $f

# Change the name and gzip
mv $f sintaxdb.fa
gzip sintaxdb.fa
