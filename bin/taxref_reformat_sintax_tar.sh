#!/bin/sh

# Handles preformatted database files suitable for sintax

# Just rename the preformatted file
# Assumes only one (gzipped) file

# Extract the fasta file without _dev in its name
tar xzf *.tgz $(tar tfz *.tgz | grep -v '_dev')

# Change the name and gzip
mv *.fasta sintaxdb.fa
gzip sintaxdb.fa
