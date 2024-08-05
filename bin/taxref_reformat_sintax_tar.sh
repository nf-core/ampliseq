#!/bin/sh

# Handles preformatted database tar files suitable for sintax
#
# This turned out to be a MISTAKE and is NOT USED, but I'm keeping the file for a while anyway.

# Extract the fasta file without _dev in its name
f=$(tar tfz *.tgz | grep fasta | grep -v '_dev')
tar xzf *.tgz $f

# Change the name and gzip
mv $f sintaxdb.fa
gzip sintaxdb.fa
