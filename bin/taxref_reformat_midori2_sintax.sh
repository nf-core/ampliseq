#!/bin/sh

# Handles the MIDORI database for SINTAX.

set -eu

# The Zenodo zip holds the files for every gene, take the preformatted SINTAX file for CO1
unzip -p *.zip 'SINTAX/uniq/*_CO1_SINTAX.fasta.gz' > sintaxdb.fa.gz
test -s sintaxdb.fa.gz
