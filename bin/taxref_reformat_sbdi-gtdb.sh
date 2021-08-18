#!/bin/sh

# Formatting script for sbdi-gtdb files.

# We get the files with numbers as names. Loop over them and find out which
# looks like assignTaxonomy and addSpecies respectively, gunzip to proper
# names.

for f in *; do
    if [ $(file -L $f | grep -c 'gzip') ]; then
        if [ $(gunzip -c $f | head -n 1 | grep -c '>.*;.*;.*;.*;.*') -eq 1 ]; then
            gunzip -c $f > assignTaxonomy.fna
        elif [ $(gunzip -c $f | head -n 1 | grep -c '>.* .*') -eq 1 ]; then
            gunzip -c $f > addSpecies.fna
        fi
    fi
done
