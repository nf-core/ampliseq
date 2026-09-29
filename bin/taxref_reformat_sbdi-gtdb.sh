#!/bin/sh

# Formatting script for sbdi-gtdb files.

# We get the files with numbers as names. Loop over them and find out which
# looks like assignTaxonomy and addSpecies respectively, gunzip to proper
# names.

for f in *; do
    if gunzip -c "$f" 2>/dev/null | head -n 1 | grep -q '>.*;.*;.*;.*;.*'; then
        # The reference repeats the domain as a second field, from SBDI's GBIF-compatible rank
        # scheme. Dropping it keeps the lineage as deep as taxlevels declares. Reverse-
        # complemented records carry a "Reversed:_" prefix on the first field only, so the
        # repeat has to be matched past it or those records keep a field the others lose and
        # every rank in them shifts. The back-references leave a release that stops repeating
        # the domain alone.
        gunzip -c "$f" \
            | sed -e '/^>/s/^>\([^;]*\);\1;/>\1;/' \
                  -e '/^>/s/^>Reversed:_\([^;]*\);\1;/>Reversed:_\1;/' \
            > assignTaxonomy.fna
    elif gunzip -c "$f" 2>/dev/null | head -n 1 | grep -q '>.* .*'; then
        gunzip -c "$f" > addSpecies.fna
    fi
done
