#!/bin/sh

# Diat.barcode (rbcL) preformatted for DADA2's assignTaxonomy().
# The staged file is named by its Dataverse id, so the single file is used whatever it is called.
# Taxonomy strings without the final ";" are completed: DADA2 would treat them as a separate group from the same taxon with it.

# Write the assignTaxonomy() fasta file: assignTaxonomy.fna
cat * | tr -d '\r' | awk 'NF' | sed '/^>/{/;$/!s/$/;/}' > assignTaxonomy.fna

# Write the addSpecies() fasta file: addSpecies.fna
# Genus and species are the two first words of the last taxonomy level; "sp.", "cf." and "aff." are not species epithets.
awk -F';' '
/^>/ { keep = 0; split($(NF-1), b, " "); if (b[2] != "" && b[2] !~ /^(sp\.|cf\.|aff\.)$/) { keep = 1; n++; print ">ID_" n " " b[1] " " b[2] } ; next }
keep { print }
' assignTaxonomy.fna > addSpecies.fna
