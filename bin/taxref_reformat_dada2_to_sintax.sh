#!/bin/sh

# Converts the DADA2 assignTaxonomy() fasta (assignTaxonomy.fna, produced by a taxref_reformat_*.sh script)
# to a gzipped fasta for VSEARCH sintax: sintaxdb.fa.gz
#
# ">Bacteria;Proteobacteria;Gammaproteobacteria;" becomes ">s1;tax=d:Bacteria,k:Proteobacteria,p:Gammaproteobacteria"
# Rank letters are assigned by position. VSEARCH accepts the nine letters d,k,p,c,o,f,g,s,t, so the letter
# does not match the rank name, and the rank names are given by --sintax_assign_taxlevels / taxlevels.
# A taxonomy string stops at its first empty rank, because sintax output is read by position.

awk '
BEGIN { split("d k p c o f g s t", letter, " ") }
/^>/ {
    n = split(substr($0, 2), rank, ";")
    tax = ""
    for (i = 1; i <= n && rank[i] != ""; i++) {
        if (i > 9) { print "ERROR: more than 9 ranks in " $0 > "/dev/stderr"; exit 1 }
        gsub(/[,:]/, "_", rank[i])
        tax = tax (i > 1 ? "," : "") letter[i] ":" rank[i]
    }
    print ">s" ++id ";tax=" tax
    next
}
{ print }
' assignTaxonomy.fna | gzip -c > sintaxdb.fa.gz
