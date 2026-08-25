#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if(length(args) != 1){
    stop("Usage: parse_dada2_taxonomy.r <ASV_tax_species.tsv>")
}

tax_file <- args[1]

OUT="tax.tsv"

# read required files
tax = read.table(tax_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = '', quote = '')

# Join the taxonomy rank columns only. This list has already missed a new non-rank column twice
# (the aggregate "confidence"/per-rank "<rank>_confidence" columns, then consolidate_dada2_taxonomy's
# "database" provenance column) -- matched by suffix pattern where possible ("_confidence", "_exact",
# covering DADA2_ADDSPECIES's "Species_exact") rather than one more literal name, so the next added
# metadata column doesn't silently leak into the taxonomy string as a bogus extra rank again.
non_rank_cols <- colnames(tax) %in% c('ASV_ID', 'sequence', 'confidence', 'database') |
    grepl('_confidence$|_exact$', colnames(tax))
r <- colnames(tax)[!non_rank_cols]
tax$taxonomy <- do.call(paste, c(tax[r], sep = ';'))

#write
print (paste("write",OUT))
write.table(tax[,c('ASV_ID', 'taxonomy')], file = OUT, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
