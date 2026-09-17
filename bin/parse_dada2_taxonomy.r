#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if(length(args) != 1){
    stop("Usage: parse_dada2_taxonomy.r <ASV_tax_species.tsv>")
}

tax_file <- args[1]

OUT="tax.tsv"

# read required files
tax = read.table(tax_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = '', quote = '')

# Join the taxonomy rank columns only, excluding known non-rank columns. Matched by suffix
# pattern where possible ("_confidence", "_exact") rather than by literal name, so a future
# added metadata column doesn't silently leak into the taxonomy string as a bogus extra rank.
non_rank_cols <- colnames(tax) %in% c('ASV_ID', 'sequence', 'confidence', 'database') |
    grepl('_confidence$|_exact$', colnames(tax))
r <- colnames(tax)[!non_rank_cols]
tax$taxonomy <- do.call(paste, c(tax[r], sep = ';'))

#write
print (paste("write",OUT))
write.table(tax[,c('ASV_ID', 'taxonomy')], file = OUT, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
