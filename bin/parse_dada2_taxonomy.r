#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if(length(args) != 1){
    stop("Usage: parse_dada2_taxonomy.r <ASV_tax_species.tsv>")
}

tax_file <- args[1]

OUT="tax.tsv"

# read required files
tax = read.table(tax_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = '', quote = '')

# Join the taxonomy rank columns only -- excludes ASV_ID, sequence, the aggregate
# "confidence" column, and the per-rank "<rank>_confidence" columns (issue #1041),
# none of which are taxonomic ranks and shouldn't end up in the QIIME2 taxonomy string
non_rank_cols <- colnames(tax) %in% c('ASV_ID', 'sequence', 'confidence') | grepl('_confidence$', colnames(tax))
r <- colnames(tax)[!non_rank_cols]
tax$taxonomy <- do.call(paste, c(tax[r], sep = ';'))

#write
print (paste("write",OUT))
write.table(tax[,c('ASV_ID', 'taxonomy')], file = OUT, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
