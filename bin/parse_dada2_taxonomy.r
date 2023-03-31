#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if(length(args) != 1){
    stop("Usage: parse_dada2_taxonomy.r <ASV_tax_species.tsv>")
}

tax_file <- args[1]

OUT="tax.tsv"

# read required files
tax = read.table(tax_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = '', quote = '')

# Join columns 2:ncol(.) - 1, the taxonomy ranks (sequence is the last)
r <- colnames(tax)[!colnames(tax) %in% c('ASV_ID', 'sequence')]
tax$taxonomy <- do.call(paste, c(tax[r], sep = ';'))

#write
print (paste("write",OUT))
write.table(tax[,c('ASV_ID', 'taxonomy')], file = OUT, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
