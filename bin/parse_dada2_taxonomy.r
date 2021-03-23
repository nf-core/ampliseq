#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if(length(args) != 2){
    stop("Usage: parse_dada2_taxonomy.r <DADA2_table.tsv> <ASV_tax_species.tsv>")
}

asv_file <- args[1]
tax_file <- args[2]

OUT="tax.tsv"

#load packages, Biostrings_2.46.0
library("Biostrings")

#read required files
tax = read.table(tax_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
asv = read.table(asv_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#merge
asvtax = merge(asv, tax, by = "sequence")

#remove sample counts
sample_columns = colnames(tax)[2:length(colnames(tax))]
df = asvtax[c("ASV_ID", sample_columns)]

#merge taxonomy & remove other columns
df$taxonomy<-do.call(paste, c(df[sample_columns], sep=";"))
for(column in sample_columns) df[column]<-NULL

#write
print (paste("write",OUT))
write.table(df, file = OUT, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")