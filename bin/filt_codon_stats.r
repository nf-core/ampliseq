#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 2){
    stop("Usage: filt_codons_stats.r <dada2_stats> <filtered_table>")
}

dada2_stats <- args[1]
coi_filt <- args[2]

old_stats <- read.table(dada2_stats, header = T, check.names = F, row.names = 1)

filt_table <- read.table(coi_filt, header = T, check.names = F, row.names = 1)

filt_counts <- as.data.frame(colSums(filt_table))
names(filt_counts) <- "codon_filtered"

new_stats <- merge(old_stats, filt_counts, by = 'row.names')

write.table(new_stats, file = "codon.filtered.stats.tsv", row.names=FALSE, sep="\t", col.names = TRUE, quote = FALSE, na = '')

