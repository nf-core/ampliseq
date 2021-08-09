#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 3){
    stop("Usage: combineTable.r <relASVTable.tsv> <repseq.fasta> <taxonomy.tsv>")
}

table <- args[1]
seq <- args[2]
taxonomy <- args[3]

OUT="combined_ASV_table.tsv"

#load packages, Biostrings_2.46.0
library("Biostrings")

#read abundance file, rename first column to ID to have that independent of DADA2/QIIME2 input
table <- read.table(file = table, sep = '\t', comment.char = "", skip=1, header=TRUE) #X.OTU.ID
colnames(table)[1] <- "ID"

#read taxonomy file, rename first column to ID to have that independent of DADA2/QIIME2 input, also remove the additional sequence column
tax <- read.table(file = taxonomy, sep = '\t', comment.char = "", header=TRUE) #Feature.ID
colnames(tax)[1] <- "ID"
tax$sequence <- NULL

#read fasta file of ASV sequences
seq <- readDNAStringSet(seq)
seq <- data.frame(ID=names(seq), sequence=paste(seq))

#check if all ids match
if(!all(seq$ID %in% tax$ID))  {paste(seq,"and",taxonomy,"dont share all IDs, this is only ok when taxa were excluded.")}
if(!all(seq$ID %in% table$ID))  {stop(paste(seq,"and",table,"dont share all IDs, exit"), call.=FALSE)}

#merge
df <- merge(tax, seq, by.x="ID", by.y="ID", all.x=FALSE, all.y=TRUE)
df <- merge(df, table, by.x="ID", by.y="ID", all=TRUE)

#write
print (paste("write",OUT))
write.table(df, file = OUT, row.names=FALSE, sep="\t")
