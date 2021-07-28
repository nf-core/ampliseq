#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 3){
    stop("Usage: combineTable.r <relASVTable.tsv> <repseq.fasta> <taxonomy.tsv>")
}

table <- args[1]
seq <- args[2]
taxonomy <- args[3]

OUT="qiime2_ASV_table.tsv"

#load packages, Biostrings_2.46.0
library("Biostrings")

#read required files
table <- read.table(file = table, sep = '\t', comment.char = "", skip=1, header=TRUE) #X.OTU.ID
tax <- read.table(file = taxonomy, sep = '\t', comment.char = "", header=TRUE) #Feature.ID
seq <- readDNAStringSet(seq)
seq <- data.frame(ID=names(seq), sequence=paste(seq))

#check if all ids match
if(!all(seq$ID %in% tax$Feature.ID))  {paste(seq,"and",taxonomy,"dont share all IDs, this is only ok when taxa were excluded.")}
if(!all(seq$ID %in% table$X.OTU.ID))  {stop(paste(seq,"and",table,"dont share all IDs, exit"), call.=FALSE)}

#merge
df <- merge(tax, seq, by.x="Feature.ID", by.y="ID", all.x=FALSE, all.y=TRUE)
df <- merge(df, table, by.x="Feature.ID", by.y="X.OTU.ID", all=TRUE)

#write
print (paste("write",OUT))
write.table(df, file = OUT, row.names=FALSE, sep="\t")
