#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if(length(args) != 2){
    stop("Usage: dada_quality.r <str; output prefix> <int; The number of records to sample from the fastq file.>")
}

OUT <- args[1]
number_of_records <- as.integer(args[2])

print(OUT)
print(number_of_records)

suppressPackageStartupMessages(library(dada2))

readfiles <- sort(list.files(".", pattern = ".fastq.gz", full.names = TRUE))
plot <- plotQualityProfile(readfiles, n = number_of_records, aggregate = TRUE)
data <- plot$data

df <- data.frame(Cycle=character(), Count=character(), Median=character(), stringsAsFactors=FALSE)
cycles <- sort(unique(data$Cycle))

#aggregate data for each sequencing cycle
for (cycle in cycles) {
    subdata <- data[data[, "Cycle"] == cycle, ]
    score <- list()
    #convert to list to calculate median
    for (j in 1:nrow(subdata)) {score <- unlist(c(score, rep(subdata$Score[j], subdata$Count[j])))}
    temp = data.frame(Cycle=cycle, Count=sum(subdata$Count), Median=median(score), stringsAsFactors=FALSE)
    df <- rbind(df, temp)
}

#write output
write.table( t(df), file = paste0(OUT,".tsv"), sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
pdf(paste0(OUT,".pdf"))
plot
dev.off()
