#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1){
    stop("Usage: metadataCategory.r <metadata.tsv>")
}

metadata <- args[1]

data = read.delim(metadata)

#remove all numeric columns
nums <- unlist(lapply(data, is.numeric))
data <- data[ , !nums]

vector <- character()
for (i in 1:ncol(data)) {

    #remove blanks or NA
    cleandata <- data[!(is.na(data[i]) | data[i]==""), ]

    #select only columns that have at least 2 of each value so that it can be used for pairwise comparisons
    noccur <- data.frame(table(cleandata[i]))
    if (nrow(unique(cleandata[i])) > 1 & nrow(unique(cleandata[i])) < nrow(cleandata[i])) {
        if ( nrow(noccur[noccur$Freq != 1,]) == nrow(noccur) ) {
            vector <- c(vector, colnames(cleandata[i]))
        }
    }
}
vector <- paste(vector, collapse=",")
cat(vector)
