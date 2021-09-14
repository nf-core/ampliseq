#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1){
    stop("Usage: metadata_all.r <metadata.tsv>")
}

metadata <- args[1]

data = read.delim(metadata)

#remove all numeric columns
nums <- unlist(lapply(data, is.numeric))
data <- data[ , !nums]

vector <- character()
for (i in 2:ncol(data)) {

    #remove blanks or NA
    cleandata <- data[!(is.na(data[i]) | data[i]==""),]

    #select only columns with multiple different values but not all unique
    if (nrow(unique(cleandata[i])) > 1 & nrow(unique(cleandata[i])) < nrow(cleandata[i])) {
        vector <- c(vector, colnames(cleandata[i]))
    }
}
vector <- paste(vector, collapse=",")
cat(vector)

