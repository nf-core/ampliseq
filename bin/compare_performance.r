#!/usr/bin/env Rscript

# compare_sequences.r

# Get params and files from the command line
args            <- commandArgs(trailingOnly=TRUE)
tag             <- args[1] # tag for each sample
alignmentFILE   <- args[2] # alignment file; expected columns: query, target, mismatches_final
obsabundFILE    <- args[3] # observed abundance*
expabundFILE    <- args[4] # expected abundance*
fbeta           <- as.numeric(args[5]) # Fbeta weight, default 2; 1=precision and recall equally weighted (=F1 score), 2=weighs recall higher than precision, 0.5=weighs recall lower than precision.
mismatch_threshold <- as.numeric(args[6])
# tab-separated file with header: first column with sequences name, following one or many columns (=samples) with numeric values (=abundance), only presence (>0)/absence are used here

# function to produce statistics
get_stats <- function(i_exp,i_obs,df,sample) {

	# if none are expected, skip!
	if( length(i_exp)==0 ) {
		print( paste("No expected seq in sample",sample, "- skipping") ); next
	} else {
		print( paste(length(i_exp),"expected seq in sample",sample) )
	}

	# stats
	TP <- intersect(i_exp, i_obs)
	FN <- setdiff(i_exp, i_obs)
	FP <- setdiff(i_obs, i_exp)
	if ( length(TP) > 0 ) {
		Fone <- ( 2 * length(TP)/length(i_obs) * length(TP)/length(i_exp) ) / ( length(TP)/length(i_obs) + length(TP)/length(i_exp) )
		Fbeta <- ( (1+fbeta^2) * length(TP)/length(i_obs) * length(TP)/length(i_exp) ) / ( (fbeta^2) * length(TP)/length(i_obs) + length(TP)/length(i_exp) )
	} else {
		Fone <- 0
		Fbeta <- 0
	}

	# save
	types <- c(
		"observed",
		"expected",
		"TP",
		"FN",
		"FP",
		"recall",
		"precision",
		"F1",
		"Fbeta",
		"fdr",
		"jaccard",
		"TPs_exp",
		"FNs_exp",
		"FPs_obs"
    )
	values <- c(
		length(i_obs),
		length(i_exp),
		length(TP),
		length(FN),
		length(FP),
		length(TP)/length(i_exp),
		length(TP)/length(i_obs),
		Fone,
		Fbeta,
		length(FP)/(length(i_obs)),
		length(TP)/(length(i_obs)+length(FN)),
		paste(head(TP, n=100), collapse=','),
		paste(head(FN, n=100), collapse=','),
		paste(head(FP, n=100), collapse=',')
	)
	ids <- rep(sample, length(types))

	df_append <- data.frame(
		sample=ids,
		type= types,
		value= values,
		stringsAsFactors=FALSE)
	df <- rbind( df, df_append )

	return(df)
}

# PREPARE INPUT

# Read sequence alignment table
alignment = read.table( alignmentFILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE)
alignment <- alignment[alignment$target !="*" & alignment$mismatch_final <= mismatch_threshold,]

# Read pipeline's ASV abundance table (either from QIIME2 [skipping first line] or from previous steps)
if (grepl("^# Constructed from biom file", readLines(obsabundFILE, n = 1))) {
	observed = read.table( obsabundFILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE, comment.char = "", skip=1 )
} else {
	observed = read.table( obsabundFILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE, comment.char = "" )
}
colnames(observed)[1] <- "ID"

# select available samples
observed_samples <- colnames(observed)[2:ncol(observed)]
print(paste( "Observed samples:", paste(observed_samples,collapse=",")))

# Read expected abundance table
exp = read.table( expabundFILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE)
colnames(exp)[1] <- "ID"
# extract samples to analyse
exp_samples <- colnames(exp)[2:ncol(exp)]
print(paste( "Expected samples:", paste( exp_samples ,collapse=",")))
SAMPLES <- exp_samples[exp_samples %in% observed_samples]
print(paste( "Investigate samples:", paste( SAMPLES ,collapse=",")))

# check if there are any samples to analyse
if (length(SAMPLES) == 0) {
	stop("ERROR - Found no samples to investigate")
}

# ANALYSE

# Initialize dataframe
df <- data.frame(
	sample=character(),
	type=character(),
	value=character(),
	stringsAsFactors=FALSE)

# for each sample
for (sample in SAMPLES) {

	# keep only non-empty columns and sample abundances
	keep_cols <- c("ID",sample)
	s_exp <- subset(exp, select = keep_cols)
	s_obs <- subset(observed, select = keep_cols)

	# keep only observed (abundance > 0)
	s_exp = s_exp[s_exp[,2] > 0,]
	s_obs = s_obs[s_obs[,2] > 0,]

	# keep only alignments where both, query and target are present
	s_alignment <- alignment[alignment$query %in% s_obs$ID & alignment$target %in% s_exp$ID,]

	# translate s_obs$ID to alignment$target if available
	obsIDs <- c()
	for (obsID in s_obs$ID){
		i_alignment <- s_alignment[s_alignment$query == obsID,]
		# collapse several targets if multiple matches
		if( nrow(i_alignment) > 0 ){
			target <- paste(i_alignment$target, collapse="-")
		} else {
			target <- obsID
		}
		obsIDs <- c(obsIDs, target)
	}

	# aggregate s_exp$ID if several alignment$query match
	expIDs <- s_exp$ID[!s_exp$ID %in% s_alignment$target]
	i_alignment <- aggregate(target ~ query, s_alignment, paste, collapse = "-")
	expIDs <- c(expIDs, i_alignment$target)

	# calculate stats
	df <- get_stats( unique(expIDs), unique(obsIDs), df, sample)
}

# Write detailed output
outfile <- "performance_per-sample.tsv"
df$tag <- rep(tag, nrow(df))
print(paste("write",outfile))
write.table(df, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep="\t")


# Make & write summary output
outfile <- "performance_summary.tsv"
df_sum <- data.frame(type=character(), median=numeric(), mean=numeric(), min=numeric(), max=numeric(), count=numeric(), stringsAsFactors=FALSE)
print("Warnings are expected at that point:")
for (type in unique(df$type)) {
	df_subset <- df[df$type == type, ]
	MEDIAN <- median( as.numeric(df_subset$value) )
	MEAN <- mean( as.numeric(df_subset$value) )
	MIN <- min( as.numeric(df_subset$value) )
	MAX <- max( as.numeric(df_subset$value) )
	COUNT <- length( as.numeric(df_subset$value) )
	df_sum[nrow(df_sum) + 1,] <- list(type, MEDIAN, MEAN, MIN, MAX, COUNT)
}
df_sum$tag <- rep(tag, nrow(df_sum))
print(paste("write",outfile))
write.table(df_sum, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep="\t")

# Plot Values
df_subset <- subset(df, type %in% c("recall","precision","F1","Fbeta","fdr","jaccard") )
df_subset$value <- as.numeric(df_subset$value)

outfile <- "performance_boxplot.svg"
print(paste("write",outfile))
svg(outfile, height = 8, width = 10)
boxplot(value~type, data=df_subset, xlab="Type", ylab="Value", ylim = c(0, 1))
invisible(dev.off())

outfile <- "performance_boxplot.png"
print(paste("write",outfile))
png(outfile, height = 400, width = 500)
boxplot(value~type, data=df_subset, xlab="Type", ylab="Value", ylim = c(0, 1))
invisible(dev.off())
