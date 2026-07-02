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

# function to produce abundance statistics
get_stats_abundance <- function(expected_abund,observed_abund,df,sample) {

	# if none are expected, skip!
	if( length(expected_abund)==0 ) {
		print( paste("No expected abundances in sample",sample, "- skipping") ); next
	} else {
		print( paste(length(expected_abund),"expected seq in sample",sample) )
	}

	# stats
	pearson_cor <- cor.test(expected_abund, observed_abund, method = "pearson")
	spearman_rho <- cor.test(expected_abund, observed_abund, method = "spearman")
	percent_dev <- abs((observed_abund - expected_abund) / expected_abund) * 100
	jensen_shannon <- function(p, q) {
		# KL divergence function (with pseudocount to handle zeros)
		kl_div <- function(x, y, pseudocount = 1e-10) {
			x <- x + pseudocount
			y <- y + pseudocount
			x <- x / sum(x)  # normalize
			y <- y / sum(y)  # normalize
			sum(x * log(x / y))
		}

		# Midpoint distribution
		m <- (p + q) / 2

		# Jensen-Shannon is the square root of the average of two KL divergences
		sqrt(0.5 * kl_div(p, m) + 0.5 * kl_div(q, m))
	}

	# save
	types <- c(
		"pearson_cor",
		"spearman_rho",
		"deviation",
		"mae",
		"rmse",
		"ps",
		"bray-curtis",
		"hellinger",
		"jensen-shannon"
	)
	values <- c(
		as.list(pearson_cor$estimate)[[1]], # Pearson's Correlation (cor)
		as.list(spearman_rho$estimate)[[1]], # Spearman's Rank Correlation (rho)
		median(percent_dev, na.rm = TRUE), # Median Percent Abundance Deviation
		mean(abs(observed_abund - expected_abund), na.rm = TRUE), # Mean Absolute Error (MAE)
		sqrt(mean((observed_abund - expected_abund)^2, na.rm = TRUE)), # Root Mean Square Error (RMSE)
		sum(pmin(observed_abund, expected_abund)) / sum(expected_abund), # Proportional Similarity (PS)
		1 - (2 * sum(pmin(observed_abund, expected_abund))) / (sum(observed_abund) + sum(expected_abund)), # Bray-Curtis Dissimilarity
		sqrt(sum((sqrt(observed_abund) - sqrt(expected_abund))^2)), # Hellinger Distance
		jensen_shannon(observed_abund, expected_abund) # Jensen-Shannon Divergence
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

# (A) PREPARE INPUT

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

# (B) ANALYSE

# Initialize dataframes
df <- data.frame(
	sample=character(),
	type=character(),
	value=character(),
	stringsAsFactors=FALSE)
abundances <- data.frame(
	sample=character(),
	ID=character(),
	query=character(),
	target=character(),
	observed_abund=numeric(),
	expected_abund=numeric(),
	stringsAsFactors=FALSE)

# for each sample
for (sample in SAMPLES) {

	# (1) FILTER - for presence/absence in that sample
	# keep only non-empty columns and sample abundances
	keep_cols <- c("ID",sample)
	s_exp <- subset(exp, select = keep_cols)
	s_obs <- subset(observed, select = keep_cols)
	# keep only observed (abundance > 0)
	s_exp = s_exp[s_exp[,2] > 0,]
	s_obs = s_obs[s_obs[,2] > 0,]
	# keep only alignments where both, query and target are present
	s_alignment <- alignment[alignment$query %in% s_obs$ID & alignment$target %in% s_exp$ID,]

	# (2) AGGREGATE - combine IDs and abundance if multiple query or targets match each other
	colnames(s_exp) <- c("ID","expected_abund")
	colnames(s_obs) <- c("ID","observed_abund")
	s_alignment <- subset(s_alignment, select = c("query","target"))
	merged <- merge(s_alignment, s_exp, by.x="target", by.y="ID", all=TRUE)
	merged <- merge(merged, s_obs, by.x="query", by.y="ID", all=TRUE)
	merged <- unique(merged)
	# retain exp without query & obs without target
	nomatch_exp <- merged[is.na(merged$query),]
	nomatch_obs <- merged[is.na(merged$target),]
	# aggregate query target by query (exp without query is lost)
	merged <- merged[order(merged$target), ]
	merged_exp <- aggregate(target ~ query, merged, paste, collapse = "-")
	merged_exp_abund <- aggregate(expected_abund ~ query, merged, sum)
	data_exp <- merge(merged_exp, merged_exp_abund, by.x="query", by.y="query", all=TRUE)
	data_exp <- unique( merge(data_exp, s_obs, by.x="query", by.y="ID", all=TRUE) )
	# aggregate query by target (obs without target is lost)
	data_exp <- data_exp[order(data_exp$query), ]
	merged_obs <- aggregate(query ~ target, data_exp, paste, collapse = "-")
	merged_obs_abund <- aggregate(observed_abund ~ target, data_exp, sum)
	data <- merge(merged_obs, merged_obs_abund, by="target", all=TRUE)
	# add exp abundance
	data_exp_abund <- subset(data_exp, select=c("target","expected_abund"))
	data_exp_abund <- unique(data_exp_abund[!is.na(data_exp_abund$target),])
	data <- merge(data, data_exp_abund, by="target", all=TRUE)
	# re-add exp without query & obs without target
	data <- rbind(data, nomatch_obs)
	data <- rbind(data, nomatch_exp)
	# add ID list
	ID_list <- ifelse(is.na(data$target), data$query, data$target)
	data$ID <- factor(ID_list, levels = ID_list)
	data$expected_abund[is.na(data$expected_abund)] <- 0
	data$observed_abund[is.na(data$observed_abund)] <- 0
	# make relative abundances
	data$expected_abund <- data$expected_abund/sum(data$expected_abund,na.rm=TRUE)
	data$observed_abund <- data$observed_abund/sum(data$observed_abund,na.rm=TRUE)
	# Sort by expected abundance
	data <- data[order(data$expected_abund, decreasing = TRUE), ]
	# append sample data to table
	data$sample <- rep(sample, nrow(data))
	data <- data[, c("sample", "ID", "query", "target", "observed_abund", "expected_abund")]
	abundances <- rbind( abundances, data )

	# (3) STATISTICS

	# calculate absence/presence stats
	expIDs <- data$target[!is.na(data$target)]
	obsIDs <- data$ID[!is.na(data$query)]
	df <- get_stats( expIDs, obsIDs, df, sample)

	# calculate abundance stats, this includes also non-matching sequences
	df <- get_stats_abundance( data$expected_abund, data$observed_abund, df, sample)

	# (4) PLOTS - pairwise comparisons

	# Side-by-Side Bar Plots
	outfile <- paste0(sample,"_abundance_barplot.svg")
	print(paste("write",outfile))
	svg(outfile, height = 8, width = 10)
	barplot(
		rbind(data$expected_abund, data$observed_abund),
		beside = TRUE,
		names.arg = data$ID,
		col = c("skyblue", "salmon"),
		legend.text = c("Expected", "Observed"),
		ylab = "Abundance",
		xlab = "Taxon",
		main = "Side-by-Side Bar Plot: Expected vs. Observed Abundance",
		cex.names = 0.7,
		las = 2
	)
	invisible(dev.off())

	# Scatter plot: Observed vs. Expected Abundance (log-log)
	data_matches <- data[data$observed_abund >0 & data$expected_abund >0,] # only on matched observed to expected
	outfile <- paste0(sample,"_scatter_loglog")
	print(paste("write",outfile))
	svg(paste0(outfile,".svg"), width = 8, height = 6)
	plot(
		log10(data_matches$expected_abund),
		log10(data_matches$observed_abund),
		xlab = "log10(Expected Abundance)",
		ylab = "log10(Observed Abundance)",
		main = "Scatter Plot: Observed vs. Expected Abundance (log-log)",
		pch = 19,
		col = "blue"
	)
	abline(a = 0, b = 1, col = "red", lty = 2)
	invisible(dev.off())

	# Rank Abundance Curves
	outfile <- paste0(sample,"_rank_abundance_curve")
	print(paste("write",outfile))
	svg(paste0(outfile,".svg"), width = 8, height = 6)
	plot(
		1:nrow(data),
		sort(data$expected_abund, decreasing = TRUE),
		type = "l",
		col = "blue",
		lwd = 2,
		xlab = "Rank",
		ylab = "Abundance",
		main = "Rank Abundance Curves",
		ylim = c(0, max(data$expected_abund, data$observed_abund))
	)
	lines(
		1:nrow(data),
		sort(data$observed_abund, decreasing = TRUE),
		col = "red",
		lwd = 2
	)
	legend(
		"topright",
		legend = c("Expected", "Observed"),
		col = c("blue", "red"),
		lwd = 2
	)
	invisible(dev.off())
}

# (5) TABLES - metrics overall

# write table
outfile <- "abundances_per-sample.tsv"
abundances$tag <- rep(tag, nrow(abundances))
print(paste("write",outfile))
write.table(abundances, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep="\t")

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

# (6) PLOT - selected metrics, overall
df_subset <- subset(df, type %in% c("recall","precision","F1","Fbeta","fdr","jaccard","pearson_cor","spearman_rho","mae","rmse","ps","bray-curtis","hellinger","jensen-shannon") )
df_subset$value <- as.numeric(df_subset$value)
outfile <- "performance_boxplot.svg"
print(paste("write",outfile))
svg(outfile, height = 8, width = 12)
par(mar=c(8, 4, 4, 2))
boxplot(value~type, data=df_subset, xlab="Type", ylab="Value", ylim = c(0, 1), las=2)
invisible(dev.off())
outfile <- "performance_boxplot.png"
print(paste("write",outfile))
png(outfile, height = 400, width = 600)
par(mar=c(8, 4, 4, 2)) #First value (8): Bottom margin (increase this if labels are still cut off); Second value (4): Left margin; Third value (4): Top margin; Fourth value (2): Right margin
boxplot(value~type, data=df_subset, xlab="Type", ylab="Value", ylim = c(0, 1), las=2)
invisible(dev.off())
