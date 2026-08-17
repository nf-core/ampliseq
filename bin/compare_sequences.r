#!/usr/bin/env Rscript

# compare_sequences.r

# Get params and files from the command line
args            <- commandArgs(trailingOnly=TRUE)
alignmentFILE   <- args[1] # expected columns: "query+target+ql+tl+qilo+qihi+tilo+tihi+gaps+mism+qstrand", see https://torognes.github.io/vsearch/misc/vsearch-userfields.7.html
obsabundFILE    <- args[2] # observed abundance*
expabundFILE    <- args[3] # expected abundance*
prefix          <- args[4] # prefix string for output files
similarity_threshold <- as.numeric(args[5]) # similarity threshold for alignmentFILE
query_or_target <- args[6] # use mismatches to query or to target?
# tab-separated file with header: first column with sequences name, following one or many columns (=samples) with numeric values (=abundance), only presence (>0)/absence are used here

# Input
if (file.exists(expabundFILE)) {
	print(paste("Compare",obsabundFILE,"to",expabundFILE))
} else {
	print(paste("Compare",obsabundFILE,"to all sequences"))
}

# PREPARE

# Read sequence alignment table
alignment = read.table( alignmentFILE, header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE)
colnames(alignment) <- c("query","target","ql","tl","qilo","qihi","tilo","tihi","gaps","mism","qstrand")

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
SAMPLES <- observed_samples

# Read expected abundance table if it exists
if (file.exists(expabundFILE)) {
	if (grepl("^# Constructed from biom file", readLines(expabundFILE, n = 1))) {
		exp = read.table( expabundFILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE, comment.char = "", skip=1 )
	} else {
		exp = read.table( expabundFILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE, comment.char = "" )
	}
	colnames(exp)[1] <- "ID"
	# extract samples to analyse
	exp_samples <- colnames(exp)[2:ncol(exp)]
	print(paste( "Expected samples:", paste( exp_samples ,collapse=",")))
	SAMPLES <- exp_samples[exp_samples %in% observed_samples]
	print(paste( "Investigate samples:", paste( SAMPLES ,collapse=",")))
	# warn if samples are missing due to incomplete overlap
	if ( !all(observed_samples %in% SAMPLES) ) {
		print(paste("WARN - sample(s) not expected (file:",expabundFILE,") and omitted:",paste(observed_samples[!observed_samples %in% SAMPLES], collapse="," )))
	}
	if ( !all(exp_samples %in% SAMPLES) ) {
		print(paste("WARN - sample(s) not observed and omitted:",paste(exp_samples[!exp_samples %in% SAMPLES], collapse="," )))
	}
	# warn about missing sequences
	alignment_target <- alignment$target[alignment$target != "*"]
	if( !all(alignment_target %in% exp$ID) ) {
		print(paste("WARN - Some expected sequences in alignment results are not in",expabundFILE,". The following are missing:",paste( unique(alignment_target[!alignment_target %in% exp$ID]),collapse=",")))
	}
	if( !all(exp$ID %in% alignment$target) ) {
		print(paste("WARN - Some expected sequences of",expabundFILE,"are not in the alignment results. The following are missing:",paste( unique(exp$ID[!exp$ID %in% alignment$target]),collapse=",")))
	}
}

# check if there are any samples to analyse
if (length(SAMPLES) == 0) {
	stop( paste0("ERROR - Found no samples to investigate. Observed samples (via sequencing data input): ", paste(observed_samples, collapse=","), ". Samples in ", expabundFILE, ": ", paste(exp_samples, collapse=",")) )
}

# Investigate alignment, required columns: query,target,ql,tl,qilo,qihi,tilo,tihi,gaps,mism
alignment$qterminalgaps <- alignment$ql - (abs(alignment$qihi-alignment$qilo)+1)
alignment$qmism <- alignment$gaps + alignment$qterminalgaps + alignment$mism
alignment$tterminalgaps <- alignment$tl - (abs(alignment$tihi-alignment$tilo)+1)
alignment$tmism <- alignment$gaps + alignment$tterminalgaps + alignment$mism
if( query_or_target=="query" ) {
	alignment$mismatch_final <- alignment$qmism
} else if( query_or_target=="target" ) {
	alignment$mismatch_final <- alignment$tmism
} else if( query_or_target=="alignment" ) {
	alignment$mismatch_final <- alignment$gaps + alignment$mism
} else {
	stop( paste("ERROR -",query_or_target,"is not valid (valid: query,target,alignment)") )
}

# order for reproducibility
alignment <- alignment[order(alignment$query, alignment$target),]

# write file
outfile <- paste0(prefix,"nucleotide-differences.tsv")
print(paste("write",outfile))
write.table(alignment, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep="\t")

# ANALYSE

# select matches above threshold
matches_above_threshold <- alignment[alignment$target !="*",]

for (sample in SAMPLES) {
	if( file.exists(expabundFILE) ) {
		# filter expected sequences in that sample (matches_above_threshold$target) with expected abundances (abundance > 0) of exp$ID
		keep_cols <- c("ID",sample)
		print(paste("Found",nrow(exp),"expected sequences overall in",expabundFILE))
		s_exp <- subset(exp, select = keep_cols)
		s_exp = s_exp[s_exp[,2] > 0,]
		print(paste("Found",nrow(s_exp),"expected sequences in sample", sample,"in",expabundFILE))
		s_matches <- matches_above_threshold[matches_above_threshold$target %in% s_exp$ID,]
	} else {
		s_matches <- matches_above_threshold
	}

	# sample will be skipped if s_matches has no data
	if( nrow(s_matches)==0 ) {
		print(paste("WARN - Skipping sample",sample,"because found no matches to expected sequences"))
		next
	} else {
		print(paste("Found",length(unique(s_matches$query)),"to be expected sequences of",length(unique(matches_above_threshold$query)),"in sample",sample))
		print(paste("Found",length(unique(s_matches$target)),"matches of",length(unique(matches_above_threshold$target)),"total in sample",sample))
	}

	# select best match (sort by mismatches and retain only first unique entry)
	s_matches <- s_matches[order(s_matches$mismatch_final, as.numeric(s_matches$mismatch_final)), ]
	s_matches <- s_matches[!duplicated(s_matches$query), ]
	print(paste("Found",length(unique(s_matches$target)),"best matches of",length(unique(matches_above_threshold$target)),"total in sample",sample))

	# filter for ASVs in that sample
	keep_cols <- c("ID",sample)
	s_observed <- subset(observed, select = keep_cols)
	# keep only observed entries in that sample (abundance > 0)
	s_observed = s_observed[s_observed[,2] > 0,]
	print(paste("Found",nrow(s_observed),"of",nrow(observed),"observed sequences in sample",sample))

	# filter alignment result by observed
	s_matches <- s_matches[s_matches$query %in% s_observed$ID,]
	print(paste("Found",nrow(s_matches),"observed sequences with match in sample",sample))
	s_below_threshold <- s_observed[!s_observed$ID %in% s_matches$query,]
	print(paste("Found",nrow(s_below_threshold),"observed sequences without match in sample",sample))

	# make barplot
	df <- as.data.frame( table(s_matches$mismatch_final) )
	df <- rbind(df, data.frame(Var1 = paste0(">=",(1-similarity_threshold)*100,"%"), Freq = length(unique(s_below_threshold$ID))))
	colnames(df) <- c("Mismatches to expected sequences","Number of sequences")

	outfile <- paste0(prefix,sample,"_nucleotide-distance_barplot.svg")
	print(paste("write",outfile))
	svg(outfile, height = 4, width = 5)
	barplot(df[,2],
		names.arg=df[,1],
		main=paste("Sample",sample,"with",length(unique(s_observed$ID)),"observed sequences\n(",length(unique(s_matches$query)),"sequences above",similarity_threshold*100,"% similarity)"),
		xlab="Mismatches to expected sequences",
		ylab="Number of sequences")
	invisible(dev.off())

	outfile <- paste0(prefix,sample,"_nucleotide-distance_barplot.png")
	print(paste("write",outfile))
	png(outfile, height = 400, width = 500)
	barplot(df[,2],
		names.arg=df[,1],
		main=paste("Sample",sample,"with",length(unique(s_observed$ID)),"observed sequences\n(",length(unique(s_matches$query)),"sequences above",similarity_threshold*100,"% similarity)"),
		xlab="Mismatches to expected sequences",
		ylab="Number of sequences")
	invisible(dev.off())

	# write file
	outfile <- paste0(prefix,sample,"_nucleotide-differences_per-sample.tsv")
	print(paste("write",outfile))
	write.table(df, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep="\t")
}
