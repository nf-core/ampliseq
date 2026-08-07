#!/usr/bin/env Rscript

# Get params and files from the command line
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
	stop("Usage: Rscript compare_profiles.r <tag> <obsFILE> <expFILE> [fbeta_val] [mode]")
}
tag <- args[1] # tag for each sample
obsFILE <- args[2] # observed abundance*, e.g. "abs-abund-table-3.tsv"
expFILE <- args[3] # expected abundance*, e.g. "expected_taxonomy.tsv"
fbeta_val <- if (length(args) >= 4) as.numeric(args[4]) else 2 # Fbeta weight, default 2; 1=precision and recall equally weighted (=F1 score), 2=weighs recall higher than precision, 0.5=weighs recall lower than precision.
mode <- if (length(args) >= 5) args[5] else "complete" # mode: "complete" for full taxonomic string, "level" for only the label at the given level

# Args check
if (!(mode %in% c("complete", "level"))) {
	stop(paste("Invalid mode:", mode, ". Must be either 'complete' or 'level'."))
}

# Load data frames from TSV files - with potential biom header
message("Loading observed taxonomies from ",obsFILE)
if (grepl("^# Constructed from biom file", readLines(obsFILE, n = 1))) {
	observed_taxonomy <- read.table(obsFILE, header = TRUE, sep = "\t", skip = 1, comment.char = "", check.names = FALSE)
} else {
	observed_taxonomy <- read.table(obsFILE, header = TRUE, sep = "\t", skip = 0, comment.char = "", check.names = FALSE)
}
message("Loading expected taxonomies from ",expFILE)
if (grepl("^# Constructed from biom file", readLines(expFILE, n = 1))) {
	expected_taxonomy <- read.table(expFILE, header = TRUE, sep = "\t", skip = 1, comment.char = "", check.names = FALSE)
} else {
	expected_taxonomy <- read.table(expFILE, header = TRUE, sep = "\t", skip = 0, comment.char = "", check.names = FALSE)
}

# Determine common samples
common_samples <- intersect(colnames(observed_taxonomy)[-1], colnames(expected_taxonomy)[-1])
if (length(common_samples) == 0) {
	stop("ERROR: No common samples found between observed and expected taxonomy files.")
}
message(paste("Processing", length(common_samples), "common samples:",paste(common_samples,collapse=",")))

# Function to get taxonomic levels based on mode (e.g., complete: L1: x1, L2: x1;x2; level: L1: x1, L2: x2)
get_taxa_by_mode <- function(taxa_vec, level, mode) {
	s <- strsplit(as.character(taxa_vec), ";")
	sapply(s, function(x) {
		if (length(x) >= level && !is.na(x[level]) && x[level] != "") {
			if (mode == "complete") {
				if (all(x[1:level] != "" & !is.na(x[1:level]))) {
					return(paste(x[1:level], collapse = ";"))
				}
			} else if (mode == "level") {
				return(x[level])
			}
		}
		NA
	})
}

# Function for classification metrics 
calculate_classification_metrics <- function(tp, fp, fn, fbeta_val) {
	recall <- if ((tp + fn) > 0) tp / (tp + fn) else NA
	precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA
	f1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) 2 * (precision * recall) / (precision + recall) else NA
	fbeta_res <- if (!is.na(precision) && !is.na(recall) && (fbeta_val^2 * precision + recall) > 0) (1 + fbeta_val^2) * precision * recall / (fbeta_val^2 * precision + recall) else NA
	fdr <- if ((tp + fp) > 0) fp / (tp + fp) else NA
	jaccard <- if ((tp + fp + fn) > 0) tp / (tp + fp + fn) else NA
	return(list(recall=recall, precision=precision, f1=f1, fbeta=fbeta_res, fdr=fdr, jaccard=jaccard))
}

# Function for correlation metrics
calculate_correlation_metrics <- function(tp_set, obs_taxa_all_lev, exp_taxa_all_lev, obs_vec, exp_vec) {
	p_cor <- NA
	s_rho <- NA
	if (length(tp_set) > 1) {
		obs_aggregated <- tapply(obs_vec, obs_taxa_all_lev, sum, NA.rm=TRUE)
		exp_aggregated <- tapply(exp_vec, exp_taxa_all_lev, sum, NA.rm=TRUE)
		obs_vals <- obs_aggregated[tp_set]
		exp_vals <- exp_aggregated[tp_set]
		if (sd(obs_vals, na.rm=TRUE) > 0 && sd(exp_vals, na.rm=TRUE) > 0) {
			p_cor <- cor(obs_vals, exp_vals, method = "pearson")
			s_rho <- cor(obs_vals, exp_vals, method = "spearman")
		}
	}
	return(list(pearson_cor=p_cor, spearman_rho=s_rho))
}

# Function for distance metrics
calculate_distance_metrics <- function(obs_taxa_at_lev, exp_taxa_at_lev, obs_taxa_all_lev, obs_vec, exp_vec) {
	all_taxa <- union(obs_taxa_at_lev, exp_taxa_at_lev)
	res <- list(bray_curtis=NA, hellinger=NA, jensen_shannon=NA, deviation=NA, mae=NA, rmse=NA, ps=NA)

	if (length(all_taxa) > 0) {
		obs_all_vals <- tapply(obs_vec, obs_taxa_all_lev, sum, NA.rm=TRUE)[all_taxa]
		exp_all_vals <- tapply(exp_vec, exp_taxa_all_lev, sum, NA.rm=TRUE)[all_taxa]
		obs_all_vals[is.na(obs_all_vals)] <- 0
		exp_all_vals[is.na(exp_all_vals)] <- 0

		sum_obs <- sum(obs_all_vals)
		sum_exp <- sum(exp_all_vals)

		if (sum_obs > 0 && sum_exp > 0) {
			obs_rel <- obs_all_vals / sum_obs
			exp_rel <- exp_all_vals / sum_exp

			res$bray_curtis <- sum(abs(obs_rel - exp_rel)) / 2
			res$hellinger <- sqrt(sum((sqrt(obs_rel) - sqrt(exp_rel))^2))

			kl_div <- function(x, y, pseudocount = 1e-10) {
				x <- x + pseudocount
				y <- y + pseudocount
				x <- x / sum(x)
				y <- y / sum(y)
				sum(x * log(x / y))
			}
			m <- 0.5 * (obs_rel + exp_rel)
			res$jensen_shannon <- sqrt(0.5 * kl_div(obs_rel, m) + 0.5 * kl_div(exp_rel, m))

			res$ps <- sum(pmin(obs_rel, exp_rel)) / sum(exp_rel)

			valid_exp_rel <- exp_rel > 0
			if (any(valid_exp_rel)) {
				perc_dev <- abs((obs_rel[valid_exp_rel] - exp_rel[valid_exp_rel]) / exp_rel[valid_exp_rel]) * 100
				res$deviation <- median(perc_dev, na.rm = TRUE)
			}
			res$mae <- mean(abs(obs_rel - exp_rel))
			res$rmse <- sqrt(mean((obs_rel - exp_rel)^2))
		}
	}
	return(res)
}

# Determine max taxonomic depth across both files
max_obs_lev <- max(sapply(strsplit(as.character(observed_taxonomy[, 1]), ";"), length))
max_exp_lev <- max(sapply(strsplit(as.character(expected_taxonomy[, 1]), ";"), length))
max_lev <- max(max_obs_lev, max_exp_lev)
message(paste("Processing", max_lev, "taxonomic levels"))

# Use a list to store results
results_list <- list()
res_counter <- 1

for (s in common_samples) {
	idx_obs <- which(colnames(observed_taxonomy) == s)
	idx_exp <- which(colnames(expected_taxonomy) == s)

	for (l in 1:max_lev) {
		level_name <- paste0("L", l)

		# Taxon strings for this level for all entries based on the specified mode
		obs_taxa_all_lev <- get_taxa_by_mode(observed_taxonomy[, 1], l, mode)
		exp_taxa_all_lev <- get_taxa_by_mode(expected_taxonomy[, 1], l, mode)

		# Observed taxa present in this sample at this level
		obs_taxa_at_lev <- unique(obs_taxa_all_lev[observed_taxonomy[, idx_obs] > 0])
		obs_taxa_at_lev <- obs_taxa_at_lev[!is.na(obs_taxa_at_lev) & obs_taxa_at_lev != ""]

		# Expected taxa present in this sample at this level
		exp_taxa_at_lev <- unique(exp_taxa_all_lev[expected_taxonomy[, idx_exp] > 0])
		exp_taxa_at_lev <- exp_taxa_at_lev[!is.na(exp_taxa_at_lev) & exp_taxa_at_lev != ""]

		# Calculate sets
		tp_set <- intersect(obs_taxa_at_lev, exp_taxa_at_lev)
		fp_set <- setdiff(obs_taxa_at_lev, exp_taxa_at_lev)
		fn_set <- setdiff(exp_taxa_at_lev, obs_taxa_at_lev)

		# Calculate counts
		tp <- length(tp_set)
		fp <- length(fp_set)
		fn <- length(fn_set)
		exp_count <- length(exp_taxa_at_lev)
		obs_count <- length(obs_taxa_at_lev)

		# Classification metrics - based on counts
		metrics_class <- calculate_classification_metrics(tp, fp, fn, fbeta_val)

		# Correlation calculations
		metrics_corr <- calculate_correlation_metrics(tp_set, obs_taxa_all_lev, exp_taxa_all_lev, observed_taxonomy[, idx_obs], expected_taxonomy[, idx_exp])

		# Distance and Error metrics
		metrics_dist <- calculate_distance_metrics(obs_taxa_at_lev, exp_taxa_at_lev, obs_taxa_all_lev, observed_taxonomy[, idx_obs], expected_taxonomy[, idx_exp])

		# Sequence of metrics as requested
		metrics_seq <- list(
			list(type = "exp", value = exp_count),
			list(type = "obs", value = obs_count),
			list(type = "TP",	value = tp),
			list(type = "FN",	value = fn),
			list(type = "FP",	value = fp),
			list(type = "recall", value = metrics_class$recall),
			list(type = "precision", value = metrics_class$precision),
			list(type = "F1", value = metrics_class$f1),
			list(type = "Fbeta", value = metrics_class$fbeta),
			list(type = "fdr", value = metrics_class$fdr),
			list(type = "jaccard", value = metrics_class$jaccard),
			list(type = "TPs_exp", value = paste(head(tp_set, 100), collapse = ",")),
			list(type = "FNs_exp", value = paste(head(fn_set, 100), collapse = ",")),
			list(type = "FPs_obs", value = paste(head(fp_set, 100), collapse = ",")),
			list(type = "bray-curtis", value = metrics_dist$bray_curtis),
			list(type = "hellinger", value = metrics_dist$hellinger),
			list(type = "jensen-shannon", value = metrics_dist$jensen_shannon),
			list(type = "pearson_cor", value = metrics_corr$pearson_cor),
			list(type = "spearman_rho", value = metrics_corr$spearman_rho),
			list(type = "deviation", value = metrics_dist$deviation),
			list(type = "mae", value = metrics_dist$mae),
			list(type = "rmse", value = metrics_dist$rmse),
			list(type = "ps", value = metrics_dist$ps)
		)
			for (m in metrics_seq) {
			results_list[[res_counter]] <- data.frame(
				sample = s,
				level = level_name,
				type = m$type,
				value = as.character(m$value),
				stringsAsFactors = FALSE
			)
			res_counter <- res_counter + 1
		}
	}
}

# Write detailed output
results <- do.call(rbind, results_list)
results$tag <- rep(tag, nrow(results))
outfile <- "profile_per-sample.tsv"
message(paste("Writing detailed output to:", outfile))
write.table(results, outfile, sep = "\t", row.names = FALSE, quote = FALSE)

# SUMMARY TABLE
outfile <- "profile_summary.tsv"
message(paste("Writing summary table to:", outfile))
df_sum <- data.frame(level=character(), type=character(), median=numeric(), mean=numeric(), standard_error=numeric(), min=numeric(), max=numeric(), count=numeric(), stringsAsFactors=FALSE)

# Convert values to numeric for summary
results_num <- results
results_num$value <- as.numeric(results_num$value)

for (lev in unique(results_num$level)) {
	for (type in unique(results_num$type)) {
		if (type %in% c("TPs_exp","FNs_exp","FPs_obs")) next #exclude non-helpful entries
		df_subset <- results_num[results_num$level == lev & results_num$type == type, ]
		if (nrow(df_subset) == 0) next
		COUNT <- sum( !is.na(df_subset$value) ) # counts all non-na values
		MEDIAN <- median( df_subset$value, na.rm = TRUE )
		MEAN <- mean( df_subset$value, na.rm = TRUE )
		SE <- sd( df_subset$value, na.rm = TRUE ) / sqrt(length(df_subset$value))
		MIN <- ifelse(COUNT==0, NA, min( df_subset$value, na.rm = TRUE ) ) # avoid -Inf when there are no values
		MAX <- ifelse(COUNT==0, NA, max( df_subset$value, na.rm = TRUE ) ) # avoid Inf when there are no values
		df_sum[nrow(df_sum) + 1,] <- list(lev, type, MEDIAN, MEAN, SE, MIN, MAX, COUNT)
	}
}
df_sum$tag <- rep(tag, nrow(df_sum))
write.table(df_sum, outfile, sep = "\t", row.names = FALSE, quote = FALSE)

# LINE PLOT PER LEVEL
plot_metrics <- c("F1", "recall", "precision")
metric_colors <- c("F1" = "black", "recall" = "red", "precision" = "blue")

all_levels <- unique(df_sum$level)
levels_num <- 1:length(all_levels)

draw_line_plot <- function() {
	plot(NULL, xlim = c(1, length(all_levels)), ylim = c(0, 1),
		xaxt = "n", xlab = "Level", ylab = "Median",
		main = "Profile Summary Metrics")
	axis(1, at = levels_num, labels = all_levels)

	for (m in plot_metrics) {
		subset_df <- subset(df_sum, type == m)
		if (nrow(subset_df) > 0) {
			y <- subset_df$median # alternative: subset_df$mean
			se <- subset_df$standard_error
			x <- 1:nrow(subset_df)
			lines(x, y, col = metric_colors[m], lwd = 2, type = "b", pch = 19)
			segments(x, y - se, x, y + se, col = metric_colors[m])
		}
	}
	legend("bottomleft", legend = names(metric_colors), col = metric_colors, lty = 1, pch = 19, lwd = 2)
}

# Plot SVG
outfile_svg <- "profile_lineplot.svg"
message(paste("Writing lineplot SVG to:", outfile_svg))
svg(outfile_svg, height = 8, width = 12)
draw_line_plot()
invisible(dev.off())

# Plot PNG
outfile_png <- "profile_lineplot.png"
message(paste("Writing lineplot PNG to:", outfile_png))
png(outfile_png, height = 400, width = 600)
draw_line_plot()
invisible(dev.off())
