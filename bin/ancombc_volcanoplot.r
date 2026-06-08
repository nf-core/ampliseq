#!/usr/bin/env Rscript
library("ggplot2")
library("jsonlite")

# Parse arguments with defaults
args <- commandArgs(trailingOnly = TRUE)
defaults <- list(
    lfc_file = "lfc.jsonl", #"lfc_slice.csv",
    qval_file = "q.jsonl", #"q_val_slice.csv",
    cutoff_pval = 0.05,
    cutoff_lfc = 1,
    n_labels = 0,
    pval_file = "p.jsonl", #"p_val_slice.csv",
    se_file = "se.jsonl", #"se_slice.csv"
    passed_ss_file = "passed_ss.jsonl", # only for ANCOM-BC2
    prefix = ""
)

get_arg <- function(i, type = "character") {
    val <- if (i <= length(args) && !is.na(args[i])) args[i] else defaults[[i]]
    if (type == "numeric") as.numeric(val) else val
}

lfc_file <- get_arg(1)
qval_file <- get_arg(2)
cutoff_pval <- get_arg(3, "numeric")
cutoff_lfc <- get_arg(4, "numeric")
n_labels <- get_arg(5, "numeric")
pval_file <- get_arg(6)
se_file <- get_arg(7)
passed_ss_file <- get_arg(8)
prefix <- get_arg(9)

# json reader function
read_jsonl_table <- function(file) {
    # Read all lines from the file, save metadata and data
    lines <- readLines(file)
    metadata <- jsonlite::fromJSON(lines[1], simplifyVector = FALSE)
    data_lines <- lines[-1]
    data_list <- lapply(data_lines, function(x) jsonlite::fromJSON(x, simplifyVector = FALSE))

    # Fill data frame with parsed data
    df <- data.frame(row.names = seq_along(data_list))
    for (field in metadata$fields) {
        col_name <- field$name
        col_type <- field$type
        values <- lapply(data_list, function(x) {
            val <- x[[col_name]]
            if (is.null(val)) NA else val
        })
        values <- unlist(values)

        # Convert to correct type
        if (col_type == "number") {
            df[[col_name]] <- as.numeric(values)
        } else if (col_type == "string") {
            df[[col_name]] <- as.character(values)
        } else {
            df[[col_name]] <- values
        }
    }

    # Rename first column to "id"
    colnames(df)[1] <- "id"

    return(df)
}

# auto-reader jsonl or csv
read_table_auto <- function(file) {
    ext <- tolower(tools::file_ext(file))
    # csv file
    if (ext == "csv") {
        cat("Read csv", file, "\n")
        return(read.csv(file, check.names = FALSE))
    }
    # json-style file
    if (ext %in% c("jsonl", "ndjson")) {
        cat("Read json-like", file, "\n")
        return(read_jsonl_table(file))
    }
}

# Read data
df_lfc  <- read_table_auto(lfc_file)
df_qval <- read_table_auto(qval_file)
df_pval <- read_table_auto(pval_file)
df_se   <- read_table_auto(se_file)
if (file.exists(passed_ss_file)) { df_pass <- read_table_auto(passed_ss_file) }

# Process each test
for (test in colnames(df_lfc)[-(1:2)]) {
    cat("Evaluating", test, "\n")

    # Merge data
    df <- merge( df_lfc[, c("id", test)], df_qval[, c("id", test)], by = "id" )
    colnames(df) <- c("id", "lfc", "qval")

    # Curate data
    df$qval[df$qval < 1e-100] <- 1e-100
    df$gene_type <- "ns"
    df$gene_type[df$lfc >= cutoff_lfc & df$qval <= cutoff_pval] <- "up"
    df$gene_type[df$lfc <= -cutoff_lfc & df$qval <= cutoff_pval] <- "down"

    # Label significant genes
    df_sign <- df[abs(df$lfc) >= cutoff_lfc & df$qval <= cutoff_pval, ]
    df_sign <- head(df_sign[order(-abs(df_sign$lfc)), ], n_labels)
    if (nrow(df_sign) > 0) {
        df_sign$label <- sapply(df_sign$id, function(x) {
            tail(strsplit(as.character(x), ";")[[1]], 1)
        })
    }

    # Create volcano plot
    p <- ggplot(df, aes(x = lfc, y = -log10(qval), fill = gene_type, size = gene_type, alpha = gene_type)) +
        geom_point(shape = 21, colour = "black") +
        geom_hline(yintercept = -log10(cutoff_pval), linetype = "dashed") +
        geom_vline(xintercept = c(-cutoff_lfc, cutoff_lfc), linetype = "dashed") +
        scale_fill_manual(values = c(up = "#4c78a8", down = "#f58518", ns = "grey"), name = "Change:") +
        scale_size_manual(values = c(up = 1, down = 1, ns = 0.5), name = "Change:") +
        scale_alpha_manual(values = c(up = 1, down = 1, ns = 0.5), name = "Change:") +
        labs(title = paste0(prefix, test), x = "log2(fold change)", y = "-log10(adjusted P-value)") +
        theme_bw()
    # Add labels
    if (n_labels > 0 && nrow(df_sign) > 0) {
        p <- p + geom_label(data = df_sign, aes(label = label), nudge_y = 2, alpha = 0.2, size = 1)
    }

    # Save plot
    svg(paste0(prefix, test, ".volcano_plot.svg"), height = 3.6, width = 3.6)
    print(p)
    dev.off()

    # Save data table
    df_full <- merge(df, df_pval[, c("id", test)], by = "id")
    df_full <- merge(df_full, df_se[, c("id", test)], by = "id")
    colnames(df_full) <- c("id", "lfc", "qval", "type", "pval", "se")
    df_full <- df_full[, c("id", "lfc", "se", "pval", "qval", "type")]
    if (file.exists(passed_ss_file)) {
        df_full <- merge(df_full, df_pass[, c("id", test)], by = "id")
        colnames(df_full) <- c("id", "lfc", "se", "pval", "qval", "type", "passed_ss")
    }
    write.table(df_full, paste0(prefix, test, ".volcano_plot.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}
