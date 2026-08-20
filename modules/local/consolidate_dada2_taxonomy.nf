process CONSOLIDATE_DADA2_TAXONOMY {
    tag "${method}"
    label 'process_low'

    conda "bioconda::bioconductor-dada2=1.38.0 conda-forge::r-base=4.5.2 conda-forge::r-digest=0.6.39 conda-forge::tbb=2022.3.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81153df5d53322e6d91b2c4c9bc4da50774fb1d101ead002fe75bb75fc6f036c/data' :
        'community.wave.seqera.io/library/bioconductor-dada2_r-base_r-digest_tbb:38acac09bac46f36' }"

    input:
    path(tax_files)
    val(method)
    val(db_key_order)
    val(outfile)

    output:
    path(outfile)        , emit: tsv
    path "versions.yml"  , emit: versions_consolidate_dada2_taxonomy, topic: versions

    script:
    """
    #!/usr/bin/env Rscript

    method <- "$method"
    db_key_order <- strsplit("$db_key_order", ",", fixed = TRUE)[[1]]
    files <- strsplit("${tax_files.join(',')}", ",", fixed = TRUE)[[1]]

    # each file's sanitized database key is the segment right before the ".tsv" extension,
    # e.g. "ASV_tax.gtdb_R07-RS207.tsv" -> "gtdb_R07-RS207" -- sanitize() (dada2_taxonomy_wf.nf)
    # already replaced every "." in a raw db_key with "_", so this is unambiguous.
    extract_db_key <- function(f) sub("^.*\\\\.([^.]+)\\\\.tsv\$", "\\\\1", basename(f))

    # process files in declared --dada_ref_taxonomy order, not channel-arrival order (which
    # varies run to run since the per-database DADA2 tasks run in parallel) -- keeps the
    # consolidated output's column order deterministic across otherwise-identical runs.
    files <- files[ order(match(sapply(files, extract_db_key), db_key_order)) ]

    tables <- lapply(files, function(f) {
        df <- read.delim(f, sep = "\\t", header = TRUE, na.strings = "", stringsAsFactors = FALSE, check.names = FALSE)
        df\$database <- extract_db_key(f)
        df
    })

    # different listed databases can resolve a different set of ranks for this data (e.g. a
    # database whose reference taxonomy never reaches species level for any ASV won't have a
    # "Species"/"species_confidence" column at all) -- align every table to the union of columns
    # seen across all of them before binding, so rbind() doesn't choke on a column-count mismatch.
    all_cols <- Reduce(union, lapply(tables, colnames))
    tables <- lapply(tables, function(df) {
        missing_cols <- setdiff(all_cols, colnames(df))
        df[missing_cols] <- NA
        df[, all_cols]
    })
    combined <- do.call(rbind, tables)

    rank_cols <- intersect(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), colnames(combined))

    if (method == "score") {
        combined\$.score <- ifelse(is.na(combined\$confidence), -Inf, combined\$confidence)
    } else if (method == "most-specific") {
        combined\$.score <- rowSums(!is.na(combined[, rank_cols, drop = FALSE]))
    } else {
        stop(paste0("Unknown consolidation method: ", method))
    }
    combined\$.tiebreak <- match(combined\$database, db_key_order)

    # per ASV_ID: highest score first, tie broken by earliest position in db_key_order
    combined <- combined[ order(combined\$ASV_ID, -combined\$.score, combined\$.tiebreak), ]
    winners <- combined[ !duplicated(combined\$ASV_ID), ]
    winners\$.score <- NULL
    winners\$.tiebreak <- NULL

    # keep original column order, with the new provenance column appended at the end
    orig_cols <- setdiff(colnames(winners), "database")
    winners <- winners[ , c(orig_cols, "database") ]

    write.table(winners, file = "$outfile", sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = "."))), "versions.yml")
    """
}
