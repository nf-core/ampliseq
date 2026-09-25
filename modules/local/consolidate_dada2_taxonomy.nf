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
    val(taxlevels_input)
    val(outfile)

    output:
    path(outfile)        , emit: tsv
    path "versions.yml"  , emit: versions_consolidate_dada2_taxonomy, topic: versions

    script:
    def taxlevels = taxlevels_input ?
        'c("' + taxlevels_input.split(",").join('","') + '")' :
        'c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")'
    """
    #!/usr/bin/env Rscript

    method <- "$method"
    db_key_order <- strsplit("$db_key_order", ",", fixed = TRUE)[[1]]
    files <- strsplit("${tax_files.join(',')}", ",", fixed = TRUE)[[1]]

    # the consolidated table feeds every downstream consumer that expects one taxonomy per ASV
    # (QIIME2 import, phyloseq/TSE, SBDI export), all of which read ranks positionally or by
    # name -- so the output carries exactly these ranks, whichever database won a given ASV.
    target_ranks <- $taxlevels

    # each file's sanitized database key is the segment right before the ".tsv" extension,
    # e.g. "ASV_tax.gtdb_R07-RS207.tsv" -> "gtdb_R07-RS207" -- sanitize() (dada2_taxonomy_wf.nf)
    # already replaced every "." in a raw db_key with "_", so this is unambiguous.
    extract_db_key <- function(f) sub("^.*\\\\.([^.]+)\\\\.tsv\$", "\\\\1", basename(f))

    # everything that isn't one of these occupies a rank position. SH and BOLD_bin are
    # identifiers rather than ranks, the same distinction bin/sbdiexport.R already makes.
    is_meta <- function(cols) cols %in% c("ASV_ID", "confidence", "sequence", "database", "SH", "BOLD_bin") | grepl("_confidence\$|_exact\$", cols)

    # rank vocabulary differs by database: a database roots at Domain or at Kingdom depending on
    # what it holds, and PR2 puts Supergroup,Division,Subdivision ahead of Class where the others
    # have Phylum. Domain fills the same slot as Kingdom and Division the same slot as Phylum, so
    # those are renamed into whichever of the pair target_ranks uses; Supergroup and Subdivision
    # have no counterpart and are dropped. See docs/usage.md for the full rationale.
    rank_synonyms <- c(Domain = "Kingdom", Kingdom = "Domain", Division = "Phylum", Phylum = "Division")

    harmonize <- function(df) {
        cols <- colnames(df)
        ranks <- cols[!is_meta(cols)]
        mapped <- ifelse(ranks %in% target_ranks, ranks, rank_synonyms[ranks])
        names(mapped) <- ranks
        mapped <- mapped[!is.na(mapped) & mapped %in% target_ranks]
        # a rank without a slot in target_ranks takes its bootstrap column with it
        dropped <- setdiff(ranks, names(mapped))
        from <- setdiff(cols, c(dropped, paste0(tolower(dropped), "_confidence")))
        to <- from
        to[match(names(mapped), from)] <- unname(mapped)
        conf_from <- paste0(tolower(names(mapped)), "_confidence")
        renamed_conf <- conf_from %in% from
        to[match(conf_from[renamed_conf], from)] <- paste0(tolower(unname(mapped)), "_confidence")[renamed_conf]
        df <- df[, from, drop = FALSE]
        colnames(df) <- to
        df
    }

    # process files in declared --dada_ref_taxonomy order, not channel-arrival order (which
    # varies run to run since the per-database DADA2 tasks run in parallel) -- keeps the
    # consolidated output's column order deterministic across otherwise-identical runs.
    files <- files[ order(match(sapply(files, extract_db_key), db_key_order)) ]

    tables <- lapply(files, function(f) {
        df <- read.delim(f, sep = "\\t", header = TRUE, na.strings = "", stringsAsFactors = FALSE, check.names = FALSE)
        df <- harmonize(df)
        df\$database <- extract_db_key(f)
        df
    })

    # a database whose reference taxonomy never reaches species level for this data won't have a
    # "Species"/"species_confidence" column at all -- align every table to the union of columns
    # seen across all of them before binding, so rbind() doesn't choke on a column-count mismatch.
    all_cols <- Reduce(union, lapply(tables, colnames))
    tables <- lapply(tables, function(df) {
        missing_cols <- setdiff(all_cols, colnames(df))
        df[missing_cols] <- NA
        df[, all_cols]
    })
    combined <- do.call(rbind, tables)

    rank_cols <- intersect(target_ranks, colnames(combined))

    if (method == "score") {
        combined\$.score <- ifelse(is.na(combined\$confidence), -Inf, combined\$confidence)
    } else if (method == "most-specific") {
        # Domain/Kingdom and Division/Phylum are two names for one slot (see rank_synonyms), and
        # target_ranks can carry both. Counting per column would hand a database that fills both a
        # free point over one that fills either, so a pair scores once.
        slots <- unique(lapply(rank_cols, function(r) sort(intersect(c(r, rank_synonyms[r]), rank_cols))))
        hits <- matrix(FALSE, nrow = nrow(combined), ncol = length(slots))
        for (i in seq_along(slots)) {
            hits[, i] <- rowSums(!is.na(combined[, slots[[i]], drop = FALSE])) > 0
        }
        combined\$.score <- rowSums(hits)
    } else {
        stop(paste0("Unknown consolidation method: ", method))
    }
    combined\$.tiebreak <- match(combined\$database, db_key_order)

    # per ASV_ID: highest score first, tie broken by earliest position in db_key_order
    combined <- combined[ order(combined\$ASV_ID, -combined\$.score, combined\$.tiebreak), ]
    winners <- combined[ !duplicated(combined\$ASV_ID), ]
    winners\$.score <- NULL
    winners\$.tiebreak <- NULL

    # same column order as a single-database dada2_taxonomy.nf table, with any identifier column
    # (SH, BOLD_bin, Species_exact) kept after the ranks and the new provenance column last
    tail_cols <- c("confidence", paste0(tolower(rank_cols), "_confidence"), "sequence", "database")
    extras <- setdiff(colnames(winners), c("ASV_ID", rank_cols, tail_cols))
    winners <- winners[ , c("ASV_ID", rank_cols, extras, intersect(tail_cols, colnames(winners))) ]

    write.table(winners, file = "$outfile", sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = "."))), "versions.yml")
    """
}
