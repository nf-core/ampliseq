process CONSOLIDATE_DADA2_TAXONOMY {
    tag "${method}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
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
    # (QIIME2 import, phyloseq/TSE), all of which read ranks positionally or by name -- so the
    # output carries exactly these ranks, whichever database won a given ASV.
    target_ranks <- $taxlevels

    # each file's sanitized database key is the segment right before the ".tsv" extension,
    # e.g. "ASV_tax.gtdb_R07-RS207.tsv" -> "gtdb_R07-RS207" -- sanitize() (dada2_taxonomy_wf.nf)
    # already replaced every "." in a raw db_key with "_", so this is unambiguous.
    extract_db_key <- function(f) sub("^.*\\\\.([^.]+)\\\\.tsv\$", "\\\\1", basename(f))

    # everything that isn't one of these occupies a rank position. SH and BOLD_bin are
    # identifiers rather than ranks, the same distinction bin/parse_dada2_taxonomy.r makes.
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
        # a synonym fills its slot only when no column of that name exists, so a database
        # holding both Domain and Kingdom keeps the one the target uses and drops the other
        mapped[!(ranks %in% target_ranks) & mapped %in% ranks] <- NA
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
        # character, so a rank that happens to be only "F" or "T" is not read as logical
        df <- read.delim(f, sep = "\\t", header = TRUE, na.strings = "", colClasses = "character", check.names = FALSE)
        conf_cols <- grep("^confidence\$|_confidence\$", colnames(df), value = TRUE)
        df[conf_cols] <- lapply(df[conf_cols], as.numeric)
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

    if (!(method %in% c("most-specific", "score"))) {
        stop(paste0("Unknown consolidation method: ", method))
    }

    # Domain/Kingdom and Division/Phylum are two names for one slot (see rank_synonyms), and
    # target_ranks can carry both, so a pair counts as one rank.
    # SH and BOLD_bin are identifiers, not ranks.
    score_cols <- setdiff(rank_cols, c("SH", "BOLD_bin"))
    slots <- unique(lapply(score_cols, function(r) sort(intersect(c(r, rank_synonyms[r]), score_cols))))
    assigned <- matrix(FALSE, nrow = nrow(combined), ncol = length(slots))
    slot_conf <- matrix(NA_real_, nrow = nrow(combined), ncol = length(slots))
    for (i in seq_along(slots)) {
        assigned[, i] <- rowSums(!is.na(combined[, slots[[i]], drop = FALSE])) > 0
        conf_cols <- intersect(paste0(tolower(slots[[i]]), "_confidence"), colnames(combined))
        if (length(conf_cols) > 0) {
            slot_conf[, i] <- suppressWarnings(apply(combined[, conf_cols, drop = FALSE], 1, max, na.rm = TRUE))
        }
    }
    # DADA2 reports a bootstrap for unassigned ranks too; only assigned ranks are compared
    slot_conf[!assigned | is.infinite(slot_conf)] <- NA
    comparable <- assigned & !is.na(slot_conf)
    depth <- apply(assigned, 1, function(x) if (any(x)) max(which(x)) else 0)
    tiebreak <- match(combined\$database, db_key_order)

    # bootstraps only fall with depth, so each database's own deepest rank would favour shallow
    # assignments; compare at the deepest rank that all remaining candidates assigned instead
    pick_by_score <- function(rows) {
        cand <- rows[depth[rows] > 0]
        if (length(cand) == 0) cand <- rows
        previous <- 0
        while (length(cand) > 1) {
            shared <- which(colSums(!comparable[cand, , drop = FALSE]) == 0)
            if (length(shared) == 0 || max(shared) <= previous) break
            previous <- max(shared)
            conf <- slot_conf[cand, previous]
            cand <- cand[conf == max(conf)]
            deeper <- cand[depth[cand] > previous]
            if (length(deeper) == 0) break
            cand <- deeper
        }
        cand <- cand[depth[cand] == max(depth[cand])]
        cand[which.min(tiebreak[cand])]
    }

    rows_by_asv <- split(seq_len(nrow(combined)), combined\$ASV_ID)
    if (method == "score") {
        winner_rows <- vapply(rows_by_asv, pick_by_score, integer(1))
    } else {
        # most slots filled first, tie broken by earliest position in db_key_order
        winner_rows <- vapply(rows_by_asv, function(rows) rows[order(-rowSums(assigned[rows, , drop = FALSE]), tiebreak[rows])][1], integer(1))
    }
    winners <- combined[winner_rows, ]

    # same column order as a single-database dada2_taxonomy.nf table, with any identifier column
    # (SH, BOLD_bin, Species_exact) kept after the ranks and the new provenance column last
    tail_cols <- c("confidence", paste0(tolower(rank_cols), "_confidence"), "sequence", "database")
    extras <- setdiff(colnames(winners), c("ASV_ID", rank_cols, tail_cols))
    winners <- winners[ , c("ASV_ID", rank_cols, extras, intersect(tail_cols, colnames(winners))) ]

    write.table(winners, file = "$outfile", sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = "."))), "versions.yml")
    """
}
