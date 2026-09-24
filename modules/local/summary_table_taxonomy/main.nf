process SUMMARY_TABLE_TAXONOMY {
    tag "${meta.classifier}.${meta.database}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a04c5424ce6fbf346430d99ae9f72d0bbb90e3a5cf4096df32fc1716f03973a4/data' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-dtplyr_pruned:a6608bc81b0e6546' }"

    input:
    tuple val(meta), path(tax_tsv), path(annotations)

    output:
    path("ampliseq.taxonomy.*.tsv.gz"), emit: tsv
    path "versions.yml"               , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def outfile = "ampliseq.taxonomy.${meta.classifier}.${meta.database}.tsv.gz"
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(readr)
        library(purrr)
        library(stringr)
    })

    RANKS <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

    tax_raw <- read_tsv("$tax_tsv", show_col_types = FALSE)

    if ("${meta.classifier}" == "QIIME2") {
        # QIIME2 gives one "k__X; p__Y; ..." string per ASV rather than named rank columns.
        prefix_to_rank <- c(k = "kingdom", p = "phylum", c = "class", o = "order", f = "family", g = "genus", s = "species")
        parse_taxon <- function(taxon) {
            out <- setNames(rep(NA_character_, length(RANKS)), RANKS)
            if (!is.na(taxon)) {
                for (part in strsplit(taxon, "; ", fixed = TRUE)[[1]]) {
                    m <- regmatches(part, regexec("^([a-z])__(.*)\$", trimws(part)))[[1]]
                    if (length(m) == 3 && m[2] %in% names(prefix_to_rank)) {
                        val <- trimws(m[3])
                        if (nzchar(val)) out[[prefix_to_rank[[m[2]]]]] <- val
                    }
                }
            }
            tibble::as_tibble(as.list(out))
        }
        parsed <- tax_raw[["Taxon"]] |> map(parse_taxon) |> list_rbind()
        tax <- bind_cols(tibble(asv_id = tax_raw[["Feature ID"]], confidence = tax_raw[["Confidence"]]), parsed)
    } else if ("${meta.classifier}" == "PPLACE") {
        # FORMAT_PPLACETAX's string carries no rank names, so ranks are assigned by position.
        parsed <- tax_raw[["taxonomy"]] |>
            strsplit(";") |>
            map(function(v) {
                v <- trimws(v)
                length(v) <- length(RANKS)
                tibble::as_tibble(as.list(setNames(v, RANKS)))
            }) |>
            list_rbind()
        tax <- bind_cols(tibble(asv_id = tax_raw[["ASV_ID"]], confidence = NA_real_), parsed)
    } else {
        # DADA2, SINTAX and VSEARCH-LCA already name their rank columns, PR2's non-standard names
        # included. The sequence and per-rank confidences are dropped: the first is recoverable from
        # asv_id, the second is DADA2-only and so not comparable across classifiers.
        tax <- tax_raw |>
            rename(asv_id = ASV_ID) |>
            select(-any_of("sequence")) |>
            select(-matches("_confidence\$")) |>
            rename_with(str_to_lower)

        # In a consolidated table (--consolidate_taxonomies), "database" is the database that won a
        # given ASV, not the table's own.
        if ("database" %in% colnames(tax)) {
            tax <- tax |> rename(source_database = database)
        }

        # DADA2_ADDSPECIES puts its exact matches in "Species_exact" and leaves "Species" out unless
        # assignTaxonomy's taxlevels had one. Many databases reach species through addSpecies alone,
        # so the promised kingdom..species schema needs both columns folded into one.
        if (!"species" %in% colnames(tax) && "species_exact" %in% colnames(tax)) {
            tax <- tax |> rename(species = species_exact)
        } else if ("species" %in% colnames(tax) && "species_exact" %in% colnames(tax)) {
            tax <- tax |> mutate(species = coalesce(species, species_exact)) |> select(-species_exact)
        } else if (!"species" %in% colnames(tax)) {
            tax <- tax |> mutate(species = NA_character_)
        }
    }

    annotations <- read_tsv("$annotations", show_col_types = FALSE)
    out <- left_join(tax, annotations, by = "asv_id")

    write_tsv(out, "$outfile")

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major", "minor")], collapse = ".")),
            paste0("    r-dplyr: ", packageVersion("dplyr"))
        ),
        "versions.yml"
    )
    """

    stub:
    def outfile = "ampliseq.taxonomy.${meta.classifier}.${meta.database}.tsv.gz"
    """
    echo -e "asv_id\\tkingdom\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies\\tconfidence" | gzip > ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: stub
    END_VERSIONS
    """
}
