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
        # Native QIIME2 classify-sklearn output: "Feature ID", "Taxon", "Confidence" -- Taxon is a
        # "k__X; p__Y; ..." rank-letter-prefixed, semicolon+space-joined string, not already-named
        # rank columns like the other classifiers.
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
        # FORMAT_PPLACETAX output: "ASV_ID", "taxonomy" -- an unranked, unprefixed semicolon-joined
        # string, best-effort positionally mapped kingdom..species (see docs/output.md caveat).
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
        # DADA2 / SINTAX / VSEARCH-LCA -- already-named rank columns (Kingdom..Species, optionally
        # non-standard names like PR2's Supergroup/Division/Subdivision, kept as-is), a single
        # most-specific-rank confidence column, optionally per-rank *_confidence columns and/or an
        # SH column (--addsh) and a sequence column. Slimmed per the confirmed summary-table design:
        # drop sequence (recoverable via asv_id) and the per-rank confidence columns (DADA2-specific
        # detail not comparable across classifiers), keep the single confidence.
        tax <- tax_raw |>
            rename(asv_id = ASV_ID) |>
            select(-any_of("sequence")) |>
            select(-matches("_confidence\$")) |>
            rename_with(str_to_lower)

        # DADA2_ADDSPECIES (--dada_addspecies) always renames its own exact-match species call to
        # "Species_exact", and only re-adds a native "Species" column (assignTaxonomy's own call)
        # when assignTaxonomy's taxlevels already included one -- many databases rely on addSpecies
        # alone for species-level calls, so without this the "species" column would be silently
        # absent for those, breaking the promised consistent kingdom..species schema. Prefer the
        # native assignTaxonomy species call when present, fall back to addSpecies' exact-match call
        # otherwise, keep the table to one species column rather than exposing both.
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
