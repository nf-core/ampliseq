process SUMMARY_TABLE_COUNTS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a04c5424ce6fbf346430d99ae9f72d0bbb90e3a5cf4096df32fc1716f03973a4/data' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-dtplyr_pruned:a6608bc81b0e6546' }"

    input:
    path(asv_table)

    output:
    path("ampliseq.counts.tsv.gz"), emit: tsv
    path "versions.yml"           , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages({
        library(dplyr)
        library(tidyr)
        library(readr)
    })

    counts <- read_tsv("$asv_table", show_col_types = FALSE) |>
        rename(asv_id = ASV_ID) |>
        pivot_longer(-asv_id, names_to = "sample", values_to = "count") |>
        filter(count > 0)

    write_tsv(counts, "ampliseq.counts.tsv.gz")

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
    """
    echo -e "asv_id\\tsample\\tcount" | gzip > ampliseq.counts.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: stub
    END_VERSIONS
    """
}
