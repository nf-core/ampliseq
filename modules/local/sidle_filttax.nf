
process SIDLE_FILTTAX {
    label 'process_single'

    conda "${projectDir}/modules/local/envs/pipesidle-0-1-0-beta.yml"
    container 'nf-core/pipesidle:0.1.0-beta'

    input:
    path(table_tofilter)
    path(table_ref)

    output:
    path("reconstructed_taxonomy.tsv"), emit: filtered
    path("reconstructed_merged.tsv")  , emit: merged
    path "versions.yml"               , emit: versions

    script:
    """
    #!/usr/bin/env Rscript

    df_tofilter <- read.table("$table_tofilter", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    colnames(df_tofilter)[1] <- "ID"

    df_ref <- read.table("$table_ref", header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 1, comment.char = "")
    colnames(df_ref)[1] <- "ID"

    df_merged <- merge(df_tofilter, df_ref, by="ID", all.x=FALSE, all.y=TRUE)
    write.table(df_merged, file = "reconstructed_merged.tsv", row.names=FALSE, sep="\t")

    df_filtered <- subset(df_tofilter, df_tofilter\$ID %in% df_ref\$ID)
    write.table(df_filtered, file = "reconstructed_taxonomy.tsv", row.names=FALSE, sep="\t")

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")) ), "versions.yml")
    """
}
