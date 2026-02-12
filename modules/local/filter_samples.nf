process FILTER_SAMPLES {
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    path(metadata, stageAs: 'input/*')
    path(table, stageAs: 'input/*')

    output:
    path("metadata.tsv"), emit: metadata
    path("table.tsv")   , emit: abundances
    path("*.log")       , emit: log, optional: true
    path "versions.yml" , emit: versions

    script:
    """
    #!/usr/bin/env Rscript

    # first column in meta has sample id
    meta <- read.table( "$metadata", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # column names are sample ids, but first column is asv id
    abund <- read.table( "$table", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    # samples that arent in both files are dropped
    meta_filtered <- meta[meta[,1] %in% colnames(abund)[2:length(colnames(abund))],]
    abund_filtered <- abund[,colnames(abund) %in% c( colnames(abund)[1], meta[,1] ) ]

    # write filtered data
    write.table(meta_filtered, file = "metadata.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep = "\t")
    write.table(abund_filtered, file = "table.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep = "\t")

    # error in case all samples were removed
    if ( nrow(meta_filtered) == 0 ) {
        stop("All samples were removed. That means no overlap between the metadata sample IDs and the abundance table sample IDs was found. Make sure that sample IDs match.")
    }

    # this is in case some samples were lost during preprocessing, i.e. samples in metadata but not in abundance table
    if ( nrow(meta) > nrow(meta_filtered) ) {
        log_message = paste("The metadata file rows were reduced from", nrow(meta), "to", nrow(meta_filtered),", because some samples were missing in the abundance table")
        write.table(log_message, file = paste0(log_message,".log"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    # this is in case some samples were not in metadata, i.e. only a subset of samples is entering downstream analysis
    if ( ncol(abund) > ncol(abund_filtered) ) {
        log_message = paste("Samples in the abundance file were reduced from", ncol(abund)-1, "to", ncol(abund_filtered)-1,", because the metadata did not contain all samples in the abundance table")
        write.table(log_message, file = paste0(log_message,".log"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    }

    # versions
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")) ), "versions.yml")
    """
}
