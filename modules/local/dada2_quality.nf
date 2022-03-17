process DADA2_QUALITY {
    tag "$meta"
    label 'process_low'

    conda (params.enable_conda ? "bioconductor-dada2=1.22.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta}_qual_stats.pdf"            , emit: pdf
    tuple val(meta), path("*_qual_stats.tsv"), emit: tsv
    path "versions.yml"                      , emit: versions
    path "*.args.txt"                        , emit: args
    path "*plotQualityProfile.txt"           , emit: warning

    script:
    def args = task.ext.args ?: ''
    def max_files = task.ext.max_files ?: '10000'
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(dada2))

    readfiles <- sort(list.files(".", pattern = ".fastq.gz", full.names = TRUE))

    #use only the first x files
    if ( length(readfiles) > $max_files ) {
        write.table('$max_files', file = paste0("WARNING Only $max_files of ",length(readfiles)," files were used for ${meta} plotQualityProfile.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        readfiles <- readfiles[1:$max_files]
    } else {
        write.table('$max_files', file = "$max_files files were used for ${meta} plotQualityProfile.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    }

    plot <- plotQualityProfile(readfiles$args)
    data <- plot\$data

    df <- data.frame(Cycle=character(), Count=character(), Median=character(), stringsAsFactors=FALSE)
    cycles <- sort(unique(data\$Cycle))

    #aggregate data for each sequencing cycle
    for (cycle in cycles) {
        subdata <- data[data[, "Cycle"] == cycle, ]
        score <- list()
        #convert to list to calculate median
        for (j in 1:nrow(subdata)) {score <- unlist(c(score, rep(subdata\$Score[j], subdata\$Count[j])))}
        temp = data.frame(Cycle=cycle, Count=sum(subdata\$Count), Median=median(score), stringsAsFactors=FALSE)
        df <- rbind(df, temp)
    }

    #write output
    write.table( t(df), file = paste0("${meta}_qual_stats",".tsv"), sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
    pdf(paste0("${meta}_qual_stats",".pdf"))
    plot
    dev.off()

    write.table('plotQualityProfile\t$args\nmax_files\t$max_files', file = "plotQualityProfile.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    """
}
