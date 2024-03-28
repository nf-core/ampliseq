process DADA2_QUALITY {
    tag "$meta"
    label 'process_low'

    conda "bioconda::bioconductor-dada2=1.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0' :
        'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    path "*_qual_stats.pdf"                  , emit: pdf
    path "*_qual_stats.svg"                  , emit: svg
    tuple val(meta), path("*_qual_stats.tsv"), emit: tsv
    path "versions.yml"                      , emit: versions
    path "*.args.txt"                        , emit: args
    path "*plotQualityProfile.txt"           , emit: warning

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(dada2))
    suppressPackageStartupMessages(library(ShortRead))

    readfiles <- sort(list.files(".", pattern = ".fastq.gz", full.names = TRUE))

    #make list of number of sequences
    readfiles_length <- countLines(readfiles) / 4
    sum_readfiles_length <- sum(readfiles_length)

    #use only the first x files when read number gets above 2147483647, read numbers above that do not fit into an INT and crash the process!
    if ( sum_readfiles_length > 2147483647 ) {
        max_files = length(which(cumsum(readfiles_length) <= 2147483647 ))
        write.table(max_files, file = paste0("WARNING Only ",max_files," of ",length(readfiles)," files and ",sum(readfiles_length[1:max_files])," of ",sum_readfiles_length," reads were used for ${meta} plotQualityProfile.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        readfiles <- readfiles[1:max_files]
    } else {
        max_files <- length(readfiles)
        write.table(max_files, file = paste0(max_files," files were used for ${prefix} plotQualityProfile.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    }

    plot <- plotQualityProfile(readfiles$args)
    data <- plot\$data

    #aggregate data for each sequencing cycle
    df <- data.frame(Cycle=character(), Count=character(), Median=character(), stringsAsFactors=FALSE)
    cycles <- sort(unique(data\$Cycle))
    for (cycle in cycles) {
        subdata <- data[data[, "Cycle"] == cycle, ]
        score <- list()
        #convert to list to calculate median
        for (j in 1:nrow(subdata)) {score <- unlist(c(score, rep(subdata\$Score[j], subdata\$Count[j])))}
        temp = data.frame(Cycle=cycle, Count=sum(subdata\$Count), Median=median(score), stringsAsFactors=FALSE)
        df <- rbind(df, temp)
    }

    #write output
    write.table( t(df), file = paste0("${prefix}_qual_stats",".tsv"), sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
    pdf(paste0("${prefix}_qual_stats",".pdf"))
    plot
    dev.off()
    svg(paste0("${prefix}_qual_stats",".svg"))
    plot
    dev.off()

    write.table(paste0('plotQualityProfile\t$args\nmax_files\t',max_files), file = "${prefix}_plotQualityProfile.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")),paste0("    ShortRead: ", packageVersion("ShortRead")) ), "versions.yml")
    """
}
