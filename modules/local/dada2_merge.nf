process DADA2_MERGE {
    label 'process_low'

    conda (params.enable_conda ? "bioconductor-dada2=1.22.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

    input:
    path(files)
    path(rds)

    output:
    path( "DADA2_stats.tsv" ), emit: dada2stats
    path( "DADA2_table.tsv" ), emit: dada2asv
    path( "ASV_table.tsv" )  , emit: asv
    path( "ASV_seqs.fasta" ) , emit: fasta
    path( "DADA2_table.rds" ), emit: rds
    path "versions.yml"      , emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    suppressPackageStartupMessages(library(digest))

    #combine stats files
    for (data in sort(list.files(".", pattern = ".stats.tsv", full.names = TRUE))) {
        if (!exists("stats")){ stats <- read.csv(data, header=TRUE, sep="\\t") }
        if (exists("stats")){
            temp <-read.csv(data, header=TRUE, sep="\\t")
            stats <-unique(rbind(stats, temp))
            rm(temp)
        }
    }
    write.table( stats, file = "DADA2_stats.tsv", sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

    #combine dada-class objects
    files <- sort(list.files(".", pattern = ".ASVtable.rds", full.names = TRUE))
    if ( length(files) == 1 ) {
        ASVtab = readRDS(files[1])
    } else {
        ASVtab <- mergeSequenceTables(tables=files, repeats = "error", orderBy = "abundance", tryRC = FALSE)
    }
    saveRDS(ASVtab, "DADA2_table.rds")

    df <- t(ASVtab)
    colnames(df) <- gsub('_1.filt.fastq.gz', '', colnames(df))
    colnames(df) <- gsub('.filt.fastq.gz', '', colnames(df))
    df <- data.frame(sequence = rownames(df), df)
    # Create an md5 sum of the sequences as ASV_ID and rearrange columns
    df\$ASV_ID <- sapply(df\$sequence, digest, algo='md5', serialize = FALSE)
    df <- df[,c(ncol(df),3:ncol(df)-1,1)]

    # file to publish
    write.table(df, file = "DADA2_table.tsv", sep = "\\t", row.names = FALSE, quote = FALSE, na = '')

    # Write fasta file with ASV sequences to file
    write.table(data.frame(s = sprintf(">%s\n%s", df\$ASV_ID, df\$sequence)), 'ASV_seqs.fasta', col.names = FALSE, row.names = FALSE, quote = FALSE, na = '')

    # Write ASV file with ASV abundances to file
    df\$sequence <- NULL
    write.table(df, file = "ASV_table.tsv", sep="\\t", row.names = FALSE, quote = FALSE, na = '')

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    """
}
