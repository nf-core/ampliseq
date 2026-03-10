process FILTER_SEQUENCES_ABUNDANCES {
    tag "${table}-${fasta}"
    label 'process_single'

    conda "bioconda::bioconductor-biostrings=2.58.0 conda-forge::r-base=4.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-biostrings:2.58.0--r40h037d062_0' :
        'biocontainers/bioconductor-biostrings:2.58.0--r40h037d062_0' }"

    input:
    path(fasta)
    path(table)

    output:
    path( "filtered_sequences.fasta" ), emit: seq
    path( "filtered_abundances.tsv" ) , emit: abund
    path( "stats_counts.tsv" )        , emit: stats_counts
    path( "stats_features.tsv" )      , emit: stats_features
    path "versions.yml"               , emit: versions

    script:
    """
    #!/usr/bin/env Rscript

    #load packages
    suppressPackageStartupMessages(library(Biostrings))

    #read abundance file, first column is ASV_ID
    table <- read.table(file = '$table', sep = '\t', comment.char = '', header=TRUE)
    colnames(table)[1] <- "ASV_ID"

    #read fasta file of ASV sequences
    seq <- readDNAStringSet("$fasta")
    seq <- data.frame(ID=names(seq), sequence=paste(seq))

    # filter
    filtered_table = table[table\$ASV_ID %in% seq\$ID,]
    filtered_seq = seq[seq\$ID %in% table\$ASV_ID,]

    # sort
    filtered_table <- filtered_table[order(filtered_table\$ASV_ID),]
    filtered_seq <- filtered_seq[order(filtered_seq\$ID),]

    #write
    write.table(filtered_table, file = "filtered_abundances.tsv", row.names=FALSE, sep="\t", col.names = TRUE, quote = FALSE, na = '')
    write.table(data.frame(s = sprintf(">%s\n%s", filtered_seq\$ID, filtered_seq\$sequence)), file = 'filtered_sequences.fasta', col.names = FALSE, row.names = FALSE, quote = FALSE, na = '')

    # log retained features
    stats <- data.frame(
        file = list("$table","$fasta"),
        input = list( nrow(table), nrow(seq) ),
        output = list( nrow(filtered_table), nrow(filtered_seq) ) )
    write.table(stats, file = "stats_features.tsv", row.names=FALSE, sep="\t")

    # log retained counts
    stats <- as.data.frame( t( rbind( colSums(table[-1]), colSums(filtered_table[-1]) ) ) )
    stats\$ID <- rownames(stats)
    colnames(stats) <- c("input","output", "sample")
    write.table(stats, file = "stats_counts.tsv", row.names=FALSE, sep="\t")

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    Biostrings: ", packageVersion("Biostrings")) ), "versions.yml")
    """
}
