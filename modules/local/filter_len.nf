process FILTER_LEN {
    tag "${fasta}"
    label 'process_low'

    conda "bioconda::bioconductor-biostrings=2.58.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-biostrings:2.58.0--r40h037d062_0' :
        'biocontainers/bioconductor-biostrings:2.58.0--r40h037d062_0' }"

    input:
    path(fasta)
    path(table)

    output:
    path( "stats.len.tsv" )      , emit: stats
    path( "ASV_table.len.tsv" )  , emit: asv, optional: true
    path( "ASV_seqs.len.fasta" ) , emit: fasta
    path( "ASV_len_orig.tsv" )   , emit: len_orig
    path( "ASV_len_filt.tsv" )   , emit: len_filt
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def min_len_asv = task.ext.min_len_asv ?: '1'
    def max_len_asv = task.ext.max_len_asv ?: '1000000'

    def read_table  = table ? "table <- read.table(file = '$table', sep = '\t', comment.char = '', header=TRUE)" : "table <- data.frame(matrix(ncol = 1, nrow = 0))"
    def asv_table_filtered  = table ? "ASV_table.len.tsv" : "empty_ASV_table.len.tsv"
    """
    #!/usr/bin/env Rscript

    #load packages
    suppressPackageStartupMessages(library(Biostrings))

    #read abundance file, first column is ASV_ID
    $read_table
    colnames(table)[1] <- "ASV_ID"

    #read fasta file of ASV sequences
    seq <- readDNAStringSet("$fasta")
    seq <- data.frame(ID=names(seq), sequence=paste(seq))

    #filter
    filtered_seq <- seq[nchar(seq\$sequence) %in% $min_len_asv:$max_len_asv,]
    list <- filtered_seq[, "ID", drop = FALSE]
    filtered_table <- merge(table, list, by.x="ASV_ID", by.y="ID", all.x=FALSE, all.y=TRUE)

    #report
    distribution_before <- table(nchar(seq\$sequence))
    distribution_before <- data.frame(Length=names(distribution_before),Counts=as.vector(distribution_before))
    distribution_after <- table(nchar(filtered_seq\$sequence))
    distribution_after <- data.frame(Length=names(distribution_after),Counts=as.vector(distribution_after))

    #write
    write.table(filtered_table, file = "$asv_table_filtered", row.names=FALSE, sep="\t", col.names = TRUE, quote = FALSE, na = '')
    write.table(data.frame(s = sprintf(">%s\n%s", filtered_seq\$ID, filtered_seq\$sequence)), 'ASV_seqs.len.fasta', col.names = FALSE, row.names = FALSE, quote = FALSE, na = '')
    write.table(distribution_before, file = "ASV_len_orig.tsv", row.names=FALSE, sep="\t", col.names = TRUE, quote = FALSE, na = '')
    write.table(distribution_after, file = "ASV_len_filt.tsv", row.names=FALSE, sep="\t", col.names = TRUE, quote = FALSE, na = '')

    #stats
    stats <- as.data.frame( t( rbind( colSums(table[-1]), colSums(filtered_table[-1]) ) ) )
    stats\$ID <- rownames(stats)
    colnames(stats) <- c("lenfilter_input","lenfilter_output", "sample")
    write.table(stats, file = "stats.len.tsv", row.names=FALSE, sep="\t")

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    Biostrings: ", packageVersion("Biostrings")) ), "versions.yml")
    """
}
