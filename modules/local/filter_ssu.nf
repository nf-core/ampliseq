process FILTER_SSU {
    tag "${fasta}"
    label 'process_low'

    conda "bioconductor::biostrings=2.58.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-biostrings:2.58.0--r40h037d062_0' :
        'quay.io/biocontainers/bioconductor-biostrings:2.58.0--r40h037d062_0' }"

    input:
    path(fasta)
    path(table)
    path(barrnap_summary)

    output:
    path( "stats.ssu.tsv" )      , emit: stats
    path( "ASV_table.ssu.tsv" )  , emit: asv
    path( "ASV_seqs.ssu.fasta" ) , emit: fasta
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def kingdom = params.filter_ssu ?: "bac,arc,mito,euk"
    """
    #!/usr/bin/env Rscript

    #load packages
    suppressPackageStartupMessages(library(Biostrings))

    kingdom <- as.list(strsplit("$kingdom", ",")[[1]])

    df = read.table("$barrnap_summary", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    # keep only ASV_ID & eval columns & sort
    df <- subset(df, select = c(ASV_ID,mito_eval,euk_eval,arc_eval,bac_eval))

    # choose kingdom (column) with lowest evalue
    df[is.na(df)] <- 1
    df\$result = colnames(df[,2:5])[apply(df[,2:5],1,which.min)]
    df\$result = gsub("_eval", "", df\$result)

    # filter ASVs
    df_filtered = subset(df, df\$result %in% kingdom)
    id_filtered = subset(df_filtered, select = c(ASV_ID))

    #error if all ASVs are removed
    if ( nrow(df_filtered) == 0 ) stop("Chosen kingdom(s) by --filter_ssu has no matches. Please choose a different kingdom (domain) or omit filtering.")

    #read abundance file, first column is ASV_ID
    table <- read.table(file = "$table", sep = '\t', comment.char = "", header=TRUE)
    colnames(table)[1] <- "ASV_ID"

    #read fasta file of ASV sequences
    seq <- readDNAStringSet("$fasta")
    seq <- data.frame(ID=names(seq), sequence=paste(seq))

    #check if all ids match
    if(!all(id_filtered\$ID %in% seq\$ID))  {stop(paste(paste(files,sep=","),"and","$fasta","dont share all IDs, exit."), call.=FALSE)}
    if(!all(id_filtered\$ID %in% table\$ASV_ID))  {stop(paste(paste(files,sep=","),"and","$table","dont share all IDs, exit"), call.=FALSE)}

    #merge
    filtered_table <- merge(table, id_filtered, by.x="ASV_ID", by.y="ASV_ID", all.x=FALSE, all.y=TRUE)
    filtered_seq <- merge(seq, id_filtered, by.x="ID", by.y="ASV_ID", all.x=FALSE, all.y=TRUE)

    #write
    write.table(filtered_table, file = "ASV_table.ssu.tsv", row.names=FALSE, sep="\t", col.names = TRUE, quote = FALSE, na = '')
    write.table(data.frame(s = sprintf(">%s\n%s", filtered_seq\$ID, filtered_seq\$sequence)), 'ASV_seqs.ssu.fasta', col.names = FALSE, row.names = FALSE, quote = FALSE, na = '')

    #stats
    stats <- as.data.frame( t( rbind( colSums(table[-1]), colSums(filtered_table[-1]) ) ) )
    stats\$ID <- rownames(stats)
    colnames(stats) <- c("ssufilter_input","ssufilter_output", "sample")
    write.table(stats, file = "stats.ssu.tsv", row.names=FALSE, sep="\t")

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    Biostrings: ", packageVersion("Biostrings")) ), "versions.yml")
    """
}
