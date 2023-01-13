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
    path(lists)

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

    #use only selected kingdom
    dir.create("./selection")
    kingdom <- as.list(strsplit("$kingdom", ",")[[1]])
    for (x in kingdom) {
        file.copy(paste0(x,".matches.txt"), paste0("./selection/",x,".matches.txt"))
    }
    files = list.files(path = "./selection", pattern="*matches.txt", full.names = TRUE)

    #error if (all) file(s) is/are empty
    if ( all(file.size(files) == 0L) ) stop("Chosen kingdom(s) by --filter_ssu has no matches. Please choose a diffferent kingdom or omit filtering.")
    files = files[file.size(files) != 0L]

    #read positive ID lists
    list = do.call(rbind, lapply(files, function(x) read.csv(x, stringsAsFactors = FALSE, header = FALSE)))
    list = unique(list)
    colnames(list)[1] <- "ID"

    #read abundance file, first column is ASV_ID
    table <- read.table(file = "$table", sep = '\t', comment.char = "", header=TRUE)
    colnames(table)[1] <- "ASV_ID"

    #read fasta file of ASV sequences
    seq <- readDNAStringSet("$fasta")
    seq <- data.frame(ID=names(seq), sequence=paste(seq))

    #check if all ids match
    if(!all(list\$ID %in% seq\$ID))  {stop(paste(paste(files,sep=","),"and","$fasta","dont share all IDs, exit."), call.=FALSE)}
    if(!all(list\$ID %in% table\$ASV_ID))  {stop(paste(paste(files,sep=","),"and","$table","dont share all IDs, exit"), call.=FALSE)}

    #merge
    filtered_table <- merge(table, list, by.x="ASV_ID", by.y="ID", all.x=FALSE, all.y=TRUE)
    filtered_seq <- merge(seq, list, by.x="ID", by.y="ID", all.x=FALSE, all.y=TRUE)

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
