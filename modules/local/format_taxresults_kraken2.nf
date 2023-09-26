process FORMAT_TAXRESULTS_KRAKEN2 {
    label 'process_low'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(report)
    tuple val(meta), path(classified_reads_assignment)
    val(taxlevels_input)

    output:
    path("*.kraken2.keys.tsv")       , emit: keys_tsv
    path("*.kraken2.complete.tsv")   , emit: complete_tsv
    path("*.kraken2.tsv")            , emit: tsv
    path("*.kraken2.into-qiime2.tsv"), emit: qiime2_tsv
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def taxlevels = taxlevels_input ? taxlevels_input : "D,P,C,O,F,G,S"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    ## report
    report <- read.table(file = "$report", header = FALSE, sep = "\t")
    colnames(report) <- c("percent","fragments_total","fragments","rank","taxid","taxon")

    taxonomy_standard_labels <- unlist(strsplit("$taxlevels",","))

    # prepare empty result data frame
    headers <- c("rank","taxid","taxon","taxonomy_rank","taxonomy_standard")
    headers <- c(headers, taxonomy_standard_labels)
    df = data.frame(matrix(nrow = 0, ncol = length(headers)))
    colnames(df) = headers

    # go through each line of report and aggregate information
    taxonomy_rank <- c("root")
    for (i in 1:nrow(report)) {
        tax <- report\$taxon[i]
        taxnospaces <- gsub("^ *", "", tax)

        # (1) extract full taxonomy with taxa ranks
        nspaces = nchar(tax)- nchar(taxnospaces)
        indent = nspaces/2 +1 #+1 is to catch also entries without any indent, e.g. "root"
        # if this is a lower taxonomic rank, add it, otherwise reset and add
        if ( indent > length(taxonomy_rank) ) {
            taxonomy_rank <- c(taxonomy_rank,paste0(report\$rank[i],"__",taxnospaces))
        } else {
            taxonomy_rank <- taxonomy_rank[1:indent-1]
            taxonomy_rank <- c(taxonomy_rank,paste0(report\$rank[i],"__",taxnospaces))
        }
        taxonomy_rank_string <- paste(taxonomy_rank,collapse=";")

        # (2) filter taxonomy_rank to only contain entries with taxonomy_standard_labels+"__"
        taxonomy_standard = list()
        taxonomy_ranks = gsub( "__.*", "", taxonomy_rank )
        for (x in taxonomy_standard_labels) {
            taxonomy_ranks_match <- match(x, taxonomy_ranks)
            if ( !is.na(taxonomy_ranks_match) ) {
                taxonomy_clean <- gsub( ".*__", "", taxonomy_rank[taxonomy_ranks_match] )
                taxonomy_standard <- c(taxonomy_standard,taxonomy_clean)
            } else {
                taxonomy_standard <- c(taxonomy_standard,"")
            }
        }
        taxonomy_standard_string <- paste(taxonomy_standard,collapse=";")
        names(taxonomy_standard) <- taxonomy_standard_labels

        # (3) propagate all in results data frame
        results <- c(
            rank=report\$rank[i],
            taxid=report\$taxid[i],
            taxon=taxnospaces,
            taxonomy_rank=taxonomy_rank_string,
            taxonomy_standard=taxonomy_standard_string,
            taxonomy_standard)
        df <- rbind(df, results)
    }
    df\$taxon_taxid <- paste0(df\$taxon," (taxid ",df\$taxid,")")
    write.table(df, file = "${prefix}.kraken2.keys.tsv", row.names=FALSE, quote=FALSE, sep="\t")

    # merge with reads
    creads <- read.table(file = "$classified_reads_assignment", header = FALSE, sep = "\t")
    colnames(creads) <- c("classified","ASV_ID","taxonomy","ASV_length","kmer_LCA")
    if ( !all( creads\$taxonomy %in% df\$taxon_taxid ) ) { stop(paste("$classified_reads_assignment","and","$report","dont share all IDs, exit"), call.=FALSE) }
    merged <- merge(creads, df, by.x="taxonomy", by.y="taxon_taxid", all.x=TRUE, all.y=FALSE)
    write.table(merged, file = "${prefix}.kraken2.complete.tsv", row.names=FALSE, quote=FALSE, sep="\t")

    # get downstream compatible table
    merged_reduced <- subset(merged, select = c("ASV_ID",taxonomy_standard_labels,"taxonomy"))
    colnames(merged_reduced) <- c("ASV_ID",taxonomy_standard_labels,"lowest_match")
    write.table(merged_reduced, file = "${prefix}.kraken2.tsv", row.names=FALSE, quote=FALSE, sep="\t")

    # get QIIME2 downstream compatible table
    qiime <- merged[c("ASV_ID", "taxonomy_standard")]
    colnames(qiime) <- c("ASV_ID", "taxonomy")
    write.table(qiime, file = "${prefix}.kraken2.into-qiime2.tsv", row.names=FALSE, quote=FALSE, sep="\t")

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")) ), "versions.yml")
    """
}
