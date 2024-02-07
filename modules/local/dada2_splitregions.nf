process DADA2_SPLITREGIONS {
    label 'process_low'

    // TODO: DADA2 not neccessary, R base sufficient!
    conda "bioconda::bioconductor-dada2=1.22.0 conda-forge::r-digest=0.6.30"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

    input:
    tuple val(meta), val(mapping)
    path(table)

    output:
    tuple val(meta), path( "DADA2_table_*.tsv" ), emit: dada2asv
    tuple val(meta), path( "ASV_table_*.tsv" )  , emit: asv
    tuple val(meta), path( "ASV_seqs_*.fasta" ) , emit: fasta
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //make groovy map to R list
    def mapping_r_list = mapping
        .toString()
        .replaceAll("':","=")
        .replaceAll(", '",",")
        .replaceAll("\\['","list(")
        .replaceAll("\\[","list(")
        .replaceAll("\\]",")")
        .replaceAll("false","'false'")
    def suffix = "region" + meta.region + "_" + meta.fw_primer + "_" + meta.rv_primer
    """
    #!/usr/bin/env Rscript

    nested_list <- $mapping_r_list
    mapping <- as.data.frame(do.call(rbind, nested_list))

    df <- read.csv("$table", header=TRUE, sep="\\t")

    # extract samples of this region
    keep <- intersect( colnames(df), c("ASV_ID", unlist(mapping\$id), "sequence" ) )
    df <- subset( df, select = keep )

    # modify sample names from .id to .sample
    for (row in 1:nrow(mapping)) {
        colnames(df) <- gsub( mapping[row,]\$id, mapping[row,]\$sample, colnames(df) )
    }

    # filter rows with only 0, occurring because many samples were removed
    df <- df[as.logical(rowSums(df[,2:(ncol(df)-1)] != 0)), ]

    # Write file with ASV abdundance and sequences to file
    write.table(df, file = "DADA2_table_${suffix}.tsv", sep = "\\t", row.names = FALSE, quote = FALSE, na = '')

    # Write fasta file with ASV sequences to file
    write.table(data.frame(s = sprintf(">%s\n%s", df\$ASV_ID, df\$sequence)), 'ASV_seqs_${suffix}.fasta', col.names = FALSE, row.names = FALSE, quote = FALSE, na = '')

    # Write ASV file with ASV abundances to file
    df\$sequence <- NULL
    write.table(df, file = "ASV_table_${suffix}.tsv", sep="\\t", row.names = FALSE, quote = FALSE, na = '')


    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    """
}
