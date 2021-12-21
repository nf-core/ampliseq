process MERGE_STATS {
    label 'process_low'

    conda (params.enable_conda ? "bioconductor-dada2=1.20.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.20.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-dada2:1.20.0--r41h399db7b_0' }"

    input:
    path('file1.tsv')
    path('file2.tsv')

    output:
    path("overall_summary.tsv") , emit: tsv

    script:
    """
    #!/usr/bin/env Rscript
    x <- read.table(\"file1.tsv\", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    y <- read.table(\"file2.tsv\", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    #merge
    df <- merge(x, y, by = "sample", all = TRUE)

    #write
    write.table(df, file = \"overall_summary.tsv\", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
    """
}
