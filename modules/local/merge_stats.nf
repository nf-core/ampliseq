process MERGE_STATS {
    label 'process_single'

    conda "bioconda::bioconductor-dada2=1.34.0 conda-forge::r-base=4.4.3 conda-forge::tbb=2020.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.34.0--r44he5774e6_2' :
        'biocontainers/bioconductor-dada2:1.34.0--r44he5774e6_2' }"

    input:
    path('file1.tsv')
    path('file2.tsv')

    output:
    path("overall_summary.tsv") , emit: tsv
    path "versions.yml"         , emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    x <- read.table(\"file1.tsv\", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    y <- read.table(\"file2.tsv\", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

    #merge
    df <- merge(x, y, by = "sample", all = TRUE)

    #write
    write.table(df, file = \"overall_summary.tsv\", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")) ), "versions.yml")
    """
}
