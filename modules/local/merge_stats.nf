// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process MERGE_STATS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconductor-dada2=1.18.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.18.0--r40h5f743cb_0"
    } else {
        container "quay.io/biocontainers/bioconductor-dada2:1.18.0--r40h5f743cb_0"
    }

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
