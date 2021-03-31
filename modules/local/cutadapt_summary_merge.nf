// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process CUTADAPT_SUMMARY_MERGE {
    tag "${files}"
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
    val(action)
    path(files)

	output:
	path("cutadapt_summary.tsv") , emit: tsv

    script:
    if (action == "merge") {
        """
        #!/usr/bin/env Rscript
        standard <- read.table(\"${files[0]}\", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        doubleprimer <- read.table(\"${files[1]}\", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        colnames(doubleprimer) <- c("sample", "cutadapt_doubleprimer_total_processed", "cutadapt_doubleprimer_reverse_complemented", "cutadapt_doubleprimer_passing_filters", "cutadapt_doubleprimer_passing_filters_percent")

        #merge
        df <- merge(standard, doubleprimer, by = "sample")

        #filter columns
        remove_columns <- c("cutadapt_doubleprimer_total_processed")
        for(column in remove_columns) df[column]<-NULL

        #write
        write.table(df, file = \"cutadapt_summary.tsv\", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
        """
    } else {
        """
        cp $files cutadapt_summary.tsv
        """
    }
}
