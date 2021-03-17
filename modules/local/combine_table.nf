// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process COMBINE_TABLE {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconductor::biostrings=2.58.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-biostrings:2.58.0--r40h037d062_0"
    } else {
        container "quay.io/biocontainers/bioconductor-biostrings:2.58.0--r40h037d062_0"
    }

	input:
	path(table)
	path(seq)
	path(tax)

	output:
	path("qiime2_ASV_table.tsv")

    script:
	"""
	combine_table.r ${table} ${seq} ${tax}
	"""
}