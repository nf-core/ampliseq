// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process CUTADAPT_SUMMARY {
    tag "${name}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

	input:
    val(name)
	tuple val(meta), path(logs)

	output:
	path("*_summary.tsv") , emit: tsv
    path "*.version.txt"  , emit: version

    script:
    def software = "python"
    def mode  = meta.single_end ? "single_end" : "paired_end"
	"""
	cutadapt_summary.py $mode *.cutadapt.log > ${name}_summary.tsv
    echo \$(python --version) > ${software}.version.txt
	"""
}