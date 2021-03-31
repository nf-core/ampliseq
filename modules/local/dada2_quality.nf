// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process DADA2_QUALITY {
    tag "$meta"
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
    tuple val(meta), path(reads)
    
    output:
    path "${meta}_qual_stats.pdf"            , emit: pdf
    tuple val(meta), path("*_qual_stats.tsv"), emit: tsv
    path "*.args.txt"                        , emit: args

    script:
    """
    dada_quality.r "${meta}_qual_stats" $options.args
    echo 'plotQualityProfile\t$options.args' > "plotQualityProfile.args.txt"
    """
}