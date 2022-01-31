process DADA2_QUALITY {
    tag "$meta"
    label 'process_low'

    conda (params.enable_conda ? "bioconductor-dada2=1.22.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    path "${meta}_qual_stats.pdf"            , emit: pdf
    tuple val(meta), path("*_qual_stats.tsv"), emit: tsv
    path "*.args.txt"                        , emit: args

    script:
    def args = task.ext.args ?: ''
    """
    dada_quality.r "${meta}_qual_stats" $args
    echo 'plotQualityProfile\t$args' > "plotQualityProfile.args.txt"
    """
}
