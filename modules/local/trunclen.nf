process TRUNCLEN {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    tuple val(meta), path(qual_stats)

    output:
    tuple val(meta), stdout

    script:
    def args = task.ext.args ?: ''
    """
    trunclen.py $qual_stats $args
    """
}
