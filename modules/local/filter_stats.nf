process FILTER_STATS {
    label 'process_low'

    conda (params.enable_conda ? "pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    path(unfiltered)
    path(filtered)

    output:
    path("count_table_filter_stats.tsv"), emit: tsv

    script:
    """
    filter_stats.py $unfiltered $filtered
    """
}
