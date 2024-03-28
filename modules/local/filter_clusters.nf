process FILTER_CLUSTERS {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::pandas=1.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5':
        'biocontainers/pandas:1.1.5' }"

    input:
    tuple val(meta), path(clusters)
    path(asv)

    output:
    path( "ASV_post_clustering_filtered.table.tsv") , emit: asv
    path( "ASV_post_clustering_filtered.fna"      ) , emit: fasta
    path( "ASV_post_clustering_filtered.stats.tsv") , emit: stats
    path( "versions.yml"                          ) , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix   = task.ext.prefix ?: "'$meta.id'"
    def clusters = "'$clusters'"
    """
    ulimit -s unlimited
    echo ${clusters} | filt_clusters.py -t ${asv} -p ${prefix} -c -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
