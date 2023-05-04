process CUTADAPT_SUMMARY {
    tag "${name}"
    label 'process_low'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    val(name)
    tuple val(meta), path(logs)

    output:
    path("*_summary.tsv") , emit: tsv
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode  = meta.single_end ? "single_end" : "paired_end"
    """
    cutadapt_summary.py $mode *.cutadapt.log > ${name}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version )
    END_VERSIONS
    """
}
