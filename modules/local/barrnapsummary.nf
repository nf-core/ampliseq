process BARRNAPSUMMARY {
    label 'process_single'

    conda (params.enable_conda ? "python=3.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'quay.io/biocontainers/python:3.9' }"

    input:
    path predictions

    output:
    path "summary.tsv" , emit: summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    summarize_barrnap.py $predictions

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version )
    END_VERSIONS
    """
}
