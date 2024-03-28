process BARRNAPSUMMARY {
    label 'process_single'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    path predictions

    output:
    path "summary.tsv" , emit: summary
    path "*warning.txt", emit: warning
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    summarize_barrnap.py $predictions

    if [[ \$(wc -l < summary.tsv ) -le 1 ]]; then
        touch WARNING_no_rRNA_found_warning.txt
    else
        touch no_warning.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version )
    END_VERSIONS
    """
}
