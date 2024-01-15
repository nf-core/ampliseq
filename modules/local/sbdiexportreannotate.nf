process SBDIEXPORTREANNOTATE {
    tag "${taxonomytable}"
    label 'process_low'

    conda "conda-forge::r-tidyverse=1.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse:1.2.1' :
        'biocontainers/r-tidyverse:1.2.1' }"

    input:
    path taxonomytable
    val  taxonomymethod
    val  dbversion
    val  cut_its
    path predictions

    output:
    path "*.tsv"       , emit: sbdiannottables
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if [[ $workflow.manifest.version == *dev ]]; then
        ampliseq_version="v$workflow.manifest.version, revision: ${workflow.scriptId.substring(0,10)}"
    else
        ampliseq_version="v$workflow.manifest.version"
    fi

    sbdiexportreannotate.R \"$dbversion\" $taxonomytable $taxonomymethod \"\$ampliseq_version\" $cut_its $predictions

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1 | sed -n 1p | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """
}
