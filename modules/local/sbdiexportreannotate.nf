process SBDIEXPORTREANNOTATE {
    tag "${taxonomytable}"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::r-tidyverse=1.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse:1.2.1' :
        'quay.io/biocontainers/r-tidyverse:1.2.1' }"

    input:
    path taxonomytable

    output:
    path "*.tsv"       , emit: sbdiannottables
    path "versions.yml", emit: versions

    script:
    """
    sbdiexportreannotate.R ${params.dada_ref_databases[params.dada_ref_taxonomy]["dbversion"]} $taxonomytable

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1 | sed -n 1p | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """
}
