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

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if [[ $workflow.manifest.version == *dev ]]; then
        ampliseq_version="v$workflow.manifest.version, revision ${workflow.scriptId.substring(0,10)}"
    else
        ampliseq_version="v$workflow.manifest.version"
    fi
    echo "This is ampliseq version: \$ampliseq_version" > testar_jt.txt

    sbdiexportreannotate.R \"${params.dada_ref_databases[params.dada_ref_taxonomy]["dbversion"]}\" $taxonomytable \"\$ampliseq_version\"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1 | sed -n 1p | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """
}
