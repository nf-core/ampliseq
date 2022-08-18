process NOVASEQ_ERR {
    tag "$meta.run"
    label 'process_medium'

    conda (params.enable_conda ? "bioconductor-dada2=1.22.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"
    input:
    tuple val(meta), path(errormodel)

    output:
    tuple val(meta), path("*.md.err.rds"), emit: errormodel
    tuple val(meta), path("*.md.err.pdf"), emit: pdf
    tuple val(meta), path("*.md.err.convergence.txt"), emit: convergence
    //path "versions.yml"                  , emit: versions

    script:
    if (!meta.single_end) {
        """
        novaseq_err_pe.r ${errormodel[0]} ${errormodel[1]} ${meta.run}
        """
    } else {
        """
        novaseq_err_se.r ${errormodel} ${meta.run}
        """
    }
}
