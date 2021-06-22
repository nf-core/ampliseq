// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process METADATA_PAIRWISE {
    tag "$metadata"
    label 'process_low'

    conda (params.enable_conda ? "bioconductor-dada2=1.18.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.18.0--r40h5f743cb_0"
    } else {
        container "quay.io/biocontainers/bioconductor-dada2:1.18.0--r40h5f743cb_0"
    }

    input:
    path(metadata)

    output:
    stdout

    script:
    """
    metadata_pairwise.r ${metadata}
    """
}