// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]

process METADATA_ALL {
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
    val(metadata_category)

    output:
    stdout

    script:
    if( !metadata_category ) {
        """
        metadata_all.r ${metadata}
        """
    } else {
        """
        printf ${metadata_category}
        """
    }
}
