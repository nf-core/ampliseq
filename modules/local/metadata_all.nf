process METADATA_ALL {
    tag "$metadata"
    label 'process_low'

    conda (params.enable_conda ? "bioconductor-dada2=1.22.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

    input:
    path(metadata)
    val(metadata_category)

    output:
    stdout

    script:
    if( !metadata_category ) {
        """
        metadata_all.r ${metadata}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            R: \$(R --version 2>&1 | sed -n 1p | sed 's/R version //' | sed 's/ (.*//')
        END_VERSIONS
        """
    } else {
        """
        printf ${metadata_category}
        """
    }
}
