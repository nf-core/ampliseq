process METADATA_PAIRWISE {
    tag "$metadata"
    label 'process_single'

    conda "bioconda::bioconductor-dada2=1.34.0 conda-forge::r-base=4.4.3 conda-forge::tbb=2020.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.34.0--r44he5774e6_2' :
        'biocontainers/bioconductor-dada2:1.34.0--r44he5774e6_2' }"

    input:
    path(metadata)

    output:
    stdout               emit: category
    path "versions.yml", emit: versions

    script:
    """
    metadata_pairwise.r ${metadata}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1 | sed -n 1p | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """
}
