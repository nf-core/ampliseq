process QIIME2_UNPACK {
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path(database)

    output:
    path("*.fna"), emit: fasta
    path("*.tax"), emit: tax
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // TODO: need to not have this be a copy.
    script:
    """
    cp $database/*.fna .
    cp $database/*.tax .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | sed -n 1p | sed 's/GNU bash, version //g')
    END_VERSIONS
    """
}
