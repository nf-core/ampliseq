process FORMAT_TAXONOMY_QIIME {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    path(database)
    
    output:
    path( "*.tax" ), emit: tax
    path( "*.fna" ), emit: fasta

    script:
    """
    ${params.qiime_ref_databases[params.qiime_ref_taxonomy]["fmtscript"]}
    """
}
