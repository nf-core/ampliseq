process FORMAT_TAXONOMY_SINTAX {
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path(database)

    output:
    path( "sintaxdb.fa.gz" ), emit: db
    path( "ref_taxonomy_sintax.txt")     , emit: ref_tax_info
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ${params.sintax_ref_databases[params.sintax_ref_taxonomy]["fmtscript"]} \\

    #Giving out information
    echo -e "--sintax_ref_taxonomy: ${params.sintax_ref_taxonomy}\\n" >ref_taxonomy_sintax.txt
    echo -e "Title: ${params.sintax_ref_databases[params.sintax_ref_taxonomy]["title"]}\\n" >>ref_taxonomy_sintax.txt
    echo -e "Citation: ${params.sintax_ref_databases[params.sintax_ref_taxonomy]["citation"]}\\n" >>ref_taxonomy_sintax.txt
    echo "All entries: ${params.sintax_ref_databases[params.sintax_ref_taxonomy]}" >>ref_taxonomy_sintax.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | sed -n 1p | grep -Eo 'version [[:alnum:].]+' | sed 's/version //')
    END_VERSIONS
    """
}
