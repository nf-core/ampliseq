
process FORMAT_TAXONOMY_QIIME {
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path(database)

    output:
    path( "*.tax" )          , emit: tax
    path( "*.fna" )          , emit: fasta
    path( "ref_taxonomy.txt"), emit: ref_tax_info
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ${params.qiime_ref_databases[params.qiime_ref_taxonomy]["fmtscript"]} \\

    #Giving out information
    echo -e "--qiime_ref_taxonomy: ${params.qiime_ref_taxonomy}\\n" >ref_taxonomy.txt
    echo -e "Title: ${params.qiime_ref_databases[params.qiime_ref_taxonomy]["title"]}\\n" >>ref_taxonomy.txt
    echo -e "Citation: ${params.qiime_ref_databases[params.qiime_ref_taxonomy]["citation"]}\\n" >>ref_taxonomy.txt
    echo "All entries: ${params.qiime_ref_databases[params.qiime_ref_taxonomy]}" >>ref_taxonomy.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | sed -n 1p | grep -Eo 'version [[:alnum:].]+' | sed 's/version //')
    END_VERSIONS
    """
}
