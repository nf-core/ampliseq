
process FORMAT_TAXONOMY_QIIME {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    path(database)

    output:
    path( "*.tax" ), emit: tax
    path( "*.fna" ), emit: fasta
    path( "ref_taxonomy.txt"), emit: ref_tax_info

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ${params.qiime_ref_databases[params.qiime_ref_taxonomy]["fmtscript"]}

    #Giving out information
    echo -e "--qiime_ref_taxonomy: ${params.qiime_ref_taxonomy}\n" >ref_taxonomy.txt
    echo -e "Title: ${params.qiime_ref_databases[params.qiime_ref_taxonomy]["title"]}\n" >>ref_taxonomy.txt
    echo -e "Citation: ${params.qiime_ref_databases[params.qiime_ref_taxonomy]["citation"]}\n" >>ref_taxonomy.txt
    echo "All entries: ${params.qiime_ref_databases[params.qiime_ref_taxonomy]}" >>ref_taxonomy.txt
    """
}
