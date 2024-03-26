process FORMAT_TAXONOMY {
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path(database)
    val(suffix)

    output:
    path( "*assignTaxonomy.fna*" ), emit: assigntax
    path( "*addSpecies.fna*")     , emit: addspecies
    path( "ref_taxonomy.*.txt")   , emit: ref_tax_info
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ${params.dada_ref_databases[params.dada_ref_taxonomy]["fmtscript"]} \\

    #Giving out information
    echo -e "--dada_ref_taxonomy: ${params.dada_ref_taxonomy}\\n" >ref_taxonomy.${suffix}.txt
    echo -e "Title: ${params.dada_ref_databases[params.dada_ref_taxonomy]["title"]}\\n" >>ref_taxonomy.${suffix}.txt
    echo -e "Citation: ${params.dada_ref_databases[params.dada_ref_taxonomy]["citation"]}\\n" >>ref_taxonomy.${suffix}.txt
    echo "All entries: ${params.dada_ref_databases[params.dada_ref_taxonomy]}" >>ref_taxonomy.${suffix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | sed -n 1p | grep -Eo 'version [[:alnum:].]+' | sed 's/version //')
    END_VERSIONS
    """
}
