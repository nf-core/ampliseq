process FORMAT_TAXONOMY {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::p7zip=15.09" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/p7zip:15.09--h2d50403_4' :
        'quay.io/biocontainers/p7zip:15.09--h2d50403_4' }"

    input:
    path(database)

    output:
    path( "*assignTaxonomy.fna*" ), emit: assigntax
    path( "*addSpecies.fna*")     , emit: addspecies
    path( "ref_taxonomy.txt")     , emit: ref_tax_info
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ${params.dada_ref_databases[params.dada_ref_taxonomy]["fmtscript"]}

    #Giving out information
    echo -e "--dada_ref_taxonomy: ${params.dada_ref_taxonomy}\n" >ref_taxonomy.txt
    echo -e "Title: ${params.dada_ref_databases[params.dada_ref_taxonomy]["title"]}\n" >>ref_taxonomy.txt
    echo -e "Citation: ${params.dada_ref_databases[params.dada_ref_taxonomy]["citation"]}\n" >>ref_taxonomy.txt
    echo "All entries: ${params.dada_ref_databases[params.dada_ref_taxonomy]}" >>ref_taxonomy.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        7za: \$(echo \$(7za --help) | sed 's/.*p7zip Version //; s/(.*//')
    END_VERSIONS
    """
}
