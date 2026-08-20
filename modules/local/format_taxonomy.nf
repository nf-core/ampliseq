process FORMAT_TAXONOMY {
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(db_key), path(database)

    output:
    tuple val(db_key), path( "*assignTaxonomy.fna*" ), emit: assigntax
    tuple val(db_key), path( "*addSpecies.fna*")     , emit: addspecies
    path( "ref_taxonomy.*.txt")                      , emit: ref_tax_info
    path "versions.yml"                              , emit: versions_format_taxonomy, topic: versions

    script:
    def suffix = db_key.replace('=','_').replace('.','_')
    """
    ${params.dada_ref_databases[db_key]["fmtscript"]} \\

    #Giving out information
    echo -e "--dada_ref_taxonomy: ${db_key}\\n" >ref_taxonomy.${suffix}.txt
    echo -e "Title: ${params.dada_ref_databases[db_key]["title"]}\\n" >>ref_taxonomy.${suffix}.txt
    echo -e "Citation: ${params.dada_ref_databases[db_key]["citation"]}\\n" >>ref_taxonomy.${suffix}.txt
    echo "All entries: ${params.dada_ref_databases[db_key]}" >>ref_taxonomy.${suffix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \"\$BASH_VERSION\")
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """
}
