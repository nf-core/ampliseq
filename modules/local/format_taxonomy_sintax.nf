process FORMAT_TAXONOMY_SINTAX {
    label 'process_single'

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

    script:
    if (params.sintax_ref_tax_custom) {
        """
        set -- \$(head -c2 "${database}" | od -An -t u1)
        if [ "\$#" -ge 2 ] && [ "\$1" = "31" ] && [ "\$2" = "139" ]; then
            cp -fL "${database}" sintaxdb.fa.gz
        else
            gzip -c "${database}" > sintaxdb.fa.gz
        fi
        echo -e "--sintax_ref_tax_custom: ${params.sintax_ref_tax_custom}\\n" >ref_taxonomy_sintax.txt
        echo -e "Title: User-supplied SINTAX reference\\n" >>ref_taxonomy_sintax.txt
        echo -e "Citation: Not specified\\n" >>ref_taxonomy_sintax.txt
        echo -e "dbversion label: user_supplied\\n" >>ref_taxonomy_sintax.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
        END_VERSIONS
        """
    } else {
        """
        ${params.sintax_ref_databases[params.sintax_ref_taxonomy]["fmtscript"]}

        echo -e "--sintax_ref_taxonomy: ${params.sintax_ref_taxonomy}\\n" >ref_taxonomy_sintax.txt
        echo -e "Title: ${params.sintax_ref_databases[params.sintax_ref_taxonomy]["title"]}\\n" >>ref_taxonomy_sintax.txt
        echo -e "Citation: ${params.sintax_ref_databases[params.sintax_ref_taxonomy]["citation"]}\\n" >>ref_taxonomy_sintax.txt
        echo "All entries: ${params.sintax_ref_databases[params.sintax_ref_taxonomy]}" >>ref_taxonomy_sintax.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
        END_VERSIONS
        """
    }
}
