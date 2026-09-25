process FORMAT_TAXONOMY_VSEARCH_LCA {
    label 'process_single'

    conda "conda-forge::bash=5.2.37 conda-forge::coreutils=9.12 conda-forge::findutils=4.10.0 conda-forge::gawk=5.4.1 conda-forge::grep=3.12 conda-forge::gzip=1.14 conda-forge::sed=4.10 conda-forge::tar=1.35 conda-forge::unzip=6.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/739a65cef3554bd727c1a6467e12d7bda4b52b94c035fcb6fa333e15c3ecb8aa/data' :
        'community.wave.seqera.io/library/bash_coreutils_findutils_gawk_pruned:63fbf75c749fd019' }"

    input:
    path(database)

    output:
    path( "sintaxdb.fa.gz" )             , emit: db
    path( "ref_taxonomy_vsearch_lca.txt"), emit: ref_tax_info
    path "versions.yml"                  , emit: versions_format_taxonomy_vsearch_lca, topic: versions

    script:
    if (params.vsearch_lca_ref_tax_custom) {
        """
        set -- \$(head -c2 "${database}" | od -An -t u1)
        if [ "\$#" -ge 2 ] && [ "\$1" = "31" ] && [ "\$2" = "139" ]; then
            cp -fL "${database}" sintaxdb.fa.gz
        else
            gzip -c "${database}" > sintaxdb.fa.gz
        fi
        echo -e "--vsearch_lca_ref_tax_custom: ${params.vsearch_lca_ref_tax_custom}\\n" >ref_taxonomy_vsearch_lca.txt
        echo -e "Title: User-supplied VSEARCH LCA reference\\n" >>ref_taxonomy_vsearch_lca.txt
        echo -e "Citation: Not specified\\n" >>ref_taxonomy_vsearch_lca.txt
        echo -e "dbversion label: user_supplied\\n" >>ref_taxonomy_vsearch_lca.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bash: \$(echo \"\$BASH_VERSION\")
        END_VERSIONS
        """
    } else {
        """
        ${params.vsearch_lca_ref_databases[params.vsearch_lca_ref_taxonomy]["fmtscript"]}

        echo -e "--vsearch_lca_ref_taxonomy: ${params.vsearch_lca_ref_taxonomy}\\n" >ref_taxonomy_vsearch_lca.txt
        echo -e "Title: ${params.vsearch_lca_ref_databases[params.vsearch_lca_ref_taxonomy]["title"]}\\n" >>ref_taxonomy_vsearch_lca.txt
        echo -e "Citation: ${params.vsearch_lca_ref_databases[params.vsearch_lca_ref_taxonomy]["citation"]}\\n" >>ref_taxonomy_vsearch_lca.txt
        echo "All entries: ${params.vsearch_lca_ref_databases[params.vsearch_lca_ref_taxonomy]}" >>ref_taxonomy_vsearch_lca.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bash: \$(echo \"\$BASH_VERSION\")
            sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
        END_VERSIONS
        """
    }
}
