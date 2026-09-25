process FORMAT_TAXONOMY_SIDLE {
    label 'process_single'

    conda "conda-forge::bash=5.2.37 conda-forge::coreutils=9.12 conda-forge::findutils=4.10.0 conda-forge::gawk=5.4.1 conda-forge::grep=3.12 conda-forge::gzip=1.14 conda-forge::sed=4.10 conda-forge::tar=1.35 conda-forge::unzip=6.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/739a65cef3554bd727c1a6467e12d7bda4b52b94c035fcb6fa333e15c3ecb8aa/data' :
        'community.wave.seqera.io/library/bash_coreutils_findutils_gawk_pruned:63fbf75c749fd019' }"

    input:
    path('database.tar.gz')
    val(suffix)

    output:
    path( "*.seq.fasta" )         , emit: seq
    path( "*.alnseq.fasta")       , emit: alnseq
    path( "*.tax.txt")            , emit: tax
    path( "ref_taxonomy.*.txt")   , emit: ref_tax_info
    path "versions.yml"           , emit: versions_format_taxonomy_sidle, topic: versions

    script:
    def derep = params.sidle_ref_databases[params.sidle_ref_taxonomy]["derep"] ?: "99"
    """
    ${params.sidle_ref_databases[params.sidle_ref_taxonomy]["fmtscript"]} ${derep} \\

    #Giving out information
    echo -e "--sidle_ref_taxonomy: ${params.sidle_ref_taxonomy}\\n" >ref_taxonomy.${suffix}.txt
    echo -e "Title: ${params.sidle_ref_databases[params.sidle_ref_taxonomy]["title"]}\\n" >>ref_taxonomy.${suffix}.txt
    echo -e "Citation: ${params.sidle_ref_databases[params.sidle_ref_taxonomy]["citation"]}\\n" >>ref_taxonomy.${suffix}.txt
    echo "All entries: ${params.sidle_ref_databases[params.sidle_ref_taxonomy]}" >>ref_taxonomy.${suffix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \"\$BASH_VERSION\")
        sed: "\$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')"
    END_VERSIONS
    """
}
