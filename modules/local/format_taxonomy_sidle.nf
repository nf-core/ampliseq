process FORMAT_TAXONOMY_SIDLE {
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path('database.tar.gz')
    val(suffix)

    output:
    path( "*.seq.fasta" )         , emit: seq
    path( "*.alnseq.fasta")       , emit: alnseq
    path( "*.tax.txt")            , emit: tax
    path( "ref_taxonomy.*.txt")   , emit: ref_tax_info
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

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
        bash: \$(bash --version | sed -n 1p | grep -Eo 'version [[:alnum:].]+' | sed 's/version //')
    END_VERSIONS
    """
}
