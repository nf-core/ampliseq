// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process FORMAT_TAXONOMY_QIIME {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    path(database)

    output:
    path( "*.tax" ), emit: tax
    path( "*.fna" ), emit: fasta
    path( "ref_taxonomy.txt"), emit: ref_tax_info

    script:
    """
    ${params.qiime_ref_databases[params.qiime_ref_taxonomy]["fmtscript"]}

    #Giving out information
    echo -e "--qiime_ref_taxonomy: ${params.qiime_ref_taxonomy}\n" >ref_taxonomy.txt
    echo -e "Title: ${params.dada_ref_databases[params.qiime_ref_taxonomy]["title"]}\n" >>ref_taxonomy.txt
    echo -e "Citation: ${params.dada_ref_databases[params.qiime_ref_taxonomy]["citation"]}\n" >>ref_taxonomy.txt
    echo "All entries: ${params.dada_ref_databases[params.qiime_ref_taxonomy]}" >>ref_taxonomy.txt
    """
}
