// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process QIIME2_ANCOM_ASV {
    tag "${table.baseName}"
    label 'process_medium'
    label 'single_cpu'
    label 'process_long'
    label 'error_ignore'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.4"

    input:
    tuple path(metadata), path(table)

    output:
    path("ancom/*")     , emit: ancom
    path "*.version.txt", emit: version

    script:
    def software     = getSoftwareName(task.process)
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    qiime composition add-pseudocount \
        --i-table ${table} \
        --o-composition-table comp-${table}
    qiime composition ancom \
        --i-table comp-${table} \
        --m-metadata-file ${metadata} \
        --m-metadata-column ${table.baseName} \
        --o-visualization comp-${table.baseName}.qzv
    qiime tools export --input-path comp-${table.baseName}.qzv \
        --output-path ancom/Category-${table.baseName}-ASV

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}
