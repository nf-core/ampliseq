// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process QIIME2_EXPORT_RELASV {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    path(table)

    output:
    path("rel-table-ASV.tsv"), emit: tsv
    path "*.version.txt"     , emit: version

    script:
    def software     = getSoftwareName(task.process)
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    #convert to relative abundances
    qiime feature-table relative-frequency \
        --i-table ${table} \
        --o-relative-frequency-table relative-table-ASV.qza

    #export to biom
    qiime tools export --input-path relative-table-ASV.qza --output-path relative-table-ASV

    #convert to tab separated text file "rel-table-ASV.tsv"
    biom convert -i relative-table-ASV/feature-table.biom \
        -o rel-table-ASV.tsv --to-tsv

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}
