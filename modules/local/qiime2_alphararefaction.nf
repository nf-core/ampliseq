// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process QIIME2_ALPHARAREFACTION {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.4"

    input:
    path(metadata)
    path(table)
    path(tree)
    path(stats)

    output:
    path("alpha-rarefaction/*"), emit: rarefaction
    path "*.version.txt"       , emit: version

    script:
    def software     = getSoftwareName(task.process)
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    maxdepth=\$(count_table_minmax_reads.py $stats maximum 2>&1)

    #check values
    if [ \"\$maxdepth\" -gt \"75000\" ]; then maxdepth=\"75000\"; fi
    if [ \"\$maxdepth\" -gt \"5000\" ]; then maxsteps=\"250\"; else maxsteps=\$((maxdepth/20)); fi
    qiime diversity alpha-rarefaction  \
        --i-table ${table}  \
        --i-phylogeny ${tree}  \
        --p-max-depth \$maxdepth  \
        --m-metadata-file ${metadata}  \
        --p-steps \$maxsteps  \
        --p-iterations 10  \
        --o-visualization alpha-rarefaction.qzv
    qiime tools export --input-path alpha-rarefaction.qzv  \
        --output-path alpha-rarefaction

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}
