// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process QIIME2_DIVERSITY_CORE {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    path(metadata)
    path(table)
    path(tree)
    path(stats)

    output:
    path("diversity_core/*_pcoa_results.qza")   , emit: pcoa
    path("diversity_core/*_vector.qza")         , emit: vector
    path("diversity_core/*_distance_matrix.qza"), emit: distance
    path "*.version.txt"                        , emit: version
    path("*rarefaction.txt")                    , emit: depth

    script:
    def software     = getSoftwareName(task.process)
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    mindepth=\$(count_table_minmax_reads.py $stats minimum 2>&1)
    if [ \"\$mindepth\" -gt \"10000\" ]; then echo \$mindepth >\"Use the sampling depth of \$mindepth for rarefaction.txt\" ; fi
    if [ \"\$mindepth\" -lt \"10000\" -a \"\$mindepth\" -gt \"5000\" ]; then echo \$mindepth >\"WARNING The sampling depth of \$mindepth is quite small for rarefaction.txt\" ; fi
    if [ \"\$mindepth\" -lt \"5000\" -a \"\$mindepth\" -gt \"1000\" ]; then echo \$mindepth >\"WARNING The sampling depth of \$mindepth is very small for rarefaction.txt\" ; fi
    if [ \"\$mindepth\" -lt \"1000\" ]; then echo \$mindepth >\"WARNING The sampling depth of \$mindepth seems too small for rarefaction.txt\" ; fi

    qiime diversity core-metrics-phylogenetic \
        --m-metadata-file ${metadata} \
        --i-phylogeny ${tree} \
        --i-table ${table} \
        --p-sampling-depth \$mindepth \
        --output-dir diversity_core \
        --p-n-jobs-or-threads ${task.cpus} \
        --verbose

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}
