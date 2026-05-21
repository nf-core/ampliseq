process QIIME2_ANCOMBC2_ASV {
    tag "${table.baseName}-${formula_in}"
    label 'process_medium'
    label 'process_long'
    label 'error_ignore'

    conda "${moduleDir}/envs/rachis-qiime2-linux-64-conda.yml"
    container "qiime2/qiime2:2026.4"

    input:
    tuple path(metadata), path(table), val(formula_in)

    output:
    path("visualizer/*")         , emit: plot
    path("differentials/*")      , emit: differentials
    path("*.qza")                , emit: qza
    path("*.qzv")                , emit: qzv
    path "versions.yml"          , emit: versions_qiime2_ancombc2_asv, topic: versions

    script:
    def args        = task.ext.args ?: ''
    def formula     = formula_in ?: "${table.baseName}"
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime composition ancombc2 \\
        --i-table "${table}" \\
        --m-metadata-file "${metadata}" \\
        $args \\
        --p-fixed-effects-formula '${formula}' \\
        --o-ancombc2-output "${formula}.differentials.qza" \\
        --p-num-processes ${task.cpus}  \\
        --verbose
    qiime tools export \\
        --input-path "${formula}.differentials.qza" \\
        --output-path "differentials/Category-${formula}-ASV"

    # Generate tabular view of ANCOMBC2 output
    qiime composition tabulate \\
        --i-data "${formula}.differentials.qza" \\
        --o-visualization "${formula}.differentials.qzv"
    qiime tools export \\
        --input-path "${formula}.differentials.qzv" \\
        --output-path "differentials/Category-${formula}-ASV"

    # Generate bar plot views of ANCOMBC2 output
    qiime composition ancombc2-visualizer \\
        --i-data "${formula}.differentials.qza" \\
        --o-visualization "${formula}.visualizer.qzv"
    # 'qiime tools export' does not produce a valid html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
