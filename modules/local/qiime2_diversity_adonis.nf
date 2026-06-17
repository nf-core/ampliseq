process QIIME2_DIVERSITY_ADONIS {
    tag "${core.baseName}-${formula}"
    label 'process_low'

    conda "${moduleDir}/envs/rachis-qiime2-linux-64-conda.yml"
    container "qiime2/qiime2:2026.4"

    input:
    tuple path(metadata), path(core), val(formula)

    output:
    path("adonis/*")     , emit: html
    path("*.qzv")        , emit: qzv
    path "versions.yml"  , emit: versions_qiime2_diversity_adonis, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    # more than 1 process is failing with QIIME2 2026.4.0!
    qiime diversity adonis \\
        --p-n-jobs 1 \\
        --i-distance-matrix ${core} \\
        --m-metadata-file ${metadata} \\
        --o-visualization ${core.baseName}_adonis.qzv \\
        $args \\
        --p-formula "${formula}"
    qiime tools export \\
        --input-path ${core.baseName}_adonis.qzv \\
        --output-path adonis/${core.baseName}-${formula}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
