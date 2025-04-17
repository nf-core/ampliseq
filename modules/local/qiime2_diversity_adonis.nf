process QIIME2_DIVERSITY_ADONIS {
    tag "${core.baseName}-${formula}"
    label 'process_low'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    tuple path(metadata), path(core), val(formula)

    output:
    path("adonis/*")     , emit: html
    path("*.qzv")        , emit: qzv
    path "versions.yml"  , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime diversity adonis \\
        --p-n-jobs $task.cpus \\
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
