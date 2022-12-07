process QIIME2_DIVERSITY_ADONIS {
    tag "${core.baseName} - ${formula}"
    label 'process_low'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2022.8"

    input:
    tuple path(metadata), path(core), val(formula)

    output:
    path("adonis/*")     , emit: html
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

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
