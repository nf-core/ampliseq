process QIIME2_DIVERSITY_ADONIS {
    tag "${core.baseName} - ${params.qiime_adonis_formula}"
    label 'process_low'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2022.2"

    input:
    tuple path(metadata), path(core)

    output:
    path("adonis/*")     , emit: html
    path "versions.yml"  , emit: versions

    script:
    def args = task.ext.args ?: ''
    def formula = params.qiime_adonis_formula
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"
    export MPLCONFIGDIR="\${PWD}/HOME"

    qiime diversity adonis \\
        --p-n-jobs $task.cpus \\
        --i-distance-matrix ${core} \\
        --m-metadata-file ${metadata} \\
        --o-visualization ${core.baseName}_adonis.qzv \\
        $args \\
        --p-formula "${formula}"
    qiime tools export --input-path ${core.baseName}_adonis.qzv \\
        --output-path adonis/${core.baseName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
    END_VERSIONS
    """
}
