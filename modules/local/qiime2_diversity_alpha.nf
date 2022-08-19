process QIIME2_DIVERSITY_ALPHA {
    tag "${core.baseName}"
    label 'process_low'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2022.2"

    input:
    tuple path(metadata), path(core)

    output:
    path("alpha_diversity/*"), emit: alpha
    path "versions.yml"      , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    qiime diversity alpha-group-significance \
        --i-alpha-diversity ${core} \
        --m-metadata-file ${metadata} \
        --o-visualization ${core.baseName}-vis.qzv
    qiime tools export --input-path ${core.baseName}-vis.qzv \
        --output-path "alpha_diversity/${core.baseName}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
    END_VERSIONS
    """
}
