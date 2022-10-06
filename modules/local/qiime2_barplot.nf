process QIIME2_BARPLOT {
    label 'process_low'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2022.8"

    input:
    path(metadata)
    path(table)
    path(taxonomy)
    val(setting)

    output:
    path("barplot${suffix}/*"), emit: folder
    path "versions.yml"       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    suffix = setting ? "_${table.baseName}" : ""
    def metadata_cmd = metadata ? "--m-metadata-file ${metadata}": ""
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    qiime taxa barplot  \
        --i-table ${table}  \
        --i-taxonomy ${taxonomy}  \
        ${metadata_cmd}  \
        --o-visualization taxa-bar-plots.qzv  \
        --verbose
    qiime tools export --input-path taxa-bar-plots.qzv  \
        --output-path barplot${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
    END_VERSIONS
    """
}
