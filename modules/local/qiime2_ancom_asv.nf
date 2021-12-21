process QIIME2_ANCOM_ASV {
    tag "${table.baseName}"
    label 'process_medium'
    label 'single_cpu'
    label 'process_long'
    label 'error_ignore'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.8"

    input:
    tuple path(metadata), path(table)

    output:
    path("ancom/*")     , emit: ancom
    path "versions.yml" , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    qiime composition add-pseudocount \
        --i-table ${table} \
        --o-composition-table comp-${table}
    qiime composition ancom \
        --i-table comp-${table} \
        --m-metadata-file ${metadata} \
        --m-metadata-column ${table.baseName} \
        --o-visualization comp-${table.baseName}.qzv
    qiime tools export --input-path comp-${table.baseName}.qzv \
        --output-path ancom/Category-${table.baseName}-ASV

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
    END_VERSIONS
    """
}
