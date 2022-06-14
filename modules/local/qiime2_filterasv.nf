process QIIME2_FILTERASV {
    tag "${category}"
    label 'process_low'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.8"

    input:
    tuple path(metadata), path(table), val(category)

    output:
    path("*.qza")       , emit: qza
    path "versions.yml" , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    qiime feature-table filter-samples \
        --i-table ${table} \
        --m-metadata-file ${metadata} \
        --p-where \"${category}<>\'\'\" \
        --o-filtered-table ${category}.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
    END_VERSIONS
    """
}
