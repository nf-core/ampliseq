process QIIME2_FEATURETABLE_GROUP {
    tag "${category}"
    label 'process_single'

    conda "${moduleDir}/envs/rachis-qiime2-linux-64-conda.yml"
    container "qiime2/qiime2:2026.4"

    input:
    tuple path(table), path(metadata), val(category)

    output:
    path("${category}.qza"), emit: qza
    path "versions.yml"    , emit: versions_qiime2_featuretable_group, topic: versions

    script:
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime feature-table filter-samples \\
        --i-table "${table}" \\
        --m-metadata-file "${metadata}" \\
        --p-where \"${category}<>\'\'\" \\
        --o-filtered-table "filtered_${category}.qza"

    qiime feature-table group \\
        --i-table "filtered_${category}.qza" \\
        --p-axis 'sample' \\
        --m-metadata-file "${metadata}" \\
        --m-metadata-column "${category}" \\
        --p-mode 'sum' \\
        --o-grouped-table "${category}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
