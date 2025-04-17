process QIIME2_FEATURETABLE_GROUP {
    tag "${category}"
    label 'process_low'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    tuple path(table), path(metadata), val(category)

    output:
    path("${category}.qza"), emit: qza
    path "versions.yml"    , emit: versions

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
