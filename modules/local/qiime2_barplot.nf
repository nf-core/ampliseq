process QIIME2_BARPLOT {
    label 'process_low'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    path(metadata)
    path(table)
    path(taxonomy)
    val(setting)

    output:
    path("barplot${suffix}/*"), emit: folder
    path "versions.yml"       , emit: versions

    script:
    suffix = setting ? "_${table.baseName}" : ""
    def metadata_cmd = metadata ? "--m-metadata-file ${metadata}": ""
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime taxa barplot  \\
        --i-table ${table}  \\
        --i-taxonomy ${taxonomy}  \\
        ${metadata_cmd}  \\
        --o-visualization taxa-bar-plots.qzv  \\
        --verbose
    qiime tools export \\
        --input-path taxa-bar-plots.qzv  \\
        --output-path barplot${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
