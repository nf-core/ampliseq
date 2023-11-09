process QIIME2_BARPLOT {
    label 'process_low'

    container "qiime2/core:2023.7"

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
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
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
