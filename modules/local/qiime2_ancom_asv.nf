process QIIME2_ANCOM_ASV {
    tag "${table.baseName}"
    label 'process_medium'
    label 'single_cpu'
    label 'process_long'
    label 'error_ignore'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    tuple path(metadata), path(table)

    output:
    path("ancom/*")     , emit: ancom
    path "versions.yml" , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime composition add-pseudocount \\
        --i-table ${table} \\
        --o-composition-table comp-${table}
    qiime composition ancom \\
        --i-table comp-${table} \\
        --m-metadata-file ${metadata} \\
        --m-metadata-column ${table.baseName} \\
        --o-visualization comp-${table.baseName}.qzv
    qiime tools export --input-path comp-${table.baseName}.qzv \\
        --output-path ancom/Category-${table.baseName}-ASV

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
