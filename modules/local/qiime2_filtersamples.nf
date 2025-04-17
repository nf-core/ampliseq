process QIIME2_FILTERSAMPLES {
    tag "${filter}"
    label 'process_low'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    tuple path(metadata), path(table), val(filter)

    output:
    path("*.qza")       , emit: qza
    path "versions.yml" , emit: versions

    script:
    def args = task.ext.args ?: "--p-where \'${filter}<>\"\"\'"
    def prefix = task.ext.prefix ?: "${filter}"
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime feature-table filter-samples \\
        --i-table ${table} \\
        --m-metadata-file ${metadata} \\
        $args \\
        --o-filtered-table ${prefix}.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
