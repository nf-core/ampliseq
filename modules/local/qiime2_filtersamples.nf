process QIIME2_FILTERSAMPLES {
    tag "${filter}"
    label 'process_low'

    container "qiime2/core:2023.7"

    input:
    tuple path(metadata), path(table), val(filter)

    output:
    path("*.qza")       , emit: qza
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
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
