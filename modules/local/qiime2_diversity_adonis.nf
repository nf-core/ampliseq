process QIIME2_DIVERSITY_ADONIS {
    tag "${core.baseName} - ${formula}"
    label 'process_low'

    container "qiime2/core:2023.7"

    input:
    tuple path(metadata), path(core), val(formula)

    output:
    path("adonis/*")     , emit: html
    path("*.qzv")        , emit: qzv
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime diversity adonis \\
        --p-n-jobs $task.cpus \\
        --i-distance-matrix ${core} \\
        --m-metadata-file ${metadata} \\
        --o-visualization ${core.baseName}_adonis.qzv \\
        $args \\
        --p-formula "${formula}"
    qiime tools export \\
        --input-path ${core.baseName}_adonis.qzv \\
        --output-path adonis/${core.baseName}-${formula}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
