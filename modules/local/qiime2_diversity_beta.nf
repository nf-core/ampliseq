process QIIME2_DIVERSITY_BETA {
    tag "${core.baseName} - ${category}"
    label 'process_low'

    container "qiime2/core:2023.7"

    input:
    tuple path(metadata), path(core), val(category)

    output:
    path("beta_diversity/*"), emit: beta
    path("*.qzv")           , emit: qzv
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime diversity beta-group-significance \\
        --i-distance-matrix ${core} \\
        --m-metadata-file ${metadata} \\
        --m-metadata-column \"${category}\" \\
        --o-visualization ${core.baseName}-${category}.qzv \\
        --p-pairwise
    qiime tools export \\
        --input-path ${core.baseName}-${category}.qzv \\
        --output-path beta_diversity/${core.baseName}-${category}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
