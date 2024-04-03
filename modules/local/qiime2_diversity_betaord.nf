process QIIME2_DIVERSITY_BETAORD {
    tag "${core.baseName}"
    label 'process_low'

    container "qiime2/core:2023.7"

    input:
    tuple path(metadata), path(core)

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
    mkdir beta_diversity

    qiime emperor plot \\
        --i-pcoa ${core} \\
        --m-metadata-file ${metadata} \\
        --o-visualization ${core.baseName}-vis.qzv
    qiime tools export \\
        --input-path ${core.baseName}-vis.qzv \
        --output-path beta_diversity/${core.baseName}-PCoA

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
