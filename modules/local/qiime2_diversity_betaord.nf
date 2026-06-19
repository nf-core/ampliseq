process QIIME2_DIVERSITY_BETAORD {
    tag "${core.baseName}"
    label 'process_low'

    conda "${moduleDir}/envs/rachis-qiime2-linux-64-conda.yml"
    container "qiime2/qiime2:2026.4"

    input:
    tuple path(metadata), path(core)

    output:
    path("beta_diversity/*"), emit: beta
    path("*.qzv")           , emit: qzv
    path "versions.yml"     , emit: versions_qiime2_diversity_betaord, topic: versions

    script:
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
