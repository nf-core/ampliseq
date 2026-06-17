process QIIME2_DIVERSITY_BETA {
    tag "${core.baseName}-${category}"
    label 'process_low'

    conda "${moduleDir}/envs/rachis-qiime2-linux-64-conda.yml"
    container "qiime2/qiime2:2026.4"

    input:
    tuple path(metadata), path(core), val(category)

    output:
    path("beta_diversity/*"), emit: beta
    path("*.qzv")           , emit: qzv
    path "versions.yml"     , emit: versions_qiime2_diversity_beta, topic: versions

    script:
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
