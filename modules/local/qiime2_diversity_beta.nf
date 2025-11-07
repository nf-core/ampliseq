process QIIME2_DIVERSITY_BETA {
    tag "${core.baseName}-${category}"
    label 'process_low'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    tuple path(metadata), path(core), val(category)

    output:
    path("beta_diversity/*"), emit: beta
    path("*.qzv")           , emit: qzv
    path "versions.yml"     , emit: versions

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
