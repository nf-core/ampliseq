process QIIME2_DIVERSITY_BETA {
    tag "${core.baseName} - ${category}"
    label 'process_low'

    container "quay.io/qiime2/core:2022.11"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple path(metadata), path(core), val(category)

    output:
    path("beta_diversity/*"), emit: beta
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

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
