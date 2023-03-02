process QIIME2_DIVERSITY_ALPHA {
    tag "${core.baseName}"
    label 'process_low'

    container "quay.io/qiime2/core:2022.11"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple path(metadata), path(core)

    output:
    path("alpha_diversity/*"), emit: alpha
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    qiime diversity alpha-group-significance \\
        --i-alpha-diversity ${core} \\
        --m-metadata-file ${metadata} \\
        --o-visualization ${core.baseName}-vis.qzv
    qiime tools export \\
        --input-path ${core.baseName}-vis.qzv \\
        --output-path "alpha_diversity/${core.baseName}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
