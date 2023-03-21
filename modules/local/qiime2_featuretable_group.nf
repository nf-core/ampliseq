process QIIME2_FEATURETABLE_GROUP {
    tag "${category}"
    label 'process_low'

    container "quay.io/qiime2/core:2022.11"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple path(table), path(metadata), val(category)

    output:
    path("${category}.qza"), emit: qza
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    qiime feature-table filter-samples \\
        --i-table "${table}" \\
        --m-metadata-file "${metadata}" \\
        --p-where \"${category}<>\'\'\" \\
        --o-filtered-table "filtered_${category}.qza"

    qiime feature-table group \\
        --i-table "filtered_${category}.qza" \\
        --p-axis 'sample' \\
        --m-metadata-file "${metadata}" \\
        --m-metadata-column "${category}" \\
        --p-mode 'sum' \\
        --o-grouped-table "${category}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
