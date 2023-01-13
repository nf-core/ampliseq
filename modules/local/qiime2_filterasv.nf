process QIIME2_FILTERASV {
    tag "${category}"
    label 'process_low'

    container "quay.io/qiime2/core:2022.8"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    tuple path(metadata), path(table), val(category)

    output:
    path("*.qza")       , emit: qza
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    qiime feature-table filter-samples \\
        --i-table ${table} \\
        --m-metadata-file ${metadata} \\
        --p-where \"${category}<>\'\'\" \\
        --o-filtered-table ${category}.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
