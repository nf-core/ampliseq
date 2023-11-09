process QIIME2_INTREE {
    tag "${meta.id}:${meta.model}"
    label 'process_low'

    container "qiime2/core:2023.7"

    input:
    tuple val(meta), path(tree)

    output:
    path("tree.qza")   , emit: qza
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime tools import \\
        --type 'Phylogeny[Rooted]' \\
        --input-path $tree \\
        --output-path tree.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
