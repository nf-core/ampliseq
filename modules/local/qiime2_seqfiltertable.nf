process QIIME2_SEQFILTERTABLE {
    tag "${repseq} filter by ${table}"
    label 'process_low'

    container "qiime2/core:2023.7"

    input:
    path(table)
    path(repseq)

    output:
    path("filtered-sequences.qza"), emit: qza
    path "versions.yml"           , emit: versions

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

    qiime feature-table filter-seqs \\
        --i-data $repseq \\
        --i-table $table \\
        --o-filtered-data filtered-sequences.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
