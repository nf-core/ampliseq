process SIDLE_ALIGN {
    tag "$meta.region"
    label 'process_medium'

    container 'nf-core/pipesidle:0.1.0-beta'

    input:
    tuple val(meta), path(kmers), path(seq)

    output:
    tuple val(meta), path("*rep-seqs_align-map.qza"), emit: aligned_map
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Sidle in QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.region}"
    def primerfw = "${meta.fw_primer}"
    def primerrv = "${meta.rv_primer}"
    """
    # https://q2-sidle.readthedocs.io/en/latest/reconstruction.html#regional-alignment
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime sidle align-regional-kmers \\
        --p-n-workers $task.cpus \\
        --i-kmers ${kmers} \\
        --i-rep-seq ${seq} \\
        --p-region ${meta.region} \\
        $args \\
        --o-regional-alignment ${prefix}_rep-seqs_align-map.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        qiime2 plugin sidle: \$( qiime sidle --version | sed 's/ (.*//' | sed 's/.*version //' )
        q2-sidle: \$( qiime sidle --version | sed 's/.*version //' | sed 's/)//' )
    END_VERSIONS
    """
}
