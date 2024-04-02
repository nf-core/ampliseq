process SIDLE_TRIM {
    tag "$meta.region,$meta.region_length"
    label 'process_single'

    container 'nf-core/pipesidle:0.1.0-beta'

    input:
    tuple val(meta), path(table), path(seq)

    output:
    tuple val(meta), path("*_table.qza")    , emit: table
    tuple val(meta), path("*_rep-seqs.qza") , emit: seq
    path "versions.yml"                     , emit: versions

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
    def length = "${meta.region_length}"
    """
    # https://q2-sidle.readthedocs.io/en/latest/read_preparation.html#dada2
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime sidle trim-dada2-posthoc \\
        --i-table ${table} \\
        --i-representative-sequences ${seq} \\
        --p-trim-length $length \\
        --o-trimmed-table ${prefix}_${length}_table.qza \\
        --o-trimmed-representative-sequences ${prefix}_${length}_rep-seqs.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        qiime2 plugin sidle: \$( qiime sidle --version | sed 's/ (.*//' | sed 's/.*version //' )
        q2-sidle: \$( qiime sidle --version | sed 's/.*version //' | sed 's/)//' )
    END_VERSIONS
    """
}
