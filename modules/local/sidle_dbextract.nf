
process SIDLE_DBEXTRACT {
    tag "$meta.region,$meta.region_length"
    label 'process_medium'

    container 'nf-core/pipesidle:0.1.0-beta'

    input:
    tuple val(meta), path(table), path(seq), path(db_seq), path(db_tax)

    output:
    tuple val(meta), path("db_*_kmers.qza"), emit: kmers
    tuple val(meta), path("db_*_map.qza")  , emit: map
    path "versions.yml"                    , emit: versions

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
    # https://q2-sidle.readthedocs.io/en/latest/database_preparation.html#prepare-a-regional-database-for-each-primer-set
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    #extract sequences
    qiime feature-classifier extract-reads \\
        --p-n-jobs $task.cpus \\
        --i-sequences $db_seq \\
        $args \\
        --p-f-primer $primerfw \\
        --p-r-primer $primerrv \\
        --o-reads db_${prefix}.qza

    #prepare to be used in alignment
    qiime sidle prepare-extracted-region \\
        --p-n-workers $task.cpus \\
        --i-sequences db_${prefix}.qza \\
        --p-region "${prefix}" \\
        --p-fwd-primer $primerfw \\
        --p-rev-primer $primerrv \\
        --p-trim-length $length \\
        --o-collapsed-kmers db_${prefix}_${length}_kmers.qza \\
        --o-kmer-map db_${prefix}_${length}_map.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        qiime2 plugin sidle: \$( qiime sidle --version | sed 's/ (.*//' | sed 's/.*version //' )
        q2-sidle: \$( qiime sidle --version | sed 's/.*version //' | sed 's/)//' )
    END_VERSIONS
    """
}
