process SIDLE_DBFILT {
    label 'process_low'

    container 'nf-core/pipesidle:0.1.0-beta'

    input:
    path(seq)
    path(tax)

    output:
    path("db_filtered_sequences.qza")     , emit: seq
    path("db_filtered_sequences_tax.qza") , emit: tax
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Sidle in QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    # https://q2-sidle.readthedocs.io/en/latest/database_preparation.html#filtering-the-database
    #pre-filtering should be very permissive!
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    # authors of SMURF recommend "--p-num-degenerates 3" for greengenes 13_8 database at 99%
    # the RESCRIPt formatted Silva 128 database is filtered to exclude sequences with more than 5 degenerates [3], [4]
    qiime rescript cull-seqs \\
        --p-n-jobs $task.cpus \\
        --i-sequences $seq \\
        $args \\
        --o-clean-sequences db_filtered_sequences.qza

    #filtering a greengenes database for features missing a phylum (p__;) or kingdom(k__;) designation.
    #CPU=1
    qiime taxa filter-seqs \\
        --i-sequences db_filtered_sequences.qza \\
        --i-taxonomy $tax \\
        $args2 \\
        --o-filtered-sequences db_filtered_sequences_tax.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        qiime2 rescript: \$( qiime rescript --version | sed 's/ (.*//' | sed 's/.*version //' )
    END_VERSIONS
    """
}
