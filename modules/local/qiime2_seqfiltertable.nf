process QIIME2_SEQFILTERTABLE {
    tag "${repseq}-filter-by-${table}"
    label 'process_low'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    path(table)
    path(repseq)

    output:
    path("filtered-sequences.qza"), emit: qza
    path "versions.yml"           , emit: versions

    script:
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
