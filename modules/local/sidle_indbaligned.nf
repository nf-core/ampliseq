process SIDLE_INDBALIGNED {
    label 'process_single'

    conda "${projectDir}/modules/local/envs/pipesidle-0-1-0-beta.yml"
    container 'nf-core/pipesidle:0.1.0-beta'

    input:
    path(seq)

    output:
    path("db_alignedsequences.qza"), emit: seq
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    # db_seq
    qiime tools import \\
        --input-path $seq \\
        --output-path db_alignedsequences.qza \\
        --type 'FeatureData[AlignedSequence]'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
