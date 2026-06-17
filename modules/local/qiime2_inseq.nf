process QIIME2_INSEQ {
    tag "${seq}"
    label 'process_low'

    conda "${moduleDir}/envs/rachis-qiime2-linux-64-conda.yml"
    container "qiime2/qiime2:2026.4"

    input:
    path(seq)

    output:
    path("rep-seqs.qza"), emit: qza
    path "versions.yml", emit: versions_qiime2_inseq, topic: versions

    script:
    """
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime tools import \\
        --input-path "$seq" \\
        --type 'FeatureData[Sequence]' \\
        --output-path rep-seqs.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
