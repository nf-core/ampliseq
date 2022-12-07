process QIIME2_INSEQ {
    tag "${seq}"
    label 'process_low'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2022.8"

    input:
    path(seq)

    output:
    path("rep-seqs.qza"), emit: qza
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
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
