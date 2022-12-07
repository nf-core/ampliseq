process QIIME2_TRAIN {
    tag "${meta.FW_primer}-${meta.RV_primer}"
    label 'process_high'
    label 'single_cpu'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2022.8"

    input:
    tuple val(meta), path(qza)

    output:
    path("*-classifier.qza"), emit: qza
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    #Train classifier
    qiime feature-classifier fit-classifier-naive-bayes \\
        --i-reference-reads ${meta.FW_primer}-${meta.RV_primer}-ref-seq.qza \\
        --i-reference-taxonomy ref-taxonomy.qza \\
        --o-classifier ${meta.FW_primer}-${meta.RV_primer}-classifier.qza \\
        --quiet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
