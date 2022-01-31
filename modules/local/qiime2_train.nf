process QIIME2_TRAIN {
    tag "${meta.FW_primer}-${meta.RV_primer}"
    label 'process_high'
    label 'single_cpu'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.8"

    input:
    tuple val(meta), path(qza)

    output:
    path("*-classifier.qza"), emit: qza
    path "versions.yml"    , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    #Train classifier
    qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads ${meta.FW_primer}-${meta.RV_primer}-ref-seq.qza \
        --i-reference-taxonomy ref-taxonomy.qza \
        --o-classifier ${meta.FW_primer}-${meta.RV_primer}-classifier.qza \
        --quiet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
    END_VERSIONS
    """
}
