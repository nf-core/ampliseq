process QIIME2_CLASSIFY {
    tag "${repseq},${trained_classifier}"
    label 'process_high'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.8"

    input:
    path(trained_classifier)
    path(repseq)

    output:
    path("taxonomy.qza"), emit: qza
    path("taxonomy.tsv"), emit: tsv
    path "versions.yml" , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    qiime feature-classifier classify-sklearn  \
        --i-classifier ${trained_classifier}  \
        --p-n-jobs ${task.cpus}  \
        --i-reads ${repseq}  \
        --o-classification taxonomy.qza  \
        --verbose
    qiime metadata tabulate  \
        --m-input-file taxonomy.qza  \
        --o-visualization taxonomy.qzv  \
        --verbose
    #produce "taxonomy/taxonomy.tsv"
    qiime tools export --input-path taxonomy.qza  \
        --output-path taxonomy
    qiime tools export --input-path taxonomy.qzv  \
        --output-path taxonomy
    cp taxonomy/taxonomy.tsv .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
    END_VERSIONS
    """
}
