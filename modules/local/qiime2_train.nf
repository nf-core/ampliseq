process QIIME2_TRAIN {
    tag "${meta.FW_primer}-${meta.RV_primer}"
    label 'process_high'
    label 'single_cpu'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    tuple val(meta), path(qza)

    output:
    path("*-classifier.qza"), emit: qza
    path "versions.yml"     , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

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
