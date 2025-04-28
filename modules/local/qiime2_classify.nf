process QIIME2_CLASSIFY {
    tag "${repseq},${trained_classifier}"
    label 'process_high'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    path(trained_classifier)
    path(repseq)

    output:
    path("taxonomy.qza"), emit: qza
    path("taxonomy.tsv"), emit: tsv
    path "versions.yml" , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime feature-classifier classify-sklearn  \\
        --i-classifier ${trained_classifier}  \\
        --p-n-jobs ${task.cpus}  \\
        --i-reads ${repseq}  \\
        --o-classification taxonomy.qza  \\
        --verbose
    qiime metadata tabulate  \\
        --m-input-file taxonomy.qza  \\
        --o-visualization taxonomy.qzv  \\
        --verbose
    #produce "taxonomy/taxonomy.tsv"
    qiime tools export \\
        --input-path taxonomy.qza  \\
        --output-path taxonomy
    qiime tools export \\
        --input-path taxonomy.qzv  \\
        --output-path taxonomy
    cp taxonomy/taxonomy.tsv .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
