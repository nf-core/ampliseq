process QIIME2_CLASSIFY {
    tag "${repseq},${trained_classifier}"
    label 'process_high'

    container "qiime2/core:2023.7"

    input:
    path(trained_classifier)
    path(repseq)

    output:
    path("taxonomy.qza"), emit: qza
    path("taxonomy.tsv"), emit: tsv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
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
