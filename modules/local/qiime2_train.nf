process QIIME2_TRAIN {
    tag "${meta.FW_primer}-${meta.RV_primer}"
    label 'process_high'
    label 'single_cpu'

    container "qiime2/core:2023.7"

    input:
    tuple val(meta), path(qza)

    output:
    path("*-classifier.qza"), emit: qza
    path "versions.yml"     , emit: versions

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
