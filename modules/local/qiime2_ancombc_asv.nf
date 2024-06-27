process QIIME2_ANCOMBC_ASV {
    tag "${table.baseName} ${formula}"
    label 'process_medium'
    label 'single_cpu'
    label 'process_long'
    label 'error_ignore'

    container "qiime2/core:2023.7"

    input:
    tuple path(metadata), path(table), val(formula)

    output:
    path("da_barplot/*")   , emit: da_barplot
    path("differentials/*"), emit: differentials
    path("*.qza")          , emit: qza
    path("*.qzv")          , emit: qzv
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args        = task.ext.args ?: ''
    def args2       = task.ext.args2 ?: ''
    def formula     = formula ?: "${table.baseName}"
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime composition ancombc \\
        --i-table "${table}" \\
        --m-metadata-file "${metadata}" \\
        $args \\
        --p-formula '${formula}' \\
        --o-differentials "${formula}.differentials.qza" \\
        --verbose
    qiime tools export \\
        --input-path "${formula}.differentials.qza" \\
        --output-path "differentials/Category-${formula}-ASV"

    # Generate tabular view of ANCOM-BC output
    qiime composition tabulate \\
        --i-data "${formula}.differentials.qza" \\
        --o-visualization "${formula}.differentials.qzv"
    qiime tools export \\
        --input-path "${formula}.differentials.qzv" \\
        --output-path "differentials/Category-${formula}-ASV"

    # Generate bar plot views of ANCOM-BC output
    qiime composition da-barplot \\
        --i-data "${formula}.differentials.qza" \\
        $args2 \\
        --o-visualization "${formula}.da_barplot.qzv"
    qiime tools export --input-path "${formula}.da_barplot.qzv" \\
        --output-path "da_barplot/Category-${formula}-ASV"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
