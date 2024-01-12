process QIIME2_EXTRACT {
    tag "${meta.FW_primer}-${meta.RV_primer}"

    container "qiime2/core:2023.7"

    input:
    tuple val(meta), path(database)

    output:
    tuple val(meta), path("*.qza"), emit: qza
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    ### Import
    qiime tools import \\
        --type \'FeatureData[Sequence]\' \\
        --input-path ${database[0]} \\
        --output-path ref-seq.qza
    qiime tools import \\
        --type \'FeatureData[Taxonomy]\' \\
        --input-format HeaderlessTSVTaxonomyFormat \\
        --input-path ${database[1]} \\
        --output-path ref-taxonomy.qza
    #Extract sequences based on primers
    qiime feature-classifier extract-reads \\
        --p-n-jobs ${task.cpus} \\
        --i-sequences ref-seq.qza \\
        --p-f-primer ${meta.FW_primer} \\
        --p-r-primer ${meta.RV_primer} \\
        $args \\
        --o-reads ${meta.FW_primer}-${meta.RV_primer}-ref-seq.qza \\
        --quiet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}
