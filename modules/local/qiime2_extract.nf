process QIIME2_EXTRACT {
    tag "${meta.FW_primer}-${meta.RV_primer}"

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    tuple val(meta), path(database)

    output:
    tuple val(meta), path("*.qza"), emit: qza
    path "versions.yml"          , emit: versions

    script:
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
