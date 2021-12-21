process QIIME2_EXTRACT {
    tag "${meta.FW_primer}-${meta.RV_primer}"
    label 'process_low'
    label 'single_cpu'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.8"

    input:
    tuple val(meta), path(database)

    output:
    tuple val(meta), path("*.qza"), emit: qza
    path "versions.yml"          , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    ### Import
    qiime tools import --type \'FeatureData[Sequence]\' \
        --input-path ${database[0]} \
        --output-path ref-seq.qza
    qiime tools import --type \'FeatureData[Taxonomy]\' \
        --input-format HeaderlessTSVTaxonomyFormat \
        --input-path ${database[1]} \
        --output-path ref-taxonomy.qza
    #Extract sequences based on primers
    qiime feature-classifier extract-reads \
        --i-sequences ref-seq.qza \
        --p-f-primer ${meta.FW_primer} \
        --p-r-primer ${meta.RV_primer} \
        --o-reads ${meta.FW_primer}-${meta.RV_primer}-ref-seq.qza \
        --quiet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
    END_VERSIONS
    """
}
