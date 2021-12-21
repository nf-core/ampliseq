process QIIME2_INTAX {
    tag "${tax}"
    label 'process_low'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.8"

    input:
    path(tax) //ASV_tax_species.tsv

    output:
    path("taxonomy.qza") , emit: qza
    path "versions.yml"  , emit: versions

    script:
    """
    parse_dada2_taxonomy.r $tax

    qiime tools import \
        --type 'FeatureData[Taxonomy]' \
        --input-format HeaderlessTSVTaxonomyFormat \
        --input-path tax.tsv \
        --output-path taxonomy.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
    END_VERSIONS
    """
}
