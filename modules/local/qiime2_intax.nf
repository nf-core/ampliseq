// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process QIIME2_INTAX {
    tag "${asv},${tax}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    path(asv) //DADA2_table.tsv
    path(tax) //ASV_tax_species.tsv
    
    output:
    path("taxonomy.qza") , emit: qza
    path "*.version.txt" , emit: version

    script:
    def software      = getSoftwareName(task.process)
    """
    parse_dada2_taxonomy.r $asv $tax

    qiime tools import \
        --type 'FeatureData[Taxonomy]' \
        --input-format HeaderlessTSVTaxonomyFormat \
        --input-path tax.tsv \
        --output-path taxonomy.qza

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}