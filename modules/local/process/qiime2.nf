// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process QIIME2_INASV {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    //tuple val(meta), path(asv)
    path(asv)
    
    output:
    //tuple val(meta), path("table.qza"), emit: qza
    path("table.qza"), emit: qza
    path "*.version.txt"       , emit: version

    script:
    def software      = getSoftwareName(task.process)
    """
    echo -n "#OTU Table" | cat - "$asv" > biom-table.txt
    biom convert -i biom-table.txt -o table.biom --table-type="OTU table" --to-hdf5
    qiime tools import \
        --input-path table.biom \
        --type 'FeatureTable[Frequency]' \
        --input-format BIOMV210Format \
        --output-path table.qza
    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}

process QIIME2_INSEQ {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    //tuple val(meta), path(seq)
    path(seq)
    
    output:
    //tuple val(meta), path("rep-seqs.qza"), emit: qza
    path("rep-seqs.qza"), emit: qza
    path "*.version.txt"       , emit: version

    script:
    def software      = getSoftwareName(task.process)
    """
    qiime tools import \
        --input-path "$seq" \
        --type 'FeatureData[Sequence]' \
        --output-path rep-seqs.qza
    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}

process QIIME2_INTAX {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "quay.io/qiime2/core:2021.2"
    } else {
        container "quay.io/qiime2/core:2021.2"
    }

    input:
    path(tax)
    
    output:
    //tuple val(meta), path("taxonomy.qza"), emit: qza
    path "*.version.txt"       , emit: version

    script:
    def software      = getSoftwareName(task.process)
    """
    #requires: Feature_ID\tTaxon, see https://forum.qiime2.org/t/dada2-taxonomy-into-qiime2/15369/2
    #qiime tools import \
    #    --type 'FeatureData[Taxonomy]' \
    #    --input-path $tax \
    #    --output-path taxonomy.qza
    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}