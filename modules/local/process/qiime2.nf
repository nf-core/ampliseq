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

process QIIME2_EXTRACT {
    tag "${meta.FW_primer}-${meta.RV_primer}"
    label 'process_low'
    label 'single_cpu'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    tuple val(meta), path(database)
    
    output:
    tuple val(meta), path("*.qza"), emit: qza
    path "*.version.txt"       , emit: version

    script:
    def software      = getSoftwareName(task.process)
    """
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

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}

process QIIME2_TRAIN {
    tag "${meta.FW_primer}-${meta.RV_primer}"
    label 'process_high'
    label 'single_cpu'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    tuple val(meta), path(qza)
    
    output:
    path("*-classifier.qza"), emit: qza
    path "*.version.txt"       , emit: version

    script:
    def software      = getSoftwareName(task.process)
    """
    #Train classifier
    qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads ${meta.FW_primer}-${meta.RV_primer}-ref-seq.qza \
        --i-reference-taxonomy ref-taxonomy.qza \
        --o-classifier ${meta.FW_primer}-${meta.RV_primer}-classifier.qza \
        --quiet

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}

process QIIME2_CLASSIFY {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    path(trained_classifier)
    path(repseq)
    
    output:
    path("taxonomy.qza"), emit: qza
    path("taxonomy.tsv"), emit: tsv
    path "*.version.txt"       , emit: version

    script:
    def software      = getSoftwareName(task.process)
    """
    qiime feature-classifier classify-sklearn  \
        --i-classifier ${trained_classifier}  \
        --p-n-jobs ${task.cpus}  \
        --i-reads ${repseq}  \
        --o-classification taxonomy.qza  \
        --verbose
    qiime metadata tabulate  \
        --m-input-file taxonomy.qza  \
        --o-visualization taxonomy.qzv  \
        --verbose
    #produce "taxonomy/taxonomy.tsv"
    qiime tools export --input-path taxonomy.qza  \
        --output-path taxonomy
    qiime tools export --input-path taxonomy.qzv  \
        --output-path taxonomy
    cp taxonomy/taxonomy.tsv .

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}