// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process QIIME2_FILTERTAXA {
    tag "taxa:${exclude_taxa};min-freq:${min_frequency};min-samples:${min_samples}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    path(table)
    path(repseq)
    path(taxonomy)
    val(min_frequency)
    val(min_samples)
    val(exclude_taxa)

    output:
    path("filtered-table.qza"), emit: asv
    path("filtered-table.tsv"), emit: tsv
    path("filtered-sequences.qza"), emit: seq
    path "*.version.txt"       , emit: version

    script:
    def software     = getSoftwareName(task.process)
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    if ! [ \"${exclude_taxa}\" = \"none\" ]; then
        #filter sequences
        qiime taxa filter-seqs \
            --i-sequences ${repseq} \
            --i-taxonomy ${taxonomy} \
            --p-exclude ${exclude_taxa} --p-mode contains \
            --o-filtered-sequences tax_filtered-sequences.qza
        #filter abundance table
        qiime taxa filter-table \
            --i-table ${table} \
            --i-taxonomy ${taxonomy} \
            --p-exclude ${exclude_taxa} --p-mode contains \
            --o-filtered-table tax_filtered-table.qza
        filtered_table="tax_filtered-table.qza"
        filtered_sequences="tax_filtered-sequences.qza"
    else
        filtered_table=${table}
        filtered_sequences=${repseq}
    fi
    qiime feature-table filter-features \
        --i-table \$filtered_table \
        --p-min-frequency ${min_frequency} \
        --p-min-samples ${min_samples} \
        --o-filtered-table filtered-table.qza

    qiime feature-table filter-seqs \
        --i-data \$filtered_sequences \
        --i-table filtered-table.qza \
        --o-filtered-data filtered-sequences.qza

    #produce raw count table in biom format "table/feature-table.biom"
    qiime tools export --input-path filtered-table.qza  \
        --output-path table
    #produce raw count table
    biom convert -i table/feature-table.biom \
        -o table/feature-table.tsv  \
        --to-tsv
    cp table/feature-table.tsv filtered-table.tsv

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}
