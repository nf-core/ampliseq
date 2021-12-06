// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process QIIME2_EXPORT_ABSOLUTE {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.4"

    input:
    path(table)
    path(repseq)
    path(taxonomy)
    val(tax_agglom_min)
    val(tax_agglom_max)

    output:
    path("rep-seq.fasta")              , emit: fasta
    path("feature-table.tsv")          , emit: tsv
    path("feature-table.biom")         , emit: biom
    path("feature-table-withtax.biom") , emit: biomtax
    path("seven_number_summary.tsv")   , emit: summary
    path("descriptive_stats.tsv")      , emit: descr
    path("abs-abund-table-*.tsv")      , emit: abundtable
    path "*.version.txt"               , emit: version

    script:
    def software     = getSoftwareName(task.process)
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    #produce raw count table in biom format "table/feature-table.biom"
    qiime tools export --input-path ${table}  \
        --output-path table
    cp table/feature-table.biom .
    
    # Export taxonomy to add to the biom file as metadata
    qiime tools export --input-path ${taxonomy}  \
        --output-path table
    # Change the header of the output file
    echo -e "#OTUID\\ttaxonomy\\tconfidence" > tax_table.tsv
    tail -n +2 table/taxonomy.tsv | sed 's/\\(.*\\);\\(.*\\)\$/\\1\\t\\2/' >> tax_table.tsv
    # Add the taxonomy to the biom file
    biom add-metadata -i feature-table.biom -o feature-table-withtax.biom --observation-metadata-fp tax_table.tsv --sc-separated taxonomy

    #produce raw count table "table/feature-table.tsv"
    biom convert -i table/feature-table.biom \
        -o feature-table.tsv  \
        --to-tsv

    #produce representative sequence fasta file "sequences.fasta"
    qiime feature-table tabulate-seqs  \
        --i-data ${repseq}  \
        --o-visualization rep-seqs.qzv
    qiime tools export --input-path rep-seqs.qzv  \
        --output-path representative_sequences
    cp representative_sequences/sequences.fasta rep-seq.fasta
    cp representative_sequences/*.tsv .

    ##on several taxa level
    array=(\$(seq ${tax_agglom_min} 1 ${tax_agglom_max}))
    for i in \${array[@]}
    do
        #collapse taxa
        qiime taxa collapse \
            --i-table ${table} \
            --i-taxonomy ${taxonomy} \
            --p-level \$i \
            --o-collapsed-table table-\$i.qza
        #export to biom
        qiime tools export --input-path table-\$i.qza \
            --output-path table-\$i
        #convert to tab separated text file
        biom convert \
            -i table-\$i/feature-table.biom \
            -o abs-abund-table-\$i.tsv --to-tsv
    done

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}
