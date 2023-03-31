process FILTER_COI_ASV {
    tag "${fasta}"
    label 'process_low'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2022.8"

    input:
    path(fasta)
    path(asv)
    path(dada2stats)

    output:
    path( "ASV_filtered.table.tsv" )  , emit: asv
    path( "ASV_filtered.fna" ) , emit: fasta
    path( "ASV_filtered.list" ) , emit: list
    path( "stats.filt.tsv" ), emit: stats
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    filt_InsBiom_asv.py -f ${fasta} -t ${asv} -p ASV
    filt_coi_stats.r ${dada2stats} ASV_filtered.table.tsv

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        R: \$(R --version | sed -n 1p | sed 's/R version //g' | sed 's/\\s.*\$//')
        END_VERSIONS
    
    """

}