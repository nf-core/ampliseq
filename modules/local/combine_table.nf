process COMBINE_TABLE {
    label 'process_low'

    conda "bioconda::bioconductor-biostrings=2.58.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-biostrings:2.58.0--r40h037d062_0' :
        'biocontainers/bioconductor-biostrings:2.58.0--r40h037d062_0' }"

    input:
    path(table)
    path(seq)
    path(tax)
    val(outfile)

    output:
    path("${outfile}")  , emit: tsv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    combine_table.r ${table} ${seq} ${tax}
    mv combined_ASV_table.tsv ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1 | sed -n 1p | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """
}
