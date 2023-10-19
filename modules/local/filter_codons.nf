process FILTER_CODONS {
    tag "${fasta}"
    label 'process_low'

    conda "conda-forge::pandas=1.1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5':
        'biocontainers/pandas:1.1.5' }"

    input:
    path(fasta)
    path(asv)

    output:
    path( "ASV_codon_filtered.table.tsv"  ) , emit: asv, optional: true
    path( "ASV_codon_filtered.fna"        ) , emit: fasta
    path( "ASV_codon_filtered.list"       ) , emit: list
    path( "codon.filtered.stats.tsv"      ) , emit: stats, optional: true
    path( "versions.yml"                  ) , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def count_table = asv ? "-t ${asv}" : ""
    def make_stats_cmd = asv ? "filt_codon_stats.py ASV_codon_filtered.table.tsv" : ""
    """
    filt_codons.py -f ${fasta} ${count_table} -p ASV_codon ${args}
    $make_stats_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
