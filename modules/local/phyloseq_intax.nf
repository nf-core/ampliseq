process PHYLOSEQ_INTAX {
    label 'process_low'

    conda "conda-forge::pandas=1.1.5 conda-forge::python=3.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5':
        'biocontainers/pandas:1.1.5' }"

    input:
    path(tax_tsv)

    output:
    path( "*.tsv" )          , emit: tsv
    path "versions.yml"      , emit: versions

    script:
    """
    reformat_tax_for_phyloseq.py $tax_tsv reformat_$tax_tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
