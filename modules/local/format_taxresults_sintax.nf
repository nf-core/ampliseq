process FORMAT_TAXRESULTS_SINTAX {
    label 'process_low'

    conda "conda-forge::python=3.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'quay.io/biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(taxtable)
    path(fastafile)
    val(outfile)
    val(taxlevels_input)

    output:
    path(outfile)      , emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def taxlevels = taxlevels_input ? taxlevels_input : "Kingdom,Phylum,Class,Order,Family,Genus,Species"
    """
    convert_sintax_output.py -i $taxtable -f $fastafile -o $outfile -t $taxlevels

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
