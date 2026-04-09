process FORMAT_TAXRESULTS_SINTAX {
    label 'process_single'

    conda "conda-forge::python=3.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(taxtable)
    path(fastafile)
    val(outfile)
    val(taxlevels_input)

    output:
    path(outfile)      , emit: tsv
    path "versions.yml", emit: versions

    script:
    def taxlevels = taxlevels_input ? taxlevels_input : "Kingdom,Phylum,Class,Order,Family,Genus,Species"
    def sintax_dbversion = params.sintax_ref_tax_custom ? 'user_supplied' : params.sintax_ref_databases[params.sintax_ref_taxonomy]["dbversion"]
    """
    convert_sintax_output.py -i $taxtable -f $fastafile -o $outfile -t $taxlevels -d \"${sintax_dbversion}\"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
