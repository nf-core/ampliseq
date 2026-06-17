process FORMAT_TAXRESULTS_VSEARCH_LCA {
    label 'process_single'

    conda "conda-forge::python=3.9.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(lcaout)
    path(fastafile)
    val(outfile)
    val(taxlevels_input)

    output:
    path(outfile)      , emit: tsv
    path "versions.yml", emit: versions_format_taxresults_vsearch_lca, topic: versions

    script:
    def taxlevels = taxlevels_input ? taxlevels_input : "Kingdom,Phylum,Class,Order,Family,Genus,Species"
    """
    convert_vsearch_lca_output.py -i $lcaout -f $fastafile -o $outfile -t $taxlevels

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
