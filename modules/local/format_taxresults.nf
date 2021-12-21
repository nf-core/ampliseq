process FORMAT_TAXRESULTS {
    label 'process_low'

    conda (params.enable_conda ? "pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    path(taxtable)
    path(taxtable_species)
    path(fastafile)

    output:
    path("ASV_tax.tsv")
    path("ASV_tax_species.tsv"), emit: tsv

    script:
    """
    add_full_sequence_to_taxfile.py $taxtable $fastafile
    add_full_sequence_to_taxfile.py $taxtable_species $fastafile
    """
}
