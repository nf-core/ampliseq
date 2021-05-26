// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process FORMAT_TAXRESULTS {
        label 'process_low'
        publishDir "${params.outdir}",
            mode: params.publish_dir_mode,
            saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'dada2', publish_id:'') }

        conda (params.enable_conda ? "pandas=1.1.5" : null)
        if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
            container "https://depot.galaxyproject.org/singularity/pandas:1.1.5"
        } else {
            container "quay.io/biocontainers/pandas:1.1.5"
        }

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